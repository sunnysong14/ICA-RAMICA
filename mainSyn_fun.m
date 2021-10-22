function [para_arr_real, Ave_amari, Std_amari, Amari, ...
    S, A, X_pn, X, RTS_est, Ainv_est] = ...
    mainSyn_fun(strDistr,para_arr,pattern_str,...
    nS,nDim,nSmp,P_noise,StrICA_lst,Seed,tag_revRAMICA)
% USAGE
%   Main function of the our experiment on BSS on synthetic data with
%   multiple noise levels, multiple ICAs, and across multiple seeds. 
%   It deals with both row and column-wise of our method assigned by the
%   input 'tag_revRAMICA'.
% 
% Liyan for AAAI'16 on 09/06/2016; 
%   o 09/12/2016 add "cRAMICA", and current version is "rRAMICA": change the 
%   codes as small as possible. 
% 
% TODO
%   Better to be able to run RAMICA's two versions togther by changing the
%   codes after deadline. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Debug for RAMICA-reverse
% If no RAMICA or not-exsit
if sum(strcmpi('ours', StrICA_lst))==0 || ~exist('tag_revRAMICA','var')
    tag_revRAMICA = 0; %@VIP to set to 0. Used to justify the condition
end

%%%% rslt dir
to_dir = ['..',filesep,'Rslt',filesep,'exp1_syn',filesep];
if ~exist(to_dir, 'file')
    mkdir(to_dir); 
end

nSeed = length(Seed);
nICA = length(StrICA_lst);

% prename of source tensor
prenameS = [strDistr, '_Par', regexprep(num2str(para_arr),'\s+','_'), ...
    '_PT', pattern_str, ...
    '_nS', num2str(nS), '_nD', num2str(nDim), '_nSp', num2str(nSmp)];


%%%% MAIN LOOP %%%%
for pn = 1 : length(P_noise)
    p_noise = P_noise(pn);
    
    % prename with noise and seed
    prename = [prenameS,'_pN',strrep(num2str(p_noise),'.','_'),'_nSd',num2str(length(Seed)),'_'];

    % init
    Amari = -ones(nICA, nSeed); 

    for ic = 1 : nICA
        strICA = lower(StrICA_lst{ic});
        
        % flnm
        name = [prename, strICA]; 
        if tag_revRAMICA %for reverse RAMICA
            name = [name, 'Rev'];
        end
        flnm = [to_dir, name];
        
        %%%% core
        if exist([flnm, '.mat'], 'file')
            load(flnm,'AmariICA',...
                'S', 'A', 'X_pn', 'X', 'RTS_est', 'Ainv_est');
            
        else
            AmariICA = -ones(1, nSeed); %init            
            
            for ss = 1 : nSeed
                seed = Seed(ss);

                %%%% Generate Source and Mixture Tensor %%%%
                [TS,S,para_arr_real] = genRandData(strDistr,para_arr,nS,nDim,nSmp,seed,pattern_str);
                
                %%%% if reverse RAMICA
                if strcmpi(strICA,'ours') && tag_revRAMICA
                    nDim_rev = nSmp;
                    nSmp_rev = nDim;
                    % reverse TS
                    TS_rev = zeros(nS,nDim_rev,nSmp_rev);
                    for iis = 1 : nS
                        TS_rev(iis,:,:) = squeeze(TS(iis,:,:))';
                    end
                    % keep and save the original and replace with the revised TS_rev. 
                    % Bec. I want the least changed code by keeping using TS.
                    TS_bfRev = TS; 
                    TS = TS_rev; %@Overwrite TS
                    
                else
                    nDim_rev = nDim;
                    nSmp_rev = nSmp;
                end                
                
                %%%% generation matrix
                A = genA(nS, seed);
                TX = ttm(tensor(TS), A, 1); 
                TX = TX.data;
                
                %%%% mixture tensor
                X = zeros(nS, nDim*nSmp);
                for i_=1 : nS
                    iTX = TX(i_,:,:);
                    X(i_,:) = iTX(:)';
                end
                
                %%%% inject data noise
                if p_noise > 0 
                    [X_pn, Noise_cov] = injctGauNoise(X, p_noise, seed);
                    % Stack back to TX_pn
                    TX_pn = zeros(nS, nDim_rev, nSmp_rev);
                    for t = 1 : nSmp_rev
                        TX_pn(:,:,t) = X_pn(:, ((t-1)*nDim_rev+1):(t*nDim_rev));
                    end
                    
                elseif p_noise == 0 
                    X_pn = X;
                    TX_pn = TX;
                    Noise_cov = 0;
                    
                else
                    error('$p_noise=\sigma^2$ should be no smaller than 0. ');
                end%if
                
                
                %%%% Run ICA Methods %%%%
                switch lower(strICA)
                    case 'fastica' % 'pow3' (default)
                        [RTS_est, ~, Ainv_est] = fastica(X);
                        
                    case 'jade'
                        Ainv_est = jadeR(X);
                        RTS_est = Ainv_est * X;
                        
                    case 'infomax'
                        Ainv_est = infomaxICA(X);
                        RTS_est = Ainv_est * X;
                        
                    case lower('ours')
                        [Ainv_est, TS_est] = mJADE(TX_pn);
                        RTS_est = TS_est;
                        
                    otherwise
                        error('Error: not defined ICA method.');                        
                end%switch                
                
                
                %%%% PF of Amari and SINR_LOSS %%%%
                if size(Ainv_est,1)<size(A,1) || sum(sum(isnan(Ainv_est)))>0 %NaN
                    amari = -1;                    
                else
                    amari = AmariIndex(A,Ainv_est);
                    %should change to use Amari_index_ISA()
                end%if
                AmariICA(ss) = amari; %assign
            end%ss
            
            
            %%%% Save Rslt %%%%
            save(flnm, 'AmariICA',...
                'strDistr', 'para_arr_real', 'pattern_str', 'nS', 'nDim', 'nSmp', ...
                'strICA', 'p_noise', 'Seed', 'AmariICA', ...
                'S', 'A', 'X_pn', 'X', 'TX_pn', 'TX', 'RTS_est', 'Ainv_est');
            %%%% extra for reverse RAMICA
            if tag_revRAMICA 
                save(flnm, 'tag_revRAMICA', 'nDim_rev', 'nSmp_rev', 'TS_bfRev','-append');
            end

        end%if-load
        
        %%%% Assign
        Amari(ic,:) = AmariICA;        
    end%ic

    
    %%%% PF Ave %%%%
    % Note to ignore the INVALID rslts
    Ave_amari = -ones(size(Amari, 1), 1);
    Std_amari = -ones(size(Amari, 1), 1);    

    for ic = 1 : size(Amari, 1)
        amari_arr = Amari(ic, :);
        pos_ = (amari_arr ~= -1);
        
        % Amari
        Ave_amari(ic) = mean(amari_arr(pos_));
        Std_amari(ic) = std(amari_arr(pos_));
    end%ic
    
end%pn

end%function
