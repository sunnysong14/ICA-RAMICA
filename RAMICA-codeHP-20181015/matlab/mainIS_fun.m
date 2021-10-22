function [Ave_amari, Std_amari, Amari, ...
    S, A, X_pn, X, RTS_est, Ainv_est] = ...
    mainIS_fun(CaseScen,nS,P_noise,StrICA_lst,Seed,CRImg_lst)
% USAGE
%   Main function of Image Separation of Caltech256 dataset on scen1 with
%   multiple noise levels, using multiple ICAs. 
% 
% INPUT
%   CaseScen: NOW only investigate scen1
%   Seed: Control both mixing matrices and chosen source images.
%   CRImg: Cell that records the row/column of the 2D image as the random
%   vector. The value is 'c' by default.
% 
% Liyan for AAAI'16 on 08/22/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSeed = length(Seed);

%%%% dir
dirData = ['..',filesep,'data',filesep,'Caltech256',filesep,'nS', num2str(nS), filesep];
to_dir = ['..', filesep, 'results', filesep, 'exp2_IS', filesep, 'nS', num2str(nS), filesep];
if ~exist(to_dir, 'file')
    mkdir(to_dir); 
end

%%%% debug: convert a single str to cell
if ischar(StrICA_lst)
    StrICA_lst = cellstr(StrICA_lst); 
end
if ischar(CRImg_lst)
    CRImg_lst = cellstr(CRImg_lst); 
end

% Caltech256 is 256x256
nRow = 64; %dimRVec
nCol = 256; %nSmp

% Special for RAMICA: if both c-&r- wise exit -> add one $StrICA_lst$
if sum(strcmpi(StrICA_lst,'ours')) && length(CRImg_lst)==2
    StrICA_lst{end+1} = 'ours'; %Add extra RAMICA at the end
end
nICA = length(StrICA_lst); %After rearrange.


%%%% Main Run %%%%
for sc = 1 : length(CaseScen)
    caseScen = CaseScen(sc);
    
    for pn = 1 : length(P_noise)
        p_noise = P_noise(pn);
        prename = ['scen', num2str(caseScen), '_nS', num2str(nS), ...
            '_pN',strrep(num2str(p_noise),'.','_'),'_nSd',num2str(nSeed),'_'];
        
        Amari = -ones(nICA, nSeed); %init
        
        for ic = 1 : nICA
            strICA = lower(StrICA_lst{ic});
            name = [prename, strICA]; 
            crImg = 'Not-used-for-this-ICA'; %overwritten by RAMICA.
            
            %%%% special for RAMICA
            if strcmpi(strICA, 'ours')
                if ic==nICA && length(CRImg_lst)==2 %Second RAMICA
                    crImg = CRImg_lst{2};
                else
                    crImg = CRImg_lst{1};
                end
                name = [name, upper(crImg)];
            end
            
            % flnm
            flnm = [to_dir, name];
            
            if exist([flnm, '.mat'], 'file')
                load(flnm,'AmariICA');

            else
                AmariICA = -ones(1, nSeed); %init
                for ss = 1 : nSeed
                    seed = Seed(ss);
                    
                    %%%% Get the Data %%%%
                    
                    %%%% load the source data in vector format
                    nmData = num2str(seed);
                    flnmData = [dirData, nmData];
                    imgdata = load(flnmData); 
                    scenstr = ['S_scen',num2str(caseScen)];
                    S_org = imgdata.(scenstr);
                    %load >> S_scen1, S_scen2, S_scen3, S_scen4, S_scen5.
                    
                    %%%% convert to data tensor: X_pn and TX_pn
                    [X_pn,X,S,A,Noise_cov] = obtainICA_IS(S_org,p_noise,seed);
                    %S: figure,for k=1:4,subplot(2,2,k),imshow(reshape(S(k, :),256,256),[]); end
                    %X_pn: figure,for k=1:4,subplot(2,2,k),imshow(reshape(X_pn(k, :),256,256),[]); end
                    
                    % TX_pn: column-wise by default
                    TX_pn = zeros(nS, nRow, nCol);
                    for t = 1 : nCol
                        TX_pn(:,:,t) = X_pn(:, ((t-1)*nRow+1):(t*nRow));
                    end
                    %Rest for row-wise
                    if strcmpi(strICA, 'ours') && strcmpi(crImg, 'r')
                        TX_pn_r = zeros(nS,nCol,nRow);
                        for kk = 1 : nS
                            TX_pn_r(kk,:,:)=squeeze(TX_pn(kk,:,:))'; %transpose
                        end
                        TX_pn = TX_pn_r; %replace
                    end
                    %figure,for k=1:4,subplot(2,2,k),imshow(squeeze(TX_pn(k,:,:)),[]);end
                    
                    
                    %%%% [core] ICA %%%%
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
                            error('Error: undefined ICA method.');
                    end%switch
                    
                    
                    %%%% Amari %%%%
                    if size(Ainv_est,1) < size(A,1)
                        amari = -1; 
                    else
                        P = Ainv_est * A;
                        amari = Amari_index_ISA(P, ones(1,size(S,1)), 'uniform', 2);
                    end%if                                      
                end%ss
                
                %%%% save
                save(flnm,'strICA','p_noise','Seed','crImg','AmariICA', ...
                    'S','A','X_pn','X','RTS_est','Ainv_est');                
            end%if
            Amari(ic,:) = AmariICA; %assign           
        end%ic
        
        
        %%%% PF Ave %%%%
        Ave_amari = -ones(size(Amari,1), 1);
        Std_amari = -ones(size(Amari,1), 1);
        
        for ic = 1 : size(Amari, 1)
            amari_arr = Amari(ic, :);
            pos_ = (amari_arr ~= -1);
            
            Ave_amari(ic) = mean(amari_arr(pos_));
            Std_amari(ic) = std(amari_arr(pos_));
        end%ic
        
    end%pn
end%sc

end%function
