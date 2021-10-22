function [TData, Data, para_arr_real] = genRandData(strDistr, para_arr, nS, nDim, nSmp, seed, pattern_str)
% Usage:
%   Generate synthetic tensor of size (nDim x nSmp x nS). 
%   [09/12/2016] It generates the [row-wise] data, and structures in row is
%   retained. Because, refer to the figure, the synthetic source is column
%   spikes and the rows repeats the pattern, and row is like the 'mother',
%   i.e. random vector.
% 
% INPUT PARAMTER
%   para_arr -- Mat of parameters of strDistr. If blank, use the default
%   values.
% 
% OUTPUT PARAMETER
%   TData -- Tensor data of size nS x nDim x nSmp.
%   Data -- Vector data of size nS x (nDim x nSmp).
%   para_arr_real -- Update the real distribution parameters. Used in "mainSyn()"
% 
% Liyan 08/30/2016 for MaICA-aaai16; [09/09/2016] Each source matrix (QxT)
% generation: force correlation among them. [09/10/2016] Synthetic
% generation: x_base^power.

% init
TData = zeros(nS, nDim, nSmp); %Note: Dif from write-up for convenience.
Data = zeros(nS, nDim*nSmp);

% contral random
rng(seed)

%% Distribution
switch lower(strDistr)
    case lower('pearsrnd') % two parameters
        meanVL = 0; stdLV = 1;
        % parameters:
        if isempty(para_arr)
            skewVL = 1; kurtLV = 1+3; %@Gaussian = 3
        elseif length(para_arr)==1
            skewVL = 1; kurtLV = para_arr+3;
        elseif length(para_arr)>=2
            skewVL = para_arr(1); kurtLV = para_arr(2)+3;
        end
        para_arr_real(1)=skewVL; para_arr_real(2)=kurtLV;
        
        S_t1 = pearsrnd(meanVL, stdLV, skewVL, kurtLV, nS, nSmp);
        %--------------------------------------------------------------   

    case lower('unif') % no parameter
        % skew = 0; % excess kurtosis = -6/5
        para_arr_real = [];
        
        S_t1 = (rand(nS, nSmp) - 0.5) * sqrt(12);
        %-------------------------------------------------------------- 

    case lower('tdistr') % one parameter
        % Variance = paraDf / (paraDf-2) (must have paraDf >= 3 in order for variance assumptions on S to be met).
        % skew = 0 (undefined if paraDf <= 3). kurtosis = 6/(paraDf-4)
        % parameters:
        if isempty(para_arr)
            df = 5;
        elseif length(para_arr) >= 1
            df = para_arr(1);
        end
        para_arr_real = df;
        
        S_t1 = trnd(df, nS, nSmp) / sqrt(df / (df - 2));
        %--------------------------------------------------------------

    case 'binomial'     % one parameter
        % mean = p; var = p(1-p); skew = (1-2p)/sqrt(p(1-p)); 
        % excess kurtosis = (1-6*p*(1-p))/(p*(1-p));
        % parameters:
        if isempty(para_arr)
            p = 0.5;
        elseif length(para_arr) >= 1
            p = para_arr(1);
        end
        para_arr_real = p;

        S_t1 = (binornd(1, p, nS, nSmp) - p) / (sqrt(p*(1-p)));
        %------------------------------------------------

    case 'exp'          % no parameter
        % skew = 2; excess kurtosis = 6;
        para_arr_real = [];
        
        S_t1 = exprnd(1, nS, nSmp) - 1;
        %------------------------------------------------

    case 'laplace'      % no parameter
        % skew = 0; excess kurtosis = 3;
        b = 1/sqrt(2);
        signs = binornd(1, 0.5, nS, nSmp)*2 - 1;
        para_arr_real = [];
        
        S_t1 = signs.*exprnd(b, nS, nSmp);
        %------------------------------------------------

    otherwise
        error(['Not defined Distribution:', strDistr]);
        %--------------------------------------------------------------
end%END OF SWITCH

% assign the first dim
TData(:,1,:) = S_t1; %@[09/12/2016] row-wise, i.e. each row is treated as a random vector.

%% Pattern

% default input
if ~exist('pattern_str', 'var')
    pattern_str = lower('plus0_1');
end

switch lower(pattern_str);      
    case lower('power')
        Mi = rand(nDim-1, 1)*2; %@Power matrix %@U(-2,2)
        for q = 2 : nDim
            qMi = Mi(q-1);
            TData(:,q,:) = sign(S_t1) .* abs(S_t1).^repmat(qMi,nS,nSmp);
        end
        %-------------------------------------------------------
    case lower('linear')
        AB = rand(nDim-1,2)*2; %@U(-2,2)
        for q = 2 : nDim
            TData(:,q,:) = AB(q-1,1)*S_t1 + AB(q-1,2);
        end
        %-------------------------------------------------------
     case lower('plus0_1') %@Linear+0.1
        for q   = 2 : nDim
            TData(:,q,:) = S_t1 + (q-1)*0.1;
        end
        %-------------------------------------------------------        
    case lower('inUse') %[09/11/2016 for RAMICA_AAAI'17]
        a_arr = 1+randn(nS,1); %@ mean=1, std=1
        b_arr = rand(nS,1);
        for q = 2 : nDim
            TData(:,q,:) = repmat(a_arr,1,nSmp).*S_t1 + repmat((q-1)*b_arr,1,nSmp);
        end
        %-------------------------------------------------------
    otherwise
        error('pattern_str not defined.');
end

%%%%%%%%%%%%%%%%%%Data Vector & Normalize %%%%%%%%%%%%%%
for is = 1 : nS
    iS = squeeze(TData(is,:,:));
    iSvec = iS(:)';
    
    % Normalization
    iSvec_nm = iSvec ./ std(iSvec);
    iS_nm = reshape(iSvec_nm, size(iS)); % recover tensor S: check OK.
    % re-assignment
    iSvec = iSvec_nm; 
    iS = iS_nm; 
    
    % assignment
    TData(is,:,:) = iS;
    Data(is,:) = iSvec;
end

end%END OF FUNCTION

%% %%%%%%%%%%%%%%%%%%%%%%%%%Script%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

clear,clc, close all
Para_cel = {'pearsrnd', 'unif', 'tdistr', 'binomial', 'exp', 'laplace'};

strDistr = 'pearsrnd';
para_arr = [];

nS = 4;
nDim = 32;
nSmp = 64;
seed = 1;

% pattern_str = 'linear'; 
% pattern_str = 'power'; 
% pattern_str = 'plus0_1'; 
pattern_str = 'inUse'; 

% main
[TData, Data, para_arr_real] = genRandData(strDistr, para_arr, nS, nDim, nSmp, seed, pattern_str);

%%%%%%%%%%%%%%%%% Visulization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for p = 1:nS
    iS = squeeze(TData(p,:,:));
    subplot(2,2,p), imshow(iS,[])
end
title_all([strDistr, ', ', pattern_str])
%}
