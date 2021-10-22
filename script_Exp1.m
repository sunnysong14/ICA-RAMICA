% USAGE 
%   Script for the first experiment, i.e. BSS on the synthetic data.
% 
% Liyan for AAAI17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc

%%%% Exp Design %%%%
Distri_cel = {'pearsrnd', 'tdistr', 'exp', 'laplace'};
para_arr = []; %use default paraS
StrICA_lst = {'infomax', 'FastICA', 'jade', 'ours'};
NS = [2 4 8 16];
NDim = [16,32,64,128];
NSmp = [16,32,64,128];
P_noise = [0:0.01:0.1,0.15,0.2];
Seed = 1:100;
% synthetic pattern
pattern_str = 'inUse'; %used in genRandData()

% rRAMICA or cRAMICA
Tag_revRAMICA = [0,1];


%%%% debug
if length(NDim) ~= length(NSmp)
    error('#NDim should equals to #NSmp.')
end


%%%% Main Run %%%%
for ndis = 1:length(Distri_cel)
    strDistr = lower(Distri_cel{ndis});
    
    for ns = 1 : length(NS)
        nS = NS(ns);
        
        for nd = 1 : length(NDim)
            nDim = NDim(nd);
            nSmp = NSmp(nd);            
            fprintf('Running %s: P=%d, (Q,T)=(%d,%d)\n', strDistr,nS,nDim, nSmp);
            
            for tg = 1:length(Tag_revRAMICA)
                tag_revRAMICA = Tag_revRAMICA(tg);

                %%%% main call
                [para_arr_real, Ave_amari, Std_amari, Amari, ...
                    S, A, X_pn, X, RTS_est, Ainv_est] = ...
                    mainSyn_fun(strDistr,para_arr,pattern_str,...
                    nS,nDim,nSmp,P_noise,StrICA_lst,Seed,tag_revRAMICA);
                
            end%tg
        end%nd
    end%ns
end%ndis
