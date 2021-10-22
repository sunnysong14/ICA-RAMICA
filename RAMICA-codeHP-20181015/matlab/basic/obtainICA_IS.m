function [X_pn, X, S, A, Noise_cov] = obtainICA_IS(S_org, p_noise, seed)
% Usage:
%   Obtain ICA data for image seperation, including sources, mixing matrix,
%   noise-free X and noisy X.
% 
% INPUT PARAMETERS:
%   S_org -- The original sources without any preprocessing.
%   p_noise -- a positive constant defining the noise powers with its value between 0
% and 1. In math, $p_noise=\sigma^2$.
% 
% OUTPUT PARAMETERS:
%   X_pn -- corrupted X with Gaussian noise level $p$.
%   X -- noise-freee observation.
%   S -- sources after preprocessing
%   A -- original mixing matrix
%   Noise_cov -- covariance matrix of noise. Used, e.g. for opt-SINR.
% 
% Liyan updated for acml16 08/09/2016, 08-10-2016: remove mat2gray() for
% non-image data. 08-11-2016 Add normalization of S back (large change on pf.)
%   Liyan for AAAI'16 on 08/22/2016, 20160825
%@Note: also used in synthetic experiments. maybe change another name.

%%
S = S_org;

% (1)convert to [0,1] for image --> not use [08/22/2016]
% S = mat2gray(S_org);
% figure, for k = 1:4, subplot(2,2,k), imshow(reshape(S(k, :),256,256),[]); end

% (2) Normalize S
% [08/11/2016 for PIMD] By doing normalization, we get better PF for
% all ICA methods, especially GIICA, PEGI, and PIMD. Also, PIMD can
% outperform others including FastICA, JADE, and others.
for s_ = 1 : size(S,1) 
    Si = S(s_, :);
    S(s_, :) = Si./std(Si);
end
% figure, for k = 1:4, subplot(2,2,k), imshow(reshape(S(k, :),256,256),[]); end

% mixing matrix
nS = size(S_org, 1);
A = genA(nS, seed);

% noise-free X
X = A * S;
% figure, for k = 1:4, subplot(2,2,k), imshow(reshape(X(k, :),256,256),[]); end

% inject noise
if p_noise > 0 
     [X_pn, Noise_cov] = injctGauNoise(X, p_noise, seed);
elseif p_noise == 0 
    X_pn = X;
    Noise_cov = 0;
else
    error('$p_noise=\sigma^2$ should be no smaller than 0. ');
end
% figure, for k = 1:4, subplot(2,2,k), imshow(reshape(X_pn(k, :),256,256),[]); end

end % END OF FUNCTION

%% Test
%{
clear,clc, close all

caseScen = 5;
p_noise = 0.1; 
seed = 1; 

[X_pn, X, S, A, Noise_cov] = obtainICA_IS(caseScen,p_noise, seed);
%}
