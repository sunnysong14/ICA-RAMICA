% Usage:
%   Inject additive Gaussian noise into X.
% 	The covariance of noise is assumed to be $p_noise*I$ where $p_noise=sigma^2$.
% 
%   X: observation matrix. The row corresponds one source, and each column is
% one sampling across all sources.
%   p_noise: a positive constant defining the noise powers with its value between 0
% and 1. In math, $p_noise=\sigma^2$
%   X_pn: the corrupted X with p_noise-level gaussian noise.
% 
% Liyan for ECML16 03-28-2016
%   last updated on 04-14-2014
%%
function [X_pn, Noise_cov] = injctGauNoise(X, p_noise, seed)

% control random
rng(seed);

% basic info
[nSrc, nSmp] = size(X);
Noise_mean = zeros(nSrc, nSmp);

% Cov(noise)
Noise_cov = p_noise*eye(nSrc);

% noisy observation
X_pn = X + mvnrnd(Noise_mean', Noise_cov)';

end  % END OF FUNCTION

%% test
%{
clear,clc
close all

% 
[X, S, A, nameSource] = genDataSyn(2);
% 
p_noise = 0.6; 

% main call
X_pn = injGausNoiseObs(X, p_noise);
%}
