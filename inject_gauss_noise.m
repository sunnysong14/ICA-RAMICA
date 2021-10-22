function [X_noisy, noise_cov] = inject_gauss_noise(X, p_noise, seed)
% Inject additive Gaussian noise into the data matrix X.
% The noise covariance is set to p_noise*I for which p_noise = sigma^2.
% 
% Parameters:
%   X: mixture matrix. Each row corresponds one source, and each column is
%      a sampling across all sources.
%   p_noise: a positive constant defining the noise powers with its value
%      between 0 and 1. In math. Note that p_noise=\sigma^2.
%   X_pn: the corrupted X with p_noise-level gaussian noise.
% 
% Liyan Song: songly@sustech.edu.cn
% Oct. 2021

rng(seed);  % control random
[nb_source, nb_sample] = size(X);
noise_means = zeros(nb_source, nb_sample);

% coverience of the noise
noise_cov = p_noise*eye(nb_source);

% noisy observation
X_noisy = X + mvnrnd(noise_means', noise_cov)';

end
