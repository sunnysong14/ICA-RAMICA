function [data_tensor, data_mat] = generate_synthetic_tensors(syn_data_name, nb_source, nb_dimension, nb_sample, seed)
% Generate synthetic tensor of size nb_dimension * nb_sample * nb_sample.
%   * It generates the "row-wise" synthetic data for which data structure
%   is shown in row.
%   * This synthetic data generation can force the correlation among
%   columns.
% 
% Returns
%   data_tensor: tensor data of the size nb_source * nb_dim * nb_sample
%   data: vector data of size nb_source * (nb_dimension * nb_sample)
% 
% Liyan Song: songly@sustech.edu.cn
% 2021-10-20

rng(seed)  % contral random
data_tensor = zeros(nb_source, nb_dimension, nb_sample);  % init
data_mat = zeros(nb_source, nb_dimension*nb_sample);  % init

switch lower(syn_data_name)
    case lower('pearsrnd')
        mean_value = 0; 
        std_value = 1;
        skew_value = 1; 
        kurt_value = 1+3;  % kurt(Gaussian) = 3        
        S_t1 = pearsrnd(mean_value, std_value, skew_value, kurt_value, nb_source, nb_sample);

    case lower('tdistr')
        % Variance = paraDf / (paraDf-2) (must have paraDf >= 3 in order
        % for variance assumptions on S to be met).
        % skew = 0 (undefined if paraDf <= 3). kurtosis = 6/(paraDf-4)
        df = 5;
        S_t1 = trnd(df, nb_source, nb_sample) / sqrt(df / (df - 2));

    case 'exp'
        % skew = 2; excess kurtosis = 6;        
        S_t1 = exprnd(1, nb_source, nb_sample) - 1;

    case 'laplace'
        b = 1/sqrt(2);
        signs = binornd(1, 0.5, nb_source, nb_sample)*2 - 1;        
        S_t1 = signs.*exprnd(b, nb_source, nb_sample);
end%switch

% design the 1st dimension
data_tensor(:, 1, :) = S_t1;  % row-wise: each row is treated as a random vector

% Pattern
a_arr = 1 + randn(nb_source, 1);  % mean=1, std=1
b_arr = rand(nb_source,1);
for q = 2 : nb_dimension
    data_tensor(:,q,:) = repmat(a_arr,1,nb_sample).*S_t1 + repmat((q-1)*b_arr,1,nb_sample);
end

% data vector & normalize
for ss = 1 : nb_source
    iS = squeeze(data_tensor(ss,:,:));
    iSvec = iS(:)';
    
    % normalization
    iSvec_nm = iSvec ./ std(iSvec);
    iS_nm = reshape(iSvec_nm, size(iS)); % recover tensor S: check OK.
    iSvec = iSvec_nm; 
    iS = iS_nm; 
    
    % assignment
    data_tensor(ss,:,:) = iS;
    data_mat(ss,:) = iSvec;
end

% visulization
is_show = false;
if is_show
    for p = 1:nb_source
        iS = squeeze(data_tensor(p,:,:));
        subplot(2,2,p), imshow(iS,[])
    end
end
end%END OF FUNCTION


