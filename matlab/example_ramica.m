% USAGE 
%   Script for the first experiment, i.e. BSS on the synthetic data.
% 
% Liyan Song: songly@sustech.edu.cn
% October, 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% synthetic experimental setup
ica_name_cell = {'infomax', 'FastICA', 'jade', 'ramica'};
syn_data_name = 'tdistr'; 
nb_source = 4; 
nb_dimension = 32; 
nb_sample = 64;
seed_array = 1:10;
noise_level = 0;  % [0:0.01:0.1,0.15,0.2]

% main run
amari_errors_runs = nan*ones(length(ica_name_cell), length(seed_array));  % init
for ic = 1 : length(ica_name_cell)
    ica_name = lower(ica_name_cell{ic});
    for ss = 1 : length(seed_array)
        seed = seed_array(ss);
        
        % generate source and mixture tensors
        [source_tensor, source_mat] = generate_synthetic_tensors(syn_data_name, nb_source, nb_dimension, nb_sample, seed);

        % mixing marix
        A = mixing_matrix(nb_source, seed);
        mixture_tensor = ttm(tensor(source_tensor), A, 1); 
        mixture_tensor = mixture_tensor.data;

        % mixture data tensor
        mixture_mat = zeros(nb_source, nb_dimension * nb_sample);
        for ii = 1 : nb_source
            i_mixture_tensor = mixture_tensor(ii,:,:);
            mixture_mat(ii,:) = i_mixture_tensor(:)';
        end

        % noisy mixture data
        if noise_level > 0 
            [mixture_mat_noisy, noise_cov] = inject_gauss_noise(mixture_mat, noise_level, seed);
            
            % stack back to mixture_tensor_noisy
            mixture_tensor_noisy = zeros(nb_source, nb_dimension, nb_sample);
            for t = 1 : nb_sample
                mixture_tensor_noisy(:,:,t) = mixture_mat_noisy(:, ((t-1)*nb_dimension+1):(t*nb_dimension));
            end
        elseif noise_level == 0 
            mixture_mat_noisy = mixture_mat;
            mixture_tensor_noisy = mixture_tensor;
            noise_cov = 0;
        else
            error('p_noise=\sigma^2 should be no smaller than 0.');
        end

        % run ICA method
        switch lower(ica_name)
            case 'fastica'
                [sources_recover, ~, A_recover] = fastica(mixture_mat);
            case 'jade'
                A_recover = jadeR(mixture_mat);
                sources_recover = A_recover * mixture_mat;
            case 'infomax'
                A_recover = infomaxICA(mixture_mat);
                sources_recover = A_recover * mixture_mat;
            case lower('ramica')
                [A_recover, source_tensor_recover] = ramica(mixture_tensor_noisy);
                sources_recover = source_tensor_recover;
            otherwise
                error('No such ICA method.');
        end%switch

        % Amari error
        if size(A_recover,1)<size(A,1) || sum(sum(isnan(A_recover)))>0 
            % if cannot fully recover all the resouces
            amari = nan;                    
        else  % if fully recovery
            amari = AmariIndex(A, A_recover);
        end
        amari_errors_runs(ic, ss) = amari;
    end%seed
end%ica

% ave amari error across seeds. 
amari_errors_ave = nanmean(amari_errors_runs, 2);
fprintf("p_noise is %0.3f, the amari errors of ICA methods are as below\n", noise_level);
amari_errors_ave
fprintf("Suceed!\n")







