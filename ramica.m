function [A_recover, source_tensor_recover, whitened_recover_matrix, whitening_matrix] = ramica(mixture_tensor)
% The proposed ramica method.
% 
% :para mixture_tensor - the mixture data tensor X. Last dimension contains
% the samplings, the 1st dimension corresponds to one component, and the
% 2nd dimension contains all elements of this component.
% 
% Liyan Song: songly@sustech.edu.cn
% Oct. 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = ndims(mixture_tensor) - 1;  % order of samples
if N~=2
    error('Error data order. Now we can only deal with tensor and N=2.'); 
end
[nb_X, nb_dimension, nb_sample] = size(mixture_tensor); % #source=#mixutre

% 1. Centering
mixture_mean_tensor = mean(mixture_tensor, N+1);
mixture_tensor = mixture_tensor - repmat(mixture_mean_tensor,[ones(1,N),nb_sample]);

% 2. Whitening
[E,D] = pcamat_ramica(mixture_tensor);
[whitened_mixture_tensor, whitening_matrix] = whiten_ramica(mixture_tensor, E, D);  % core

% 3. RAMICA Recover ============================
% 3-1. compute sample cumulant tensor
R = eye(nb_X);
CMT = -ones(nb_X, nb_X, nb_X, nb_X);  % init
for i1 = 1 : nb_X
    i1X = squeeze(whitened_mixture_tensor(i1,:,:));  % This should be 256x250 with random vectors in columns.    
    for i2 = 1 : nb_X
        i2X = squeeze(whitened_mixture_tensor(i2,:,:));
        for i3 = 1 : nb_X
            i3X = squeeze(whitened_mixture_tensor(i3,:,:));
            for i4 = 1 : nb_X
                i4X = squeeze(whitened_mixture_tensor(i4,:,:));                
                CMT(i1,i2,i3,i4) = mean(mean(i1X.*i2X.*i3X.*i4X)) ...
                    - R(i1,i2)*R(i3,i4) ...
                    - R(i1,i3)*R(i2,i4) ...
                    - R(i1,i4)*R(i2,i3) ;
            end%i4
        end%i3
    end%i2
end%i1

% 3-2. Cumulant basis matrices
% compute a basis of the space of symmetric nX*nX matrices
nCMB = (nb_X*(nb_X+1))/2;  % Dim. of the space of real symm matrices
CBM = zeros(nb_X,nb_X,nCMB);  % Holds the basis.   
icm = 0;  % index to the elements of the basis
Id = eye(nb_X);  % convenience
for im = 1:nb_X
    vi = Id(:,im);  % the ith basis vetor of R^nX
    icm = icm + 1;
    CBM(:,:,icm) = vi*vi';
    for jm = 1:im-1
        vj = Id(:,jm); % the jth basis vetor of R^nX
        icm = icm + 1;
        CBM(:,:,icm) = sqrt(0.5) * (vi*vj'+vj*vi');
    end%jm
end%im
% >>CBM(:,:,i) is the i-th element of an orthonormal basis for the space of
% nX*nX symmetric matrices. 

% 3-3. Images of basis cumulant matrices with cumulant tensor
CMM = zeros(nb_X,nb_X,nCMB); %init
for nn = 1:nCMB
    Eij = CBM(:,:, nn);
    
    for i1 = 1 : nb_X
        for i2 = 1 : nb_X
            CMM(i1, i2, nn) = sum(sum(squeeze(CMT(i1,i2,:,:)).*Eij));
        end%i2
    end%i1
end%nn

% 3-4. update CMM-->CM, 
% rearrange tensor of cumulant matrices to a wide matrix
CM = zeros(nb_X, nb_X*nCMB);
for nn = 1 : nCMB
    CM(:, (nn-1)*nb_X+1:((nn-1)*nb_X+nb_X)) = squeeze(CMM(:,:,nn));
end%nn

% 3-5. Jacobi Joint Diagonalization of CM
WX = eye(nb_X); % whitening has been done outside.
verbose	= 0;
whitened_recover_matrix = Jacobi(CM, WX, nb_X, nb_X, nb_sample, nCMB, verbose);

% 4. Estimated Source
source_tensor_recover = zeros(size(mixture_tensor));
for tt = 1 : size(source_tensor_recover,3)
    source_tensor_recover(:,:,tt) = whitened_recover_matrix * whitened_mixture_tensor(:,:,tt);
end

% the recover (demix) matrix
A_recover = whitened_recover_matrix * whitening_matrix;

end
