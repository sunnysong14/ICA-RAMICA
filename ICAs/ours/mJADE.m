function [Ainv_est, TS_est, Uinv_est, White] = mJADE(TX)
% USAGE
%   JADE for matrix ICA. I.E. the sources are 2D images. We assume that
%   #mixtures = #sources.
% 
% INPUT
%   TX: The mixture data tensor X. Last dimension contains the samplings,
%   the first dimension corresponds to one component, and the second
%   dimension contains all elements of this component.
% 
% OUTPUT
%   Ainv_est: Demixing matrix
%   Uinv_est: Whitented Demixing matrix: Uinv_est * TX_white = Sources
%   TS_est: The estimated tensor source
%   White: The whitening matrix
% 
% Liyan for AAAI on 07/18/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = ndims(TX) - 1; %order of samples
if N~=2
    error('Error data order. Now we can only deal with tensor and N=2.'); 
end
[nX,dimXi,nSpl] = size(TX); % #source=#mixutre


%%%% 1. Centering %%%%
TXmean = mean(TX, N+1);
TX = TX - repmat(TXmean,[ones(1,N),nSpl]);


%%%% 2. Whitening %%%%
[E,D] = pcamatLY(TX);

%%%% Debug: PCA may reduces dim(E) or dim(D) causing errors in mWhitenv().
if min(size(E)) < max(size(E))
    Ainv_est = -1; %assign -1 so that Amari will not be computed
    TS_est = -1;
    Uinv_est = -1;
    White = -1;
    return
end

[TX_white, White] = mWhitenv(TX, E, D); %core


%%%% 3. MJADE Recover %%%%

%%%% 1) Sample cumulant tensor %%%%
R = eye(nX);
CMT = -ones(nX, nX, nX, nX); %init

for i1 = 1 : nX
    i1X = squeeze(TX_white(i1,:,:)); %should be 256x250 with random vectors in columns.
    
    for i2 = 1 : nX
        i2X = squeeze(TX_white(i2,:,:));
        
        for i3 = 1 : nX
            i3X = squeeze(TX_white(i3,:,:));
            
            for i4 = 1 : nX
                i4X = squeeze(TX_white(i4,:,:));
                
                CMT(i1,i2,i3,i4) = mean(mean(i1X.*i2X.*i3X.*i4X)) ...
                    - R(i1,i2)*R(i3,i4) ...
                    - R(i1,i3)*R(i2,i4) ...
                    - R(i1,i4)*R(i2,i3) ;
            end%i4
        end%i3
    end%i2
end%i1

%%%% 2) Cumulant basis matrices %%%%

% % % case-1: identify cumulant basis matrices
% % Compute basis matrices {E^{ij}}
% % nCMB = nX^2; % #(cumulant basis: E^{ij})
% % CBM = zeros(nX, nX, nCMB); % cumulant matrices
% % % Form E^{ij}
% % nn = 1;
% % while nn <= nCMB
% %     for i1 = 1:nX
% %         for i2 = 1:nX
% %             CBM(i1,i2,nn) = 1;
% %             nn = nn+1;
% %         end
% %     end
% % end


% case-2: JADE's cumulant basis matrices <-- Use This One Now.
% compute a basis of the space of symmetric nX*nX matrices
nCMB = (nX*(nX+1))/2; % Dim. of the space of real symm matrices
CBM = zeros(nX,nX,nCMB); % Holds the basis.   
icm = 0; % index to the elements of the basis
Id = eye(nX); % convenience
for im = 1:nX
    vi = Id(:,im); %the ith basis vetor of R^nX
    icm = icm + 1;
    CBM(:,:,icm) = vi*vi';
    
    for jm = 1:im-1
        vj = Id(:,jm); % the jth basis vetor of R^nX
        icm = icm + 1;
        CBM(:,:,icm) = sqrt(0.5) * (vi*vj'+vj*vi');
    end%jm
end%im
%>>CBM(:,:,i) is the i-th element of an orthonormal basis for the space of
%nX*nX symmetric matrices.

% Note: the two cases have the same demixing matrix in a sense that
% multipling a permutation matrix. 


%%%% 3) Images of basis cumulant matrices with cumulant tensor %%%%
CMM = zeros(nX,nX,nCMB); %init
for nn = 1:nCMB
    Eij = CBM(:,:, nn);
    
    for i1 = 1 : nX
        for i2 = 1 : nX
            CMM(i1, i2, nn) = sum(sum(squeeze(CMT(i1,i2,:,:)).*Eij));
        end%i2
    end%i1
end%nn


%%%% 4) If use the above, we need to update CMM-->CM %%%%
% Rearrange tensor of cumulant matrices to a wide matrix
CM = zeros(nX, nX*nCMB);
for nn = 1 : nCMB
    CM(:, (nn-1)*nX+1:((nn-1)*nX+nX)) = squeeze(CMM(:,:,nn));
end%nn


%%%% 5) Jacobi Joint Diagonalization of CM %%%%
WX = eye(nX); % Whitening has been done outside.
verbose	= 0;
Uinv_est = Jacobi(CM, WX, nX, nX, nSpl, nCMB, verbose);


%%%% 4. Estimated Source %%%%
TS_est = zeros(size(TX));
for tt = 1 : size(TS_est,3)
    TS_est(:,:,tt) = Uinv_est * TX_white(:,:,tt);
end%tt


%%%% Demix Matrix %%%%
Ainv_est = Uinv_est * White;

end%function
