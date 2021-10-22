function [newVectors, whiteningMatrix, dewhiteningMatrix] = mWhitenv ...
    (TX, E, D, s_verbose)
% Log:
%   LY extends this funciton for tensor TX <-- vectors, and is used to
%   separate images. 
%   TX: tensor with the size(TX)-->[dimXi, nX,nSample], and after
%   centering.
% 
% For AAAI from 17/07/2016
% 
%WHITENV - Whitenv TX.
%
% [newVectors, whiteningMatrix, dewhiteningMatrix] = ...
%                               whitenv(TX, E, D, verbose);
%
% Whitens the data (row TX) and reduces dimension. Returns
% the whitened TX (row TX), whitening and dewhitening matrices.
%
% ARGUMENTS
%
% TX       Data in row TX.
% E             Eigenvector matrix from function 'pcamat'
% D             Diagonal eigenvalue matrix from function 'pcamat'
% verbose       Optional. Default is 'on'
%
% EXAMPLE
%       [E, D] = pcamat(TX);
%       [nv, wm, dwm] = whitenv(TX, E, D);
%
%
% This function is needed by FASTICA and FASTICAG
%
%   See also PCAMAT

% @(#)$Id: whitenv.m,v 1.3 2003/10/12 09:04:43 jarmo Exp $

% ========================================================
% Default value for 'verbose'
if nargin < 4, s_verbose = 'off'; end

% Check the optional parameter verbose;
switch lower(s_verbose)
 case 'on'
  b_verbose = 1;
 case 'off'
  b_verbose = 0;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end

% ========================================================
% In some cases, rounding errors in Matlab cause negative
% eigenvalues (elements in the diagonal of D). Since it
% is difficult to know when this happens, it is difficult
% to correct it automatically. Therefore an error is 
% signalled and the correction is left to the user.
if any (diag (D) < 0),
    error (sprintf (['[ %d ] negative eigenvalues computed from the' ...
		   ' covariance matrix.\nThese are due to rounding' ...
		   ' errors in Matlab (the correct eigenvalues are\n' ...
		   'probably very small).\nTo correct the situation,' ...
		   ' please reduce the number of dimensions in the' ...
		   ' data\nby using the ''lastEig'' argument in' ...
		   ' function FASTICA, or ''Reduce dim.'' button\nin' ...
		   ' the graphical user interface.'], ...
           sum (diag (D) < 0)));
end

% ========================================================
% Calculate the whitening and dewhitening matrices (these handle
% dimensionality simultaneously).
whiteningMatrix = sqrt(D)\E'; %inv (sqrt (D)) * E';%<-- the original code.
dewhiteningMatrix = E * sqrt (D);

% Project to the eigenvectors of the covariance matrix.
% Whiten the samples and reduce dimension simultaneously.
if b_verbose, fprintf ('Whitening...\n'); end

%%%%%%%%%%%%%%%%%%%%%% LY: tensor whiening %%%%%%%%%%%%%%%%%%%%%%%%%%%
N = ndims(TX) - 1; % The order of samples.
numSpl = size(TX,3); % dimXi: dim of random vector; nMix: #(mixture)

if N == 1 % if "TX" is vectors, use pcamat()'s original codes.
    % The correctness has been justified by obtaining the same pca results
    % of this function and pcamat() for matrix X. 07/17/2016.
    newVectors =  whiteningMatrix * TX;
    
elseif N == 2  % if "TX" is real matrices, use my definition.
    TX_white = zeros(size(TX));
    for j = 1 : numSpl
        iTX = TX(:,:,j);
        TX_white(:,:,j) = whiteningMatrix * iTX;
    end
    newVectors = TX_white;
end
%%%%%%%%%%%%%%%%%%%%%% LY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ========================================================
% Just some security...
if ~isreal(newVectors)
  error ('Whitened TX have imaginary values.');
end

% Print some information to user
if b_verbose
  fprintf ('Check: covariance differs from identity by [ %g ].\n', ...
    max (max (abs (cov (newVectors', 1) - eye (size (newVectors, 1))))));
end
