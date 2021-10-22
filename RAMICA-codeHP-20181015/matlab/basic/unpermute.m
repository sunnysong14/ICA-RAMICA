function W_new = unpermute(A, W)
% function W_new = unpermute(A, W)
% inputs:  A -- actual mixing matrix
%          W -- demixing matrix
% output:  a permutation of the rows of W such that W(i, :) * A is
% approximately the ith indicator function
    d = size(A, 1);
    W_new = zeros(d, d);
    indices = 1:d;
    newspots = zeros(d, 1);
    % eps = 1e-4;  % eps is actually a matlab default defined to be some
    % small number.
    
    PMat = abs(W * A) + 10*eps* ones(d, d); % adding ones to make sure no terms are zero in our "permutation matrix"
    
    for i = 1:d
        max_val = max(max(PMat));
        [r, c] = find(PMat >= max_val-eps, 1);
        W_new(c, :) = W(r, :);
        
        % Make certain we do not reuse any row or column
        % PMat
        % r
        % d
        PMat(r, :) = zeros(1, d);
        PMat(:, c) = zeros(d, 1);
    end
end