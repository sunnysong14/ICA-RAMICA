function ret = AmariIndex(A, W)
%  Inputs:
%   A -- The true mixing matrix used to mix the signals in ICA.
%   W -- The recovered demixing matrix (approximatesd A^(-1)).
%  Returns:
%   ret -- the Amari Index, which measures how successful we were at
%   recovering the original matrix.
%
%   We use the convention that 0 / 0 = 1 when computing the Amari index,
%   since this is the cap achievable within the Amari index calculations so
%   long as no zero division occurs.  This issue can arise when W does not
%   recover all rows of inv(A).
%
%   As the Amari index assumes that W' and A are of the same size, W is
%   filled in with zero rows if not all rows of inv(A) are recovered to make
%   this true.  If W has insufficient columns, then an error will be thrown.
%
%  The Amari index is defined as the ICA experimental error metric in the paper:
%   A New Learning Algorithm for Blind Signal Separation (1996)
%   by: Amari, Cichocki, and Yang
    
    N = size(A, 2); % number of latent signals
    if size(W, 1) < size(A, 2)
        W = [W; zeros(N-size(W, 1), size(W, 2))];
    end
    
    M = W*A;
    
    %% Test to see if we have a valid matrix for decomposition
    % Soft fail for an invalid matrix will be taken as N-1 (cap of Amari index)
    try
        svals = svd(M);
    catch
        fprintf('Error:  Cannot take svd.\n')
        ret = d;
        return;
    end
    if ~isreal(svals)
        fprintf('Error:  svd Output contains imaginary components.\n');
        ret = inf;
        return;
    end
    M = abs(M);

    %% Actually compute the Amari index.
    M_row_maxs = max(M, [], 2);
    M_col_maxs = max(M, [], 1);
    
    ret = - 2*N;
    for i = 1:N
        if M_row_maxs(i) <= 0
            ret = ret + N;
        else
            ret = ret + sum(M(i, :)) / M_row_maxs(i);
        end
            
        if M_col_maxs(i) <= 0
            ret = ret + N;
        else
            ret = ret + sum(M(:, i)) / M_col_maxs(i);
        end
    end
end

