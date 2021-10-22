function A = mixing_matrix(nb_source, seed)
% Setup the mixing matrix A.
% 
% Liyan Song: songly@sustech.edu.cn

rng(seed);  % contrl random

% random with ||column||=1 & plus identity matrix
A = rand(nb_source, nb_source);
for cc = 1 : nb_source
    A(:, cc) = A(:, cc) / norm(A(:,cc));
end
A = A + eye(nb_source);
end
