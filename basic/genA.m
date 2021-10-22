function A = genA(nS, seed)
% Usage: 
%   Generate mixing matrix $A$
% 
% nS: #source
% seed: random seed
% A: the generated mixing data
% 
% Liyan Dec.25 2015 for IJCAI16
%   last updated 04-14-2016

%%
% contrl random
rng(seed);

% generate mixing matrix A: random with ||column||=1 & plus identity matrix
A = rand(nS, nS);
for c_ = 1 : nS
    A(:, c_) = A(:, c_) / norm(A(:,c_));
end
A = A + eye(nS);

end  % End OF FUNCTION
%% test
%{

nS = 4; 
seed = 1;

A = genA(nS, seed);
%}
