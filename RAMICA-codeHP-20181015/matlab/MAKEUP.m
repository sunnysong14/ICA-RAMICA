function MAKEUP()
% USAGE
%   It settles downs all the PATH required.
% 
% Liyan 2015/05/18/ Mon.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc

% add code path
addpath(genpath(['.', filesep]));    

% screen
fprintf('Exp environment set up accomplished. Have fun.\n') 
fprintf('Running in the file:\n\t%s.\n',pwd);
end
