%================================================================
% Codes to add all files under this file to MATLAB path:
%   Liyan Dec.4th 2015 Friday
%================================================================

% cd to this directory

% add recursive path
addpath(pwd) 			%<-- Add the tensor toolbox to the MATLAB path

% Also add the met directory
cd met; 
addpath(pwd)                    %<-- [OPTIONAL]

% Save for future MATLAB sessions
savepath

% remove path 
% rmpath(genpath(pwd))          %<-- remove path
%================================================================