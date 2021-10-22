function configurate()
% Configurate the folders and paths for this project. 
% Please run this scirpt at first.
% 
% Liyan Song: songly@sustech.edu.cn
% Oct.2021

clear, clc

% add code path
addpath(genpath(['.', filesep]));    

% screen
fprintf('Experimental environment has been set up. Have fun.\n') 
fprintf('Running in the file:\n\t%s.\n',pwd);

end
