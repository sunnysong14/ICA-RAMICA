% USAGE 
%   Script for the second experiment, i.e. bline image separation
% 
% Liyan for AAAI17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc

CaseScen = 1; %only this scenario
P_noise = [0:0.01:0.1, 0.15, 0.2];
StrICA_lst = {'infomax','FastICA','jade','ours'};
CRImg_lst = {'c','r'};
nS = 4;
Seed = 1:100;

%%%% Main Run %%%%
[Ave_amari, Std_amari, Amari, ...
    S, A, X_pn, X, RTS_est, Ainv_est] = ...
    mainIS_fun(CaseScen,nS,P_noise,StrICA_lst,Seed,CRImg_lst);
