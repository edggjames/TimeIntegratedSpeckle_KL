clc
close all hidden
clearvars

% linewidth and fontsize
LW = 1.5;
fs = 15;

%% add functions to file_path
addpath(pwd,"Functions\")

%% NB
% run script generate_standard_data.m prior to running this script

% load standard data
load Data\comp_data

%% generate vector of exposure times (Song et al. 2016)
log_vec = 10.^(-3:0.1:3);

T = t_c.*log_vec;
T(T==0) = [];
num_T_vals = length(T);
K_sq = zeros(1,num_T_vals);

%% generate images calculate speckle contrast values
for i = 1:num_T_vals
    % generate new seed images per T_val
    M_1 = generateM(dim_1,dim_2);
    M_2 = generateM(dim_1,dim_2);
    I_int = integrateLASCASpeckleOverTau(M_1,M_2,T(i),t_c);
    K_sq(i) = calcGlobalK(I_int)^2;
    disp(i)
end
%% 
save Data\Song_data