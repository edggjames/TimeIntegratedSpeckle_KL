clc
close all
clearvars

%% add functions to file_path
addpath(pwd,"Functions\")

%% generate standard data set for testing
% specify image dimensions
dim_1 = 600; % y in paper
dim_2 = 600; % x in paper
% then generate two M matrices
M_1 = generateM(dim_1,dim_2);
M_2 = generateM(dim_1,dim_2);
% generate matrix H, the pupil plane in the frequency domain
radius = 100;
H = generateH(dim_1,dim_2,radius);
t_c = 370e-6;
T = 5*t_c;
m = 1e6;
N_s = 3; % must be odd, usually 3,5,7,9, 11 - dependent on speckle size
delta = eps;

save Data\comp_data