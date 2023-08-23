% Author: Edward James, PhD Student, UCL, November 2018.
% e.james.14@ucl.ac.uk

% script to load precomputed ACF data from AOM-DCS and fit to a
% modulated DCS model of solution to CDE in semi-inifinite geometry using
% autocorrelationFitAOTCW.m

clearvars
clc
close all

%% add functions to file_path
addpath(pwd,"Functions\")

%% To set line widths and font size for plots
fs=20;
LW=1.5;

%% assign lower and upper bounds for tau (i.e autocorrelation time delay)
lower_bound = 1e-7;
upper_bound = 1e-2;
% refractive index of biological tissue (from Dong et al. 2013)
n = 1.34; 
% wavelength of laser in cm in vaccum
lambda = 784e-7; 
% in-plane 1D distance between source and detector in cm
rho = 2.05; 
% set optical properties of sample
mu_a = 0.10; % cm^-1
mu_s_p = 11.275; %cm^-1
% specify ranges for 3 parameters to recover (as per Dong et al. 2013)
beta_range = 0.1:0.01:1.0; % unitless coherence factor
alpha_Db_range = 0.4e-8:0.1e-8:2e-8; % cm^2/s
% acoustic frequency
f_a = 2e6;

% load experimental data
data = importdata(strcat("Data\AOM_experimental_data.SIN"));

% Run semi-infinite geometry DCS fit function with US modulation
[tau,~,beta_brownian,alpha_Db_brownian,A] = ...
    autocorrelationFitAOTCWf0(data,lower_bound,upper_bound,n,lambda,rho,mu_a,...
    mu_s_p,beta_range,alpha_Db_range,f_a,1);
title({'SDS: 20.5 mm - Duration: 10 s - Source Power: 40 mW - Intralipid Dilution: 1:19';['$\mu_a:'...
    ,num2str(mu_a),'cm^{-1} - \mu_s'':',num2str(mu_s_p),'cm^{-1}$ - PZT - 2 MHz - 3.92 $V_{RMS}$']},...
    'Interpreter','Latex')

%% save data
save Data\AOM_data tau beta_brownian alpha_Db_brownian A
