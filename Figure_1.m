clc
clearvars
close all

%% Figure 1 from Edward James, Samuel Powell, and Peter Munro, "Simulation
% of statistically accurate time-integrated dynamic speckle patterns in
% biomedical optics," Opt. Lett. 46, 4390-4393 (2021)

% script to develop taking the weighted sum of N statistically indepdenent
% speckle patterns, where the weighting factors are the eigenvalues of the
% g_1 function

%% add functions to file_path
addpath(pwd,"Ancillary_scripts\Functions\")

%% NB
% run script generate_standard_data.m prior to running this script

%% load data
load Ancillary_scripts\Data\comp_data
clear delta N_s radius T M_1 M_2
% linewidth
LW = 1.5;

%% specify parameters
N = 500; % number of eigen values
T = 1*t_c;
alpha = 0.9;

%% calculate eigen values
t2 = linspace(0,T,N);
t1 = t2;
[n,k] = meshgrid(t2,t1);
diff = abs(k-n);
K = alpha*exp(-diff/t_c) + (1-alpha);
lambda = eigs(K,N)/N;

%% show plots
figure('units','normalized','outerposition',[0 0 1 1])

fs = 24;
subplot(1,2,1)
imagesc(t1*1000,t2*1000,K)
ax = gca;
ax.FontSize = fs;
colormap bone
xlabel('$t_1$ (ms)','FontWeight','bold','Interpreter','Latex','Fontsize',fs+3)
ylabel('$t_2$ (ms)','FontWeight','bold','Interpreter','Latex','Fontsize',fs+3)
title('$(a)$','FontWeight','bold','Interpreter','Latex','Fontsize',fs+5)
axis square
colorbar
xticks([0 0.10 0.20 0.30])
yticks([0 0.10 0.20 0.30])

subplot(1,2,2)
scatter(1:N,lambda,'k.')
ax = gca;
ax.FontSize = fs;
set(gca,'YScale','log');
xlabel('$n$','FontWeight','bold','Interpreter','Latex','Fontsize',fs+3)
ylabel('$\lambda_n$','FontWeight','bold','Interpreter','Latex','Fontsize',fs+3)
title('$(b)$','FontWeight','bold','Interpreter','Latex','Fontsize',fs+5)
axis square
xlim([0 N])
box on