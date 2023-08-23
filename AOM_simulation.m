% Author: Edward James, PhD Student, UCL, December 2021.
% e.james.14@ucl.ac.uk

% script to simulate 2D-TIDSP of DCS image with and without AOM

clearvars
clc
close all

%% NB
% run script Fit_D2_DCS_AOT_AOT_CW_precomp.m prior to running this script

%% add functions to file_path
addpath(pwd,"Ancillary_scripts\Functions\")

%% To set line widths and font size for plots
fs=20;
LW=1.5;

%% assign variables
% refractive index of biological tissue (from Dong et al. 2013)
n = 1.34; 
% wavelength of laser in cm in vaccum
lambda = 784e-7; 
% in-plane 1D distance between source and detector in cm
rho = 2.05; 
% set optical properties of sample
mu_a = 0.10; % cm^-1
mu_s_p = 11.275; %cm^-1
% acoustic frequency
f_a = 2e6;

k_0 = 2*pi*n/lambda; % wavevector magnitude of the incident light field
R_eff = -1.440*n^-2 + 0.710*n^-1 + 0.668 + 0.0636*n; % effective reflection
% coefficient to account for the index mismatch between air and tissue
z_0 = 1/mu_s_p; % the collimated source is usually approximated as an 
% isotropic positive source located at depth z_0 into the medium
z_b = (2*z_0*(1+R_eff))/(3*(1-R_eff)); % the boundary condition requirement 
% leads to a signal size of zero (i.e. for fluence rate in the case of a NIR
% diffuse reflectance measurement or for the temporal autocorrelation function
% in the case of the DCS measurement) at z = -z_b, which is generally called
% the extrapolated zero-boundary condition
r_1 = sqrt((z_0^2)+(rho^2)); %3D distance b/w detector and positive isotropic 
% imaging source
r_2 = sqrt(((2*z_b+z_0)^2)+(rho^2)); %3D distance b/w detector and negative
% isotropic imaging source at z = -(z_0+2z_b)
omega = 2*pi*f_a; % acoustic radial frequency in rads/s

% define function
k_D_1 = @(alpha_DB,tau) sqrt((3*mu_a*mu_s_p)+(mu_s_p^2*k_0^2*(6*alpha_DB)*tau));

% load data - NB this is experimental data from 
load Ancillary_scripts\Data\AOM_data.mat

%% unmodulated g_1
g_1_DCS = @(tau) abs(...
    ((r_2*exp(-1*k_D_1(alpha_Db_brownian,tau)*r_1) - r_1*exp(-1*k_D_1(alpha_Db_brownian,tau)*r_2)))./ ...
    ((r_2*exp(-1*k_D_1(alpha_Db_brownian,0)*r_1)   - r_1*exp(-1*k_D_1(alpha_Db_brownian,0)*r_2)))...
    ).^2;

%% modulated g_1
g_1_AOM_DCS = @(tau) abs( (1 - ...
    A/2 + A/2*cos(omega.*tau) ... % first harmonic only 
    ).*...
    ((r_2*exp(-1*k_D_1(alpha_Db_brownian,tau)*r_1) - r_1*exp(-1*k_D_1(alpha_Db_brownian,tau)*r_2)))./ ...
    ((r_2*exp(-1*k_D_1(alpha_Db_brownian,0)*r_1)   - r_1*exp(-1*k_D_1(alpha_Db_brownian,0)*r_2)))...
    ).^2;

%% plot

figure('units','normalized','outerposition',[0 0 1 1])
plot(tau,g_1_DCS(tau),'r','LineWidth',LW)
set(gca,'xscale','log')
hold on
plot(tau,g_1_AOM_DCS(tau),'g','LineWidth',LW)
hold off
ax = gca;
ax.FontSize = fs;

xlabel('$\tau$ (s)','FontSize',fs+10,'FontWeight','bold','Interpreter','Latex')
ylabel('$g_1(\tau)$','FontSize',fs+10,'FontWeight','bold','Interpreter','Latex')
legend('DCS - Brownian fit','AOM-DCS - Brownian fit',...
    'location','NorthEast','FontWeight','bold','Fontsize',fs+7,'Interpreter','Latex')
ylim([0 1])
lower_bound = 5e-7;
upper_bound = 5e-4;
xlim([lower_bound upper_bound])
hold off
box on

%% now get eigenvalues 
dim_1 = 600;
dim_2 = 600; 
radius = 100;
H = generateH(dim_1,dim_2,radius);
T = 10e-6; % ??
N = 500; % number of eigen values
alpha = 1;
t2 = linspace(0,T,N);
t1 = t2;
[n,k] = meshgrid(t2,t1);
diff = abs(k-n);

%% check exposure time
exp_time_check = 1/(T*omega);

% calculate eigen values for DCS
K_DCS = alpha*g_1_DCS(diff) + (1-alpha);
lambda_DCS = eigs(K_DCS,N)/N;

% calculate eigen values for AOM DCS
K_AOM_DCS = alpha*g_1_AOM_DCS(diff) + (1-alpha);
lambda_AOM_DCS = eigs(K_AOM_DCS,N)/N;


%% plot eigenvalues
% plot K matrices too

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
fs = 24;
imagesc(t1*1e6,t2*1e6,K_DCS)
ax = gca;
ax.FontSize = fs;
colormap bone
xlabel('$t_1$ ($\mu$s)','FontWeight','bold','Interpreter','Latex','Fontsize',fs+3)
ylabel('$t_2$ ($\mu$s)','FontWeight','bold','Interpreter','Latex','Fontsize',fs+3)
axis square
colorbar
title('(a)','Interpreter', 'latex',...
    'FontWeight','bold','Interpreter','Latex','Fontsize',fs+5)

subplot(2,2,2)
fs = 24;
imagesc(t1*1e6,t2*1e6,K_AOM_DCS)
ax = gca;
ax.FontSize = fs;
colormap bone
xlabel('$t_1$ ($\mu$s)','FontWeight','bold','Interpreter','Latex','Fontsize',fs+3)
ylabel('$t_2$ ($\mu$s)','FontWeight','bold','Interpreter','Latex','Fontsize',fs+3)
axis square
colorbar
title('(b)','Interpreter', 'latex',...
    'FontWeight','bold','Interpreter','Latex','Fontsize',fs+5)

subplot(2,1,2)
plot(1:N,lambda_DCS,'r--','LineWidth',LW)
hold on
plot(1:N,lambda_AOM_DCS,'g-','LineWidth',LW)
hold off
ax = gca;
ax.FontSize = fs;
set(gca,'YScale','log');
xlabel('$n$','FontWeight','bold','Interpreter','Latex','Fontsize',fs+3)
ylabel('$\lambda_n$','FontWeight','bold','Interpreter','Latex','Fontsize',fs+3)
title('(c)','FontWeight','bold','Interpreter','Latex','Fontsize',fs+5)
xlim([0 N])
box on
% add legend
legend('DCS','AOM-DCS',...
    'location','NorthEast','FontWeight','bold','Fontsize',fs+3,'Interpreter','Latex')

%% simulate images for each 
W_DCS     = zeros(dim_1,dim_2);
W_AOM_DCS = zeros(dim_1,dim_2);

for j = 1:N
    % calculate random field - uncorrelated and independent - this is b_n
    U = generateM(dim_1,dim_2);
    % calculate intensity - factor in CTF here
    I = calcIOneImage(U,H);
    % convert to unit mean
    I = I./mean(I(:));
    % convert to mean lambda_n
    W_DCS     = W_DCS     + I*lambda_DCS(j);
    W_AOM_DCS = W_AOM_DCS + I*lambda_AOM_DCS(j);
    disp(j)
end

%% show images and histograms
figure('units','normalized','outerposition',[0 0 1 1])
data_1 = W_DCS./mean(W_DCS(:));
subplot(1,3,1)
imshow(data_1(1:150,1:150),[])
title('(a)','FontWeight','bold','Interpreter','Latex','Fontsize',fs+5)

subplot(1,3,2)
data_2 = W_AOM_DCS./mean(W_AOM_DCS(:));
imshow(data_2(1:150,1:150),[])
title('(b)','FontWeight','bold','Interpreter','Latex','Fontsize',fs+5)

subplot(1,3,3)
histogram(data_1,'BinMethod','fd','DisplayStyle','stairs',...
    'LineWidth',LW,'EdgeColor','r','Normalization','pdf','LineStyle','--')
hold on
histogram(data_2,'BinMethod','fd','DisplayStyle','stairs',...
    'LineWidth',LW,'EdgeColor','g','Normalization','pdf')
hold off
ax = gca;
ax.FontSize = fs; 
title('(c)','FontWeight','bold','Interpreter','Latex','Fontsize',fs+5)

axis square
xlabel('$\mathbf{W}/\overline{\mathbf{W}}$','FontSize',fs,'FontWeight','bold','Interpreter','Latex')
ylabel('PDF','FontSize',fs,'FontWeight','bold','Interpreter','Latex')
xlim([0 3])
legend('DCS','AOM-DCS',...
    'FontWeight','bold','Fontsize',fs,'Interpreter','Latex')

%% calculate speckle contrast 
SK_DCS = std(data_1(:))/mean(data_1(:));
SK_AOM_DCS = std(data_2(:))/mean(data_2(:));

%% does it fit to theory of A
test_A = SK_DCS/SK_AOM_DCS - 1;
