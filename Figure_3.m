% show validation on all three cases
clc
clearvars
close all

%% Figure 3 from Edward James, Samuel Powell, and Peter Munro, "Simulation
% of statistically accurate time-integrated dynamic speckle patterns in
% biomedical optics," Opt. Lett. 46, 4390-4393 (2021)

%% NB
% run scripts
%  a)generate_standard_data.m,
%  b)Case_3_calcK_and_Plot.m
%  c)Case_4_calcK_and_Plot.m
%  d)Case_5_calcK_and_Plot.m
% prior to running this script

%% add functions to file_path
addpath(pwd,"Ancillary_scripts\Functions\")

%% NB
% run script generate_standard_data.m prior to running this script

% linewidth and fontsize
LW = 1.5;
fs = 15;

%% load data from case 3
load Ancillary_scripts\Data\case_3_data_T_vec.mat
case_3_K_model = K_model_plot;
case_3_K_sim   = k_vec_temp;
t_c_fit_3      = t_c_fit;
beta_fit_3     = beta_fit;
alpha_fit_3    = alpha_fit;

%% load data from case 4
load Ancillary_scripts\Data\case_4_data_T_vec.mat
case_4_K_model = K_model_plot;
case_4_K_sim   = k_vec_temp;
t_c_fit_4      = t_c_fit;
beta_fit_4     = beta_fit;
alpha_fit_4    = alpha_fit;

%% load data from case 5
load Ancillary_scripts\Data\case_5_data_T_vec.mat
case_5_K_model = K_model_plot;
case_5_K_sim   = k_vec_temp;
t_c_fit_5      = t_c_fit;
beta_fit_5     = beta_fit;
alpha_fit_5    = alpha_fit;

%% plot
figure('units','normalized','outerposition',[0 0 1 1])
set(gca, 'XScale', 'log')
ax = gca;
ax.FontSize = fs; 
hold on
% show simulated data and model fits
scatter(T/t_c,case_4_K_sim,'g','LineWidth',LW)
semilogx(T_fine/t_c,case_4_K_model,'g-.','LineWidth',LW)
scatter(T/t_c,case_3_K_sim,'rd','LineWidth',LW)
semilogx(T_fine/t_c,case_3_K_model,'r:','LineWidth',LW)
scatter(T/t_c,case_5_K_sim,'bs','LineWidth',LW)
semilogx(T_fine/t_c,case_5_K_model,'b--','LineWidth',LW)
hold off

legend('$p=0.5$ - simulated data',...
    ['Model fit - $\tau_c$ = ',num2str(t_c_fit_4*1e3,'%.2f'),' ms - $\beta$ = ',...
    num2str(beta_fit_4,'%.2f'),' - $\alpha$ = ',num2str(alpha_fit_4,'%.2f')],...
    '$p=1$ - simulated data',...
    ['Model fit - $\tau_c$ = ',num2str(t_c_fit_3*1e3,'%.2f'),' ms - $\beta$ = ',...
    num2str(beta_fit_3,'%.2f'),' - $\alpha$ = ',num2str(alpha_fit_3,'%.2f')],...
    '$p=2$ - simulated data',...
    ['Model fit - $\tau_c$ = ',num2str(t_c_fit_5*1e3,'%.2f'),' ms - $\beta$ = ',...
    num2str(beta_fit_5,'%.2f'),' - $\alpha$ = ',num2str(alpha_fit_5,'%.2f')],...
    'location','southwest',...
    'Interpreter','Latex','FontWeight','bold','Fontsize',fs+8)
xlabel('$T/\tau_c$','FontWeight','bold','Interpreter','Latex','Fontsize',fs+8)
ylabel('Global speckle contrast, $K_p$','Interpreter','Latex','Fontsize',fs+8)
ylim([0 1])
xlim([T(1)/t_c T(end)/t_c])
box on 