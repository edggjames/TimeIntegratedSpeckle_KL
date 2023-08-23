clc
close all
clearvars

%% Figure S2 from Edward James, Samuel Powell, and Peter Munro, "Simulation
% of statistically accurate time-integrated dynamic speckle patterns in
% biomedical optics," Opt. Lett. 46, 4390-4393 (2021)

%% NB
% run script generate_standard_data.m prior to running this script

%% add functions to file_path
addpath(pwd,"Ancillary_scripts\Functions\")

% load_data = true;
load_data = false;

%% To set line widths and font size for plots
fs=20;
LW=1.5;

load Ancillary_scripts\Data\comp_data
clear delta M_1 M_2 radius H dim_1 dim_2 m N_s

log_vec = 10.^(-2:0.1:3);
T = t_c.*log_vec;
num_Ts = length(T);
M = Mone(T,t_c);
rGood = defrGood(T,t_c);
N = 1000;
alpha = 1;

if load_data == false
    
    e_vec_1 = zeros(1,num_Ts);
    e_vec_2 = zeros(1,num_Ts);
    e_vec_3 = zeros(1,num_Ts);
    e_vec_4 = zeros(1,num_Ts);
    e_vec_5 = zeros(1,num_Ts);
    e_vec_6 = zeros(1,num_Ts);
    e_vec_7 = zeros(1,num_Ts);
    
    eig_sum = zeros(1,num_Ts);
    
    %%  loop through T vals
    for i = 1:num_Ts
        
        T_temp = T(i);
        t2 = linspace(0,T_temp,N);
        t1 = t2;
        [n,k] = meshgrid(t2,t1);
        diff = (k-n);
        K = alpha*exp(-abs(diff)/t_c) + (1-alpha);  
        
        e = eigs(K,N)/N;
        
        eig_sum(i) = sum(e);
        e_vec_1(i) = e(1);
        e_vec_2(i) = e(2);
        e_vec_3(i) = e(3);
        e_vec_4(i) = e(4);
        e_vec_5(i) = e(5);
        disp(i)
    end
    save Ancillary_scripts\Data\figure_S1_data
else
    load Ancillary_scripts\Data\figure_S1_data
end

%% figure for letters paper
figure('units','normalized','outerposition',[0 0 1 1])
loglog(rGood,e_vec_1,'k','Linewidth',LW);
hold on
loglog(rGood,e_vec_2,'k','Linewidth',LW);
loglog(rGood,e_vec_3,'k','Linewidth',LW);
loglog(rGood,e_vec_4,'k','Linewidth',LW);
loglog(rGood,e_vec_5,'k','Linewidth',LW);
ax = gca;
ax.FontSize = fs; 
box on
hold off
xlabel('$T/\tau_c$','FontWeight','bold','Interpreter','Latex','Fontsize',fs+10)
ylabel('$\lambda_n$','FontWeight','bold','Interpreter','Latex','Fontsize',fs+10)
xlim([min(rGood) max(rGood)])
ylim([1e-4 1.2])
% add labels 
x = 2.5e-2;
text(x,7.0e-1,'$n=1$','Interpreter','Latex','Fontsize',fs+5,'HorizontalAlignment','right')
text(x,7.0e-3,'$2$','Interpreter','Latex','Fontsize',fs+5)
text(x,1.7e-3,'$3$','Interpreter','Latex','Fontsize',fs+5)
text(x,7.8e-4,'$4$','Interpreter','Latex','Fontsize',fs+5)
text(x,4.2e-4,'$5$','Interpreter','Latex','Fontsize',fs+5)