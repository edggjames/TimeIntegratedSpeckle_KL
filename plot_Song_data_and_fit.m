clc
close all hidden
clearvars

%% NB
% run script generate_data_Song_et_al.m prior to running this script

%% load generated data and plot
load Ancillary_scripts\Data\Song_data.mat

% define M for Lorentzian spectrum, i.e. neg exp autocorrelation
r = T./t_c;
M = (1./r + 0.5*r.^-2.*(exp(-2*r)-1)).^-1;

figure('units','normalized','outerposition',[0 0 1 1])
scatter(r,sqrt(K_sq),'k+','Linewidth',LW)
set(gca, 'XScale', 'log')

hold on
xlabel('$T/\tau_c$','FontSize',fs,'FontWeight','bold','Interpreter','Latex')
ylabel('Global speckle contrast, $K$','FontSize',fs,'FontWeight','bold','Interpreter','Latex')

%% develop forward model from Parthasarathy et al 2008 - "Robust flow 
% measurement with multi-exposure speckle imaging."

obj_func = @(x) sum((K_sq-parthasarathyMethod(T,t_c,x)).^2);
x0 = 0;
x = fminsearch(obj_func,x0);
beta = x;
K_sq_model = parthasarathyMethod(T,t_c,1);
semilogx(r,sqrt(K_sq_model),'k','Linewidth',LW)
xlim([min(r) r(56)])
ylim([0 1.1])
box on
legend('Simulated Data','Expected Data','location','SouthWest',...
    'FontSize',fs)
ax = gca;
ax.FontSize = fs; 


function K_sq_model = parthasarathyMethod(T,t_c,beta)
% equation 11 from paper with rho = 1
W = T./t_c;
K_sq_model = beta*(exp(-2*W)-1+2*W)./(2*W.^2);
end


