clc
clearvars
close all

%% add functions to file_path
addpath(pwd,"Functions\")

%% NB
% run script generate_standard_data.m prior to running this script

% load_data = true;
load_data = false;

%% load data
load Data\comp_data
clear delta N_s radius T M_1 M_2
% linewidth and fontsize
LW = 1.5;
fs = 15;

%% specify parameters

N = 500; % number of eigen values
alpha = 0.90;
log_vec = 10.^(-3:0.1:3);

T = log_vec*t_c;
num_Ts = length(T);
T_fine = linspace(T(1),T(end),m);
W_stack = zeros(num_Ts,dim_1,dim_2);

%% calculate M values
M = Mtwo(T,t_c);

if load_data == false
    for i = 1:num_Ts
        % calculate eigen values
        t2 = linspace(0,T(i),N);
        t1 = t2;
        [n,k] = meshgrid(t2,t1);
        diff = abs(k-n);
        K = alpha*exp(-(diff/t_c).^2) + (1-alpha);
        lambda = eigs(K,N)/N;
        
        % loop through and take weighted sum of images
        W = zeros(dim_1,dim_2);
        for j = 1:N
            U = generateM(dim_1,dim_2);
            % add CTF here
            I = calcIOneImage(U,H);
            I = I./mean(I(:));
            W = W + I*lambda(j);
            disp([i,j])
        end
        W_stack(i,:,:) = W;
    end
    save case_5_data_T_vec
elseif load_data == true
    load case_5_data_T_vec
end

%% add code to reanalyse K using a different technique if need be and also get more metrics
N_s_vec = 1;
num_N_s = length(N_s_vec);
k_vec_temp = zeros(1,num_Ts);
for i = 1:num_Ts
    k_vec_temp(i) = calcGlobalK(squeeze(W_stack(i,:,:)));
    disp(i)
end

%% calculate K_model fit
for m = 1:num_N_s
    N_s = N_s_vec(m);
    obj_fun = @(x) sum((k_vec_temp-1./sqrt(Mtwo(T,x(1),x(2),x(3)))).^2);
    x0 = [1e-3,1,1];
    lb = [0,0,0];
    ub = [inf,1,1];
    % assign termination tolerance on fitted variables
    options = optimset('TolX',1e-11,'Display','off');
    x = fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],options);
    t_c_fit   = x(1);
    beta_fit  = x(2);
    alpha_fit = x(3); 
    
    M_fine = Mtwo(T_fine,t_c_fit,beta_fit,alpha_fit);
    K_model_plot = 1./sqrt(M_fine);
    
    %% plot results
    figure('units','normalized','outerposition',[0 0 1 1])
    scatter(T,k_vec_temp,'b+','LineWidth',1.5)
    set(gca, 'XScale', 'log')
    hold on
    semilogx(T_fine,K_model_plot,'b','LineWidth',1)
    % show asymptotes
    yline(sqrt(beta_fit),'k-.','LineWidth',1.5)
    yline(sqrt(beta_fit)*(1-alpha_fit),'k--','LineWidth',1.5)
    legend('Simulated',...
        ['Model Fit - $\tau_c$ = ',num2str(t_c_fit*1e3,'%.2f'),' ms - $\beta$ = ',...
        num2str(beta_fit,'%.2f'),' - $\alpha$ = ',num2str(alpha_fit,'%.2f')],...
        '$K = \sqrt{\beta}$','$K = \sqrt{\beta}(1-\alpha)$',...
        'location','best',...
        'Interpreter','Latex')
    xlabel('Camera Integration Time, $T$ (seconds)','Interpreter','Latex')
    ylabel('Global Speckle Contrast, $K$','Interpreter','Latex')
    grid minor
    title (['Case 5 - $g_{1d}(\tau) = \exp\left(-\left(\frac{\tau}{\tau_c}\right)^2\right)$ - $\tau_c$ = ',...
    num2str(t_c*1e3),' ms - CTF Modelled - $\alpha$ = ',num2str(alpha)],'Interpreter','Latex')
    xlim([T(1) T(end)])
    ylim([0 1])
end
save case_5_data_T_vec