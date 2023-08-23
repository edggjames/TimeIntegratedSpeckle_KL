clc
clearvars
close all

%% Figure 2 from Edward James, Samuel Powell, and Peter Munro, "Simulation
% of statistically accurate time-integrated dynamic speckle patterns in
% biomedical optics," Opt. Lett. 46, 4390-4393 (2021)

% load_data = true;
load_data = false;

%% add functions to file_path
addpath(pwd,"Ancillary_scripts\Functions\")

%% NB
% run script generate_standard_data.m prior to running this script

%% load data
load Ancillary_scripts\Data\comp_data
clear delta N_s radius T M_1 M_2
% linewidth and fontsize
LW = 1.5;
fs = 15;

%% specify parameters
N = 500; % number of eigen values
log_vec = [0.1,1,5,25];

T = t_c.*log_vec;
num_Ts = length(T);
thresh = 1e-5; % for Goodman exact PDFs
res = m;
alpha = 0.90;
W_stack = zeros(num_Ts,dim_1,dim_2);
beta = 0.960612858930412; % from Case_3_calcK_and_Plot.m
rGood = defrGood(T,t_c);
P_states = 1;
Mone = Mone(T,t_c,beta,alpha);
Mone = Mone*P_states;
I_vec = linspace(0,20,res);
case_number = 3;

if load_data == false
    for i = 1:num_Ts
        W = zeros(dim_1,dim_2);
        % loop through all polarisation states
        for m = 1:P_states
            % calculate eigen values
            t2 = linspace(0,T(i),N);
            t1 = t2;
            [n,k] = meshgrid(t2,t1);
            diff = abs(k-n);
            K = alpha*exp(-diff/t_c) + (1-alpha);
            lambda = eigs(K,N)/N;
            
            for j = 1:N
                % calculate random field - uncorrelated and independent - this is b_n
                U = generateM(dim_1,dim_2);
                % calculate intensity - factor in CTF here
                I = calcIOneImage(U,H);
                % convert to unit mean
                I = I./mean(I(:));
                % convert to mean lambda_n
                W = W + I*lambda(j);
                disp([i,m,j])
            end
        end
        % assign
        W_stack(i,:,:) = W;
    end
    save Ancillary_scripts\Data\4_plot_data
elseif load_data
    load Ancillary_scripts\Data\4_plot_data
end

%% call function for exact solution
[pW_all,W_bar] = exactSolutionCases345(case_number,T,N,t_c,thresh,res,I_vec,alpha);

%% show all images and histograms
close all
figure('units','normalized','outerposition',[0 0 1 1])
title_vec = {'(a)','(c)','(e)','(g)','(b)','(d)','(f)','(h)'};
x1 = 1;
y1 = 1;
w = 150;
h = 150;

for i = 1:num_Ts
    data = squeeze(W_stack(i,:,:));
    data = data./mean(data(:));
    % show image
    subplot(2,4,i)
    imshow(data,[])
    xlim([x1, x1+w])
    ylim([y1, y1+h])
    ax = gca;
    ax.FontSize = fs; 
    title(title_vec(i),'FontWeight','bold','Interpreter','Latex','Fontsize',fs+5)
    % show histogram
    subplot(2,4,i+4)
    shade = 0.6;
    histogram(data,'BinMethod','auto','DisplayStyle','bar',...
        'EdgeColor','none','Normalization','pdf','FaceColor',[shade shade shade],...
        'FaceAlpha',shade)
    hold on
    plot(I_vec/W_bar(i),pW_all(i,:)*W_bar(i),'k-','LineWidth',LW)
    I_bar = W_bar(i);
    GDF = getGDF(Mone(i),I_bar,I_vec);
    plot(I_vec/I_bar,I_bar*GDF,'k--','LineWidth',LW)
    ax = gca;
    ax.FontSize = fs; 
    hold off
    title(title_vec(i+4),'FontWeight','bold','Interpreter','Latex','Fontsize',fs+2)
    xlabel('$\mathbf{W}/\overline{\mathbf{W}}$','Interpreter','Latex','FontWeight','bold','Interpreter','Latex','Fontsize',fs)
    ylabel('PDF','FontWeight','bold','Fontsize',fs,'Interpreter','Latex')
    axis square
    xlim([0 3])
    if i == 1
        legend('Simulated data','Exact expected','Approximate expected',...
            'FontWeight','bold','Fontsize',fs-4.5,'Interpreter','Latex')
    end
    max_val = max([max(pW_all(i,:)*W_bar(i)),max(I_bar*GDF)]);
    max_val = ceil(max_val*10)/10;
    ylim([0 max_val])
    box on 
end