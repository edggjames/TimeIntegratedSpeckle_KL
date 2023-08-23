clc
clearvars
close all

% load_data = true;
load_data = false;

%% add functions to file_path
addpath(pwd,"Ancillary_scripts\Functions\")

%% load data
load Ancillary_scripts\Data\comp_data
clear delta radius T M_1 M_2 m N_s t_c
% linewidth and fontsize
LW = 1.5;

%% load seed image
seed_image = phantom('Modified Shepp-Logan',dim_1);
labels = unique(seed_image(:));
n_labels = length(labels);

%% define alpha
alpha = [1,1,1,1,1,0.8,0.8];

%% define p
p = [1,1,1,1,1,1,0.5];

%% define tau_c
t_c = [0.05,0.10,0.20,0.50,1.00,1.00,1.00]*1e-3;

%% produce a label template for each label
labels_stencil = zeros(n_labels,dim_1,dim_2);

for i = 1:n_labels
    % make template for this label
    label_temp = zeros(dim_1,dim_2);
    label_temp(seed_image==labels(i)) = 1;
    % assign
    labels_stencil(i,:,:) = label_temp;
    clear label_temp
end

%%
if load_data == false
    
    %% then simulate 6 images (one for each label)
    N = 1000; % number of eigen values
    n_images = 30;
    W = zeros(n_images,dim_1,dim_2);
    T = 20e-3;
    t2 = linspace(0,T,N);
    t1 = t2;
    [n,k] = meshgrid(t2,t1);
    diff = abs(k-n);
    
    %% calculate matrix of eigenvalues
    lambda = zeros(n_labels,N);
    % loop through labels
    tic
    for j = 1:n_images
        for i = 1:n_labels
            % calculate eigen values
            K = alpha(i)*exp(-(diff/t_c(i)).^(p(i))) + (1-alpha(i));
            lambda(i,:) = eigs(K,N)/N;
            disp([j,i])
        end
    end
    time_1 = toc;
    time_1 = time_1/n_images;
    %%
    tic
    for k = 1:n_images
        % add N fully developed speckle images
        for j = 1:N
            % calculate random field - uncorrelated and independent - this is b_n
            U = generateM(dim_1,dim_2);
            % calculate intensity - factor in CTF here
            I = calcIOneImage(U,H);
            % convert to unit mean
            I = I./mean(I(:));
            % convert to mean lambda_n for each label
            for i = 1:n_labels
                W(k,:,:) = squeeze(W(k,:,:)) + I*lambda(i,j).*squeeze(labels_stencil(i,:,:));
            end
            disp([k,j])            
        end
    end
    time_2 = toc;
    time_2 = time_2/n_images;
    
    save Ancillary_scripts\Data\figure_4_data
elseif load_data
    load Ancillary_scripts\Data\figure_4_data
end

%% calculate spatial speckle contrast
n_images_K = n_images;
[K_im] = methodLSTCA(W,n_images_K);

%% plot Figure 4
fs = 20;
idx = 1;
x1 = 250;
w  = 200;
y1 = 150;
h  = 200;
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
imshow(seed_image,[])
ax = gca;
ax.FontSize = fs-idx;
colormap(gray(n_labels));
c = colorbar;
c.Ticks = linspace(0+1/(2*n_labels),1-1/(2*n_labels),n_labels);
c.TickLabels = ["1" "2" "3" "4" "5" "6" "7"];
box on
title('(a)','FontWeight','bold','Interpreter','Latex','Fontsize',fs)

subplot(1,3,2)
imshow(squeeze(W(1,:,:)),[])
ax = gca;
ax.FontSize = fs-idx;
colorbar
hold on
hold on
rectangle('Position',[x1 y1 w h],'EdgeColor','w','LineStyle','--')
box on
title('(b)','FontWeight','bold','Interpreter','Latex','Fontsize',fs)

subplot(1,3,3)
imshow(K_im,[])
ax = gca;
ax.FontSize = fs-idx;
colorbar
title('(c)','FontWeight','bold','Interpreter','Latex','Fontsize',fs)
box on

%% plot Figure S2
figure('units','normalized','outerposition',[0 0 1 1])
imshow(squeeze(W(1,:,:)),[])
ax = gca;
ax.FontSize = fs-idx;
% colorbar
box on
xlim([x1 x1+w-1])
ylim([y1 y1+h-1])