function [tau,g_2,beta_brownian,alpha_Db_brownian,A,fit_brownian_unmod,fit_brownian_mod] ...
    = autocorrelationFitAOTCWf0(data,lower_bound,upper_bound,n,lambda,rho,mu_a,...
    mu_s_p,beta_range,alpha_Db_range,f_a,~)

% Function to fit measured g_2 CW AOT data to clinically relevant semi-infinite
% geometry solution of correlation diffusion equation (CDE). Uses same code
% from autocorrelationFit.m.

% Then simply modulates this by (1 + M*cos(w_a*t)) - therefore M is a
% coefficient representing modulation depth

% uses fundamental frequency only.

% Does Brownian model only - no random flow model at the moment.

% Inputs:
    % 1)  data - structure with two fields: autocorrelation g_2 'data' and
    %     experimental parameters in 'textdata'
    % 2)  lower bound - lower bound on tau
    % 3)  upper bound - upper bound on tau
    % 4)  n - refractive index of medium under examination
    % 5)  lambda - wavelength of laser (cm)
    % 6)  rho - 1D source detector spacing
    % 7)  mu_a - absorption coefficient (cm^-1)
    % 8)  mu_s_p - reduced scattering coefficient (cm^-1)
    % 9)  beta_range - range of beta values from which to randomly pick a
    %     starting point for optimisation
    % 10) alpha_Db_range - range of alpha_Db values from which to randomly
    %     pick a starting point for optimisation (cm^2/s)
    % 11) f_a - acoustic frequency in Hz
    % 12) ~ - if present then plots 3 figures, otherwise does not

% Outputs:
    % 1) tau - truncated tau data set
    % 2) g_2 - truncated g_2 data set
    % 3) beta_brownian - beta value from brownian fit model
    % 4) alpha_Db_brownian - alpha_Db value (cm^2/s)
    % 6) A - modulation depth
    % 7) fit_brownian_unmod - unmodulated output of brownian fit model
    % 8) fit_brownian_mod - modulated output of brownian fit model
    
% Plots 3 figures:
    % 1) Figure 1 - entire g_2 data set
    % 2) Figure 2 - truncated g_2 data set
    % 3) Figure 3 - raw data, brownian fit (modulated and unmodulated), 
    %    3 extracted parameters
    
% Integration time is displayed in command window. 
    
% Author: Edward James, PhD Student, UCL, November 2018.
% e.james.14@ucl.ac.uk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To set line widths and font size for plots
fs=20;
LW=1.5;

%extract and assign relevant fields
g_2 = data.data(:,2);
tau = data.data(:,1);
if nargin == 12
    % display integration time in seconds
    disp('Integration time in seconds:')
    data.textdata{4} %#ok<NOPRT>
end

% truncate at lower and minimum bounds for tau
[~, lower_index] = min( abs( tau-lower_bound ));
[~, upper_index] = min( abs( tau-upper_bound ));
tau = tau(lower_index:upper_index);
g_2 = g_2(lower_index:upper_index);

%  assign known parameters
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

% assign initial guesses for each of 3 variables:
beta = datasample(beta_range,1);
alpha_Db = datasample(alpha_Db_range,1);
A = 0; % amplitude modulation (i.e. at fundamental frequency)

% Brownian model

% specify objective function to minimise, i.e. the sum of square differences
% between measured data and analytical model, let x = [beta, alpha_Db]
% NB we are using the normalised elecric field autocorrelation function
% here, g_1, i.e. [G_1(tau)/G_1(0)]

k_D_1 = @(alpha_DB,tau) sqrt((3*mu_a*mu_s_p)+(mu_s_p^2*k_0^2*(6*alpha_DB)*tau));
% minimise the SSD between the measured and modelled functions
obj_fun_1 = @(x) sum((g_2 -(1 + (x(1)*abs( (1 - ...
    x(3)/2 + x(3)/2*cos(omega.*tau) ... % fundamental frequency 
    ) .*...   
    ((r_2*exp(-1*k_D_1(x(2),tau)*r_1) - r_1*exp(-1*k_D_1(x(2),tau)*r_2)))./ ...
    ((r_2*exp(-1*k_D_1(x(2),0  )*r_1) - r_1*exp(-1*k_D_1(x(2),0)*r_2)))...
    ).^2))).^2);

% assign start points and termination tolerance on fitted variables
x0 = [beta, alpha_Db, A];
options = optimset('TolX',1e-11,'MaxFunEvals',length(x0)*300);

% find x which minimises this objective function, also return residual at
% termination
[x,~,~,~] = fminsearch(obj_fun_1,x0,options);

%extract these optimised values
beta_brownian = x(1);
alpha_Db_brownian = x(2);
A = x(3);

% For case of diffusive motion, Db is the effective Brownian diffusion
% coefficient of the tissue scatterers.
% alpha is a unitless factor, which represents the fraction of light-scattering
% events from moving scatterers, units of area/time (here cm^2/s)

% calculate vector of fitted data for later plotting 
% unmodulated signal
fit_brownian_unmod = 1+(beta_brownian*abs(...
    ((r_2*exp(-1*k_D_1(alpha_Db_brownian,tau)*r_1) - r_1*exp(-1*k_D_1(alpha_Db_brownian,tau)*r_2)))./ ...
    ((r_2*exp(-1*k_D_1(alpha_Db_brownian,0)*r_1)   - r_1*exp(-1*k_D_1(alpha_Db_brownian,0)*r_2)))...
    ).^2);

% modulated signal
fit_brownian_mod = 1+(beta_brownian*abs( (1 - ...
    A/2 + A/2*cos(omega.*tau) ... % first harmonic only 
    ).*...
    ((r_2*exp(-1*k_D_1(alpha_Db_brownian,tau)*r_1) - r_1*exp(-1*k_D_1(alpha_Db_brownian,tau)*r_2)))./ ...
    ((r_2*exp(-1*k_D_1(alpha_Db_brownian,0)*r_1)   - r_1*exp(-1*k_D_1(alpha_Db_brownian,0)*r_2)))...
    ).^2);

if nargin == 12
    % Plot these two fitted models
    figure('units','normalized','outerposition',[0 0 1 1])
    scatter(tau,g_2,'b.')
    set(gca,'xscale','log')
    hold on
    plot(tau,fit_brownian_unmod,'r','LineWidth',LW)
    plot(tau,fit_brownian_mod,'g','LineWidth',LW)

    ax = gca;
    ax.FontSize = fs; 

    xlabel('$\tau$ (s)','FontSize',fs+10,'FontWeight','bold','Interpreter','Latex')
    ylabel('$g_2(\tau)$','FontSize',fs+10,'FontWeight','bold','Interpreter','Latex')
    xlim([lower_bound upper_bound])
    legend('Autocorrelated data','DCS - Brownian fit','AOM-DCS - Brownian fit',...
        'location','NorthEast','FontWeight','bold','Fontsize',fs+7,'Interpreter','Latex')
    ylim([0.99 1.3])
    hold off
    box on
end

end

