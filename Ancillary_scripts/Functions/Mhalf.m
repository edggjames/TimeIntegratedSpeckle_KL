function [M_half,M_half_asymp] = Mhalf(T,t_c,beta,alpha)
%Function to produce M values for n = 1/2

if ~exist('beta','var')
    % if beta is not parsed set to default value of 1
    beta = 1;
%     disp('Setting beta to default value of 1 in Mhalf')
end
rGood = defrGood(T,t_c);
if ~exist('alpha','var')
    % if alpha is not parsed set to default value of 1
%     disp('Setting alpha to default value of 1 in Mhalf')
    term_1 = 3*exp(-2*sqrt(rGood));
    term_2 = 1+2*sqrt(rGood)+4/3*rGood;
    term_3 = -3 + 2*rGood;
    term_4 = 2*rGood.^2;
    M_half = 1/beta*((term_1.*term_2+term_3)./term_4).^-1;
    M_half_asymp = T./(beta*t_c);
else
    % call function with no alpha or beta for first term and invert
    term_1 = (alpha^2)*Mhalf(T,t_c).^-1;
    % then calculate other terms
    term_2 = 16*alpha*(1-alpha);
    term_3 = (3*exp(-sqrt(rGood)).*(1+sqrt(rGood)+rGood/3)-3+rGood/2)./(rGood.^2);
    term_4 = (1-alpha)^2;
    M_half = 1/beta*(term_1+term_2.*term_3+term_4).^-1;
    M_half_asymp = T./(beta*t_c);
end

end

