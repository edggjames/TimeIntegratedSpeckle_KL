function [M_two,M_two_asymp] = Mtwo(T,t_c,beta,alpha)
%Function to produce M values for n = 2

if ~exist('beta','var')
    % if beta is not parsed set to default value of 1
    beta = 1;
end
rGood = defrGood(T,t_c);
if ~exist('alpha','var')
    % if alpha is not parsed set to default value of 1
    term_1 = exp(-2*rGood.^2)-1;
    term_2 = erf(rGood*sqrt(2)).*rGood*sqrt(2*pi);
    term_3 = 2*rGood.^2;
    M_twp  = (term_1+term_2)./(term_3);
    M_two  = 1/beta*(M_twp).^-1;
else
    % call function with no alpha or beta for first term and invert
    term_1 = (alpha^2)*Mtwo(T,t_c).^-1;
    % then calculate other terms
    term_2 = 2*alpha*(1-alpha);
    term_3 = (exp(-rGood.^2)-1+erf(rGood).*rGood*sqrt(pi))./(rGood.^2);
    term_4 = (1-alpha)^2;
    M_two  = 1/beta*(term_1+term_2.*term_3+term_4).^-1;
    M_two_asymp = 1/beta*rGood.*sqrt(2/pi);
end

end

