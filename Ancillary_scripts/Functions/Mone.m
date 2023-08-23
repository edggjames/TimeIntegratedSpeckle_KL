function [M_one,M_one_asymp] = Mone(T,t_c,beta,alpha)
%Function to produce M values for n = 1

if ~exist('beta','var')
    % if beta is not parsed set to default value of 1
    beta = 1;
end
rGood = defrGood(T,t_c);
if ~exist('alpha','var')
    % if alpha is not parsed set to default value of 1
    term_1 = exp(-2*rGood);
    term_2 = term_1-1+(2*rGood);
    M_one = 1/beta*(term_2./(2*rGood.^2)).^-1;
    M_one_asymp = rGood/beta;
else
    % call function with no alpha or beta for first term and invert
    term_1 = (alpha^2)*Mone(T,t_c).^-1;
    % then calculate other terms
    term_2 = 4*alpha*(1-alpha);
    term_3 = (exp(-rGood)-1+rGood)./(rGood.^2);
    term_4 = (1-alpha)^2;
    M_one = 1/beta*(term_1+term_2.*term_3+term_4).^-1;
    M_one_asymp = rGood/beta;
end

end

