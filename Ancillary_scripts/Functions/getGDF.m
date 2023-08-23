function [GDF] = getGDF(M,I_bar,I_vec)
% Function to generate Gamma density function of order M

GDF = (M/I_bar)^M*I_vec.^(M-1).*exp(-M*I_vec/I_bar)/gamma(M);
end

