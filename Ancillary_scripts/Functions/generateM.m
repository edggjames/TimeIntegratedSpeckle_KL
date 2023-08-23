function [M,omega] = generateM(dim_1,dim_2)
% function to generate M matrix from Song et al. 2016 

omega = rand(dim_1,dim_2)*2*pi-pi;
M = exp(-1i*omega);
end

