function [H] = generateH(dim_1,dim_2,radius)
% Generate CTF as pupil plane in the frequency domain. 
% Defined as circular pupil in centre of image

% add options for other shapes (e.g. rectangle)

[Y,X] = meshgrid(1:dim_1,1:dim_2);
y_centre = dim_1/2;
x_centre = dim_2/2;
H = sqrt(((X-x_centre).^2) + (Y-y_centre).^2) <= radius;
end

