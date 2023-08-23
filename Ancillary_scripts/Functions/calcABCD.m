function [a,b,c,d] = calcABCD(M_1,M_2,H)
% function to calculate the values of a, b, c and d for two uncorrelated
% speckle patterns

% H is the CTF of the instrument, if passed in the appropriate calculation
% is done. 

% calculate FFTs
A = fftshift(fft2(M_1));
B = fftshift(fft2(M_2));
if nargin == 3
    A = ifftshift(ifft2(A.*H));
    B = ifftshift(ifft2(B.*H));    
    % disp('CTF in calcABCD')
else
    % warning('No CTF in calcABCD')
end
% take real and imaginary components
a = real(A);
b = imag(A);
c = real(B);
d = imag(B);
end

