function [I] = calcIOneImage(M,H)
% Function to calculate I from one seed image

% calculate intensity distribution at the image plane
if nargin == 1
    I = abs(fftshift(fft2(M))).^2;
elseif nargin == 2
    % factor in H too
    I = abs(ifftshift(ifft2(H.*fftshift(fft2(M))))).^2;
end

end