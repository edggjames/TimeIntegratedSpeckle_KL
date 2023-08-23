function [K] = calcGlobalK(I)
% Calculate global speckle contrast of an image

K = std(I(:))/mean(I(:));
end

