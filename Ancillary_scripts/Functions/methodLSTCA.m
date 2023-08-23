function [K_im] = methodLSTCA(W_stack,n_images)
% Function to implement LSTCA on stack of images

sigma = squeeze(std(W_stack(1:n_images,:,:),0,1));
mu    = squeeze(mean(W_stack(1:n_images,:,:),1));
K_im  = sigma./mu;

end

