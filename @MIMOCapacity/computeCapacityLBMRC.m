function lb_nats = computeCapacityLBMRC(G, x_max, sigma_w)
% computeCapacityLBMRC computes a lower bound on capacity of y=G*x+w by
% reducing performing maximal ratio combining at the receiver to reduce the
% received constellation to 1D and applying SmithCapacity.
%
%   lb_nats = computeCapacityLBMRC(G, x_max, sigma_w)
%
% lb_nats is a lower bound on capacity.
%
% G is the channel matrix.
% x_max (scalar) is the maximum value of x.
% sigma_w (scalar) is the the standard deviation of the noise w (along each
%   dimension). 

% Don't bother running simplifyChannelMatrix.
% The simplifications would not speed up this method.  

[~, n_t] = size(G);

Gx_max = G * repmat(x_max, n_t, 1);

mrc_xmax = sqrt(Gx_max' * Gx_max);

% Since sigma_w is the same along each dimension, the Gaussian noise has
% circular symmetry, so the standard deviation of the noise along the line
% from 0 to Gx_max is still sigma_w.  

lb_nats = MIMOCapacity.computeSmithScalarCapacity(mrc_xmax, sigma_w);

end

