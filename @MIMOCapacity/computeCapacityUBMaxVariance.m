function ub_nats = computeCapacityUBMaxVariance(G, x_max, sigma_w)
% computeCapacityUBMaxVariance computes an upper bound on the capacity of
% y=G*x+w for x in [0,x_max] by using a Gaussian-distributed RV for each
% dimension of x, where the RV's variance is the maximum possible variance
% given the constraint that the RV is in [0,x_max].  
%
% Algorithm described in lab book #4, p77-79.  
%
%   ub_nats = computeCapacityUBMaxVariance(G, x_max, sigma_w)
%
% ub_nats is the calculated upper bound on capacity.
%
% G is the channel matrix.
% x_max (scalar) is the maximum possible value for each dimension of x.
% sigma_w (scalar) is the standard deviation of the AWGN for each receiver.

[n_r, ~] = size(G);

% maximum variance and standard deviation for each dimension of x.
max_var = (x_max/2)^2;

% Covariance of G*x, w, and y.
cov_Gx = max_var * (G*G');
cov_w = sigma_w^2 * eye(n_r);
cov_y = cov_Gx + cov_w;

h_y = MIMOCapacity.calculateDiffEntropyOfMVGaussian(cov_y);
h_y_given_x = MIMOCapacity.calculateDiffEntropyOfMVGaussian(cov_w);

ub_nats = h_y - h_y_given_x;

end

