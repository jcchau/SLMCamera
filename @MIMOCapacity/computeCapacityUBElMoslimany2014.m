function ub_nats = computeCapacityUBElMoslimany2014(G, x_max, sigma_w)
% computeCapacityUBElMoslimany2014 computes an upper bound on the capacity
% of y=G*x+w for x in [0,x_max] using the method described in
% ElMoslimany2014.
%
%   ub_nats = computeCapacityUBElMoslimany2014(G, x_max, sigma_w)
%
% ub_nats is the upper bound on capacity in nats.
%
% G is the channel matrix.
% x_max (scalar) is the maximum value of x.
% sigma_w (scalar) is the standard deviation of the noise w.

% Don't bother running simplifyChannelMatrix.
% The simplifications would not noticeably speed up this method.  

[~, n_t] = size(G);

Gx_max = G * repmat(x_max, n_t, 1);

C_for_each_dim_of_Gx = arrayfun( ...
    @(xm) MIMOCapacity.computeSmithScalarCapacity(xm, sigma_w), Gx_max);

ub_nats = sum(C_for_each_dim_of_Gx);

end

