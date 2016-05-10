function lb_nats = computeCapacityLBUnifGx(G, x_max, sigma_w)
%COMPUTECAPACITYLBUNIFGX Summary of this function goes here
%   Detailed explanation goes here

[n_r, ~] = size(G);

sigma_w = repmat(sigma_w, n_r, 1);

max_nbins = 1e7; % A tradeoff between speed and precision

lb_nats = MIMOCapacityLBUnifGx.calculateMutualInfoWithUnifGx( ...
    G, x_max, sigma_w, max_nbins);

end

