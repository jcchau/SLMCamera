function lb_nats = computeCapacityLBUnifGx(G, x_max, sigma_w)
%COMPUTECAPACITYLBUNIFGX Summary of this function goes here
%   Detailed explanation goes here

[n_r, n_t] = size(G);

sigma_w = repmat(sigma_w, n_r, 1);

max_nbins_linprog = 6e2; % A tradeoff between speed and precision: 6e4
max_nbins_mem = 1e9/8/5;

% TODO: move this max_nbins_linprog feature in to MIMOCapacityLBUnifGx.
Gx_max = G*repmat(x_max, n_t, 1);

expected_nbins_linprog = ceil(prod( ...
    max_nbins_linprog^(1/n_r) .* (12*sigma_w+Gx_max) ./ Gx_max ));

max_nbins = min(max_nbins_mem, expected_nbins_linprog);

lb_nats = MIMOCapacityLBUnifGx.calculateMutualInfoWithUnifGx( ...
    G, x_max, sigma_w, max_nbins);

end

