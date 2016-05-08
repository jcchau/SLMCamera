function [lb_nats_conservative, lb_nats] = ...
    computeCapacityLBUnifX(G, x_max, sigma_w)
%COMPUTECAPACITYLBUNIFX Summary of this function goes here
%   Detailed explanation goes here

[n_r, ~] = size(G);

variance_noise_out = repmat(sigma_w^2, n_r, 1);

max_nbins = 1e9/8/5;
min_trials = 1e9;

[lb_nats, lb_nats_conservative] = ...
    MIMOCapacityLBUnifX.calculateCapacityLBForUniformInput( ...
    G, x_max, variance_noise_out, max_nbins, min_trials);

end

