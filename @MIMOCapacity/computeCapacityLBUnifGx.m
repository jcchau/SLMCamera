function lb_nats = computeCapacityLBUnifGx(G, x_max, sigma_w)
%COMPUTECAPACITYLBUNIFGX Summary of this function goes here
%   Detailed explanation goes here

%% Prepare the channel matrix

G = MIMOCapacity.simplifyChannelMatrix(G);

% Handle special case where G is 0.
if(isequal(G,0))
    lb_nats = 0;
    return
end

% Apply Q' transform
Q = MIMOCapacity.computeQTTransform(G);
G = Q' * G;

%% Compute the lower bound on capacity

n_r = size(G, 1);

sigma_w = repmat(sigma_w, n_r, 1);

max_nbins_linprog = 6e4; % A tradeoff between speed and precision: 6e4
max_nbins_mem = 1e9/8/5;

% TODO: move this max_nbins_linprog feature in to MIMOCapacityLBUnifGx.
[umin, umax] = MIMOCapacity.computeUExtremes(G, x_max);
usize = umax-umin;

expected_nbins_linprog = ceil(prod( ...
    max_nbins_linprog^(1/n_r) .* (12*sigma_w+usize) ./ usize ));

max_nbins = min(max_nbins_mem, expected_nbins_linprog);

lb_nats = MIMOCapacityLBUnifGx.calculateMutualInfoWithUnifGx( ...
    G, x_max, sigma_w, max_nbins);

end

