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

% Maximum number of bins to use for the PMF.
% Primarily constrained by the amount of memory available.  
% Here, we're assuming 1GB, with each PMF bin taking up 8 bytes, and
% expecting up to 5 matrices of the size of the PMF.  
% Can increase max_nbins when we transition to the SCC computers with more
% memory.  
max_nbins = 1e9/8/5;

lb_nats = MIMOCapacityLBUnifGx.calculateMutualInfoWithUnifGx( ...
    G, x_max, sigma_w, max_nbins);

end

