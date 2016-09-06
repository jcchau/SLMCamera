function testComputeCapacityLBUnifXFullRankHighSNR(tc)
% testComputeCapacityLBUnifXFullRankComparison uses a randomly-generated
% full column rank channel matrix to verify that the lower-bound on
% capacity calculated via the full column rank branch (using
% MIMOCapacityLBUnifX.approximateUnifXLBForNegligibleNoiseAndFullColRankG)
% is approximately equal to the lower bound on capacity calculated using
% the PMF-based method when the SNR is high.  
%
% In this comparison, we treat the PMF-based answer as the correct capacity
% when we constrain x to be uniformly distributed. 

% To ensure full column rank for this randomly-generated channel matrix,
% have n_r>=n_t.  
% Although it's possible that the randomly-generated channel matrix would
% have linearly-dependent columns, it is a nearly-impossible probability.  
n_r = randi(10);
n_t = randi(n_r);

G = rand(n_r, n_t);

sigma_w = rand();

%% Choose xmax so that the SNR is very high

weakest_gain = min(G(:));
% This way, the ratio of the weakest peak signal value to sigma_w would be
% at least 100.  
x_max = 100 * sigma_w / weakest_gain / rand();

%% Compute the UnifX lower bound on capacity using the method under test
% Here, the full-rank method is used because G has full column rank and
% would be full-rank after the Q'-transform.
[lb_nats_conservative, lb_nats] = ...
    MIMOCapacity.computeCapacityLBUnifX(G, x_max, sigma_w);

tc.verifyEqual(lb_nats_conservative, lb_nats, ...
    ['The full-rank method for the UnifX LB should not calculate a ' ...
    'separate/different lb_nats_conservative.']);

%% Compute the UnifX lower bound using the PMF method

G = MIMOCapacity.simplifyChannelMatrix(G);

% Handle the special case where G is 0.
if(isequal(G,0))
    lb_nats_conservative = 0;
    lb_nats = 0;
    return
end

% Apply the Q' transform
Q = MIMOCapacity.computeQTTransform(G);
G = Q' * G;

variance_noise_out = repmat(sigma_w^2, size(G,1), 1);
max_nbins = 1e9/8/5;
min_trials = 1e9;

[lb_nats_expected, lb_nats_conservative_expected] = ...
    MIMOCapacityLBUnifX.calculateCapacityLBForUniformInput( ...
    G, x_max, variance_noise_out, max_nbins, min_trials);

%% Compare results

% Here, we must consider the differences between lb_nats and
% lb_nats_conservative as computed by
% MIMOCapacityLBUnifX.calculateCapacityLBForUniformInput.  
% lb_nats is likely a closer estimate of the true lower bound on capacity
% (given the uniformly-distributed x constraint), but has a possibility of
% over-estimating the lower bound due to its assumption that the
% probability is uniformly-distributed within each bin.  
% lb_nats_conservative, on the other hand, is designed to avoid
% overestimating the lower bound, but likely is too conservative, severely
% underestimating the capacity (especially when the generated PMF does not
% have enough trials to be smooth).  

tc.verifyEqual(lb_nats, lb_nats_expected, 'RelTol', 0.02, ...
    'lb_nats should be approximately equal lb_nats_expected.');
tc.verifyGreaterThanOrEqual(lb_nats, lb_nats_conservative_expected, ...
    'lb_nats should not be less than lb_nats_conservative_expected.');

end
