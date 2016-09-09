function testComputeCapacityLBUnifXFullRankLowSNR(tc)
% testComputeCapacityLBUnifXFullRankRankLowSNR uses a randomly-generated
% full column rank channel matrix to verify that the lower-bound on
% capacity calculated via the full column rank branch (using
% MIMOCapacityLBUnifX.approximateUnifXLBForNegligibleNoiseAndFullColRankG)
% is approximately equal to the lower bound on capacity calculated using
% the PMF-based method when the SNR is *LOW*.  
%
% In this comparison, we treat the PMF-based answer as the correct capacity
% when we constrain x to be uniformly distributed. 

% To ensure full column rank for this randomly-generated channel matrix,
% have n_r>=n_t.  
% Although it's possible that the randomly-generated channel matrix would
% have linearly-dependent columns, it is a nearly-impossible probability.  
% Don't want n_r to be too big to keep the PMF-based approach accurate.
n_r = randi(5);
n_t = randi(n_r);

G = rand(n_r, n_t);

sigma_w = rand();
% Choose x_max so that SNR is low enough that the noise is NOT negligible.
avg_rx_gain = mean(sum(G,2)); % likely > 1 with enough columns
x_max = sigma_w / avg_rx_gain * 3*rand();

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

[lb_nats_expected, ~] = ...
    MIMOCapacityLBUnifX.calculateCapacityLBForUniformInput( ...
    G, x_max, variance_noise_out, max_nbins, min_trials);

%% Compare results

% Since noise is not negligible, we expect the calculated lb_nats
% calculated by the method under test to be less than lb_nats_expected.

tc.verifyLessThan(lb_nats, lb_nats_expected, ...
    ['lb_nats should be not be significantly greater than ' ...
    'lb_nats_expected.']);

end
