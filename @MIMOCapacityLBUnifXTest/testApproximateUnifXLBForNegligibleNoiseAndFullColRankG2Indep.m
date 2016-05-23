function testApproximateUnifXLBForNegligibleNoiseAndFullColRankG2Indep(tc)
% testApproximateUnifXLBForNegligibleNoiseAndFullColRankG2Indep tests
% approximateUnifXLBForNegligibleNoiseAndFullColRankG using 2 independent
% channels.

n_t = 2;
g = rand(n_t,1) + 0.1;
G = diag(g);
sigma_w = rand(); % sigma_w <= 1
x_max = 1e3 ./ rand(); % x_max >= 1e3

[lb_approx, Q, h_x, h_QTransposeGx, h_QTransposew] = ...
    MIMOCapacityLBUnifX. ...
    approximateUnifXLBForNegligibleNoiseAndFullColRankG( ...
    G, x_max, sigma_w);

% Check Q
tc.verifyEqual(Q, eye(n_t), ...
    'Q should be an identity matrix for these independent channel.');

% Check that h_x is the sum of the differential entropy of each x
expected_h_x = 2*log(x_max);
tc.verifyEqual(h_x, expected_h_x, 'AbsTol', 1e-15, ...
    'h_x should be the sum of the differential entropy of each x.');

% Check that h_QTransposeGx is the sum of the differential entropy of each
% G*x.
expected_h_Gx = sum(log(g.*x_max));
tc.verifyEqual(h_QTransposeGx, expected_h_Gx, 'AbsTol', 1e-14, ...
    ['h_QTransposeGx should equal the sum of the ' ...
    'differential entropy of each G*x.']);

% Check h_QTransposew
expected_h_w = ...
    MIMOCapacity.calculateDiffEntropyOfGaussian(repmat(sigma_w^2, n_t, 1));
tc.verifyEqual(h_QTransposew, expected_h_w, 'AbsTol', 1e-14, ...
    ['h_QTransposew should equal the sum of the differential entropy ' ...
    'of each w.']);

% check the computed lower bound
expected_lb = expected_h_Gx - expected_h_w;
tc.verifyEqual(lb_approx, expected_lb, 'AbsTol', 1e-12, ...
    'lb_approx should equal h_Gx - h_w.');

%% Check against the PMF-based result

variance_noise_out = repmat(sigma_w^2, n_t, 1);
max_nbins = 1e5; % 316 bins per dimension
min_trials = 1e8;

[lb_pmf, ~, ~, ~, ~, h_y, ~] = ...
    MIMOCapacityLBUnifX.calculateCapacityLBForUniformInput( ...
    G, x_max, variance_noise_out, max_nbins, min_trials);

% h(y) should equal h(z) since we're not removing any dimensions (Q=1).
% Hence, h(y) should approximately equal h(G*x).
tc.verifyEqual(h_QTransposeGx, h_y, 'RelTol', 1e-3, ...
    'h_QTransposeGx should be approximately equal h_y.');

tc.verifyEqual(lb_approx, lb_pmf, 'RelTol', 1e-3, ...
    ['The lower bound on capacity from both methods should be ' ... 
    'approximately equal.']);

end

