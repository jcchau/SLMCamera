function testApproximateUnifXLBForNegligibleNoiseAndFullColRankGScalar( ...
    tc)
% testApproximateUnifXLBForNegligibleNoiseAndFullColRankGScalar tests
% approximateUnifXLBForNegligibleNoiseAndFullColRankG using a scalar G.

G = rand() + 0.1;
sigma_w = rand(); % sigma_w <= 1
x_max = 1e3 ./ rand(); % x_max >= 1e3

[lb_approx, Q, h_x, h_QTransposeGx, h_QTransposew] = ...
    MIMOCapacityLBUnifX. ...
    approximateUnifXLBForNegligibleNoiseAndFullColRankG( ...
    G, x_max, sigma_w);

% Check Q
tc.verifyEqual(Q, 1, ...
    'Q should be 1 for this scalar channel.');

% Check h_x equals ln(x_max)
tc.verifyEqual(h_x, log(x_max), 'AbsTol', 1e-15, ...
    'h_x should be log(x_max).');

% Check h_QTransposeGx is ln(G*x_max);
expected_h_QTransposeGx = log(G*x_max);
tc.verifyEqual(h_QTransposeGx, expected_h_QTransposeGx, ...
    'AbsTol', 1e-15, 'h_QTransposeGx should equal log(G*x_max).');

% Check h_QTransposew
expected_h_QTransposew = ...
    MIMOCapacity.calculateDiffEntropyOfGaussian(sigma_w^2);
tc.verifyEqual(h_QTransposew, expected_h_QTransposew, 'AbsTol', 1e-15, ...
    ['h_QTransposew should equal the differential entropy of a ' ...
    'Gaussian RV with variance sigma_w^2.']);

% check the computed lower bound
expected_lb = expected_h_QTransposeGx - expected_h_QTransposew;
tc.verifyEqual(lb_approx, expected_lb, 'AbsTol', 1e-12, ...
    'lb_approx should equal h_QTransposeGx - h_QTransposew.');

%% Check against the PMF-based result

variance_noise_out = sigma_w^2;
max_nbins = 1e5;
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

