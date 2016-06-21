function testApproximateUnifXLBForNegligibleNoiseAndFullColRankG3x2(tc)
% testApproximateUnifXLBForNegligibleNoiseAndFullColRankG3x2 test
% approximateUnifXLBForNegligibleNoiseAndFullColRankG using a random 3x2 G.
%
% For comparison, we use calculateCapacityLBForUniformInput both with and
% without transform Q.

n_r = 3;
n_t = 2;
G = rand(n_r, n_t);
sigma_w = rand();
x_max = 1e3 ./ rand();

% method under test
[lb_approx, Q, h_x, h_QTransposeGx, h_QTransposew] = ...
    MIMOCapacityLBUnifX. ...
    approximateUnifXLBForNegligibleNoiseAndFullColRankG( ...
    G, x_max, sigma_w);

%% check with calculateCapacityLBForUniformInput using transformed channel

A = Q'*G;
variance_noise_out_QT = repmat(sigma_w^2, n_t, 1);
max_nbins_QT = 1e5; % 316 bins per dimension
min_trials_QT = 1e8;

[lb_pmf_QT, ~, ~, ~, ~, h_QTy, ~] = ...
    MIMOCapacityLBUnifX.calculateCapacityLBForUniformInput( ...
    A, x_max, variance_noise_out_QT, max_nbins_QT, min_trials_QT);

% h(Q'*y) should equal h(z)
tc.verifyEqual(h_QTransposeGx, h_QTy, 'RelTol', 1e-3, ...
    ['h_QTransposeGx should be approximately equal h(Q''*y) ' ...
    'computed using the PMF.']);

% And the lower bound on capacity should be approximately the same.
tc.verifyEqual(lb_approx, lb_pmf_QT, 'RelTol', 1e-3, ...
    ['The I(x;z) should approximately equal I(x;Q''*y) ' ...
    'computed using the PMF.']);

%% check with calculateCapacityLBForUniformInput without transformed ch.

variance_noise_out = repmat(sigma_w^2, n_r, 1);
max_nbins = 125e6; % 500 bins per dimension
min_trials = 1e3 * max_nbins;

[lb_pmf, ~, ~, ~, ~, h_y, ~] = ...
    MIMOCapacityLBUnifX.calculateCapacityLBForUniformInput( ...
    G, x_max, variance_noise_out, max_nbins, min_trials);

% h(y) should differ from h(Q'*y), so there is no point in comparing h(y)
% against h(Q'*G*x) since we don't expect them to be approximately equal.

% The lower bound on capacity should be approximately the same.
tc.verifyEqual(lb_approx, lb_pmf, 'RelTol', 1e-3, ...
    'I(x;z) should approximately equal I(x;y) computed using the PMF.');

end

