function testGenerateReceivedPmfForUniformInput2DCorrelatedGNoNoise(tc)
% testGenerateReceivedPmfForUniformInput2DCorrelatedGNoNoise verifies that
% the returned output is correct even when there is correlation.

G = [ 1 1; 1 1 ];
x_max = 2;
variance_noise_out = [0; 0];
bins_per_dimension = 4;
min_trials = 1e6 * log(4^2)^2; % 7.7e6

%% run the method under test

[pmf, trials, y_min, y_max, hits] = ...
    MIMOCapacityOld.generateReceivedPmfForUniformInput( ...
    G, x_max, variance_noise_out, bins_per_dimension, min_trials);

%% preliminary checks

% verify that y_min and y_max are 5 noise standard deviations out.
tc.verifyEqual(y_min, [0; 0], 'y_min should be a column vector of 0.');
tc.verifyEqual(y_max, [4; 4], 'y_max should be a column vector of 4.');

% check that pmf = hits/trial
tc.verifyEqual(pmf, hits/trials, 'pmf should be the rate of hits');

tc.verifyEqual(size(pmf), [bins_per_dimension, bins_per_dimension], ...
    'Size of pmf does not match bins_per_dimension.');

tc.verifyGreaterThanOrEqual(trials, min_trials, ...
    'trials should be at least min_trials.');

total_hits = sum(sum(hits));
tc.verifyEqual(total_hits, trials, ...
    ['With no noise, the total hits in all bins should equal the ' ...
    'number of trials.']);

%% compute the expected pmf
% Non-diagonal elements are 0 because the two receivers always receive the
% same value (no noise).
% The pdf of each single receiver is the convolution of the pdfs for two
% uniformly distributed RVs from 0 to 2.

expected_pmf = [ ...
    1/8, 0, 0, 0; ...
    0, 3/8, 0, 0; ...
    0, 0, 3/8, 0; ...
    0, 0, 0, 1/8];

%% compute the expected standard deviation

stddev_pmf = sqrt(expected_pmf .* (1-expected_pmf) ./ trials);

%% check that the pmf from the method under test is within 4.*stddev_pmf.

min_pmf = expected_pmf - 4.*stddev_pmf;
min_pmf(min_pmf<0) = 0;
max_pmf = expected_pmf + 4.*stddev_pmf;

tc.verifyGreaterThanOrEqual(pmf, min_pmf, ...
    'pmf should be at least min_pmf.');
tc.verifyLessThanOrEqual(pmf, max_pmf, ...
    'pmf should not be greater than max_pmf.');

end

