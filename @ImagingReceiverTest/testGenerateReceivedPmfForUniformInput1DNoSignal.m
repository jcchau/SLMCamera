function testGenerateReceivedPmfForUniformInput1DNoSignal(tc)
% testGenerateReceivedPmfForUniformInput1DNoSignal verifies that
% generateReceivedPmfForUniformInput works when x_max = 0.

G = rand()./rand();
x_max = 0;
variance_noise_out = rand()./rand();
bins_per_dimension = randi(100);

min_trials = poissrnd(1e6 * log(bins_per_dimension)^2) + 1; ...
    % lambda < 2.13e7; min_trials >= 1.

trials_per_batch = randi(2^randi(15)); % random int from 1 to 32768

% Prevent this from taking too many batches (since that would take too
% long).
max_num_batches = 1e5;
if(min_trials/trials_per_batch > max_num_batches)
    min_trials = randi(max_num_batches * trials_per_batch);
end

%% run the method under test

[pmf, trials, min_noise_to_delta_ratio, y_min, y_max, hits] = ...
    ImagingReceiver.generateReceivedPmfForUniformInput( ...
    G, x_max, variance_noise_out, bins_per_dimension, ...
    min_trials, trials_per_batch);

%% preliminary checks

tc.verifyEqual(y_min, -5 * sqrt(variance_noise_out), ...
    'y_min should be -5 times the noise standard deviation.');
tc.verifyEqual(y_max, 5 * sqrt(variance_noise_out), ...
    'y_max should be 5 times the noise standard deviation.');

tc.verifyEqual(pmf, hits./trials, 'pmf should be the rate of hits.');

tc.verifyLength(pmf, bins_per_dimension, ...
    'pmf should have length bins_per_dimension.');

tc.verifyGreaterThanOrEqual(trials, min_trials, ...
    'trials should be at least min_trials.');

tc.verifyLessThanOrEqual(sum(hits), trials, ...
    'There should not be more hits than trials.');

missed_hits = trials - sum(hits);
expected_missed_hits = 2 * normcdf(-5, 0, 1) * trials;
tc.verifyLessThanOrEqual(missed_hits, ceil(10 * expected_missed_hits), ...
    sprintf(['Too many misses: %d hits are expected to miss the bins ' ...
    'but %d actually missed.'], ...
    expected_missed_hits, missed_hits));

%% check min_noise_to_delta_ratio
% expected value calculated according to p.191 of lab book 3.

delta = (y_max - y_min) / bins_per_dimension;
expected_min_noise_to_delta_ratio = sqrt(variance_noise_out)/delta;
tc.verifyEqual(min_noise_to_delta_ratio, ...
    expected_min_noise_to_delta_ratio, ...
    sprintf('min_noise_to_delta_ratio should be equal %f.', ...
    expected_min_noise_to_delta_ratio));

%% compute the expected pmf
% the expected pmf should just be that of a Gaussian pdf representing the
% noise.

sigma = sqrt(variance_noise_out);
bin_boundaries = linspace(y_min, y_max, bins_per_dimension+1);

expected_pmf = normcdf(bin_boundaries(2:end), 0, sigma) - ...
    normcdf(bin_boundaries(1:end-1), 0, sigma);

%% compute the expected standard deviation of resulting pmf of each bin
% equation 7 from p.196 in lab book 3

stddev_pmf = sqrt(expected_pmf .* (1-expected_pmf) ./ trials);

%% check that the pmf is within 4 standard deviations of expected_pmf

min_pmf = expected_pmf - 4*stddev_pmf;
min_pmf(min_pmf<0) = 0;

max_pmf = expected_pmf + 4*stddev_pmf;

tc.verifyLessThanOrEqual(pmf, max_pmf(:), ...
    'pmf should not be greater than max_pmf.');
tc.verifyGreaterThanOrEqual(pmf, min_pmf(:), ...
    'pmf should not be less than min_pmf.');

end

