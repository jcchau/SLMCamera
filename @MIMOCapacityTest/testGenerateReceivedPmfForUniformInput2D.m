function testGenerateReceivedPmfForUniformInput2D(tc)
% testGenerateReceivedPmfForUniformInput2D 

G = [ 1, 0; 0, 1 ];
x_max = 2;
variance_noise_out = [1; 1];
bins_per_dimension = 12;
min_trials = 1e6 * log(12^2)^2; % 2.5e7

%% run the method under test

[pmf, trials, y_min, y_max, hits] = ...
    MIMOCapacity.generateReceivedPmfForUniformInput( ...
    G, x_max, variance_noise_out, bins_per_dimension, min_trials);

%% preliminary checks

% verify that y_min and y_max are 5 noise standard deviations out.
tc.verifyEqual(y_min, [-5; -5], 'y_min should be a column vector of -5.');
tc.verifyEqual(y_max, [7; 7], 'y_max should be a column vector of 7.');

% check that pmf = hits/trial
tc.verifyEqual(pmf, hits/trials, 'pmf should be the rate of hits');

tc.verifyEqual(size(pmf), [bins_per_dimension, bins_per_dimension], ...
    'Size of pmf does not match bins_per_dimension.');

tc.verifyGreaterThanOrEqual(trials, min_trials, ...
    'trials should be at least min_trials.');

total_hits = sum(sum(hits));

tc.verifyLessThanOrEqual(total_hits, trials, ...
    'There should not be more hits than trials.');

missed_hits = trials - total_hits;
expected_hit_rate = (1 - 2*normcdf(-5, 0, 1))^2;
expected_missed_hits = (1 - expected_hit_rate) * trials;
tc.verifyLessThanOrEqual(missed_hits, ceil(10 * expected_missed_hits), ...
    sprintf(['Too many misses: %d hits are expected to miss the bins ' ...
    'but %d actually missed.'], ...
    expected_missed_hits, missed_hits));

%% compute the expected PMF in one dimension

sigma_noise = sqrt(variance_noise_out);

bin_boundaries = -5:7;
pdf_y_1d = (normcdf(bin_boundaries, 0, sigma_noise(1)) - ...
    normcdf(bin_boundaries-x_max, 0, sigma_noise(1))) ./ x_max;

% min and max pdf values for each 1d pdf
min_pdf_y_1d = min(pdf_y_1d(1:end-1), pdf_y_1d(2:end));
max_pdf_y_1d = max(pdf_y_1d(1:end-1), pdf_y_1d(2:end));

% turn both into column vectors (so we know whether it's a column or row
% vector)
min_pdf_y_1d = min_pdf_y_1d(:);
max_pdf_y_1d = max_pdf_y_1d(:);

% Since both dimensions of the pdf are independent, the 2d pdf is just the
% product of the corresponding 1d pdf for each dimension.
% This min and max pdf is also the min and max pmf because each bin is unit
% size.  
min_pmf = min_pdf_y_1d * min_pdf_y_1d';
max_pmf = max_pdf_y_1d * max_pdf_y_1d';

%% calculate the expected standard deviation of the returned pmf

% Assume that all of the pmfs are < 0.5, so we can assume that the standard
% deviation calculated using max_pmf would be the maximum pf
if(any(max_pmf > 0.5))
    error('We assume that all of the pmf values should be < 0.5.');
end

stddev_pmf = sqrt(max_pmf.*(1-max_pmf) ./ trials);

%% verify that the calculated pmf is not too far beyond what we expect
% accept up to 4 standard deviations beyond the min_pmf and max_pmf
% calculated above.

min_pmf = min_pmf - 4 .* stddev_pmf;
min_pmf(min_pmf<0) = 0;

max_pmf = max_pmf + 4 .* stddev_pmf;

tc.verifyGreaterThanOrEqual(pmf, min_pmf, ...
    'pmf should be at least min_pmf.');
tc.verifyLessThanOrEqual(pmf, max_pmf, ...
    'pmf should not be greater than max_pmf.');

% y_centers = -4.5:6.5;
% surf(y_centers, y_centers, pmf);

end

