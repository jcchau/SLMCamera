function testGenerateReceivedPmfForUniformInput3D(tc)
% testGenerateReceivedPmfForUniformInput3D verifies that
% generateReceivedPmfForUniformInput works for 3 dimensions.  
% If it works for 3D, probably safe to assume that it also works for even
% more receivers.  

G = [ 1 0 0; 0 1 0; 0 0 1 ];
x_max = 2;
variance_noise_out = [ 1; 1; 1 ];
bins_per_dimension = 12; % each bin is unit big
min_trials = 1e6 * log(12^3)^2; % 5.6e7
trials_per_batch = 45*2^15;

%% run the method under test

[pmf, trials, y_min, y_max, hits] = ...
    MIMOCapacityOld.generateReceivedPmfForUniformInput( ...
    G, x_max, variance_noise_out, bins_per_dimension, ...
    min_trials, trials_per_batch);

%% preliminary checks

tc.verifyEqual(y_min, [-5; -5; -5], ...
    'y_min should be -5 for every receiver.');
tc.verifyEqual(y_max, [7; 7; 7], ...
    'y_max should be 7 for every receiver.');

tc.verifyEqual(pmf, hits/trials, 'pmf should be the rate of hits.');

tc.verifyEqual(size(pmf), repmat(bins_per_dimension, 1, 3), ...
    'Size of pmf does not match bins_per_dimension');

tc.verifyGreaterThanOrEqual(trials, min_trials, ...
    'trials should be at least min_trials.');

total_hits = sum(hits(:));

tc.verifyLessThanOrEqual(total_hits, trials, ...
    'There should not be more hits than trials.');

missed_hits = trials - total_hits;
expected_hit_rate = (1 - 2*normcdf(-5, 0, 1))^3;
expected_missed_hits = (1 - expected_hit_rate) * trials;

% Whether a trial lands in any of the bins is a Bernoulli RV (with variance
% p*(1-p).  The variance of trials trials is trials*p*(1-p).  
stddev_missed_hits = sqrt(...
    trials * expected_hit_rate * (1-expected_hit_rate));
% Accept up to 4 standard deviations additional misses.
max_missed_hits = expected_missed_hits + 4*stddev_missed_hits;

tc.verifyLessThanOrEqual(missed_hits, max_missed_hits, ...
    sprintf(['Too many misses: %d hits are expected to miss the bins ' ...
    'but %d actually missed.'], ...
    expected_missed_hits, missed_hits));

%% compute the expected PMF

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

% Since both dimensions of the pdf are independent, the 3d pdf is just the
% product of the corresponding 1d pdf for each dimension.
% This min and max pdf is also the min and max pmf because each bin is unit
% size.  
min_pmf = zeros(repmat(bins_per_dimension, 1, 3));
max_pmf = zeros(repmat(bins_per_dimension, 1, 3));
for ii = 1:bins_per_dimension
    for jj = 1: bins_per_dimension
        for kk = 1:bins_per_dimension
            min_pmf(ii,jj,kk) = min_pdf_y_1d(ii) * min_pdf_y_1d(jj) * ...
                min_pdf_y_1d(kk);
            max_pmf(ii,jj,kk) = max_pdf_y_1d(ii) * max_pdf_y_1d(jj) * ...
                max_pdf_y_1d(kk);
        end
    end
end

%% calculate the expected standard deviation of the returned pmf

% Assume that all of the pmfs are < 0.5, so we can assume that the standard
% deviation calculated using max_pmf would be the maximum pf
if(any(max_pmf > 0.5))
    error('We assume that all of the pmf values should be < 0.5.');
end

stddev_pmf = sqrt(max_pmf.*(1-max_pmf) ./ trials);

%% verify that the calculated pmf is not too far beyond what we expect
% Accept up to 5 standard deviations beyond the min_pmf and max_pmf
% calculated above.
% Use 5 standard deviations because we did encounter a case that just
% seemed to be bad luck with 4 standard deviations.

min_pmf = min_pmf - 5 .* stddev_pmf;
min_pmf(min_pmf<0) = 0;

max_pmf = max_pmf + 5 .* stddev_pmf;

tc.verifyGreaterThanOrEqual(pmf, min_pmf, ...
    'pmf should be at least min_pmf.');
tc.verifyLessThanOrEqual(pmf, max_pmf, ...
    'pmf should not be greater than max_pmf.');

end

