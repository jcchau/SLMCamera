function testGenerateReceivedPmfForUniformInput1DNoNoise(tc)
% testGenerateReceivedPmfForUniformInput1DNoNoise verifies that
% MIMOCapacityLBUnifX.generateReceivedPmfForUniformInput works even when
% variance_noise_out = 0.

G = rand()./rand();
x_max = rand()./rand();
variance_noise_out = 0;
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

[pmf, trials, y_min, y_max, hits] = ...
    MIMOCapacityLBUnifX.generateReceivedPmfForUniformInput( ...
    G, x_max, variance_noise_out, bins_per_dimension, ...
    min_trials, trials_per_batch);

%% preliminary checks

tc.verifyEqual(y_min, 0, ...
    'y_min should be 0.');
Gxmax = G*x_max;
tc.verifyEqual(y_max, Gxmax, ...
    sprintf('y_max should be G*x_max = %f.', Gxmax));

tc.verifyEqual(pmf, hits./trials, 'pmf should be the rate of hits.');

tc.verifyLength(pmf, bins_per_dimension, ...
    'pmf should have length bins_per_dimension.');

tc.verifyGreaterThanOrEqual(trials, min_trials, ...
    'trials should be at least min_trials.');

tc.verifyEqual(sum(hits), trials, ...
    'Without noise, all of the trials should land in a bin.');

%% the expected pmf

bin_boundaries = linspace(y_min, y_max, bins_per_dimension+1);

delta = (y_max - y_min) / bins_per_dimension;
peak_pmf = delta / Gxmax;

expected_pmf = zeros(bins_per_dimension, 1);

for ii = 1:bins_per_dimension
    if(bin_boundaries(ii+1) < 0 || bin_boundaries(ii) > Gxmax)
        % expected_pmf should be (and already is) 0.
    elseif(bin_boundaries(ii) >= 0 && bin_boundaries(ii+1) <= Gxmax)
        expected_pmf = peak_pmf;
    elseif(bin_boundaries(ii) < 0 && bin_boundaries(ii+1) > 0)
        % the right edge of the bin is > 0 and that part has a pdf equal
        % Gxmax.  
        expected_pmf = bin_boundaries(ii+1) / Gxmax;
    elseif(bin_boundaries(ii) < Gxmax && bin_boundaries(ii+1) > Gxmax)
        % the left part of the bin intersects the part of the pdf that
        % equals Gxmax.
        expected_pmf = (Gxmax - bin_boundaries(ii)) / Gxmax;
    else
        error('The above cases should have covered every possibility.');
    end % if-elseif-else ladder
end % for ii

%% and the expected standard deviation of the pmf
% from equation 7 of p.196 of lab book 3

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

