function testGenerateReceivedPmfForUniformInput1D(tc)
% testGenerateReceivedPmfForUniformInput1D verifies that
% generateReceivedPmfForUniformInput returns (approximately) the correct
% result for a SISO system (one input and one output).  
%
% Note that many of the expected values were calculated by hand and are
% hard-coded, so don't change the initial parameters (otherwise, these
% hard-coded expected results will need to be recalculated).

G = 0.5;
x_max = 6;
variance_noise_out = 4;
bins_per_dimension = 23;
min_trials = 1e7;

%% run the method under test

[pmf, trials, y_min, y_max, hits] = ...
    MIMOCapacityOld.generateReceivedPmfForUniformInput( ...
    G, x_max, variance_noise_out, bins_per_dimension, min_trials);

%% preliminary checks

% verify that y_min and y_max are 5 noise standard deviations out.
tc.verifyEqual(y_min, -10, 'y_min should be -10.');
tc.verifyEqual(y_max, 13, 'y_max should be 13.');

tc.verifyEqual(pmf, hits/trials, 'pmf should be the rate of hits');

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

%% compute the expected pmf

% the boundaries around each bin
y = -10:13;
% the pdf for each value of y (calculated via the flip and slide method of
% convolution), using the uniform pdf as a window of height 1/3 and length
% 3.
pdf_y = (normcdf(y, 0, 2) - normcdf(y-3, 0, 2)) ./ 3;

% the pmf should equal delta (1) times the mean pdf in that bin.  
% the mean pdf is between the minimum pdf and the maximum pdf.
% these are row vectors.
min_pmf_y = min(pdf_y(1:end-1), pdf_y(2:end));
max_pmf_y = max(pdf_y(1:end-1), pdf_y(2:end));

% there is a peak in the pdf of y in the 12th bin (at 1.5, between y=1 and
% y=2) that's not captured above.  
max_pmf_y(12) = (normcdf(1.5, 0, 2) - normcdf(-1.5, 0, 2)) ./ 3;

%% calculate the expected standard deviation of the returned pmf

% with these test parameters, all of the pmfs are < 0.5.
% using equation 7 from p.196 of lab book 3, the standard deviation of the
% error in the pmf equals sqrt(p_i.*(1-p_i)./n), where p_i is the pmf of
% bin i.  
% We choose max_pmf_y since for pmf < 0.5, a larger pmf yields a larger
% expected standard deviation for the simulated pmf.  

max_stddev_pmf = sqrt(max_pmf_y .* (1-max_pmf_y) ./ trials);

%% verify that the calculated pmf is not too far beyond what we expect
% That the calculated pmf from the method under test is not more than
% 4*max_stddev_pmf beyond the bounds of min_pmf_y and max_pmf_y.

% expand the bounds min_pmf_y and max_pmf_y by 4*max_stddev_pmf
min_pmf_y = min_pmf_y - 4*max_stddev_pmf;
min_pmf_y(min_pmf_y<0) = 0;

max_pmf_y = max_pmf_y + 4*max_stddev_pmf;

tc.verifyLessThanOrEqual(pmf, max_pmf_y(:), ...
    'pmf should not exceed max_pmf_y.');
tc.verifyGreaterThanOrEqual(pmf, min_pmf_y(:), ...
    'pmf should not be less than min_pmf_y.');

% ycenters = -9.5:1:12.5;
% plot(ycenters, pmf, 'bo-');
% hold on
% plot(ycenters, max_pmf_y, 'r-')
% plot(ycenters, min_pmf_y, 'g-')
% grid minor

end

