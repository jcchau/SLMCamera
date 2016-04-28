function testCalculateCapacityForUniformInputVarianceH(tc)
% testCalculateCapacityForUniformInputVarianceH checks that the h_y
% in calculateCapacityForUniformInput is within 5 standard deviations (as
% specified by variance_H) of the expected result.  

mu = 0;
variance = 0.1;
sigma = sqrt(variance);
a = -5 * sigma;
b = 5 * sigma;
nbins = 1e9/8/5;
trials = 1e9;

[nats, min_nats, variance_H, nbins, trials, h_y, pmf] = ...
    MIMOCapacityOld.calculateCapacityForUniformInput(0, 1, variance, ...
    nbins, trials);

[h_y_from_true_pmf, true_pmf] = ...
    MIMOCapacity.calculateDiffEntropyOfClippedNormal( ...
    a, b, nbins, mu, sigma);

h_noclip = MIMOCapacity.calculateDiffEntropyOfGaussian(variance);

stddev_h = sqrt(variance_H);

tc.verifyEqual(h_y, h_y_from_true_pmf, 'AbsTol', 5*stddev_h, ...
    ['The calculated h_y should most likely be within 5 of the ' ...
    'calculated standard deviations of h_y calculated from the ' ...
    'true PMF.']);

y = linspace(a, b, nbins);
plot(y, pmf, 'b');
hold on
plot(y, true_pmf, 'r');
legend('pmf', 'true pmf');
xlabel('y');

end

