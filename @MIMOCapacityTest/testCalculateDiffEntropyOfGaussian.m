function testCalculateDiffEntropyOfGaussian(tc)
% testCalculateDiffEntropyOfGaussian checks that
% calculateDiffEntropyOfGaussian returns the expected result for a random
% number of independent Gaussian random variables, each with random
% variance.  

% number of dimensions
d = randi(5);

sigma = rand(d, 1) ./ rand(d, 1);

de_test = MIMOCapacity.calculateDiffEntropyOfGaussian(sigma.^2);

% check using ln(sigma * sqrt(2*pi*e)) as the differential entropy of a
% single Gaussian random variable.  

expected_de_each = log(sigma .* sqrt(2*pi*exp(1)));
% Total differential entropy should be the sum of the differential entropy
% for each Gaussian random variable (assuming that they're independent).
expected_de = sum(expected_de_each);

tc.verifyEqual(de_test, expected_de, ...
    ['The differential entropy from calculateDiffEntropyOfGaussian ' ...
    'should match the expected differential entropy.']);

end

