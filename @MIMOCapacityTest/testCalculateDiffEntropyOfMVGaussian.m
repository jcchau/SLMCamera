function testCalculateDiffEntropyOfMVGaussian(tc)
% testCalculateDiffEntropyOfMVGaussian tests
% MIMOCapacity.calaculateDiffEntropyOfMVGaussian against Theorem 8.6.5 in
% Cover2005.

n_r = randi(10);
n_t = randi(10);

G = rand(n_r, n_t);
K = G * G';

nats_expected = 0.5 * log((2*pi*exp(1))^n_r * det(K));

nats_test = MIMOCapacity.calculateDiffEntropyOfMVGaussian(K);

tc.verifyEqual(nats_test, nats_expected, ...
    'Calculated differential entropy does not match.');

end

