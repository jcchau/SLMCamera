function nats = calculateDiffEntropyOfMVGaussian(K)
% calculateDiffEntropyOfMVGaussian calculates the differential entropy of a
% multivariate Gaussian random variable given covariance matrix K as
% specified in Cover2005, Theorem 8.6.5.
%
%   nats = calculateDiffEntropyOfMVGaussian(K)
%
% nats is the calculated differential entropy in nats.
%
% K is the covariance matrix of the Gaussian random variable.  

n = size(K, 1);

% nats = 0.5*log((2*pi*exp(1))^n * det(K));
% From lab book #4, p.87:
nats = 0.5*(log((2*pi)^n * det(K)) + n);

end

