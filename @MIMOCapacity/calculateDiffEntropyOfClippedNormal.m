function [de, pmf] = calculateDiffEntropyOfClippedNormal( ...
    a, b, nbins, mu, sigma)
% calculateDiffEntropyOfClippedNormal calculates the differential entropy
% of a clipped normal distribution.
%
%   DE = calculateDiffEntropyOfClippedNormal(A, B, NBINS, MU, SIGMA)
%
% DE is the calculated differential entropy.
% PMF is the generated PMF used to calculate the differential entropy.
%   (Provided for debugging.)
%
% A is the lower cutoff of the clipped normal distribution.
% B is the upper cutoff ot the clipped normal distribution.
% NBINS is number of bins used to approximate the PMF.  
% MU is the mean of the normal distribution (before clipping).  Default is
%   0.
% SIGMA is the standard deviation of the normal distribution (before
%   clipping).  Default is 1. 
%
% Unlike what's called the "truncated normal distribution", we do not
% normalize the probabilities so that the total probability of the
% remaining range of values is 1.  

if(nargin < 5)
    sigma = 1;
end
if(nargin < 4)
    mu = 0;
end

if(a>b)
    error('Parameter A must be less than B.');
end

% bin boundaries
bb = linspace(a, b, nbins+1);

cdf_at_boundaries = normcdf(bb, mu, sigma);

pmf = cdf_at_boundaries(2:end) - cdf_at_boundaries(1:end-1);

delta = (b-a)/nbins;

de = MIMOCapacity.calculateDiffEntropyFromPmf(pmf, delta);

end

