function de = calculateDiffEntropyOfGaussian(variance)
% calculateDiffEntropyOfGaussian calculates the differential entropy of a
% (multivariate) Gaussian random variable.  
%
% If more than one variance is provided (through parameter VARIANCE), the
% Gaussian random variables are assumed to be independent. 
%
%   DE = calculateDifEntropyOfGaussian(VARIANCE)
%
% DE is the differential entropy in nats.  
% VARIANCE is a vector containing the variance of each Gaussian random
%   variable.  

if(~isvector(variance))
    error('Parameter VARIANCE should be a vector.');
end

% From p.100 of lab book #3, the differential entropy of multiple
% independent Gaussian random variables is:
% 0.5 * ln( (2*pi*e)^{n_r} * prod_{i=1}^{n_r}(sigma^2) )
de = 0.5 * sum(log((2*pi .* variance)) + 1);

end

