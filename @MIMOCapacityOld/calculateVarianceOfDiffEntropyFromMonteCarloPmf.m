function variance = calculateVarianceOfDiffEntropyFromMonteCarloPmf( ...
    pmf, num_trials)
% calculateVarianceOfDiffEntropyFromMonteCarloPmf calculates/estimates the
% variance of the error in the differential entropy calculated using a PMF
% generated through a Monte Carlo simulation.  
%
%   VARIANCE = calculateVarianceOfDiffEntropyFromMonteCarloPmf( ...
%       PMF, NUM_TRIALS)
%
% VARIANCE (scalar) is the computed variance of the calculated differential
%   entropy of the random variable described by the PMF.
%
% PMF (matrix) is the probability mass function generated through Monte
%   Carlo simulations. 
% NUM_TRIALS (scalar) is the number of trials used to generate the PMF.  

% unroll pmf and drop any values that are zero.
% Justification: lab book #3, p. 183-184.  
pmf = pmf(pmf>0);

% Computed according to p.199 of lab book #3 (Imaging Receivers &
% Photodetector Arrays).
variance = sum( pmf.*(1-pmf) .* (-log(pmf)-1).^2 ./ num_trials );

end

