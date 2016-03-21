function de = calculateDiffEntropyFromPmf(pmf, bin_size)
% calculateDiffEntropyFromPmf computes the differential entropy of a random
% variable given its probability mass function (PMF) and the size of the
% bin used to compute the PMF.  
%
%   [DE, DE_BITS] = calculateDiffEntropyFromPmf(PMF, BIN_SIZE)
%
% DE (scalar) is the differential entropy of the random variable in nats.  
%
% PMF (vector or multi-dimension matrix) is the probability mass function
%   of the random variable.  PMF is unrolled in this method (as PMF(:)), so
%   the specific shape of the matrix does not matter here. 
% BIN_SIZE (vector) if BIN_SIZE is NOT a scalar, each element will be
%   treated as the "delta" for a dimension of the bin and the bin is
%   assumed to be rectangular.
% BIN_SIZE (scalar) is the size (multi-dimension volume) of each bin of
%   random variable values used to calculate the PMF (so that each PMF
%   value is the probability that the random variable has the range of
%   values of a single bin).  All PMF bins are assumed to be the same size.
%   BIN_SIZE can be computed as the product of the "delta" for each
%   dimension of the bin (for rectangular bins).  
%
% This method calculates differential entropy using by first computing the
% discrete entropy of the random variable described by the PMF, and then
% applying theorem 8.3.1 of Cover2005 to convert the discrete entropy into
% differential entropy.  This conversion is described in p.121 of lab book
% #3 (Imaging Receivers & Photodetector Arrays).  

if(isvector(bin_size))
    % Also covers the case where bin_size is a scalar.
    log_delta = sum(log(bin_size));
else
    error('BIN_SIZE must be a scalar or a vector.');
end

% When the pmf is zero, that bin's contribution to the entropy is zero, so
% only use pmf(pmf>0).  
% Also unroll pmf.
pmf = pmf(pmf>0);

% Compute the discrete entropy (in nats).
entropy_discrete = -sum(pmf .* log(pmf));

% compute the differential entropy (applying theorem 8.3.1 from Cover2005).
de = entropy_discrete + log_delta;

end

