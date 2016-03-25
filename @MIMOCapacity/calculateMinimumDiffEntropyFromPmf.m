function de = calculateMinimumDiffEntropyFromPmf(pmf, bin_size)
% calculateMinimumDiffEntropyFromPmf computes the minimum differential
% entropy of a random variable given its probability mass function (PMF)
% and the size of the bin used to compute the PMF.  
%
% This method assumes that the PMF is centered so that the peak probability
% density function (PDF) occurs at the values of the random variable that
% correspond to the center bin(s) of the matrix PMF, that the PDF descends
% monotonically from this peak, and that the PDF is concave (negative
% second derivative) (or flat) within the bins around this peak.  
%
% This method differs from calculateDiffEntropyFromPmf because that method
% assumes that the PDF is uniform within each bin (which yields the maximum
% differential entropy).  
%
%   DE = calculateMinimumDiffEntropyFromPmf(PMF, BIN_SIZE)
%
% DE is the minimum differential entropy of the random variable in nats.  
%
% PMF (vector or multi-dimension matrix) is the probability mass function
%   of the random variable.  Each dimension of the PMF matrix should
%   correspond to a dimension of the random variable.  (PMF is NOT
%   unrolled, so the shape of PMF is important.)
% BIN_SIZE is a column vector containing the size of each bin used for the
%   PMF.  

%% input validation

if(~iscolumn(bin_size))
    error('Parameter BIN_SIZE must be a column vector.');
end

num_dimensions = length(bin_size);
if(~( num_dimensions==1 && isvector(pmf) || ...
        (ndims(pmf) == num_dimensions) ))
    error(['The implied dimension of the random variable don''t match ' ...
        'in parameters PMF and BIN_SIZE.']);
end

% Check if pmf is empty because this method assumes that each dimension has
% at least size 1.  
if(numel(pmf) == 0)
    % pmf is empty
    error('Parameter PMF is empty.');
end

%% average and minimum PDF

bin_volume = prod(bin_size);

% average pdf
pavg = pmf ./ bin_volume;

% min pdf
pmin = MIMOCapacity.minimumNeighbor(pavg, num_dimensions);

%% maximum pdf
% For the bins at (or immediately surrounding) the center of the
% distribution, the max possible PDF value is not the max pavg of the
% neighboring bins, but larger.  
% Method maximumNeighbor does not account for this, so pmax needs to be
% adjusted for these central bins.  

pmax = MIMOCapacity.maximumNeighbor(pavg);

% Determine which bins need to be adjusted to have the peak PMF
index_central_bins = cell(1, num_dimensions);

for d = 1:num_dimensions
    num_bins_along_d = size(pmf, d);
    
    if(mod(num_bins_along_d, 2) == 0)
        % this dimension has an even number of bins
        index_central_bins(d) = {[num_bins_along_d/2, ...
            num_bins_along_d/2+1]};
    else
        % this dimension has an odd number of bins
        index_central_bins(d) = {(num_bins_along_d+1)/2};
    end % if(mod(bin_size, 2) == 0)
end % for d = 1:num_dimensions

% From p.17 of lab book #4 (Imaging Receivers & Photodetector Arrays).
pmax(index_central_bins{:}) = 2 .* pavg(index_central_bins{:}) - ...
    pmin(index_central_bins{:});

%% use mean, min, and max PDF to compute the minimum entropy
% Algorithm from p. 15-17 of lab book #4.

alpha = (pavg - pmin) ./ (pmax - pmin);

% When (pmax-pmin) is 0, alpha is NaN.  In this case, the proportion of the
% bin that is pmax vs pmin doesn't matter because pmax == pmin, so just set
% alpha to 1 to avoid propagating NaN through to the calculated
% differential entropy.  
alpha(pmax==pmin) = 1;

% FREE MEMORY: remove pavg from memory since it's no longer needed
clear pavg

% Drop zero values from pmax and pmin since zero-probability events don't
% contribute to entropy.
% And also remove the corresponding entries from alpha so they still match.
% (This also unrolls the matrices.)
alpha_pmin = alpha(pmin>0);
alpha_pmax = alpha(pmax>0);
clear alpha % FREE MEMORY: alpha no longer needed

pmin = pmin(pmin>0);
pmax = pmax(pmax>0);

de = -bin_volume .* ( ...
    sum(alpha_pmax .* pmax .* log(pmax)) + ...
    sum((1-alpha_pmin) .* pmin .* log(pmin)) );

end
