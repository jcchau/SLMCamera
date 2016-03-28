function nbins = fillMaxNBins(max_nbins, dimensions)
% fillMaxNBins returns the maximum number of bins along each dimension
% given a maximum number of bins while trying to avoid having a
% disproportionate number of bins in each dimension.  
%
%   NBINS = fillMaxNBins(MAX_NBINS, DIMENSIONS)
%
% NBINS is a row vector of the number of bins along each dimension.
%
% MAX_NBINS is the maximum total number of bins.
% DIMENSIONS is the number of dimensions.
%
% Essentially starts with nbins = max_nbins^(1/dimensions) along each
% dimension (a (hyper-)cube), then tries to incrementally expand the nbins
% in some dimensions without decreasing nbins in any dimension and without
% exceeding max_nbins.  

nbins = zeros(1, dimensions);

for d = dimensions:-1:1
    
    nbins(d) = floor(max_nbins^(1/d));
    
    max_nbins = max_nbins / nbins(d);
    
end % for d

end

