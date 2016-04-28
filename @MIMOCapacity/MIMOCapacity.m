classdef MIMOCapacity
    %MIMOCAPACITY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        de = calculateDiffEntropyFromPmf(pmf, bin_size)
        
        li = convertToLinearIndex(weights, subs)
        weights = convertToLinearIndexWeights(mat_size)
        msi = convertLinearToSubscriptIndex(weights, li)
        
        out = minimumNeighbor(in, num_dimensions)
        out = maximumNeighbor(in)
        
        de = calculateMinimumDiffEntropyFromPmf(pmf, bin_size)
        
        de = calculateDiffEntropyOfGaussian(variance)
        
        nbins = fillMaxNBins(max_nbins, dimensions)
        
        [de, pmf] = calculateDiffEntropyOfClippedNormal( ...
            a, b, nbins, mu, sigma)
        
        out = removeZeroRowsAndCols(in)
        
        msi = convertPointToSubscriptIndex(y, ymin, delta, nbins)
        
        pmf = generateUniformPmfForGx(G, x_max, y_min, delta, nbins)
        
        [pmf, reachable] = computeUniformPmfForGx(G, x_max, ...
            y_min, delta, nbins)
        
        [pmf, delta] = computeReceivedPmfViaUnifThenConv( ...
            G, x_max, sigma_w, ns, nbins)
        
    end % methods(Static)
    
end

