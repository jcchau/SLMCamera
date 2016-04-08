classdef MIMOCapacity
    %MIMOCAPACITY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        de = calculateDiffEntropyFromPmf(pmf, bin_size)
        
        [pmf, trials, y_min, y_max, hits] = ...
            generateReceivedPmfForUniformInput( ...
            G, x_max, variance_noise_out, bins_per_dimension, ...
            min_trials, trials_per_batch)
        
        li = convertToLinearIndex(weights, subs)
        weights = convertToLinearIndexWeights(mat_size)
        
        out = minimumNeighbor(in, num_dimensions)
        out = maximumNeighbor(in)
        
        de = calculateMinimumDiffEntropyFromPmf(pmf, bin_size)
        
        de = calculateDiffEntropyOfGaussian(variance)
        
        nbins = fillMaxNBins(max_nbins, dimensions)
        
        [nats, min_nats, variance_H, nbins, trials, h_y, pmf] = ...
            calculateCapacityForUniformInput( ...
            G, x_max, variance_noise_out, max_nbins, min_trials)
        
        variance = calculateVarianceOfDiffEntropyFromMonteCarloPmf( ...
            pmf, num_trials)
        
        [de, pmf] = calculateDiffEntropyOfClippedNormal( ...
            a, b, nbins, mu, sigma)
        
        out = removeZeroRowsAndCols(in)
        
        msi = convertPointToSubscriptIndex(y, ymin, delta, nbins)
        
        pmf = generateTransformedUniformPmf( ...
            G, x_max, nbins, ymin, ymax)
    end % methods(Static)
    
end

