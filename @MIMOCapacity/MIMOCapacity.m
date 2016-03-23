classdef MIMOCapacity
    %MIMOCAPACITY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        de = calculateDiffEntropyFromPmf(pmf, bin_size)
        [pmf, trials, min_noise_to_delta_ratio, y_min, y_max, hits] = ...
            generateReceivedPmfForUniformInput( ...
            G, x_max, variance_noise_out, bins_per_dimension, ...
            min_trials, trials_per_batch)
        
        li = convertToLinearIndex(weights, subs)
        weights = convertToLinearIndexWeights(mat_size)
        
        out = minimumNeighbor(in, num_dimensions)
        out = maximumNeighbor(in)
        
        de = calculateMinimumDiffEntropyFromPmf(pmf, bin_size)
    end % methods(Static)
    
end

