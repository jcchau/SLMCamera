classdef MIMOCapacity
    %MIMOCAPACITY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        de = calculateDiffEntropyOfOutputForUniformInput( ...
            S, SH, x_max, variance_bg, variance_t)
        [pmf, trials, min_noise_to_delta_ratio, y_min, y_max, hits] = ...
            generateReceivedPmfForUniformInput( ...
            G, x_max, variance_noise_out, bins_per_dimension, ...
            min_trials, trials_per_batch)
        
        li = convertToLinearIndex(weights, subs)
        weights = convertToLinearIndexWeights(mat_size)
    end % methods(Static)
    
end

