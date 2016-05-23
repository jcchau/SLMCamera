classdef MIMOCapacityLBUnifX
    % MIMOCapacityLBUnifX holds methods and implementations for calculating
    % a lower bound for the MIMO capacity of VLC channels using a
    % uniformly-distributed input x.
    
    properties
    end
    
    methods(Static)
        
        [nats, min_nats, variance_H, nbins, trials, h_y, pmf] = ...
            calculateCapacityLBForUniformInput( ...
            G, x_max, variance_noise_out, max_nbins, min_trials)
        [lb, Q, h_x, h_QTransposeGx, h_QTransposew] = ...
            approximateUnifXLBForNegligibleNoiseAndFullColRankG( ...
            G, x_max, sigma_w)
        
        variance = calculateVarianceOfDiffEntropyFromMonteCarloPmf( ...
            pmf, num_trials)
        
        [pmf, trials, y_min, y_max, hits] = ...
            generateReceivedPmfForUniformInput( ...
            G, x_max, variance_noise_out, bins_per_dimension, ...
            min_trials, trials_per_batch)
        
    end
    
end

