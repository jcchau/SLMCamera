classdef MIMOCapacityOld
    % MIMOCapacityOld holds old methods and implementations for calculating
    % the MIMO capacity of VLC channels.  
    %
    % The methods in this class should probably NOT be used since they
    % produce incorrect results.  In particular, the following methods
    % assume incorrectly that the entropy of y is maximized using a
    % uniformly distributed x when y = G*x + w and x is bounded:
    %   - calculateCapacityForUniformInput
    %   - calculateDiffEntropyOfOutputForUniformInput
    %   - calculateVarianceOfDiffEntropyFromMonteCarloPmf
    %   - generateReceivedPmfForUniformInput
    
    properties
    end
    
    methods(Static)
        
        [nats, min_nats, variance_H, nbins, trials, h_y, pmf] = ...
            calculateCapacityForUniformInput( ...
            G, x_max, variance_noise_out, max_nbins, min_trials)
        
        variance = calculateVarianceOfDiffEntropyFromMonteCarloPmf( ...
            pmf, num_trials)
        
        [pmf, trials, y_min, y_max, hits] = ...
            generateReceivedPmfForUniformInput( ...
            G, x_max, variance_noise_out, bins_per_dimension, ...
            min_trials, trials_per_batch)
        
    end
    
end

