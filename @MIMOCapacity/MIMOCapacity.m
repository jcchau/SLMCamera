classdef MIMOCapacity
    %MIMOCAPACITY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        
        C = computeSmithScalarCapacity(x_max, sigma_w)
        
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
        
        [out, nz_rows, nz_cols] = removeZeroRowsAndCols(in)
        
        msi = convertPointToSubscriptIndex(y, ymin, delta, nbins)
        
        ub_nats = computeCapacityUBElMoslimany2014(G, x_max, sigma_w)
        ub_nats = computeCapacityUBMaxVariance(G, x_max, sigma_w)
        lb_nats = computeCapacityLBMRC(G, x_max, sigma_w)
        [lb_nats_conservative, lb_nats] = ...
            computeCapacityLBUnifX(G, x_max, sigma_w)
        lb_nats = computeCapacityLBUnifGx(G, x_max, sigma_w)
        
        nats = calculateDiffEntropyOfMVGaussian(K)
        
        % To remove extra dimensions from the channel matrix.
        [Q, G_B] = simplifyChannelMatrix(G)
        Q = computeQTTransform(G)
        
    end % methods(Static)
    
end

