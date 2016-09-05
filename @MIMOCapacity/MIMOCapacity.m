classdef MIMOCapacity
    %MIMOCAPACITY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        
        % If there is at least MinOrthogonalPartFactor of the original
        % column left after orthognonalization, then this column is
        % linearly independent of the previous columns.
        % 1e-4 was chosen on p. 127 of lab book 4 with the reasoning that
        % if two channels have the same noise, but one channel is 1e-4 as
        % strong as the other channel, then the former's contribution to
        % the capacity would be negligible compared to the latter's.
        % Otherwise, consider this column to be linearly dependent on the
        % previous columns.
        MinOrthogonalPartFactor = 1e-4;
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
        G_B = simplifyChannelMatrix(G)
        Q = computeQTTransform(G)
        
        [umin, umax] = computeUExtremes(F, xmax)
        
    end % methods(Static)
    
end

