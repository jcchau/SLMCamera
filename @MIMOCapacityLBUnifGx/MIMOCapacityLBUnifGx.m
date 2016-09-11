classdef MIMOCapacityLBUnifGx
    %MIMOCAPACITYLBUNIFGX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        [mi_nats, nbins, h_y, pmf] = calculateMutualInfoWithUnifGx( ...
            G, x_max, sigma_w, max_nbins, ns)
        [pmf, delta] = computeReceivedPmfViaUnifThenConv( ...
            G, x_max, sigma_w, ns, nbins)
        [pmf, reachable] = computeUniformPmfForGx(G, x_max, ...
            y_min, delta, nbins)
        pmf = generateUniformPmfForGx(G, x_max, y_min, delta, nbins)
        pmf = fillUniformPmfForGx(G, x_max, y_min, delta, nbins)
    end % methods(Static)
    
end

