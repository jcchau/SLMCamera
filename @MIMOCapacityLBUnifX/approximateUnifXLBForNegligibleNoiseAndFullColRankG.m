function [lb, Q, h_x, h_QTransposeGx, h_QTransposew] = ...
    approximateUnifXLBForNegligibleNoiseAndFullColRankG( ...
    G, x_max, sigma_w)
% approximateUnifXLBForNegligibleNoiseAndFullColRankG approximates a lower
% bound on the capacity using a uniformly distributed x by assuming that
% the noise is negligible, given that the channel matrix G is full column
% rank.  
%
%   [lb, Q, h_x, h_QTransposeGx, h_QTransposew] = ...
%       approximateUnifXLBForNegligibleNoiseAndFullColRankG( ...
%           G, x_max, sigma_w)
%
% lb (scalar) is I(x;y) (as the lower bound on capacity) in nats computed
%   using the assumption that the noise is negligible.
% Q, h_x, h_QTransposeGx, and h_QTransposew are provided for debugging.
%
% G is the channel matrix.
% x_max (scalar) is the maximum transmit signal value for each transmitter.
% sigma_w (scalar) is the standard deviation of the noise w of each
%   receiver element.  
%
% Implemented as described in lab book #4 p.108-111.

[~, n_t] = size(G);

% Check that G is full column rank.
if(rank(G) < n_t)
    error(['MIMOCapacityLBUnifX:approximateUnifXLBForNegligibleNoise' ...
        'AndFullColRankG:GNotFullColRank'], ...
        'Parameter G must have full column rank.');
end

% "economy size" QR-decomposition of G.
% Replaced with Gram-Schmidth orthogonalization to keep A=Q'*G positive.
% We need Q s.t. A is positive for
% testApproximateUnifXLBForNegligibleNoiseAndFullColRankG3x2.
Q = MIMOCapacityLBUnifX.gramSchmidt(G);

A = Q'*G;

% Skip checking that the noise is negligible since the lower bound on
% capacity would still be valid even with significant noise.  

% Derived in lab book #4, p.111-112.
h_x = n_t * log(x_max);

% h(Q'*G*x)
h_QTransposeGx = h_x + log(abs(det(A)));

% h(Q'*w) is the differential entropy of n_t independent Gaussian RVs, each
% with variance sigma_w^2.  
K_QTransposew = sigma_w^2 .* eye(n_t);
h_QTransposew = MIMOCapacity.calculateDiffEntropyOfMVGaussian( ...
    K_QTransposew);

% I(x;y) is approximately h(Q'*G*x) - h(Q'*w).
lb = h_QTransposeGx - h_QTransposew;

end

