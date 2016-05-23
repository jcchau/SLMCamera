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

% "economy size" QR-decomposition of G
[Q, ~] = qr(G, 0);

A = Q'*G;

% Sanity check that the noise can be considered negligible.
Ax_max = A * repmat(x_max,n_t,1);
if(any(abs(Ax_max) < 10*sigma_w))
    error(['MIMOCapacityLBUnifX:approximateUnifXLBForNegligibleNoise' ...
        'AndFullColRankG:SignificantNoise'], ...
        ['At least one element of Q''*G*x is not much greater than ' ...
        'sigma_w. It''s not reasonable to assume that the noise is ' ...
        'negligible.']);
end
% Though even if the sanity check passes, the noise can still be considered
% significant.  Fortunately, the lower bound on capacity would still be
% valid, just less tight (because we end up underestimating h(Q'*G*x+Q'*w)
% as h(Q'*G*x)).

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

