function [mi_nats, nbins, h_y, pmf] = calculateMutualInfoWithUnifGx( ...
    G, x_max, sigma_w, max_nbins, ns)
% calculateMutualInfoWithUnifGx calculates I(x;y) for y = G*x + w with a
% uniformly-distributed G*x.  
%
%   [MI_NATS, NBINS, H_Y, PMF] = calculateMutualInfoWithUnifGx( ...
%       G, X_MAX, SIGMA_W, MAX_NBINS, NS)
%
% MI_NATS (scalar) is the calculated mutual information in nats.
% NBINS (row vector) is the number of bins along each dimension of y.  
% H_Y (scalar) is the calculated differential entropy of y.
% PMF (n_r dimension matrix of size NBINS) is the PMF of y used to
%   calculate the mutual information.  
%
% (H_Y and PMF are provided for debugging and testing purposes.)
%
% G is the channel matrix (of gains from each transmitter to each receiver
%   element).  G is a n_r by n_t matrix, where n_r is the number of
%   receiver elements and n_t is the number of transmitter elements.  
% X_MAX (scalar or n_t-element column vector, positive) is the maximum
%   value for x (from each transmitter).  If x_max is provided as a scalar,
%   it will be replicated into a n_t-element column vector.
% SIGMA_W (n_r-element column vector) is the standard deviation of the
%   independent white Gaussian noise (w) of each receiver element.  
% MAX_NBINS (scalar) is the maximum number of bins that may be used for the
%   PMF.  The PMF is stored as a matrix of double-precision floats (where
%   each bin takes 8 bytes).  
% NS (scalar, default: 6) determines the domain of the PMF generated.
%   Along each dimension d of y, where Gx_max = G*x_max, the PMF is
%   computed for: -ns*sigma_w(d) <= y(d) <= Gx_max(d) + ns*sigma_w(d).

% To avoid confusion, all vectors are stored as column vectors in this
% method; transpose them to get row vectors as necessary.

%% default parameters

% Default ns
if(nargin < 5)
    ns = 6;
end

%% process input

% G
if(~ismatrix(G))
    error('Parameter G must be a 2D matrix.')
end
[n_r, n_t] = size(G);

% x_max
if(isscalar(x_max))
    x_max = repmat(x_max, n_t, 1);
end
if(~isequal(size(x_max), [n_t, 1]))
    error('Parameter x_max must be a n_t-element column vector.');
end
if(any(x_max < 0))
    error('Parameter x_max must be non-negative.');
end

% sigma_w
if(~isequal(size(sigma_w), [n_r, 1]))
    error('Parameter sigma_w must be a n_r-element column vector.');
end
if(any(sigma_w < 0))
    error('Parameter sigma_w must be non-negative.');
end

% max_nbins
if(~isscalar(max_nbins) || max_nbins<1)
    error('Parameter MAX_NBINS must be a scalar whole number.');
end
nbins = (MIMOCapacity.fillMaxNBins(max_nbins, n_r))';

% ns
if(~isscalar(ns) || ns<0)
    error('Parameter NS must be a non-negative scalar number.');
end

%% generate the PMF of y

[pmf, delta] = MIMOCapacity.computeReceivedPmfViaUnifThenConv( ...
    G, x_max, sigma_w, ns, nbins);

%% compute the outputs

% h(y|x)
h_y_given_x = MIMOCapacity.calculateDiffEntropyOfGaussian(sigma_w.^2);

h_y = MIMOCapacity.calculateDiffEntropyFromPmf(pmf, delta);

mi_nats = h_y - h_y_given_x;

end

