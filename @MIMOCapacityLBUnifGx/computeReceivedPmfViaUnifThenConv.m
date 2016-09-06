function [pmf, delta] = computeReceivedPmfViaUnifThenConv( ...
    G, x_max, sigma_w, ns, nbins)
% computeReceivedPmfViaUnifThenConv computes the PMF (of y = G*x + w, with
% bounded values of x and independent AWGN w, for uniformly-distributed
% G*x) by first computing the uniformly-distributed PMF of G*x, and then
% convolving this PMF against the noise PMF.
%
% Channel model: y = G*x + w.
%
%   pmf = computeReceivedPmfViaUnifThenConv(G, x_max, sigma_w, nbins)
%
% pmf (n_r-dimension matrix of size nbins) is the PMF of y given a
%   uniformly-distributed G*x. 
% delta (n_r-element column vector) is the size of each bin in the PMF.
%
% G (n_r by n_t matrix) is the channel matrix. 
% x_max (n_t-element column vector) is a vector of the maximum values of x.
%   For each of the n_t transmitters, 0 <= x <= x_max.  If x_max is
%   provided as a scalar, it will be replicated into a n_t-element column
%   vector.
% sigma_w (n_r-element column vector) is the standard deviation of w.  
% ns (scalar) determines the domain of the PMF generated.  Along each
%   dimension d of y, where umin and umax are the min and max values of the
%   noise-free component of the received signal, the PMF is computed for
%   umin-ns*sigma_w(d) <= y(d) <= umax + ns*sigma_w(d).
% nbins (n_r-element column vector) is the number of bins along each
%   dimension of the output pmf.  

% To avoid confusion, all vectors are stored as column vectors in this
% method; transpose them to get row vectors as necessary.

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

% nbins
if(~isequal(size(nbins), [n_r, 1]))
    error('Parameter nbins must be a n_r-element column vector.');
end
if(any(nbins < 1))
    error('Parameter nbins must consist of whole numbers.');
end

%% figure out how big everything should be

[umin, umax] = MIMOCapacity.computeUExtremes(G, x_max);
y_min = umin - ns .* sigma_w;
y_max = umax + ns .* sigma_w;
delta = (y_max - y_min) ./ nbins;

if(any(delta == 0))
    % Channel matrix G should be simplified by the calling function that
    % computes the entropy to eliminate dimensions of y that have zero
    % variance.  Eliminating these dimensions should not affect the
    % entropy, but keeping these dimensions (receivers) would needlessly
    % increase computation time and cause problems due to zero bin sizes
    % (using non-zero bin sizes along these dimensions would reduce the
    % precision of the PMF and the calculated entropy).
    error(['delta is zero for one dimension of the PMF. ' ...
        'Please eliminate receivers that receive nothing.']);
end

% pre-allocate the pmf to be returned
% The last singleton dimension is so that MATLAB does not make a square pmf
% if n_r == 1.
pmf = zeros([nbins' 1]);

% Only have computeUniformPmfForGx compute the pmf for bins in the
% rectangular prism between umin and umax because the computation time for
% that method is proportional to the number of bins (and because that
% method is relatively slow).  

% Index of 0 and Gx_max in the pmf for this method
% (computeReceivedPmfViaUnifThenConv).
i_unif_min = (MIMOCapacity.convertPointToSubscriptIndex( ...
    umin', y_min', delta', nbins'))';
i_unif_max = (MIMOCapacity.convertPointToSubscriptIndex( ...
    umax', y_min', delta', nbins'))';

% Parameters for method computeUniformPmfForGx so that the bins in that
% method align with the bins in this method and so we don't process more
% bins than we need to in computeUniformPmfForGx.
unif_y_min = y_min + (i_unif_min-1) .* delta;
unif_nbins = i_unif_max - i_unif_min + 1;

subpmf_indices = cell(n_r, 1);
for d = 1:n_r
    subpmf_indices(d) = {i_unif_min(d):i_unif_max(d)};
end % d = 1:n_r

%% compute the uniform PMF for G*x

pmf(subpmf_indices{:}) = MIMOCapacityLBUnifGx.computeUniformPmfForGx( ...
    G, x_max, ...
    unif_y_min', delta', unif_nbins');

%% add the Gaussian noise to the PMF

% To ensure that the PMF of w does not need to be zero-padded when
% convolved against the non-zero probabilities of the uniform PMF of G*x,
% the PMF of w needs to be w_nbins long.  
w_nbins = nbins + i_unif_max - i_unif_min;

% Ensure that the length of the PMF of w is odd so that conv (or convn) can
% clip off the same number of samples on each end (so that the result of
% the convolution remains centered and so that the ends of the convolved
% PMF are clipped symmetrically).
% Knowing whether w_nbins is even or odd also makes testing easier (by
% removing some variability).
even_w_nbins = mod(w_nbins, 2) == 0;
w_nbins(even_w_nbins) = w_nbins(even_w_nbins) + 1;

% Now for each dimension of w:
% - Generate the PMF of w along that dimension, symmetric about w(d) = 0
% with w_nbins, using the same bin size specified in delta.
% - Perform the convolution of the PMF of y (so far) against the PMF of w
% along dimension d, and use the result as the PMF of y so far for the
% later dimensions.  
for d = 1:n_r
    w_max = w_nbins(d) * delta(d) / 2;
    w_min = -w_max;
    bin_boundaries = (linspace(w_min, w_max, w_nbins(d)+1))';
    cdf_w_at_bb = normcdf(bin_boundaries, 0, sigma_w(d));
    pmf_w = cdf_w_at_bb(2:end) - cdf_w_at_bb(1:end-1);
    
    % reshape pmf_w to be along dimension d
    shape_pmf_w = ones(1, n_r+1); % row vector; add last singleton dim in case n_r is 1.
    shape_pmf_w(d) = w_nbins(d);
    pmf_w = reshape(pmf_w, shape_pmf_w);
    
    % Perform the convolution (keeping the same size because pmf is already
    % the desired size).
    pmf = convn(pmf, pmf_w, 'same');
end % for d = 1:n_r

end

