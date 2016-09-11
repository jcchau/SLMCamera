function pmf = fillUniformPmfForGx(G, x_max, y_min, delta, nbins)
% fillUniformPmfForGx generates the uniformly-distributed PMF that covers
% the possible values of G*x where x in each dimension is in the range [0,
% x_max].
% 
% In contrast to generateUniformPmfForGx, this method recursively scans
% every discretized value of x instead of randomly generating x.
%
% Uses the channel model: y = G*x + w.
%
%   PMF = generateTransformedUniformPmf( ...
%       G, X_MAX, Y_MIN, DELTA, NBINS)
%
% PMF (n_r-dimension matrix) is the generated PMF.
%
% G (n_r*n_t matrix) is the channel matrix.  
% X_MAX (scalar or n_t-element column vector) is the maximum signal that
%   can be transmitted by each transmitter.
% Y_MIN (n_r-element row vector) is the minimum value (respectively)
%   covered covered by the PMF along each dimension. 
% DELTA (n_r-element row vector) is the size (in terms of y) of each bin of
%   the PMF.
% NBINS (n_r-element row vector) is the number of bins along each dimension
%   of y.  NBINS should consist of whole numbers.  
%
% Note that y_min, delta, and nbins should be specified to be large enough
% to cover G*x (without the noise from w).  When the PMF of G*x is
% convolved against the PMF of the noise w, the matrix would be
% correspondingly expanded (using the same bin sizes).

%% check inputs and get dimensions

% G
if(~ismatrix(G))
    error('Parameter G must be 2D matrix.');
end
[n_r, n_t] = size(G);

% x_max
if(isscalar(x_max))
    x_max = repmat(x_max, n_t, 1);
end
if(~isequal(size(x_max), [n_t, 1]))
    error(['Parameter X_MAX must be a n_t-element column vector ' ...
        '(or a scalar).']);
end
if(any(x_max<0))
    error('Parameter X_MAX should be non-negative.');
end

% y_min
if(~isequal(size(y_min), [1, n_r]))
    error('Parameter Y_MIN must be a n_r element row vector.');
end

% delta
if(~isequal(size(delta), [1, n_r]))
    error('Parameter delta must be a n_r-element row vector.');
end
if(any(delta <= 0))
    % delta should not be equal zero; otherwise, we'd have duplicate bins
    % covering the same values (assuming that we have more than one bin in
    % that dimension).
    error('Parameter delta must be positive.');
end

% nbins
if(~isequal(size(nbins), [1, n_r]))
    error('Parameter NBINS must be a n_r element row vector.');
end
if(any(nbins<1))
    error('PARAMETER NBINS should be a whole number.');
end

%% check that y_min and y_max include every possible signal w/o noise value
[umin, umax] = MIMOCapacity.computeUExtremes(G, x_max);
if(any(y_min'>umin))
    warning(['y_min is larger than umin: ' ...
        'the PMF does not include all possible values of ' ...
        'the signal without noise.']);
end
y_max = y_min + nbins .* delta; % row vector
if(any(y_max'<umax))
    warning(['y_max is smaller than umax: ' ...
        'the PMF does not include all possible values of ' ...
        'the signal without noise.']);
end

%% initial setup

% The number of trials to run together through vectorization.
% Select a value that's
%  - big enough to fully/efficiently utilize the CPU(s),
% - is divisible in many ways, and
% - is not too big (to keep memory use small and to avoid extending the
%   computation time too long). 
trials_per_batch = 1474560; % 45*2^15; 11.25 MiB per double vector.

% Allocate matrix to indicate whether a bin has been reached.
% The last dimension of 1 in the argument of zeros ensures that zeros does
% not return a square matrix when n_r is 1 (and doesn't do anything
% otherwise).
hits = false([nbins, 1]);

% pre-compute the weights needed to convert a row of index_x into a linear
% index (needed to index elements in arbitrary-dimension matrices).
indexing_weights = MIMOCapacity.convertToLinearIndexWeights(nbins);

%% include x_max into G
% rand() gives uniformly distributed random values in the interval (0,1).
% However, we want x to be uniformly distributed in the interval (0, x_max)
% for each dimension.  
% We generate many many values for x; rather than multiply each one by
% x_max to get the appropriate range of random values, multiple each column
% of G by the x_max for that transmitter.  This should yield the same
% result as multiplying each column of x by the corresponding column of
% x_max before multiplying by G.  

% x_max a column vector; transpose it to get a row vector.
G = bsxfun(@times, G, x_max');

% Now, x essential runs from [0,1].  

%% For each transmitter, the number of steps we want in our sweep

% The change in y (in units of bins) that a unit change in w would produce.
numBinsPerUnitX = bsxfun(@rdivide, G, delta');

% Number of steps that we want to take from x=0 to x=1 inclusive (to not
% skip over any bins in y). (row vector)
stepsX = ceil(max(numBinsPerUnitX, [], 1)) + 1;

deltaX = 1./(stepsX-1);

%% Go through each of the steps
% Trying to avoid expensive recursive function calls.

totalSteps = prod(stepsX);

% The process of converting from ii=0:totalSteps-1 to ix is the same as
% converting from linear matrix indices to matrix subscript indices (except
% that here, we want to end with zero-indexed indices).  
tS_to_ix_weights = MIMOCapacity.convertToLinearIndexWeights(stepsX);

% TODO: vectorize this for loop (but be careful about running out of
% memory).
for ii = 1:totalSteps
    
    % Compute the corresponding index for the step of x.
    % ix is zero-indexed so that x=(ix.*deltaX)' (to get x as a column
    % vector).
    ix = MIMOCapacity.convertLinearToSubscriptIndex( ...
        tS_to_ix_weights, ii) - 1;
    
    x = (ix .* deltaX)';
    
    y = G * x;
    
    si_y = MIMOCapacity.convertPointToSubscriptIndex(y', ...
        y_min, delta, nbins);
    li_y = MIMOCapacity.convertToLinearIndex(indexing_weights, si_y);
    
    % Mark the corresponding bin as hit.
    hits(li_y) = true;
    
end % for ii

%% return a uniform PMF that reaches the same bins

% Probability for each bin hit.
unif_prob_per_bin = 1/nnz(hits);

% Generate a PMF that is unif_prob_per_bin where G*x reaches and 0
% elsewhere.
pmf = zeros([nbins, 1]);
pmf(hits) = unif_prob_per_bin;

end
