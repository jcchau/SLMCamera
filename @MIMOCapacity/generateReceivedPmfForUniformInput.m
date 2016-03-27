function [pmf, trials, y_min, y_max, hits] = ...
    generateReceivedPmfForUniformInput( ...
    G, x_max, variance_noise_out, bins_per_dimension, ...
    min_trials, trials_per_batch)
% generateReceivedPmfForUniformInput computes the PMF of the received
% signal for a MIMO system with uniformly distributed inputs.  
%
% Assume channel:
% y = G*x + w
%
%   [PMF, TRIALS, Y_MIN, Y_MAX, HITS] = ...
%       generateReceivedPmfForUniformInput( ...
%       G, X_MAX, VARIANCE_NOISE_OUT, BINS_PER_DIMENSION, ...
%       MIN_TRIALS, TRIALS_PER_BATCH)
%
% PMF is a n_r-dimension matrix, where each dimension has size
%   BINS_PER_DIMENSION.  PMF is the computed PMF of y.  In each dimension
%   of y, each bin_i covers values of y in
%   (y_min + (i-1)*delta) < y <= (y_min + i*delta). 
% TRIALS is a scalar for the number of trials simulated to compute the PMF.
% Y_MIN AND Y_MAX are n_r-element column vectors of the minimum and maximum
%   values of y that are included in any bin of the calculated PMF. 
% HITS is a matrix with the same size as PMF.  HITS is a count of how many
%   times a trial yields a y that lands in each bin of the PMF.  
%
% G is the channel matrix (from each transmitter to each receiver element).
%   G is a n_r by n_t matrix, where n_r is the number of reciever elements
%   and n_t is the number of transmitter elements.  
% X_MAX is the maximum value for each x (from each transmitter).  X_MAX is
%   a scalar, positive number.  
% VARIANCE_NOISE_OUT is a n_r-element column vector representing the
%   variance of the independent AWGN received by each receiver.  
% BINS_PER_DIMENSION a scalar positive integer that is the number of bins
%   to use per dimension of the (PMF); i.e., into how many bins do we
%   discretize the signal value of each receiver element.  
% MIN_TRIALS is a scalar for the minimum number of Monte Carlo trials that
%   should be used to compute the PMF of y.  
% TRIALS_PER_BATCH (optional, default=4096) is the number of trials to run
%   together for vectorization.  
%
% WARNING: This function may require a lot of memory.  Approximately
% 2 * (8 bytes * BINS_PER_DIMENSION^n_r):
% - Both PMF and HITS are each (8 bytes * BINS_PER_DIMENSION^n_r) bytes.

%% set default for optional argument

% for trials_per_batch
if(nargin < 6)
    % Select a default that's
    % - big enough to fully/efficiently utilize the CPU(s),
    % - is divisible in many ways, and
    % - is not too big (to keep memory use small and to avoid extending the
    %   computation time too long). 
    % 45*2^15 performed better than 2^20 or 2^21 when tested using
    % MIMOCapacityTest.testGenerateReceivedPmfForUniformInput23; tried 2^13
    % to 2^25.
    trials_per_batch = 1474560; % 45*2^15; 11.25 MiB. 
end

%% validate inputs and determine sizes

if(~ismatrix(G))
    error('Parameter G must be a 2-dimensional matrix.');
end
[n_r, n_t] = size(G);

if(~isscalar(x_max) || x_max<0)
    error('Parameter X_MAX must be scalar and non-negative.');
end

if(~isequal(size(variance_noise_out), [n_r, 1]))
    error(['Parameter VARIANCE_NOISE_OUT must be a ' ...
        'n_r-element column vector.']);
end
if(any(variance_noise_out < 0))
    error(['Variance must be non-negative in ' ...
        'parameter VARIANCE_NOISE_OUT.']);
end

if(~isscalar(bins_per_dimension))
    error('Parameter BINS_PER_DIMENSION must be scalar.');
end
if(~isscalar(min_trials))
    error('Parameter MIN_TRIALS must be scalar.');
end
if(~isscalar(trials_per_batch))
    error('Parameter TRIALS_PER_BATCH must be scalar.');
end

if(min_trials < 1)
    error('Parameter MIN_TRIALS must be at least 1.');
end
if(trials_per_batch < 1)
    error('Parameter TRIALS_PER_BATCH must be at least 1.');
end

%% initial setup of variables

% noise standard deviation for each receiver element.
stddev_noise_out = sqrt(variance_noise_out);

% Determine the maximum value of G*x (the maximum received signals when
% excluding noise).
Gx_max = G * repmat(x_max, n_t, 1);

% we clip the PMF of y from y_min to y_max in our appoximation
y_min = -5 * stddev_noise_out; % n_r element column matrix
y_max = Gx_max + 5*stddev_noise_out; % n_r element column matrix
delta = (y_max-y_min) ./ bins_per_dimension;

% Allocate a n_r dimension matrix to store the hit counts.
% The last dimension of 1 in the argument of zeros ensures that zeros does
% not return a square matrix when n_r is 1 (and doesn't do anything
% otherwise).
hits = zeros([repmat(bins_per_dimension, 1, n_r), 1]);
% In each dimension of y, each bin_i covers values of y in
% (y_min + (i-1)*delta) < y <= (y_min + i*delta).  
% As calculated in p.134 of my lab book #3, the index for each dimension
% index_y = ceil((y-y_min)/delta).

% We will calculate PMF by dividing hits by the total number of trials into
% batches.
batches = ceil( min_trials / trials_per_batch);
trials = batches * trials_per_batch;

% pre-compute the weights needed to convert a row of index_x into a linear
% index (needed to index elements in arbitrary-dimension matrices).
indexing_weights = MIMOCapacity.convertToLinearIndexWeights(size(hits));

%% calculate the pmf of y through Monte Carlo simulation

wb = waitbar(0, 'Starting...');

for ii = 1:batches
    
    waitbar(ii/batches, wb, sprintf('Batch %u of %u', ii, batches));
    % so the bar graphic shows batch ii completed (even though it just
    % started), but batches is so big, not worth the extra computation to
    % make this waitbar appear more accurate.
    
    % In this vectorization, each row is a trial in the batch.  
    
    x = x_max * rand(trials_per_batch, n_t);
    
    % Each row of w is the AWGN for a trial, where each column is the AWGN
    % for each dimension of y (i.e., for each receiver element).  
    w = normrnd(0, repmat(stddev_noise_out', trials_per_batch, 1), ...
        trials_per_batch, n_r);
    
    y = x * G' + w;
    
    % Convert y into corresponding matrix-subscript indices for hits.
    % In each dimension of y, index_y = ceil((y-y_min)/delta).  
    index_y = ceil(bsxfun(@rdivide, ...
        bsxfun(@minus, y, y_min'), ...
        delta'));

    % Discard any row of index_y for which any value is outside of the
    % range of valid indices for matrix hits: [1, bins_per_dimension].
    index_y = index_y( ...
        all(index_y >= 1, 2) & ...
        all(index_y <= bins_per_dimension, 2), :);
    
    % And convert the matrix subscript indices to linear indices.
    li_y = MIMOCapacity.convertToLinearIndex( ...
        indexing_weights, index_y);
    
    % Tally each hit in matrix hits
    % NOTE: if we do "for jj = li_y", it does NOT iterate through all of
    % the elements of li_y, but instead just does one iteration in the
    % loop, so it's necessary to do "for jj = 1:length(li_y)" instead.
    for jj = 1:length(li_y)
        index_bin = li_y(jj);
        hits(index_bin) = hits(index_bin) + 1;
    end % for jj = 1:length(li_y)
    
end % for ii = 1:batches

close(wb)

% Compute the PMF.
% Note that if any trials yielded a y outside of bins in hits, then the
% total of all pmf below < 1. 
pmf = hits ./ trials;

end