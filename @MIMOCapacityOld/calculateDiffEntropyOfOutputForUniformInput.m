function de = calculateDiffEntropyOfOutputForUniformInput( ...
    S, SH, x_max, variance_bg, variance_t)
% calculateDiffEntropyOfOutputForUniformInput calculates the differential
% entropy of the output of an imaging SLM VLC receiver given shot noise
% (due to background illumination) and thermal noise. 
%
% Assume channel:
% y = S(H*x + w_bg) + w_sd + w_t
% Where we ignore the signal-dependent shot noise w_sd.
%
%   DE = calculateDiffEntropyOfOutputForUniformInput( ...
%       S, SH, X_MAX, VARIANCE_BG, VARIANCE_T)
%
% DE is the calculated differential entropy in nats.  
% S is the SLM slection matrix.  
% SH is the product of S*H, where H is a matrix representing the channel
%   gain from each transmitter through each pixel (to a photodetector). 
% X_MAX is the maximum value for each x (from each transmitter). 
% VARIANCE_BG is a scalar value for the variance of shot noise of each SLM
%   pixel (w_bg).  w_bg is treated as i.i.d. AWGN. 
% VARIANCE_T is a scalar value for the variance of the thermal noise (w_t).
%   w_t is i.i.d. AWGN.  
%
% Note that the memory (and compute time) requirement for this method
% scales exponentially with n_r (the number of receiver elements).  
% n_r should not exceed 5 (recommend not exceed 4).  

% TODO: separate out the SLMCOR-specific code so that we can also use this
% for traditional imaging VLC receivers.  

% The sampling resolution for each dimension (n_r dimensions) for the PMF
% matrix.
PMF_BINS_PER_DIMENSION = 100;
% Use PMF_BINS_PER_DIMENSION = 100 for test runs
% 256 bins is 8 bits of precision for each receiver

% We run ceil(AVERAGE_TRIALS_PER_BIN * PMF_BINS_PER_DIMENSION^n_r /
% TRIALS_PER_BATCH) batches of trials.  
% Increasing AVERAGE_TRIALS_PER_BIN decreases the statistical distance
% between the simulated PMF and the actual PMF.
% We specify average trials per bin instead of total trials because
% otherwise, as the number of bins increase, the number of trials per bin
% would decrease, increasing the standard error of the probability for that
% bin.  
% 10^6 trials would take approximately 0.29 seconds on ECE-MCL-04.  Tested
% using the following code:
%   tic;
%   trials = 10^6;
%   SH = rand(5,5);
%   x = rand(trials, 5);
%   w = normrnd(0,1,trials,5);
%   Y = x * SH' + w;
%   toc

AVERAGE_TRIALS_PER_BIN = 10^4; % temporarily reduce to shorten run time
%AVERAGE_TRIALS_PER_BIN = 10^6;

% Do trials in batches to facilitate vectorization for computational
% efficiency.
TRIALS_PER_BATCH = 4096;

%% get dimensions and validate inputs
if(~ismatrix(S))
    error('Parameter S must be a 2-dimensional matrix.');
end
if(~ismatrix(SH))
    error('Parameter SH must be a 2-dimensional matrix.');
end

% n_r: number of photodetectors
% n_t: number of transmitters
% n_s: number of SLM pixels
[n_r, n_t] = size(SH);
n_s = size(S, 2);

% check that S and SH have the same number of rows (i.e., same number of
% receivers)
if(size(S,1) ~= n_r)
    error('Parameters S and SH must have the same number of rows.');
end

if(~isscalar(x_max))
    error('Parameter X_MAX must be scalar.');
end
if(~isscalar(variance_bg))
    error('Parameter VARIANCE_BG must be scalar.');
end
if(~isscalar(variance_t))
    error('Parameter VARIANCE_T must be scalar.');
end

if(x_max < 0)
    error('Parameter X_MAX must be non-negative.');
end
if(variance_bg < 0 || variance_t < 0)
    error('Parameters VARIANCE_BG and VARIANCE_T must be non-negative.');
end

% TODO: revisit memory to see how big each matrix can actually be.

% Avoid using more than 50 GiB for the PMF matrix
% (Save another 50 GiB for other variables.)
% Doubles in MATLAB are 8 bytes.
MAX_MEM_FOR_PMF = 40*2^30;
if(PMF_BINS_PER_DIMENSION^n_r * 8 > MAX_MEM_FOR_PMF)
    error(['n_r=%d is too large and would require more than ' ...
        'MAX_MEM_FOR_PMF bytes of memory.'], ...
        n_r)
end

%% allocate variables to store hits counts to calculate PMF
% In each dimension, extend 5 standard deviations of noise beyond the
% extremes of SH*x.

% calculate the Gaussian noise variance at the output
variance_noise_out = S * repmat(variance_bg, n_s, 1) + variance_t;
stddev_noise_out = sqrt(variance_noise_out);

% determine the maximum values for each dimension of SH*x (i.e., y without
% noise)
SHx_max = SH * repmat(x_max, n_t, 1);

% we clip the PMF of y from y_min to y_max in our appoximation
y_min = -5 * stddev_noise_out; % n_r element column matrix
y_max = SHx_max + 5*stddev_noise_out; % n_r element column matrix
delta = (y_max-y_min) ./ PMF_BINS_PER_DIMENSION;

% Allocate a n_r dimension matrix to store the hit counts
hits = zeros(repmat(PMF_BINS_PER_DIMENSION, 1, n_r));
% In each dimension of y, each bin_i covers values of y in
% (y_min + (i-1)*delta) < y <= (y_min + i*delta).  
% As calculated in p.134 of my lab book #3, the index for each dimension
% index_y = ceil((y-y_min)/delta).  

% calculate PMF by dividing hits by the total number of trials
batches = ceil( ...
    AVERAGE_TRIALS_PER_BIN * PMF_BINS_PER_DIMENSION^n_r / ...
    TRIALS_PER_BATCH);
trials = batches * TRIALS_PER_BATCH;

%% calculate the pmf of y through Monte Carlo simulation

% pre-compute the weights needed to convert a row of index_x into a linear
% index (needed to index elements in arbitrary-dimension matrices).

wb = waitbar(0, 'Starting...');

% pre-compute the weights needed to convert from subscript to linear
% indexing.
indexing_weights = ImagingReceiver.convertToLinearIndexWeights(size(hits));

for ii = 1:batches
    
    waitbar(ii/batches, wb, sprintf('Batch %u of %u', ii, batches));
    % so the bar graphic shows batch ii completed (even though it just
    % started), but batches is so big, not worth the extra computation to
    % make this waitbar appear more accurate.
    
    % In this vectorization, each row is a trial in the batch.  
    
    x = rand(TRIALS_PER_BATCH, n_t);
    
    % Each row of w is the AWGN for a trial, where each column is the AWGN
    % for each dimension of y (i.e., for each receiver element).  
    w = normrnd(0, repmat(stddev_noise_out', TRIALS_PER_BATCH, 1), ...
        TRIALS_PER_BATCH, n_r);
    
    y = x * SH' + w;
    
    % Convert y into corresponding indices for hits.
    % In each dimension of y, index_y = ceil((y-y_min)/delta).  
    index_y = ceil(bsxfun(@rdivide, ...
        bsxfun(@minus, y, y_min'), ...
        delta'));

    % Discard any row of index_y for which any value is outside of the
    % range of valid indices for matrix hits: [1, PMF_BINS_PER_DIMENSION].
    index_y = index_y( ...
        all(index_y >= 1, 2) & ...
        all(index_y <= PMF_BINS_PER_DIMENSION, 2), :);
    
    % And convert the matrix subscript indices to linear indices.
    li_y = ImagingReceiver.convertToLinearIndex( ...
        indexing_weights, index_y);
    
    % Tally each hit in matrix hits
    for jj = li_y
        hits(jj) = hits(jj) + 1;
    end % for jj = li_y
    
end % for ii = 1:batches

close(wb)

% Compute the PMF.
% Note that if any trials yielded a y outside of bins in hits, then the
% total of all pmf below < 1. 
pmf = hits ./ trials;

%% calculate the discrete entropy of y using the pmf

% Eliminate zero entries in pmf.
% Their contribution to ent_quantized should be 0, but MATLAB calculates it
% as NaN since we have 0*log(0) = 0*-Inf.
pmf = pmf(pmf>0);

% ent_quantized is in nats since we're using the natural log
ent_quantized = -sum(pmf(:) .* log(pmf(:)));

%% account for the factor described in Theorem 8.3.1 in Cover2005
% quantization for the pmf of y is in steps of delta

de = ent_quantized + sum(log(delta));

end
