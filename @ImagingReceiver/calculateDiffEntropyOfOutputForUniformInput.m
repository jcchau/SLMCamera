function diffent = calculateDiffEntropyOfOutputForUniformInput( ...
    S, SH, x_max, variance_bg, variance_t)
% calculateDiffEntropyOfOutputForUniformInput calculates the differential
% entropy of the output of an imaging SLM VLC receiver given shot noise
% (due to background illumination) and thermal noise. 
%
% Assume channel:
% y = S(H*x + w_bg) + w_sd + w_t
% Where we ignore the signal-dependent shot noise w_sd.
%
%   DIFFENT = calculateDiffEntropyOfOutputForUniformInput( ...
%       S, SH, X_MAX, VARIANCE_BG, VARIANCE_T)
%
% DIFFENT is the calculated differential entropy in nats.  
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

% The sampling resolution for each dimension (n_r dimensions) for the PMF
% matrix.
PMF_SAMPLES_PER_DIMENSION = 100;

% The sampling resolution for each dimension of x (n_t) dimensions.
X_SAMPLES_PER_DIMENSION = 100;

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

% Avoid using more than 50 GiB for the PMF matrix
% (Save another 50 GiB for other variables.)
MAX_MEM_FOR_PMF = 50*2^30;
if(PMF_SAMPLES_PER_DIMENSION^n_r * 8 > MAX_MEM_FOR_PMF)
    error('n_r=%d is too large and would require more than MAX_MEM_FOR_PMF bytes of memory.', ...
        n_r)
end

%% calculate the pmf of y through numerical convolution

% calculate the Gaussian noise variance at the output
variance_noise_out = S * repmat(variance_bg, n_s, 1) + variance_t;
stddev_noise_out = sqrt(variance_noise_out);

% determine the maximum values for each SH*x (i.e., y without noise)
SHx_max = SH * repmat(x_max, n_t, 1);

% we clip the PMF of y from y_min to y_max in our appoximation
y_min = -5 * stddev_noise_out; % n_r element column matrix
y_max = SHx_max + 5*stddev_noise_out; % n_r element column matrix
delta = (y_max-y_min) ./ PMF_SAMPLES_PER_DIMENSION;

% Allocate a n_r dimension matrix to store the PMF approximating the PDF of
% y.  
% The PMF (we divide by the number of sample points at the end).  
% Pre-allocate a n_r dimension matrix of doubles where each dimension has
% 100 elements. 
pmf = zeros(repmat(PMF_SAMPLES_PER_DIMENSION, 1, n_r));

% Compute the sample points for each dimension of x.
% X_SAMPLE_PER_DIMENSION samples for each dimension.
delta_x = x_max/X_SAMPLES_PER_DIMENSION;
x_sample_points = 0:delta_x:x_max;
% Center the first X_SAMPLE_PER_DIMENSION sample points (the previous
% command gave us X_SAMPLES_PER_DIMENSION+1 points).
x_sample_points = x_sample_points(1:X_SAMPLES_PER_DIMENSION) + delta_x/2;

% generate a grid for pmf
y_values = cell(1, n_r);
for ii = 1:n_r
    y_values{ii} = ...
        linspace(y_min(ii), y_max(ii), PMF_SAMPLES_PER_DIMENSION);    
end % for ii = 1:n_r
% DON'T generate a grid since the resulting grid would be 3x as large as
% variable pmf (and our algorithm is constrained by memory):
% [grid_y{1:n_r}] = ndgrid(y_values{:});

wb = waitbar(0, 'Starting...');
wb_max = X_SAMPLES_PER_DIMENSION^n_t;

for ii = 1:X_SAMPLES_PER_DIMENSION^n_t
    % map ii to a n_t dimension sample point, x
    [sample_index{1:n_t}] = ...
        ind2sub(repmat(X_SAMPLES_PER_DIMENSION, 1, n_t), ii);
    x = zeros(n_t, 1);
    for jj = 1:n_t
        x(jj) = x_sample_points(sample_index{jj});
    end % for jj = 1:n_t
    
    waitbar((ii-1)/wb_max, wb, ...
        sprintf('For x at point %d of %d', ii, wb_max))
    
    for jj = 1:PMF_SAMPLES_PER_DIMENSION^n_r
        [pmf_index{1:n_r}] = ...
            ind2sub(repmat(PMF_SAMPLES_PER_DIMENSION, 1, n_r), jj);
        
        % calculate the n_r-dimension y value corresponding to this PMF
        % sample.
        y = zeros(1, n_r);
        for kk = 1:n_r
            y(kk) = y_values{kk}(pmf_index{kk});
        end % for kk = 1:n_r
        
        % for mvnpdf, y (n_r-dim row vector) is where the pdf is evaluated,
        % (SH*x)' is a n_r-dim row vector representing the mean of the
        % joint-Gaussian RV, and variance_noise_out' is a n_r-dim row
        % vector for the variance along each of the n_r dimensions.  
        % Each of the n_r dimensions are independently distributed.  
        pmf(pmf_index{:}) = pmf(pmf_index{:}) + ...
            mvnpdf(y, (SH*x)', variance_noise_out');
    end % for jj = 1:PMF_SAMPLES_PER_DIMENSION^n_r
end % for ii

close(wb)

% Divide by the number of x samples used to account for the pdf of x (which
% we didn't factor in before to save a little computation).
pmf = pmf ./ X_SAMPLES_PER_DIMENSION^n_t;

%% calculate the entropy of y using the pmf

% ent_quantized is in nats
ent_quantized = -sum(pmf(:) .* log(pmf(:)));

%% account for the factor described in Theorem 8.3.1 in Cover2005
% quantization for the pmf of y is in steps of delta

diffent = ent_quantized + sum(log(delta));

end

