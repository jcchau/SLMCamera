function [nats, min_nats, variance_H, nbins, trials] = ...
    calculateCapacityForUniformInput( ...
    G, x_max, variance_noise_out, max_nbins, min_trials)
% calculateCapacityForUniformInput calculates the capacity of a MIMO
% communication channel with uniformly distributed inputs and with AWGN at
% the output.  
%
% Assume channel:
% y = G*x + w
%
%   [NATS, MIN_NATS, VARIANCE_H, NBINS, TRIALS] = ...
%       calculateCapacityForUniformInput( ...
%       G, X_MAX, VARIANCE_NOISE_OUT, MAX_NBINS, MIN_TRIALS)
%
% NATS (scalar) is the upper bound on the capacity (assuming that the
%   generated PMF is correct).
% MIN_NATS (scalar) is the lower bound on the capacity (assuming that the
%   generated PMF is correct). 
% VARIANCE_H (scalar) is the variance in the error of entropy of y (due to
%   the variances in the PMF because the PMF is generated through a Monte
%   Carlo simulation).  
% NBINS (row vector) is number of bins along each dimension of y.  
% TRIALS (scalar) is the number of trials simulated to compute the PMF.
%
% G is the channel matrix (from each transmitter to each receiver element).
%   G is a n_r by n_r matrix, where n_r is the number of received elements
%   and n_t is the number of transmitter elements.  
% X_MAX (scalar, positive) is the maximum value for each x (from each
%   transmitter). 
% VARIANCE_NOISE_OUT (n_r-element column vecotr) is the variance of the
%   independent white Gaussian noise (w) of each receiver.  
% MAX_NBINS (scalar) is the maximum number of bins that may be used for the
%   PMF.  The PMF is stored as a matrix of double-precision floats (where
%   each bin takes 8 bytes).  During the operation of this method, there
%   may be up to 5 matrices of doubles of the same size as the PMF, so it
%   is recommended to set MAX_NBINS < MAX_MEM / 8 / 5, where MAX_MEM is the
%   maximum memory available for this method.  Note that this method also
%   needs to store other variables (though they are much smaller than the
%   PMF) in memory.
% MIN_TRIALS is the minimum number of trials to use for
%   generateReceivedPmfForUniformInput.
%
% Information is measured in units of nats (instead of bits) unless
% otherwise specified.  
%
% This method is designed to make calculating the capacity from the channel
% matrix, the maximum input value, and the noise at the output simple.
% However, better performance may be 

%% validate inputs and determine sizes

if(~ismatrix(G))
    error('Parameter G must be a 2-dimensional matrix.');
end
[n_r, n_t] = size(G);

if(~isscalar(max_nbins) || max_nbins<1)
    error('Parameter MAX_NBINS must be a scalar whole number.');
end
nbins = MIMOCapacity.fillMaxNBins(max_nbins, n_r);

%% generate the PMF

[pmf, trials, y_min, y_max, ~] = ...
    MIMOCapacity.generateReceivedPmfForUniformInput( ...
    G, x_max, variance_noise_out, nbins, min_trials);

bin_size = (y_max - y_min) ./ nbins';

%% compute the outputs

h_y_given_x = MIMOCapacity.calculateDiffEntropyOfGaussian( ...
    variance_noise_out);

% Computed according to p.199 of lab book #3 (Imaging Receivers &
% Photodetector Arrays).
variance_H = sum( pmf(:).*(1-pmf(:)) .* (-log(pmf(:))-1).^2 ./ n );

h_y = MIMOCapacity.calculateDiffEntropyFromPmf(pmf, bin_size);
min_h_y = MIMOCapacity.calculateMinimumDiffEntropyFromPmf(pmf, bin_size);

nats = h_y - h_y_given_x;
min_nats = min_h_y - h_y_given_x;

end
