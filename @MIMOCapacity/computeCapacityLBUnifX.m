function [lb_nats_conservative, lb_nats] = ...
    computeCapacityLBUnifX(G, x_max, sigma_w)
%COMPUTECAPACITYLBUNIFX Summary of this function goes here
%   Detailed explanation goes here
%
% G is the channel matrix.
% x_max (scalar) is the maximum value of x.
% sigma_w (scalar) is the standard deviation of the noise w.

%% Prepare the channel matrix

G = MIMOCapacity.simplifyChannelMatrix(G);

% Handle the special case where G is 0.
if(isequal(G,0))
    lb_nats_conservative = 0;
    lb_nats = 0;
    return
end

% Apply the Q' transform
G = Q' * G';
% After the Q'-transform, the channel matrix has full row rank.  

[n_r, n_t] = size(G);

%% Compute the lower bound on capacity

if(n_r == n_t)
    % If n_t == n_r, then we also have full column rank and the channel
    % matrix is full rank (square and invertible).  
    
    % This lower bound would be tighter when the SNR is higher (so the
    % noise does not contribute significantly to h(z) for the channel
    % z=G*x+w.  (Here, we are using the post-Q'-transform G and assuming
    % that the noise is isotropic).  
    % However, the lower bound on capacity is still valid even if the SNR
    % is not very high; in this case, the lower bound would be looser.  
    
    lb_nats = ...
            MIMOCapacityLBUnifX. ...
            approximateUnifXLBForNegligibleNoiseAndFullColRankG( ...
            G, x_max, sigma_w);
    lb_nats_conservative = lb_nats;
    
else % if(n_r == n_t)
    % In this case, n_t > n_r, so we do not have full column rank and we
    % cannot use the more computationally-efficient algorithm above.  
    % Instead, we do it the slow way by calculating the PMF.  
    % We still use the channel matrix after the Q'-transform here to get
    % rid of extra (flat) dimension in the received signals that would make
    % the bins too large along those dimensions.  
    
    variance_noise_out = repmat(sigma_w^2, n_r, 1);

    max_nbins = 1e9/8/5;
    min_trials = 1e9;

    [lb_nats, lb_nats_conservative] = ...
        MIMOCapacityLBUnifX.calculateCapacityLBForUniformInput( ...
        G, x_max, variance_noise_out, max_nbins, min_trials);
end % else of if(n_r == n_t)

end

