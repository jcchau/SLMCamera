function plotCapacityWithUnifGxVsSNR(G, xmax_to_sigmaw_ratio, max_nbins)
% plotCapacityWithUnifGxVsSNR plots a loose approximation of the channel
% capacity (using a uniformly-distributed G*x) as a function of SNR.
%
% G is the channel matrix.
% XMAX_TO_SIGMAW_RATIO is the ratio of x_max to sigma_w (to set the SNR).
%   x_max is fixed at 1 and sigma_w is adjusted accordingly.  
% max_nbins is the maximum number of bins.  

[n_r, n_t] = size(G);

x_max = 1;
list_sigma_w = x_max ./ xmax_to_sigmaw_ratio;

wb = waitbar(0, '0');
for ii = 1:length(list_sigma_w)
    sigma_w = repmat(list_sigma_w(ii), n_r, 1);
    
    mi(ii) = MIMOCapacity.calculateMutualInfoWithUnifGx( ...
        G, x_max, sigma_w, max_nbins);
    
    waitbar(ii/length(list_sigma_w), wb, sprintf('%d', ii));
end % for sigma_w
close(wb)

semilogx(xmax_to_sigmaw_ratio, mi);
xlabel('x_{max} to \sigma_w ratio');
ylabel('I(x;y) (in nats)');
title('I(x;y) for y=G*x+w with uniformly distributed G*x');

end

