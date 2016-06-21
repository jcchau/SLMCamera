function [lb_nats_conservative, lb_nats] = ...
    computeCapacityLBUnifX(G, x_max, sigma_w)
%COMPUTECAPACITYLBUNIFX Summary of this function goes here
%   Detailed explanation goes here
%
% G is the channel matrix.
% x_max (scalar) is the maximum value of x.
% sigma_w (scalar) is the standard deviation of the noise w.

[n_r, ~] = size(G);

if(x_max >= 100*sigma_w)
    try
        lb_nats = ...
            MIMOCapacityLBUnifX. ...
            approximateUnifXLBForNegligibleNoiseAndFullColRankG( ...
            G, x_max, sigma_w);
        lb_nats_conservative = lb_nats;
        return
    catch me
        switch me.identifier
            case {['MIMOCapacityLBUnifX:approximateUnifXLBFor' ...
                    'NegligibleNoiseAndFullColRankG:GNotFullColRank'], ...
                    ['MIMOCapacityLBUnifX:approximateUnifXLBFor' ...
                    'NegligibleNoiseAndFullColRankG:SignificantNoise']}
                % Ok. Just use the slower method below instead.
            otherwise
                rethrow(me)
        end % switch me.identifier
    end % try-catch
end % if(x_max >= 100*sigma_w)

variance_noise_out = repmat(sigma_w^2, n_r, 1);

max_nbins = 1e9/8/5;
min_trials = 1e9;

[lb_nats, lb_nats_conservative] = ...
    MIMOCapacityLBUnifX.calculateCapacityLBForUniformInput( ...
    G, x_max, variance_noise_out, max_nbins, min_trials);

end

