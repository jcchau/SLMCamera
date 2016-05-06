function [ub_ElMos, ub_MVar, lb_MRC, lb_UnifX] = ...
    plotCapacityBounds(G, x_max_to_sigma_w_ratio)
%PLOTCAPACITYBOUNDS Summary of this function goes here
%   Detailed explanation goes here

fprintf('ElMos...\n');
ub_ElMos = arrayfun( ...
    @(xsr) MIMOCapacity.computeCapacityUBElMoslimany2014(G, xsr, 1), ...
    x_max_to_sigma_w_ratio);

fprintf('MVar...\n');
ub_MVar = arrayfun( ...
    @(xsr) MIMOCapacity.computeCapacityUBMaxVariance(G, xsr, 1), ...
    x_max_to_sigma_w_ratio);

fprintf('MRC...\n');
lb_MRC = arrayfun( ...
    @(xsr) MIMOCapacity.computeCapacityLBMRC(G, xsr, 1), ...
    x_max_to_sigma_w_ratio);

fprintf('UnifX...\n');
lb_UnifX = arrayfun( ...
    @(xsr) MIMOCapacity.computeCapacityLBUnifX(G, xsr, 1), ...
    x_max_to_sigma_w_ratio);

%% Plot the result
figx = 20*log10(x_max_to_sigma_w_ratio);

figure();
hold on
plot(figx, ub_ElMos, 'bs-');
plot(figx, ub_MVar, 'r^-');
plot(figx, lb_MRC, 'g+-');
plot(figx, lb_UnifX, 'co-');

legend('ub ElMos', 'ub MVar', 'lb MRC', 'lb UnifX', ...
    'Location', 'NorthWest');
xlabel('20*log_{10}(x_{max}/\sigma_w) dB')
ylabel('I(y;x) nats')
title('Capacity bounds for channel matrix G')

end

