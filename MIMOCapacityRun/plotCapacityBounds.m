function [ub_ElMos, ub_MVar, lb_MRC, lb_UnifX] = ...
    plotCapacityBounds(G, x_max_to_sigma_w_ratio)
%PLOTCAPACITYBOUNDS Summary of this function goes here
%   Detailed explanation goes here

fprintf('ElMos...\n');
ub_ElMos = arrayfun( ...
    @(xsr) MIMOCapacity.computeCapacityUBElMoslimany2014(G, xsr, 1), ...
    x_max_to_sigma_w_ratio);

fprintf('MVar...');
ub_MVar = arrayfun( ...
    @(xsr) MIMOCapacity.computeCapacityUBMaxVariance(G, xsr, 1), ...
    x_max_to_sigma_w_ratio);

fprintf('MRC...');
lb_MRC = arrayfun( ...
    @(xsr) MIMOCapacity.computeCapacityLBMRC(G, xsr, 1), ...
    x_max_to_sigma_w_ratio);

fprintf('UnifX...');
lb_UnifX = arrayfun( ...
    @(xsr) MIMOCapacity.computeCapacityLBUnifX(G, xsr, 1), ...
    x_max_to_sigma_w_ratio);

figx = 20*log10(x_max_to_sigma_w_ratio);

figure();
semilogx(figx, ub_ElMos, 'b+');
semilogx(figx, ub_MVar, 'rv');
semilogx(figx, lb_MRC, 'g^');
semilogx(figx, lb_UnifX, 'co');

legend('ub ElMos', 'ub MVar', 'lb MRC', 'lb UnifX');
xlabel('20*log_{10}(x_{max}/\sigma_w)')
ylabel('I(y;x) nats')
end

