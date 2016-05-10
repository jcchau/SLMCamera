% analyzeData.m

ntrials = 100;

%% Settings from runTrial

% The ratios of x_max to sigma_w in which to run evaluate capacity.
snr_dB = [60, 80, 100];
snr = 10.^(snr_dB/20);
nsnr = length(snr);

% for the traditional imaging VLC receivers
npix_rows = [ 1, 2, 3, 10, 100 ];
npix_cols = [ 2, 2, 3, 10, 100 ];
npix = npix_rows .* npix_cols;
nres = length(npix);

%% load the data
[trial_rank_dmd, trial_cn_dmd, ...
    trial_ub_ElMos_dmd, trial_lb_MRC_dmd, trial_lb_UnifX_dmd, ...
    trial_lb_UnifGx_dmd, ...
    trial_rank, trial_cn, ...
    trial_ub_ElMos, trial_lb_MRC, trial_lb_UnifX, ...
    trial_lb_UnifGx] = ...
    aggregateResults(ntrials);

%% Average rank & median condition number

fprintf('\n**Average Rank**\n');
mean_rank_dmd = mean(trial_rank_dmd);
mean_rank = mean(trial_rank);
fprintf('                DMD: %f\n', mean_rank_dmd);
for ii = 1:nres
    fprintf('%3dx%3d traditional: %f\n', ...
        npix_cols(ii), npix_rows(ii), mean_rank(ii));
end

fprintf('\n**Median Condition Number**\n');
median_cn_dmd = median(trial_cn_dmd);
median_cn = median(trial_cn);
fprintf('                DMD: %f\n', median_cn_dmd);
for ii = 1:nres
    fprintf('%3dx%3d traditional: %f\n', ...
        npix_cols(ii), npix_rows(ii), median_cn(ii));
end

%% Capacity bounds
% Convert to bits

% DMD
mean_ub_ElMos_dmd = mean(trial_ub_ElMos_dmd) ./ log(2);
mean_lb_MRC_dmd = mean(trial_lb_MRC_dmd) ./ log(2);
mean_lb_UnifGx_dmd = mean(trial_lb_UnifGx_dmd) ./ log(2);
mean_lb_UnifX_dmd = mean(trial_lb_UnifX_dmd) ./ log(2);
mean_lb_dmd = max(mean_lb_MRC_dmd, ...
    max(mean_lb_UnifGx_dmd, mean_lb_UnifX_dmd));

figure(1);
plot(snr_dB, mean_ub_ElMos_dmd, 's-', ...
    snr_dB, mean_lb_MRC_dmd, 'x-', ...
    snr_dB, mean_lb_UnifGx_dmd, '>-', ...
    snr_dB, mean_lb_UnifX_dmd, 'o-', ...
    snr_dB, mean_lb_dmd, ':');
xlabel('20*log_{10}(x_{max}/\sigma_w) dB');
ylabel('I(x;y) bits');
grid on;
legend('UB ElMos', 'LB MRC', 'LB UnifGx', 'LB UnifX', 'LB', ...
    'Location', 'NorthWest');
title('Average capacity bounds for the SLM VLC Receiver');

% Traditional imaging
mean_ub_ElMos = mean(trial_ub_ElMos) ./ log(2);
mean_lb_MRC = mean(trial_lb_MRC) ./ log(2);
mean_lb_UnifGx = mean(trial_lb_UnifGx) ./ log(2);
mean_lb_UnifX = mean(trial_lb_UnifX) ./ log(2);

mean_ub_ElMos = reshape(mean_ub_ElMos, nres, nsnr);
mean_lb_MRC = reshape(mean_lb_MRC, nres, nsnr);
mean_lb_UnifGx = reshape(mean_lb_UnifGx, nres, nsnr);
mean_lb_UnifX = reshape(mean_lb_UnifX, nres, nsnr);

% Drop UnifGx and UnifX for npix>20.
mean_lb = zeros(nres, nsnr);
mean_lb(npix<=20,:) = max(mean_lb_MRC(npix<=20,:), ...
    max(mean_lb_UnifGx(npix<=20,:), mean_lb_UnifX(npix<=20,:)));

% Plot each traditional imaging VLC receiver
for ires = 1:nres
    figure(ires+1); % the DMD plot is figure 1
    if(npix(ires) <= 20)
        plot(snr_dB, mean_ub_ElMos(ires,:), 's-', ...
            snr_dB, mean_lb_MRC(ires,:), 'x-', ...
            snr_dB, mean_lb_UnifGx(ires,:), '>-', ...
            snr_dB, mean_lb_UnifX(ires,:), 'o-', ...
            snr_dB, mean_lb(ires,:), ':');
        legend('UB ElMos', 'LB MRC', 'LB UnifGx', 'LB UnifX', 'LB', ...
            'Location', 'NorthWest');
    else
        % Omit UnifGx and UnifX because the bins are spread out over too
        % many dimensions.  As a result, the bins are not small enough to
        % neglect the effects of the boundary bins (on increasing the
        % differential entropy).  
        plot(snr_dB, mean_ub_ElMos(ires,:), 's-', ...
            snr_dB, mean_lb_MRC(ires,:), 'x-', ...
            snr_dB, mean_lb(ires,:), ':');
        legend('UB ElMos', 'LB MRC', 'LB', ...
            'Location', 'NorthWest');
    end
    xlabel('20*log_{10}(x_{max}/\sigma_w) dB');
    ylabel('I(x;y) bits');
    grid on;
    title(sprintf(['Average capacity bounds for the %dx%d ' ...
        'traditional imaging VLC receiver'], ...
        npix_cols(ires), npix_rows(ires)));
end % for ires

% Plot the receivers against each other
plotmarkers = '.ox+*vs^p<h>'; % omit d to reserve for the DMD plot
num_pm = length(plotmarkers);

la = cell(2+nres,1);

h = figure(nres+2); % previous figures go to nres+1
clf
hold on
plot(snr_dB, mean_ub_ElMos_dmd, 'd-');
la{1} = 'DMD UB';
for ires = 1:nres
    plot(snr_dB, mean_ub_ElMos(ires,:), ...
        [plotmarkers(mod(ires-1,num_pm)+1), '-']);
    la{ires+1} = sprintf('%dx%d UB', npix_cols(ires), npix_rows(ires));
end % for ires

% reset the colors for the lower bound
set(gca, 'ColorOrderIndex', 1);

plot(snr_dB, mean_lb_dmd, 'd:');
la{nres+2} = 'DMD LB';
for ires = 1:nres
    plot(snr_dB, mean_lb(ires,:), ...
        [plotmarkers(mod(ires-1,num_pm)+1), ':']);
    la{ires+nres+2} = sprintf('%dx%d LB', ...
        npix_cols(ires), npix_rows(ires));
end % for ires

legend(la{:}, 'Location', 'NorthWest');
xlabel('20*log_{10}(x_{max}/\sigma_w) dB');
ylabel('I(x;y) bits');
grid on;
title('Comparison of average capacity bounds');