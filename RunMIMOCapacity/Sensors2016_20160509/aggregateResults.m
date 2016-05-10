function [trial_rank_dmd, trial_cn_dmd, ...
    trial_ub_ElMos_dmd, trial_lb_MRC_dmd, trial_lb_UnifX_dmd, ...
    trial_lb_UnifGx_dmd, ...
    trial_rank, trial_cn, ...
    trial_ub_ElMos, trial_lb_MRC, trial_lb_UnifX, ...
    trial_lb_UnifGx] = ...
    aggregateResults(ntrials)

S = load('data/trial1.mat', 'snr_dB', 'npix');
nsnr = length(S.snr_dB);
nres = length(S.npix);
clear S

%% preallocate space

% DMD-based
trial_rank_dmd = zeros(ntrials, 1);
trial_cn_dmd = zeros(ntrials, 1);
trial_ub_ElMos_dmd = zeros(ntrials, nsnr);
trial_lb_MRC_dmd = zeros(ntrials, nsnr);
trial_lb_UnifX_dmd = zeros(ntrials, nsnr);

% traditional imaging
trial_rank = zeros(ntrials, nres);
trial_cn = zeros(ntrials, nres);
trial_ub_ElMos = zeros(ntrials, nres, nsnr);
trial_lb_MRC = zeros(ntrials, nres, nsnr);
trial_lb_UnifX = zeros(ntrials, nres, nsnr);

%% load data
for ii = 1:ntrials
    mat_fname = sprintf('data/trial%d.mat', ii);
    S = load(mat_fname, 'trial_rank_dmd', 'trial_cn_dmd', ...
        'trial_ub_ElMos_dmd', 'trial_lb_MRC_dmd', 'trial_lb_UnifX_dmd', ...
        'trial_rank', 'trial_cn', ...
        'trial_ub_ElMos', 'trial_lb_MRC', 'trial_lb_UnifX', ...
        'trial_lb_UnifGx_dmd', 'trial_lb_UnifGx');
    
    trial_rank_dmd(ii) = S.trial_rank_dmd;
    trial_cn_dmd(ii) = S.trial_cn_dmd;
    trial_ub_ElMos_dmd(ii,:) = reshape(S.trial_ub_ElMos_dmd, 1, nsnr);
    trial_lb_MRC_dmd(ii,:) = reshape(S.trial_lb_MRC_dmd, 1, nsnr);
    trial_lb_UnifX_dmd(ii,:) = reshape(S.trial_lb_UnifX_dmd, 1, nsnr);
    trial_lb_UnifGx_dmd(ii,:) = reshape(S.trial_lb_UnifGx_dmd, 1, nsnr);
    
    trial_rank(ii,:) = reshape(S.trial_rank, 1, nres);
    trial_cn(ii,:) = reshape(S.trial_cn, 1, nres);
    trial_ub_ElMos(ii,:,:) = reshape(S.trial_ub_ElMos, 1, nres, nsnr);
    trial_lb_MRC(ii,:,:) = reshape(S.trial_lb_MRC, 1, nres, nsnr);
    trial_lb_UnifX(ii,:,:) = reshape(S.trial_lb_UnifX, 1, nres, nsnr);
    trial_lb_UnifGx(ii,:,:) = reshape(S.trial_lb_UnifGx, 1, nres, nsnr);
end % for ii

end % function