function runTrial(ii_s, srand_s, savedir)
%RUNTRIAL Summary of this function goes here
%   Detailed explanation goes here

%% load trial parameters
i_trial = sscanf(ii_s, '%d');
srand = sscanf(srand_s, '%d');
rng(srand);

%% common parameters
% copied from MadCom2016_20151110.m

r_aperture = 25.4e-3;
lenspoint = [0,0,0];
lens_to_array_distance = 49.9e-3;
zenith_angle = 0;
azimuth = 0;
tilt = 0;
array_width = 14.52e-3;
array_height = 8.16e-3;

% transmitter parameters
txheight = 2.5;
fov_width = array_width * txheight / lens_to_array_distance;
fov_height = array_height * txheight / lens_to_array_distance;
%tx_side = min(fovwidth, fovheight)/10;
    % 40.881 mm
tx_side = 40e-3;
tx_template = [ -tx_side/2, -tx_side/2, txheight; ...
    -tx_side/2, tx_side/2, txheight; ...
    tx_side/2, tx_side/2, txheight; ...
    tx_side/2, -tx_side/2, txheight];

nverticestx = size(tx_template, 1);

% The ratios of x_max to sigma_w in which to run evaluate capacity.
snr_dB = [-20, 0, 20] + 40;
snr = 10.^(snr_dB/20);

%% trial parameters
% specific to this trial

offsetx_a = fov_width .* (rand()-0.5);
offsety_a = fov_height .* (rand()-0.5);

offsetx_b = fov_width .* (rand()-0.5);
offsety_b = fov_height .* (rand()-0.5);

txa = Polygon(tx_template + ...
    [ repmat(offsetx_a,nverticestx,1), ...
    repmat(offsety_a,nverticestx,1), zeros(nverticestx,1) ]);
txb = Polygon(tx_template + ...
    [ repmat(offsetx_b,nverticestx,1), ...
    repmat(offsety_b,nverticestx,1), zeros(nverticestx,1) ]);

%% For the DMD prototype

fprintf('For the DMD...\n');

% DMD prototype parameters
nrows_dmd = 1080;
ncols_dmd = 1920;
element_width_dmd = array_width / ncols_dmd;
element_height_dmd = array_height / nrows_dmd;

ir_dmd = ImagingReceiver(r_aperture, ...
    lenspoint, lens_to_array_distance, zenith_angle, azimuth, tilt, ...
    element_width_dmd, element_height_dmd, ...
    nrows_dmd, ncols_dmd);

% calculate
pxgain_dmd_a = ir_dmd.calculateTransmitterToPixelGain(txa);
pxgain_dmd_b = ir_dmd.calculateTransmitterToPixelGain(txb);

% compute for the DMD receiver which pixels to which PD

pxdirection = nan(nrows_dmd, ncols_dmd);
pxdirection(pxgain_dmd_a > pxgain_dmd_b) = 0;
pxdirection(pxgain_dmd_a < pxgain_dmd_b) = 1;

% randomly assign the remaining pixels with 50% probability either way.
equalcells = pxgain_dmd_a == pxgain_dmd_b;
pxdirection(equalcells) = randi(2, sum(equalcells(:)), 1) -1;

% the S selection matrix
% mirror is directed towards the first photodetector iff pxdirection is
% 0; and is directed towards the second photodetector otherwise.
S_dmd = [ ~pxdirection(:)'; pxdirection(:)' ];

% The resulting channel matrix, where G=SH and H is the gain from each
% transmitter through each SLM pixel.
G_dmd = S_dmd * [ pxgain_dmd_a(:), pxgain_dmd_b(:) ];

% metrics
fprintf('Computing metrics for the DMD...\n');

trial_rank_dmd = rank(G_dmd);
trial_cn_dmd = cond(G_dmd);

trial_ub_ElMos_dmd = zeros(1, length(snr));
trial_lb_MRC_dmd = zeros(1, length(snr));
trial_lb_UnifX_dmd = zeros(1, length(snr));

G_dmd = MIMOCapacity.removeZeroRowsAndCols(G_dmd);

for isnr = 1:length(snr)
    trial_ub_ElMos_dmd(isnr) = ...
        MIMOCapacity.computeCapacityUBElMoslimany2014(G_dmd, snr(isnr), 1);
    trial_lb_MRC_dmd(isnr) = ...
        MIMOCapacity.computeCapacityLBMRC(G_dmd, snr(isnr), 1);
    trial_lb_UnifX_dmd(isnr) = ...
        MIMOCapacity.computeCapacityLBUnifX(G_dmd, snr(isnr), 1);
end % for isnr

%% For the traditional imaging VLC receiver

npix_rows = [ 1, 2, 3, 10, 100, 400, 600, 720, 1080 ];
npix_cols = [ 2, 2, 2, 10, 100, 600, 800, 1280, 1920 ];
npix = npix_rows .* npix_cols;

trial_rank = zeros(length(npix), 1);
trial_cn = zeros(length(npix), 1);
trial_ub_ElMos = zeros(length(npix), length(snr));
trial_lb_MRC = zeros(length(npix), length(snr));
trial_lb_UnifX = zeros(length(npix), length(snr));

for ipix = 1:length(npix)
    
    fprintf('For %dx%d...\n', npix_cols(ipix), npix_rows(ipix));

    % npix-pixel imaging receiver
    nrows = npix_rows(ipix);
    ncols = npix_cols(ipix);
    element_width = array_width / ncols;
    element_height = array_height / nrows;

    ir = ImagingReceiver(r_aperture, ...
        lenspoint, lens_to_array_distance, zenith_angle, azimuth, tilt, ...
        element_width, element_height, ...
        nrows, ncols);

    % calculate pixel gains
    pxgain_a = ir.calculateTransmitterToPixelGain(txa);
    pxgain_b = ir.calculateTransmitterToPixelGain(txb);
    
    % compute the final channel matrices
    H = [ pxgain_a(:), pxgain_b(:) ];
    
    trial_rank(ipix) = rank(H);
    trial_cn(ipix) = cond(H);
    
    % remove zero rows and columns for capacity calculations
    H = MIMOCapacity.removeZeroRowsAndCols(H);
    nrows_nonzero = size(H,1);
    fprintf('nrows_nonzero=%d\n', nrows_nonzero);
    
    for isnr = 1:length(snr)
        fprintf('Computing metrics for %dx%d, snr=%ddB...\n', ...
            npix_cols(ipix), npix_rows(ipix), snr_dB(isnr));
        trial_ub_ElMos(ipix,isnr) = ...
            MIMOCapacity.computeCapacityUBElMoslimany2014(H, snr(isnr), 1);
        trial_lb_MRC(ipix,isnr) = ...
            MIMOCapacity.computeCapacityLBMRC(H, snr(isnr), 1);
        if(nrows_nonzero <= 12)
            trial_lb_UnifX(ipix,isnr) = ...
                MIMOCapacity.computeCapacityLBUnifX(H, snr(isnr), 1);
        else
            % Too many dimensions.  Assuming 1e9/8/5 bins, we'd average
            % approximately 4 bins per dimension if we had 12 dimensions.  
            trial_lb_UnifX(ipix,isnr) = NaN;
        end
    end % for isnr
    
end % for ipix

%% Save the results

fprintf('Saving...\n');

fname = sprintf('%s/trial%d.mat', savedir, i_trial);

% Warning: Variable 'S_dmd' cannot be saved to a MAT-file whose version is
% older than 7.3.
% To save this variable, use the -v7.3 switch.
save(fname, '-v7.3', 'snr_dB', 'npix', ...
    'trial_rank_dmd', 'trial_cn_dmd', ...
    'trial_ub_ElMos_dmd', 'trial_lb_MRC_dmd', 'trial_lb_UnifX_dmd', ...
    'trial_rank', 'trial_cn', ...
    'trial_ub_ElMos', 'trial_lb_MRC', 'trial_lb_UnifX')

fprintf('Done.\n');

end

