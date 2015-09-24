function testCalculateTransmitterToPixelGainOverAllPixels(tc)
% testCalculateTransmitterToPixelGainOverAllPixels verifies that, when the
% transmitter is entirely within the ImagingReceiver's field of view, that
% the total gain across all ImagingReceiver pixels (approximately) equals
% the calculated optical gain from the transmitter to the aperature. 

poly_tx_x = 10e-3 .* [-1; -1; 1; 1];
poly_tx_y = 10e-3 .* [-1; 1; 1; -1];
poly_tx_z = [1.996; 1.996; 2.004; 2.0004]; % slight upward tilt
poly_tx = Polygon(poly_tx_x, poly_tx_y, poly_tx_z);

lenspoint = [0,0,0];
% So l = 2 and A_tr*cos(theta) = (20e-3)^2
l = 2;
a_tr_cos_theta = (20e-3)^2;

r_aperture = 25e-3;
% so A_r is
a_r = pi * (r_aperture)^2;

% make the camera have a 30 degree (half-angle) field of view.
lens_to_array_distance = 1e-2;
array_side_length = 2*lens_to_array_distance * tan(pi/6);

% Allow up to 5 degree zenith angle (whole transmitter should
% still remain in FoV) 
zenith_angle = 5 * pi/180 * rand();
% And use random azimuth and tilt
azimuth = 2*pi*rand();
tilt = 2*pi*rand();

% With this setup, expect a magnification of
magnification_approx = lens_to_array_distance ./ 2;
% where 2 is the distance to the transmitter.

% Let's have approximately 20 pixels wide cover the transmitter
pixels_per_mm = 20 / (magnification_approx * 20e-3);
% where 20e-3 is the width of the transmitter polygon.

nrows = ceil(array_side_length * pixels_per_mm);
ncols = nrows;

element_width = array_side_length / ncols;
element_height = array_side_length / nrows;

% calculate expected total gain (where A_tr == A_t)
total_gain_expected = a_tr_cos_theta/poly_tx.area() * ...
    a_r / pi * cos(zenith_angle) / l^2;

% make the ImagingReceiver object.
ir = ImagingReceiver(r_aperture, lenspoint, ...
    lens_to_array_distance, zenith_angle, azimuth, tilt, ...
    element_width, element_height, nrows, ncols);

gains = ir.calculateTransmitterToPixelGain(poly_tx);

tc.verifyEqual(size(gains), [nrows, ncols], ...
    'The returned gains matrix was not the correct size.');
tc.verifyGreaterThanOrEqual(gains, 0, ...
    'None of the gains should be negative.');

total_gain_test = sum(gains(:));

% NOTE: the approximately 6% error we get doesn't instill too much
% confidence.  Do more tests. 
tc.verifyEqual(total_gain_test, total_gain_expected, ...
    'RelTol', 0.06, ...
    'The test total gain doesn''t match the expected total gain.');
end % function testCalculateTransmitterToPixelGainOverAllPixels