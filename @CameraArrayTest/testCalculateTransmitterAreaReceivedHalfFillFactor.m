function testCalculateTransmitterAreaReceivedHalfFillFactor(tc)
    % testCalculateTransmitterAreaReceivedHalfFillFactor verifies
    % that a high-resolution camera with 50% fill factor sees 50%
    % of the area of a transmitter that's entirely within its field
    % of view.  

    lenspoint = [0,0,-1];

    % poly_tx is located directly over lenspoint (2m above).
    % A camera looking directly up with at least a 10.1 degree
    % (half-angle) field of view should be able to see the whole
    % thing. 
    poly_tx = Polygon( ...
        [-0.25, -0.25, 1; ...
        -0.25, 0.25, 1; ...
        0.25, 0.25, 1; ...
        0.25, -0.25, 1]);

    % make the camera have a 30 degree (half-angle) field of view.
    lens_to_array_distance = 1e-2;
    array_side_length = 2*lens_to_array_distance * tan(pi/6);

    % Allow up to 5 degree zenith angle (whole transmitter should
    % still remain in FoV) 
    zenith_angle = 5 * pi/180 * rand();
    % And use random azimuth and tilt
    azimuth = 2*pi*rand();
    tilt = 2*pi*rand();

    % With this setup, expect a magnification of of 1e-2/2 = 0.005.
    % So out poly_tx with 0.5m side length forms an image with
    % 2.5mm side length.  
    % Let's say we want at least 200 pixels in 2.5mm (so the
    % percent of the transmitter area seen is determined by fill
    % factor and is less affected by pixel alignment).
    linear_res_on_tx = 200;
    min_linear_resolution = ceil(array_side_length * ...
        linear_res_on_tx / 2.5e-3);
    nrows = randi([min_linear_resolution, ...
        ceil(1.4*min_linear_resolution)]);
    ncols = randi([min_linear_resolution, ...
        ceil(1.4*min_linear_resolution)]);
    % Chose ceil(1.4*min_linear_resolution) as the max linear
    % resolution because min_linear_resolution is large and we
    % don't want to more than double the already long minimum
    % compute time.

    element_width = array_side_length/ncols;
    element_height = array_side_length/nrows;

    % pixel_template with 50% fill factor
    pixel_template = Polygon( ...
        [0,0,0; ...
        0, element_height, 0; ...
        element_width, 0, 0]);

    % the CameraArray
    ca = CameraArray(lenspoint, lens_to_array_distance, ...
        zenith_angle, azimuth, tilt, ...
        element_width, element_height, nrows, ncols, ...
        pixel_template);

    % get the output of the method under test
    a_tr = ca.calculateTransmitterAreaReceived(poly_tx);

    area_seen_test = sum(sum(a_tr));
    area_seen_expected = poly_tx.area()/2;

    acceptable_reltol = 1-(1-1/linear_res_on_tx)^2;
    tc.verifyEqual(area_seen_test, area_seen_expected, ...
        'RelTol', acceptable_reltol, ...
        'The camera should see half of the transmitter polygon.');
end % function testCalculateTransmitterAreaReceivedHalfFillFactor