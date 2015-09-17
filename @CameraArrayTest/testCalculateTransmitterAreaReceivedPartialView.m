function testCalculateTransmitterAreaReceivedPartialView(tc)
    % testCalculateTransmitterAreaReceivedPartialView verifies that
    % the total area received as calculated by
    % calculateTransmitterAreaReceived equals the area of the
    % transmitter within the field of view when the transmitter
    % polygon is only partially in the field of view.  

    % Design the system for 1/100 magnification, placing the
    % transmitter at z=1.
    % Have the camera see everything within -0.5<=x<=0.5,
    % -0.5<=y<=0.5, z=1. (Nothing more or less on the z=1 plane.)

    lenspoint = [0,0,0];
    lens_to_array_distance = 1e-2;
    zenith_angle = 0;
    azimuth = 0;
    tilt = 0;

    array_side_length = lens_to_array_distance;
    nrows = randi([1,100]);
    ncols = randi([1,100]);
    element_width = array_side_length/ncols;
    element_height = array_side_length/nrows;

    ca = CameraArray(lenspoint, lens_to_array_distance, ...
        zenith_angle, azimuth, tilt, ...
        element_width, element_height, nrows, ncols);

    % place a transmitter polygon on z=1 and randomly shift it out
    % of the FoV.
    % Don't make this transmitter polygon square to better expose
    % the potential problems of using
    % RectangularArray.listIntersectingElements (which may return
    % some candidates that don't really intersect with the
    % specified polygon).
    poly_matrix_tx = [-0.5, -0.5, 1; ...
        -0.5*rand(), 0.5*rand(), 1; ...
        0.5, 0.5, 1; ...
        0.5*rand(), -0.5*rand(), 1];
    xshift = 2*(rand()-0.5);
    yshift = 2*(rand()-0.5);
    poly_matrix_tx(:,1) = poly_matrix_tx(:,1) + xshift;
    poly_matrix_tx(:,2) = poly_matrix_tx(:,2) + yshift;
    poly_tx = Polygon(poly_matrix_tx);

    %% calculate the expected tx area seen by the camera

    % create a 2D polygon defining the FoV at z=1.
    poly_2D_fov = [-0.5, -0.5; ...
        -0.5, 0.5; ...
        0.5, 0.5; ...
        0.5, -0.5];
    [x_in_fov, y_in_fov] = polybool('intersection', ...
        poly_matrix_tx(:,1), poly_matrix_tx(:,2), ...
        poly_2D_fov(:,1), poly_2D_fov(:,2));
    area_tx_in_fov = polyarea(x_in_fov, y_in_fov);

    %% get and check a_tr
    a_tr = ca.calculateTransmitterAreaReceived(poly_tx);

    tc.verifyEqual(size(a_tr), [nrows, ncols], ...
        'a_tr does not have the expected matrix size.');

    total_area_test = sum(sum(a_tr));

    tc.verifyEqual(total_area_test, area_tx_in_fov, ...
        'AbsTol', 1e-12, ...
        'The transmitter area received calculated by the calculateTransmitterAreaReceived does not match the transmitter area in the field of view.');

end % function testCalculateTransmitterAreaReceivedPartialView