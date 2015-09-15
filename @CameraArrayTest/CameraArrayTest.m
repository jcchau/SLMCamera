classdef CameraArrayTest < matlab.unittest.TestCase
    % CAMERAARRAYTEST Unit tests for CameraArray.
    
    properties
    end
    
    methods(Test)
        %% constructor
        function testConstructorHappy(tc)
            %% generate parameters
            [lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template] = ...
                CameraArrayTest.genCameraArrayParameters();
            
            %% construct the object
            obj = CameraArray(lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, ...
                nrows, ncols, pixel_template);
            
            %% verify the properties
            
            tc.verifyEqual(obj.lens_to_array_distance, ...
                lens_to_array_distance);
            
            [localZ(1,1), localZ(1,2), localZ(1,3)] = ...
                rotateTo(zenith_angle, azimuth, tilt, 0, 0, 1);
            [localX(1,1), localX(1,2), localX(1,3)] = ...
                rotateTo(zenith_angle, azimuth, tilt, 1, 0, 0);
            [localY(1,1), localY(1,2), localY(1,3)] = ...
                rotateTo(zenith_angle, azimuth, tilt, 0, 1, 0);
            
            centerpoint = lenspoint - lens_to_array_distance * localZ;
            tc.verifyLessThan(norm(obj.pixel_array.centerpoint - ...
                centerpoint), 1e-10);
            
            tc.verifyLessThan(norm(...
                obj.pixel_array.plane_axes.horizontal - localX), 1e-10);
            tc.verifyLessThan(...
                norm(obj.pixel_array.plane_axes.vertical - localY), 1e-10);
            
            tc.verifyEqual(obj.pixel_array.element_width, element_width);
            tc.verifyEqual(obj.pixel_array.element_height, element_height);
            
            tc.verifyEqual(obj.pixel_array.nrows, nrows);
            tc.verifyEqual(obj.pixel_array.ncols, ncols);
            
            tc.verifyEqual(obj.pixel_array.polygon_template.toMatrix(), ...
                pixel_template.toMatrix());
        end % testConstructorHappy
        
        function testConstructorDefaultPixelTemplate(tc)
            % testConstructorDefaultPolygonTemplate verifies that the
            % constructor still works when the optional pixel_template
            % argument is not provided.
            
            [lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ~] = ...
                CameraArrayTest.genCameraArrayParameters();
            
            obj = CameraArray(lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, ...
                nrows, ncols);
            
            %% verify the properties
            % same tests as in testConstructorHappy
            
            tc.verifyEqual(obj.lens_to_array_distance, ...
                lens_to_array_distance);
            
            [localZ(1,1), localZ(1,2), localZ(1,3)] = ...
                rotateTo(zenith_angle, azimuth, tilt, 0, 0, 1);
            [localX(1,1), localX(1,2), localX(1,3)] = ...
                rotateTo(zenith_angle, azimuth, tilt, 1, 0, 0);
            [localY(1,1), localY(1,2), localY(1,3)] = ...
                rotateTo(zenith_angle, azimuth, tilt, 0, 1, 0);
            
            centerpoint = lenspoint - lens_to_array_distance * localZ;
            tc.verifyLessThan(norm(obj.pixel_array.centerpoint - ...
                centerpoint), 1e-10);
            
            tc.verifyLessThan(norm(...
                obj.pixel_array.plane_axes.horizontal - localX), 1e-10);
            tc.verifyLessThan(...
                norm(obj.pixel_array.plane_axes.vertical - localY), 1e-10);
            
            tc.verifyEqual(obj.pixel_array.element_width, element_width);
            tc.verifyEqual(obj.pixel_array.element_height, element_height);
            
            tc.verifyEqual(obj.pixel_array.nrows, nrows);
            tc.verifyEqual(obj.pixel_array.ncols, ncols);
            
            % check that the pixel_template fills the pixel
            pixel_template_matrix_expected = ...
                [0, 0, 0; ...
                0, element_height, 0; ...
                element_width, element_height, 0; ...
                element_width, 0, 0];
            tc.verifyEqual(obj.pixel_array.polygon_template.toMatrix(), ...
                pixel_template_matrix_expected);
            
        end % function testConstructorDefaultPixelTemplate
        
        %% property access methods
        function testSetLensToArrayDistanceHappy(tc)
            [lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template] = ...
                CameraArrayTest.genCameraArrayParameters();
            
            obj = CameraArray(lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, ...
                nrows, ncols, pixel_template);
            
            new_lens_to_array_distance = 0.1*rand();
            
            obj.lens_to_array_distance = new_lens_to_array_distance;
            
            tc.verifyEqual(obj.lens_to_array_distance, ...
                new_lens_to_array_distance);
        end % function testSet_lens_to_array_distance
        
        function testSetLensToArrayDistanceVector(tc)
            [lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template] = ...
                CameraArrayTest.genCameraArrayParameters();
            
            obj = CameraArray(lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, ...
                nrows, ncols, pixel_template);
            
            new_lens_to_array_distance = 0.1*rand([2,1]);
            
            try
                obj.lens_to_array_distance = new_lens_to_array_distance;
                tc.verifyFail(...
                    'The set property access method should have generated an error.');
            catch
            end
        end % function testSetLensToArrayDistanceVector
        
        %% preclipPolygon
        function testPreclipPolygonSimple(tc)
            [~, lens_to_array_distance, ...
                ~, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template] = ...
                CameraArrayTest.genCameraArrayParameters();
            
            % set the lenspoint and zenith_angle so that I can personally
            % calculate the expected result.
            lenspoint = [0,0,0];
            zenith_angle = 0;
            
            obj = CameraArray(lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, ...
                nrows, ncols, pixel_template);
            
            % generate the polygon that crosses the xy-plane
            % intersects the xy-plane at [2.5, +/-1, 0]
            polygon_orig_matrix = [2, -1, 0.5; ...
                2, 1, 0.5; ...
                3, 1, -0.5; ...
                3, -1, -0.5];
            
            polygon_orig = Polygon(polygon_orig_matrix);
            
            poly_plane = polygon_orig.getPlane();
            
            [polyout, isvalid] = obj.preclipPolygon(polygon_orig);
            
            tc.verifyTrue(isvalid, ...
                'The clipped polygon should be valid in this case.');
            
            polyout_matrix = polyout.toMatrix();
            
            nvectors_out = size(polyout_matrix, 1);
            
            for ii = 1:nvectors_out
                tc.verifyGreaterThan(polyout_matrix(ii,3), 0, ...
                    'The clipped polygon should not have vertices below the xy-plane.');
            end % for ii
            
            if(nvectors_out > 4)
                % in case an extra vertex was added to close the polygon
                tc.verifyEqual(polyout_matrix(1,:), polyout_matrix(end,:));
                
                % remove the extra vertex
                polyout_matrix = polyout_matrix(1:end-1, :);
                nvectors_out = size(polyout_matrix, 1);
            end % if(nvectors_out > 4)
            
            tc.verifyEqual(nvectors_out, 4, ...
                'The clipped polygon should only have 4 vertices');
            
            % Ensure that we still have the two vertices that should be
            % kept
            tc.verifyTrue(...
                any(ismember(polyout_matrix, [2, -1, 0.5], 'rows')), ...
                'preclipPolygon should have kept vertex [2, -1, 0.5].');
            tc.verifyTrue(...
                any(ismember(polyout_matrix, [2, 1, 0.5], 'rows')), ...
                'preclipPolygon should have kept vertex [2, 1, 0.5].');
            
            % TODO: more tests that the polygons are (approximately) eqaul
            % are needed
            
        end % function testPreclipPolygonSimple
        
        function testPreclipPolygonDiscardAll(tc)
            % testPreclipPolygonDiscardAll checks that the returned isvalid
            % is false when the entire polygon is discarded.
            
            [~, lens_to_array_distance, ...
                ~, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template] = ...
                CameraArrayTest.genCameraArrayParameters();
            
            % set the lenspoint and zenith_angle so that I can personally
            % calculate the expected result.
            lenspoint = [0,0,0];
            zenith_angle = 0;
            
            obj = CameraArray(lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, ...
                nrows, ncols, pixel_template);
            
            % create a polygon with all vertices on the side of the plane
            % to be removed in clipping.
            polyin_matrix = [1, 0, -1; 1 1, -1; 2, 0, -1];
            polyin = Polygon(polyin_matrix);
            
            [~, isvalid] = obj.preclipPolygon(polyin);
            
            tc.verifyFalse(isvalid, 'isvalid should be false.');            
        end % function testPreclipPolygonDiscardAll
        
        %% calculateTransmitterAreaReceived
        
        function testCalculateTransmitterAreaReceivedWholeTx(tc)
            % testCalculateTransmitterAreaReceivedWholeTx verifies that if
            % the whole transmitter is in the field-of-view of the receiver
            % with 100% fill factor, then the total area of the transmitter
            % received equals the area of the transmitter polygon.  
            %
            % Since the method calcualteTransmitterAreaReceived is rather
            % complex, I use this approach to test the method (rather than
            % to compare the test result a_tr against a known-correct
            % a_tr) for now.
            
            %% create the CameraArray object ca
            % Aim a CameraArray with a wide FoV in the [1,1,1] direction to
            % cover most of the first octant.  
            lenspoint = [-1, -1, -1];
            lens_to_array_distance = 1e-2;  % 1cm
            
            % vector [1,1,1] points like the hypotenuse of a right triangle
            % with a sqrt(2) base and a height of 1.  The sqrt(2) base is
            % due to the distance of 1 along both x and y directions.
            zenith_angle = pi/2 - atan(1/sqrt(2));
            % so that projected onto the x-y plane, the camera points 45
            % degrees between the +x and +y axis.
            azimuth = 3*pi/4;
            tilt = 0;
            
            % The array_width and array_height determine the angle of view.
            % Making them twice the lens_to_array_distance ensures a 45
            % degree angle-of-view in both the horizontal and vertical
            % direction.
            array_width = 2*lens_to_array_distance;
            array_height = 2*lens_to_array_distance;
            
            % randomly choose the number of rows and columns; this should
            % not significantly affect the result.
            nrows = randi([1,3000]);
            ncols = randi([1,3000]);
            
            % and then scale the pixel sizes to fill the array size
            element_width = array_width / ncols;
            element_height = array_height / nrows;
            
            ca = CameraArray(lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, ...
                nrows, ncols);
            
            % check that the camera array is pointed the correct direction
            expected_array_normal = [1/sqrt(3), 1/sqrt(3), 1/sqrt(3)];
            test_array_normal = ca.pixel_array.plane_axes.normal();
            array_normal_error = ...
                norm(test_array_normal - expected_array_normal);
            tc.verifyLessThan(array_normal_error, 1e-10, ...
                'The CameraArray is not aimed as expected.');
            
            %% put a transmitter polygon entirely in the FoV
            % a rectangle
            % Make sure it's oriented so that the camera sees at least some
            % of the face (and not just see the edge)
            poly_tx = Polygon(...
                [1, 1, 1; ...
                1, 1.2, 1; ...
                1.2, 1.2, 1.1; ...
                1.2, 1, 1.1]);
            
            % a rectangle with one side 0.2 long and the other side is the
            % hypotenuse of a right triangle with other sides 0.2 and 0.1
            % long.
            area_tx_expected = 0.2 * sqrt(0.2^2 + 0.1^2);
            
            % check that poly_tx.area() matches
            area_tx_test = poly_tx.area();
            tc.verifyEqual(area_tx_test, area_tx_expected, ...
                'RelTol', 1e-10, ...
                'The poly_tx.area() did not return the expected area.');
            
            %% use calculateTransmitterAreaReceived to get a_tr and check
            a_tr = ca.calculateTransmitterAreaReceived(poly_tx);
            
            % check that a_tr is the expected size
            size_a_tr = size(a_tr);
            tc.verifyEqual(size_a_tr, [nrows, ncols]);
            
            % check the total area
            total_area_test = sum(sum(a_tr));
            tc.verifyEqual(total_area_test, area_tx_expected, ...
                'RelTol', 1e-10, ...
                'The total area in a_tr does not equal the area of the transmitter polygon.');
        end % function testCalculateTransmitterAreaReceivedWholeTx
    end % methods(Test)
    
    methods(Static)
        function [lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template] = genCameraArrayParameters()
            % genCameraArrayParameters generates a set of parameters for
            % constructing a CameraArray object.
            
            lenspoint = [ 10*(rand(1,2)-0.5), 1.5*rand() ];
            lens_to_array_distance = 0.1*rand();
            zenith_angle = pi * rand();
            azimuth = 2*pi * rand();
            tilt = 2*pi * rand();
            element_width = 1e-5 * rand();
            element_height = 1e-5 * rand();
            nrows = randi(2000);
            ncols = randi(2000);
            normalized_pixel_template = [...
                0, 0; ...
                0, 1; ...
                0.7, 0.8; ...
                0.9, 0];
            pixel_template = Polygon(...
                element_width * normalized_pixel_template(:,1), ...
                element_height * normalized_pixel_template(:,2), ...
                zeros(size(normalized_pixel_template,1), 1));
        end % function genCameraArrayParameters
    end % methods(Static)
    
end

