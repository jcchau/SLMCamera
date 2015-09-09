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
            
            polyout = obj.preclipPolygon(polygon_orig);
            
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

