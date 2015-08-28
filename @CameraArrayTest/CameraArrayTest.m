classdef CameraArrayTest < matlab.unittest.TestCase
    % CAMERAARRAYTEST Unit tests for CameraArray.
    
    properties
    end
    
    methods(Test)
        %% constructor
        function testConstructorHappy(tc)
            %% generate parameters
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
    end % methods(Test)
    
end

