classdef CameraArray < handle
    % CAMERAARRAY Represents a unified camera (with both the pixel array
    % and the point representing the lens)
    %   Provides methods to move the camera as a connected unit.
    %
    % This class is a Handle class so that it can be modified (moved and
    % re-oriented) by class methods without needing to create and replace
    % the entire CameraArray object every time.  
    
    properties
        % the distance from the array to the point representing the lens        
        lens_to_array_distance
        
        % the rectangular array object representing the pixel array
        pixel_array
    end
    
    methods
        function obj = CameraArray(lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, ...
                nrows, ncols, pixel_template)
            % CameraArray constructs a CameraArray object.
            %
            % OBJ = CameraArray(LENSPOINT, LENS_TO_ARRAY_DISTANCE, ...
            %	ZENITH_ANGLE, AZIMUTH, TILT, ...
            %	ELEMENT_WIDTH, ELEMENT_HEIGHT, ...
            %	NROWS, NCOLS, PIXEL_TEMPLATE)
            %
            % LENSPOINT is the point representing the lens in the camera
            %   system.
            % LENS_TO_ARRAY_DISTANCE is the distance between LENSPOINT and
            %   the plane of the pixel array
            % ZENITH_ANGLE, AZIMUTH, and TILT specify the orientation of
            %   the camera (about LENSPOINT).
            % ELEMENT_WIDTH and ELEMENT_HEIGHT specify the dimension of
            %   each element in the array.
            % NROWS and NCOLS specify the number of rows and columns in the
            %   array.  
            % PIXEL_TEMPLATE is analogous to POLYGON_TEMPLATE in
            %   RectangularArray and specifies the shape of the pixel in
            %   each element.  
            
            if(nargin > 0)
                obj.lens_to_array_distance = lens_to_array_distance;
                
                plane_axes = CameraArray.genRotatedAxes(zenith_angle, ...
                    azimuth, tilt);
                
                array_center = lenspoint - ...
                    lens_to_array_distance * plane_axes.normal();
                
                if(nargin >= 9) % if polygon_template is provided
                    obj.pixel_array = RectangularArray(...
                        array_center, plane_axes, ...
                        element_width, element_height, ...
                        nrows, ncols, pixel_template);
                else
                    obj.pixel_array = RectangularArray(...
                        array_center, plane_axes, ...
                        element_width, element_height, ...
                        nrows, ncols);
                end % if(nargin >= 9)
            end % if(nargin > 0)
        end % function CameraArray
        
        %% Property access methods
        function set.lens_to_array_distance(obj, newval)
            if(~isscalar(newval))
                error('Property lens_to_array_distance must be scalar.');
            end
            if(~isreal(newval) || newval<0)
                error('Property lens_to_array_distance must be a non-negative real number.');
            end
            
            obj.lens_to_array_distance = newval;
        end % function set.lens_to_array_distance
        
        function set.pixel_array(obj, newval)
            if(~isa(newval, 'RectangularArray'))
                error('Property pixel_array must be a RectangularArray object.');
            end
            
            obj.pixel_array = newval;
        end % function set.pixel_array
        
        %% Methods to move the CameraArray
        function obj = setOrientation(obj, zenith_angle, azimuth, tilt)
            % setOrientation changes the orientation of the CameraArray
            % object while maintaining the position of the lenspoint.  
            %
            % OBJ = setOrientation(OBJ, ZENITH_ANGLE, AZIMUTH, TILT)
            %
            % ZENITH_ANGLE, AZIMUTH, and TILT define the new orientation.
            
            % calculate the lenspoint around which to rotate the camera
            lenspoint = obj.lenspoint();
            
            % Update the plane axes and the array's centerpoint to the new
            % orientation.
            obj.pixel_array.plane_axes = CameraArray.genRotatedAxes(...
                zenith_angle, azimuth, tilt);
            obj.pixel_array.centerpoint = lenspoint - ...
                obj.lens_to_array_distance * ...
                obj.pixel_array.plane_axes.normal();
        end % function setOrientation
        
        function obj = setLensPoint(obj, lenspoint)
            % setLensPoint changes the location of the CameraArray object
            % while maintaining the orientation.
            %
            % OBJ = setLensPoint(OBJ, LENSPOINT)
            %
            % LENSPOINT is a 3-element row vector representing the 3D
            % coordinates of the new point where the lens is located.  
            
            obj.pixel_array.centerpoint = lenspoint - ...
                obj.lens_to_array_distance * ...
                obj.pixel_array.plane_axes.normal();
        end % function setLensPoint
        
        %% General methods
        
        function lp = lenspoint(obj)
            % lenspoint Returns a 3-column row vector representing the
            % coordinates of the lenspoint.  
            lp = obj.pixel_array.centerpoint + ...
                obj.lens_to_array_distance * ...
                obj.pixel_array.plane_axes.normal();
        end % function lenspoint
        
        polyout = preclipPolygon(obj, polyin)
        
        a_tr = calculateTransmitterAreaReceived(obj, transmitter_polygon);
        
    end % methods
    
    methods(Static)
        function oa = genRotatedAxes(zenith_angle, azimuth, tilt)
            [h1, h2, h3] = rotateTo(zenith_angle, azimuth, tilt, ...
                1, 0, 0);
            [v1, v2, v3] = rotateTo(zenith_angle, azimuth, tilt, ...
                0, 1, 0);
            oa = OrthogonalAxes([h1, h2, h3], [v1, v2, v3]);
        end % function genRotatedAxes
    end % methods(Static)
end

