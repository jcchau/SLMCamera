classdef ImagingReceiver < CameraArray
    % ImagingReceiver represents an imaging VLC receiver. 
    
    properties
        % The lens/aperture radius.
        r_aperture
    end
    
    methods
        function obj = ImagingReceiver(r_aperture, varargin)
            % ImagingReceiver constructs an ImagingReceiver object.
            %
            % OBJ = ImagingReceiver(R_APERTURE, ...
            %   LENSPOINT, LENS_TO_ARRAY_DISTANCE, ...
            %	ZENITH_ANGLE, AZIMUTH, TILT, ...
            %	ELEMENT_WIDTH, ELEMENT_HEIGHT, ...
            %	NROWS, NCOLS, PIXEL_TEMPLATE)
            %
            % R_APERTURE is the radius of the lens/aperture.  It is used to
            %   calculate how much light is gathered by the receiver.
            % LENSPOINT, LENS_TO_ARRAY_DISTANCE, ZENITH_ANGLE, AZIMUTH,
            %   TILT, ELEMENT_WIDTH, ELEMENT_HEIGHT, NROWS, NCOLS, and
            %   PIXEL_TEMPLATE are the arguments for the CameraArray
            %   constructor.  
            %
            % See also CameraArray.CameraArray.
            
            if(nargin == 0)
                args_camera_array = {};
            elseif(nargin <= 1)
                % Verify that we're not invoking CameraArray's no-argument
                % constructor.  
                % We don't want to invoke the no-argument constructor for
                % superclass CameraArray unless we're using the no-argument
                % class for this class too.
                error('ImagingReceiver:ImagingReceiver:notEnoughInputs', ...
                    'Not enough input arguments.');
            else
                args_camera_array = varargin;
            end
            obj@CameraArray(args_camera_array{:});
            
            if(nargin > 0)
                obj.r_aperture = r_aperture;
            end
        end % function ImagingReceiver
        
        gains = calculateTransmitterToPixelGain(obj, transmitter_polygon)
    end % methods
    
    methods(Static)
        diffent = calculateDiffEntropyOfOutputForUniformInput( ...
            S, SH, x_max, variance_bg, variance_t)
    end % methods(Static)
end
