classdef ImagingReceiverTest < matlab.unittest.TestCase
    %IMAGINGRECEIVERTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Test)
        %% Constructor
        function testConstructorHappy(tc)
            [r_aperture, lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template] = ...
                ImagingReceiverTest.genImagingReceiverParams();
            
            obj = ImagingReceiver(r_aperture, ...
                lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template);
            
            ImagingReceiverTest.verifyProperties(tc, obj, ...
                r_aperture, lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template);
        end % function testConstructorHappy

        function testConstructorNoArgument(~)
            ImagingReceiver();
        end % function testConstructorNoArgument
        
        function testConstructorDefaultPixelTemplate(tc)
            [r_aperture, lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                ~] = ...
                ImagingReceiverTest.genImagingReceiverParams();
            
            obj = ImagingReceiver(r_aperture, ...
                lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols);
            
            pixel_template_expected = ...
                CameraArrayTest.defaultPixelTemplate( ...
                element_width, element_height);
            
            ImagingReceiverTest.verifyProperties(tc, obj, ...
                r_aperture, lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template_expected);
        end % function testConstructorDefaultPixelTemplate
        
        function testConstructorNoSuperconstructorArgs(tc)
            % testConstructorNoSuperconstructorArgs tries to call the
            % constructor with no arguments for the superconstructor
            % (CameraArray)
            
            r_aperture = rand()/rand();
            
            tc.verifyError(@() ImagingReceiver(r_aperture), ...
                'ImagingReceiver:ImagingReceiver:notEnoughInputs', ...
                'Constructor should not accept just one argument.');
        end % function testConstructorNoSuperconstructorArgs
        
        %% calculateTransmitterToPixelGain
        testCalculateTransmitterToPixelGainOverAllPixels(tc)
        
        
    end % methods(Test)
    
    methods(Static)
        function verifyProperties(tc, obj, ...
                r_aperture, lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template)
            
            tc.verifyEqual(obj.r_aperture, r_aperture);
            
            CameraArrayTest.verifyProperties(tc, obj, ...
                lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template);
        end % function verifyProperties
        
        function [r_aperture, lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template] = genImagingReceiverParams()
            
            r_aperture = rand()/rand();
            
            [lenspoint, lens_to_array_distance, ...
                zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols, ...
                pixel_template] = ...
                CameraArrayTest.genCameraArrayParameters();
        end % function genImagingReceiverParams()
    end % methods(Static)
end

