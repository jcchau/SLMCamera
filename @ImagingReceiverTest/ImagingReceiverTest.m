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
        
        function testCalculateTransmitterToPixelGainPerfectAlignment(tc)
            % In this test case, the transmitter is only seen by one pixel
            % on the receiver such that theta and phi are exactly 0.  This
            % allows us to precisely compute the expected gain for
            % verification.  
            %
            % The only error (difference) should be due to
            % precision/rounding errors.  
            
            tx_height = 2;
            tx_side = 10e-3;
            poly_tx = Polygon( ...
                [-tx_side/2, -tx_side/2, tx_height; ...
                -tx_side/2, tx_side/2, tx_height; ...
                tx_side/2, tx_side/2, tx_height; ...
                tx_side/2, -tx_side/2, tx_height]);
            
            r_aperture = 12.7e-3;   % 1 inch diameter
            lenspoint = [0,0,0];
            lens_to_array_distance = 20e-3;
            
            % set element_width and element_height to equal the width of
            % the transmitter image
            element_width = tx_side * lens_to_array_distance / tx_height;
            element_height = element_width;
            
            % surround the central pixel with pixels that should have zero
            % gain so we can check that too.
            nrows = 3;
            ncols = 3;
            
            zenith_angle = 0;
            azimuth = 0;
            tilt = 0;
            
            % calculate the expected gain: A_r / pi / l^2
            gain_expected = r_aperture^2 / tx_height^2;
            
            ir = ImagingReceiver(r_aperture, lenspoint, ...
                lens_to_array_distance, zenith_angle, azimuth, tilt, ...
                element_width, element_height, nrows, ncols);
            
            gains = ir.calculateTransmitterToPixelGain(poly_tx);
            
            tc.verifyEqual(size(gains), [nrows, ncols], ...
                'The returned matrix gains is not the correct size.');
            
            tc.verifyGreaterThanOrEqual(gains, 0, ...
                'None of the gains should be negative.');
            tc.verifyLessThanOrEqual(gains, 1, ...
                'None of the gains greater than unity.');
            
            tc.verifyEqual(gains(2,2), gain_expected, 'RelTol', 1e-10);
            
            % Check that aside from cell (2,2), all gains are zero
            matrix_expected = zeros(nrows, ncols);
            matrix_expected(2,2) = gains(2,2);
            tc.verifyEqual(gains, matrix_expected, ...
                'The gains for all other elements should be exactly zero.')
        end % function testCalculateTransmitterToPixelGainPerfectAlignment
        
        function testConvertToLinearIndex(tc)
            % Compares the result of convertToLinearIndex against sub2ind
            % with a random-dimension random-size matrix, and a random
            % number of random indices to convert.
            % Simultaneously tests method convertToLinearIndexWeights.
            
            % up to 4 dimensions
            num_dims = randi(4);
            
            % each dimension up to 11 long
            mat_size = randi(11, 1, num_dims);
            
            % up to 5 indices
            num_indices = randi(5);
            
            % the indices
            subs_cell = cell(1, num_dims);
            subs_mat = zeros(num_indices, num_dims);
            for ii = 1:num_dims
                subs_cell{ii} = randi(mat_size(ii), num_indices, 1);
                subs_mat(:, ii) = subs_cell{ii};
            end % for ii = 1:dimensions
            
            weights = ImagingReceiver.convertToLinearIndexWeights( ...
                mat_size);
            li_test = ImagingReceiver.convertToLinearIndex( ...
                weights, subs_mat);
            
            if(num_dims == 1)
                % sub2ind doesn't work for 1D vectors; in this case, the
                % matrix subscript indexing is equal to the linear
                % indexing.
                li_expected = subs_cell{1};
            else
                li_expected = sub2ind(mat_size, subs_cell{:});
            end
            
            tc.verifyTrue(isequal(li_test, li_expected), ...
                ['convertToLinearIndex did not return the same ' ...
                'results as sub2ind.'])
        end
        
        testGenerateReceivedPmfForUniformInput1D(tc)
        testGenerateReceivedPmfForUniformInput1DNoSignal(tc)
        testGenerateReceivedPmfForUniformInput1DNoNoise(tc)
        testGenerateReceivedPmfForUniformInput2D(tc)
        testGenerateReceivedPmfForUniformInput2DCorrelatedGNoNoise(tc)
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

