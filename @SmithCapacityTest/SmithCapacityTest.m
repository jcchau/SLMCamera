classdef SmithCapacityTest < matlab.unittest.TestCase
    %SMITHCAPACITYTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Test)
        function testCheckCorollary1ForA1_1(tc)
            % Verifies that checkCorollary1 returns true for A=1.1.
            % According to Smith1969 p.51, n=2 is optimal for A<=1.6.
            
            A=1.1;
            poi = [-A, A];
            voi = [0.5, 0.5];
            I_Fo = SmithCapacity.I([-A, A], [0.5, 0.5]);
            
            %x = linspace(poi(1), poi(end), 100);
            %figure();
            %plot(x, SmithCapacity.i(x, poi, voi))
            %hold on
            %plot([x(1), x(end)], repmat(I_Fo, 1, 2), 'r')
            
            r = SmithCapacity.checkCorollary1(A, poi, voi, I_Fo);
            tc.verifyTrue(r, 'poi and voi should be optimal for A=1.1.');
        end % function testCheckCorollary1ForA1_1
        
        function testCheckCorollary1ForA1_6(tc)
            % Verifies that checkCorollary1 returns true for A=1.6.
            % According to Smith1969 p.51, n=2 is optimal for A<=1.6.
            
            A=1.6;
            poi = [-A, A];
            voi = [0.5, 0.5];
            I_Fo = SmithCapacity.I([-A, A], [0.5, 0.5]);
            
            %x = linspace(poi(1), poi(end), 100);
            %figure();
            %plot(x, SmithCapacity.i(x, poi, voi))
            %hold on
            %plot([x(1), x(end)], repmat(I_Fo, 1, 2), 'r')
            
            r = SmithCapacity.checkCorollary1(A, poi, voi, I_Fo);
            tc.verifyTrue(r, 'poi and voi should be optimal for A=1.6.');
        end % function testCheckCorollary1ForA1_6
        
        function testCheckCorollary1ForA1_7(tc)
            % Verifies that checkCorollary1 returns true for A=1.7.
            % According to Smith1969 p.51, n=2 is not optimal for A>=1.7.
            
            A=1.7;
            poi = [-A, A];
            voi = [0.5, 0.5];
            I_Fo = SmithCapacity.I([-A, A], [0.5, 0.5]);
            
            %x = linspace(poi(1), poi(end), 100);
            %figure();
            %plot(x, SmithCapacity.i(x, poi, voi))
            %hold on
            %plot([x(1), x(end)], repmat(I_Fo, 1, 2), 'r')
            
            r = SmithCapacity.checkCorollary1(A, poi, voi, I_Fo);
            tc.verifyFalse(r, 'n=2 should not be optimal for A=1.7.');
        end % function testCheckCorollary1ForA1_7
        
        testComputeCapacityOnlyAmplitudeConALE1_6(tc)
        testComputeCapacityOnlyAmplitudeConA1_7(tc)
        testComputeCapacityOnlyAmplitudeConA6(tc)
        
        testComputeCapacityOnlyAmplitudeConShort(tc)
        testComputeCapacityOnlyAmplitudeConShortA100(tc)
        
        testSmithCapacity_i(tc)
        
        testSmithCapacity_p_Y(tc)
    end % methods(Test)
    
end

