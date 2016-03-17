classdef MIMOCapacityTest < matlab.unittest.TestCase
    % MIMOCapacityTest A test class for MIMOCapacity.
    
    properties
    end
    
    methods(Test)
        testConvertToLinearIndex(tc)
        
        testGenerateReceivedPmfForUniformInput1D(tc)
        testGenerateReceivedPmfForUniformInput1DNoSignal(tc)
        testGenerateReceivedPmfForUniformInput1DNoNoise(tc)
        testGenerateReceivedPmfForUniformInput2D(tc)
        testGenerateReceivedPmfForUniformInput2DCorrelatedGNoNoise(tc)
        testGenerateReceivedPmfForUniformInput3D(tc)
    end % methods(Test)
    
end

