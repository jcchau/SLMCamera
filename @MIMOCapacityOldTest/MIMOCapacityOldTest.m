classdef MIMOCapacityOldTest < matlab.unittest.TestCase
    % MIMOCapacityOldTest are tests belonging to MIMOCapacityOld.
    
    properties
    end
    
    methods(Test)
        
        testCalculateCapacityForUniformInputVarianceH(tc)
        
        testGenerateReceivedPmfForUniformInput1D(tc)
        testGenerateReceivedPmfForUniformInput1DNoSignal(tc)
        testGenerateReceivedPmfForUniformInput1DNoNoise(tc)
        testGenerateReceivedPmfForUniformInput2D(tc)
        testGenerateReceivedPmfForUniformInput2DCorrelatedGNoNoise(tc)
        testGenerateReceivedPmfForUniformInput3D(tc)
        
    end
    
end

