classdef MIMOCapacityLBUnifXTest < matlab.unittest.TestCase
    % MIMOCapacityLBUnifXTest are tests belonging to MIMOCapacityLBUnifX.
    
    properties
    end
    
    methods(Test)
        
        testCalculateCapacityLBForUniformInputVarianceH(tc)
        
        testGenerateReceivedPmfForUniformInput1D(tc)
        testGenerateReceivedPmfForUniformInput1DNoSignal(tc)
        testGenerateReceivedPmfForUniformInput1DNoNoise(tc)
        testGenerateReceivedPmfForUniformInput2D(tc)
        testGenerateReceivedPmfForUniformInput2DCorrelatedGNoNoise(tc)
        testGenerateReceivedPmfForUniformInput3D(tc)
        
    end
    
end

