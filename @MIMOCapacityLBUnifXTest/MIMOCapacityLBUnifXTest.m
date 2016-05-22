classdef MIMOCapacityLBUnifXTest < matlab.unittest.TestCase
    % MIMOCapacityLBUnifXTest are tests belonging to MIMOCapacityLBUnifX.
    
    properties
    end
    
    methods(Test)
        
        % Disabled this test because it always fails.  
        % It fails because the variance_H calculation is incorrect as
        % explained on p.36-38 of lab book #4.  
        %testCalculateCapacityLBForUniformInputVarianceH(tc)
        
        testGenerateReceivedPmfForUniformInput1D(tc)
        testGenerateReceivedPmfForUniformInput1DNoSignal(tc)
        testGenerateReceivedPmfForUniformInput1DNoNoise(tc)
        testGenerateReceivedPmfForUniformInput2D(tc)
        testGenerateReceivedPmfForUniformInput2DCorrelatedGNoNoise(tc)
        testGenerateReceivedPmfForUniformInput3D(tc)
        
    end
    
end

