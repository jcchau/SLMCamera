classdef MIMOCapacityLBUnifGxTest < matlab.unittest.TestCase
    %MIMOCAPACITYLBUNIFGXTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Test)
        testComputeReceivedPmfViaUnifThenConv2D(tc)
        testComputeReceivedPmfViaUnifThenConv2DExpandedSigmaW(tc)
        
        testComputeUniformPmfForGx3DGEye(tc)
        testComputeUniformPmfForGxGFromP27(tc)
        testComputeUniformPmfForGxIdenticalRowsOfG(tc)
        
        testGenerateUniformPmfForGxGFromP27(tc)
        testFillUniformPmfForGxGFromP27(tc)
        
        testCalculateMutualInfoWithUnifGxUnitaryT(tc)
        
        testFillUniformPmfForGxAgainstComputeUniformPmfForGx(tc)
    end
    
end

