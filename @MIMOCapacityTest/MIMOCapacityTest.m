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
        
        testCalculateDiffEntropyFromPmf1DUniform(tc)
        
        testMinimumNeighborEmptyIn(tc)
        testMinimumNeighbor1DRow(tc)
        testMinimumNeighbor1DCol(tc)
        testMinimumNeighborSingletonDimension(tc)
        testMinimumNeighbor2DManual(tc)
        testMinimumNeighbor2DRandom(tc)
        testMinimumNeighbor3D(tc)
        
        testMaximumNeighborEmptyIn(tc)
        testMaximumNeighbor1DCol(tc)
        testMaximumNeighbor1DRow(tc)
        testMaximumNeighbor2DManual(tc)
        testMaximumNeighbor2DRandom(tc)
        testMaximumNeighbor3D(tc)
    end % methods(Test)
    
end

