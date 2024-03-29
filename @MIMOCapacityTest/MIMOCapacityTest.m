classdef MIMOCapacityTest < matlab.unittest.TestCase
    % MIMOCapacityTest A test class for MIMOCapacity.
    
    properties
    end
    
    methods(Test)
        testConvertToLinearIndex(tc)
        
        testCalculateDiffEntropyFromPmf1DUniform(tc)
        testCalculateDiffEntropyFromPmf2DGaussian(tc)
        
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
        
        testCalculateMinimumDiffEntropyFromPmfConvergesUniform(tc)
        testCalculateMinimumDiffEntropyFromPmfConvergesGaussianEven(tc)
        testCalculateMinimumDiffEntropyFromPmfConvergesGaussianOdd(tc)
        testCalculateMinimumDiffEntropyFromPmf2D(tc)
        
        combinedTestCalculateMinAndMaxDiffEntropy3DConverges(tc)
        
        testCalculateDiffEntropyOfGaussian(tc)
        
        testConvertPointToSubscriptIndex(tc)
        
        testConvertLinearToSubscriptIndexInverse(tc)
        testConvertLinearToSubscriptIndexInd2sub(tc)
        
        testRemoveZeroRowsAndCols(tc)
        testRemoveZeroRowsAndColsAllZeros(tc)
        
        testCalculateDiffEntropyOfMVGaussian(tc)
        
        % Test simplifyChannelMatrix and computeQTTransform
        testComputeQTTransformNoDependence(tc)
        testComputeQTTransformDependentColMoreTx(tc)
        testComputeQTTransformDependentRx(tc)
        testSimplifyChannelMatrixStepA(tc)
        testSimplifyChannelMatrixStepB(tc)
        
        testComputeUExtremesEmpirically(tc)
        
        testComputeCapacityLBUnifXFullRankHighSNR(tc)
        testComputeCapacityLBUnifXFullRankLowSNR(tc)
        
        testComputeCapacityLBUnifGxAgainstUnifX(tc)
        
    end % methods(Test)
    
end

