function testRemoveZeroRowsAndColsAllZeros(tc)
% testRemoveZeroRowsAndColsAllZeros checks
% MIMOCapacity.removeZerosRowsAndCols using a matrix of zeros.

G = zeros(randi(10), randi(10));

test_out = MIMOCapacity.removeZeroRowsAndCols(G);

expected_out = 0;

tc.verifyEqual(test_out, expected_out, ...
    'test_out does not match expected_out');

end

