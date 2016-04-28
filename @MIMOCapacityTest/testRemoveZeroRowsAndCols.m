function testRemoveZeroRowsAndCols(tc)
% testRemoveZeroRowsAndCols checks MIMOCapacity.removeZerosRowsAndCols
% using a pre-defined matrix.

G = rand(10);
zero_rows = [1, 3, 4, 5, 9];
zero_cols = [1, 2, 3, 7];

G(zero_rows, :) = 0;
G(:, zero_cols) = 0;

test_out = MIMOCapacity.removeZeroRowsAndCols(G);

nonzero_rows = [2, 6, 7, 8, 10];
nonzero_cols = [4, 5, 6, 8, 9, 10];
expected_out = G(nonzero_rows, nonzero_cols);

tc.verifyEqual(test_out, expected_out, ...
    'test_out does not match expected_out');

end

