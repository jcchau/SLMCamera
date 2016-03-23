function testMinimumNeighbor1DCol(tc)
% testMinimumNeighbor1D verifies that MIMOCapacity.minimumNeighbor using a
% random 1D (column vector) input.

matlen = randi(5);

in = rand(matlen, 1);

out = MIMOCapacity.minimumNeighbor(in, 1);

%% compute the expected result
expected_out = zeros(size(in));

expected_out(2:end-1) = min(in(1:end-2), in(3:end));
expected_out(2:end-1) = min(expected_out(2:end-1), in(2:end-1));

%% check

tc.verifyEqual(out, expected_out, ...
    'minimumNeighbor did not return the expected output.');

end

