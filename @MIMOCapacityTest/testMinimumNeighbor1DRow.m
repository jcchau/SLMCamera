function testMinimumNeighbor1DRow(tc)
% testMinimumNeighbor1DRow verifies that MIMOCapacity.minimumNeighbor using
% a random 1D (row vector) input.

matlen = randi(5);

in = rand(1, matlen);

out = MIMOCapacity.minimumNeighbor(in, 1);

%% compute the expected result
expected_out = zeros(size(in));

expected_out(2:end-1) = min(in(1:end-2), in(3:end));
expected_out(2:end-1) = min(expected_out(2:end-1), in(2:end-1));

%% check

tc.verifyEqual(out, expected_out, ...
    'minimumNeighbor did not return the expected output.');

end

