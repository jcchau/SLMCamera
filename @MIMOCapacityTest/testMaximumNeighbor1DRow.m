function testMaximumNeighbor1DRow(tc)
% testMaximumNeighbor1DRow verifies that MIMOCapacity.maximumNeighbor using
% a random 1D (row vector) input.

matlen = randi(5);

in = rand(1, matlen);

out = MIMOCapacity.maximumNeighbor(in);

%% compute the expected result
expected_out = zeros(size(in));

% Check each element against the next
expected_out(1:end-1) = max(in(1:end-1), in(2:end));

% And then check each maximum against the previous element
expected_out(2:end) = max(expected_out(2:end), in(1:end-1));

%% check

tc.verifyEqual(out, expected_out, ...
    'minimumNeighbor did not return the expected output.');

end

