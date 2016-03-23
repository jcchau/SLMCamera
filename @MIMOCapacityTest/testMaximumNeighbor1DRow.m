function testMaximumNeighbor1DRow(tc)
% testMaximumNeighbor1DRow verifies that MIMOCapacity.maximumNeighbor using
% a random 1D (row vector) input.

matlen = randi(5);

in = rand(1, matlen);

out = MIMOCapacity.maximumNeighbor(in);

%% compute the expected result

% start out with each element as its maximum
expected_out = in;

% Then replace each "maximum" with the next element if it's bigger
expected_out(1:end-1) = max(expected_out(1:end-1), in(2:end));

% Or its previous element if it's bigger.
expected_out(2:end) = max(expected_out(2:end), in(1:end-1));

%% check

tc.verifyEqual(out, expected_out, ...
    'minimumNeighbor did not return the expected output.');

end

