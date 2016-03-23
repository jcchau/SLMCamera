function testMinimumNeighbor2DManual(tc)
% testMinimumNeighbor2DManual checks that MIMOCapacity.minimumNeighbor
% works correctly for a 2D example.
% Expected result was verified manually.

A = [ ...
    0.8147    0.0975    0.1576    0.1419    0.6557    0.7577; ...
    0.9058    0.2785    0.9706    0.4218    0.0357    0.7431; ...
    0.1270    0.5469    0.9572    0.9157    0.8491    0.3922; ...
    0.9134    0.9575    0.4854    0.7922    0.9340    0.6555; ...
    0.6324    0.9649    0.8003    0.9595    0.6787    0.1712];

out = MIMOCapacity.minimumNeighbor(A, 2);

expected_out = [ ...
    0         0         0         0         0         0; ...
    0    0.0975    0.0975    0.0357    0.0357         0; ...
    0    0.1270    0.2785    0.0357    0.0357         0; ...
    0    0.1270    0.4854    0.4854    0.1712         0; ...
    0         0         0         0         0         0];

tc.verifyEqual(out, expected_out, ...
    'The output did not match the manually-verified expected output.');

end

