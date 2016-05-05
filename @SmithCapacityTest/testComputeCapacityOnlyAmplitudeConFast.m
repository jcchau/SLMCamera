function testcomputeCapacityOnlyAmplitudeConFast(tc)
% testcomputeCapacityOnlyAmplitudeConFast ensures that
% computeCapacityOnlyAmplitudeConFast returns the same result as
% computeCapacityOnlyAmplitudeCon, but runs faster.  

A = 6 * rand();

%% Methods under test, where the long version gives the "expected" result.

tic
[C1, poi1, voi1] = SmithCapacity.computeCapacityOnlyAmplitudeCon(A);
t1 = toc;

tic
[C2, poi2, voi2] = SmithCapacity.computeCapacityOnlyAmplitudeConFast(A);
t2 = toc;

%% Verify

tc.verifyEqual(C2, C1, 'AbsTol', 1e-10, 'C');
tc.verifyEqual(poi2, poi1, 'AbsTol', 1e-5, 'poi');
tc.verifyEqual(voi2, voi1, 'AbsTol', 1e-5, 'voi');

% Allow for a 0.5 second tolerance.
tc.verifyLessThanOrEqual(t2, t1, 'time');

end

