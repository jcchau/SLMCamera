function testComputeCapacityOnlyAmplitudeConALE1_6(tc)
% testComputeCapacityOnlyAmplitudeConALE1_6 test
% computeCapacityOnlyAmplitudeCon for A <= 1.6, where the optimal
% distribution of X should have 2 points of increase at -A and A of equal
% probability.  

A = 1.6 * rand();

[C, poi, voi] = SmithCapacity.computeCapacityOnlyAmplitudeCon(A);

tc.verifyEqual(poi, [-A; A], 'AbsTol', 1e-6, 'poi');
tc.verifyEqual(voi, [0.5; 0.5], 'AbsTol', 1e-6, 'voi');

C_expected = SmithCapacity.I([-A, A], [0.5, 0.5]);

tc.verifyEqual(C, C_expected, 'AbsTol', 1e-6, 'C');

end

