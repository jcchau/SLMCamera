function testComputeCapacityOnlyAmplitudeConA1_7(tc)
% testComputeCapacityOnlyAmplitudeConALE1_7 tests
% computeCapacityOnlyAmplitudeCon for A = 1.7, where the optimal
% distribution of X should have 3 points of increase at -A, 0, and A, where
% the points of increase at A have equal probability.  

A = 1.7;

[C, poi, voi] = SmithCapacity.computeCapacityOnlyAmplitudeCon(A);

tc.verifyEqual(poi, [-A; 0; A], 'AbsTol', 1e-6, 'poi');

tc.verifyLength(voi, 3, 'length(voi)');
tc.verifyEqual(voi(1), voi(3), 'AbsTol', 1e-6, 'voi(1) ~= voi(3)');
tc.verifyEqual(sum(voi), 1, 'AbsTol', 1e-12, ...
    'Total probability constraint violated.');

% C should be greater than if n=2
I_n2 = SmithCapacity.I([-A, A], [0.5, 0.5]);
tc.verifyGreaterThan(C, I_n2, 'C');

end

