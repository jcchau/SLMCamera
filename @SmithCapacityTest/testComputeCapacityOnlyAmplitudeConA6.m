function testComputeCapacityOnlyAmplitudeConA6(tc)
% testComputeCapacityOnlyAmplitudeConALE6 tests
% computeCapacityOnlyAmplitudeCon for A = 6.0.  

% According to Smith 1969:
% - 6 points of increase,
% - poi and voi symmetric about 0.
% - voi grows as the point-of-interest gets farther from 0.  

A = 6.0;

[C, poi, voi] = SmithCapacity.computeCapacityOnlyAmplitudeCon(A);

tc.verifyLength(poi, 6, 'length(poi)');
tc.verifyLength(voi, 6, 'length(voi)');

tc.verifyEqual(sum(voi), 1, 'AbsTol', 1e-12, ...
    'Total probability constraint violated.');

% check symmetry
tc.verifyEqual(poi(3:-1:1), -poi(4:6), 'AbsTol', 1e-6, 'poi symmetry.');
tc.verifyEqual(voi(3:-1:1), voi(4:6), 'AbsTol', 1e-6, 'voi symmetry.');

% check voi ascending
for ii=0:1
    tc.verifyGreaterThan(voi(5+ii), voi(4+ii), ...
        'voi should increase with poi''s distance from 0.');
end

% check that poi(1) and poi(end) are -A and A respectively.
tc.verifyEqual(poi([1,end]), [-A; A], 'AbsTol', 1e-6, 'poi extremes');

% C
I_expected = SmithCapacity.I(poi, voi);
tc.verifyEqual(C, I_expected, 'C');

end

