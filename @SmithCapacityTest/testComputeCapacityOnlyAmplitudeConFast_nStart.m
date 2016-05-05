function testComputeCapacityOnlyAmplitudeConFast_nStart(tc)
% testComputeCapacityOnlyAmplitudeConFast_nStart checks the nStart feature
% of computeCapacityOnlyAmplitudeConFast. 

A = 6 * rand();

if(A>=3)
    nStart = 4;
else
    nStart = 2;
end

%% Methods under test, where the long version gives the "expected" result.

[C1, poi1, voi1] = SmithCapacity.computeCapacityOnlyAmplitudeConFast(A);
[C2, poi2, voi2] = SmithCapacity.computeCapacityOnlyAmplitudeConFast( ...
    A, nStart);

%% Verify

tc.verifyEqual(C2, C1, 'AbsTol', 1e-10, 'C');
tc.verifyEqual(poi2, poi1, 'AbsTol', 1e-5, 'poi');
tc.verifyEqual(voi2, voi1, 'AbsTol', 1e-5, 'voi');

end

