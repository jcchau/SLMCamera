function testComputeCapacityOnlyAmplitudeConShortA100(tc)
% testComputeCapacityOnlyAmplitudeConShortA100 tests
% computeCapacityOnlyAmplitudeConShort for A = 100.  
% Want to ensure that this method under test will at least work for a
% signal-to-noise ratio of up to 100.  

A = 100;

[C, poi, voi] = SmithCapacity.computeCapacityOnlyAmplitudeConShort(A);

n = length(poi);
tc.verifyLength(voi, n, 'length(voi)');

tc.verifyEqual(sum(voi), 1, 'AbsTol', 1e-12, ...
    'Total probability constraint violated.');

% check symmetry
if(mod(n,2) == 0)
    negside_i = n/2:-1:1;
    posside_i = n/2+1:2*n;
else
    negside_i = floor(n/2):-1:1;
    center_i = floor(n/2) + 1;
    posside_i = floor(n/2)+2:2*n;
    tc.verifyEqual(poi(center_i), 0, 'AbsTol', 1e-5, 'poi symmetry.');
end
tc.verifyEqual(poi(negside_i), -poi(posside_i), ...
    'AbsTol', 1e-5, 'poi symmetry.');
tc.verifyEqual(voi(negside_i), voi(posside_i), ...
    'AbsTol', 1e-6, 'voi symmetry.');

% check voi ascending
for ii = length(posside_i)-1;
    tc.verifyGreaterThan(voi(posside_i(ii+1)), voi(posside_i(ii)), ...
        'voi should increase with poi''s distance from 0.');
end

% check that poi(1) and poi(end) are -A and A respectively.
tc.verifyEqual(poi([1,end]), [-A; A], 'AbsTol', 1e-6, 'poi extremes');

% C
I_expected = SmithCapacity.I(poi, voi);
tc.verifyEqual(C, I_expected, 'C');

end

