function testComputeCapacityOnlyAmplitudeConFastA32(tc)
% testcomputeCapacityOnlyAmplitudeConFastA32 tests
% computeCapacityOnlyAmplitudeConFast for A = 32.  
% Want to ensure that this method under test will at least work for a
% signal-to-noise ratio of up to 32 (a little larger than 30dB =
% 10^(30/20)).  
%
% To avoid taking excessively long on weaker computers, only go up to A=12
% unless we have at least 12 parpool workers.  
%
% This was originally attempted for A=100.  However, even on the SCC with
% 12 MATLAB workers we only made it to n=87 before reaching the time limit.
% And then, resuming at n=87, we only got to n=94 before timing out again.
% Finally, starting at n=100, we reached n=114 before timing out.

pp = gcp();
if(isempty(pp) || pp.NumWorkers < 12)
    % Getting to A=100 takes a very long time.  Only try if we're running
    % with a lot of parpool workers.  
    A = 12;
    nStart = 12; % skip to the n=12 to be fast
    tol_poi = 1e-3;
    tol_voi = 1e-5;
    
    [C, poi, voi] = SmithCapacity.computeCapacityOnlyAmplitudeConFast( ...
        A, nStart);
else
    A = 32;
    tol_poi = 1e-2;
    tol_voi = 1e-3;
    
    % Don't use nStart to be thorough.
    [C, poi, voi] = SmithCapacity.computeCapacityOnlyAmplitudeConFast(A);
end

% Sort points of increase in ascending order
[poi, sortindex] = sort(poi);
voi = voi(sortindex); % does not change shape

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
    posside_i = floor(n/2)+2:n;
    tc.verifyEqual(poi(center_i), 0, 'AbsTol', tol_poi, ...
        'poi (center) symmetry.');
end
tc.verifyEqual(poi(negside_i), -poi(posside_i), ...
    'AbsTol', tol_poi, 'poi symmetry.');
tc.verifyEqual(voi(negside_i), voi(posside_i), ...
    'AbsTol', tol_voi, 'voi symmetry.');

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

