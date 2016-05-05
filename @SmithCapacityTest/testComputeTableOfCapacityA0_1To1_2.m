function testComputeTableOfCapacityA0_1To1_2(tc)
%TESTCOMPUTETABLEOFCAPACITYA0_1TO6 checks for A=0.1:0.1:1.2
%
% I test this particular case because I got this error before:
% n=3 fmincon... checking...
% Error using SmithCapacity.computeCapacityOnlyAmplitudeConFast (line 173)
% n=3 is far too large for A=0.9.
% 
% Error in SmithCapacity.computeTableOfCapacity (line 26)
%     [C(ii), poi, ~] = ...


A = 0.1:0.1:1.2;
npoints = length(A);

[C, n] = SmithCapacity.computeTableOfCapacity(A);

tc.verifyLength(C, npoints);
tc.verifyLength(n, npoints);

tc.verifyGreaterThanOrEqual(n, 2, 'n should always be >= 2.');

for ii = 2:npoints
    tc.verifyGreaterThan(C(ii), C(ii-1), 'C should be ascending.');
    tc.verifyGreaterThanOrEqual(n(ii), n(ii-1), 'n should not decrease.');
end % for ii

end

