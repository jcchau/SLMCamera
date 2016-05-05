function testComputeTableOfCapacity(tc)

Amax = 6.*rand();
npoints = randi(10)+1;

A = linspace(Amax/npoints, Amax, npoints);

[C, n] = SmithCapacity.computeTableOfCapacity(A);

tc.verifyLength(C, npoints);
tc.verifyLength(n, npoints);

tc.verifyGreaterThanOrEqual(n, 2, 'n should always be >= 2.');

for ii = 2:npoints
    tc.verifyGreaterThan(C(ii), C(ii-1), 'C should be ascending.');
    tc.verifyGreaterThanOrEqual(n(ii), n(ii-1), 'n should not decrease.');
end % for ii

end

