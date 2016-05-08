function testComputeUniformPmfForGx3DGEye(tc)
% testComputeUniformPmfForGx3DGEye tests computeUniformPmfForGx using G =
% eye(3).

G = eye(3);
x_max = 10;
delta = repmat(1, 1, 3);
y_min = zeros(1,3) - 2*delta;
nbins = x_max' ./ delta + 4;

[pmf, reachable] = MIMOCapacityLBUnifGx.computeUniformPmfForGx( ...
    G, x_max, y_min, delta, nbins);

reachable_expected = false(nbins);
reachable_expected(2:end-1, 2:end-1, 2:end-1) = true;

tc.verifyEqual(reachable, reachable_expected, 'reachable');

pmf_expected = zeros(nbins);
pmf_expected(reachable_expected) = 1/nnz(reachable_expected);

tc.verifyEqual(pmf, pmf_expected, 'pmf');

end

