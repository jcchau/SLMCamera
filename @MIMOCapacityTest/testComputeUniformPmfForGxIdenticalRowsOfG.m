function testComputeUniformPmfForGxIdenticalRowsOfG(tc)
% TESTCOMPUTEUNIFORMPMFFORGXIDENTICALROWSOFG checks that
% computeUniformPmfForGx outputs the correct results when rows of G are
% identical. 

G = [ 0.5, 0.5; 0.5, 0.5 ];
x_max = [1; 1];
y_min = [-.2, -.2];
delta = [0.1, 0.1];
nbins = [14, 14];

[pmf, reachable] = MIMOCapacity.computeUniformPmfForGx(G, x_max, ...
    y_min, delta, nbins);

reachable_expected = false(nbins);
for ii = 2:13
    reachable_expected(ii,ii) = true;
end
for ii = 3:13
    % and the bins that touch by a corner
    reachable_expected(ii-1, ii) = true;
    reachable_expected(ii, ii-1) = true;
end

tc.verifyEqual(reachable, reachable_expected, 'reachable');

pmf_expected = zeros(nbins);
pmf_expected(reachable_expected) = 1/nnz(reachable_expected);

tc.verifyEqual(pmf, pmf_expected, 'pmf');

end

