function testFillUniformPmfForGxAgainstComputeUniformPmfForGx(tc)
% Compares the result of fillUniformPmfForGx against the result of
% computeUniformPmfForGx.

n_r = randi(6);
n_t = randi(6);
G = rand(n_r, n_t);
G = MIMOCapacity.simplifyChannelMatrix(G);
Q = MIMOCapacity.computeQTTransform(G);

% Apply Q' transform, and use this new G.
G = Q' * G;
[n_r, n_t] = size(G);

xmax = 1./rand(n_t,1);

% Pick a value so that computeUniformPmfForGx does not take too long, but
% want this to be a large value.
max_nbins = 6e4;

% nbins is a col vector
nbins = (MIMOCapacity.fillMaxNBins(max_nbins, n_r))';

% umin, umax are col vectors
[umin, umax] = MIMOCapacity.computeUExtremes(G, xmax);

% Pick ymin and ymax to cover a slightly larger range than umin and umax.
% (Do this so we test the case where not everything is aligned with a bin
% boundary: a more typical case.)
ymin = umin - rand(n_r,1) .* (umax-umin)/100;
ymax = umax + rand(n_r,1) .* (umax-umin)/100;

delta = (ymax-ymin)./nbins;

%% Method under test

pmf = MIMOCapacityLBUnifGx.fillUniformPmfForGx( ...
    G, xmax, ymin', delta', nbins');

tc.verifyEqual(sum(pmf(:)), 1, 'AbsTol', 1e-12);

mode_pmf = mode(pmf(pmf>0));
tc.verifyTrue(all(pmf(pmf>0)==mode_pmf), ...
    'Expect all pmf(pmf>0) to have the same value.');

%% computeUniformPmfForGx for comparison

pmf_expected = MIMOCapacityLBUnifGx.computeUniformPmfForGx( ...
    G, xmax, ymin', delta', nbins');

tc.verifyEqual(sum(pmf_expected(:)), 1, 'AbsTol', 1e-12);

mode_pmf_expected = mode(pmf_expected(pmf_expected>0));
tc.verifyTrue(all(pmf_expected(pmf_expected>0)==mode_pmf_expected), ...
    'Expect all pmf_expected(pmf_expected>0) to have the same value.');

%% Verify the result through comparison

pmf_differs = xor(pmf>0, pmf_expected>0);

tc.verifyEqual(pmf, pmf_expected, 'AbsTol', 1e-14, ...
    sprintf('Of %d bins in the PMF, %d differ.', ...
    prod(nbins), nnz(pmf_differs)));

end

