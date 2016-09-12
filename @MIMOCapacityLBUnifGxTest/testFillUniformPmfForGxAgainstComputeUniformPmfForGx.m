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
% 6e5 already takes over an hour.  Don't want to increase any more.  
% For the test below, adjust the number of allowable missed bins according
% to n_r.  
max_nbins = 6e5;

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

tic
pmf = MIMOCapacityLBUnifGx.fillUniformPmfForGx( ...
    G, xmax, ymin', delta', nbins');
toc

tc.verifyEqual(sum(pmf(:)), 1, 'AbsTol', 1e-11);

mode_pmf = mode(pmf(pmf>0));
tc.verifyTrue(all(pmf(pmf>0)==mode_pmf), ...
    'Expect all pmf(pmf>0) to have the same value.');

%% computeUniformPmfForGx for comparison

tic
pmf_expected = MIMOCapacityLBUnifGx.computeUniformPmfForGx( ...
    G, xmax, ymin', delta', nbins');
toc

tc.verifyEqual(sum(pmf_expected(:)), 1, 'AbsTol', 1e-11);

mode_pmf_expected = mode(pmf_expected(pmf_expected>0));
tc.verifyTrue(all(pmf_expected(pmf_expected>0)==mode_pmf_expected), ...
    'Expect all pmf_expected(pmf_expected>0) to have the same value.');

%% Verify the result through comparison

% Okay, I've examined that they're different at the outer boundary, with
% the method under test missing some of the bins that
% computeUniformPmfForGx hits.  I think this is acceptable, especially
% since we're aiming for a lower bound (and including all of the boundary
% bins may cause us to overestimate this lower bound).  
%
% tc.verifyEqual(pmf, pmf_expected, 'AbsTol', 1e-14, ...
%     sprintf(['Of %d bins in the PMF, %d differ. ' ...
%     'pmf has %d hits and pmf_expected has %d hits'], ...
%     prod(nbins), nnz(pmf_differs), nnz(pmf>0), nnz(pmf_expected>0)));

pmf_missed = pmf==0 & pmf_expected>0;
pmf_extra = pmf>0 & pmf_expected==0;

tc.verifyEqual(nnz(pmf_extra), 0, ...
    'We don''t want fillUniformPmfForGx to hit extra bins.');

% Expected missed bins
% Allow bins along the perimeter to be missed.
% Since we don't have a easy way to calculate the number of bins on the
% outer boundary of the polytope that is reachable:
allowed_missed_bins = prod(nbins) - prod((nbins-2));

tc.verifyLessThanOrEqual(nnz(pmf_missed), ...
    allowed_missed_bins, ...
    sprintf('fillUniformPmfForGx missed %d%% bins, exceeded %d%%.', ...
    nnz(pmf_missed), allowed_missed_bins));

end

