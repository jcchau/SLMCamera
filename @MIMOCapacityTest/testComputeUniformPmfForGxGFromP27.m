function testComputeUniformPmfForGxGFromP27(tc)
% testComputeUniformPmfForGxGFromP27 tests
% MIMOCapacity.computeUniformPmfForGx using the channel matrix G from
% p.26-27 of lab book #4 and x_max = 1.  

G = [1, 1, 0.5;
    0.5, 1, 1];

[n_r, n_t] = size(G);
x_max = ones(n_t, 1);
y_min = zeros(1, n_r);
nbins = [5, 5];
y_max = [2.5, 2.5];
delta = (y_max - y_min) ./ nbins;

[pmf, reachable] = MIMOCapacity.computeUniformPmfForGx(G, x_max, ...
    y_min, delta, nbins);

reachable_expected = ...
    [ 1 1 1 0 0;
    1 1 1 1 0;
    1 1 1 1 1;
    0 1 1 1 1;
    0 0 1 1 1];
pmf_expected = zeros(nbins);
pmf_expected(reachable_expected==1) = 1/sum(reachable_expected(:));

tc.verifyEqual(reachable, (reachable_expected==1), 'reachable');
tc.verifyEqual(pmf, pmf_expected, 'pmf');

end

