function testGenerateUniformPmfForGxGFromP27(tc)
% testComputeUniformPmfForGxGFromP27 tests
% MIMOCapacity.computeUniformPmfForGx using the channel matrix G from
% p.26-27 of lab book #4 and x_max = 1.  

G = [1, 1, 0.5;
    0.5, 1, 1];

[n_r, n_t] = size(G);
x_max = ones(n_t, 1);
y_min = zeros(1, n_r);
nbins = [5, 5];
y_max = (G*x_max)';
delta = (y_max - y_min) ./ nbins;

pmf = MIMOCapacity.generateUniformPmfForGx(G, x_max, y_min, delta, nbins);

pmf_expected = [ ...
    1 1 0 0 0;
    1 1 1 0 0;
    0 1 1 1 0;
    0 0 1 1 1;
    0 0 0 1 1];
pmf_expected(pmf_expected==1) = 1/nnz(pmf_expected);

tc.verifyEqual(pmf, pmf_expected, 'pmf');

end

