function testComputeReceivedPmfViaUnifThenConv2D(tc)
% testComputeReceivedPmfViaUnifThenConv2D tests
% MIMOCapacity.computeReceivedPmfViaUnifThenConv with a basic 2D example.

G = [ 1, 0; 0, 1 ];
x_max = 2;
sigma_w = [1; 1];
ns = 5;
nbins = [12; 12];

%% run the method under test

[pmf, delta] = MIMOCapacity.computeReceivedPmfViaUnifThenConv( ...
    G, x_max, sigma_w, ns, nbins);

%% preliminary checks

delta_expected = [1; 1];
tc.verifyEqual(delta, delta_expected, 'delta');

tc.verifyLessThanOrEqual(sum(pmf(:)), 1, ...
    'Total probability in pmf should not exceed 1.');

%% check pmf against the expected pmf

pmf_Gx_1d = zeros(1, 12);
pmf_Gx_1d(5:7) = 1;
pmf_Gx_1d = pmf_Gx_1d ./ nnz(pmf_Gx_1d);

% odd-length normpmf (of ample length)
bb = -30.5:30.5;
normcdf_at_bb = normcdf(bb, 0, sigma_w(1));
normpmf = normcdf_at_bb(2:end) - normcdf_at_bb(1:end-1);

pmf_y_1d = conv(pmf_Gx_1d, normpmf, 'same');

% figure(1)
% stem(pmf_Gx_1d, 'r');
% hold on
% stem(normpmf, 'g');
% stem(pmf_y_1d, 'b');

pmf_expected = pmf_y_1d' * pmf_y_1d;

tc.verifyEqual(pmf, pmf_expected, 'AbsTol', 1e-12, ...
    'pmf does not match the expected pmf.');


end

