function testCalculateMinimumDiffEntropyFromPmf2D(tc)
% testCalculateMinimumDiffEntropyFromPmf2D tests
% calculateMinimumDiffEntropyFromPmf2D using a 2D Gaussian RV (2
% independent dimensions each with random standard deviation).  

sigma = rand(1,2) ./ rand(1,2);
y_max = 6 .* sigma;
y_min = -y_max;

% From p.100 of lab book #3.
de_true = 0.5 * sum(log((2*pi*exp(1) .* sigma.^2)));

% A progression of nbins that should steadily decrease the difference
% between de_true and de_min.
list_nbins = [ ...
    5, 10;
    7, 10;
    9, 10;
    25, 10;
    25, 16;
    25, 36;
    49, 36;
    49, 64;
    81, 64;
    81, 100;
    121, 100; 
    121, 144;
    169, 144;
    169, 196;
    225, 196;
    225, 256;
    289, 256;
    289, 324;
    361, 324;
    361, 400];
list_nbins = [ list_nbins;
    ((21:31).^2)', ((22:32).^2)';
    1447, 1448;
    2047, 2048;
    2895, 2896;
    4095, 4096;
    5793, 5792];

length_list_nbins = size(list_nbins,1);

de_min = zeros(length_list_nbins, 1);

for ibin = 1:length_list_nbins
    
    nbins = list_nbins(ibin, :);
    bin_size = (y_max - y_min) ./ nbins;
    
    % Bin boundaries for each of the two dimensions
    bin_boundaries_1 = linspace(y_min(1), y_max(1), nbins(1)+1);
    bin_boundaries_2 = linspace(y_min(2), y_max(2), nbins(2)+1);
    
    % compute the pmf along each dimension
    cdf_at_bounds_1 = normcdf(bin_boundaries_1, 0, sigma(1));
    pmf_1 = cdf_at_bounds_1(2:end) - cdf_at_bounds_1(1:end-1);
    cdf_at_bounds_2 = normcdf(bin_boundaries_2, 0, sigma(2));
    pmf_2 = cdf_at_bounds_2(2:end) - cdf_at_bounds_2(1:end-1);
    
    % combine into single pmf with both dimensions independent
    pmf = pmf_1(:) * pmf_2(:)';
    
    % run the method under test
    de_min(ibin) = MIMOCapacity.calculateMinimumDiffEntropyFromPmf(pmf, ...
        bin_size');
    
end % for ibin

%% do checks

tc.verifyLessThanOrEqual(de_min, de_true, ...
    ['The minimum differential entropy should be no larger than the ' ...
    'actual differential entropy.']);

diff_true_min = de_true - de_min;
tc.verifyLessThan(diff_true_min(2:end), diff_true_min(1:end-1), ...
    ['As the number of bins increase, de_min should converge ' ...
    'toward de_true.']);

% From the plot on p.19 in lab book #4, I expect the diff_true_min to be
% well below 1e-4 for 4095x4096 bins.
tc.verifyLessThan(diff_true_min(end), 1e-4, ...
    sprintf('The tolerance should be < 1e-4 for %d bins.', ...
    prod(list_nbins(end,:))));

loglog(prod(list_nbins, 2), diff_true_min, '+-');
xlabel('Number of bins')
ylabel('de_{true} - de_{min}')
title('Difference between the true differential entropy and the minimum.')

end

