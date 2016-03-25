function testCalculateDiffEntropyFromPmf2DGaussian(tc)
% testCalculateDiffEntropyFromPmf2DGaussian checks that the result of
% calculateDiffEntropyFromPmf converges towards the true result for a 2D
% Gaussian RV.  

sigma = rand(1,2) ./ rand(1,2);
y_max = 6 .* sigma;
y_min = -y_max;

% From p.100 of lab book #3.
% de_true = 0.5 * sum(log((2*pi*exp(1) .* sigma.^2)));
de_true = 0.5 * sum(log((2*pi .* sigma.^2)) + 1);

% Exponentially growing number of bins (starting from 2).
% Randomly subtracts 1 to test both even and odd numbers of bins
length_list_nbins = 12;
list_nbins = repmat((2.^(1:length_list_nbins))', 1, 2) + ...
    randi([-1, 0], length_list_nbins, 2);

% calculateDiffEntropyFromPmf computes the maximum entropy given the
% average PDF for each bin (provided in the PMF).
de_max = zeros(length_list_nbins, 1);

for ibin = 1:length_list_nbins
    
    nbins = list_nbins(ibin, :);
    bin_size = (y_max - y_min) ./ nbins; % a row vector
    
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
    
    % randomly make bin_size either the area of each bin or a column vector
    % to test both cases in calculateDiffEntropyFromPmf
    if(rand()>0.5)
        bin_size = prod(bin_size);
    else
        bin_size = bin_size(:);
    end
    
    % run the method under test
    de_max(ibin) = MIMOCapacity.calculateDiffEntropyFromPmf(pmf, bin_size);
        
end % for ibin

%% do checks

tc.verifyGreaterThanOrEqual(de_max, de_true, ...
    ['The maximum differential entropy should be no smaller than the ' ...
    'actual differential entropy.']);

diff_true_max = de_max - de_true;
tc.verifyLessThan(diff_true_max(2:end), diff_true_max(1:end-1), ...
    ['As the number of bins increase, de_max should converge ' ...
    'toward de_true.']);

% Expect that for a large number of bins, that the error is small
tc.verifyLessThan(diff_true_max(end), 1e-5, ...
    sprintf('The tolerance should be < 1e-5 for %d bins.', ...
    prod(list_nbins(end,:))));

loglog(prod(list_nbins, 2), diff_true_max, 'b+');
hold on
loglog(prod(list_nbins, 2), abs(diff_true_max), 'rx-');

end

