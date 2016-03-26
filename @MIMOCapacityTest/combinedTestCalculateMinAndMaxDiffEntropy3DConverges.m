function combinedTestCalculateMinAndMaxDiffEntropy3DConverges(tc)
% combinedTestCalculateMinAndMaxDiffEntropy3DConverges verifies that the
% differential entropy calculated by calculateDiffEntropyFromPmf and
% calculateMinimumDiffEntropyFromPmf bounds and converges toward the true
% differential entropy for 3 independent Gaussian random variables (3D).  

sigma = rand(1,3) ./ rand(1,3);
y_max = 6 .* sigma;
y_min = -y_max;

% From p.100 of lab book #3.
de_true = 0.5 * sum(log((2*pi .* sigma.^2)) + 1);

% Allow up to 1.8 GB, each double takes 8 bytes, and allow for up to 5
% matrices of equal size to matrix pmf.  
max_num_bins = 1.8e9 / 8 / 5; % 45e6 (about 356^3)

% Generate a progression of increasing bin counts
list_nbins = [10, 10, 10];
while(prod(list_nbins(end,:)) <= max_num_bins)
    dim_to_increase = randi(3);
    new_nbins = list_nbins(end,:);
    new_nbins(dim_to_increase) = ...
        randi(round([1.2, 1.5] .* new_nbins(dim_to_increase)));
    list_nbins = [ list_nbins; new_nbins ];
end % while(prod(list_nbins(end,:)) <= max_num_bins)
% The last row of list_nbins exceeds max_num_bins
list_nbins = list_nbins(1:end-1, :);

length_list_nbins = size(list_nbins, 1);

% To hold the results
de_min = zeros(length_list_nbins, 1);
de_max = zeros(length_list_nbins, 1);

for ibin = 1:length_list_nbins
    
    nbins = list_nbins(ibin, :);
    bin_size = (y_max - y_min) ./ nbins; % a row vector
    
    % bin boundaries (row vectors)
    bb1 = linspace(y_min(1), y_max(1), nbins(1)+1);
    bb2 = linspace(y_min(2), y_max(2), nbins(2)+1);
    bb3 = linspace(y_min(3), y_max(3), nbins(3)+1);
    
    % compute the pmf along each dimension (also row vectors)
    cdf1 = normcdf(bb1, 0, sigma(1));
    cdf2 = normcdf(bb2, 0, sigma(2));
    cdf3 = normcdf(bb3, 0, sigma(3));
    pmf1 = cdf1(2:end) - cdf1(1:end-1);
    pmf2 = cdf2(2:end) - cdf2(1:end-1);
    pmf3 = cdf3(2:end) - cdf3(1:end-1);
    
    % reshape pmf3 to be along the 3rd dimension
    pmf3 = reshape(pmf3, 1, 1, nbins(3));
    
    % compute the 3D PMF
    pmf = bsxfun(@times, pmf1' * pmf2, pmf3);
    
    % convert bin_size into a column vector for the methods under test
    bin_size = bin_size(:);
    
    % run both methods under test
    de_min(ibin) = MIMOCapacity.calculateMinimumDiffEntropyFromPmf(pmf, ...
        bin_size);
    de_max(ibin) = MIMOCapacity.calculateDiffEntropyFromPmf(pmf, bin_size);
    
end % for ibin

%% do checks

% That they're on the correct side of the true differential entropy
tc.verifyGreaterThanOrEqual(de_max, de_true, ...
    ['The maximum differential entropy should be no smaller than the ' ...
    'actual differential entropy.']);
tc.verifyLessThanOrEqual(de_min, de_true, ...
    ['The minimum differential entropy should be no larger than the  ' ...
    'actual differential entropy.']);

% That the difference between the minimum and maximum bound shrinks as the
% number of bins increase
diff_max_min = de_max - de_min;
tc.verifyLessThan(diff_max_min(2:end), diff_max_min(1:end-1), ...
    ['As the number of bins increase, the difference between de_max ' ...
    'and de_min should tend towards zero.']);

% Also check that de_max and de_min both individually converges toward the
% true differential entropy (instead of just one of them converging).
diff_max_true = de_max - de_true;
tc.verifyLessThan(diff_max_true(2:end), diff_max_true(1:end-1), ...
    ['As the number of bins increase, the difference between de_max ' ...
    'and de_true should tend towards zero.']);

diff_true_min = de_true - de_min;
% diff_true_min may rise locally (due to changes between the algorithm for
% even number of bins and odd number of bins), but in the long term,
% diff_true_min should tend toward zero as the number of bins increase.
tc.verifyLessThan(diff_true_min(12:end), diff_true_min(1:end-11), ...
    ['As the number of bins increase, the difference between de_true ' ...
    'and de_min should tend towards zero.']);

% Plot the result for visual check.
% figure();
% total_nbins = prod(list_nbins, 2);
% loglog(total_nbins, diff_max_min, 'b+-');
% hold on
% loglog(total_nbins, diff_max_true, 'g+-');
% loglog(total_nbins, diff_true_min, 'r+-');
% legend('max-min', 'max-true', 'true-min');

end

