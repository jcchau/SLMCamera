function testCalculateMinimumDiffEntropyFromPmfConvergesGaussianOdd(tc)
% testCalculateMinimumDiffEntropyFromPmfConvergesGaussianOdd verifies using
% a Gaussian RV that as the number of bins increase, the difference between
% the true and minimum differential entropy approaches zero.
%
% Here, we only use odd numbers of bins (so that the difference between
% using even and odd number of bins doesn't interfere with our simple/naive
% test of convergence.
%
% Also verifies that the true differential entropy is always at least as
% large as the calculated minimum differential entropy.  

y_min = -5;
y_max = 5;
sigma = 1;

% Differential entropy of a Gaussian-distributed RV with standard deviation
% sigma.  
de_true = log(sigma * sqrt(2*pi*exp(1)));

% 10 is the smallest number of even bins that satisfies the requirement
% that the pdf is concave for the bins around the peak.
% 5 is the smallest number of odd bins that satisfiest this condition.
list_nbins = 10 * (1:10).^2;
% Change to odd bin counts
list_nbins = sort([ 5, 7, list_nbins-1, list_nbins+1 ], 2, 'ascend');

% array to store the calculated minimum differential entropy
de_min = zeros(length(list_nbins), 1);

for ibin = 1:length(list_nbins)
    
    nbins = list_nbins(ibin);
    bin_size = (y_max - y_min) ./ nbins;
    bin_boundaries = linspace(y_min, y_max, nbins+1);
    
    % compute the pmf
    cdf_at_boundaries = normcdf(bin_boundaries, 0, sigma);
    pmf = cdf_at_boundaries(2:end) - cdf_at_boundaries(1:end-1);
    
    % run the method under test
    de_min(ibin) = MIMOCapacity.calculateMinimumDiffEntropyFromPmf(pmf, ...
        bin_size);
    
end % for ibin

%% do checks

tc.verifyLessThanOrEqual(de_min, de_true, ...
    ['The minimum differential entropy should be no larger than the ' ...
    'actual differential entropy.']);

diff_true_min = de_true - de_min;
tc.verifyLessThan(diff_true_min(2:end), diff_true_min(1:end-1), ...
    ['As the number of bins increase, de_min should converge ' ...
    'toward de_true.']);

% figure();
% plot(list_nbins, repmat(de_true, length(list_nbins), 1), 'b-');
% hold on;
% plot(list_nbins, de_min, 'ro-');

end

