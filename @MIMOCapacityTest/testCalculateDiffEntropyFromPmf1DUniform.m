function testCalculateDiffEntropyFromPmf1DUniform(tc)
% testCalculateDiffEntropyFromPmf1DUniform tests method
% calculateDiffEntropyFromPmf using a 1D uniformly-distributed RV.

y_min = -rand() / rand();
y_max = rand() / rand();

unif_min = rand() * (y_max - y_min) + y_min;
unif_max = rand() * (y_max - unif_min) + unif_min;

unif_pdf = 1 / (unif_max - unif_min);

nbins = randi([1e3, 1e4]);

bin_size = (y_max - y_min) / nbins;

bin_boundaries = linspace(y_min, y_max, nbins+1);

pmf = zeros(nbins, 1);
for ii = 1:nbins
    if(bin_boundaries(ii) >= unif_min && bin_boundaries(ii+1) <= unif_max)
        pmf(ii) = unif_pdf * bin_size;
    elseif(bin_boundaries(ii) < unif_min && ...
            bin_boundaries(ii+1) > unif_min)
        pmf(ii) = (bin_boundaries(ii+1) - unif_min) / ...
            (unif_max - unif_min);
    elseif(bin_boundaries(ii) < unif_max && ...
            bin_boundaries(ii+1) > unif_max)
        pmf(ii) = (unif_max - bin_boundaries(ii)) / (unif_max- unif_min);
    end
end

expected_de = log(unif_max - unif_min);

%% method under test

de = MIMOCapacity.calculateDiffEntropyFromPmf(pmf, bin_size);

%% check

tc.verifyGreaterThanOrEqual(de, expected_de, ...
    [ 'The differential entropy calculated from the PMF may be greater' ...
    ' than the actual differential entropy, but should not be less.']);

maximum_de = log(unif_max - unif_min + 2*bin_size);
tc.verifyLessThan(de, maximum_de, ...
    ['The differential entropy calculated from the PMF, ' ...
    'even considering the additional entropy from the bins that ' ...
    'contain the extreme values of the uniformly distributed random ' ...
    'variable, should be less than if the uniform distribution ' ...
    'spanned two extra bins.']);

end