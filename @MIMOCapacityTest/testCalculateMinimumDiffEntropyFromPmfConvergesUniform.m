function testCalculateMinimumDiffEntropyFromPmfConvergesUniform(tc)
% testCalculateMinimumDiffEntropyFromPmfConvergesUniform verifies that as
% the number of bins increase, the difference between the true differential
% entropy and the minimum differential entropy converges to (or approaches)
% zero.
%
% Also verifies that that the calculated minimum differential entropy is
% less than (or equal to) the true differential entropy.
%
% Uses a uniform random variable for this test.

y_min = -1;
y_max = 1;

% Need to center the PMF so that the "peak" of the PDF is in the center of
% the PMF.  
unif_max = rand();
unif_min = -unif_max;

probability_density = 1/(unif_max - unif_min);

% true differential entropy of a uniform random variable in nats
de_true = log(unif_max - unif_min);

% List of nbins to consider
list_nbins = 2.^(0:10);

de_min = zeros(length(list_nbins), 1);

% Indicates whether the conditions for the method-under-test is satisfied.
mut_condition_satisfied = false(length(list_nbins), 1);

for ibin = 1:length(list_nbins)
    
    nbins = list_nbins(ibin);
    bin_size = (y_max - y_min) / nbins;
    bin_boundaries = linspace(y_min, y_max, nbins+1);
    
    pmf = zeros(nbins, 1);
    for ii = 1:nbins
        if(bin_boundaries(ii) >= unif_min && ...
                bin_boundaries(ii+1) <= unif_max)
            % the bin is entirely within the range of values of the uniform
            % random variable.
            pmf(ii) = probability_density * bin_size;
            mut_condition_satisfied(ibin) = true;
        elseif(bin_boundaries(ii) <= unif_min && ...
                bin_boundaries(ii+1) >= unif_max)
            % the bin contains the whole range of possible values of the
            % uniform random variable.
            pmf(ii) = 1;
        elseif(bin_boundaries(ii) < unif_min && ...
                bin_boundaries(ii+1) > unif_min)
            % bin includes left end of the pdf
            pmf(ii) = probability_density * ...
                (bin_boundaries(ii+1) - unif_min);
        elseif(bin_boundaries(ii) < unif_max && ...
                bin_boundaries(ii+1) > unif_max)
            % bin includes right end of the pdf
            pmf(ii) = probability_density * ...
                (unif_max - bin_boundaries(ii));
        end
    end % for ii = 1:nbins
    
    % run the method under test: calculateMinimumDiffEntropyFromPmf
    de_min(ibin) = MIMOCapacity.calculateMinimumDiffEntropyFromPmf(pmf, ...
        bin_size);
    
    % check that the error is never more than if we chopped a bin off of
    % each end of the pdf.
    de_chopped = log(max(unif_max - unif_min - 2*bin_size, 0));
    tc.verifyGreaterThan(de_min(ibin), de_chopped, ...
        ['The minimum entropy calculated should not be lower than if ' ...
        'the uniform PDF was truncated by two bins.']);
    
end % for ibin

%% checks
% However, calculateMinimumDiffEntropy relies on the assumption that the
% PDF is concave in the bins around the peak.  So only run this test when
% at least one bin is completely within the range of possible values for
% the uniform random variable. 

if(any(mut_condition_satisfied))
    
    % Verify that de_min is always less than or equal de_true.  
    tc.verifyLessThanOrEqual(de_min(mut_condition_satisfied), ...
        de_true + 1e-12, ...
        ['The minimum bound on the differential entropy should always ' ...
        'be less than or equal the true differential entropy.']);

    % Check that the error is decreasing.
    de_error = abs(de_true - de_min(mut_condition_satisfied));
    tc.verifyLessThanOrEqual(de_error(2:end), ...
        de_error(1:end-1) + 1e-12, ...
        ['The error should decrease as the number of bins increase ' ...
        '(until de_min equals de_true).']);

end % if(any(mut_condition_satisfied))

end

