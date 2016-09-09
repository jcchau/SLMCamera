function testComputeCapacityLBUnifGxAgainstUnifX(tc)
% testComputeCapacityLBUnifGxAgainstUnifX tests
% MIMOCapacity.computeCapacityLBUnifGx against
% MIMOCapacity.computeCapacityLBUnifX.  When the channel matrix G has full
% column rank, the two should approximately compute the same lower bound on
% capacity.  When the channel matrix does NOT have full column rank, the
% UnifGx lower bound should be larger than the UnifX lower bound.  
%
% Note however that when the SNR is low, the UnifX lower bound may compute
% a lower lower bound.  So for this test, we use high SNR.  

n_r = randi(6);
n_t = randi(6);

G = rand(n_r, n_t);

% Determine if the channel matrix has full column rank.  
% (Since G is randomly generated)
is_full_col_rank = n_r >= n_t;

% Pick a high SNR (using the same approach as in
% testComputeCapacityLBUnifXFullRankHighSNR.  
sigma_w = rand();
weakest_gain = min(G(:));
x_max = 100 * sigma_w / weakest_gain / rand();

%% Run the method under test

lb_nats = MIMOCapacity.computeCapacityLBUnifGx(G, x_max, sigma_w);

%% And for comparison

[~, lb_nats_expected] = ...
    MIMOCapacity.computeCapacityLBUnifX(G, x_max, sigma_w);

%% Check the result

if(is_full_col_rank)
    % Expect lb_nats to be approximately the same as lb_nats_expected
    % (which would be identical to lb_nats_conservative_expected).
    tc.verifyEqual(lb_nats, lb_nats_expected, 'RelTol', 0.03, ...
        ['Expected the UnifGx and then UnifX lower bounds to be ' ...
        'approximately the same.']);
else
    % Expect lb_nats to be larger than lb_nats_expected.
    tc.verifyGreaterThan(lb_nats, lb_nats_expected, ...
        ['Expected the UnifGx LB to be larger than the UnifX lower '...
        'bound since G does not have full column rank.']);
end % if(is_full_col_rank)

end

