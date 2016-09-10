function testApproximateUnifXLBForNegligibleNoiseAndFullColRankGUnitaryT(tc)
% Verifies that the result of
% approximateUnifXLBForNegligibleNoiseAndFullColRankG does not change when
% a unitary transform is applied to the channel matrix G.
%
% This is more of a sanity check.  It should not fail unless
% approximateUnifXLBForNegligibleNoiseAndFullColRankG is really wrong.  

n_r = randi(20);
n_t = n_r;
G = rand(n_r, n_t);

sigma_w = rand();
% Not necessarily high SNR
x_max = sigma_w / min(G(:)) / rand();

% Random unitary transform matrix
U = orth(rand(n_r)-0.5);

%% On the original G
[lb_A, ~, h_x_A, h_QTransposeGx_A, h_QTransposew_A] = ...
    MIMOCapacityLBUnifX.approximateUnifXLBForNegligibleNoiseAndFullColRankG( ...
    G, x_max, sigma_w);

%% On the U'*G
[lb_B, ~, h_x_B, h_QTransposeGx_B, h_QTransposew_B] = ...
    MIMOCapacityLBUnifX.approximateUnifXLBForNegligibleNoiseAndFullColRankG( ...
    U'*G, x_max, sigma_w);

%% Check that the results are equal

% The unitary transform on the channel matrix should not change mutual
% information.
tc.verifyEqual(lb_B, lb_A, 'AbsTol', 1e-15);

% The unitary transform on the channel has no effect on h(x).
tc.verifyEqual(h_x_B, h_x_A, 'AbsTol', 1e-15);

% h(Q'*G*x) == h(Q'*U'*G*x) should both equal h(G*x).
tc.verifyEqual(h_QTransposeGx_B, h_QTransposeGx_A, 'AbsTol', 1e-15);

% Again, the unitary transform on G has no effect on h(Q'*w).
tc.verifyEqual(h_QTransposew_B, h_QTransposew_A, 'AbsTol', 1e-15);


end

