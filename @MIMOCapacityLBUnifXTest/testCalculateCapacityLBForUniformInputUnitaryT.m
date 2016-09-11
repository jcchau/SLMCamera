function testCalculateCapacityLBForUniformInputUnitaryT(tc)
% Verifies that the result of calculateCapacityLBForUniformInput does not
% change when an unitary transform is applied to G (when x_max is the same
% for all transmitters).

n_r = randi(20);
n_t = n_r;
G = rand(n_r, n_t);

sigma_w = rand();
% Not necessarily high SNR
x_max = sigma_w / min(G(:)) / rand();

% Random unitary transform matrix
U = orth(rand(n_r)-0.5);

variance_noise_out = repmat(sigma_w^2, n_r, 1);
max_nbins = 1e9/8/5;
min_trials = 1e9;

%% On the original G

[nats_A, min_nats_A, ~, ~, ~, h_y_A, ~] = ...
    MIMOCapacityLBUnifX.calculateCapacityLBForUniformInput( ...
    G, x_max, variance_noise_out, max_nbins, min_trials);

%% On the U'*G

[nats_B, min_nats_B, ~, ~, ~, h_y_B, ~] = ...
    MIMOCapacityLBUnifX.calculateCapacityLBForUniformInput( ...
    U'*G, x_max, variance_noise_out, max_nbins, min_trials);

%% Compare the results

tc.verifyEqual(nats_B, nats_A, 'RelTol', 0.01);

% Kinda expecting more variance in min_nats, but we'll see.
% min_nats might be zero, so use 0.02*nats_A to emulate 'RelTol'.
tc.verifyEqual(min_nats_B, min_nats_A, 'AbsTol', 0.02*nats_A);

tc.verifyEqual(h_y_B, h_y_A, 'RelTol', 0.01);

end

