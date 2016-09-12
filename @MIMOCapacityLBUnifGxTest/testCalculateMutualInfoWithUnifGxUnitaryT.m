function testCalculateMutualInfoWithUnifGxUnitaryT(tc)
% Tests calculateMutualInfoWithUnifGx to ensure that a unitary transform on
% the channel matrix does not affect the computed capacity lower bound when
% the noise is isotropic (identical across all of the receivers).

% Here we ensure that G has full *row* rank so we don't have to apply the
% Q'-transform.  Otherwise, we would need to apply the Q'-transform to
% ensure that the received constellation is not "flat" in any direction
% (which would mean that the PMF bin sizes would necessarily be too large
% along that dimension).
n_t = randi(6);
n_r = randi(n_t);
G = rand(n_r, n_t);

sigma_w = rand();
% Not necessarily high SNR
x_max = sigma_w / mean(G(:)) / rand();

% Random unitary transform matrix
U = orth(rand(n_r)-0.5);

%% Prepare parameters

n_r = size(G, 1);

sigma_w = repmat(sigma_w, n_r, 1);

max_nbins_linprog = 6e4; % A tradeoff between speed and precision: 6e4
max_nbins_mem = 1e9/8/5;

%% For the original G

% TODO: move this max_nbins_linprog feature in to MIMOCapacityLBUnifGx.
[uminA, umaxA] = MIMOCapacity.computeUExtremes(G, x_max);
usizeA = umaxA-uminA;

expected_nbins_linprogA = ceil(prod( ...
    max_nbins_linprog^(1/n_r) .* (12*sigma_w+usizeA) ./ usizeA ));

max_nbinsA = min(max_nbins_mem, expected_nbins_linprogA);

lb_natsA = MIMOCapacityLBUnifGx.calculateMutualInfoWithUnifGx( ...
    G, x_max, sigma_w, max_nbinsA);

%% For U'*G

[uminB, umaxB] = MIMOCapacity.computeUExtremes(U'*G, x_max);
usizeB = umaxB-uminB;

expected_nbins_linprogB = ceil(prod( ...
    max_nbins_linprog^(1/n_r) .* (12*sigma_w+usizeB) ./ usizeB ));

max_nbinsB = min(max_nbins_mem, expected_nbins_linprogB);

lb_natsB = MIMOCapacityLBUnifGx.calculateMutualInfoWithUnifGx( ...
    U'*G, x_max, sigma_w, max_nbinsB);

%% Compare

tc.verifyEqual(lb_natsB, lb_natsA, 'RelTol', 0.01);

end

