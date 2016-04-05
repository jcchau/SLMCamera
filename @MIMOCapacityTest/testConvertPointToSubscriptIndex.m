function testConvertPointToSubscriptIndex(tc)
% testConvertPointToSubscriptIndex

npoints = randi(10);
num_dims = randi(5);

y_min = rand(1, num_dims);

nbins = randi(10, 1, num_dims);

delta = 1./nbins;

expected_indices = zeros(npoints, num_dims);
y = zeros(npoints, num_dims);
for jj = 1:num_dims
    expected_indices(:,jj) = randi(nbins(jj), npoints, 1);
    y(:,jj) = (expected_indices(:,jj)-1) .* delta(jj) + ...
        repmat(y_min(jj), npoints, 1);
end

% Add some random value to y that shouldn't change the resulting index.
y_rand = bsxfun(@times, rand(npoints, num_dims), delta);
y = y + y_rand;


msi = MIMOCapacity.convertPointToSubscriptIndex(y, y_min, delta, nbins);


tc.verifyEqual(msi, expected_indices);

end

