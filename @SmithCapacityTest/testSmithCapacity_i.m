function testSmithCapacity_i(tc)
% testSmithCapacity_i tests SmithCapacity.i against a slower/simpler
% implementation.  

n = randi([2, 10]);
A = rand() * 2*n;

poi = rand(n, 1);
% poi is in ascending order from -A to A. 
poi = sort(2 .* (poi - 0.5) .* A);

voi = rand(n, 1);
% Normalize probabilites to sum to 1.
voi = voi ./ sum(voi);

numsamples = randi(10);
x = 2.*(rand(numsamples, 1)-0.5) .* A;

%% run method under test

r_test = SmithCapacity.i(x, poi, voi);

%% compute expected result

expected = zeros(size(x));
for ii = 1:numel(x)
    % marginal information density
    expected(ii) = integral( ...
        @(y) SmithCapacity.aloga_over_b( ...
        SmithCapacity.p_N(y-x(ii)), ...
        SmithCapacity.p_Y(y, poi, voi)), ...
        -Inf, Inf);
end % for ii

%% verify result

tc.verifyEqual(r_test, expected, 'AbsTol', 1e-6, ...
    'test results do not match expected results.');

end

