function testSmithCapacity_p_Y(tc)
% testSmithCapacity_p_Y tests SmithCapacity.p_Y against a slower/simpler
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
y = 2.*(rand(numsamples, 1)-0.5) .* A;

%% run method under test

pdf_test = SmithCapacity.p_Y(y, poi, voi);

%% compute expected result

pdf_expected = zeros(size(y));
for ii = 1:length(poi)
    pdf_expected = pdf_expected + voi(ii)*normpdf(y, poi(ii), 1);
end % for ii

%% verify result

tc.verifyEqual(pdf_test, pdf_expected, ...
    'test results do not match expected results.');

end

