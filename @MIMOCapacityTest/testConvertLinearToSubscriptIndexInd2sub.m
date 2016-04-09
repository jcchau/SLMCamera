function testConvertLinearToSubscriptIndexInd2sub(tc)
% testConvertLinearToSubscriptIndexInd2sub compares the result of
% convertLinearToSubscriptIndex against that of ind2sub.

% number of dimensions)
nd = randi(5);

% size of matrix to be indexed
matsize = randi(7, 1, nd);

% pick an element
li = randi(prod(matsize));

weights = MIMOCapacity.convertToLinearIndexWeights(matsize);
ms_test = MIMOCapacity.convertLinearToSubscriptIndex(weights, li);

[ms{1:nd}] = ind2sub(matsize, li);

ms_expected = cell2mat(ms);

tc.verifyEqual(ms_test, ms_expected);

end

