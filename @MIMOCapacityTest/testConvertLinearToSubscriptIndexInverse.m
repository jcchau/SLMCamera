function testConvertLinearToSubscriptIndexInverse(tc)
% testConvertLinearToSubscriptIndexInverse checks that
% MIMOCapacity.convertLinearToSubscriptIndex is the inverse of
% convertToLinearIndex.

% number of dimensions)
nd = randi(5);

% size of matrix to be indexed
matsize = randi(7, 1, nd);

% pick an element
ms_orig = arrayfun(@(x) randi(x), matsize);

weights = MIMOCapacity.convertToLinearIndexWeights(matsize);
li = MIMOCapacity.convertToLinearIndex(weights, ms_orig);
ms_test = MIMOCapacity.convertLinearToSubscriptIndex(weights, li);

tc.verifyEqual(ms_test, ms_orig, ...
    'ms_test should be the same subscript index as ms_orig.');

end

