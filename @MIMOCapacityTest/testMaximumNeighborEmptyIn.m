function testMaximumNeighborEmptyIn(tc)
% testMaximumNeighborEmptyIn verifies that MIMOCapacity.maximumNeighbor
% returns an empty matrix if provided with an empty input matrix.  

nd = randi([2,5]);

% create an empty nd-dimension matrix
in = zeros(zeros(1, nd));

out = MIMOCapacity.maximumNeighbor(in);

tc.verifyTrue(isempty(out), ...
    'minimumNeighbor did not return an empty matrix.');

tc.verifyEqual(ndims(out), nd, ...
    ['The number of dimensions of the output differs from the number ' ...
    'of dimensions of the input.']);

end

