function testMinimumNeighborSingletonDimension(tc)
% testMinimumNeighborSingletonDimension verifies that minimumNeighbor
% returns all zeros when there is a singleton dimension (since every cell
% is adjacent to the outer boundary.  

num_dimensions = randi(5);

singleton_dim = randi(num_dimensions);

matsize = randi(30, [1, num_dimensions]);
matsize(singleton_dim) = 1;

in = rand(matsize);

out = MIMOCapacity.minimumNeighbor(in, num_dimensions);

tc.verifyEqual(size(out), size(in), ...
    'Input and output size should match.');

tc.verifyEqual(out, zeros(matsize), ...
    'Every element of out should equal zero.');

end

