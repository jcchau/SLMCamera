function testMinimumNeighbor2DRandom(tc)
% testMinimumNeighbor2D checks that MIMOCapacity.minimumNeighbor works
% correctly for a random 2D example.  
% Expected result was generated automatically using a different algorithm
% than the method under test. 

nrows = randi(6);
ncols = randi(6);

A = rand(nrows, ncols);

out = MIMOCapacity.minimumNeighbor(A, 2);

expected_out = zeros(size(A));
for row = 2:size(A,1)-1
    for col = 2:size(A,2)-1
        expected_out(row,col) = min(min(A(row-1:row+1, col-1:col+1)));
    end % for col
end % for row

tc.verifyEqual(out, expected_out, ...
    ['The output did not match the independently-generated expected ' ...
    'output.']);

end

