function testMinimumNeighbor3D(tc)
% testMinimumNeighbor3D verifies that MIMOCapacity.minimumNeighbor returns
% the correct result for a randomly-generated 3D matrix input.
% The expected result is computed using a different algorithm than what's
% used in minimumNeighbor.  

matsize = randi(7, 1, 3);
A = rand(matsize);

out = MIMOCapacity.minimumNeighbor(A, 3);

expected_out = zeros(matsize);

for ii = 2:matsize(1)-1
    for jj = 2:matsize(2)-1
        for kk = 2:matsize(3)-1
            expected_out(ii,jj,kk) = min(min(min( ...
                A(ii-1:ii+1, jj-1:jj+1, kk-1:kk+1))));
        end % for kk
    end % for jj
end % for ii

tc.verifyEqual(out, expected_out, ...
    ['The output did not match the independently-generated expected ' ...
    'result.']);

end

