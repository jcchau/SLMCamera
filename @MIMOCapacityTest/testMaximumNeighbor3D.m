function testMaximumNeighbor3D(tc)
% testMaximumNeighbor3D tests MIMOCapacity.maximumNeighbor using a 3D
% matrix input.  
%
% Like in testMaximumNeighbor2DRandom, only the middle of the output is
% thoroughly checked.  The outer edges/faces are checked less rigorously.  

matsize = randi([2, 7], 1, 3);
A = rand(matsize);

out = MIMOCapacity.maximumNeighbor(A);

%% Expected middle of output

expected_out = zeros(matsize);

for ii = 2:matsize(1)-1
    for jj = 2:matsize(2)-1
        for kk = 2:matsize(3)-1
            expected_out(ii,jj,kk) = max(max(max( ...
                A(ii-1:ii+1, jj-1:jj+1, kk-1:kk+1))));
        end % for kk
    end % for jj
end % for ii

tc.verifyEqual(out(2:matsize(1)-1, 2:matsize(2)-1, 2:matsize(3)-1), ...
    expected_out(2:matsize(1)-1, 2:matsize(2)-1, 2:matsize(3)-1), ...
    'The middle of the output did not match the expected output.');

%% additional checks

% that no output is smaller than the input
tc.verifyGreaterThanOrEqual(out, A, ...
    'No output elements should be smaller than the corresponding input.');

% that no output is larger than the maximum input
tc.verifyLessThanOrEqual(out, max(A(:)), ...
    'No output should be larger than the largest input.');

% Check each of the 6 faces:
% That the outer faces of the output are not smaller than the adjacent
% parallel face one layer in.

tc.verifyGreaterThanOrEqual(out(1, :, :), A(2, :, :));
tc.verifyGreaterThanOrEqual(out(end, :, :), A(end-1, :, :));

tc.verifyGreaterThanOrEqual(out(:, 1, :), A(:, 2, :));
tc.verifyGreaterThanOrEqual(out(:, end, :), A(:, end-1, :));

tc.verifyGreaterThanOrEqual(out(:, :, 1), A(:, :, 2));
tc.verifyGreaterThanOrEqual(out(:, :, end), A(:, :, end-1));

end
