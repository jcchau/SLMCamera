function testMaximumNeighbor2DRandom(tc)
% testMaximumNeighbor2DRandom verifies MIMOCapacity.maximumNeighbor using a
% randomly generated 2D matrix.
%
% Note that this test method does not thoroughly check the edges or corners
% of the output since they're somewhat complicated to check (and it didn't
% seem worth the effort to reimplement method maximumNeighbor, but
% completely differently to verify the edges).

nrows = randi([2,7]);
ncols = randi([2,7]);

A = rand(nrows, ncols);

out = MIMOCapacity.maximumNeighbor(A);

%% expected output

% for the middle
expected_out = zeros(nrows, ncols);
for row = 2:size(A,1)-1
    for col = 2:size(A,2)-1
        expected_out(row,col) = max(max(A(row-1:row+1, col-1:col+1)));
    end % for col
end % for row

%% check

% the center elements
tc.verifyEqual(out(2:nrows-1, 2:ncols-1), ...
    expected_out(2:nrows-1, 2:ncols-1), ...
    ['The middle of the output did not match the middle of the ' ...
    'expected output.']);

% that no output is smaller than the input
tc.verifyGreaterThanOrEqual(out, A, ...
    'The output should not be smaller than the input.');

% that no output is larger than the maximum input
tc.verifyLessThanOrEqual(out, max(A(:)), ...
    'The output should not be larger than the largest input.');

% and check that the edges are at least as large as the adjacent parallel
% edges of the input.

tc.verifyGreaterThanOrEqual(out(1,:), A(2,:));
tc.verifyGreaterThanOrEqual(out(end,:), A(end-1,:));
tc.verifyGreaterThanOrEqual(out(:,1), A(:,2));
tc.verifyGreaterThanOrEqual(out(:,end), A(:,end-1));

end

