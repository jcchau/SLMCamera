function testComputeQTTransformDependentRx(tc)
% Tests method MIMOCapacity.computeQTTransform to ensure that any
% linearly-dependent rows of the channel matrix are properly eliminated by
% the Q' transform.  

numrows = randi(20);
numcols = randi(20);

num_indep_rows = randi(numrows);

%% Generate matrix G with size [numrows, numcols] and rank rankG.
G = zeros(numrows, numcols);
G(1:num_indep_rows, :) = rand(num_indep_rows, numcols);

% Fill in remaining rows (num_indep_rows+1:end) of G
num_dependent_rows = numrows - num_indep_rows;
A = 2 * (rand(num_dependent_rows, numrows) - 0.5);
G(num_indep_rows+1:end, :) = A * G;

% Check expected rank
rankG_expected = min(num_indep_rows, numcols);
tc.assumeEqual(rank(G), rankG_expected, ...
    sprintf('Expected G to have rank %d.', rankG_expected));

%% Method under test
Q = MIMOCapacity.computeQTTransform(G);

%% Check results

% Check the size
tc.verifyEqual(size(Q), [numrows, rankG_expected], ...
    'Q is the wrong size.');

% Verify each column of Q is orthogonal
tc.verifyEqual(rank(Q), size(Q,2), ...
    'Columns of Q are not orthogonal (linearly independent).');

% Check that each column is normalized
for ii = 1:size(Q,2)
    tc.verifyEqual(norm(Q(:,ii)), 1, 'AbsTol', 1e-15, ...
        sprintf('Column %d is not normalized.', ii));
end % for ii

% Check that each column of G can be represented as a linear combination of
% the columns of Q.
for ii = 1:numcols
    g = G(:,ii);
    inBasisQ = Q' * g;
    backToOrigBasis = Q * inBasisQ;
    tc.verifyEqual(backToOrigBasis, g, 'RelTol', 1e-4, ...
        sprintf(['Column %d of G cannot be represented as a linear ' ...
        'combination of the columns of Q.', ii]));
end

end

