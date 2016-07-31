function testComputeQTTransformDependentColMoreTx(tc)
% Tests method MIMOCapacity.computeQTTransform using a random channel
% matrix with more columns than rows (guaranteeing that some columns are
% linearly dependent).

numcols = randi(20);
numrows = randi(numcols);
G = rand(numrows, numcols);

tc.assumeEqual(rank(G), numrows, ...
    'Such a random channel matrix should have full row rank.');

% Method under test
Q = MIMOCapacity.computeQTTransform(G);

% Check the size
% Number of rows and number of columns should equal the rank of G.  
tc.verifyEqual(size(Q), [numrows, numrows]);

% Check that each column is orthogonal
tc.verifyEqual(rank(Q), numrows);

% Check that each column is normalized
for ii = 1:size(Q,2)
    tc.verifyEqual(norm(Q(:,ii)), 1, 'AbsTol', 1e-15, ...
        sprintf('Column %d is not normalized.', ii));
end

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

