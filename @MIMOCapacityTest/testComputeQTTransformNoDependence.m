function testComputeQTTransformNoDependence(tc)
% Tests method MIMOCapacity.computeQTTransform using a random channel
% matrix with linearly-independent columns (i.e., full column rank).  
%
% This test is modeled after MIMOCapacityLBUnifXTest.gramSchmidtTest.

% Select channel matrix dimensions to make full column rank likely.  
numrows = randi(20);
numcols = randi(numrows);
G = rand(numrows, numcols);

tc.assumeEqual(rank(G), numcols, ...
    ['Test assumes that channel matrix G has no linear dependence ' ...
    'between columns.']);

% Method under test
Q = MIMOCapacity.computeQTTransform(G);

% Check size
tc.verifyEqual(size(Q), [numrows, numcols], ...
    'Q is the wrong size.');

% Check that each column is orthogonal.
for ii = 1:numcols
    for jj = 1:ii-1
        tc.verifyEqual(dot(Q(:,ii), Q(:,jj)), 0, 'AbsTol', 1e-15, ...
            sprintf('Columns %d and %d are not orthogonal.', ii, jj));
    end
end

% Check that each column is normalized.
for ii = 1:numcols
    tc.verifyEqual(norm(Q(:,ii)), 1, 'AbsTol', 1e-15, ...
        sprintf('Column %d is not orthogonal.', ii));
end

% Check that each column of G can be represented as a linear combination of
% the columns of Q.
for ii = 1:numcols
    g = G(:,ii);
    inBasisQ = Q' * g;
    backToOrigBasis = Q * inBasisQ;
    tc.verifyEqual(backToOrigBasis, g, 'AbsTol', 1e-14, ...
        sprintf( ...
        'Column %d of G cannot be represented as the columns of Q.', ii));
end

end

