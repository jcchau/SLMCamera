function gramSchmidtTest(tc)
% gramSchmidtTest checks method gramSchmidt using a random V

numrows = randi(20);
numcols = randi(numrows);
V = rand(numrows, numcols);

tc.assumeEqual(rank(V), numcols, ...
    'Method gramSchmidt requires full column rank.');

W = MIMOCapacityLBUnifX.gramSchmidt(V);

% Check size.
tc.verifyEqual(size(W), [numrows, numcols], ...
    'W is the wrong size.');

% Check that each column is orthogonal.
for ii = 1:numcols
    for jj = 1:ii-1
        tc.verifyEqual(dot(W(:,ii), W(:,jj)), 0, 'AbsTol', 1e-15, ...
            sprintf('Columns %d and %d are not orthogonal.', ii, jj));
    end
end

% Check that each column is normalized.
for ii = 1:numcols
    tc.verifyEqual(norm(W(:,ii)), 1, 'AbsTol', 1e-15, ...
        sprintf('Column %d is not orthogonal.', ii));
end

% Check that each column of V can be represented as a linear combination of
% the columns of W.
for ii = 1:numcols
    v = V(:,ii);
    inBasisW = W' * v;
    backToOrigBasis = W * inBasisW;
    tc.verifyEqual(backToOrigBasis, v, 'AbsTol', 1e-14, ...
        sprintf( ...
        'Column %d of V cannot be represented as the columns of W.', ii));
end

end

