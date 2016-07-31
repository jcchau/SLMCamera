function testSimplifyChannelMatrixStepB(tc)
% Tests Step B in MIMOCapacity.simplifyChannelMatrix to ensure that columns
% of the channel matrix that are multiples of each other are properly
% combined.  

numrows = randi(20);
numcols = randi(20);

G_B_expected = rand(numrows, numcols);

%% Split some columns to generate columns that are multiples of each other

num_cols_to_split = randi(numcols);

G = zeros(numrows, numcols+num_cols_to_split);
G(:,1:numcols) = G_B_expected;

for ii = 1:num_cols_to_split
    col_to_split = randi(numcols);
    fraction = rand();
    
    G(:,numcols+ii) = fraction .* G(:,col_to_split);
    G(:,col_to_split) = (1-fraction) .* G(:,col_to_split);
end % for ii

% Step B of MIMOCapacity.simplifyChannelMatrix should recombine these split
% columns to get back to G_B_expected.  

%% Method under test
[Q, G_B] = MIMOCapacity.simplifyChannelMatrix(G);

%% Check G_B

tc.verifyEqual(G_B, G_B_expected, 'AbsTol', 1e-15, ...
    'The calculated G_B does not match the expected G_B.');

%% Quick checks on Q

rank_G_B_expected = rank(G_B_expected);
tc.verifyEqual(size(Q,2), rank_G_B_expected, ...
    ['The Q^T Transform should yield only as many receivers as the ' ...
    'rank of G_B.']);

% Check that the calculated Q' Transform is invertible.
tc.verifyEqual(Q*(Q'*G_B), G_B, 'AbsTol', 1e-15, ...
    'Multiplying by Q should undo the Q^T transform.');

end

