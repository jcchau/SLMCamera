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
G_B = MIMOCapacity.simplifyChannelMatrix(G);

%% Check G_B

tc.verifyEqual(G_B, G_B_expected, 'AbsTol', 1e-15, ...
    'The calculated G_B does not match the expected G_B.');

end

