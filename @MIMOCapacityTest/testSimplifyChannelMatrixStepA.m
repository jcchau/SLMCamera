function testSimplifyChannelMatrixStepA(tc)
% Basic check that MIMOCapacity.simplifyChannelMatrix removes rows and
% columns of zeros.
%
% Not designed to be a thorough test of MIMOCapacity.removeZeroRowsAndCols.
% More intendened to catch if MIMOCapacity.simplifyChannelMatrix forgets to
% call removeZeroRowsAndCols.  

numrows = randi(20);
numcols = randi(20);

G = rand(numrows, numcols);

num_zero_rows = 0;
for ii = 1:numrows
    if(rand() < 2/numrows)
        G(ii,:) = 0;
        num_zero_rows = num_zero_rows + 1;
    end
end

num_zero_cols = 0;
for ii = 1:numcols
    if(rand() < 2/numcols)
        G(:,ii) = 0;
        num_zero_cols = num_zero_cols + 1;
    end
end

%% Method under test
G_B = MIMOCapacity.simplifyChannelMatrix(G);

%% Check that removeZeroRowsAndCols was run

G_A_expected = MIMOCapacity.removeZeroRowsAndCols(G);

if(numrows-num_zero_rows > 1)
    % If we have at least 2 non-zero rows, then for the random G, we can
    % assume that G_B equals G_A because the columns would likely not be
    % scaled copies of each other.  
    % Otherwise, with just one non-zero row, the columns are certainly
    % scaled copies of each other.  
    
    if(numcols > num_zero_cols)
        % For simplicity, ignore the special case where we have no non-zero
        % values in the channel matrix G. 
        
        tc.verifyEqual(size(G_B), ...
            [numrows-num_zero_rows, numcols-num_zero_cols], ...
            'G_B is not the expected size.');
    end

    % Check against G_A_expected.
    tc.verifyEqual(G_B, G_A_expected, ...
        ['G_B should equal G_A_expected since we have at least 2 ' ...
        'non-zero rows.']);

elseif(numrows-num_zero_rows == 1)
	tc.verifyEqual(G_B, sum(G_A_expected, 2), 'AbsTol', 1e-13, ...
        ['Since we have only one non-zero row, the method under test ' ...
        'should sum all of the rows together in Step B.']);
end % if-elseif

end

