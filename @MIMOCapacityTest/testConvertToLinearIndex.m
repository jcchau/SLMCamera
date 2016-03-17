function testConvertToLinearIndex(tc)
    % Compares the result of convertToLinearIndex against sub2ind with a
    % random-dimension random-size matrix, and a random number of random
    % indices to convert.
    % Simultaneously tests method convertToLinearIndexWeights.

    % up to 4 dimensions
    num_dims = randi(4);

    % each dimension up to 11 long
    mat_size = randi(11, 1, num_dims);

    % up to 5 indices
    num_indices = randi(5);

    % the indices
    subs_cell = cell(1, num_dims);
    subs_mat = zeros(num_indices, num_dims);
    for ii = 1:num_dims
        subs_cell{ii} = randi(mat_size(ii), num_indices, 1);
        subs_mat(:, ii) = subs_cell{ii};
    end % for ii = 1:dimensions

    weights = MIMOCapacity.convertToLinearIndexWeights( ...
        mat_size);
    li_test = MIMOCapacity.convertToLinearIndex( ...
        weights, subs_mat);

    if(num_dims == 1)
        % sub2ind doesn't work for 1D vectors; in this case, the
        % matrix subscript indexing is equal to the linear
        % indexing.
        li_expected = subs_cell{1};
    else
        li_expected = sub2ind(mat_size, subs_cell{:});
    end

    tc.verifyTrue(isequal(li_test, li_expected), ...
        ['convertToLinearIndex did not return the same ' ...
        'results as sub2ind.'])
end

