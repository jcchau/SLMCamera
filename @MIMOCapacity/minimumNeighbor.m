function out = minimumNeighbor(in, num_dimensions)
% minimumNeighbor returns a matrix OUT (of the same size as input matrix
% IN) that contains for each cell, the minimum of the values of that cell
% or its adjacent neighbors (including neighbors that only touch by a
% corner).  
%
% The input matrix IN is assumed to only contain non-negative real numbers,
% and the outside of matrix IN is assumed to consist of zeros (i.e., matrix
% IN is surrounded by zeros).  
%
%   OUT = minimumNeighbor(IN, NUM_DIMENSIONS)
%
% IN is a arbitrary-dimension matrix of non-negative real numbers.  
% NUM_DIMENSIONS is the number of dimensions of matrix IN.  This parameter
%   is necessary because MATLAB is otherwise unable to distinguish between
%   a matrix of a certain size, and the matrix with trailing singleton
%   dimensions.  This method needs to treat these cases differently since
%   the first and last elements of each dimension are set to zero.  
%
% OUT is a non-negative matrix if the same dimensions as IN, consisting of
%   the minimum values of each cell's (in matrix IN) neighbors (including
%   the cell itself it it contains the minimum value).  

if(isempty(in))
    % If the input is empty, return an empty matrix (with the same number
    % of dimensions to avoid causing potential problems).
    out = in;
else

    %% set out equal in
    % So we won't need to code the first dimension differently to keep a
    % matrix cell's value the same if its value is the smallest of its
    % neighbors.

    out = in;

    %% update along each dimension to get the smallest neighboring value

    % loop through each dimension
    for d = 1:num_dimensions

        if(num_dimensions == 1)
            % Special case for row vectors because size(in,d) is always 1
            % (instead of the only dimension).  
            last_index = length(in);
        else
            last_index = size(in,d);
        end % end if-else
        
        % Use a cell-array as the index to enable vectorization of this
        % algorithm.  
        % This matindex selects everything along every dimension. Inside
        % the for ii loop below, it replaces the index for one of the
        % dimensions with a specific value to select a "(hyper)plane" of
        % bins.
        % References:
        % - http://www.mathworks.com/matlabcentral/answers/58825-how-can-i-dynamically-assign-access-elements-in-a-matrix-with-arbitrary-dimensions
        % - https://www.mathworks.com/matlabcentral/newsreader/view_thread/239004

        matindex = repmat({':'}, 1, num_dimensions);

        % In dimension d, select indices 2:(end-1).
        % Since the "end" keyword only works inside an index expression
        % (and since we're creating this cell of indices outside of an
        % index expression), we compute the equivalent value of "end" using
        % method size.
        matindex(d) = {2:last_index-1};

        % Also create a matrix index that selects indices 1:(end-2) in
        % dimension d.
        matindex_prev = repmat({':'}, 1, num_dimensions);
        matindex_prev(d) = {1:last_index-2};

        % And a matrix that selects indices 3:end in dimension d.
        matindex_next = repmat({':'}, 1, num_dimensions);
        matindex_next(d) = {3:last_index};

        % Set out to the minimum of either its own or its neighbor's (in
        % dimension d) value.  
        % Using out for the neighbor instead of matrix "in" allows us to
        % in effect also consider cells that only share a corner as
        % neighbors.
        out(matindex{:}) = min(out(matindex{:}), out(matindex_prev{:}));
        out(matindex{:}) = min(out(matindex{:}), out(matindex_next{:}));

    end % for d = 1:ndims

    %% set the border cells of out to zero
    % Since we assume that values are zero outside of the border of out,
    % and zero is the minimum possible value (assuming that matrix "in"
    % only contains non-negative real numbers), the "minimum neighbor"
    % value of each cell along the border is zero.

    for d = 1:num_dimensions
        
        if(num_dimensions == 1)
            % Special case for row vectors because size(in,d) is always 1
            % (instead of the only dimension).  
            last_index = length(in);
        else
            last_index = size(in,d);
        end % end if-else
        
        matindex = repmat({':'}, 1, num_dimensions);
        matindex(d) = {[1, last_index]};

        out(matindex{:}) = 0;
        
    end % for d = 1:num_dimensions

end % else (of if(isempty(in)))

end

