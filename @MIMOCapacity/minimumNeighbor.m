function out = minimumNeighbor(in)
% minimumNeighbor returns a matrix OUT (of the same size as input matrix
% IN) that contains for each cell, the minimum of the values of that cell
% or its adjacent neighbors (including neighbors that only touch by a
% corner).  
%
% The input matrix IN is assumed to only contain non-negative real numbers,
% and the outside of matrix IN is assumed to consist of zeros (i.e., matrix
% IN is surrounded by zeros).  
%
%   OUT = minimumNeighbor(IN)
%
% IN is a arbitrary-dimension matrix of non-negative real numbers.  
% OUT is a non-negative matrix if the same dimensions as IN, consisting of
% the minimum values of each cell's (in matrix IN) neighbors (including the
% cell itself it it contains the minimum value).  

if(isempty(in))
    out = [];
else

    %% set out equal in
    % So we won't need to code the first dimension differently to keep a
    % matrix cell's value the same if its value is the smallest of its
    % neighbors.

    out = in;

    %% update out along each dimension to get the smallest neighboring value

    num_dimensions = ndims(in);

    % loop through each dimension
    for d = 1:num_dimensions

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
        matindex(d) = {2:(size(in,d)-1)};

        % Also create a matrix index that selects indices 1:(end-2) in
        % dimension d.
        matindex_prev = repmat({':'}, 1, num_dimensions);
        matindex_prev(d) = {1:(size(in,d)-2)};

        % And a matrix that selects indices 3:end in dimension d.
        matindex_next = repmat({':'}, 1, num_dimensions);
        matindex_next(d) = {3:size(in,d)};

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
        matindex = repmat({':'}, 1, num_dimensions);
        matindex(d) = {[1, size(in,d)]};

        out(matindex{:}) = 0;
    end % for d = 1:num_dimensions

end % else (of if(isempty(in)))

end

