function out = maximumNeighbor(in)
% maximumNeighbor returns a matrix of the largest value of the neighboring
% matrix elements.  
%
% This method assumes that for the matrix elements at the outer boundary of
% the matrix, that the "outside" of the matrix does not contain larger
% elements.  
%
% This method is used by MIMOCapacity.calculateMinimumDiffEntropyFromPmf;
% however, this method does not insert the peak PMF or PDF value at the
% center of the distribution; calculateMinimumDiffEntropyFromPmf should do
% that separately.  
%
%   OUT = maximumNeighbor(IN)
%
% OUT is a matrix (of the same size as IN) that contains the maximum value
%   of the immediately adjacent matrix elements (including elements that
%   only touch by a corner) in matrix IN; or if the corresponding matrix
%   element in IN is largest, that value.  
%
% IN is an arbitrary-dimension matrix of real numbers.  

num_dimensions = ndims(in);

% start with out = in
% This way, we won't have to use different code for the first dimension.  
out = in;

%% Compute the maximum neighboring value 

for d = 1:num_dimensions
    
    last_index = size(in, d);
    
    %% For the elements at the end along dimension d
    % Do this before the elements in the middle along dimension d. 
    % Otherwise, matrix out in index 2 along dimension d, may have already
    % copied a larger value from index 3 along dimension d.  And then if we
    % copy the maximum value between index 1 and 2 along dimension d into
    % index 1 along dimension d, then we may end up with some values from
    % index 3 (as the maximum neighbor for elements in index 1).  
    % This problem does not occur the other way around since the elements
    % at the end don't have elements further beyond the end that may
    % erroneously copied (e.g., don't need to worry about values from index
    % 0 along dimension d being copied into index 2 along dimension d).  
    
    if(last_index >= 2) % last_index is the length along dimension d
        
        % Matrix indices for index 1 and 2 along dimension d
        matindex_end = repmat({':'}, 1, num_dimensions);
        matindex_end(d) = {1};
        matindex_adjacent = repmat({':'}, 1, num_dimensions);
        matindex_adjacent(d) = {2};
        
        % Store in index 1 along dimension d of out the maximum of the
        % values of either index 1 or 2 along dimension d. 
        out(matindex_end{:}) = max(out(matindex_end{:}), ...
            out(matindex_adjacent{:}));
        
        % Matrix indices for the last and second-to-last indices along
        % dimension d.  
        matindex_end(d) = {last_index};
        matindex_adjacent(d) = {last_index-1};
        
        % Store in index "end" along dimension d of matrix out the maximum
        % of the values of either index "end" or "end-1" along dimension d.
        out(matindex_end{:}) = max(out(matindex_end{:}), ...
            out(matindex_adjacent{:}));
        
    end % if(last_index >= 2)
    
    %% For the elements in the middle along dimension d
    
    % Matrix index for everything except the first and last along dimension
    % d.
    matindex = repmat({':'}, 1, num_dimensions);
    matindex(d) = {2:last_index-1};
    
    % Index for everything one step before matindex along dimension d
    matindex_prev = repmat({':'}, 1, num_dimensions);
    matindex_prev(d) = {1:last_index-2};
    
    % Index for everything one step after matindex along dimension d
    matindex_next = repmat({':'}, 1, num_dimensions);
    matindex_next(d) = {3:last_index};
    
    % Note that in the above index calculations, no elements would be
    % selected if size(in,d)<3.  k
    
    out(matindex{:}) = max(out(matindex{:}), out(matindex_prev{:}));
    out(matindex{:}) = max(out(matindex{:}), out(matindex_next{:}));
        
end % for d = 1:num_dimensions

end

