function li = convertToLinearIndex(weights, subs)
% convertToLinearIndex converts a matrix subscripts into linear indices.
% Should provide the same results as sub2ind, but accepts the each set of
% matrix subscripts as a row in SUBS instead of requiring them to be comma
% separated parameters.
%
% LI = convertToLinearIndex(WEIGHTS, SUBS)
%
% LI is a row vector of corresponding linear indices.
% WEIGHTS is a row-vector of the weights used to convert SUBS into LI.
%   WEIGHTS = convertToLinearIndex(MAT_SIZE)
%   where MAT_SIZE is the matrix size (as a row vector) as specified by the
%   MATLAB size(...) method.  
% SUBS is a 2-D matrix, where each column corresponds to a dimension and
%   each row is a set of matrix subscript indices.  
%
% WARNING: this method is not designed to handle empty matrices.
%
% WARNING: This method assumes that subscript indices SUBS are valid (and
% for efficiency, does not check).  

% Reference:
% http://stackoverflow.com/questions/10146082/indexing-of-unknown-dimensional-matrix

% convert to zero-indexed
subs = subs - 1;

% Compute the linear index as a row vector.
% Add one to every linear index to convert from 0-indexing back to
% 1-indexing.
li = subs * weights' + 1;

end