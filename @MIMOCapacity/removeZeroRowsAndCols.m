function [out, nz_rows, nz_cols] = removeZeroRowsAndCols(in)
% removeZeroRowsAndCols removes any rows or columns that contain all zeros
% from the provided matrix.  
%
% If no rows or no columns remain, returns a 1x1 matrix with a single zero
% (to maintain compatibility with methods that do not expect empty
% matrices).  
%
% By removing rows and columns that consist of only zero from the channel
% matrix, the computation and memory requirements (of computing and storing
% the PMF) can be greatly reduced in computing the PMF and the differential
% entropy of the recieved signal from the PMF without altering the computed
% differential entropy.  
%
%   [OUT, NZ_ROWS, NZ_COLS] = removeZeroRowsAndCols(IN)
%
% OUT is a 2D matrix that is IN with all rows and columns that contain only
%   zeros removed.  
% NZ_ROWS is a logical column vector indicating which rows have non-zero
%   values.  
% NZ_COLS is a logical row vector indicating which columns have non-zero
%   values.
%
% IN is a 2D matrix.  

if(~ismatrix(in))
    error('This method expects parameter IN to be a 2D matrix.');
end

% identify rows and columns with any non-zero elements
nz_rows = any(in~=0, 2);
nz_cols = any(in~=0, 1);

% output a matrix containing only those rows and columns from the input.
out = in(nz_rows, nz_cols);

% in case there are no rows or no columns left
if(isempty(out))
    out = 0;
end

end

