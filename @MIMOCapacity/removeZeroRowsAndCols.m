function out = removeZeroRowsAndCols(in)
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
%   OUT = removeZeroRowsAndCols(IN)
%
% OUT is a 2D matrix that is IN with all rows and columns that contain only
%   zeros removed.  
%
% IN is a 2D matrix.  

if(~ismatrix(in))
    error('This method expects parameter IN to be a 2D matrix.');
end

%% remove any rows of zeros

% pre-allocate space for the matrix without zero rows
in_without_zero_rows = zeros(size(in));

% counter for number of non-zero rows
nrows = 0;

for r = 1:size(in, 1)
    if(any(in(r, :) ~= 0))
        nrows = nrows + 1;
        in_without_zero_rows(nrows, :) = in(r, :);
    end
end

% trim any extra pre-allocated rows
in_without_zero_rows = in_without_zero_rows(1:nrows, :);

%% remove any columns of zeros

% pre-allocate space for the matrix without zero columns
out = zeros(size(in_without_zero_rows));

% counter for number of non-zero columns
ncols = 0;

for c = 1:size(in_without_zero_rows, 2)
    if(any(in_without_zero_rows(:, c) ~= 0))
        ncols = ncols + 1;
        out(:, ncols) = in_without_zero_rows(:, c);
    end
end

% trim any extra pre-allocated columns
out = out(:, 1:ncols);

%% in case there are no rows or no columns left

if(nrows == 0 || ncols == 0)
    out = 0;
end

end

