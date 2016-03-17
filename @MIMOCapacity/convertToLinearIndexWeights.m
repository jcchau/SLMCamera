function weights = convertToLinearIndexWeights(mat_size)
% convertToLinearIndexWeights computes the weights parameter for method
% convertToLinearIndex.
%
% WEIGHTS = convertToLinearIndexWeights(MAT_SIZE)
%
% WEIGHTS is the WEIGHTS parameter for method convertToLienarIndex.
%   WEIGHTS is a row-vector.  
% MAT_SIZE is the size of the matrix for which the indices should be
%   generated.  This should be the output of the size(...) method.  (A
%   row-vector).
%
% We separate the calculation of WEIGHTS to avoid repeatedly re-computing
% it every time we convert indices (since it remains unchanged as long as
% the matrix size remains unchanged).

if(any(mat_size == 0))
    error(['convertToLinearIndexWeights does not handle zero-sized ' ...
        'matrices.']);
end

weights = cumprod([1, mat_size(1:end-1)]);

end