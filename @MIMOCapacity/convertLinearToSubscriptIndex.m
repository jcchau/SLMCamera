function msi = convertLinearToSubscriptIndex(weights, li)
% convertLinearToSubscriptIndex converts linear indices to matrix subscript
% indices given the size of the matrix.
%
%   MSI = convertLinearToSubscriptIndex(MATSIZE, LI)
%
% MSI is the corresponding matrix subscript index (as a row vector) for the
%   given linear index.  
%
% WEIGHTS is a row vector of the weights computed by
%   convertToLinearIndexWeights using the size of the matrix to be indexed.
%   (WEIGHTS, as generated by convertToLinearIndexWeights has as many
%   elements as the resulting MSI).  
% LI is the linear index of the element to be indexed.  
%
% This method should return the same matrix subscript index as ind2sub.
% However, this method returns the matrix subscript index as a row vector
% instead of separate comma-separated outputs. 
%
% Used the MATLAB implementation of ind2sub as a reference while
% implementing this method.

% number of dimensions of the matrix
nd = length(weights);

% preallocate the output
msi = zeros(1, nd);

% change li from 1-indexed to 0-indexed;
li = li - 1;

% Work from the last to the first dimension
for d = nd:-1:1
    
    % compute the zero-indexed linear index for all of the previous
    % dimensions (as if we didn't have dimension d).
    lipd = rem(li, weights(d));
    
    % The zero-indexed subscript index for dimension d is: 
    % floor(li / weights(d)).
    msi(d) = (li - lipd) / weights(d);
    
    % And then repeat with the next iteration of the loop, using the linear
    % index with dimension d removed.  
    li = lipd;
    
end % for d

% convert the 0-indexed subscript indices back to 1-indexed subscript
% indices.
msi = msi + 1;

end

