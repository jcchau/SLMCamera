function msi = convertPointToSubscriptIndex(y, ymin, delta, nbins)
% convertPointToSubscriptIndex converts a point in multidimension space to
% a corresponding matrix subscript index.  
%
% In each dimension of y, index_y = ceil((y-y_min)/delta).  
% This algorithm was copied from
% MIMOCapacity.generateReceivedPmfForUniformInput, where it has been
% tested.
%
%   MSI = convertPointToSubscriptIndex(Y, YMIN, DELTA, NBINS)
%
% MSI is the matrix subscript index of each point in Y.  Each row in MSI is
%   a set of subscript indices for that point, and each row corresponds to
%   a point represented by a row in Y.
%
% Y is a matrix in which each row represents a point to convert.
% YMIN (row vector) is the point that corresponds to the zero matrix
%   subscript index. 
% DELTA (row vector) is the extent of a bin along each dimension of Y.  
% NBINS (row vector) is the size of the matrix being indexed.  Subscript
%   indices for rows of Y that are not covered by the matrix are discarded.
%
% Note that since index_y = ceil((y-y_min)/delta), if any dimension of y is
% exactly y_min, y would not land in any bin (the corresponding matrix
% subscript index would be exactly zero).  
%
% In this case, for convenience, if y(d)==y_min(d) for a point y for
% dimension d, then the matrix subscript index for that dimension is set to
% 1.  

% An alternate implementation could be:
% index_y = floor((y-y_min)/delta) + 1,
% but whether one implementation is used or the other should not have a
% significant effect for y with continuous PDFs.  

%% check input

if(~ismatrix(y))
    error('Parameter Y must be a matrix.');
end
% Note that ndims is the name of a MATLAB method, so call this something
% else.
num_dim_y = size(y, 2);

if(~isequal(size(ymin), [1, num_dim_y]))
    error(['Parameter YMIN should be a row vector with the same ' ...
        'number of columns as Y.']);
end

if(~isequal(size(delta), [1, num_dim_y]))
    error(['Parameter DELTA should be a row vector with the same ' ...
        'number of columns as Y.']);
end

if(~isequal(size(nbins), [1, num_dim_y]))
    error(['Parameter NBINS should be a row vector with the same ' ...
        'number of columns as Y.']);
end

%% do conversion

% Convert y into corresponding matrix-subscript indices for hits.
% In each dimension of y, index_y = ceil((y-y_min)/delta).  
msi = ceil(bsxfun(@rdivide, ...
    bsxfun(@minus, y, ymin), ...
    delta));

% Handle cases when y==y_min for any dimension.
% When indexing PMFs, this should not have any significant effect for
% continuous CDFs since the probability of y==y_min should be zero.  
% However, being able to set y==y_min is convenient.  
msi(bsxfun(@eq, y, ymin)) = 1;

% Discard any row of index_y for which any value is outside of the
% range of valid indices for matrix hits: [1, nbins(d)] for each dimension
% d. 
valid_msi = all(msi >= 1, 2) & ...
    all(msi <= repmat(nbins,size(msi,1),1), 2);
msi = msi(valid_msi, :);

end