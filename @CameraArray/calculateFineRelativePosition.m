function [l, cos_theta, cos_phi] = ...
        calculateFineRelativePosition(obj, pixel_indices, ...
        transmitter_plane)
% calculateFineRelativePosition calculates the distance and
% relative angles between the pixels of the ImagingReceiver and
% the transmitter.  
%
% [L, COS_THETA, COS_PHI] = calculateFineRelativePosition( ...
%   OBJ, PIXEL_INDICES, TRANSMITTER_PLANE)
%
% OBJ is the ImagingReceiver object.
% PIXEL_INDICES are the MATLAB matrix indices of the pixels for
%   which the L, THETA, and PHI outputs should be calculated.
%   PIXEL_INDICES should be provided as a column vector
%   containing the linear indices of the pixel.  
%   The linear indices returned by find(a_tr>0) may be directly
%   used as this argument.  
%   The ind2sub method will be used to convert these linear
%   indices to convert PIXEL_INDICES to (row, column) indices.
%   sub2ind can be used to convert (row, column) indices to the
%   format required for PIXEL_INDICES.  
% TRANSMITTER_PLANE is the Plane object representing the plane
%   on which the transmitter polygon lies.  This can be
%   obtained by calling the getPlane method of the
%   transmitter's Polygon object.  
%
% L, COS_THETA, and COS_PHI are each column vectors where each
%   row corresponds to the result for the element specified by
%   the same row in PIXEL_INDICES.
%
% L is the distance between OBJ.lenspoint() and the projected
%   center of the specified pixel onto the TRANSMITTER_PLANE.  
% COS_THETA is the cosine of the angle between the normal of
%   the TRANSMITTER_PLANE and the line that connects the
%   lenspoint to the projected center of the specified pixel on
%   the TRANSMITTER_PLANE.  If TRANSMITTER_PLANE.normal points
%   away from OBJ.lenspoint, then the opposite normal vector
%   for the TRANSMITTER_PLANE will be used instead.  
% COS_PHI is the cosine of the angle between the direction the
%   camera is pointing and the line from the lenspoint to the
%   projected center of the specified pixel.  
%
% This method returns COS_THETA and COS_PHI instead of THETA
% and PHI because it saves a step to compute and because it's
% what we'll use to compute the gains from each transmitter to
% each pixel anyway.  
%
% See also find, ind2sub, sub2ind,
% CameraArray.calculateTransmitterAreaReceived.

[row, col] = ind2sub( ...
    [obj.pixel_array.nrows, obj.pixel_array.ncols], ...
    pixel_indices);

pixel_center = obj.pixel_array.getElementCenter(row, col);

replicated_lenspoint = ...
    repmat(obj.lenspoint(), size(pixel_center, 1), 1);

projection_direction = replicated_lenspoint - pixel_center;

[projected_center, ray_intersects, is_parallel] = ...
    transmitter_plane.intersectRay( ...
    replicated_lenspoint, projection_direction);

% check for ray projection errors
if(~all(ray_intersects))
    error('CameraArray:BadPixelProjection', ...
        'The pixel projection misses the transmitter plane.');
    % Error here because the returned values for L, COS_THETA,
    % and COS_PHI would be incorrect/invalid.  
end
if(any(is_parallel))
    warning('CameraArray:ParallelPixelProjection', ...
        'The pixel projection ray is parallel to TRANSMITTER_PLANE.  L may be incorrect.');
    % Used a warning here because cos_theta would be nearly
    % zero for this pixel, resulting in (nearly) zero gain
    % regardless of what value of L is computed.
end

% the scalar distance from lenspoint to the projected center
l = sqrt(sum((projected_center - replicated_lenspoint).^2, 2));

% normalize projection_direction so that we can use the
% dot-product to determine the cosine of theta and phi.
magnitude_projection_direction = ...
    sqrt(sum(projection_direction.^2, 2));
normalized_projection_direction = projection_direction ./ ...
    repmat(magnitude_projection_direction, 1, 3);

% Since all vectors below are normalized, the dot-product of
% two vectors below is equal to the cosine of the angle between
% them.  
% Take the absolute value of this cosine to effectively ensure
% that we're choosing the smaller angle (since each plane can
% be defined with either normal vector pointing in opposite
% directions).
cos_theta = abs(normalized_projection_direction * ...
    transmitter_plane.normal()');
cos_phi = abs(normalized_projection_direction * ...
    obj.pixel_array.plane_axes.normal()');

end
