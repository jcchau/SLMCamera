function a_tr = calculateTransmitterAreaReceived(obj, transmitter_polygon)
% CALCULATETRANSMITTERAREARECEIVED Calculates the area of the transmitter
% that is within the angle-of-view of each of the CameraArray's pixels.
%
% A_TR = calculateTransmitterAreaReceived(OBJ, TRANSMITTER_POLYGON)
%
% TRANSMITTER_POLYGON is a Polygon object representing the transmitter.
% A_TR is a NROWS by NCOLS matrix, where each element of the matrix
%   represent the area in the TRANSMITTER_POLYGON that is seen by the
%   corresponding pixel in the CameraArray.  (To clarify, the fraction of
%   the TRANSMITTER_POLYGON seen by the camera pixels would be
%   A_TR/TRANSMITTER_POLYGON.area().

% pre-clip the transmitter_polygon to avoid tracing parallel and
% non-intersecting rays

[transmitter_polygon, isvalid] = obj.preclipPolygon(transmitter_polygon);

% Default values for the returned a_tr
a_tr = zeros(obj.pixel_array.nrows, obj.pixel_array.ncols);
% those pixels that aren't listed in pixel_indices don't see the
% transmitter (and hence, see zero area of the transmitter).

if(isvalid)
    % The following operations rely on the assumption that
    % transmitter_polygon is a valid polygon that's entirely in front of
    % the camera's lenspoint.  Otherwise, a_tr with all zeros is already
    % correct.  

    polygon_matrix_tx = transmitter_polygon.toMatrix();
    nvertices = size(polygon_matrix_tx, 1);

    ray_direction = polygon_matrix_tx - ...
        repmat(obj.lenspoint(), nvertices, 1);

    image_plane = Plane(obj.pixel_array.centerpoint, ...
        obj.pixel_array.plane_axes.normal());

    % we don't need the returned ray_intersects or is_parallel values since
    % the vertices of the clipped polygon all project onto the image_plane.
    [polygon_matrix_image, ~, ~] = ...
        image_plane.intersectRay(polygon_matrix_tx, ray_direction);

    % get a candidate list of all CameraArray pixels that could possibly
    % intersect with image_polygon.
    % Indices are (row,column).  
    polygon_image = Polygon(polygon_matrix_image);
    pixel_indices = ...
        obj.pixel_array.listIntersectingElements(polygon_image);

    ncandidates = size(pixel_indices, 1);

    %% determine the area of the transmitter seen by each pixel
    
    % Determine on what plane does the transmitter polygon lie.
    % And for transforming projected pixel polygons onto the same
    % coordinate system, create an OrthogonalAxes object.
    [a_tx, b_tx, origin_tx, axisA_tx, axisB_tx] = ...
        transmitter_polygon.to2D();
    tx_axes = OrthogonalAxes(axisA_tx, axisB_tx);
    plane_tx = Plane(origin_tx, tx_axes.normal());
    
    % reorder a_tx and b_tx so that the vertices are clockwise (to avoid
    % the warning from polybool).
    [a_tx, b_tx] = poly2cw(a_tx, b_tx);
    
    % determine whether we want to keep vertices in the + or -
    % plane_tx.normal() direction.
    tx_to_lenspoint_displacement_along_normal = ...
        dot(obj.lenspoint() - origin_tx, plane_tx.normal());
    
    if(tx_to_lenspoint_displacement_along_normal ~= 0)
        % If (tx_to_lenspoint_displacement_along_normal == 0), then a_tr
        % should be 0 as it's already set.  
        
        if(tx_to_lenspoint_displacement_along_normal > 0)
            % lenspoint is on the +normal side of plane_tx
            clip_plane = Plane(obj.lenspoint, plane_tx.normal());
        else
            % lenspoint is on the -normal side of plane_tx
            clip_plane = Plane(obj.lenspoint, -plane_tx.normal());
        end
        
        % For each candidate pixel polygon,
        %   - clip the pixel polygon to remove any vertices behind the
        %   plane parallel to the transmitter polygon that goes through
        %   lenspoint.
        %   - project the clipped pixel polygon through the lenspoint onto
        %   plane_tx.
        %   - intersect the projected pixel polygon with the transmitter
        %   polygon and set the corresponding a_tr element to the area of
        %   this intersection. 
        for ii = 1:ncandidates
            pixel_polygon = obj.pixel_array.getPolygon(...
                pixel_indices(ii,1), pixel_indices(ii,2));
            
            [pixel_polygon, isvalid] = clip_plane.clipPolygon(...
                pixel_polygon);

            if(isvalid)
                % If ~isvalid, then
                % a_tr[pixel_indices(ii,1), % pixel_indices(ii,2)] = 0.
                
                % project the pixel back to plane_tx
                ray_points = pixel_polygon.toMatrix();
                ray_directions = ...
                    repmat(obj.lenspoint(), size(ray_points, 1), 1) - ...
                    ray_points;
                
                pixel_on_plane_tx = plane_tx.intersectRay(...
                    ray_points, ray_directions);
                
                % Need to convert the projected pixel polygon to 2D on the
                % transmitter's plane to do polygon intersection. 
                % (We already have the transmitter's polygon as 2D on this
                % plane.)
                [a_rx, b_rx, ~] = tx_axes.transformPointsFromGlobal(...
                    origin_tx, pixel_on_plane_tx);
                % re-order polygon vertices for clockwise order to avoid
                % warning from polybool
                [a_rx, b_rx] = poly2cw(a_rx, b_rx);
                
                % find the intersection
                [a_intersection, b_intersection] = ...
                    polybool('intersection', a_rx, b_rx, a_tx, b_tx);

                % finally, determine the area
                % where pixel_indices(ii,1) is the row of the pixel and
                % pixel_indices(ii,2) is the column of the pixel.
                a_tr(pixel_indices(ii,1), pixel_indices(ii,2)) = ...
                    polyarea(a_intersection, b_intersection);
                
            end % if(isvalid)
        end % for ii
    
    end % if(tx_to_lenspoint_displacement_along_normal ~= 0)
    
end % if(isvalid)

end
