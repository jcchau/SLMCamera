function polyout = preclipPolygon(obj, polyin)
% PRECLIPPOLYGON Pre-clips a polygon, ensuring that the polygon remains
% within clip_half_angle.  
%
% This function is useful before projecting a polygon through the lenspoint
% onto the image plane to ensure that all remaining vertices form rays
% through the lenspoint that intersect the image plane.  
%
% POLYOUT = preclipPolygon(OBJ, POLYIN, CLIP_HALF_ANGLE)
%
% OBJ is the CameraArray object.
% POLYIN is the polygon representing the object in the scene to be imaged.
% POLYOUT is the clipped polygon.

% determine the plane to use to clip the polygon
% clip at 1e-3 the lens_to_array_distance in front of the lenspoint.
clip_plane_normal = obj.pixel_array.plane_axes.normal();
clip_plane_point = obj.pixel_array.centerpoint + ...
    1.001 * obj.lens_to_array_distance * ...
    clip_plane_normal;
clip_plane = Plane(clip_plane_point, clip_plane_normal);

% then perform the clipping
polyout = clip_plane.clipPolygon(polyin);

end

