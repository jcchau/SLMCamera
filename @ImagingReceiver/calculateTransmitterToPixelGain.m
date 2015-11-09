function gains = calculateTransmitterToPixelGain(obj, transmitter_polygon)
% calculateTransmitterToPixelGain calculates the optical channel gain from
% a diffuse (Lambertian) transmitter to each pixel of the ImagingReceiver.
%
% GAINS = calculateTransmitterToPixelGain(OBJ, TRANSMITTER_POLYGON)
%
% OBJ is the ImagingReceiver object.
% TRANSMITTER_POLYGON is Polygon object representing the transmitter's
%   emitting surface.
%
% GAINS is an OBJ.pixel_array.nrows row by OBJ.pixel_aray.ncols column
%   matrix.  Element GAINS(row,col) contains the calculated optical gain
%   from the transmitter to the pixel at (row,col).
%
% This method uses the gain formula on p. 49 of Book 3 of my lab book,
% "Imaging Receivers & Photodetector Arrays":
%  P_r/P_t = (A_tr*A_r)/(pi*A_t * l^2) * cos(theta) * cos(phi) = gain

if(~isa(obj, 'ImagingReceiver'))
    error('OBJ must be an ImagingReceiver object.');
end
if(~isa(transmitter_polygon, 'Polygon'))
    error('transmitter_polygon must be a Polygon object.');
end

% area of the transmitter polygon
a_t = transmitter_polygon.area();

% area of the receiver's lens/aperture
a_r = pi * obj.r_aperture^2;

% c_precomputed is defined as below and is precomputed once for efficiency.
c_precomputed = a_r/(pi*a_t);

%% calculate the remaining factors in the gain formula

a_tr = obj.calculateTransmitterAreaReceived(transmitter_polygon);
% !!!QUIRK!!!
% When the argument for method find is a row vector, find returns a row
% vector instead of the regular column vector.
% When indices_receiving_pixels is a row vector, it triggers an error in
% obj.calculateFineRelativePosition.
% Furthermore, even when indices_receiving_pixels is converted to a column
% vector, the fact that a_tr is a row vector triggers an error on the last
% command to calculate the variable gains.
% To work around this problem, convert a_tr to a column vector if it is a
% row vector.
if(isrow(a_tr))
    a_tr = a_tr(:);
end

indices_receiving_pixels = find(a_tr > 0);

[l_squared, cos_theta, cos_phi] = obj.calculateFineRelativePosition( ...
    indices_receiving_pixels, transmitter_polygon.getPlane());


matrix_size = [obj.pixel_array.nrows, obj.pixel_array.ncols];

% compute gains
gains = zeros(matrix_size); % default value for pixels where (a_tr == 0).
gains(indices_receiving_pixels) = c_precomputed .* ...
    a_tr(indices_receiving_pixels) ./ l_squared .* cos_theta .* cos_phi;

end

