function detyz=bctodetyz(geometry, bc)

% This function takes beam center as an input and gives detector
% position as an output. The detector position is the xyz coordinates
% of the detector center in mm

horzpx=geometry.detector.det_pixels_horz;
vertpx=geometry.detector.det_pixels_vert;

horzmm=geometry.detector.det_size_horz;
vertmm=geometry.detector.det_size_vert;

horz_offset=round(horzpx/2);
vert_offset=round(vertpx/2);

pxpermm_horz=horzpx/horzmm;
pxpermm_vert=vertpx/vertmm;

detyz=[(bc(1)-horz_offset)/pxpermm_horz (bc(2)-vert_offset)/pxpermm_vert];

end