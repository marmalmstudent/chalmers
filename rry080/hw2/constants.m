R0_near = 10446;  % [m] Range distance to top edge of image
dR = 4;  % [m] Pixel size in range direction
dA = 0.43;  % [m] Pixel size in azimuth direction
theta_a = pi/60.0;  % degree Antenna beamwidth in azimuth direction
D_a = 1.3;  % [m] Antenna length in azimuth direction
lbda = 0.0566;  % Radar wavelength
v = 131;  % [ms^−1] Speed of antenna, e.g. the aircraft
PRF_v = 2.32;  % [m^−1] First guess