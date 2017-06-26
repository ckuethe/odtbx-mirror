function [site_Pos] = Site_Position(lat_GD, theta_LST, h_ELLP)
% SITE_POSITION Calculation of observation site position vector
%
%    [site_Pos] = Site_Position(lat_GD, theta_LST, h_ELLP) Calculates an 
%    observation site's position vector in the Geocentric Equatorial 
%    Coordinate System, IJK. This can act as an approximation for the ECI 
%    coordinate vectors. If the site longitude is input instead of the 
%    local sidereal time, the site position vector in the International 
%    Terrestrial Reference Frame (ITRF) frame is given instead.
%
% INPUTS:
% lat_GD - Geodetic Latitude (degrees)
% theta_LST - Either the local sidereal time or the longitude (degrees)
% h_ELLP - Altitude of observation site (km)
%
% OUTPUTS:
% site_Pos - The position vector of the observation site (km)
%
% CALLS:
% None

% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/02/2015                Original
%   Ryan Willmot    07/02/2015                Added sun functionality

% ODTBX: Orbit Determination Toolbox
% 
% Copyright (c) 2003-2015 United States Government as represented by the
% administrator of the National Aeronautics and Space Administration. All
% Other Rights Reserved.
% 
% This file is distributed "as is", without any warranty, as part of the
% ODTBX. ODTBX is free software; you can redistribute it and/or modify it
% under the terms of the NASA Open Source Agreement, version 1.3 or later.
% 
% You should have received a copy of the NASA Open Source Agreement along
% with this program (in a file named License.txt); if not, write to the 
% NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.

global radius
global Type

% Case for Earth-centered orbit
% Radius of curvature in the meridian
if strcmpi(Type,'Earth')
    e = 0.081819221456;    %Eccentricity of Earth
    c_E = radius / sqrt(1 - (e^2)*(sind(lat_GD)^2));    
    s_E = c_E * (1 - e^2);

% Case of Sun-centered orbit
else
    c_E = 0;
    s_E = 0;
end

r_d = (c_E + h_ELLP)*cosd(lat_GD);    %(km)
r_k = (s_E + h_ELLP)*sind(lat_GD);    %(km)

site_Pos = [r_d*cosd(theta_LST); r_d*sind(theta_LST); r_k];    %(km)

end