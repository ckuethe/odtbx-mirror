function [theta_LST] = LST(date, long)
% LST Local Sidereal Time Calculation
%
%    [theta_LST] = LST(date,long) Takes in a gregorian calendar date, 
%    time, and longitude of a location and returns the local sideral 
%    time in degrees (theta_LST). 
%
% INPUTS:
% date - Datevector of observation (UT1)
% long - longitude (Degrees - negative for West, positive for East)
%
% OUTPUTS:
% theta_LST - Local Sidereal Time (degrees)
%
% CALLS:
% None

% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/02/2015                Original

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

jd = JDate(date);
t_UT1 = (jd - 2451545.0) / 36525;

theta_GMST = 67310.54841 + (876600*3600 + 8640184.812866)*t_UT1...
    + 0.093104*(t_UT1^2) - (6.5*10^-6)*(t_UT1^3);

while abs(theta_GMST) > 86400
    if theta_GMST > 0
        theta_GMST = theta_GMST - 86400;
    end
    if theta_GMST < 0
        theta_GMST = theta_GMST + 86400;
    end
end

% Conversion from seconds to degrees
theta_GMST = theta_GMST / 240;

if theta_GMST < 0
    theta_GMST = 360 + theta_GMST;
end

theta_LST = theta_GMST + long;
if theta_LST < 0
    theta_LST = 360 + theta_LST;
elseif theta_LST > 360
    theta_LST = theta_LST - 360;
else
end
