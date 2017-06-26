function [] = FileFormatter_Angles(epoch, times, meas, sigma, long, lat, alt, filename)
% FILEFORMATTER_ANGLES Write angle observations to data file with formatting
%
%    [filename] = FileFormatter_Angles(epoch, times, meas, sigma, lat, long, alt)
%    Takes the appropriate angles-only observational data and writes it to
%    a data file following the correct file formatting procedure outlined
%    in the IOD Application documentation. This is provided for convenience
%    only and is not required for the application's operation.
% 
% INPUTS:
% epoch - Datenumber of epoch for observations
% times - Array of datenumbers to be added to the epoch. These represent
%         the exact observation times.
% meas - Angle observations in a nx2 matrix [RA Declination] (degrees)
% sigma - Matrix of uncertainty values for observations. nx3 where the
%         first column is RA sigma (degrees), second column is Dec sigma 
%         (degrees) and the third column is the correlation between the 
%         RA and Dec values.
% long - Array of longitude or RA values for the observer's location at
%        each observation time.  (Degrees - Negative for West, Positive for East)
% lat - Array of latitude or Dec values for the observer's location at each
%       observation time. (Degrees - Negative for South, Positive for North)
% alt - Altitude (range) values for observer's location at each observation
%       time. (km)
% filename - A string of a filename that the data should be written to.
%            Must be .txt or .dat
%
% OUTPUTS:
% None
%
% CALLS:
% None

% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    07/15/2015                Original

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

tspan = times + epoch;
tspan = datevec(tspan);

fid = fopen(filename, 'w');
for i = 1:length(lat)
    fprintf(fid, '%d %0.2d %0.2d ',tspan(i,1:3));
    fprintf(fid, '%0.2d:%0.2d:%0.3f ',tspan(i,4:6));
    fprintf(fid, '%f %0.6e %f %0.6e %f ', meas(i,1), sigma(i,1), meas(i,2), sigma(i,2), sigma(i,3));
    fprintf(fid, '%0.10f %f %f\n', lat(i), long(i), alt(i));
end
fclose(fid);
end