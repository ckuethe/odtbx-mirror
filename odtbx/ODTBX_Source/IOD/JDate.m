function [JD] = JDate(time)
% JDATE Julian Date calculation
%
%    [JD] = JDate(time) Calculates the Julian Date for a given Gregorian
%    calendar date. Use the function juliandate from the Aerospace Toolbox
%    in MATLAB if available. This function serves as a backup in case the
%    toolbox cannot be accessed.
%
% INPUTS:
% time - nx6 matrix of n datevectors for n observation times (datevector)
%
% OUTPUTS:
% JD - 1xn vector of Julian Dates
%
% CALLS:
% None

% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    07/27/2015                Original

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

for i = 1:size(time,1)
    yr(i) = time(i,1);
    mo(i) = time(i,2);
    d(i) = time(i,3);
    h(i) = time(i,4);
    min(i) = time(i,5);
    s(i) = time(i,6);
    JD(i) = yr(i)*367 - floor((7*(yr(i) + floor((mo(i)+9)/12)))/4) + ...
        floor(275*mo(i)/9) + d(i) + 1721013.5 + ((s(i)/60 + min(i))/60 + h(i)) / 24;
end

end