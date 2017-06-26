function [LOS] = LOS_Vectors(meas)
% LOS_VECTORS Calculation of Line-Of-Site Vectors
% 
%    [LOS] = LOS_Vectors(meas) Returns the line-of-site unit vectors 
%    given topocentric right ascension and declination observations.
%
% INPUTS:
% meas - nx2 matrix containing n observations with the first column holding
%        right ascention values and the second column holding declination
%        values [RA Dec] (degrees)
%
% OUTPUTS:
% LOS - 3xn matrix of LOS vectors with each column, n, corresponding to
%       each observation given in the input
%
% CALLS:
% None

% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/03/2015                Original

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

LOS = [];
for i = 1:length(meas)
    RA = meas(i, 1);
    dec = meas(i, 2);
    LOS = [LOS; [cosd(dec)*cosd(RA) cosd(dec)*sind(RA) sind(dec)]];
end

LOS = LOS.';

end

