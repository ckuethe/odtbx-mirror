function[Q] = Double_r_Accuracy(M32, M21, n, time)
% DOUBLE_R_ACCURACY Determination of the accurcacy of a computed orbit
% using the Double-r IOD method.
%
%    [Q] = Double_r_Accuracy(M32, M21, n, time) determines how well a 
%    computed orbit (using the double_r method) fits the observational
%    measurements by comparing the times between measurements according to 
%    the computed orbit, with the true measured times.
%
% INPUTS:
% M32 - Change in mean anomaly from observation 2 to 3    (rad)
% M21 - Change in mean anomaly from observation 1 to 2    (rad)
% n - Mean motion
% time - nx6 matrix of n datevectors for n observation times (datevector)
% 
% OUTPUTS:
% Q - Overall accuracy of the fit
% 
% CALLS:
% None
%
% WARNING: This accurcay calculation assumes that consecutive measurements
% are taken on the same pass

% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/19/2015                Original
 
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

% True time interval between measurements
t1 = etime(time(1,:),time(2,:));
t3 = etime(time(3,:),time(2,:));

%% Vallado's Accuracy method
F1 = t1 + M21/n;
F2 = t3 - M32/n;
Q = sqrt(F1^2 + F2^2);


