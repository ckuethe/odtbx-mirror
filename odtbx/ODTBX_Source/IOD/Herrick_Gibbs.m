function [v2] = Herrick_Gibbs(r, time)
% HERRICK_GIBBS Herrick-Gibbs method of velocity determination given three postion vectors
% 
%    [v2] = Herrick_Gibbs(r, time)Uses the Herrick-Gibbs method to return 
%    an estimation of the velocity vector for the middle observation time 
%    out of three observations. 
%
% This method is best used with position vectors from a single pass of the
% object over a particular ground station.
%
% INPUTS:
% r - 3x3 matrix containing orbiting object's position vectors at three
%     sequential observation times. (km)
% time - 3x6 matrix of datevectors for the three observation times
%
% OUTPUTS:
% v2 - Estimation of object's velocity vector at second time (km/s)
%
% CALLS:
% None
%
% WARNING: To use the Herrick-Gibbs method, check that the position vectors are
% coplanar to a certain tolerance (<3deg) and that the separation between
% the vectors is small (< 5deg). Failure to do this will cause erroneous
% estimations. This function assumes that these checks are done before
% passing arguments through, but will check again and provide error
% messages to indicate potential errors.

% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/08/2015                Original

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

% Gravitational Parameter
global muglobal;
mu = muglobal;   % (km^3/s^2)

r1 = r(:,1);
r2 = r(:,2);
r3 = r(:,3);

Z12 = cross(r1, r2);
Z23 = cross(r2, r3);
Z31 = cross(r3, r1);

% Check that position vectors are within coplanar tolerance (~3deg)
a_COP = 90 - acosd(dot(Z23,r(:,1))/(norm(Z23)*norm(r(:,1))));
if abs(a_COP) > 3
    fprintf(1,['\nWARNING: Position Vectors are not coplanar to 3 degree', ... 
    ' tolerance!\nEstimations  using Gibbs or Herrick-Gibbs', ...
    ' may be off by a large amount.\n']); 
    fprintf(1,'Coplanar angle is: %f deg\n',a_COP); 
end

% Check that angle between position vectors is small enough
alpha_12 = acosd(dot(r(:,1),r(:,2))/(norm(r(:,1))*norm(r(:,2))));
alpha_23 = acosd(dot(r(:,2),r(:,3))/(norm(r(:,2))*norm(r(:,3))));

if (abs(alpha_12) > 5 && abs(alpha_23) > 5)
    fprintf(1, ['\nWARNING: Angle between position vectors is greater than',...
        ' 3 degrees. This will result in erroneous results when using',...
        ' Herrick-Gibbs Method!\n']);
    fprintf(1,'Angle between first and second vector is: %f deg\n',alpha_12);
    fprintf(1,'Angle between second and third vector is: %f deg\n',alpha_23);
else
end

t21 = etime(time(2,:), time(1,:));
t32 = etime(time(3,:), time(2,:));
t31 = etime(time(3,:), time(1,:));

% Velocity Vector estimation using Taylor Series expansion
v2 = -t32*(1/(t21*t31)+mu/(12*norm(r(:,1))^3))*r(:,1) + ...
    (t32 - t21)*(1/(t21*t32)+mu/(12*norm(r(:,2))^3))*r(:,2)+...
    t21*(1/(t32*t31)+mu/(12*norm(r(:,3))^3))*r(:,3);    % (km/s)
end