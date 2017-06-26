function [v2, e, P] = Gibbs(r)
% GIBBS Gibbs method of velocity determination given three position vectors
%
%    [v2, e, P] = Gibbs(r) Uses the Gibbs method to return an estimation 
%    of the velocity vector for the middle observation time out of three
%    observations. The function also estimates the eccentricity, e, and
%    semilatus rectum, P.
% 
% INPUTS:
% r - 3x3 matrix containing orbiting object's position vectors at three
%     sequential observation times. (km)
%
% OUTPUTS:
% v2 - Estimation of object's velocity vector at second time (km/TU)
% e - Eccentricity estimate of the orbit
% P - Semiparameter of the orbit (km)
%
% CALLS:
% None
%
% WARNING: To use the Gibbs method, check that the position vectors are
% coplanar to a certain tolerance (<3deg) and that the separation between
% the vectors is large (>1deg). Failure to do this will cause erroneous
% estimations. This function assumes that these checks are done before
% passing arguments through, but will check again and provide error
% messages to indicate potential errors.

% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/05/2015                Original

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

alpha_12 = acosd(dot(r(:,1),r(:,2))/(norm(r(:,1))*norm(r(:,2))));
alpha_23 = acosd(dot(r(:,2),r(:,3))/(norm(r(:,2))*norm(r(:,3))));

if (abs(alpha_12) > 1 && abs(alpha_23) > 1)
else
    fprintf(1, ['\nWARNING: Angle between position vectors is smaller than',...
        ' 1 degree. This will result in erroneous results when using Gibbs',...
        ' Method!\n']);
    fprintf(1,'Angle between first and second vector is: %f deg\n',alpha_12);
    fprintf(1,'Angle between second and third vector is: %f deg\n',alpha_23);
end

N = norm(r1)*Z23 + norm(r2)*Z31 + norm(r3)*Z12;
D = Z12 + Z23 + Z31;
S = (norm(r2) - norm(r3))*r1 + (norm(r3) - norm(r1))*r2 + (norm(r1) - norm(r2))*r3;
B = cross(D, r2);

Lg = sqrt(mu/(norm(N)*norm(D)));

% Direction of eccentricity (e_hat) used to find the true anomaly
W_hat = N / norm(N);
Q_hat = S / norm(S);
e_hat = cross(Q_hat, W_hat);

v2 = (Lg/norm(r2))*B + Lg*S;    % (km/s)
e = norm(S)/norm(D);
P = norm(N)*e / norm(S);    % (km)
end