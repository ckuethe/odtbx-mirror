function [sat_Pos, sat_Vel] = Laplace(meas, time, rsite)
% LAPLACE Laplace's method of angles-only initial orbit determination
%
%    [sat_Pos, sat_Vel] = Laplace(meas, time, rsite) Utilizes Laplace's 
%    method of initial orbit determination to provide an intertial position
%    and velocity vector for an orbiting body. 
% 
% INPUTS:
% meas - Angle observations in a nx2 matrix [RA declination] (degrees)
% time - nx6 matrix of n datevectors for n observation times (datevector)
% rsite - 3xn matrix of site position vectors in the intertial frame (km)
%
% WARNING: Input units must be degrees for meas, datevector format
% for time, and km for rsite.
%
% OUTPUTS:
% sat_Pos - Orbiting object's position vector for second observation time (km)
% sat_Vel - Orbiting object's velocity vector for the second observation
%           time (km/s)
%
% CALLS:
% 1) LOS_Vectors
% 2) GeneralLagrange

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

disp('%%%%%%%%%% LAPLACE METHOD %%%%%%%%%%');

% Angular velocity of site location relative to IJK frame
w = [0 0 0.000072921158553]';    % (rad/s)

% Gravitational Parameter
global muglobal
mu = muglobal;    % (km^3/s^2)

% Obtains the LOS unit vectors from the given RA and dec observations
LOS = LOS_Vectors(meas);

%% Lagrange Interpolation
% Lagrange interpolation formula to estimate unit position vector at any
% time. Accuracy of estimation increases with a larger number of
% observations.
try
    % will fail for more than 3 obs:
    L_gen = SimpleLagrange(time, LOS);    % (1, s^-1, s^-2)
catch
    % requires Symbolic Toolbox:
    L_gen = GeneralLagrange(time, LOS);    % (1, s^-1, s^-2)
end

% Calculation of average site velocity and accelaration vectors
if length(rsite) > 1
     r2site = GeneralLagrange(time, rsite);
     r2site_dot = r2site(:,2);
     r2site_2dot = r2site(:,3);
     r2site = r2site(:,1);
else
    % If the observations occur from one site in intertial frame
    r2site = rsite;
    r2site_dot = cross(w, rsite);
    r2site_2dot = cross(w, r2site_dot);
end

%% Position Magnitude Estimation 
% Setup of matrix determinants that will be needed for iterative
% calculation of radius from middle site to object
D = 2*det(L_gen);    % (s^-3)
D1 = det([L_gen(:,1) L_gen(:,2) r2site_2dot]);    % (km/s^3)
D2 = det([L_gen(:,1) L_gen(:,2) r2site]);    % (km/s)

C = dot(LOS(:,2), rsite(:,2));    % (km)

% Eighth Order polynomial coefficients for finding the middle radius.
% Polynomial is of the form r^8 + c1*r^6 + c2*r^3 + c3 = 0
c1 = 4*C*D1/D - (4*D1^2)/D^2 - dot(r2site,r2site);    % (km^2)
c2 = mu*(4*C*D2/D - 8*D1*D2/D2);    % (km^5)
c3 = -4*mu^2*D2^2/D^2;    %(km^8)

% Now the roots of the eighth order polynomial must be solved for. 
r2 = roots([1 0 c1 0 0 c2 0 0 c3]);

% Now extract only the real roots
r2 = r2(r2 == real(r2));

% The positive of the two real roots represents the first iteration
% solution. Now the guess for the f and g series coefficient u can be
% updated. 
if r2(1) > 0
   r2 = r2(1);
elseif r2(2) > 0
    r2 = r2(2);
else
    error('No positive root found. No solution\n');
end

% Determine the magnitude of the slant-range vector
rho = -2*D1/D - 2*mu*D2/(D*r2^3);

sat_Pos = rho*LOS(:,2) + r2site;    % (km)

%% Velocity Vector Calculation
% Repeat of the process to calculate the velocity vector associated with 
% the middle observation

D3 = det([L_gen(:,1), r2site_2dot, L_gen(:,3)]);
D4 = det([L_gen(:,1), rsite(:,2), L_gen(:,3)]);

% Range rate
rho_dot = -D3 / D - mu*D4/(D*r2^3);

% Velocity
sat_Vel = rho_dot*LOS(:,2) + rho*L_gen(:,2) + r2site_dot;
end
