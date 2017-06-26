function [KOE, M32, M21, n, R2, R3, E32, v2] = Double_r_OrbitElements(R1mag, R2mag, rsite, LOS, time)
% DOUBLE_R_ORBITELEMENTS Orbital elements calculation for Double-r method
%
%    [KOE, M32, M21, n, R2, R3, E32] = Double_r_OrbitElementS(R1mag, R2mag,
%    rsite, LOS, time) Computes the orbital elements for the Double-r 
%    initial orbit determination method.
%
%    This method is based off of Goddard Trajectory Determination System 
%    (GTDS) Mathematical Theory, Revision 1 
%    (1989, Long, Cappellari, Velez, Fuchs).
%
% INPUTS:
% R1mag - Magnitude of orbiting object's first position vector (km)
% R2mag - Magnitude of orbiting object's second position vector (km)
% rsite - 3x3 matrix of position vectors for observation sites (km)
% LOS - 3x3 matrix of line of site unit vectors for observations
% time - nx6 matrix of n datevectors for n observation times (datevector)
%
% OUTPUTS:
% KOE - Structure containing Keplerian Orbital Elements
%   KOE.sma - Semi-major axis (km)
%   KOE.incl - Inclination (rad)
%   KOE.ecc - Eccentricity
%   KOE.tran - True anomaly (rad)
% M32 - Change in mean anomaly from observation 2 to 3 (rad)
% M21 - Change in mean anomaly from observation 1 to 2 (rad)
% n - True anomaly at observation 2 (rad)
% R1 - Inertial position vector at observation 1 time (km)
% R2 - Interial position vector at observation 2 time (km)
% E32 - Change in eccentric anomaly between 3 and 2 (rad)
% v2 - Velocity at second observation time (km/s)
%
% CALLS:
% 1) kepel (ODTBX)
% 2) Gibbs
% 3) Herrick_Gibbs

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

% Constants
global muglobal
global radius
global Type
Re = radius;
mu = muglobal;

% Time differences between observations
t1 = etime(time(1,:),time(2,:));
t3 = etime(time(3,:),time(2,:));

% Magnitude of site position vectors
rs1mag = norm(rsite(:,1));
rs2mag = norm(rsite(:,2));
rs3mag = norm(rsite(:,3));

C1 = 2*dot(LOS(:,1),rsite(:,1));
C2 = 2*dot(LOS(:,2),rsite(:,2));
rho1 = 0.5*(-C1+sqrt(C1^2 - 4*(rs1mag^2 - R1mag^2)));
rho2 = 0.5*(-C2+sqrt(C2^2 - 4*(rs2mag^2 - R2mag^2)));

% Intertial Position vectors at first and second times
R1 = rsite(:,1) + rho1*LOS(:,1);
R2 = rsite(:,2) + rho2*LOS(:,2);
r1mag = norm(R1);
r2mag = norm(R2);
k = cross(R1,R2)/(r1mag*r2mag);

rho3 = abs(dot(rsite(:,3),k) / dot(LOS(:,3),k));

R3 = rsite(:,3) + rho3*LOS(:,3);
r3mag = norm(R3);
r = [R1 R2 R3];

% Vector parallel to the angular momentum unit vector
W = cross(R1,R2)/(r1mag*r2mag);
Wz = W(3);

% Checking direction of orbit (direct or retrograde)  
Tmin = 2*pi*sqrt(Re^3/mu);
if Tmin < t3 && t3 < -t1
    s = sign(cross(R2,R3));
elseif -t1 < Tmin && -t1<t3
    s = sign(Wz);
else
    s = 1;    % Default assumption that orbit is direct
end

% Calculating the difference in the true anomalies
cosnu21 = dot(R2,R1)/(r2mag*r1mag);
sinnu21 = s*sqrt(1-cosnu21^2);
nu21 = acosd(cosnu21);
cosnu31 = dot(R3,R1)/(r3mag*r1mag);
sinnu31 = s*sqrt(1-cosnu31^2);
nu31 = atand(sinnu31/cosnu31);
if cosnu31 < 0 && sinnu31 > 0
    nu31 = nu31 + 180;
elseif cosnu31 < 0 && sinnu31 < 0
    nu31 = nu31 + 180;
end
cosnu32 = dot(R3,R2)/(r3mag*r2mag);
sinnu32 = s*sqrt(1-cosnu32^2);
nu32 = atand(sinnu32/cosnu32);
if cosnu32 < 0 && sinnu32 > 0
    nu32 = nu32 + 180;
elseif cosnu32 < 0 && sinnu32 < 0
    nu32 = nu32 + 180;
end

% Calculation of orbit inclination, i
i = acos(Wz/sinnu21);    % (deg)
if strcmpi(Type,'Earth')
    
    % Calculation of semilatus rectum, p
    Ct1 = r1mag*sinnu31/(r2mag*sinnu32);
    Ct3 = r1mag*sinnu21/(r3mag*sinnu32);
    c1 = r2mag*sinnu32/(r1mag*sinnu31);
    c3 = r2mag*sinnu21/(r3mag*sinnu31);
    
    % Temp value for v2 to be send back up and recalculated in Double_r
    v2 = 0;
    if abs(nu31) >= abs(nu32)
        p = (c1*r1mag + c3*r3mag - r2mag) / (c1 + c3 -1);    % (km)

    else
        p = (r1mag + Ct3*r3mag - Ct1*r2mag) / (1 + Ct3 - Ct1);   % (km)
    end
else
    % Gibbs/Herrick Gibbs and Kepel to get p
    
    % Check the angular separation to determine whether to use Gibbs or
    % Herrick-Gibbs method. Uses Gibbs for angles greater than 3 degrees
    % and Herrick-Gibbs for angles less than 3 degrees
    alpha_12 = acosd(dot(r(:,1),r(:,2))/(norm(r(:,1))*norm(r(:,2))));
    alpha_23 = acosd(dot(r(:,2),r(:,3))/(norm(r(:,2))*norm(r(:,3))));

    % Use Gibbs method if separation angles are large enough

    if (abs(alpha_12) > 3 && abs(alpha_23) > 3)
        % Gibbs returns the velocity, eccentricity, and semiparameter.
        [v2, e, p] = Gibbs(r);
    else
        % Use Herrick-Gibbs method if separation angles are small
        v2 = Herrick_Gibbs(r,time);
        KOE = kepel(r(:,2),v2,mu);
        p = KOE.sma*(1-KOE.ecc^2);
    end
end

% Determining the eccentricity of the orbit, e
ecosnu1 = p/r1mag - 1;
ecosnu2 = p/r2mag - 1;
ecosnu3 = p/r3mag - 1;

if abs(sinnu21) > abs(sinnu32)
    esinnu2 = (-ecosnu2*cosnu21 + ecosnu1)/sinnu21;
else
    esinnu2 = (ecosnu2*cosnu32 - ecosnu3)/sinnu32;
end

e = sqrt(ecosnu2^2 + esinnu2^2);

% True Anomaly
nu1 = acos(ecosnu1/e);
nu2 = acos(ecosnu2/e);
nu3 = acos(ecosnu3/e);

% Semimajor axis, a
a = p / (1 - e^2);    % (km)

%%%%%%%%%%%%%NEED TO ADD IN OTHER CASES (HYPERBOLIC,PARABOLIC)%%%%%%%%%%
% For Elliptical Orbit
if e < 1
    % Mean motion
    n = sqrt(mu/a^3);    % (rad/s)
    
    % Orbital period 
    T = 2*pi/(60*n);    % (min)
    
    % Perigee height
    hp = a*(1-e) - Re;
    
    % Apogee height
    ha = hp + 2*a*e;
else
    n = 1;
%     fprintf(1,'Non-elliptical orbit!\n');
end

Se = (r2mag/p)*sqrt(1-e^2)*esinnu2;
Ce = (r2mag/p)*(e^2 + ecosnu2);

% Calculation of differences in eccentric anomalies
sinE32 = (r3mag/sqrt(a*p))*sinnu32-(r3mag/p)*(1-cosnu32)*Se;
cosE32 = 1 - (r3mag*r2mag/(a*p))*(1-cosnu32);
sinE21 = (r1mag/sqrt(a*p))*sinnu21+(r1mag/p)*(1-cosnu21)*Se;
cosE21 = 1-(r2mag*r1mag/(a*p))*(1-cosnu21);

if strcmpi(Type,'Earth')
    E32 = atan(sinE32/cosE32);
    E21 = atan(sinE21/cosE21);
    
    % Calculation of differences in mean anomalies
    M32 = E32 + 2*Se*(sin(E32/2))^2 - Ce*sin(E32);
    M21 = E21 - 2*Se*(sin(E21/2))^2 - Ce*sin(E21);
else
    E32 = acos(cosE32);
    E21 = acos(cosE21);
end

 % Calculation of differences in mean anomalies
M32 = E32 + 2*Se*(sin(E32/2))^2 - Ce*sin(E32);
M21 = E21 - 2*Se*(sin(E21/2))^2 - Ce*sin(E21);

KOE.sma = a;
KOE.ecc = e;
KOE.incl = i;
KOE.tran = nu2;
end