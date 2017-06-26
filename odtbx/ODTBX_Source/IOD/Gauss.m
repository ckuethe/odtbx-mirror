function [sat_Pos, sat_Vel, Iterations, r1, r3] = Gauss(meas, time, rsite)
% GAUSS Gauss method of angles-only initial orbit determination
%
%    [sat_Pos, sat_Vel, Iterations] = Gauss(meas, time, rsite) Uses Gauss' 
%    initial orbit determination technique in order to provide an intertial
%    position and velocity vector for an orbiting body.
%
% INPUTS
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
% Iterations -  Number of iterations to reach convergence
% r1 - Orbiting object's position vector for first observation time. Needed
%      for the initial range estimates of Double-r method. (km)
% r3 - Orbiting object's position vector for third observation time. Needed
%      for the initial range estimates of Double-r method. (km)
%
% CALLS:
% 1) LOS_Vectors
% 2) Gibbs
% 3) Herrick-Gibbs
% 4) fgcalc

% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/04/2015                Original

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

disp('%%%%%%%%%% GAUSS METHOD %%%%%%%%%%');
% Gravitational Parameter
global muglobal;
mu = muglobal;   % (km^3/s^2)

Ltemp = [];
for i = 1:3
    Ltemp = [Ltemp; [-cosd(meas(i,1))*cosd(meas(i,2)); sind(meas(i,1))*cosd(meas(i,2)); sind(meas(i,2))]];
    
end

LOS = LOS_Vectors(meas);

% Change in observation times
t1 = etime(time(1,:), time(2,:));
t3 = etime(time(3,:), time(2,:));

% Rewriting coefficients
a1 = t3 / (t3 -t1);
a1u = t3*((t3-t1)^2 - t3^2) / (6*(t3 - t1));
a3 = -t1/ (t3-t1);
a3u = -t1*((t3-t1)^2 - t1^2) / (6*(t3 - t1));

M = inv(LOS) * rsite;

d1 = M(2,1)*a1 - M(2,2) + M(2,3)*a3;
d2 = M(2,1)*a1u + M(2,3)*a3u;

C = dot(LOS(:,2), rsite(:,2));

% The coefficients of the eighth order polynomial that must be solved in
% order to find the middle radius magnitude can be condensed for 
% simplicity. Given the following form: r^8 + q1*r^6 + q2*r^3 + q3 = 0;
q1 = -(d1^2 + 2*C*d1 + dot(rsite(:,2),rsite(:,2)));
q2 = -2*mu*(C*d2 + d1*d2);
q3 = - mu^2*d2^2;

% The roots of the eighth order polynomial are solved using MATLAB's 
% roots() command 
r2 = roots([1 0 q1 0 0 q2 0 0 q3]);

% Extract only the real roots
r2 = r2(r2 == real(r2));

% The positive of the two real roots represents the first iteration
% solution. Now the guess for the f and g series coefficient u can be
% updated. 
if r2(1) > 0
   r2 = r2(1);
elseif r2 (2) > 0
    r2 = r2(2);
else
    disp(muglobal)
    error('No positive, real roots found. No Solution.\n');
    
end

% Update guess for the f and g series ceofficient
u = mu / r2^3;

% Calculation of the original coefficients which defined the linear
% independence of each position vector
c1 = a1 + a1u*u;
c2 = -1;
c3 = a3 + a3u*u;
c = [c1; c2; c3];

% Calculate the initial guess of the slant-ranges, rho
cp = -1*(M*c);
rhonot = [cp(1)/c1; cp(2)/c2; cp(3)/c3];    % (ER)
rho = rhonot;
rho_next = [0; 0; 0];

% Initial position vectors estimate
r = [rhonot(1)*LOS(:,1) + rsite(:,1) rhonot(2)*LOS(:,2) + rsite(:,2) rhonot(3)*LOS(:,3) + rsite(:,3)];    % (ER)

% Iteration count to keep track of how many iterations until convergence
i = 0;

% Set the error tolerance for the iterative calculations of rho values
Error = rhonot(1)*1e-6;

% Counter for number of times the error tolerance is increased
error_count = 0;

%% Iterative Convergance of Solution
% Iterates until slant-range estimate stops changing significantly
while (abs(rho(1)-rho_next(1)) > Error && abs(rho(2)-rho_next(2)) > Error && abs(rho(3)-rho_next(3)) > Error)
%     disp(abs(rho(1)-rho_next(1)))
    % Check the angular separation to determine whether to use Gibbs or
    % Herrick-Gibbs method. Uses Gibbs for angles greater than 3 degrees
    % and Herrick-Gibbs for angles less than 3 degrees
    alpha_12 = acosd(dot(r(:,1),r(:,2))/(norm(r(:,1))*norm(r(:,2))));
    alpha_23 = acosd(dot(r(:,2),r(:,3))/(norm(r(:,2))*norm(r(:,3))));

%     Use Gibbs method if separation angles are large enough
    if (abs(alpha_12) > 3 && abs(alpha_23) > 3)
        % Gibbs returns the velocity, eccentricity, and semiparameter.
        [v2, e, p] = Gibbs(r);
        
    else
        % Use Herrick-Gibbs if separation angles are small
        v2 = Herrick_Gibbs(r,time);
        KOE = kepel(r(:,2),v2,mu);
        p = KOE.sma*(1-KOE.ecc^2);
    end

%     % Using Lambert's Problem with Minimum Energy Solution
%     v2 = Lambert_Min_Energy(r(:,2),r(:,3));
%     disp('Lambert v2')
%     disp(v2);
%     KOE = kepel(r(:,2),v2,mu);
%     p = KOE.sma*(1-KOE.ecc^2);
    
    %% Calculation of eccentricity and semiparameter using ELORB techniques

    % Angular momentum vector
    h = cross(r(:,2),v2);

    % Node vector
    n = cross([0;0;1],h);

    % Eccentricity vector
    e = ((dot(v2,v2)-mu/norm(r(:,2)))*r(:,2) - dot(r(:,2),v2)*v2)/mu;

    % Specific mechanical energy (sme)
    sme = dot(v2,v2)/2 - mu/norm(r(:,2));

    if e == 1.0
        % Semiparameter
        p = dot(h,h)/mu;
    else
        % Semimajor Axis
        a = -mu/(2*sme);
        p = a*(1-dot(e,e));
    end
   
    % Calculation of f and g coefficients
    f1 = 1 - (norm(r(:,1)) / p)*(1-cosd(-alpha_12));
    f3 = 1 - (norm(r(:,3)) / p)*(1-cosd(alpha_23));
    g1 = (norm(r(:,1))*norm(r(:,2))*sind(-alpha_12))/sqrt(mu*p);
    g3 = (norm(r(:,3))*norm(r(:,2))*sind(alpha_23))/sqrt(mu*p);
    c1 = g3 / (f1*g3 - f3*g1);
    c3 = -g1 / (f1*g3 - f3*g1);
    
    % Recalculation of slant range distances
    c = [c1; c2; c3];
    cp = M*-1*c;
    rho_temp = [cp(1)/c1; cp(2)/c2; cp(3)/c3];
   
    if i == 1
        rho = rho;
        rho_next = rho_temp;
    else
        rho = rho_next;
        rho_next = rho_temp;
    end

    i = i + 1;
    fprintf(1,'Iteration: %f\n',i);
    r = [rho_next(1)*LOS(:,1) + rsite(:,1) rho_next(2)*LOS(:,2) + rsite(:,2) rho_next(3)*LOS(:,3) + rsite(:,3)];
    if i >= 500
       Error = Error * 10;
       error_count = error_count + 1;
       i = 0;
       r = [rhonot(1)*LOS(:,1) + rsite(:,1) rhonot(2)*LOS(:,2) + rsite(:,2) rhonot(3)*LOS(:,3) + rsite(:,3)]; 
       fprintf(1,'Increasing error tolerance by factor of 10.\n');
    end
    if error_count == 6
       beep;
       errordlg(['Solution diverged. Stopped after 3000 iterations'],'Computation Error','modal');
       fprintf(1,'WARNING: Solution diverged. Stopped at 3000 iterations.\n');
       sat_Pos = NaN;
       sat_Vel = NaN;
       return;
    end
end

% Updated position vectors estimate
sat_Pos = r(:,2);    % (km)
sat_Vel = v2;    % (km/s)
r1 = r(:,1);
r3 = r(:,3);
Iterations = i;
end