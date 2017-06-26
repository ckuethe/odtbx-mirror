function [sat_Pos, sat_Vel, Qo, Iterations] = Double_r(meas, time, rsite)
% DOUBLE_R Double-r method of angles-only initial orbit determination
%
%    [Sat_Pos, sat_Vel, Q, Iterations] = Double_r(meas, time, rsite) Uses
%    the Double-r initial orbit determination technique to provide an
%    intertial position and velocity vector for an orbiting body.
%
% INPUTS:
% meas - Angle observations in a nx2 matrix [RA Declination] (degrees)
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
% Q - Accuracy estimation of solution
% Iterations - Number of iterations to reach convergence
%
% CALLS:
% 1) LOS_Vectors
% 2) Double_r_OrbitElements
% 3) Double_r_Accuracy

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

%% Preliminary Orbit Search
disp('%%%%%%%%%% DOUBLE-R METHOD %%%%%%%%%%');

global muglobal
global Type
global IODSettings
% Constants
mu = muglobal;    % (km^3/s^2)

Qmin = 0.20;    % Upper bound on accuracy tolerance

% Default initial range estimates
[r1o, r2o] = Double_r_RangeEstimates(meas, time, rsite);

% Can adjust accracy tolerances for different situations as needed
if strcmpi(Type, 'Earth')
    Qmin = 0.20;
elseif strcmpi(Type, 'Sun')
    Qmin = 0.20;
else
    Qmin = 0.20;
end
r1initial = r1o;
r2initial = r2o;

% Magnitude of site position vectors
rs1mag = norm(rsite(:,1));
rs2mag = norm(rsite(:,2));
rs3mag = norm(rsite(:,3));

% Line of site vectors for each observation
LOS = LOS_Vectors(meas);

% Difference in observation times
tau1 = etime(time(1,:),time(2,:));
tau3 = etime(time(3,:),time(2,:));

% Set the inital accuracy with the r1o and r2o default values, Qo.
% Iterative guesses of r1o and r2o will produce new accuracy values which
% will be compared to this value.
[KOE, M32, M21, n, R2, R3, E32] = Double_r_OrbitElements(r1o, r2o, rsite, LOS, time);
F1 = tau1 + M21/n;
F2 = tau3 - M32/n;

Qo = Double_r_Accuracy(M32, M21, n, time);    % Initial Accuracy
fprintf(1,'Initial Accuracy:%f\n',Qo);

sat_Pos = R2;

% Iteration counter
Iterations = 1;

% Number of time Qmin has been increased
Resets = 0;

% Setting the max number of allowable iterations
if ~isfield(IODSettings, 'Iterations')
    % Default value
    maxIter = 500;
else
    % User-defined iteration max
    maxIter = IODSettings.Iterations;
end

%% Loop here until accuracy of estimate is less than tolerance Qmin
while (Qo > Qmin)
    fprintf(1,'Iteration: %f\n',Iterations);
    % Set the initial change in magnitude of the estimate
     delta_r1 = 0.00001*r1o;    % Default step size of 10^-5 km
     delta_r2 = 0.00001*r2o;
    
    % Repeat Calculations for F1
    [KOE, M32, M21, n, R2temp, R3temp, E32temp] = Double_r_OrbitElements(r1o+delta_r1, r2o, rsite, LOS, time);
    F1new = tau1 + M21/n;
    F2new = tau3 - M32/n;
    dF1dr1 = (F1new - F1)/delta_r1;
    dF2dr1 = (F2new - F2)/delta_r1;
    
    % Repeat Calculations for F2
    [KOE, M32, M21, n, R2temp, R3temp, E32temp] = Double_r_OrbitElements(r1o, r2o+delta_r2, rsite, LOS, time);
    F1new = tau1 + M21/n;
    F2new = tau3 - M32/n;
    dF1dr2 = (F1new - F1)/delta_r2;
    dF2dr2 = (F2new - F2)/delta_r2;
    
    delta = dF1dr1*dF2dr2 - dF2dr1*dF1dr2;
    delta1 = dF2dr2*F1 - dF1dr2*F2;
    delta2 = dF1dr1*F2 - dF2dr1*F1;
    
    delta_r1 = -delta1/delta;
    delta_r2 = -delta2/delta;
    
    % Update range estimates
    r1o = r1o + delta_r1;
    r2o = r2o + delta_r2;

    % Check to see if the accuracy has improved
    [KOE, M32, M21, n, R2temp, R3temp, E32temp, v2] = Double_r_OrbitElements(r1o, r2o, rsite, LOS, time);
    Q = Double_r_Accuracy(M32, M21, n, time);
    if Q < Qo && imag(Q) == 0
       Qo = Q; 
       sat_Pos = R2temp;
       R3 = R3temp;
       E32 = E32temp;
       fprintf(1,'tau1: %f\n',M21/n);
       fprintf(1,'Qo: %f\n',Qo);
    end

    % Update the F-values for the next pass
    F1 = tau1 + M21/n;
    F2 = tau3 - M32/n;

    Iterations = Iterations + 1;
    

    % Divergence counter value can be increased, but it is unlikely to
    % converge after 500 iterations. The error tolerance 
    if Iterations >= maxIter
        beep;
        errordlg(['Solution diverged. Reached maximum iterations.'],'Computation Error','modal');
        fprintf(1,'WARNING: Solution diverged. Stopped at %d iterations.\n', maxIter);
        sat_Pos = NaN;
        sat_Vel = NaN;
        Qo = NaN;
        Iterations = maxIter;
        return;
    end
end

f = 1.0 - (KOE.sma/r2o)*(1.0 - cos(E32));
g = tau3 - sqrt(KOE.sma^3/mu)*(E32 - sin(E32));

% Calculation of final position and velocity vectors
if strcmpi(Type,'Sun')
    % V2 value sent up from Gibbs or Herrick-Gibbs method
    sat_Vel = v2;
elseif strcmpi(Type,'Earth')
    sat_Vel = (R3 - f*sat_Pos) / g;
else
    sat_Vel = v2;
end
end

   