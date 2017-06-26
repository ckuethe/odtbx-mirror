function [r1o, r2o] = Double_r_RangeEstimates(meas, time, rsite)
% DOUBLE_R_RANGEESTIMATES Estimates starting range values for Double-r method
%
%    [r1o, r2o] = Double_r_range_estimates(meas, time, rsite) Finds a pair
%    of range values, r1o and r2o, that generate a specified accuracy for a
%    given orbit. These initial range values can be used as the starting
%    values for further differential correction techniques. This estimator
%    uses the cone masking technique outlined in GTDS.
%
% INPUTS:
% meas - Angle observations in a nx2 matrix [RA Declination] (degrees)
% time - nx6 matrix of n datevectors for n observation times (datevector)
% rsite - 3xn matrix of site position vectors in the intertial frame (km)
%
% OUTPUTS:
% r1o - Range estimate at first observation time (km)
% r2o - Range estimate at second observation time (km)
% 
% CALLS:
% 1) Double_r_OrbitElements
% 2) Double_r_Accuracy
%
% WARNING: This accurcay calculation assumes that consecutive measurements
% are taken on the same pass

% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/19/2015                Original
%   Ryan Willmot    07/14/2015                Fix issues with heliocentric
%   Ryan Willmot    07/17/2015                GUI user-defined tolerances

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

% ALGORITHMIC ESTIMATION OF UPPER AND LOWER BOUNDS ON SPACECRAFT HEIGHTS
global Type
global IODSettings
global radius
% Magnitude of site position vectors
rs1mag = norm(rsite(:,1));
rs2mag = norm(rsite(:,2));
rs3mag = norm(rsite(:,3));

% Initial range estimate
if strcmpi(Type, 'Sun')
    % Initial range estimate is radius of Sun, then iteration moves guess
    % outward toward actual solution
    if ~isfield(IODSettings, 'r2o')
        % Default values
        r1o = 0.75*149597871;    % Starting value of 0.75*1AU
        r2o = r1o;
    else
        % User-defined initial estimates
        r1o = IODSettings.r2o*149597871;
        r2o = r1o;
    end
    
    % Initial height estimate at time t2
    if ~isfield(IODSettings, 'h2')
        h2 = radius;    % (km)
    else
        h2 = IODSettings.h2*149597871;
    end

    if ~isfield(IODSettings, 'Qmin')
        % Default Value
        Qmin = 3000;
    else
        Qmin = IODSettings.Qmin;
    end
    
    K = 5;
else
    if ~isfield(IODSettings, 'r2o')
        % Default values
        r1o = radius;
        r2o = r1o; 
    else
        
        % User-defined initial estimates
        r1o = IODSettings.r2o;
        r2o = r1o;
    end
    
    % Initial height estimate at time t2
    if ~isfield(IODSettings, 'h2')
        h2 = radius;    % (km)
    else
        h2 = IODSettings.h2;
    end
    
    % Default value
    K = 1.25;
    if ~isfield(IODSettings, 'Qmin')
        % Default value
        Qmin = 150;
    else
        Qmin = IODSettings.Qmin;
    end
end

% Line of site vectors for each observation
LOS = LOS_Vectors(meas);

% Difference in observation times
tau1 = etime(time(1,:),time(2,:));
tau3 = etime(time(3,:),time(2,:));

% Need to set the inital accuracy with the r1o and r2o default values.
% Iterative guesses of r1o and r2o will produce new accuracy values which
% will be compared to this value.
[KOE, M32, M21, n, R2, R3, E32] = Double_r_OrbitElements(r1o, r2o, rsite, LOS, time);
F1 = tau1 + M21/n;
F2 = tau3 - M32/n;

Qo = Double_r_Accuracy(M32, M21, n, time);
fprintf(1,'Q1:%f\n',Qo);

% Level Set and iteration counter
L = 0;

% Set max levels
if ~isfield(IODSettings, 'Levels')
    Lmax = 20;
else
    Lmax = IODSettings.Levels;
end

% Current number of Levels checked
Lcount = 0;

% h2 reduction count
h2reduc = 0; 

% h2 reduction factor
if ~isfield(IODSettings, 'h2reducFactor')
    h2reducFactor = 2;
else
    h2reducFactor = IODSettings.h2reducFactor;
end

%% Forming Initial Guess
while Qo > Qmin
    
    % Formation of integer exponents using diamond pattern in Figure 9-3 
    % of GTDS. This is based on a cone-masking technique.
    if L == 0
        E1 = 0;
        E2 = E1;
    else
        E1 = zeros(1, L*4);
        E2 = E1; 
        E1(1) = L;
        for j = 2:L*2
            E1(j) = E1(j-1) - 1;
            E2(j) = L - abs(E1(j));
        end 
        k = 1;
        for j = L*2+1:L*4
             E1(j) = -E1(k);
             E2(j) = -E2(k);
             k = k + 1;
        end
    end 
    J = max(5/6, 1/K);
    
    % Formation of R2mag guesses
    R2magnewpos = r2o+h2*K.^E2;    % Length = 4*L
    R1magnewpos = r1o+h2*K.^E1;    % Length = 4*L
    R2magnewneg = r2o-h2*K.^E2;
    R1magnewneg = r1o-h2*K.^E1;
   
    for i =1:length(E1)
        for j = 1:length(E2)
            [KOE, M32, M21, n, R2, R3, E32] = Double_r_OrbitElements(R1magnewpos(i), R2magnewpos(j), rsite, LOS, time);
            Q = Double_r_Accuracy(M32, M21, n, time);
            if Q < Qo && imag(Q) == 0
               
                Qtemp = Qo;
                Qo = Q;
                
                % Shift search pattern to newly found values
                L = 0;
                r1o = R1magnewpos(i);
                r2o = R2magnewpos(j);
                fprintf(1,'Qo:%f\n',Qo);
                if Qo < Qtemp / 2
                    % Reducing step size
                    h2 = h2 / (1 / h2reducFactor^-1);
                    h2reduc = h2reduc + 1;
                end
                if abs(Qtemp - Q) < 5 && (strcmpi(Type, 'Earth')||strcmpi(Type, 'Other'))
                    h2 = h2*2;
                elseif abs(Qtemp - Q) < 50 && strcmpi(Type, 'Sun')
                    h2 = h2*2;
                end
                if h2reduc > 5 && abs(Qtemp - Q) < 1e-5
                    h2 = h2 * (1 / h2reducFactor^-1)^h2reduc;
                    h2reduc = 0;
                end
            end
            
%             % Negative checks only needed if initial range guess could be
%             % larger than true range values
%             if R1magnewneg(i) > 0 && R2magnewneg(j) > 0
%                 [KOE, M32, M21, n, R2, R3, E32] = Double_r_OrbitElements(R1magnewneg(i), R2magnewneg(j), rsite, LOS, time);
%                 Q = Double_r_Accuracy(M32, M21, n, time);
%                 if Q < Qo && imag(Q) == 0
%                     Qtemp = Qo;
%                     Qo = Q;
%                     L = 0;
%                     r1o = R1magnewneg(i);
%                     r2o = R2magnewneg(j);
%                     fprintf(1,'Qo:%f\n',Qo);
%                     if Qo < Qtemp / 2
%                         % Reducing step size
%                         h2 = h2 / 2;
%                         h2reduc = h2reduc + 1;
%                     end
%                     if abs(Qtemp - Q) < 5
%                         h2 = h2*2;
%                     end
%                     if h2reduc > 5 && abs(Qtemp - Q) < 1e-5
%                         h2 = h2 * 2;
%                         h2reduc = 0;
%                     end
%                 end
%             end
%             if R1magnewneg(i) > 0
%                 [KOE, M32, M21, n, R2, R3, E32] = Double_r_OrbitElements(R1magnewneg(i), R2magnewpos(j), rsite, LOS, time);
%                 Q = Double_r_Accuracy(M32, M21, n, time);
%                 if Q < Qo && imag(Q) == 0
%                     Qtemp = Qo;
%                     Qo = Q;
%                     L = 0;
%                     r1o = R1magnewneg(i);
%                     r2o = R2magnewpos(j);
%                     fprintf(1,'Qo:%f\n',Qo);
%                     if Qo < Qtemp / 2
%                         % Reducing step size
%                         h2 = h2 / 2;
%                         h2reduc = h2reduc + 1;
%                     end
%                     if abs(Qtemp - Q) < 5
%                         h2 = h2*2;
%                     end
%                     if h2reduc > 5 && abs(Qtemp - Q) < 1e-5
%                         h2 = h2 * 2;
%                         h2reduc = 0;
%                     end
%                 end
%             end
%             if R2magnewneg(i) < 0
%                 [KOE, M32, M21, n, R2, R3, E32] = Double_r_OrbitElements(R1magnewpos(i), R2magnewneg(j), rsite, LOS, time);
%                 Q = Double_r_Accuracy(M32, M21, n, time);
%                 if Q < Qo && imag(Q) == 0
%                     Qtemp = Qo;
%                     Qo = Q;
%                     L = 0;
%                     r1o = R1magnewpos(i);
%                     r2o = R2magnewneg(j);
%                     fprintf(1,'Qo:%f\n',Qo);
%                     if Qo < Qtemp / 2
%                         % Reducing step size
%                         h2 = h2 / 2;
%                         h2reduc = h2reduc + 1;
%                     end
%                     if abs(Qtemp - Q) < 5
%                         h2 = h2*2;
%                     end
%                     if h2reduc > 5 && abs(Qtemp - Q) < 1e-5
%                         h2 = h2 * 2;
%                         h2reduc = 0;
%                     end
%                 end
%             end

        end
    end
    L = L + 1;
    Lcount = Lcount + 1;
    
    % Case if no feasible range estimates are found using cone-masking
    % technique
    if L > Lmax
        fprintf(1,'ERROR: No solution found. Either passed it or did not reach it.\n');
        fprintf(1,'Using range estimates from Gauss method.\n');
        [r2o, sat_Vel, Iterations, r1o] = Gauss(meas, time, rsite);
        r1o = norm(r1o);
        r2o = norm(r2o);
        return;  
    end
end

end