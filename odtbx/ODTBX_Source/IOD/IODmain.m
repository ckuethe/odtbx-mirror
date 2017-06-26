function [r2, v2, P, Iterations] = IODmain(filename, calcMethod, muType, DoublerSettings, mu, R, cosMatrix)
% IODMAIN Main file for Initial Orbit Determination Application. 
% Interact with this file during integration with another program.
%
%    [r2, v2, P] = IODmain(filename, calcMethod, cosMatrix) Computes a
%    position and velocity vector for a given orbit with observational data
%    contained in the string filename. Note that this function only handles
%    angles-only IOD methods (Laplace, Gauss, Double-r). The function also
%    returns a covariance matrix, P.
%
% This function is another interface for the IOD application. Where the GUI
% cannot interact with files and programs outside of the application, this
% function can.
% 
% WARNING: This function only works for angles-only IOD processing. The
% first three inputs are mandatory, the last four are optional. Default
% values will be used if only three inputs are given.
%
% INPUTS:
% filename - Full filename of observational data file
%            See documentation for proper data file formatting.
% calcMethod - Calculation method for angles-only IOD. Options are 'Laplace',
%          'Gauss', and 'Double-r'.
% muType - Signifies heliocentric or geocentric orbit. Options are 'Sun',
%          'Earth', or 'Other' for user-defined gravitational body.
% DoublerSettings - Structure for Double-r tolerances and iteration
%                   settings
%            .r2o - Initial range estimate (km for geocentric, AU for
%                   heliocentric)
%             .h2 - Range step size (km for geocentric, AU for
%                   heliocentric)
%           .Qmin - Accuracy tolerance on initial range estimation
%  .h2reducFactor - Step size reduction factor
%         .Levels - Number of search levels for cone-maskign technique
%     .Iterations - Allowable iterations for Double-r method
% 
% mu - User-defined gravitational parameter. For future functionality of
%      orbits around bodies other than the Sun or Earth
% R - User-defined radius for gravitational body (km).
% cosMatrix - Direction cosine matrix for conversion from Earth-Centered
%             Earth-Fixed frame to Earth-Centered Intertial frame. If no
%             argument is given, LST approximations are used.
%
% OUTPUTS:
% r2 - Estimation of orbiting object's position vector at the second time
%     (km)
% v2 - Estimation of orbiting object's velocity vector at the second time
%     (km/s)
% P - Covariance matrix
% Iterations - Number of iterations required to reach convergence
% 
% CALLS:
% 1) Site_Position
% 2) LST
% 3) Laplace
% 4) Gauss
% 5) Double_r
% 6) Casotto
% 7) kepel    (ODTBX)
% 8) sigpt

% REVISION HISTORY:
%   Author          Date (MM/DD/YYYY)         Comment
%   Ryan Willmot    06/23/2015                Original

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

if nargin < 3
    beep;
    errordlg('Not enough input arguments.','Input Error', 'modal')
    r2 = NaN;
    v2 = NaN;
    P = NaN;
    Iterations = NaN;
    return;
end
if nargin >= 4
    global IODSettings
    IODSettings = DoublerSettings;
end
if nargin == 3
    clearvars -global
end

% Global variables used by many functions in IOD application
global muglobal;
global radius;
global Type;
Type = muType;
if strcmpi(muType, 'Sun')
    if nargin < 5
        muglobal = 132712400180.6;    % km^3/s^2
    else
        muglobal = mu;
    end
    radius = 695800;    % km
elseif strcmpi(muType, 'Earth')
    % Default Earth's Gravitational Parameter and radius
    if nargin < 5
        muglobal = 398600.4418;    % km^3/s^2
    else
        muglobal = mu;
    end
    radius = 6378.137;    % km
elseif strcmpi(muType, 'Other')
    if nargin < 6
        beep;
        errordlg(['Not enough input arguments. Must enter a gravitational'...
            ' parameter and radius value for muType "Other".'...
            ' If not using Double-r method, input an empty structure for DoublerSettings.'], 'Input Error', 'modal');
        r2 = NaN;
        v2 = NaN;
        P = NaN;
        Iterations = NaN;
        return;
    else
        muglobal = mu;
        radius = R;
    end
else
    beep;
    errordlg('Invalid muType entered. Options are "Sun","Earth" or "Other".', 'Input Error', 'modal')
    disp(muType)
    r2 = NaN;
    v2 = NaN;
    P = NaN;
    Iterations = NaN;
    return;
end

% Check to see if the input file exists. For now, it is assumed that the
% input file is in the same directory as IODmain.m
if ~exist(filename, 'file')
    fprintf(1,'WARNING: File does not exist (%s)\n',filename);
    beep;
    errordlg('WARNING: File does not exist.');
    r2 = NaN;
    v2 = NaN;
    P = NaN;
    return;
end

% Check to see if file is formattted correctly before performing any
% calculations
fid = fopen(filename);
nObs = 0;
while ~feof(fid)
    strTime = fscanf(fid, '%c', 23);
    dataTemp = fgets(fid);
    dataTemp = str2num(dataTemp);

    % After reading in the date information from the file, there should
    % be exactly 8 values left for the angles-only data files
    if length(dataTemp) ~= 8
        fclose(fid);
        beep;
        errordlg(['File not formatted properly. Please check the',...
            ' application documentation for data file formatting',...
            ' guidelines.'], 'File Error','modal');
        fprintf(1, ['WARNING: File not formatted properly. Please',...
            ' check the application documentation for data file',...
            ' formatting guidelines.\n']);
        r2 = NaN;
        v2 = NaN;
        P = NaN;
        return;
    end

    % Counting the number of lines/observations
    nObs = nObs + 1; 
end 
frewind(fid);

% Initialization of all variable vectors
time = [];
data = [];
RA = [];
Dec = [];
RAsigma = [];
Decsigma = [];
RADecCorr = [];
lat = [];
long = [];
alt = [];
Px = [];


%% Reading in all data from input file
while ~feof(fid)
    strTime = fscanf(fid, '%c', 23);
    time = [time; datevec(strTime, 'yyyy mm dd HH:MM:SS.FFF')];
    dataTemp = fgets(fid);
    dataTemp = str2num(dataTemp);
    RA = [RA; dataTemp(1)];
    RAsigma = [RAsigma; dataTemp(2)];
    Dec = [Dec; dataTemp(3)];
    Decsigma = [Decsigma; dataTemp(4)];
    RADecCorr = [RADecCorr; dataTemp(5)];
    lat = [lat; dataTemp(6)];
    long = [long; dataTemp(7)];
    alt = [alt; dataTemp(8)];
end

data = [RA Dec];

% Formation of domain covariance matrix, Px
for i = 1:nObs
    block = [RAsigma(i)^2 RADecCorr(i)*RAsigma(i)*Decsigma(i);
        RADecCorr(i)*RAsigma(i)*Decsigma(i) Decsigma(i)^2];
    Px = blkdiag(Px, block);
end

% Check to see if Px is positive-definite
[a, b] = chol(Px);

% Should be all the data in the file. Perform Error check to make sure

% This conditional only checks to make sure the input file is of the
% correct length. 
if ~feof(fid)
    set(handles.Error, 'String', ['ERROR: Input file not formatted',...
        ' correctly. Please exit and re-format the file.']);
    r2 = NaN;
    v2 = NaN;
    P = NaN;
    return;
end
fclose(fid);

%% Formation of site position vectors using either given cosine matrix or LST approximations
rsite = [];
if nargin == 7
    rsiteECEF = [];
    for i = 1:nObs
        rsiteECEF = [rsiteECEF Site_Position(lat(i), long(i), alt(i))];
    end
    rsite = cosMatrix * rsiteECEF;
elseif strcmpi(muType, 'Sun') || strcmpi(muType, 'Other')
    % Do not need LST approx. when given exact heliocentric state vector of
    % the observer or state vector relative to different orbiting body
    for i = 1:nObs
        rsite = [rsite Site_Position(lat(i), long(i), alt(i))];
    end
else
    LST_values = [];
    for i = 1:nObs
        LST_values = [LST_values, LST(time(i,:),long(i))];
        rsite = [rsite Site_Position(lat(i),LST_values(i),alt(i))];
    end
end

%% Calculation of Covariance Matrix, P

% Domain mean
x = [];
for i = 1:nObs
    x = [x; data(i,:)'];
end

% Function handles for Laplace, Gauss, and Double-r
L = @Laplace;
G = @Gauss;
D = @Double_r;

%% Calculation of observation position and velocity vectors
if strcmpi(calcMethod, 'Laplace')
    if b
        [r2, v2] = Laplace(data, time, rsite);
        P = NaN;
        Iterations = NaN;
    else
        [y, P] = sigpt(x, Px, L, time, rsite);
        r2 = y(1:3);
        v2 = y(4:6);
    end
    KOE = kepel(r2, v2, muglobal);
elseif strcmpi(calcMethod, 'Gauss')
    if b
        [r2, v2, Iterations, r1, r3] = Gauss(data, time, rsite);
        P = NaN;
    else
        [y, P] = sigpt(x, Px, G, time, rsite);
        r2 = y(1:3);
        v2 = y(4:6);
        Iterations = NaN;
      
    end
    KOE = kepel(r2, v2, muglobal);

elseif strcmpi(calcMethod, 'Double-r')
    if b
        [r2, v2, Q, Iterations] = Double_r(data, time, rsite);
        P = NaN;

    else
        [y, P] = sigpt(x, Px, D, time, rsite);
        r2 = y(1:3);
        v2 = y(4:6);
        Iterations = NaN;
    end
    KOE = kepel(r2, v2, muglobal);
    
elseif strcmpi(calcMethod, 'Casotto')
    if b
        [r2, v2] = Casotto(data, time, rsite);
        P = NaN;
    else
        [y, P] = sigpt(x, Px, D, time, rsite);
        r2 = y(1:3);
        v2 = y(4:6);
    end
else
    beep;
    errordlg(['Invalid calculation method entered. Please choose',...
        ' Laplace, Gauss or Double-r.'],'Input Error','modal');
    fprintf(1, ['WARNING: Invalid calculation method entered.',...
        ' Please choose Laplace, Gauss or Double-r.\n']);
    r2 = NaN;
    v2 = NaN;
    P = NaN;
    return;
end
end