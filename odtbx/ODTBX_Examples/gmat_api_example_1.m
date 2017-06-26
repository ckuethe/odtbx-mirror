% GMAT_API_EXAMPLE_1 Demonstrates using the GMAT Java API
% to calculate dynamics and measurements for use in a squential estimator

clear variables; clear java; clc; %#ok<CLJAVA>
java.lang.System.gc();

currDir = pwd;

% Initialize GMAT without loading a script
[myMod, gmatBinPath, result] = load_gmat();

% Get the instance of the GMAT File Manager
fileManager = gmat.FileManager.Instance();

%% Common setup for Dyn and Dat

% Get the SolarSystem object from GMAT
% Contains data on solar system bodies
ss = myMod.GetDefaultSolarSystem();

%% Set up Dynamics in GMAT

% Prepare the force model to be used for dynamics
fm = gmat.ODEModel('ODEModel', 'DefaultProp_ForceModel');

fm.SetStringParameter('CentralBody', 'Earth');
fm.SetStringParameter('PrimaryBodies', '{Earth}');
fm.SetStringParameter('Drag', 'None');
fm.SetOnOffParameter('SRP', 'Off');
fm.SetOnOffParameter('RelativisticCorrection', 'Off');
fm.SetStringParameter('ErrorControl', 'RSSStep');

gravityField = gmat.GravityField('Earth', 'Earth');
gravityField.SetIntegerParameter('Degree', 4);
gravityField.SetIntegerParameter('Order', 4);
gravityField.SetStringParameter('PotentialFile', 'JGM2.cof');
gravityField.SetStringParameter('EarthTideModel', 'None');

% EOP file is handled behind the scenes in GMAT script
eopFileName = fileManager.GetFullPathname('EOP_FILE');
eopFile = gmat.EopFile(eopFileName);
gravityField.SetEopFile(eopFile);

% Inform Java that GMAT will handle deleting gravityField instead of the 
% Java garbage collector
gravityField.setSwigOwnership(false);

fm.AddForce(gravityField);

% Create a state for derivatives
% Must be of size (n + (n x n))
% Contains both n d/dt values, and n x n Jacobian values.
state = gmat.GmatState(6+6^2);
fm.SetSolarSystem(ss); % Set solar system pointer in force model
fm.SetState(state); % Provide force model with the state placeholder


sat = gmat.Spacecraft('Sat'); % Create new Spacecraft
% Create PropagationStateManager to manage calculation of derivatives
propManager = gmat.PropagationStateManager();
propManager.SetObject(sat); % Add Spacecraft to PropagationStateManager
propManager.SetProperty('AMatrix', sat); % Want to calculate Jacobian 
propManager.BuildState(); % Fills in state data for PropagationStateManager


% Tell force model to use propmanager
fm.SetPropStateManager(propManager); 
fm.UpdateInitialData(); % Update model to catch changes that were added
fm.BuildModelFromMap(); % Sets up the models in the force model

fm.Initialize(); % Initialize the force model

%% Set up Measurements in GMAT

prop = myMod.CreatePropSetup('prop');

% Get Earth from solar system
earth = ss.GetBody('Earth');

% Create an EarthFixed and EarthMJ2000Eq coordinate system
ef = gmat.CoordinateSystem.CreateLocalCoordinateSystem('EarthFixed','BodyFixed',earth,[],[],earth,ss);
j2k = gmat.CoordinateSystem.CreateLocalCoordinateSystem('EarthMJ2000Eq','MJ2000Eq',earth,[],[],earth,ss);

% Create a Spacecraft for measurements
simsat = gmat.Spacecraft('SimSat');
simsat.SetSolarSystem(ss);
% Set the J2000 body for the Spacecraft
% WARNING: The J2000 body must be identical for all objects in a GMAT run
simsat.SetJ2000BodyName('Earth');
simsat.SetJ2000Body(earth);
simsat.SetInternalCoordSystem(j2k); % Want to use inertial coordinates
simsat.SetRefObject(j2k,gmat.ObjectType.COORDINATE_SYSTEM, 'EarthMJ2000Eq');
simsat.SetStringParameter('Attitude', 'CoordinateSystemFixed');
simsat.SetRealParameter('X', 70000);
simsat.SetRealParameter('Y', 20000);
simsat.SetRealParameter('Z', 5000);

% Create a ground station to send and receive signals
gds = gmat.GroundStation('GDS');
gds.SetSolarSystem(ss);
% Set the J2000 body for the GroundStation
% WARNING: The J2000 body must be identical for all objects in a GMAT run
gds.SetJ2000BodyName('Earth');
gds.SetJ2000Body(earth);
gds.SetInternalCoordSystem(ef); % Want to use Earth fixed coordinates
gds.SetRefObject(earth,gmat.ObjectType.CELESTIAL_BODY, 'EarthFixed');
% Select if the location of the ground station is in 'Cartesian' or in
% 'Spherical' (Lat, Lon, Alt) coordinates
gds.SetStringParameter('StateType', 'Spherical'); 
% If the StateType is Spherical, the HorizonReference needs to be specified
% Options are 'Sphere' or 'Ellipsoid' which accounts for flattening
gds.SetStringParameter('HorizonReference', 'Ellipsoid');
gds.SetRealParameter('Location1', 0.0); % Latitude of Ground station
gds.SetRealParameter('Location2', 0.0); % Longitude of Ground station
gds.SetRealParameter('Location3', 0.0); % Altitude of Ground station


ant1  = gmat.Antenna('Antenna','Antenna1'); % Create station antenna
tmit  = gmat.Transmitter('Transmitter','Transmitter1'); % Create transmitter
rec   = gmat.Receiver('Receiver','Receiver1'); % Create receiver
tmit.SetRealParameter('Frequency',2067.5); % Set transmitter frequency

ant2  = gmat.Antenna('Antenna','Antenna2'); % Create spacecraft antenna
tpond = gmat.Transponder('Transponder','Transponder1'); % Create transponder
tpond.SetStringParameter('TurnAroundRatio','240/221'); % Set TurnAroundRatio

% Use Antenna1 for Transmitter1 and Receiver1
tmit.SetStringParameter('PrimaryAntenna','Antenna1');
rec.SetStringParameter('PrimaryAntenna','Antenna1');
tmit.SetRefObject(ant1,gmat.ObjectType.HARDWARE,'Antenna1');
rec.SetRefObject(ant1,gmat.ObjectType.HARDWARE,'Antenna1');

% Use Antenna2 for Transponder1
tpond.SetStringParameter('PrimaryAntenna','Antenna2');
tpond.SetRefObject(ant2,gmat.ObjectType.HARDWARE,'Antenna2');

% Add Antenna2 and Transponder1 to spacecraft
simsat.SetStringParameter('AddHardware','Antenna2');
simsat.SetStringParameter('AddHardware','Transponder1');
simsat.SetRefObject(ant2,gmat.ObjectType.HARDWARE,'Antenna2');
simsat.SetRefObject(tpond,gmat.ObjectType.HARDWARE,'Transponder1');

% Add Antenna1, Transmitter1, and Receiver1 to station
gds.SetStringParameter('AddHardware','Antenna1');
gds.SetStringParameter('AddHardware','Transmitter1');
gds.SetStringParameter('AddHardware','Receiver1');
gds.SetRefObject(ant1,gmat.ObjectType.HARDWARE,'Antenna1');
gds.SetRefObject(tmit,gmat.ObjectType.HARDWARE,'Transmitter1');
gds.SetRefObject(rec,gmat.ObjectType.HARDWARE,'Receiver1');

% Define range measurements and error model
tem = gmat.ErrorModel('TheErrorModel');
% Specify these measurements are range measurements in km
tem.SetStringParameter('Type','Range');
tem.SetRealParameter('NoiseSigma', 0.050); % Standard deviation of noise
tem.SetRealParameter('Bias',0); % Bias in measurement

% Define doppler range rate measurements and error model
tem2 = gmat.ErrorModel('TheErrorModel2');
% Specify these measurements are doppler range rate measurements
tem2.SetStringParameter('Type','RangeRate');
tem2.SetRealParameter('NoiseSigma', 5e-5); % Standard deviation of noise
tem2.SetRealParameter('Bias',0); % Bias in measurement

% Add ErrorModels to the ground station
gds.SetStringParameter('ErrorModels','TheErrorModel',0);
gds.SetStringParameter('ErrorModels','TheErrorModel2',1);
gds.SetRefObject(tem,gmat.ObjectType.ERROR_MODEL,'TheErrorModel');
gds.SetRefObject(tem2,gmat.ObjectType.ERROR_MODEL,'TheErrorModel2');

% Create a TrackingFileSet to manage the observations
tfs = gmat.TrackingFileSet('Ranger');
tfs.SetStringParameter('FileName','TrkFile_API_GN.gmd'); % in GMAT output directory
tfs.SetBooleanParameter('UseLightTime',false);
tfs.SetBooleanParameter('UseRelativityCorrection',false);
tfs.SetBooleanParameter('UseETminusTAI',false);
tfs.SetSolarSystem(ss);

% Define signal paths and measurement type(s)
% 2-way measurements are used here along the path GDS -> SimSat -> GDS
% Add range measurements to TrackingFileSet
tfs.SetStringParameter('AddTrackingConfig', '{{GDS,SimSat,GDS}, Range}',0);
% Add doppler range rate measurements to TrackingFileSet
tfs.SetStringParameter('AddTrackingConfig', '{{GDS,SimSat,GDS}, RangeRate}',1);
% Add ground station to TrackingFileSet
tfs.SetRefObject(gds,gmat.ObjectType.SPACE_POINT,'GDS');
% Add spacecraft to TrackingFileSet
tfs.SetRefObject(simsat,gmat.ObjectType.SPACECRAFT,'SimSat');
tfs.SetPropagator(prop); % Tell TrackingFileSet the propagator to use

%% Initialize measurement components

gds.Initialize();
simsat.Initialize();
tfs.Initialize();

%% Set initial conditions

% Can get parameters from GMAT
rEarth = earth.GetEquatorialRadius();

x0 = [7100;0;1300;0;7.35;1]; % Initial state
tSpan = 0:100:90*60; % Time span

P0 = diag([1e-6, 1e-6, 1e-6, 1e-12, 1e-12, 1e-12]); % Initial covariance

%%

odeOpts = odeset('reltol', 1e-6, 'abstol', 1e-6); % Numerical integration tolerances
odtbxOpts = setOdtbxOptions('OdeSolvOpts', odeOpts);

% Call sequential estimator
% Need to send force model to dyn function
% Need to send TrackingFileSet and Spacecraft to dat function
[t,xhat,P,e,dy,Pa,Pv,Pw,Phata,Phatv,Phatw,SigSA,Pdy,Pdyt] = ...
    estseq(@gmat_dyn,@gmat_dat,tSpan,x0,P0,odtbxOpts,fm,{tfs,simsat});


%% Plot results
X = xhat{1};

figure(1); clf;
[XX,YY,ZZ] = sphere;
surf(XX*rEarth, YY*rEarth, ZZ*rEarth, 'FaceColor', 'b', 'FaceAlpha', 0.25, 'EdgeAlpha', 0.15);
hold on

plot3(X(1,:), X(2,:), X(3,:), 'r', 'LineWidth', 2)
axis equal


estval(t,e,P,scrunch(Pa+Pv+Pw),gcf);