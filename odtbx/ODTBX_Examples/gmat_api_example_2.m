% GMAT_API_EXAMPLE_2 Demonstrates using the GMAT Java API
% to calculate dynamics and measurements for use in a squential
% estimator. This version is the same as GMAT_API_ODTBX_filter,
% however most of the initialization is shifted over to the GMAT script

clear variables; clear java; clc; %#ok<CLJAVA>
java.lang.System.gc();

currDir = pwd;

scriptFileName = fullfile(currDir, 'gmat_est.script'); % GMAT script to use

% Initialize GMAT and load script
[myMod, gmatBinPath, result] = load_gmat(scriptFileName);

%% Common setup for Dyn and Dat

% Get the SolarSystem object from GMAT
% Contains data on solar system bodies
ss = myMod.GetDefaultSolarSystem();

%% Set up Dynamics in GMAT

% Prepare the force model to be used for dynamics
fm = myMod.GetODEModel('DefaultProp_ForceModel');
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

%% Set up PropSetup for Measurements in GMAT

prop = myMod.CreatePropSetup('prop');

tfs = myMod.GetConfiguredObject('simData');
tfs = gmat.TrackingFileSet.dynamic_cast(tfs);
tfs.SetPropagator(prop);

%% Run the mission
% Must cd if the ROOT_PATH in gmat_startup_file.txt is a relative path
cd(gmatBinPath)
runResult = myMod.RunMission(); % Initialization for estimation is handled here
cd(currDir)

if runResult < 1
    % No point continuing if RunMission() failed
    error('Error in Moderator::RunMission\nError code:\n%d', runResult);
end

%% Get references to the Spacecraft and GroundStation
simsat = myMod.GetInternalObject('SimSat');
tfs = myMod.GetInternalObject('simData');

% Downcast from GmatBase to correct class
simsat = gmat.Spacecraft.dynamic_cast(simsat);
tfs = gmat.TrackingFileSet.dynamic_cast(tfs);

%% Set initial conditions

% Can get parameters from GMAT
ss = myMod.GetDefaultSolarSystem();
earth = ss.GetBody('Earth');
rEarth = earth.GetEquatorialRadius();

x0 = [7100;0;1300;0;7.35;1]; % Initial state
tSpan = 0:100:90*60; % Time span

P0 = diag([1e-6, 1e-6, 1e-6, 1e-12, 1e-12, 1e-12]); % Initial covariance

%% ODTBX performs estimation

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