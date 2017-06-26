% ODTBX startup.m file

% ODTBX: Orbit Determination Toolbox
% 
% Copyright (c) 2003-2014 United States Government as represented by the
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

% This file sets the paths for ODTBX use.  Read the sections below to
% properly configure your setup.  (End-users who need their own startup.m
% customizations should append them to the end of this file.)
%
% This file serves three different purposes:
% 1) It is pulled from source control for the nightly regression tests.
%        The regression tests assume a particular directory layout for the
%        clean JAT build and ODTBX tests.  The regression driver will set
%        the basePath variable in the environment.  If basePath is defined,
%        this file assumes it is part of a regression test and ignores
%        other setup types.
% 2) It is the basis for the installer to place in the end-user's installation.
%        The installer will modify this file as it is placed on the 
%        user's system.
% 3) It is used by developers to describe their particular setup.
%        Each developer will have a different set of paths for the various
%        components and will have to modify this file accordingly.
%
% Each setup describes the following for the rest of the file to process:
%     odtbxPath         The path to the ODTBX folder
%     odtbxMicePath     The path to the MICE folder
%     jatSrcPath        JAT source and data files (and maybe the bytecode)
%     jatBytecodePath   (optional) JAT bytecode path for some builds

%
%% DEVELOPER SETUP:
%
% Developers, you will have to modify this section before you start development
% on ODTBX.  Don't modify other parts of this file unless you're sure how
% the changes will affect all developers, the installer, and the
% regression test setup.
%
% DO NOT CHECK IN THESE CHANGES!
% They are for your own development working copy and their presence may 
% break the installer and regression tests.  The rest of this file keys off
% of these either being commented out or set.
%
% Also, if you are changing between multiple installations on the same
% machine, uncommenting the following lines may help avoid conflicts with
% existing installations, old paths, and other ODTBX development trees.
% These will blow away the MATLAB path to a clean state and populate the
% new path out of this file.
%
% clear java         % Clears Java classes
% restoredefaultpath % Restores the default Matlab path
% javaclasspath({})  % Clears dynamic java class path (leaves static path)
%
%
% Developers, 
% modify these paths for your system.  They must be absolute paths, and 
% don't include the trailing slash.
%
% The path to the top-level ODTBX folder (containing odtbx, vendor, ...)
myBasePath = '/Users/ravidavi/Development/odtbx-git';
%
% The path to the ODTBX folder (the directory containing ODTBX_Source)
% e.g., the clone of ssh://git.code.sf.net/p/odtbx/git
% UNCOMMENT THIS if you are use the ODTBX repository instead of an
% installed version of ODTBX.
%odtbxPath = [myBasePath, '/odtbx'];
%
% The path to a folder containing an installed version of ODTBX. Change
% this if you want to use JAT & Mice from an installed version of
% ODTBX, while keeping the actual ODTBX code from odtbxPath. Some devs do
% this so that they don't have to build JAT.
myInstPath = myBasePath;
%
% The path to a specific MICE folder (the top-level folder of a MICE 
% distribution, usually named "mice", and it contains directories like:
% data, doc, include, lib, src, etc.
odtbxMicePath = [myInstPath, '/vendor/MacIntelMice_64/mice'];
%
% The path to the JAT folder that contains the source and data files
% (the directory above jat/)
% Some developers may choose to build their bytecode in these folders as well.
jatSrcPath =  [myInstPath, '/vendor/Jat'];
%
% (Optional) The path to the JAT folder that contains the Java bytecode
% (*.class files).  Most developers won't need this if they want to place
% their bytecode in the same folders as their source location, above.
% Only uncomment and set this if you know how you are compiling JAT.
%jatBytecodePath = '';

%
%% DEVELOPER SETUP ENDS  (That's it for a standard developer setup.)
%

%
%% GMAT WITH ODTBX:
%
% The path to the top-level folder of a GMAT distribution, usually named
% "GMAT". It contains directories like: bin, matlab, plugins, etc.
%gmatPath = '/absolute/path/to/GMAT';

%
%% CHECK STARTUP TYPE:
%

% First, check if we've got regression test going on or not.
startuptype = 0; % unknown
if exist('basePath','var') && ~isempty(basePath)
    startuptype = 1; % regression (trumps all others)
elseif exist('odtbxPath','var') && ~isempty(odtbxPath)
    startuptype = 2; % dev setup, above
end
% Note, startuptype = 3 is for an installed env

%
%% INSTALLER SETUP:
%
% This block is used by the installer script and assumes that the entire
% developer block above is still commented out and that no basePath is
% being passed in.
% Installer directory layout:
% (Keep this coordinated with the installer .xml files!)
% baseIstPath
%   ODTBX_Examples
%   ODTBX_Source
%       JAT_Adapters
%   Regression_Validation
%   Jat
%       jat (all source, data, and bytecode files are assumed here)
%   mice (the appropriate Mice arch version is installed here)
%       doc
%       src
%       lib
%       kernels, (etc)
%   gmat (the appropriate GMAT arch version is installed here)
%       bin
%       data
%       docs
%       matlab, (etc)
%       
% Note, the installer may not install all ODTBX directories.  See below.

% If we don't know what type of run this is, check for an installed run.
if startuptype == 0

    % This will be replaced by the installer:
    baseIstPath = '$INSTALL_PATH' ;
    
    % Check to see if the installer replaced baseIstPath, above.
    if length(baseIstPath) == 13 && baseIstPath(1) == '$'
        % Nope, the installer hasn't replaced it.  This isn't a user
        % install.
        startuptype = 0;
    else
        % Yep, the installer replaced it.        
        startuptype = 3; % installed run
        
        odtbxPath = baseIstPath;
        odtbxMicePath = fullfile(baseIstPath,'mice');
        jatSrcPath = fullfile(baseIstPath,'Jat');
        jatBytecodePath = jatSrcPath;
        
        % If gmatPath isn't already specified, then search for GMAT either
        % inside or next to the ODTBX installation directory. Note that
        % this currently doesn't consider multiple GMAT installations.
        if ~exist('gmatPath', 'var')
            % Search inside ODTBX directory
            reldir = baseIstPath; %
            gmatdir = dir(fullfile(reldir, 'GMAT*'));
            
            % Search next to ODTBX directory
            if isempty(gmatdir)
                reldir = fullfile(baseIstPath, '..');
                gmatdir = dir(fullfile(reldir, 'GMAT*'));
            end
            
            % Additional default search paths can be added here
            
            % If GMAT was potentially found, then save its full path.
            % This path will be validated below.
            if ~isempty(gmatdir)
                gmatPath = fullfile(reldir, gmatdir(1).name);
            end
            
            clear reldir gmatdir; % Cleanup
        end
    end
    
    clear baseIstPath;
end

%
%% INSTALLER SETUP ENDS
%

% Environment type check and error message:
if startuptype == 0
    error('ODTBX:startup','ODTBX startup.m: Unable to determine the ODTBX startup environment. \nIf you are a user please contact the ODTBX development team.  This is likely an installer issue.\nIf you are a developer you will need to modify startup.m appropriately.  See the file comments.');
end

%
%% REGRESSION SETUP:
%
% This block is used by the regression test setup.  Its assumes that
% basePath is being passed in. 
% Regression directory layout:
% (Keep this coordinated with regressionTest.sh!)
%
% odtbx (basePath)
%   ODTBX_Source
%   Regression_Validation
% vendor
%   Jat
%     jat (source, perhaps data, but JAT has frequent hardcoded assumptions)
%     maven
%       target
%           classes (bytecode, data files must be here too for JAT's assumptions)
%   LinuxMice_64
%     mice
%       doc, src, lib, kernels, etc...
%
if startuptype == 1
    odtbxPath = basePath;
    odtbxMicePath = fullfile(basePath,'..','vendor','LinuxMice_64','mice');
    jatSrcPath = fullfile(basePath,'..','vendor','Jat');
    jatBytecodePath = fullfile(basePath,'..','vendor','Jat','maven','target','classes');
    gmatPath = fullfile(basePath,'..','vendor','GMAT_Linux64','GMAT');
end
%
%% REGRESSION SETUP END
%

%
% Finally, set all the paths for any setup case.
%

% Set Matlab paths for OD Toolbox
addpath(fullfile(odtbxPath,'ODTBX_Source'));
addpath(fullfile(odtbxPath,'ODTBX_Source','JAT_Adapters'));
addpath(fullfile(odtbxPath,'ODTBX_Source','GMAT_Adapters'));
addpath(fullfile(odtbxPath,'ODTBX_Source','Integrators'));
addpath(fullfile(odtbxPath,'ODTBX_Source','meas_sched_gui'));
addpath(fullfile(odtbxPath,'ODTBX_Source','GPS_Time_Utilities'));
addpath(fullfile(odtbxPath,'ODTBX_Source','IOD'));
addpath(fullfile(odtbxPath,'ODTBX_Examples'));
addpath(fullfile(odtbxPath,'HTML'));

% Set these Matlab paths for OD Toolbox
% However, an end-user install may not provide these directories.
% Check for their existence when it is an end-user install.
d = fullfile(odtbxPath,'Regression_Validation');
if startuptype ~= 3 || isdir(d)
    addpath(d);
end
d = fullfile(odtbxPath,'Regression_Validation','Release2');
if startuptype ~= 3 || isdir(d)
    addpath(d);
end
d = fullfile(odtbxPath,'Regression_Validation','JAT_Adapters_Tests');
if startuptype ~= 3 || isdir(d)
    addpath(d);
end
d = fullfile(odtbxPath,'Regression_Validation','DataFiles');
if startuptype ~= 3 || isdir(d)
    addpath(d);
end

% Make sure MICE version is correct for MATLAB version
miceversiontxt = fileread(fullfile(odtbxMicePath,'doc','version.txt'));
miceversionstart = strfind(miceversiontxt, 'V.N');
miceversionnum = str2num(miceversiontxt(miceversionstart+3:miceversionstart+6));
mlver = ver('MATLAB');
mlver = str2num(mlver.Version);
if((mlver >= 9.2) && (miceversionnum < 66))
    warning('ODTBX:MICEVersion', ['Using MATLAB R2017a or greater, but ' ...
            'MICE version is not 66. Please download MICE v0066 or greater ' ...
            'from <a href="http://sourceforge.net/projects/odtbx/files">ODTBX on SourceForge</a>.']);
elseif((mlver < 9.2) && (miceversionnum >= 66))
    warning('ODTBX:MICEVersion', ['Using MATLAB R2016b or earlier, but ' ...
            'MICE version is 66 or greater. Please use ODTBX installer to ' ...
            'install MICE v0064']);
end
clear miceversiontxt miceversionstart miceversionnum mlver;

% Set Matlab paths for MICE:
addpath(fullfile(odtbxMicePath,'doc'));
addpath(fullfile(odtbxMicePath,'doc','html'));
addpath(fullfile(odtbxMicePath,'doc','html','mice'));
addpath(fullfile(odtbxMicePath,'src','mice'));
addpath(fullfile(odtbxMicePath,'lib'));
addpath(fullfile(odtbxMicePath,'kernels'));

% Set Matlab paths for JAT
addpath(jatSrcPath);
addpath(fullfile(jatSrcPath,'jat'));

% Set Java paths for JAT 
javaaddpath(jatSrcPath);
javaaddpath(fullfile(jatSrcPath,'jat'));
if exist('jatBytecodePath','var') && ~isempty(jatBytecodePath)
    javaaddpath(jatBytecodePath);
end

% Set Matlab paths for GMAT
if exist('gmatPath', 'var')
    fprintf('Adding GMAT to MATLAB path: %s\n', gmatPath);
    addpath(fullfile(gmatPath,'bin'));
    addpath(fullfile(gmatPath,'matlab','libCInterface'));
end

if(startuptype == 1)
    % Recursively generate a path string containing the ODTBX_Data directory
    % and everything below it.
    p = genpath([basePath filesep 'ODTBX_Data']);
else
    % Recursively generate a path string containing the ODTBX_Data directory
    % and everything below it.
    p = genpath([odtbxPath filesep 'ODTBX_Data']);
end

% genpath returns a path separator after each directory, 
% find out where each path separator appears in p, 
% then use this to check each path, 
% test it and add it to the MATLAB path if it passes the test.
if ~isempty(p)
    pseps = strfind(p,pathsep);
    for i = 1:length(pseps)
        % get a single path to check
        if i == 1
            p2check = p(1:pseps(1));
        else
            p2check = p(pseps(i-1)+1:pseps(i));
        end
        
        % The test: any path that contains a file separator with a leading
        % '.', as in .../.git/..., in it should be ignored (a.k.a the
        % standard UNIX convention).
        if isempty(strfind(p2check,[filesep '.']))
            % no "/." or "\.", etc. found so it is good:
            addpath(p2check);
        end
    end
end

% Clear the remnants out of the workspace to leave only the path setttings.
clear d i p p2check pseps ;
clear basePath jatSrcPath jatBytecodePath miceTopPath miceArch gmatPath startuptype odtbxMicePath odtbxPath myBasePath myInstPath;

%
%% End-User Customization Goes Below Here:
% (Note this file may be overwritten and these changes lost when 
% uninstalling or re-installing ODTBX.)
%
