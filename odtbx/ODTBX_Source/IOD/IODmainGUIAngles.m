function varargout = IODmainGUIAngles(varargin)
% IODMAINGUIANGLES MATLAB code for IOD angles-only GUI
%      IODmainGUIAngles, by itself, creates a new IODmainGUIAngles or raises the existing
%      singleton*.
%
%      H = IODmainGUIAngles returns the handle to a new IODmainGUIAngles or the handle to
%      the existing singleton*.
%
%      IODmainGUIAngles('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IODmainGUIAngles.M with the given input arguments.
%
%      IODmainGUIAngles('Property','Value',...) creates a new IODmainGUIAngles or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IODmainGUIAngles_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IODmainGUIAngles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IODmainGUIAngles

% Last Modified by GUIDE v2.5 31-Jul-2015 15:12:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IODmainGUIAngles_OpeningFcn, ...
                   'gui_OutputFcn',  @IODmainGUIAngles_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before IODmainGUIAngles is made visible.
function IODmainGUIAngles_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IODmainGUIAngles (see VARARGIN)

% Choose default command line output for IODmainGUIAngles
handles.output = hObject;
% Ensure that the GUI window sizes properly for all MATLAB versions
set(0,'ScreenPixelsPerInch',96);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IODmainGUIAngles wait for user response (see UIRESUME)
% uiwait(handles.IODAngles);


% --- Outputs from this function are returned to the command line.
function varargout = IODmainGUIAngles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function calculationMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calculationMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
calcMethod = 'Laplace';
handles.calcMethod = calcMethod;
guidata(hObject, handles);


% --- Executes when selected object is changed in calculationMethod.
function calculationMethod_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in calculationMethod 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

calcMethod = get(eventdata.NewValue,'String');
handles.calcMethod = calcMethod;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function CalculateButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CalculateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Default time unit
timeUnit = 's';
handles.timeunit = timeUnit;

% Default length unit
lengthUnit = 'km';
handles.lengthunit = lengthUnit;

guidata(hObject, handles);


% --- Executes on button press in CalculateButton.
function CalculateButton_Callback(hObject, eventdata, handles)
% hObject    handle to CalculateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Data Extraction from user's data file
h = findobj('Tag', 'IODMain');
if ~isempty(h)
    g1data = guidata(h);
end

% Clears all previous values in GUI
set(handles.accuracy, 'String','');
set(handles.sma,'String','');
set(handles.ecc, 'String', '');
set(handles.incl, 'String', '');
set(handles.tran, 'String', '');
set(handles.raan, 'String', '');
set(handles.argp, 'String', '');
set(handles.iterationsText, 'String', '');
set(handles.text6, 'String', ''); % Position Vector
set(handles.text7, 'String', ''); % Velocity Vector
set(handles.Error, 'String', 'Calculating...');
pause(1);

global muglobal;
global radius;
global Type;
global IODSettings;
if strcmpi(Type, 'Sun')
    muglobal = 132712400180.6;    % km^3/s^2
    mu = muglobal;
    radius = 695800;    % km
elseif strcmpi(Type, 'Earth')
    muglobal = 398600.4418;    % km^3/s^2
    mu = muglobal;
    radius = 6378.137;
else
    muglobal = g1data.mu;
    mu = muglobal;
    radius = g1data.Radius;
end

% Determining which method to run the input data through

if strcmp(handles.calcMethod, 'Laplace')
    % Initiate calculation with Laplace function
    % Return units should always be km and km/s
    [r2, v2, P] = IODmain(g1data.filename, 'Laplace', Type, IODSettings, muglobal, radius);
    KOE = kepel(r2, v2, mu);
    KOE.period = 2*pi*sqrt(KOE.sma^3/mu);
    set(handles.iterationsText, 'String', 'N/A for Laplace');
elseif strcmp(handles.calcMethod, 'Gauss')
    % Initiate calculation with Gauss function
    % Return units should always be km and km/s coming out
    [r2, v2, P, Iterations] = IODmain(g1data.filename, 'Gauss', Type, IODSettings, muglobal, radius);
    KOE = kepel(r2, v2, mu);
    KOE.period = 2*pi*sqrt(KOE.sma^3/mu);
    set(handles.iterationsText, 'String', num2str(Iterations));
else
    % Initiate calculation with Double-r iteration
    [r2, v2, P, Iterations] = IODmain(g1data.filename, 'Double-r', Type, IODSettings, muglobal, radius);
    KOE = kepel(r2, v2, mu);
    KOE.period = 2*pi*sqrt(KOE.sma^3/mu);
    set(handles.iterationsText, 'String', num2str(Iterations));
end


handles.r2 = r2;     % Always km at this line
handles.v2 = v2;     % Always km/s at this line

% Setting Keplerian Orbital Elements text fields in GUI
if strcmpi(handles.lengthunit, 'ER')
    set(handles.sma, 'String', num2str(KOE.sma/6378.137));
    set(handles.smaUnits, 'String', handles.lengthunit);
elseif strcmpi(handles.lengthunit, 'AU')
    set(handles.sma, 'String', num2str(KOE.sma/149597871));
    set(handles.smaUnits, 'String', handles.lengthunit);
else
    set(handles.sma, 'String', num2str(KOE.sma));
    set(handles.smaUnits, 'String', handles.lengthunit);
end

set(handles.ecc, 'String', num2str(KOE.ecc));
set(handles.incl, 'String', num2str(KOE.incl*180/pi));
set(handles.tran, 'String', num2str(KOE.tran*180/pi));
set(handles.raan, 'String', num2str(KOE.raan*180/pi));
set(handles.argp, 'String', num2str(KOE.argp*180/pi));

% Setting position and velocity vector text fields in GUI
if strcmpi(handles.lengthunit, 'km')
    lconversion = 1;
    if strcmpi(handles.timeunit, 'TU')
        timeconversion = 806.81112382429;
    elseif strcmpi(handles.timeunit, 'Day')
        timeconversion = 86400;
    else
        timeconversion = 1;
    end
elseif strcmpi(handles.lengthunit, 'ER')
    lconversion = 6378.137;
    if strcmpi(handles.timeunit, 'TU')
        timeconversion = 806.81112382429;
    elseif strcmpi(handles.timeunit, 'Day')
        timeconversion = 86400;
    else
        timeconversion = 1;
    end
else
    lconversion = 149597871;
    if strcmpi(handles.timeunit, 'TU')
        timeconversion = 806.81112382429;
    elseif strcmpi(handles.timeunit, 'Day')
        timeconversion = 86400;
    else
        timeconversion = 1;
    end
end

% Check for imaginary solutions
for i = 1:3
    if imag(r2(i)) ~= 0 || imag(v2(i)) ~= 0
        beep;
        errordlg('Error, no real solution found. Recommend adjusting advanced settings.', 'Imaginary Solution', 'modal');
        break;
    end
end

if strcmpi(handles.lengthunit, 'km') && strcmpi(Type, 'Sun')
    set(handles.text6, 'String', num2str(round(handles.r2'/lconversion)));
    set(handles.accuracy, 'String', 'Rounded position vector to nearest integer.')
else
    set(handles.text6, 'String', num2str(handles.r2'/lconversion));
    set(handles.accuracy, 'String','');
end

set(handles.sma, 'String', num2str(KOE.sma/lconversion));
set(handles.text7, 'String', num2str(v2'*timeconversion / lconversion));
set(handles.posUnits, 'String', handles.lengthunit);
set(handles.velUnits, 'String', strcat(handles.lengthunit, '/', handles.timeunit));
set(handles.smaUnits, 'String', handles.lengthunit);
set(handles.Error, 'String', 'Calculation Complete');


handles.KOE = KOE;    % Length units should be km, angles in rad

%% Creating the visualization using EarthOrbitPlot.m from ODTBX
if KOE.ecc < 1 && KOE.ecc > 0 && handles.visual == 1
    
    % First set the number of times to be propagated    
    t = 10:KOE.period/1000:KOE.period;
    KOEf = kepprop2b(KOE, t, mu);    %(ODTBX function)

    Trajectory = [];
    for i = 1:length(KOEf.sma)
        temp.sma = KOEf.sma(i);
        temp.ecc = KOEf.ecc(i);
        temp.incl = KOEf.incl(i);
        temp.raan = KOEf.raan(i);
        temp.argp = KOEf.argp(i);
        temp.tran = KOEf.tran(i);
        Y = kep2cart(temp, mu);    %(ODTBX function)

        % Resizing output matrix for convenience
        Y = [Y(1) Y(4); Y(2) Y(5); Y(3) Y(6)];

        % Creation of Interial Trajectory Matrix for EarthOrbitPlot
        Trajectory = [Trajectory Y(:,1)];
    end
    

    OrbitPlot(Trajectory);
    set(handles.orbitTrajectory,'XColor','k');
    set(handles.orbitTrajectory,'YColor', 'k');
    set(handles.orbitTrajectory,'ZColor', 'k');
    h = findobj('Tag', 'IODAngles');
    g2data = guidata(h);
    set(g2data.IODAngles, 'Color', [0.941 0.941 0.941]);
    
elseif KOE.ecc >= 1 && handles.visual == 1
    beep;
    errordlg('ERROR: Cannot currently plot parabolic or hyperbolic orbits', 'Plotting Error','modal');
else
    cla(handles.orbitTrajectory);
end

guidata(hObject, handles);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function visualizationPanel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to visualizationPanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Default visualization setting: 1 for on, 0 for off
visual = 1;
handles.visual = visual;
guidata(hObject, handles);


% --- Executes when selected object is changed in visualizationPanel.
function visualizationPanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in visualizationPanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
temp = get(eventdata.NewValue,'String');
if strcmpi(temp, 'On')
    visual = 1;
else
    visual = 0;
end
handles.visual = visual;
guidata(hObject, handles);


% --- Executes on selection change in lengthUnit.
function lengthUnit_Callback(hObject, eventdata, handles)
% hObject    handle to lengthUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lengthUnit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lengthUnit
global Type;
contents = cellstr(get(hObject,'String'));
lengthUnit = contents{get(hObject, 'Value')};
smaValue = get(handles.sma, 'String');
if isempty(smaValue)
    if strcmpi(lengthUnit, 'Select an option:')
        handles.lengthunit = 'km';
    else
        handles.lengthunit = lengthUnit;
    end
else
    if strcmpi(lengthUnit, 'km')
        handles.lengthunit = 'km';
        lconversion = 1;
        if strcmpi(handles.timeunit, 'TU')
            timeconversion = 806.81112382429;
        elseif strcmpi(handles.timeunit, 'Day')
            timeconversion = 86400;
        else
            timeconversion = 1;
        end
    elseif strcmpi(lengthUnit, 'ER')
        handles.lengthunit = 'ER';
        lconversion = 6378.137;
        if strcmpi(handles.timeunit, 'TU')
            timeconversion = 806.81112382429;
        elseif strcmpi(handles.timeunit, 'Day')
            timeconversion = 86400;
        else
            timeconversion = 1;
        end
    else
        handles.lengthunit = 'AU';
        lconversion = 149597871;
        if strcmpi(handles.timeunit, 'TU')
            timeconversion = 806.81112382429;
        elseif strcmpi(handles.timeunit, 'Day')
            timeconversion = 86400;
        else
            timeconversion = 1;
        end
    end
    if strcmpi(lengthUnit, 'km') && strcmpi(Type, 'Sun')
        set(handles.text6, 'String', num2str(round(handles.r2'/lconversion)));
        set(handles.accuracy, 'String', 'Rounded position vector to nearest integer.')
    else
        set(handles.text6, 'String', num2str(handles.r2'/lconversion));
        set(handles.accuracy, 'String','');
    end
    
    set(handles.sma, 'String', num2str(handles.KOE.sma/lconversion));
    set(handles.text7, 'String', num2str(handles.v2'*timeconversion / lconversion));
    set(handles.posUnits, 'String', handles.lengthunit);
    set(handles.velUnits, 'String', strcat(handles.lengthunit, '/', handles.timeunit));
    set(handles.smaUnits, 'String', handles.lengthunit);
    
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function lengthUnit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lengthUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.lengthunit = 'km';
guidata(hObject, handles)

% --- Executes on selection change in timeUnit.
function timeUnit_Callback(hObject, eventdata, handles)
% hObject    handle to timeUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns timeUnit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from timeUnit
global Type;
contents = cellstr(get(hObject,'String'));
timeunit = contents{get(hObject, 'Value')};
smaValue = get(handles.sma, 'String');

if isempty(smaValue)
    if strcmpi(timeunit, 'Select an option:')
        handles.timeunit = 's';
    else
        handles.timeunit = timeunit;
    end
else
    if strcmpi(timeunit, 's')
        handles.timeunit = 's';
        timeconversion = 1;
        if strcmpi(handles.lengthunit, 'km')
            lconversion = 1;
        elseif strcmpi(handles.lengthunit, 'ER')
            lconversion = 6378.137;
        else
            lconversion = 149597871;
        end
    elseif strcmpi(timeunit, 'Day')
        handles.timeunit = 'Day';
        timeconversion = 86400;
        if strcmpi(handles.lengthunit, 'km')
            lconversion = 1;
        elseif strcmpi(handles.lengthunit, 'ER')
            lconversion = 6378.137;
        else
            lconversion = 149597871;
        end
    else
        handles.timeunit = 'TU';
        timeconversion = 806.81112382429;
        if strcmpi(handles.lengthunit, 'km')
            lconversion = 1;
        elseif strcmpi(handles.lengthunit, 'ER')
            lconversion = 6378.137;
        else
            lconversion = 149597871;
        end
    end
    if strcmpi(handles.lengthunit, 'km') && strcmpi(Type, 'Sun')
        set(handles.text6, 'String', num2str(round(handles.r2'/lconversion)));
        set(handles.accuracy, 'String', 'Rounded position vector to nearest integer.')
    else
        set(handles.text6, 'String', num2str(handles.r2'/lconversion));
        set(handles.accuracy, 'String','');
    end
    set(handles.sma, 'String', num2str(handles.KOE.sma/lconversion));
    set(handles.text7, 'String', num2str(handles.v2'*timeconversion / lconversion));
    set(handles.posUnits, 'String', handles.lengthunit);
    set(handles.velUnits, 'String', strcat(handles.lengthunit, '/', handles.timeunit));
    set(handles.smaUnits, 'String', handles.lengthunit);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function timeUnit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.timeunit = 's';
guidata(hObject, handles);


% --- Executes on button press in settingsButton.
function settingsButton_Callback(hObject, eventdata, handles)
% hObject    handle to settingsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
IODmainGUISettings


% --- Executes during object deletion, before destroying properties.
function IODAngles_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to IODAngles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clearvars -global;
try
    close(findobj('Tag', 'settingsGUI'));
catch ME
    if findobj('Tag','settingsGUI')
        close(findobj('Tag','settingsGUI'));
    end
end


% --- Executes during object creation, after setting all properties.
function settingsButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to settingsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Error_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'String','');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function uipanel10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

guidata(hObject, handles);
