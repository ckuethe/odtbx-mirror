function varargout = IODmainGUISettings(varargin)
% IODMAINGUISETTINGS MATLAB code for IOD angles-only settings GUI
%      IODMAINGUISETTINGS, by itself, creates a new IODMAINGUISETTINGS or raises the existing
%      singleton*.
%
%      H = IODMAINGUISETTINGS returns the handle to a new IODMAINGUISETTINGS or the handle to
%      the existing singleton*.
%
%      IODMAINGUISETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IODMAINGUISETTINGS.M with the given input arguments.
%
%      IODMAINGUISETTINGS('Property','Value',...) creates a new IODMAINGUISETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IODmainGUISettings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IODmainGUISettings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IODmainGUISettings

% Last Modified by GUIDE v2.5 03-Aug-2015 14:39:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IODmainGUISettings_OpeningFcn, ...
                   'gui_OutputFcn',  @IODmainGUISettings_OutputFcn, ...
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


% --- Executes just before IODmainGUISettings is made visible.
function IODmainGUISettings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IODmainGUISettings (see VARARGIN)

% Choose default command line output for IODmainGUISettings
handles.output = hObject;
% Ensure that the GUI window sizes properly for all MATLAB versions
set(0,'ScreenPixelsPerInch',96);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IODmainGUISettings wait for user response (see UIRESUME)
% uiwait(handles.settingsGUI);


% --- Outputs from this function are returned to the command line.
function varargout = IODmainGUISettings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in okButton.
function okButton_Callback(hObject, eventdata, handles)
% hObject    handle to okButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Do a bunch of stuff to send information to the methods
global IODSettings;
if handles.Qminvalue <= 0 || isnan(str2double(get(handles.Qmin,'String')))
    beep;
    errordlg('Qmin must be a positive number greater than zero.', 'Input Error','modal');
elseif handles.h2 <= 0 || isnan(str2double(get(handles.rangeStepSize,'String')))
    beep;
    errordlg('Step size must be a positive number greater than zero.', 'Input Error', 'modal');
elseif handles.h2reducFactor <= 0 || isnan(str2double(get(handles.stepSizeReducFactor,'String')))
    beep;
    errordlg('Reduction Factor must be a positive number greater than zero.', 'Input Error', 'modal');
elseif handles.Iterations <= 0 || isnan(str2double(get(handles.doublerIterations,'String')))
    beep;
    errordlg('Number of iterations must be a positive number greater than zero.', 'Input Error', 'modal');
elseif handles.r2o <= 0 || isnan(str2double(get(handles.rangeEstimate,'String')))
    beep;
    errordlg('Initial Range Estimate must be a positive number greater than zero.', 'Input Error', 'modal');
elseif handles.Levels <= 0 || isnan(str2double(get(handles.searchLevels,'String')))
    beep;
    errordlg('Number of Search Levels must be a positive number greater than zero.', 'Input Error', 'modal');
else
    
    IODSettings.r2o = handles.r2o;
    IODSettings.Levels = handles.Levels;
    IODSettings.Iterations = handles.Iterations;
    IODSettings.Qmin = handles.Qminvalue;
    IODSettings.h2reducFactor = handles.h2reducFactor;
    IODSettings.h2 = handles.h2;
    close(handles.settingsGUI);
end
x= 0;

% --- Executes on button press in applyButton.
function applyButton_Callback(hObject, eventdata, handles)
% hObject    handle to applyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global IODSettings;
if handles.Qminvalue <= 0 || isnan(str2double(get(handles.Qmin,'String')))
    beep;
    errordlg('Qmin must be a positive number greater than zero.', 'Input Error','modal');
elseif handles.h2 <= 0 || isnan(str2double(get(handles.rangeStepSize,'String')))
    beep;
    errordlg('Step size must be a positive number greater than zero.', 'Input Error', 'modal');
elseif handles.h2reducFactor <= 0 || isnan(str2double(get(handles.stepSizeReducFactor,'String')))
    beep;
    errordlg('Reduction Factor must be a positive number greater than zero.', 'Input Error', 'modal');
elseif handles.Iterations <= 0 || isnan(str2double(get(handles.doublerIterations,'String')))
    beep;
    errordlg('Number of iterations must be a positive number greater than zero.', 'Input Error', 'modal');
elseif handles.r2o <= 0 || isnan(str2double(get(handles.rangeEstimate,'String')))
    beep;
    errordlg('Initial Range Estimate must be a positive number greater than zero.', 'Input Error', 'modal');
elseif handles.Levels <= 0 || isnan(str2double(get(handles.searchLevels,'String')))
    beep;
    errordlg('Number of Search Levels must be a positive number greater than zero.', 'Input Error', 'modal');
else
    
    IODSettings.r2o = handles.r2o;
    IODSettings.Levels = handles.Levels;
    IODSettings.Iterations = handles.Iterations;
    IODSettings.Qmin = handles.Qminvalue;
    IODSettings.h2reducFactor = handles.h2reducFactor;
    IODSettings.h2 = handles.h2;
end
guidata(hObject, handles)



% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.settingsGUI);



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



function rangeStepSize_Callback(hObject, eventdata, handles)
% hObject    handle to rangeStepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rangeStepSize as text
%        str2double(get(hObject,'String')) returns contents of rangeStepSize as a double
h2 = str2double(get(hObject,'String'));
handles.h2 = h2;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function rangeStepSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rangeStepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global Type;
global IODSettings;
global radius;

if ~isfield(IODSettings, 'h2')
    if strcmpi(Type, 'Earth')
        h2 = 6378.137;
    elseif strcmpi(Type, 'Sun')
        h2 = 0.0047;
    else
        h2 = radius;
    end
else
    h2 = IODSettings.h2;
end
set(hObject, 'String', num2str(h2));
handles.h2 = h2;
guidata(hObject, handles)


function Qmin_Callback(hObject, eventdata, handles)
% hObject    handle to Qmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Qmin as text
%        str2double(get(hObject,'String')) returns contents of Qmin as a double
Qmin = str2double(get(hObject,'String'));
handles.Qminvalue = Qmin;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Qmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Qmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global Type;
global IODSettings;

if ~isfield(IODSettings, 'Qmin')
    if strcmpi(Type, 'Earth')
        Qmin = 150;
    elseif strcmpi(Type, 'Sun')
        Qmin = 3000;
    else
        Qmin = 100;
    end
else
    Qmin = IODSettings.Qmin;
end

set(hObject, 'String', num2str(Qmin));
handles.Qminvalue = Qmin;
guidata(hObject, handles)


function doublerIterations_Callback(hObject, eventdata, handles)
% hObject    handle to doublerIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of doublerIterations as text
%        str2double(get(hObject,'String')) returns contents of doublerIterations as a double
Iterations = str2double(get(hObject,'String'));
handles.Iterations = Iterations;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function doublerIterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to doublerIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global IODSettings;
if ~isfield(IODSettings, 'Iterations')
    Iterations = 500;
else
    Iterations = IODSettings.Iterations;
end
set(hObject, 'String', num2str(Iterations));
handles.Iterations = Iterations;
guidata(hObject, handles)


function stepSizeReducFactor_Callback(hObject, eventdata, handles)
% hObject    handle to stepSizeReducFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stepSizeReducFactor as text
%        str2double(get(hObject,'String')) returns contents of stepSizeReducFactor as a double
h2reducFactor = str2double(get(hObject,'String'));
handles.h2reducFactor = h2reducFactor;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function stepSizeReducFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stepSizeReducFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global IODSettings;
if ~isfield(IODSettings, 'h2reducFactor')
    h2reducFactor = 2;
else
    h2reducFactor = IODSettings.h2reducFactor;
end
set(hObject, 'String', num2str(h2reducFactor));
handles.h2reducFactor = h2reducFactor;
guidata(hObject, handles)



function searchLevels_Callback(hObject, eventdata, handles)
% hObject    handle to searchLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of searchLevels as text
%        str2double(get(hObject,'String')) returns contents of searchLevels as a double
handles.Levels = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function searchLevels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to searchLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global IODSettings;
if ~isfield(IODSettings, 'Levels')
    Levels = 20;
else
    Levels = IODSettings.Levels;
end

set(hObject, 'String', num2str(Levels));
handles.Levels = Levels;
guidata(hObject, handles);


function rangeEstimate_Callback(hObject, eventdata, handles)
% hObject    handle to rangeEstimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rangeEstimate as text
%        str2double(get(hObject,'String')) returns contents of rangeEstimate as a double
r2o = str2double(get(hObject, 'String'));
handles.r2o = r2o;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function rangeEstimate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rangeEstimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global radius;
global Type;
global IODSettings;
if ~isfield(IODSettings, 'r2o')
    if strcmpi(Type, 'Earth')
        r2o = 6378.137;
    elseif strcmpi(Type, 'Sun')
        r2o = 0.75;
    else
        r2o = radius;
        disp(radius)
        disp('hello')
    end
else
    r2o = IODSettings.r2o;
end
set(hObject, 'String', num2str(r2o));
handles.r2o = r2o;
% b = uicontrol;
% s = sprintf('Range estimate at second observation time\nthat is used as the starting value for the\n cone-masking technique in Double_r_range_estimates.m');
% b.TooltipString = s;
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function rangeEstimateText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rangeEstimateText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global Type;
if strcmpi(Type, 'Sun')
    set(hObject, 'String', 'Initial Range Estimate (AU)');
else 
    set(hObject, 'String', 'Initial Range Estimate (km)');
end


% --- Executes on mouse press over figure background.
function settingsGUI_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to settingsGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse motion over figure - except title and menu.
function settingsGUI_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to settingsGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function uitoggletool1_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function text5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global Type;
if strcmpi(Type, 'Sun')
    set(hObject, 'String', 'Range Estimator Step Size (AU)');
else
    set(hObject, 'String', 'Range Estimator Step Size (km)');
end
