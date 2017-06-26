function varargout = IODmainGUI(varargin)
% IODMAINGUI MATLAB code for main IOD GUI
%      IODMAINGUI, by itself, creates a new IODMAINGUI or raises the existing
%      singleton*.
%
%      H = IODMAINGUI returns the handle to a new IODMAINGUI or the handle to
%      the existing singleton*.
%
%      IODMAINGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IODMAINGUI.M with the given input arguments.
%
%      IODMAINGUI('Property','Value',...) creates a new IODMAINGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IODmainGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IODmainGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IODmainGUI

% Last Modified by GUIDE v2.5 03-Aug-2015 11:09:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IODmainGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @IODmainGUI_OutputFcn, ...
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


% --- Executes just before IODmainGUI is made visible.
function IODmainGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IODmainGUI (see VARARGIN)

% Choose default command line output for IODmainGUI
handles.output = hObject;
% Ensure that the GUI window sizes properly for all MATLAB versions
set(0,'ScreenPixelsPerInch',96);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IODmainGUI wait for user response (see UIRESUME)
% uiwait(handles.IODMain);


% --- Outputs from this function are returned to the command line.
function varargout = IODmainGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function IODMain_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to IODMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
contents = cellstr(get(hObject,'String'));
datatype = contents{get(hObject, 'Value')};
handles.datatype = datatype;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
datatype = contents{get(hObject, 'Value')};
handles.datatype = datatype;
guidata(hObject, handles);

% --- Executes on mouse press over figure background.
function IODMain_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to IODMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function upload_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upload_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% setappdata(hObject,'filename', 0);
filename = 0;
handles.filename = filename;
guidata(hObject, handles);

% --- Executes on button press in upload_file.
function upload_file_Callback(hObject, eventdata, handles)
% hObject    handle to upload_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname] = uigetfile('*.*', 'File Selection');
if filename ~=0
    set(handles.text3, 'String', strcat(pathname,filename));
    set(handles.error_msg, 'String', '');
    filename = strcat(pathname, filename);
else
    set(handles.text3, 'String', '');
end
handles.filename = filename;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function nasa_logo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nasa_logo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
I = imread('NASA-Logo.jpg');
imhandle = imshow(I);
% Hint: place code in OpeningFcn to populate nasa_logo


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Type
contents = get(handles.orbitType,'String');
if strcmpi(contents(get(handles.orbitType,'Value')),'Heliocentric')
    Type = 'Sun';
elseif strcmpi(contents(get(handles.orbitType,'Value')),'Geocentric')
    Type = 'Earth';
else
    Type = 'Other';
end

if handles.filename == 0
    error = 'Please upload a file to continue!';
    set(handles.error_msg, 'String', error);
elseif strcmpi(Type, 'Other') && (isempty(get(handles.inputRadius,'String')) || isempty(get(handles.inputMu,'String')))
    beep;
    errordlg('Please enter a radius and gravitational parameter value to continue!','Input error','modal');
elseif strcmpi(Type, 'Other') && (isnan(str2double(get(handles.inputRadius,'String'))) ||...
        isnan(str2double(get(handles.inputMu,'String'))) || (str2double(get(handles.inputRadius,'String')) <= 0) ||...
        (str2double(get(handles.inputMu,'String')) <= 0))
    beep;
    errordlg('Radius and gravitational parameter must be positive numbers greater than zero.', 'Input Error', 'modal');
else
       
    % Selects the appropriate second GUI to pull up based on user datatype
    % input
    if strcmp(handles.datatype, 'Range, Azimuth and Elevation')
    elseif strcmp(handles.datatype, 'Angles-only')
        % Check to see if file is formatted correctly before passing to
        % calculation and visualization GUI
        fid = fopen(handles.filename);
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
                errordlg(['File not formatted properly. Please check the ',...
                    'application documentation for data file formatting guidelines.'], 'File Error','modal');
                return;
            end
            
            % Counting the number of lines/observations
            nObs = nObs + 1; 
        end 
        fclose(fid);
        handles.Radius = str2double(get(handles.inputRadius, 'String'));
        handles.mu = str2double(get(handles.inputMu, 'String'));
        guidata(hObject, handles);
        IODmainGUIAngles;    
    else
        % Case of Range and Range-Rate processing
    end
end
    


% --- Executes on selection change in orbitType.
function orbitType_Callback(hObject, eventdata, handles)
% hObject    handle to orbitType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns orbitType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from orbitType
contents = cellstr(get(hObject,'String'));
orbittype = contents{get(hObject, 'Value')};
global Type;
if strcmpi(orbittype, 'Geocentric')
    Type = 'Earth';
    set(handles.text6, 'Visible', 'Off');
    set(handles.text8, 'Visible', 'Off');
    set(handles.inputRadius, 'Visible', 'Off');
    set(handles.inputMu, 'Visible', 'Off');
elseif strcmpi(orbittype, 'Heliocentric')
    Type = 'Sun';
    set(handles.text6, 'Visible', 'Off');
    set(handles.text8, 'Visible', 'Off');
    set(handles.inputRadius, 'Visible', 'Off');
    set(handles.inputMu, 'Visible', 'Off');
else
    Type = 'Other';
    set(handles.text6, 'Visible', 'On');
    set(handles.text8, 'Visible', 'On');
    set(handles.inputRadius, 'Visible', 'On');
    set(handles.inputMu, 'Visible', 'On');
    beep;
    errordlg(['Please enter the radius of the gravitational body and a gravitational parameter'...
        ' in the boxes that have appeared below the "Orbit Type" panel.'],'Input Required','modal');
    
end
handles.orbittype = Type;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function orbitType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to orbitType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    
end
global Type;
Type = 'Earth';
handles.orbittype = Type;
guidata(hObject, handles);

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function IODMain_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to IODMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clearvars -global;
try
    close(findobj('Tag', 'settingsGUI'));
    close(findobj('Tag', 'IODAngles'));
catch ME
    if findobj('Tag','settingsGUI')
        close(findobj('Tag','settingsGUI'));
    end
    if findobj('Tag', 'IODAngles')
        close(findobj('Tag', 'IODAngles'));
    end
end



% --- Executes on button press in FileFormatGuidelinesButton.
function FileFormatGuidelinesButton_Callback(hObject, eventdata, handles)
% hObject    handle to FileFormatGuidelinesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try 
    open(fullfile(fileparts(pwd),'Documentation','IOD_File_Formatting_Guidelines.html'))
catch ME
    try
        open('IOD_File_Formatting_Guidelines.html');
    catch ME
        beep;
        errordlg('Unable to open file formatting guidelines.', 'Error', 'modal');
    end
end



function inputMu_Callback(hObject, eventdata, handles)
% hObject    handle to inputMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputMu as text
%        str2double(get(hObject,'String')) returns contents of inputMu as a double


% --- Executes during object creation, after setting all properties.
function inputMu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputMu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function inputRadius_Callback(hObject, eventdata, handles)
% hObject    handle to inputRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inputRadius as text
%        str2double(get(hObject,'String')) returns contents of inputRadius as a double


% --- Executes during object creation, after setting all properties.
function inputRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
