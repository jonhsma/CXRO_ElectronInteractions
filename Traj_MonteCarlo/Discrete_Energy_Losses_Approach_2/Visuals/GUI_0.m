function varargout = GUI_0(varargin)
%GUI_0 M-file for GUI_0.fig
%      GUI_0, by itself, creates a new GUI_0 or raises the existing
%      singleton*.
%
%      H = GUI_0 returns the handle to a new GUI_0 or the handle to
%      the existing singleton*.
%
%      GUI_0('Property','Value',...) creates a new GUI_0 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to GUI_0_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GUI_0('CALLBACK') and GUI_0('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GUI_0.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_0

% Last Modified by GUIDE v2.5 15-Nov-2018 11:37:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_0_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_0_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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
end

% --- Executes just before GUI_0 is made visible.
function GUI_0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for GUI_0
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

    % Initialize the scan archive
    global scanArchiveGUI
    if nargin >= 4
        scanArchive     =    varargin{1};
    end

    % Check if scanArchive is available
    if exist('scanArchive','var')
        activate(hObject, eventdata, handles)    
        scanArchiveGUI  =   scanArchive;    
    else
        deactivate(hObject, eventdata, handles)
    end
    
end

% UIWAIT makes GUI_0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_0_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%% Nevigation stuff

% --- Executes on button press in pushbutton1.
function buttonVisualize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    plotTraj(hObject, eventdata, handles);
end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
    contents = cellstr(get(hObject,'String'));
    command = contents{get(hObject,'Value')};
    eval(command);
    
end



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
end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
    if isnan(str2double(get(hObject,'String')))||...
        str2double(get(hObject,'String')) < 1
        set(hObject,'String','1')
    end
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    set(hObject,'String','Trial Number');
end
end


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
    if isnan(str2double(get(hObject,'String')))||...
        str2double(get(hObject,'String')) < 0
        set(hObject,'String','0')
    end
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    set(hObject,'String','Incidence Number');
end
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton1.
function pushbutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    plotHierachy(hObject, eventdata, handles);
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
    if isnan(str2double(get(hObject,'String')))||...
        str2double(get(hObject,'String')) < 1
        set(hObject,'String','1')
    end
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


% --- Executes on key press with focus on pushbutton1 and none of its controls.
function pushbutton1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
    plotHierachy(hObject, eventdata, handles);
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.edit3,'String',num2str(max(str2double(get(handles.edit3,'string')-1),1)))
    buttonVisualize_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.edit3,'String',num2str(str2double(get(handles.edit3,'string'))+1))
    buttonVisualize_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.edit1,'String',num2str(max(str2double(get(handles.edit1,'string'))-1,1)))
    buttonVisualize_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.edit1,'String',num2str(str2double(get(handles.edit1,'string'))+1))
    buttonVisualize_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.edit2,'String',num2str(max(str2double(get(handles.edit2,'string'))-1,0)))
    buttonVisualize_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    set(handles.edit2,'String',num2str(str2double(get(handles.edit2,'string'))+1))
    buttonVisualize_Callback(hObject, eventdata, handles)
end

%% File Input
% --- Executes on button press in pushbutton9.
function loadButton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if isempty(handles)
        disp('test')
        GUI_0;
    end
    loadArchive(hObject, eventdata,handles);
end

function result = activate(hObject, eventdata, handles)
    result = set(handles.pushbutton1,'string','Visualize/Hierachy',...
        'enable','on',...
        'ForegroundColor',[0 0 0]);
end

function result = deactivate(hObject, eventdata, handles)
    result = set(handles.pushbutton1,'string','No valid database',...
        'enable','off',...
        'ForegroundColor', [0.6 0.6 0.6]);
end



%% Custom functions
function y = plotTraj(hObject, eventdata, handles)
global scanArchiveGUI
    try
        currEnergy  = round(str2double(get(handles.edit3,'string')));
        currTrial   = round(str2double(get(handles.edit1,'string')));
        currInc     = round(str2double(get(handles.edit2,'string')));
        
        hold off
        if size(scanArchiveGUI,1)> 1
            if currInc > 0
                testTraj    =   scanArchiveGUI{currEnergy,currTrial}.incidences{currInc};
                displayTrajectory(testTraj,'energy',[],[],80,jet(100))
            else
                testTraj    =   scanArchiveGUI{currEnergy,currTrial}.incidences{1};
                for ii = 2:length(scanArchiveGUI{currEnergy,currTrial}.incidences)
                    testTraj = [testTraj scanArchiveGUI{currEnergy,currTrial}.incidences{ii}];
                end
                displayTrajectory(testTraj,'energy',...
                    scanArchiveGUI{currEnergy,currTrial}.acid_xyz,...
                    scanArchiveGUI{currEnergy,currTrial}.activation_xyz,...
                    80,jet(100))
            end
        else
            if currInc > 0
                testTraj    =   scanArchiveGUI{currEnergy}{currTrial}.incidences{currInc};
                displayTrajectory(testTraj,'energy',[],[],80,jet(100))
            else
                testTraj    =   scanArchiveGUI{currEnergy}{currTrial}.incidences{1};
                for ii = 2:length(scanArchiveGUI{currEnergy,currTrial}.incidences)
                    testTraj = [testTraj scanArchiveGUI{currEnergy}{currTrial}.incidences{ii}];
                end
                displayTrajectory(testTraj,'energy',...
                    scanArchiveGUI{currEnergy}{currTrial}.acid_xyz,...
                    scanArchiveGUI{currEnergy}{currTrial}.activation_xyz,...
                    80,jet(100))
            end
        end        
        colormap(jet)
        caxis([0 80])
        daspect([1 1 1])
        colorbar;
        rotate3d on;
        
        xlabel('x(nm)')
        ylabel('y(nm)')
        zlabel('z(nm)')
        
        
        y = 1;
    catch exception
        y = exception;
    end
end

function y = plotHierachy(hObject, eventdata, handles)
global scanArchiveGUI
    try
        currEnergy  = round(str2double(get(handles.edit3,'string')));
        currTrial   = round(str2double(get(handles.edit1,'string')));
        currInc     = round(str2double(get(handles.edit2,'string')));

        if size(scanArchiveGUI,1)> 1
            if currInc > 0
                testTraj    =   {scanArchiveGUI{currEnergy,currTrial}.incidences{currInc}};
            else
                testTraj    =   scanArchiveGUI{currEnergy,currTrial}.incidences;
            end
        else
            if currInc > 0
                testTraj    =   {scanArchiveGUI{currEnergy}{currTrial}.incidences{currInc}};
            else
                testTraj    =   scanArchiveGUI{currEnergy}{currTrial}.incidences;
            end
        end

        hold off
        for jj = 1:length(testTraj)
            displayTrajHierachy(testTraj{jj})
            hold on
        end
        daspect([1 1 1])
        xlabel('x(nm)')
        ylabel('y(nm)')
        zlabel('z(nm)')
        rotate3d on;        
        y = 1;
    catch exception
        y = exception;
    end
end

function result = loadArchive(hObject, eventdata, handles)
global scanArchiveGUI

    [fileName, filePath]    =    uigetfile('*.mat');
    fullFilePath = fullfile(filePath,fileName);
    if exist(fullFilePath,'file')
        set(hObject,'string','Loading...',...
            'enable','off',...
            'ForegroundColor', [0.6 0.6 0.6]);
        drawnow;
        load(fullFilePath)
        set(hObject,'string','Load database',...
            'enable','on',...
            'ForegroundColor', [0 0 0])
        activate(hObject,eventdata, handles)
        scanArchiveGUI  = scanArchive;
        result = 1;
    else
        warndlg('File does not exist');
        if ~exist('scanArvhiceGUI','var')
            deactivate(hObject,eventdata, handles)
            result = -1;
        else
            results = 1;
        end
    end
    
end
