function varargout = ImageAnalysisGUI(varargin)
% IMAGEANALYSISGUI M-file for ImageAnalysisGUI.fig
%      IMAGEANALYSISGUI, by itself, creates a new IMAGEANALYSISGUI or raises the existing
%      singleton*.
%
%      H = IMAGEANALYSISGUI returns the handle to a new IMAGEANALYSISGUI or the handle to
%      the existing singleton*.
%
%      IMAGEANALYSISGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGEANALYSISGUI.M with the given input arguments.
%
%      IMAGEANALYSISGUI('Property','Value',...) creates a new IMAGEANALYSISGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImageAnalysisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImageAnalysisGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImageAnalysisGUI

% Last Modified by GUIDE v2.5 05-Oct-2009 10:45:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ImageAnalysisGUI_OpeningFcn, ...
    'gui_OutputFcn',  @ImageAnalysisGUI_OutputFcn, ...
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


% --- Executes just before ImageAnalysisGUI is made visible.
function ImageAnalysisGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImageAnalysisGUI (see VARARGIN)

% Choose default command line output for ImageAnalysisGUI
handles.output = hObject;

%initialization
global RAW_IMAGE_FOLDER;
handles.RAW_IMAGE_FOLDER = RAW_IMAGE_FOLDER;

set(handles.projection_bgroup,'SelectionChangeFcn',@calculateProjection);
set(handles.projectionMeanVsMax_bgroup,'SelectionChangeFcn',@calculateProjection);
set(handles.combinedVsSeparate_bgroup,'SelectionChangeFcn',@plotOnSelectionChange);
set(handles.viewChannel_bgroup,'SelectionChangeFcn',@plotOnSelectionChange);
set(handles.units_bgroup,'SelectionChangeFcn',@plotOnSelectionChange);

%init conditions
set(handles.gamma1_edit,'String','1');
set(handles.lowCut1_edit,'String','0');
set(handles.highCut1_edit,'String','1');
set(handles.gamma2_edit,'String','1');
set(handles.lowCut2_edit,'String','0');
set(handles.highCut2_edit,'String','1');

%set(handles.combinedVsSeparate_bgroup,'SelectionChangeFcn',@(hObject,eventData)combinedVsSeparate_SelectionChangeFcn);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImageAnalysisGUI wait for user response (see UIRESUME)
% uiwait(handles.MainPanel);


% --- Outputs from this function are returned to the command line.
function varargout = ImageAnalysisGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in crop_button.
function crop_button_Callback(hObject, eventdata, handles)
% hObject    handle to crop_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
point1 = [];
point2 = [];
while isequal(point1,point2)
    k = waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
end
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions
yRange = round([p1(1) p1(1)+offset(1)]);
xRange = round([p1(2) p1(2)+offset(2)]);

setappdata(handles.MainPanel,'xRange',xRange);
setappdata(handles.MainPanel,'yRange',yRange);
set(handles.xCropRegion_edit,'String',num2str(xRange));
set(handles.yCropRegion_edit,'String',num2str(yRange));

%crop (if not xz projection)
channel1Image = getappdata(handles.MainPanel,'curImage_ch1');
channel2Image = getappdata(handles.MainPanel,'curImage_ch2');
if ~strcmp(get(get(handles.projection_bgroup,'SelectedObject'),'Tag'),'projXZ_radio')    
    xRange = getappdata(handles.MainPanel,'xRange');
    yRange = getappdata(handles.MainPanel,'yRange');    
    if ~isempty(xRange)
        channel1Image = channel1Image(xRange(1):xRange(2),yRange(1):yRange(2));
        channel2Image = channel2Image(xRange(1):xRange(2),yRange(1):yRange(2));
    end
end
setappdata(handles.MainPanel,'curImage_ch1',channel1Image);
setappdata(handles.MainPanel,'curImage_ch2',channel2Image);

plotImage(handles);


% --- Executes on button press in resetCrop_button.
function resetCrop_button_Callback(hObject, eventdata, handles)
% hObject    handle to resetCrop_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

xRange = [1 getappdata(handles.MainPanel,'Xpixels')];
yRange = [1 getappdata(handles.MainPanel,'Ypixels')];

setappdata(handles.MainPanel,'xRange',xRange);
setappdata(handles.MainPanel,'yRange',yRange);
set(handles.xCropRegion_edit,'String',num2str(xRange));
set(handles.yCropRegion_edit,'String',num2str(yRange));
calculateProjection(hObject,eventdata);
plotImage(handles);

function xCropRegion_edit_Callback(hObject, eventdata, handles)
% hObject    handle to xCropRegion_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xCropRegion_edit as text
%        str2double(get(hObject,'String')) returns contents of xCropRegion_edit as a double
xRange = str2num(get(hObject,'String'));
setappdata(handles.MainPanel,'xRange',xRange);
plotImage(handles);

% --- Executes during object creation, after setting all properties.
function xCropRegion_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xCropRegion_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yCropRegion_edit_Callback(hObject, eventdata, handles)
% hObject    handle to yCropRegion_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yCropRegion_edit as text
%        str2double(get(hObject,'String')) returns contents of yCropRegion_edit as a double
yRange = str2num(get(hObject,'String'));
setappdata(handles.MainPanel,'yRange',yRange);
plotImage(handles);

% --- Executes during object creation, after setting all properties.
function yCropRegion_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yCropRegion_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zRange1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to zRange1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zRange1_edit as text
%        str2double(get(hObject,'String')) returns contents of zRange1_edit as a double
calculateProjection(hObject,eventdata);

% --- Executes during object creation, after setting all properties.
function zRange1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zRange1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gamma1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to gamma1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma1_edit as text
%        str2double(get(hObject,'String')) returns contents of gamma1_edit as a double
calculateProjection(hObject,eventdata);

% --- Executes during object creation, after setting all properties.
function gamma1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lowCut1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to lowCut1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowCut1_edit as text
%        str2double(get(hObject,'String')) returns contents of lowCut1_edit as a double
calculateProjection(hObject,eventdata);

% --- Executes during object creation, after setting all properties.
function lowCut1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowCut1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function highCut1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to highCut1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of highCut1_edit as text
%        str2double(get(hObject,'String')) returns contents of highCut1_edit as a double
calculateProjection(hObject,eventdata);

% --- Executes during object creation, after setting all properties.
function highCut1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to highCut1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveImage_button.
function saveImage_button_Callback(hObject, eventdata, handles)
% hObject    handle to saveImage_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.RAW_IMAGE_FOLDER), DefaultName = [handles.RAW_IMAGE_FOLDER '/processed'];
else DefaultName = '~'; end
onscreenImage = getappdata(handles.MainPanel,'onscreenImage');
onscreenImage1 = getappdata(handles.MainPanel,'onscreenImage1');
onscreenImage2 = getappdata(handles.MainPanel,'onscreenImage2');

[FileName,PathName,FilterIndex] = uiputfile('*.tif','SaveFileAs',DefaultName);

if ~isempty(FilterIndex) %if file selected
    if ~isempty(onscreenImage) %just 1
        fname = [PathName FileName];
        imwrite(onscreenImage,fname,'tif');
    elseif ~isempty(onscreenImage1) %2 images        
        fname = [PathName FileName];
        [pathstr, name, ext] = fileparts(fname); 
        fname1 = [pathstr '/' name '_ch1' ext];
        fname2 = [pathstr '/' name '_ch2' ext];
        imwrite(onscreenImage1,fname1,'tif');
        imwrite(onscreenImage2,fname2,'tif');
    end
    
end

% --- Executes on button press in saveMovie_button.
function saveMovie_button_Callback(hObject, eventdata, handles)
% hObject    handle to saveMovie_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in saveMatrix_button.
function saveMatrix_button_Callback(hObject, eventdata, handles)
% hObject    handle to saveMatrix_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function zFrame_slider_Callback(hObject, eventdata, handles)
% hObject    handle to zFrame_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if ~getappdata(handles.MainPanel,'drawInProgress') %if drawing, do nothing
    curZ = round(get(hObject,'Value'));
    set(handles.zFrame_edit,'String',num2str(curZ));
    Image1Raw = getappdata(handles.MainPanel,'Image1Raw');
    setappdata(handles.MainPanel,'curImage_ch1',Image1Raw(:,:,curZ));
    if getappdata(handles.MainPanel,'NumChannels') > 1
        Image2Raw = getappdata(handles.MainPanel,'Image2Raw');
        setappdata(handles.MainPanel,'curImage_ch2',Image2Raw(:,:,curZ));
    end
    plotImage(handles);
end

% --- Executes during object creation, after setting all properties.
function zFrame_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zFrame_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function zFrame_edit_Callback(hObject, eventdata, handles)
% hObject    handle to zFrame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zFrame_edit as text
%        str2double(get(hObject,'String')) returns contents of zFrame_edit as a double


% --- Executes during object creation, after setting all properties.
function zFrame_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zFrame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function myPixelX_edit_Callback(hObject, eventdata, handles)
% hObject    handle to myPixelX_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of myPixelX_edit as text
%        str2double(get(hObject,'String')) returns contents of myPixelX_edit as a double


% --- Executes during object creation, after setting all properties.
function myPixelX_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to myPixelX_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function myPixelY_edit_Callback(hObject, eventdata, handles)
% hObject    handle to myPixelY_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of myPixelY_edit as text
%        str2double(get(hObject,'String')) returns contents of myPixelY_edit as a double


% --- Executes during object creation, after setting all properties.
function myPixelY_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to myPixelY_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function myPixelZ_edit_Callback(hObject, eventdata, handles)
% hObject    handle to myPixelZ_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of myPixelZ_edit as text
%        str2double(get(hObject,'String')) returns contents of myPixelZ_edit as a double


% --- Executes during object creation, after setting all properties.
function myPixelZ_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to myPixelZ_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function myPixelI1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to myPixelI1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of myPixelI1_edit as text
%        str2double(get(hObject,'String')) returns contents of myPixelI1_edit as a double


% --- Executes during object creation, after setting all properties.
function myPixelI1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to myPixelI1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selectLocation_button.
function selectLocation_button_Callback(hObject, eventdata, handles)
% hObject    handle to selectLocation_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(get(handles.projection_bgroup,'SelectedObject'),'Tag'),'projXZ_radio')   %if XZ    
    [x,z] = ginput(1);
    x = round(x); z = round(z);
    set(handles.myPixelX_edit,'String',num2str(x));
    set(handles.myPixelZ_edit,'String',num2str(z));
    set(handles.myPixelY_edit,'String','');
    
    channel1Image = getappdata(handles.MainPanel,'curImage_ch1');
    channel2Image = getappdata(handles.MainPanel,'curImage_ch2');

    set(handles.myPixelI1_edit,'String',num2str(channel1Image(z,x)));
    if ~isempty(channel2Image)
        set(handles.myPixelI2_edit,'String',num2str(channel2Image(z,x)));
    end
else    
    [x,y] = ginput(1);
    x = round(x); y = round(y);
    set(handles.myPixelX_edit,'String',num2str(x));
    set(handles.myPixelY_edit,'String',num2str(y));
    set(handles.myPixelZ_edit,'String','');
    
    channel1Image = getappdata(handles.MainPanel,'curImage_ch1');
    channel2Image = getappdata(handles.MainPanel,'curImage_ch2');

    set(handles.myPixelI1_edit,'String',num2str(channel1Image(y,x)));
    if ~isempty(channel2Image)
        set(handles.myPixelI2_edit,'String',num2str(channel2Image(y,x)));
    end
end


% --- Executes on selection change in ch1Color_popup.
function ch1Color_popup_Callback(hObject, eventdata, handles)
% hObject    handle to ch1Color_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ch1Color_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ch1Color_popup
plotImage(handles);

% --- Executes during object creation, after setting all properties.
function ch1Color_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1Color_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ch2Color_popup.
function ch2Color_popup_Callback(hObject, eventdata, handles)
% hObject    handle to ch2Color_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ch2Color_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ch2Color_popup
plotImage(handles);

% --- Executes during object creation, after setting all properties.
function ch2Color_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch2Color_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in applySubtraction_button.
function applySubtraction_button_Callback(hObject, eventdata, handles)
% hObject    handle to applySubtraction_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function subtractionScaleFactor_edit_Callback(hObject, eventdata, handles)
% hObject    handle to subtractionScaleFactor_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subtractionScaleFactor_edit as text
%        str2double(get(hObject,'String')) returns contents of subtractionScaleFactor_edit as a double


% --- Executes during object creation, after setting all properties.
function subtractionScaleFactor_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subtractionScaleFactor_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in resetSubtraction_button.
function resetSubtraction_button_Callback(hObject, eventdata, handles)
% hObject    handle to resetSubtraction_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function plotOnSelectionChange(hObject, eventdata)
handles = guidata(hObject);
plotImage(handles);

% --------------------------------------------------------------------
function calculateProjection(hObject, eventdata)
% hObject    handle to projection_bgroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
handles = guidata(hObject);
switch get(get(handles.projection_bgroup,'SelectedObject'),'Tag')
    case 'projStack_radio'
        set(handles.projectionMeanVsMax_bgroup,'visible','off');
        set(handles.zFrame_slider,'visible','on'); set(handles.zFrame_edit,'visible','on');
        curZ = round(get(handles.zFrame_slider,'Value'));
        Image1Raw = getappdata(handles.MainPanel,'Image1Raw');
        setappdata(handles.MainPanel,'curImage_ch1',Image1Raw(:,:,curZ));
        if getappdata(handles.MainPanel,'NumChannels') > 1
            Image2Raw = getappdata(handles.MainPanel,'Image2Raw');
            setappdata(handles.MainPanel,'curImage_ch2',Image2Raw(:,:,curZ));
        end
    case 'projXY_radio'
        set(handles.projectionMeanVsMax_bgroup,'visible','on');
        set(handles.zFrame_slider,'visible','off'); set(handles.zFrame_edit,'visible','off');
        Image1Raw = getappdata(handles.MainPanel,'Image1Raw');
        Zrange = str2num(get(handles.zRange1_edit,'String'));
        Image1Raw = Image1Raw(:,:,Zrange(1):Zrange(2));
        switch get(get(handles.projectionMeanVsMax_bgroup,'SelectedObject'),'Tag')            
            case 'projMax_radio'
                setappdata(handles.MainPanel,'curImage_ch1',max(Image1Raw,[],3));
            case 'projMean_radio'
                setappdata(handles.MainPanel,'curImage_ch1',mean(Image1Raw,3));
        end
        if getappdata(handles.MainPanel,'NumChannels') > 1
            Image2Raw = getappdata(handles.MainPanel,'Image2Raw');
            Zrange = str2num(get(handles.zRange2_edit,'String'));
            Image2Raw = Image2Raw(:,:,Zrange(1):Zrange(2));
            switch get(get(handles.projectionMeanVsMax_bgroup,'SelectedObject'),'Tag')
                case 'projMax_radio'
                    setappdata(handles.MainPanel,'curImage_ch2',max(Image2Raw,[],3));
                case 'projMean_radio'
                    setappdata(handles.MainPanel,'curImage_ch2',mean(Image2Raw,3));
            end
        end
    case 'projXZ_radio'
        set(handles.projectionMeanVsMax_bgroup,'visible','on');
        set(handles.zFrame_slider,'visible','off'); set(handles.zFrame_edit,'visible','off');
        Image1Raw = getappdata(handles.MainPanel,'Image1Raw');
        Zrange = str2num(get(handles.zRange1_edit,'String'));
        Image1Raw = Image1Raw(:,:,Zrange(1):Zrange(2));
        switch get(get(handles.projectionMeanVsMax_bgroup,'SelectedObject'),'Tag')
            case 'projMax_radio'
                setappdata(handles.MainPanel,'curImage_ch1',squeeze(max(Image1Raw,[],2))');
            case 'projMean_radio'
                setappdata(handles.MainPanel,'curImage_ch1',squeeze(mean(Image1Raw,2))');
        end
        if getappdata(handles.MainPanel,'NumChannels') > 1
            Image2Raw = getappdata(handles.MainPanel,'Image2Raw');
            Zrange = str2num(get(handles.zRange2_edit,'String'));
            Image2Raw = Image2Raw(:,:,Zrange(1):Zrange(2));
            switch get(get(handles.projectionMeanVsMax_bgroup,'SelectedObject'),'Tag')
                case 'projMax_radio'
                    setappdata(handles.MainPanel,'curImage_ch2',squeeze(max(Image2Raw,[],2))');
                case 'projMean_radio'
                    setappdata(handles.MainPanel,'curImage_ch2',squeeze(mean(Image2Raw,2))');
            end
        end
end
plotImage(handles);


% --------------------------------------------------------------------
function loadConfocalImage_item_Callback(hObject, eventdata, handles)
% hObject    handle to loadConfocalImage_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.RAW_IMAGE_FOLDER)
    DefaultFolder = handles.RAW_IMAGE_FOLDER;
else
    DefaultFolder = '~';
end

[FileName,PathName,FilterIndex] = uigetfile('*.ics','Select Image File',DefaultFolder);
if FilterIndex %if something was selected
    [pathstr, name, ext] = fileparts([PathName FileName]);
    handles.ICSFileName = [pathstr '/' name ext];
    handles.IDSFileName = [pathstr '/' name '.ids'];
    set(handles.MainPanel,'Name','ImageAnalysisGUI (loading)');
    drawnow;
    readICSFile(handles);
    readIDSFIle(handles);
    set(handles.MainPanel,'Name','ImageAnalysisGUI');
    drawnow;
    Image1Raw = getappdata(handles.MainPanel,'Image1Raw');
    Image2Raw = getappdata(handles.MainPanel,'Image2Raw');
    setappdata(handles.MainPanel,'curImage_ch1', Image1Raw(:,:,1));
    %set Zrange
    Zframes = getappdata(handles.MainPanel,'Zframes');
    set(handles.zRange1_edit,'String',num2str([1,Zframes]));    
    if getappdata(handles.MainPanel,'NumChannels') > 1
        setappdata(handles.MainPanel,'curImage_ch2', Image2Raw(:,:,1)); 
        %set Zrange
        set(handles.zRange2_edit,'String',num2str([1,Zframes])); 
    end
    %plot
    plotImage(handles);
end

% --------------------------------------------------------------------
function quit_item_Callback(hObject, eventdata, handles)
% hObject    handle to quit_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.MainPanel);

function readICSFile(handles)
% get nchannels, x, y, z size of zstack
icsfid = fopen(handles.ICSFileName, 'r');
TextLine = fgetl(icsfid);
TargetTextLine1 = 'layout';
TargetTextLine2 = 'sizes';
while (isempty(strfind(TextLine, TargetTextLine1)) | isempty(strfind(TextLine, TargetTextLine2)))
    TextLine = fgetl(icsfid);
    if (TextLine == -1)
        break;
    end
end
ImageParameters = sscanf(TextLine(14:length(TextLine)), '%f %f %f %f %f');
setappdata(handles.MainPanel,'NumChannels',ImageParameters(2));
setappdata(handles.MainPanel,'Xpixels',ImageParameters(3));
setappdata(handles.MainPanel,'Ypixels',ImageParameters(4));
setappdata(handles.MainPanel,'Zframes',ImageParameters(5));
fclose(icsfid);
% get x,y,z,scale (microns/pixel)
icsfid = fopen(handles.ICSFileName, 'r');
TextLine = fgetl(icsfid);
TargetTextLine1 = 'parameter';
TargetTextLine2 = 'scale';
while (isempty(strfind(TextLine, TargetTextLine1)) | isempty(strfind(TextLine, TargetTextLine2)))
    TextLine = fgetl(icsfid);
    if (TextLine == -1)
        break;
    end
end
ImageParameters = sscanf(TextLine(17:length(TextLine)), '%f %f %f %f %f');
setappdata(handles.MainPanel,'Xscale',ImageParameters(3));
setappdata(handles.MainPanel,'Yscale',ImageParameters(4));
setappdata(handles.MainPanel,'Zscale',abs(ImageParameters(5)));

%initialize stuff
Zframes = getappdata(handles.MainPanel,'Zframes');
set(handles.zFrame_slider,'Min',1);
set(handles.zFrame_slider,'Max',Zframes);
set(handles.zFrame_slider,'SliderStep',[1 1]./Zframes);
set(handles.zFrame_slider,'Value',1);
set(handles.zFrame_edit,'String','1');

fclose(icsfid);

function readICSFile_Andor(handles)
% get nchannels, x, y, z size of zstack
icsfid = fopen(handles.ICSFileName, 'r');
TextLine = fgetl(icsfid);
TargetTextLine1 = 'layout';
TargetTextLine2 = 'sizes';
while (isempty(strfind(TextLine, TargetTextLine1)) | isempty(strfind(TextLine, TargetTextLine2)))
    TextLine = fgetl(icsfid);
    if (TextLine == -1)
        break;
    end
end
ImageParameters = sscanf(TextLine(14:length(TextLine)), '%f %f %f %f %f');
setappdata(handles.MainPanel,'NumChannels',1);
setappdata(handles.MainPanel,'Xpixels',ImageParameters(2));
setappdata(handles.MainPanel,'Ypixels',ImageParameters(3));
setappdata(handles.MainPanel,'Zframes',ImageParameters(4));
fclose(icsfid);
% get x,y,z,scale (microns/pixel)
setappdata(handles.MainPanel,'Xscale',1);
setappdata(handles.MainPanel,'Yscale',1);
setappdata(handles.MainPanel,'Zscale',1);

%initialize stuff
Zframes = getappdata(handles.MainPanel,'Zframes');
set(handles.zFrame_slider,'Min',1);
set(handles.zFrame_slider,'Max',Zframes);
set(handles.zFrame_slider,'SliderStep',[1 1]./Zframes);
set(handles.zFrame_slider,'Value',1);
set(handles.zFrame_edit,'String','1');


function readIDSFIle(handles)
%channel 1
fid = fopen(handles.IDSFileName, 'r');
fseek(fid, 4, 'bof');
Xpixels = getappdata(handles.MainPanel,'Xpixels');
Ypixels = getappdata(handles.MainPanel,'Ypixels');
Zframes = getappdata(handles.MainPanel,'Zframes');
im = zeros(Xpixels,Ypixels,Zframes);
for r = 1:Zframes
    im(:, :, r) = fread(fid, [Xpixels Ypixels], 'uint16', 2*(getappdata(handles.MainPanel,'NumChannels')-1), 'l');
end
im = im./max(max(max(im)));
setappdata(handles.MainPanel,'Image1Raw',im);
fclose(fid);

%channel 2
if getappdata(handles.MainPanel,'NumChannels') > 1
    fid = fopen(handles.IDSFileName, 'r');
    fseek(fid, 2, 'bof');
    im2 = zeros(Xpixels,Ypixels,Zframes);
    for r = 1:Zframes
        im2(:, :, r) = fread(fid, [Xpixels Ypixels], 'uint16', 2*(getappdata(handles.MainPanel,'NumChannels')-1), 'l');
    end
    im2 = im2./max(max(max(im2)));
    setappdata(handles.MainPanel,'Image2Raw',im2);
    fclose(fid);
end

function readIDSFIle_Andor(handles)
%channel 1
fid = fopen(handles.IDSFileName, 'r');
%fseek(fid, 4, 'bof');
Xpixels = getappdata(handles.MainPanel,'Xpixels');
Ypixels = getappdata(handles.MainPanel,'Ypixels');
Zframes = getappdata(handles.MainPanel,'Zframes');
im = zeros(Xpixels,Ypixels,Zframes);
for r = 1:Zframes
    im(:, :, r) = fread(fid, [Xpixels Ypixels], 'uint16', 0, 'l');
end
im = im./max(max(max(im)));
setappdata(handles.MainPanel,'Image1Raw',im);
fclose(fid);


function processCurImages(handles)
gamma1 = str2num(get(handles.gamma1_edit,'String'));
lowCut1 = str2num(get(handles.lowCut1_edit,'String'));
highCut1 = str2num(get(handles.highCut1_edit,'String'));
gamma2 = str2num(get(handles.gamma2_edit,'String'));
lowCut2 = str2num(get(handles.lowCut2_edit,'String'));
highCut2 = str2num(get(handles.highCut2_edit,'String'));

channel1Image = getappdata(handles.MainPanel,'curImage_ch1');
channel2Image = getappdata(handles.MainPanel,'curImage_ch2');

channel1Image = channel1Image./max(max(channel1Image)); %normalize
channel1Image = channel1Image.^gamma1;
channel1Image = channel1Image./max(max(channel1Image)); %normalize
[r,c] = size(channel1Image);
channel1ImageFlat = reshape(channel1Image,1,r*c);
sdev = std(channel1ImageFlat);
m = mean(channel1ImageFlat);
lowVal = m - lowCut1*sdev;
highVal = m + highCut1*sdev;
channel1Image(channel1Image < lowVal) = lowVal;
channel1Image(channel1Image > highVal) = highVal;
channel1Image = channel1Image - lowVal;
channel1Image = channel1Image./max(max(channel1Image)); %normalize

%channel1Image(channel1Image<lowCut1) = 0;
%channel1Image(channel1Image>highCut1) = 1;
setappdata(handles.MainPanel,'curImage_ch1',channel1Image);

channel2Image = channel2Image./max(max(channel2Image)); %normalize
channel2Image = channel2Image.^gamma2;
channel2Image = channel2Image./max(max(channel2Image)); %normalize
channel2Image(channel2Image<lowCut2) = 0;
channel2Image(channel2Image>highCut2) = 1;
[r,c] = size(channel2Image);
channel2ImageFlat = reshape(channel2Image,1,r*c);
sdev = std(channel2ImageFlat);
m = mean(channel2ImageFlat);
lowVal = m - lowCut2*sdev;
highVal = m + highCut2*sdev;
channel2Image(channel2Image < lowVal) = lowVal;
channel2Image(channel2Image > highVal) = highVal;
channel2Image = channel2Image - lowVal;
channel2Image = channel2Image./max(max(channel2Image)); %normalize

setappdata(handles.MainPanel,'curImage_ch2',channel2Image);

function plotImage(handles)
setappdata(handles.MainPanel,'drawInProgress',1);
processCurImages(handles);

%check which type of view
switch get(get(handles.combinedVsSeparate_bgroup,'SelectedObject'),'Tag')
    case 'viewCombined_radio'
        set(handles.viewChannel_bgroup,'visible','off');
        %clear separate axes
        cla(handles.imageCh1_axis); set(handles.imageCh1_axis,'visible','off');
        cla(handles.imageCh2_axis); set(handles.imageCh2_axis,'visible','off');
        
        set(handles.imageCombined_axis,'visible','on');
        axes(handles.imageCombined_axis);
        %colormap('gray');
        channel1Image = getappdata(handles.MainPanel,'curImage_ch1');
        channel2Image = getappdata(handles.MainPanel,'curImage_ch2');
        [r,c] = size(channel1Image);
        %get color for each channel
        ch1Color = get(handles.ch1Color_popup,'Value');
        if ch1Color > 3, ch1Color = 1; end;
        ch2Color = get(handles.ch2Color_popup,'Value');
        if ch2Color > 3, ch1Color = 2; end;
        %make composite image
        curImage = zeros(r,c,3);
        curImage(:,:,ch1Color) = channel1Image;
        curImage(:,:,ch2Color) = channel2Image;
        %plot it
        image(curImage);
        setappdata(handles.MainPanel,'onscreenImage',curImage);
        setappdata(handles.MainPanel,'onscreenImage1',[]);
        setappdata(handles.MainPanel,'onscreenImage2',[]);
        %set axis lables
        switch get(get(handles.units_bgroup,'SelectedObject'),'Tag')
            case 'pixelUnits_radio'
                tickLocations = get(handles.imageCombined_axis,'Xtick');
                set(handles.imageCombined_axis,'XtickLabel',num2str(tickLocations'));
                tickLocations = get(handles.imageCombined_axis,'Ytick');
                set(handles.imageCombined_axis,'YtickLabel',num2str(tickLocations'));
            case 'micronsUnits_radio'
                tickLocations = get(handles.imageCombined_axis,'Xtick');
                tickLabels = round(tickLocations.*getappdata(handles.MainPanel,'Xscale'));
                set(handles.imageCombined_axis,'XtickLabel',num2str(tickLabels'));
                switch get(get(handles.projection_bgroup,'SelectedObject'),'Tag')
                    case 'projXZ_radio'
                        tickLocations = get(handles.imageCombined_axis,'Ytick');
                        tickLabels = round(tickLocations.*getappdata(handles.MainPanel,'Zscale'));
                        set(handles.imageCombined_axis,'YtickLabel',num2str(tickLabels'));
                    otherwise
                        tickLocations = get(handles.imageCombined_axis,'Ytick');
                        tickLabels = round(tickLocations.*getappdata(handles.MainPanel,'Yscale'));
                        set(handles.imageCombined_axis,'YtickLabel',num2str(tickLabels'));
                end
        end
    case 'viewSingle_radio'
        set(handles.viewChannel_bgroup,'visible','on');
        %clear separate axes
        cla(handles.imageCh1_axis); set(handles.imageCh1_axis,'visible','off');
        cla(handles.imageCh2_axis); set(handles.imageCh2_axis,'visible','off');
        
        switch get(get(handles.viewChannel_bgroup,'SelectedObject'),'Tag')
            case 'viewCh1_radio'
                curImage = round(getappdata(handles.MainPanel,'curImage_ch1')*255);
                curColor = get(handles.ch1Color_popup,'Value');
            case 'viewCh2_radio'
                curImage = round(getappdata(handles.MainPanel,'curImage_ch2')*255);
                curColor = get(handles.ch2Color_popup,'Value');
        end
        set(handles.imageCombined_axis,'visible','on');
        axes(handles.imageCombined_axis);
        [r,c] = size(curImage);
        
        if curColor <= 3 %r,g,b
            colorImage = zeros(r,c,3);
            colorImage(:,:,curColor) = curImage./255;
            image(colorImage);
            setappdata(handles.MainPanel,'onscreenImage',colorImage);
        elseif curColor == 4 %greyscale
            colormap('gray');
            image(curImage);
            setappdata(handles.MainPanel,'onscreenImage',curImage);
        elseif  curColor == 5 %pseudocolor
            colormap('jet');
            image(curImage);
            setappdata(handles.MainPanel,'onscreenImage',curImage);
        end
        setappdata(handles.MainPanel,'onscreenImage1',[]);
        setappdata(handles.MainPanel,'onscreenImage2',[]);
        %set axis lables
        switch get(get(handles.units_bgroup,'SelectedObject'),'Tag')
            case 'pixelUnits_radio'
                tickLocations = get(handles.imageCombined_axis,'Xtick');
                set(handles.imageCombined_axis,'XtickLabel',num2str(tickLocations'));
                tickLocations = get(handles.imageCombined_axis,'Ytick');
                set(handles.imageCombined_axis,'YtickLabel',num2str(tickLocations'));
            case 'micronsUnits_radio'
                tickLocations = get(handles.imageCombined_axis,'Xtick');
                tickLabels = round(tickLocations.*getappdata(handles.MainPanel,'Xscale'));
                set(handles.imageCombined_axis,'XtickLabel',num2str(tickLabels'));
                switch get(get(handles.projection_bgroup,'SelectedObject'),'Tag')
                    case 'projXZ_radio'
                        tickLocations = get(handles.imageCombined_axis,'Ytick');
                        tickLabels = round(tickLocations.*getappdata(handles.MainPanel,'Zscale'));
                        set(handles.imageCombined_axis,'YtickLabel',num2str(tickLabels'));
                    otherwise
                        tickLocations = get(handles.imageCombined_axis,'Ytick');
                        tickLabels = round(tickLocations.*getappdata(handles.MainPanel,'Yscale'));
                        set(handles.imageCombined_axis,'YtickLabel',num2str(tickLabels'));
                end
        end        
    case 'viewSeparate_radio'
        set(handles.viewChannel_bgroup,'visible','off');
        %clear combined axis
        cla(handles.imageCombined_axis); set(handles.imageCombined_axis,'visible','off');
        
        %channel 1
        axes(handles.imageCh1_axis);
        set(handles.imageCh1_axis,'visible','on');
        channel1Image = round(getappdata(handles.MainPanel,'curImage_ch1')*255);
        [r,c] = size(channel1Image);
        curColor = get(handles.ch1Color_popup,'Value');
        if curColor <= 3 %r,g,b
            colorImage = zeros(r,c,3);
            colorImage(:,:,curColor) = channel1Image./255;
            image(colorImage);
            setappdata(handles.MainPanel,'onscreenImage1',colorImage);
        elseif curColor == 4 %greyscale
            colormap('gray');
            image(channel1Image);
            setappdata(handles.MainPanel,'onscreenImage1',channel1Image);
        elseif  curColor == 5 %pseudocolor
            colormap('hsv');
            image(channel1Image);
            setappdata(handles.MainPanel,'onscreenImage1',channel1Image);
        end
        %channel 2
        if getappdata(handles.MainPanel,'NumChannels') > 1
            axes(handles.imageCh2_axis);
            set(handles.imageCh1_axis,'visible','on');
            channel2Image = round(getappdata(handles.MainPanel,'curImage_ch2')*255);
            
            curColor = get(handles.ch2Color_popup,'Value');
            if curColor <= 3 %r,g,b
                colorImage = zeros(r,c,3);
                colorImage(:,:,curColor) = channel2Image./255;
                image(colorImage);
                setappdata(handles.MainPanel,'onscreenImage2',colorImage);
            elseif curColor == 4 %greyscale
                colormap('gray');
                image(channel2Image);
                setappdata(handles.MainPanel,'onscreenImage2',channel2Image);
            elseif  curColor == 5 %pseudocolor
                colormap('hsv');
                image(channel2Image);
                setappdata(handles.MainPanel,'onscreenImage2',channel2Image);
            end
        end
        setappdata(handles.MainPanel,'onscreenImage',[]);
        %set axis lables
        switch get(get(handles.units_bgroup,'SelectedObject'),'Tag')
            case 'pixelUnits_radio'
                tickLocations = get(handles.imageCh1_axis,'Xtick');
                set(handles.imageCh1_axis,'XtickLabel',num2str(tickLocations'));
                tickLocations = get(handles.imageCh1_axis,'Ytick');
                set(handles.imageCh1_axis,'YtickLabel',num2str(tickLocations'));
                tickLocations = get(handles.imageCh2_axis,'Xtick');
                set(handles.imageCh2_axis,'XtickLabel',num2str(tickLocations'));
                tickLocations = get(handles.imageCh2_axis,'Ytick');
                set(handles.imageCh2_axis,'YtickLabel',num2str(tickLocations'));
            case 'micronsUnits_radio'
                tickLocations = get(handles.imageCh1_axis,'Xtick');
                tickLabels = round(tickLocations.*getappdata(handles.MainPanel,'Xscale'));
                set(handles.imageCh1_axis,'XtickLabel',num2str(tickLabels'));
                tickLocations = get(handles.imageCh2_axis,'Xtick');
                tickLabels = round(tickLocations.*getappdata(handles.MainPanel,'Xscale'));
                set(handles.imageCh2_axis,'XtickLabel',num2str(tickLabels'));
                switch get(get(handles.projection_bgroup,'SelectedObject'),'Tag')
                    case 'projXZ_radio'
                        tickLocations = get(handles.imageCh1_axis,'Ytick');
                        tickLabels = round(tickLocations.*getappdata(handles.MainPanel,'Zscale'));
                        set(handles.imageCh1_axis,'YtickLabel',num2str(tickLabels'));
                        tickLocations = get(handles.imageCh2_axis,'Ytick');
                        tickLabels = round(tickLocations.*getappdata(handles.MainPanel,'Zscale'));
                        set(handles.imageCh2_axis,'YtickLabel',num2str(tickLabels'));
                    otherwise
                        tickLocations = get(handles.imageCh1_axis,'Ytick');
                        tickLabels = round(tickLocations.*getappdata(handles.MainPanel,'Yscale'));
                        set(handles.imageCh1_axis,'YtickLabel',num2str(tickLabels'));
                        tickLocations = get(handles.imageCh2_axis,'Ytick');
                        tickLabels = round(tickLocations.*getappdata(handles.MainPanel,'Yscale'));
                        set(handles.imageCh2_axis,'YtickLabel',num2str(tickLabels'));
                end
        end 
end

setappdata(handles.MainPanel,'drawInProgress',0);


% --------------------------------------------------------------------
function projection_bgroup_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to projection_bgroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function zRange2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to zRange2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zRange2_edit as text
%        str2double(get(hObject,'String')) returns contents of zRange2_edit as a double
calculateProjection(hObject,eventdata);

% --- Executes during object creation, after setting all properties.
function zRange2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zRange2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gamma2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to gamma2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma2_edit as text
%        str2double(get(hObject,'String')) returns contents of gamma2_edit as a double
calculateProjection(hObject,eventdata);

% --- Executes during object creation, after setting all properties.
function gamma2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lowCut2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to lowCut2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowCut2_edit as text
%        str2double(get(hObject,'String')) returns contents of lowCut2_edit as a double
calculateProjection(hObject,eventdata);

% --- Executes during object creation, after setting all properties.
function lowCut2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowCut2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function highCut2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to highCut2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of highCut2_edit as text
%        str2double(get(hObject,'String')) returns contents of highCut2_edit as a double
calculateProjection(hObject,eventdata);

% --- Executes during object creation, after setting all properties.
function highCut2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to highCut2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function myPixelI2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to myPixelI2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of myPixelI2_edit as text
%        str2double(get(hObject,'String')) returns contents of myPixelI2_edit as a double


% --- Executes during object creation, after setting all properties.
function myPixelI2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to myPixelI2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function loadAndorImage_item_Callback(hObject, eventdata, handles)
% hObject    handle to loadAndorImage_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.RAW_IMAGE_FOLDER)
    DefaultFolder = handles.RAW_IMAGE_FOLDER;
else
    DefaultFolder = '~';
end

[FileName,PathName,FilterIndex] = uigetfile('*.ics','Select Image File',DefaultFolder);
if FilterIndex %if something was selected
    [pathstr, name, ext] = fileparts([PathName FileName]);
    handles.ICSFileName = [pathstr '/' name ext];
    handles.IDSFileName = [pathstr '/' name '.ids'];
    set(handles.MainPanel,'Name','ImageAnalysisGUI (loading)');
    drawnow;
    readICSFile_Andor(handles);
    readIDSFIle_Andor(handles);
    set(handles.MainPanel,'Name','ImageAnalysisGUI');
    drawnow;
    Image1Raw = getappdata(handles.MainPanel,'Image1Raw');
    setappdata(handles.MainPanel,'curImage_ch1', Image1Raw(:,:,1));
    %set Zrange
    Zframes = getappdata(handles.MainPanel,'Zframes');
    set(handles.zRange1_edit,'String',num2str([1,Zframes]));    
     %plot
    plotImage(handles);
end