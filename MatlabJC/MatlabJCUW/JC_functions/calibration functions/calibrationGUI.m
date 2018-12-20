function varargout = calibrationGUI(varargin)
% CALIBRATIONGUI M-file for calibrationGUI.fig
%      CALIBRATIONGUI, by itself, creates a new CALIBRATIONGUI or raises the existing
%      singleton*.
%
%      H = CALIBRATIONGUI returns the handle to a new CALIBRATIONGUI or the handle to
%      the existing singleton*.
%
%      CALIBRATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALIBRATIONGUI.M with the given input arguments.
%
%      CALIBRATIONGUI('Property','Value',...) creates a new CALIBRATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calibrationGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calibrationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calibrationGUI

% Last Modified by GUIDE v2.5 01-Mar-2010 08:50:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @calibrationGUI_OpeningFcn, ...
    'gui_OutputFcn',  @calibrationGUI_OutputFcn, ...
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


% --- Executes just before calibrationGUI is made visible.
function calibrationGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calibrationGUI (see VARARGIN)

% Choose default command line output for calibrationGUI

[s,r] = system('df | grep /Volumes/Groups | grep ^afp | awk ''{print $6}''') ; % this line finds the most recent "/Volumes/Group"

handles.output = hObject;

handles.calibrationRootDir = [r(1:end-1),'/Lab Staff/Calibration/'] ;

handles.SpectraRange = [370:720]; %in nm
handles.hc_over_lambda = 6.63E-34*3E8 ./ (handles.SpectraRange .* 1E-9);

%display colors
handles.colors{1} = [1 0 0]; %red
handles.colors{2} = [0 0.5 0]; %green
handles.colors{3} = [0 0 1]; %blue

%load default stuff
handles = loadDefaults(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes calibrationGUI wait for user response (see UIRESUME)
% uiwait(handles.background);


% --- Outputs from this function are returned to the command line.
function varargout = calibrationGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function File_menu_Callback(hObject, eventdata, handles)
% hObject    handle to File_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Load_LED_Spectra_item_Callback(hObject, eventdata, handles)
% hObject    handle to Load_LED_Spectra_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, path] = uigetfile('*.mat','LED Spectra File');
handles.LEDSPectraFile = [path fname];
handles = loadLEDSpectra(handles);
guidata(hObject, handles);

function handles = loadLEDSpectra(handles)
load(handles.LEDSPectraFile,'Spectra');
handles.LEDSpectra_X = Spectra(:,1);
handles.LEDSpectra_red = Spectra(:,2);
handles.LEDSpectra_green = Spectra(:,3);
handles.LEDSpectra_blue = Spectra(:,4);

plotLEDSpectra(handles);


% --------------------------------------------------------------------
function Load_PR_Spectra_item_Callback(hObject, eventdata, handles)
% hObject    handle to Load_PR_Spectra_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, path] = uigetfile('*.mat','Photoreceptor Spectra File');
handles.PRSpectraFile = [path fname];
handles = loadPRSpectra(handles);
guidata(hObject, handles);


function handles = loadPRSpectra(handles)
load(handles.PRSpectraFile ,'Spectra');
handles.PRSpectra_X = Spectra(:,1);
handles.PRSpectra_rod = Spectra(:,2);
handles.PRSpectra_scone = Spectra(:,3);
handles.PRSpectra_mcone = Spectra(:,4);

if size(Spectra,2) == 5
    handles.PRSpectra_lcone = Spectra(:,5);
else
    handles.PRSpectra_lcone = ones(size(Spectra(:,1))).*NaN;
end

plotPRSpectra(handles);

function plotPRSpectra(handles)

h = plot(handles.PR_Spectra_axis,handles.PRSpectra_X,handles.PRSpectra_rod,'k',...
    handles.PRSpectra_X,handles.PRSpectra_scone,'b',...
    handles.PRSpectra_X,handles.PRSpectra_mcone,'g',...
    handles.PRSpectra_X,handles.PRSpectra_lcone,'r');

axis(handles.PR_Spectra_axis,[370 720 0 1.2]);
xlabel(handles.PR_Spectra_axis,'Wavelength (nm)');
ylabel(handles.PR_Spectra_axis,'Relative Amplitude');
legend(h,{'rod','s-cone','m-cone','l-cone'});


function ret = mystr2num(s)
if isempty(s), ret = nan;
else
    ret = str2num(s);
end


function handles = calibrateLEDSpectra(hObject,handles)

red1Value = mystr2num(get(handles.red1_input,'String'));
red2Value = mystr2num(get(handles.red2_input,'String'));
red3Value = mystr2num(get(handles.red3_input,'String'));

green1Value = mystr2num(get(handles.green1_input,'String'));
green2Value = mystr2num(get(handles.green2_input,'String'));
green3Value = mystr2num(get(handles.green3_input,'String'));

blue1Value = mystr2num(get(handles.blue1_input,'String'));
blue2Value = mystr2num(get(handles.blue2_input,'String'));
blue3Value = mystr2num(get(handles.blue3_input,'String'));
blue4Value = mystr2num(get(handles.blue4_input,'String'));

%Note: we want to put light meter spectrum here

%1e-9 to convert power values off the meter into watts
handles.LEDSpectra_red_calibrated(:,1) = (handles.LEDSpectra_red.* handles.hc_over_lambda').*(red1Value*1E-9/sum(handles.LEDSpectra_red.* handles.hc_over_lambda'))./handles.hc_over_lambda' ;
handles.LEDSpectra_red_calibrated(:,2) = (handles.LEDSpectra_red.* handles.hc_over_lambda').*(red2Value*1E-9/sum(handles.LEDSpectra_red.* handles.hc_over_lambda'))./handles.hc_over_lambda' ;
handles.LEDSpectra_red_calibrated(:,3) = (handles.LEDSpectra_red.* handles.hc_over_lambda').*(red3Value*1E-9/sum(handles.LEDSpectra_red.* handles.hc_over_lambda'))./handles.hc_over_lambda' ;

handles.LEDSpectra_green_calibrated(:,1) = (handles.LEDSpectra_green.* handles.hc_over_lambda').*(green1Value*1E-9/sum(handles.LEDSpectra_green.* handles.hc_over_lambda'))./handles.hc_over_lambda' ;
handles.LEDSpectra_green_calibrated(:,2) = (handles.LEDSpectra_green.* handles.hc_over_lambda').*(green2Value*1E-9/sum(handles.LEDSpectra_green.* handles.hc_over_lambda'))./handles.hc_over_lambda' ;
handles.LEDSpectra_green_calibrated(:,3) = (handles.LEDSpectra_green.* handles.hc_over_lambda').*(green3Value*1E-9/sum(handles.LEDSpectra_green.* handles.hc_over_lambda'))./handles.hc_over_lambda' ;

handles.LEDSpectra_blue_calibrated(:,1) = (handles.LEDSpectra_blue.* handles.hc_over_lambda').*(blue1Value*1E-9/sum(handles.LEDSpectra_blue.* handles.hc_over_lambda'))./handles.hc_over_lambda' ;
handles.LEDSpectra_blue_calibrated(:,2) = (handles.LEDSpectra_blue.* handles.hc_over_lambda').*(blue2Value*1E-9/sum(handles.LEDSpectra_blue.* handles.hc_over_lambda'))./handles.hc_over_lambda' ;
handles.LEDSpectra_blue_calibrated(:,3) = (handles.LEDSpectra_blue.* handles.hc_over_lambda').*(blue3Value*1E-9/sum(handles.LEDSpectra_blue.* handles.hc_over_lambda'))./handles.hc_over_lambda' ;
handles.LEDSpectra_blue_calibrated(:,4) = (handles.LEDSpectra_blue.* handles.hc_over_lambda').*(blue4Value*1E-9/sum(handles.LEDSpectra_blue.* handles.hc_over_lambda'))./handles.hc_over_lambda' ;

%Filter paper attenuation
if get(handles.Use_filterPaperFactors_checkbox,'Value')
    handles.FilterPaperAttenuation.red = mystr2num(get(handles.filter_attenuation_red_edit,'String'));
    handles.FilterPaperAttenuation.green = mystr2num(get(handles.filter_attenuation_green_edit,'String'));
    handles.FilterPaperAttenuation.blue = mystr2num(get(handles.filter_attenuation_blue_edit,'String'));
    
    handles.LEDSpectra_red_calibrated = handles.LEDSpectra_red_calibrated ./ handles.FilterPaperAttenuation.red;
    handles.LEDSpectra_green_calibrated = handles.LEDSpectra_green_calibrated ./ handles.FilterPaperAttenuation.green;
    handles.LEDSpectra_blue_calibrated = handles.LEDSpectra_blue_calibrated ./ handles.FilterPaperAttenuation.blue;
end


function handles = createConversionMatrix(hObject,handles)

%concatenate calibrated spectra [r1 r2 r3 g1 g2 g3 b1 b2 b3 b4]
allLEDSpectra_calibrated_layer = [handles.LEDSpectra_red_calibrated handles.LEDSpectra_green_calibrated handles.LEDSpectra_blue_calibrated];

%divide everything by spot size: Values are now photons / (sec*V*um^2)
allLEDSpectra_calibrated_layer = allLEDSpectra_calibrated_layer ./ handles.spotSize;
allLEDSpectra_calibrated_layerNDF = allLEDSpectra_calibrated_layer.* repmat(handles.NDFSpectra,1,10);

%get collecting area
%collecting area
CollectingArea_rod = mystr2num(get(handles.rod_collecting_area_edit,'String'));
CollectingArea_cone = mystr2num(get(handles.cone_collecting_area_edit,'String'));

%multiply by PR spectra for each PR to get isomerizations, correct for
%collection area as well
isom_rod = handles.PRSpectra_rod'*allLEDSpectra_calibrated_layer .* CollectingArea_rod;
isom_scone = handles.PRSpectra_scone'*allLEDSpectra_calibrated_layer .* CollectingArea_cone;
isom_mcone = handles.PRSpectra_mcone'*allLEDSpectra_calibrated_layer .* CollectingArea_cone;
isom_lcone = handles.PRSpectra_lcone'*allLEDSpectra_calibrated_layer .* CollectingArea_cone;

isom_rodNDF = handles.PRSpectra_rod'*allLEDSpectra_calibrated_layerNDF .* CollectingArea_rod;
isom_sconeNDF = handles.PRSpectra_scone'*allLEDSpectra_calibrated_layerNDF .* CollectingArea_cone;
isom_mconeNDF = handles.PRSpectra_mcone'*allLEDSpectra_calibrated_layerNDF .* CollectingArea_cone;
isom_lconeNDF = handles.PRSpectra_lcone'*allLEDSpectra_calibrated_layerNDF .* CollectingArea_cone;

%multiply by NDFs to make a conversion table (V to isom) for each PR.
%Table will be size NDFs x LED switches*colors (ex. 3x10)
L = length(handles.NDFs);
for i=1:L
    if handles.NDFs(i) == 0
        handles.conversionMatrix_rod(i,:) = isom_rod;
        handles.conversionMatrix_scone(i,:) = isom_scone;
        handles.conversionMatrix_mcone(i,:) = isom_mcone;
        handles.conversionMatrix_lcone(i,:) = isom_lcone; 
    else
        handles.conversionMatrix_rod(i,:) = isom_rodNDF.*10^-handles.NDFs(i);
        handles.conversionMatrix_scone(i,:) = isom_sconeNDF.*10^-handles.NDFs(i);
        handles.conversionMatrix_mcone(i,:) = isom_mconeNDF.*10^-handles.NDFs(i);
        handles.conversionMatrix_lcone(i,:) = isom_lconeNDF.*10^-handles.NDFs(i);
    end
end

function plotLEDSpectra(handles)

plot(handles.LED_Spectra_axis,handles.LEDSpectra_X,handles.LEDSpectra_red,'r',...
    handles.LEDSpectra_X,handles.LEDSpectra_green,'g',...
    handles.LEDSpectra_X,handles.LEDSpectra_blue,'b')

set(handles.LED_Spectra_axis,'xlim',[handles.SpectraRange(1),handles.SpectraRange(end)]);
%axis(handles.LED_Spectra_axis,[370 720 0 1.2]);
xlabel(handles.LED_Spectra_axis,'Wavelength (nm)');
ylabel(handles.LED_Spectra_axis,'Energy (uW/V)');

function red1_input_Callback(hObject, eventdata, handles)
% hObject    handle to red1_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of red1_input as text
%        str2double(get(hObject,'String')) returns contents of red1_input
%        as a double
handles = clearTables(handles);



% --- Executes during object creation, after setting all properties.
function red1_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red1_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function blue2_input_Callback(hObject, eventdata, handles)
% hObject    handle to blue2_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blue2_input as text
%        str2double(get(hObject,'String')) returns contents of blue2_input as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function blue2_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue2_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function blue3_input_Callback(hObject, eventdata, handles)
% hObject    handle to blue3_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blue3_input as text
%        str2double(get(hObject,'String')) returns contents of blue3_input as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function blue3_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue3_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function blue1_input_Callback(hObject, eventdata, handles)
% hObject    handle to blue1_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blue1_input as text
%        str2double(get(hObject,'String')) returns contents of blue1_input as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function blue1_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue1_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function green2_input_Callback(hObject, eventdata, handles)
% hObject    handle to green2_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of green2_input as text
%        str2double(get(hObject,'String')) returns contents of green2_input as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function green2_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green2_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function green3_input_Callback(hObject, eventdata, handles)
% hObject    handle to green3_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of green3_input as text
%        str2double(get(hObject,'String')) returns contents of green3_input as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function green3_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green3_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function green1_input_Callback(hObject, eventdata, handles)
% hObject    handle to green1_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of green1_input as text
%        str2double(get(hObject,'String')) returns contents of green1_input as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function green1_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to green1_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function red2_input_Callback(hObject, eventdata, handles)
% hObject    handle to red2_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of red2_input as text
%        str2double(get(hObject,'String')) returns contents of red2_input as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function red2_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red2_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function red3_input_Callback(hObject, eventdata, handles)
% hObject    handle to red3_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of red3_input as text
%        str2double(get(hObject,'String')) returns contents of red3_input as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function red3_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to red3_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Load_NDF_Spectra_item_Callback(hObject, eventdata, handles)
% hObject    handle to Load_NDF_Spectra_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in colorSelector_menu.
function colorSelector_menu_Callback(hObject, eventdata, handles)
% hObject    handle to colorSelector_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns colorSelector_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colorSelector_menu

%recreate conversion matrix
handles = calibrateLEDSpectra(hObject,handles);
handles = createConversionMatrix(hObject,handles);
guidata(hObject, handles);

D_rod = get(handles.rod_table,'Data');
%D_scone(2:end,:) = D(2:end,:);
D_rod = doConversionCalculation(handles,D_rod,handles.conversionMatrix_rod,0);
set(handles.rod_table,'Data',D_rod);

D_scone = get(handles.scone_table,'Data');
%D_scone(2:end,:) = D(2:end,:);
D_scone = doConversionCalculation(handles,D_scone,handles.conversionMatrix_scone,0);
set(handles.scone_table,'Data',D_scone);

D_mcone = get(handles.mcone_table,'Data');
%D_mcone(2:end,:) = D(2:end,:);
D_mcone = doConversionCalculation(handles,D_mcone,handles.conversionMatrix_mcone,0);
set(handles.mcone_table,'Data',D_mcone);

D_lcone = get(handles.lcone_table,'Data');
%D_lcone(2:end,:) = D(2:end,:);
D_lcone = doConversionCalculation(handles,D_lcone,handles.conversionMatrix_lcone,0);
set(handles.lcone_table,'Data',D_lcone);

set(handles.rod_table,'ForegroundColor',handles.colors{get(hObject,'Value')});
set(handles.scone_table,'ForegroundColor',handles.colors{get(hObject,'Value')});
set(handles.mcone_table,'ForegroundColor',handles.colors{get(hObject,'Value')});
set(handles.lcone_table,'ForegroundColor',handles.colors{get(hObject,'Value')});



% --- Executes during object creation, after setting all properties.
function colorSelector_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colorSelector_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function date_input_Callback(hObject, eventdata, handles)
% hObject    handle to date_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of date_input as text
%        str2double(get(hObject,'String')) returns contents of date_input as a double

date = get(hObject,'String') ;
if length(date)~=8 | ~strcmp(date([3,6]),'//')
    errordlg('date must be entered as mm/dd/yy')
end
    
% --- Executes during object creation, after setting all properties.
function date_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to date_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String',datestr(datenum(date),2));


% --- Executes on button press in Use_filterPaperFactors_checkbox.
function Use_filterPaperFactors_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Use_filterPaperFactors_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Use_filterPaperFactors_checkbox
handles = clearTables(handles);

% --- Executes when entered data in editable cell(s) in Results_table.
function Results_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to Results_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on Results_table and none of its controls.
function Results_table_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Results_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in rig_selector_menu.
function rig_selector_menu_Callback(hObject, eventdata, handles)
% hObject    handle to rig_selector_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns rig_selector_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rig_selector_menu

rigNumber = get(hObject,'Value');
%rigLetters = get(hObject,'String');
rigLetters = {'rig','rigB - two-photon','rigC - suction','rigE - confocal','rigF - old slice'} ;
rigLetter = rigLetters{rigNumber};

handles.currentRigDir = [handles.calibrationRootDir,rigLetter,'/'] ;

%calibration text file
handles.calibrationFileName = [handles.currentRigDir,rigLetter,'calibrationValues'];

%spot size
handles.spotSize = load([handles.currentRigDir,[[rigLetters{rigNumber}(1:4)],'spotSize.txt']]);
set(handles.spotSize_edit,'String',num2str(handles.spotSize));

%LED Spectra
handles.LEDSPectraFile = [handles.currentRigDir,[rigLetters{rigNumber}(1:4),'LEDSpectra.mat']];
handles = loadLEDSpectra(handles);

%clear tables
handles = clearTables(handles);

guidata(hObject, handles);
%keyboard;

function handles = loadDefaults(handles)

%NDFs 
handles.NDFs = [0 2 4];
set(handles.ND_filters_edit,'String',num2str(handles.NDFs));
load([handles.calibrationRootDir,'BlueNDFSpectra']); %all LEDs will use the blue NDF spectra
handles.NDFSpectra = BlueNDFSpectra(:,2) ;


%PR Spectra
handles.PRSpectraFile = [handles.calibrationRootDir, 'photoreceptor_spectra/macaque spectra/macaque_all.mat'];
handles = loadPRSpectra(handles);

%collecting area
load([handles.calibrationRootDir, 'CollectingArea.mat'],'CollectingArea');
handles.CollectingArea = CollectingArea;
set(handles.rod_collecting_area_edit,'String',num2str(handles.CollectingArea.rod));
set(handles.cone_collecting_area_edit,'String',num2str(handles.CollectingArea.cone));

if get(handles.Use_filterPaperFactors_checkbox,'Value')
    load([handles.calibrationRootDir 'FilterPaperAttenuation.mat'],'FilterPaperAttenuation');
    handles.FilterPaperAttenuation = FilterPaperAttenuation;
    
    set(handles.filter_attenuation_red_edit,'String',num2str(handles.FilterPaperAttenuation.red));
    set(handles.filter_attenuation_green_edit,'String',num2str(handles.FilterPaperAttenuation.green));
    set(handles.filter_attenuation_blue_edit,'String',num2str(handles.FilterPaperAttenuation.blue));
end

clearTables(handles);

% --- Executes during object creation, after setting all properties.
function rig_selector_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rig_selector_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = datefixer(handles)
% change date dtring into number
date = get(handles.date_input,'string') ;
handles.dateNum = str2num(date([1:2,4:5,7:8])) ;


% --- Executes on button press in save_calibration_button.
function save_calibration_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_calibration_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get calibration inputs
red1Value = mystr2num(get(handles.red1_input,'String'));
red2Value = mystr2num(get(handles.red2_input,'String'));
red3Value = mystr2num(get(handles.red3_input,'String'));

green1Value = mystr2num(get(handles.green1_input,'String'));
green2Value = mystr2num(get(handles.green2_input,'String'));
green3Value = mystr2num(get(handles.green3_input,'String'));

blue1Value = mystr2num(get(handles.blue1_input,'String'));
blue2Value = mystr2num(get(handles.blue2_input,'String'));
blue3Value = mystr2num(get(handles.blue3_input,'String'));
blue4Value = mystr2num(get(handles.blue4_input,'String'));

MonitorRedValue = mystr2num(get(handles.MonitorRed_input,'String'));
MonitorGreenValue = mystr2num(get(handles.MonitorGreen_input,'String'));
MonitorBlueValue = mystr2num(get(handles.MonitorBlue_input,'String'));
monitorBlckValue = mystr2num(get(handles.monitorBlck_input,'String'));

notes = get(handles.notes_input,'String'); 

calibInput = [red1Value ,red2Value ,red3Value ,green1Value ,green2Value,green3Value,blue1Value,blue2Value,blue3Value,blue4Value,MonitorRedValue,MonitorGreenValue,MonitorBlueValue,monitorBlckValue,notes] ;

handles = datefixer(handles) ; % change date into number

dlmwrite(handles.calibrationFileName,[handles.dateNum, calibInput],...
    '-append','delimiter','\t','precision','%.2f') ;

CalibValues = dlmread(handles.calibrationFileName) ;

if size(CalibValues,1)>1
    
    CalibChange = abs(CalibValues(end,2:14)-CalibValues(end-1,2:14))./CalibValues(end-1,2:14) ;
    
    if any(CalibChange>.1) ;
        
        warnings = find(CalibChange>.1) ;
        
        Messages = {'red switch 1 calibration has changed more than 10% since last measured',...
            'red switch 2 calibration has changed more than 10% since last measured',...
            'red switch 3 calibration has changed more than 10% since last measured',...
            'green switch 1 calibration has changed more than 10% since last measured',...
            'green switch 2 calibration has changed more than 10% since last measured',...
            'green switch 3 calibration has changed more than 10% since last measured',...
            'blue switch 1 calibration has changed more than 10% since last measured',...
            'blue switch 2 calibration has changed more than 10% since last measured',...
            'blue switch 3 calibration has changed more than 10% since last measured',...
            'blue switch 4 calibration has changed more than 10% since last measured',...
            'monitor red calibration has changed more than 10% since last measured',...
            'monitor green calibration has changed more than 10% since last measured',...
            'monitor blue calibration has changed more than 10% since last measured',...
            'monitor blck/wht calibration has changed more than 10% since last measured'} ;
        
        finalMessage = cell(1,length(warnings)) ;
        for a=1:length(warnings) ; % for each warning
            finalMessage{a} = Messages{warnings(a)} ;
        end
        
        errordlg(finalMessage,'calibration warning') ; % display these
        
    end
end

function calib_file_button_Callback(hObject, eventdata, handles)
% hObject    handle to calib_file_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calib_file_button as text
%        str2double(get(hObject,'String')) returns contents of calib_file_button as a double

CalibValues = dlmread(handles.calibrationFileName);

temp = get(0,'screensize');
screenW = temp(3);
screenH = temp(4);

h = figure;
pos = round([.2*screenW .8*screenH .6*screenW  .4*screenH]);
set(h,'position',pos,'menubar','none','name','Calibration Table','numbertitle','off');

column_headings = {'date','red1','red2','red3','green1','green2','green3','blue1','blue2','blue3','blue4'};
t = uitable('parent',h,'ColumnName',column_headings,'Data',CalibValues,'units','normalized');
pos = [0.01 0.01 .98 .98];
set(t,'position',pos);

% --- Executes during object creation, after setting all properties.
function calib_file_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calib_file_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_calibration_button.
function load_calibration_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_calibration_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = datefixer(handles) ;

CalibValues = dlmread(handles.calibrationFileName) ;
    
if any(CalibValues(:,1) == handles.dateNum)
    ValuesRequested = find(CalibValues(:,1) == handles.dateNum,1,'last') ;
    toDisplay = CalibValues(ValuesRequested,2:end) ;

    set(handles.red1_input,'String',toDisplay(1));
    set(handles.red2_input,'String',toDisplay(2));
    set(handles.red3_input,'String',toDisplay(3));
    set(handles.green1_input,'String',toDisplay(4));
    set(handles.green2_input,'String',toDisplay(5));
    set(handles.green3_input,'String',toDisplay(6));
    set(handles.blue1_input,'String',toDisplay(7));
    set(handles.blue2_input,'String',toDisplay(8));
    set(handles.blue3_input,'String',toDisplay(9));
    set(handles.blue4_input,'String',toDisplay(10));
    set(handles.MonitorRed_input,'String',toDisplay(11));
    set(handles.MonitorGreen_input,'String',toDisplay(12));
    set(handles.MonitorBlue_input,'String',toDisplay(13));
    set(handles.monitorBlck_input,'String',toDispay(14));
    set(handles.notes_input,'String',toDisplay(15));
else
    errordlg('no calibration values recorded for this date on this rig')
end

handles = clearTables(handles);

function filter_attenuation_red_edit_Callback(hObject, eventdata, handles)
% hObject    handle to filter_attenuation_red_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filter_attenuation_red_edit as text
%        str2double(get(hObject,'String')) returns contents of filter_attenuation_red_edit as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function filter_attenuation_red_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_attenuation_red_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filter_attenuation_green_edit_Callback(hObject, eventdata, handles)
% hObject    handle to filter_attenuation_green_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filter_attenuation_green_edit as text
%        str2double(get(hObject,'String')) returns contents of filter_attenuation_green_edit as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function filter_attenuation_green_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_attenuation_green_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filter_attenuation_blue_edit_Callback(hObject, eventdata, handles)
% hObject    handle to filter_attenuation_blue_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filter_attenuation_blue_edit as text
%        str2double(get(hObject,'String')) returns contents of filter_attenuation_blue_edit as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function filter_attenuation_blue_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_attenuation_blue_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spotSize_edit_Callback(hObject, eventdata, handles)
% hObject    handle to spotSize_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spotSize_edit as text
%        str2double(get(hObject,'String')) returns contents of spotSize_edit as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function spotSize_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spotSize_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rod_collecting_area_edit_Callback(hObject, eventdata, handles)
% hObject    handle to rod_collecting_area_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rod_collecting_area_edit as text
%        str2double(get(hObject,'String')) returns contents of rod_collecting_area_edit as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function rod_collecting_area_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rod_collecting_area_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cone_collecting_area_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cone_collecting_area_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cone_collecting_area_edit as text
%        str2double(get(hObject,'String')) returns contents of cone_collecting_area_edit as a double
handles = clearTables(handles);

% --- Executes during object creation, after setting all properties.
function cone_collecting_area_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cone_collecting_area_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ND_filters_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ND_filters_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ND_filters_edit as text
%        str2double(get(hObject,'String')) returns contents of ND_filters_edit as a double

handles.NDFs = str2num(get(hObject,'String'));
handles = clearTables(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ND_filters_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ND_filters_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = clearTables(handles)
D = cell(4,10);
set(handles.rod_table,'Data',D);
set(handles.scone_table,'Data',D);
set(handles.mcone_table,'Data',D);
set(handles.lcone_table,'Data',D);


% --- Executes when entered data in editable cell(s) in rod_table.
function rod_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to rod_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%recreate conversion matrix
handles = calibrateLEDSpectra(hObject,handles);
handles = createConversionMatrix(hObject,handles);
guidata(hObject, handles);

%set display color
set(handles.rod_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});
set(handles.scone_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});
set(handles.mcone_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});
set(handles.lcone_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});

Ind = eventdata.Indices;
row = Ind(1);
col = Ind(2);
if row == 1 %P* supplied, calculate voltage, switch, and NDF
    calc_isom = 0;
else %Voltage, switch, or NDF edited, calculate P*
    calc_isom = 1;
end

D = get(hObject,'Data');

D = doConversionCalculation(handles,D,handles.conversionMatrix_rod,calc_isom);
set(hObject,'Data',D);

D_scone = get(handles.scone_table,'Data');
D_scone(2:end,:) = D(2:end,:);
D_scone = doConversionCalculation(handles,D_scone,handles.conversionMatrix_scone,1);
set(handles.scone_table,'Data',D_scone);

D_mcone = get(handles.mcone_table,'Data');
D_mcone(2:end,:) = D(2:end,:);
D_mcone = doConversionCalculation(handles,D_mcone,handles.conversionMatrix_mcone,1);
set(handles.mcone_table,'Data',D_mcone);

D_lcone = get(handles.lcone_table,'Data');
D_lcone(2:end,:) = D(2:end,:);
D_lcone = doConversionCalculation(handles,D_lcone,handles.conversionMatrix_lcone,1);
set(handles.lcone_table,'Data',D_lcone);


function D = doConversionCalculation(handles,D,conversionMatrix,calc_isom)
redSwitches = 1:3;
greenSwitches = 4:6;
blueSwitches = 7:10;

cols = size(D,2);
color = get(handles.colorSelector_menu,'Value');
switch color
    case 1, switches = redSwitches;
    case 2, switches = greenSwitches;
    case 3, switches = blueSwitches;
end
Nswitches = length(switches);

if calc_isom
    for c=1:cols
        V = D{2,c};        
        LEDswitch = switches(D{3,c});
        NDF = D{4,c};
        
        if ~isempty(V) && ~isempty(LEDswitch) && ~isempty(NDF)
            NDF_index = find(handles.NDFs==NDF);
            D{1,c} = V.*conversionMatrix(NDF_index,LEDswitch);
        end
    end
else
    Ncols = 0;
    choiceMatrix_sum = zeros(length(handles.NDFs),Nswitches);%NDFs x switches
    for c=1:cols
        Pstar = D{1,c};
        if ~isempty(Pstar)
            Ncols = Ncols + 1;
            choiceMatrix{c} = zeros(length(handles.NDFs),Nswitches); %NDFs x switches
            
            for s=1:Nswitches %switches
                for n=1:length(handles.NDFs)
                    LEDswitch = switches(s);                    
                    NDF_index = n;
                    choiceMatrix{c}(n,s) = Pstar./conversionMatrix(NDF_index,LEDswitch);                    
                end
            end
            %can't have V>10, so set them to Inf
            choiceMatrix{c}(find(choiceMatrix{c}>=10)) = Inf;
            choiceMatrix_sum = choiceMatrix_sum + abs(choiceMatrix{c})-1;
            minV(c) = min(min(abs(choiceMatrix{c}-1)));
        end
    end
    
    if get(handles.linkSwitch_checkbox,'Value')              
        minV_sum = min(min(abs(choiceMatrix_sum-Ncols)));
        [bestR,bestC] = find(abs(choiceMatrix_sum-Ncols) == minV_sum,1); %find V nearest 1                
        for c=1:Ncols
            D{2,c} = choiceMatrix{c}(bestR,bestC);
            D{3,c} = bestC;
            D{4,c} = handles.NDFs(bestR);
        end
    else
         for c=1:Ncols
            [bestR,bestC] = find(abs(choiceMatrix{c}-1) == minV(c),1); %find V nearest 1      
            D{2,c} = choiceMatrix{c}(bestR,bestC);
            D{3,c} = bestC;
            D{4,c} = handles.NDFs(bestR);
         end        
    end
end

% --- Executes on button press in linkSwitch_checkbox.
function linkSwitch_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to linkSwitch_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of linkSwitch_checkbox


% --- Executes when entered data in editable cell(s) in mcone_table.
function mcone_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to mcone_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%recreate conversion matrix
handles = calibrateLEDSpectra(hObject,handles);
handles = createConversionMatrix(hObject,handles);
guidata(hObject, handles);

%set display color
set(handles.rod_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});
set(handles.scone_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});
set(handles.mcone_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});
set(handles.lcone_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});

Ind = eventdata.Indices;
row = Ind(1);
col = Ind(2);
if row == 1 %P* supplied, calculate voltage, switch, and NDF
    calc_isom = 0;
else %Voltage, switch, or NDF edited, calculate P*
    calc_isom = 1;
end

D = get(hObject,'Data');

D = doConversionCalculation(handles,D,handles.conversionMatrix_mcone,calc_isom);
set(hObject,'Data',D);

D_scone = get(handles.scone_table,'Data');
D_scone(2:end,:) = D(2:end,:);
D_scone = doConversionCalculation(handles,D_scone,handles.conversionMatrix_scone,1);
set(handles.scone_table,'Data',D_scone);

D_rod = get(handles.rod_table,'Data');
D_rod(2:end,:) = D(2:end,:);
D_rod = doConversionCalculation(handles,D_rod,handles.conversionMatrix_rod,1);
set(handles.rod_table,'Data',D_rod);

D_lcone = get(handles.lcone_table,'Data');
D_lcone(2:end,:) = D(2:end,:);
D_lcone = doConversionCalculation(handles,D_lcone,handles.conversionMatrix_lcone,1);
set(handles.lcone_table,'Data',D_lcone);



% --- Executes when entered data in editable cell(s) in scone_table.
function scone_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to scone_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%recreate conversion matrix
handles = calibrateLEDSpectra(hObject,handles);
handles = createConversionMatrix(hObject,handles);
guidata(hObject, handles);

%set display color
set(handles.rod_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});
set(handles.scone_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});
set(handles.mcone_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});
set(handles.lcone_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});

Ind = eventdata.Indices;
row = Ind(1);
col = Ind(2);
if row == 1 %P* supplied, calculate voltage, switch, and NDF
    calc_isom = 0;
else %Voltage, switch, or NDF edited, calculate P*
    calc_isom = 1;
end

D = get(hObject,'Data');

D = doConversionCalculation(handles,D,handles.conversionMatrix_scone,calc_isom);
set(hObject,'Data',D);

D_rod = get(handles.scone_table,'Data');
D_rod(2:end,:) = D(2:end,:);
D_rod = doConversionCalculation(handles,D_rod,handles.conversionMatrix_rod,1);
set(handles.rod_table,'Data',D_rod);

D_mcone = get(handles.mcone_table,'Data');
D_mcone(2:end,:) = D(2:end,:);
D_mcone = doConversionCalculation(handles,D_mcone,handles.conversionMatrix_mcone,1);
set(handles.mcone_table,'Data',D_mcone);

D_lcone = get(handles.lcone_table,'Data');
D_lcone(2:end,:) = D(2:end,:);
D_lcone = doConversionCalculation(handles,D_lcone,handles.conversionMatrix_lcone,1);
set(handles.lcone_table,'Data',D_lcone);

% --- Executes when entered data in editable cell(s) in lcone_table.
function lcone_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to lcone_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%recreate conversion matrix
handles = calibrateLEDSpectra(hObject,handles);
handles = createConversionMatrix(hObject,handles);
guidata(hObject, handles);

%set display color
set(handles.rod_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});
set(handles.scone_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});
set(handles.mcone_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});
set(handles.lcone_table,'ForegroundColor',handles.colors{get(handles.colorSelector_menu,'Value')});

Ind = eventdata.Indices;
row = Ind(1);
col = Ind(2);
if row == 1 %P* supplied, calculate voltage, switch, and NDF
    calc_isom = 0;
else %Voltage, switch, or NDF edited, calculate P*
    calc_isom = 1;
end

D = get(hObject,'Data');

D = doConversionCalculation(handles,D,handles.conversionMatrix_lcone,calc_isom);
set(hObject,'Data',D);

D_scone = get(handles.scone_table,'Data');
D_scone(2:end,:) = D(2:end,:);
D_scone = doConversionCalculation(handles,D_scone,handles.conversionMatrix_scone,1);
set(handles.scone_table,'Data',D_scone);

D_mcone = get(handles.mcone_table,'Data');
D_mcone(2:end,:) = D(2:end,:);
D_mcone = doConversionCalculation(handles,D_mcone,handles.conversionMatrix_mcone,1);
set(handles.mcone_table,'Data',D_mcone);

D_rod = get(handles.rod_table,'Data');
D_rod(2:end,:) = D(2:end,:);
D_rod = doConversionCalculation(handles,D_rod,handles.conversionMatrix_rod,1);
set(handles.rod_table,'Data',D_rod);



function blue4_input_Callback(hObject, eventdata, handles)
% hObject    handle to blue4_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blue4_input as text
%        str2double(get(hObject,'String')) returns contents of blue4_input as a double
handles = clearTables(handles);


% --- Executes during object creation, after setting all properties.
function blue4_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blue4_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in throughRGC_button.
function throughRGC_button_Callback(hObject, eventdata, handles)
% hObject    handle to throughRGC_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of throughRGC_button

set(handles.cone_collecting_area_edit,'String',handles.CollectingArea.cone_throughRetina);

    


% --- Executes on button press in PRdirect_button.
function PRdirect_button_Callback(hObject, eventdata, handles)
% hObject    handle to PRdirect_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PRdirect_button

set(handles.cone_collecting_area_edit,'String',handles.CollectingArea.cone);



function MonitorGreen_input_Callback(hObject, eventdata, handles)
% hObject    handle to MonitorGreen_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MonitorGreen_input as text
%        str2double(get(hObject,'String')) returns contents of MonitorGreen_input as a double


% --- Executes during object creation, after setting all properties.
function MonitorGreen_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MonitorGreen_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MonitorRed_input_Callback(hObject, eventdata, handles)
% hObject    handle to MonitorRed_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MonitorRed_input as text
%        str2double(get(hObject,'String')) returns contents of MonitorRed_input as a double


% --- Executes during object creation, after setting all properties.
function MonitorRed_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MonitorRed_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MonitorBlue_input_Callback(hObject, eventdata, handles)
% hObject    handle to MonitorBlue_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MonitorBlue_input as text
%        str2double(get(hObject,'String')) returns contents of MonitorBlue_input as a double


% --- Executes during object creation, after setting all properties.
function MonitorBlue_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MonitorBlue_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function notes_input_Callback(hObject, eventdata, handles)
% hObject    handle to notes_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of notes_input as text
%        str2double(get(hObject,'String')) returns contents of notes_input as a double


% --- Executes during object creation, after setting all properties.
function notes_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to notes_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function monitorBlck_input_Callback(hObject, eventdata, handles)
% hObject    handle to monitorBlck_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of monitorBlck_input as text
%        str2double(get(hObject,'String')) returns contents of monitorBlck_input as a double


% --- Executes during object creation, after setting all properties.
function monitorBlck_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to monitorBlck_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
