function playtime(z,numfile,struct,strind)
% playtime(z,numfile,struct)
%
% CHANGES TO THIS VERSION
% Version 3.3  Keeps users default for label and decimate during the
% parsing of the whole cell.
% Version 3.2
% Fixed bug.  If you have saved CellInfo, and you enter playtime and decide to return to epochs
% without saving that epochcondition or familycondition than it use to save that condition anyway. 
% Now it doesn't, you have to actually hit save analysis to add a condition to the CellInfo structure.
%
% New addition:  Threshold setting in the histograms pull-down menu.  Can look at the max or
% min histogram, set a threshold, and you will be returned a vector of epochs that fall
% above the threshold.  A new interface pops up, that lets you look at each of these epochs 
% individually.  Reject the ones you want to reject.  When you close the window, the vector of
% all the epochs that you want to reject will be sent to playtime, and a button will appear
% in the upper right 'Deselect More', this will activate your selections, so they no longer appear
% in the epoch list on the left.
% Version 3.0
% Biggest change is the implementation of CellInfo.  Single and SinCondition no longer 
% exist, instead we have two different structures within CellInfo, EpochCondition and 
% FamilyCondition.  For more info on the CellInfo structure see the end of this help file
% or type help StructCellInfo at the command line.  Also familystructure has been removed, this
% is because the transitions between epochs and playtime or familysearch and playtime are
% now the same, ie. the CellInfo structure is always created in either epochs or 
% familysearch, if it doesn't already exist, and this structure is passed as a variable to 
% playtime.  Both epochs and familysearch load the epochnumbers they have selected in the 
% corresponding structure, Epoch or Family, and playtime just adds to these structures.
% The CellInfo is either saved to the same file it was in to begin with, or if a structure
% file has not been created yet, it loads CellInfo to the command line, and the user saves
% it to a file using the command
%
% >> save newfilename CellInfo
%
% Data loading is left up to the user, commands LoadCIData andRmCIData are used.
% Returning to Epochs after saving analysis works, continues to add to EpochCondition or 
% FamilyCondition.  When you exit, it will have either updated the original structure file
% with all additional Condition structures, or will have updated the CellInfo sitting in 
% the command window with all of the additional Condition structures.
%
% Version 2 Basically,  the only major difference between this and the earlier version is 
% the command familystructure, in the finalsave case.  This calls up the script familystructure,
% which checks to see if familyplay was run and put the structure SinCondition into the 
% workspace.  If so, it loads this structure into the structure SingleCondition as the 
% field, family.  Then it loads the epoch data into the structure.  This can be removed
% at a later time by using RmCIData, and re-entered again using LoadCIData.
%
% Last of 3 interfaces, can be entered from epochs or directly.  To load directly
% must have either a file created by epochs or a structure in the SingleCondition 
% format (see below)
% 
% To load the last EpochCondition of a structure, three inputs are necessary 
% >> playtime('init','matfilename',1)
% 'init' starts the interface, and 'matfilename' is the file with the structure in it  
% (the .mat part is not necessary), and 1 tells it is an EpochCondition
%
% To load the last FamilyCondition of a structure, three inputs are necessary 
% >> playtime('init','matfilename',2)
% 'init' starts the interface, and 'matfilename' is the name of the file created by Epochs 
% (the .mat part is not necessary), and 2 telss it is a FamilyCondition
%
% To load a previously analyzed structure, where you want to look at a Condition other
% than the latest one, use the same general syntax, but add one more input.  This is to 
% give the index number of the appropriate condition
%
% >> playtime('init','structurename',1,2)  % for EpochCondition(2)
% >> playtime('init','structurename',2,3)  % for FamilyCondition(3)
%
% This is the breakdown of the structure returned at the end of playtime:
%
%
%		CellInfo.CellFile 				File Name of original data file, incl. path
%		CellInfo.Comment 				Filled in by user
%		CellInfo.Label 					Filled in by user
%		CellInfo.CellType 				Filled in by user
%		CellInfo.Rig 					Filled in by user
%		CellInfo.OutputScaleFactor		Filled in by user
%		CellInfo.OutputConfiguration	Filled in by user
%		CellInfo.NDFConfiguration		Filled in by user
%		CellInfo.UserInfo				Can be structure, user defined
%
%	SingleCondition becomes either CellInfo.EpochCondition or CellInfo.FamilyCondition

%		SingleCondition.Segnum			Segment number analyzed	    
%		SingleCondition.Epochnumbers  	List of epoch numbers as selected in Epochs
%  		SingleCondition.ExcludeEpochs		List of 0's and 1's, length of Epochnumbers, 
%										0's deselected during analysis
%		SingleCondition.Comment			Filled in by user
%		SingleCondition.Plotpref		plotpref, same as in analyzethis		
%
%		SingleCondition.SearchCrit		datanum2, same as in analyzethis
%		SingleCondition.SearchPara		datanum1, same as in analyzethis
%
%	These are unique to FamilyCondition
%		SingleCondition.EpochNumbers	index #'s for epochs in chosen families
%		SingleCondition.FamilyFlag		size of EpochNumbers, gives each family a number
%		SingleCondition.FamilyStep		values of step parameter for epcohs in EpochNumbers
%		SingleCondition.FamilyStepGuide	name of parameter that steps (ie. amp)
%		SingleCondition.FamilyCues		values of cue parameter for families selected
% 		SingleCondition.FamilyCueGuide	name of parameter(s) that change from family to family
%
%
%  To save this file, use the syntax:
%  >> save filename SingleCondition
%  You will probably want to clear the workspace once you have saved the structure.
%  To look at the structure once it has been saved and cleared, use the syntax
%  >> load filename
%  The command whos will give you the structure name and other info, printing CellInfo
%  at the command line will list the fields.
%
%  To add more analysis tools please see analysisplay.m (can be accessed as a regular file,
%  or by printing <help analysisplay> at the command line.

% Other Variables used in playtime
%  fileptr, segnum, and epochnumber are all used in ITC functions
% 
%  structtype is 1 for EpochCondition, 2 for FamilyCondition
%  strindflag is the Condition #
%  newflag lets the program know if this is the first time this file has been analyzed, so
%		the appropriate data can be added into CellInfo.
%  listepochs is a vector of the epoch numbers selected in Epochs
%  searchcrit is the criteria used to choose the epochs (boxes checked in Epochs - datanum2)
%  searchpara is the parameters (actual values of criteria) used to select epochs - datnum1
%  newaxx and newaxy are the x and y coordinates of the graph when the user changes the axes
%  analine is the temporary vector with analysis data in it
%  groupepochs is the matrix (cell) with all of the epoch data currently being used in it
%  stdev is the standard deviation vector 
%  tempepochs is used to convert the epochnumbers in a previously saved structure
%  Created January 2001 MKMK
%  2nd version completed June 2001 MKMK
%  Oct 2004 MKMK Changed ITCFinishAnalysis to use evalin in order to get around
%  wierd bug in Matlab version 7
%  Dec 2004 MKMK Changed how NDFConfiguration and OutputScaleFactor are
%  saved - now a number instead of a string
%  May 2005 MKMK Turns out it wasn't matlab; Fred fixed the bug in
%  ITCFinishAnalysis, so we don't have to use evalin anymore.

persistent fp fileptr plotpref segnum otherseg structname listepochs epochnumber newaxx newaxy average average2 groupepochs groupepochs2 stdev epochplot epochplot2 histoplot structtype strindflag newflag CellInfo condense rejectplots newrejectlist saved
% Removed tempepochs
if nargin < 1
	z='init';
	prompt = {'Name of structure file','1 = EpochCondition, 2 = FamilyCondition','Condition #'};
	wintitle = 'Structure info';
	def = {'cell032901c1','1','1'};
	info = inputdlg(prompt,wintitle,1,def);
	numfile = info{1};
	structtype = info{2};
	structtype = str2num(structtype);
	strindflag = info{3};
	strindflag = str2num(strindflag);
	newflag = 0;
	load(numfile)	
	structname = numfile;
end

if nargin > 1
	%numfile
	load(numfile)
	structname = numfile;
end

if nargin == 3
	structtype = struct;
 	strindflag = 1;
	newflag = 1;
end

if nargin == 4
	structtype = struct;
	strindflag = strind;
	newflag = 0;
end


	% If this is the first time analyzing data from this cell, then Epochs or 
	% Familysearch has created CellInfo with only two fields, CellFile and the Condition
	% field.  The program knows if this is the first time by whether there are three or
	% four inputs, and sets newflag accordingly, newflag = 1 will prompt a request for 
	% more info at the end of playtime, to fill in the rest of the fields of CellInfo.  
switch z

case 'init'
saved = 0;
fp = CellInfo.CellFile;
if structtype == 1
	shortsin = CellInfo.EpochCondition(strindflag);
	tempepochs = shortsin.EpochNumbers;
	plotpref = shortsin.PlotPref;
	segnum = shortsin.SegNum;
	if ~isempty(shortsin.ExcludeEpochs)
		exclude = shortsin.ExcludeEpochs;
		finalepochs = tempepochs(find(exclude==0));
	else 
		finalepochs = tempepochs;
	end
elseif structtype == 2
	shortfam = CellInfo.FamilyCondition(strindflag);
	conden = ~isnan(shortfam.EpochNumbers);
	tempepochs = shortfam.EpochNumbers(conden);
	plotpref = shortfam.PlotPref;
	segnum = shortfam.SegNum;
	if ~isempty(shortfam.ExcludeEpochs)
		exclude = shortfam.ExcludeEpochs;
		finalepochs = tempepochs(find(exclude==0));
	else
		finalepochs = tempepochs;
	end
end	
listepochs = finalepochs;
%disp('Entering Playtime')
%  If there is no file, go back to begining.
if isempty(fp)
	playtime exit
	icaphoton
else
	[fileptr, ecode] = ITCInitializeAnalysis(2e6, fp);
	%disp('initialized');
	if ecode~=0
		playtime error
	end
end
%  These list the epochs that were chosen in Epochs screen.  The one called str is used by
%  the listbox, and will change when user wants to remove an epoch from the data to be 
%  analyzed.  The one called originalstr is divided in half, if it is long, so that all
%  numbers are displayed, and it is a permanent record of the epochs chosen originally from
%  the Epochs interface.
%  If we are looking at something analyzed previously, str and originalstr will be different
%  already.  
str = listepochs;
%  MAY NEED TO CHANGE THIS IF I HAVE EXCLUDE EMPTY AS THE DEFAULT IN THE STRUCTURE, DEPENDS
%  ON HOW I GET MORE STRUCTURES IN
if exist('exclude') 
	if ~isempty(exclude)
		originalstr = tempepochs;
		origilong = length(originalstr);
	end
else
	originalstr = listepochs;
	origilong = length(listepochs);
end
if origilong<50
	oristr = originalstr(1:end);
	oristr2 = [];
else 
	oristr = originalstr(1:50);
	oristr2 = originalstr(51:end);
end
load playtime
h0 = figure('CloseRequestFcn','playtime close',...
	'Color',[0.8 0.8 0.8], ...
	'Colormap',mat0, ...
	'Position',[60 40 800 650], ...
	'Tag','Fig1');
%  This is the drop down menu for additional analysis functions.  Add a uimenu and an 
%  appropriate callback.
h1 = uimenu('Parent',h0, ...
	'Label','select_analysis', ...
	'Tag','playone');
	 uimenu(h1,'Label','Average','Callback','playtime average');
	 uimenu(h1,'Label','Std Deviaton','Callback','playtime deviation');
h1 = uimenu('Parent',h0, ...
	'Label','histograms', ...
	'Tag','playone');
	 h1a = uimenu(h1,'Label','Set Threshold');
	 uimenu(h1a,'Label','Maximum','Callback','playtime threshold','Userdata','playtime maximum');
	 uimenu(h1a,'Label','Minimum','Callback','playtime threshold','Userdata','playtime minimum');
	 h1b = uimenu(h1,'Label','Minimum','Callback','playtime minimum');
	 h1c = uimenu(h1,'Label','Maximum','Callback','playtime maximum');
	 h1d = uimenu(h1,'Label','Std Deviation','Callback','playtime stdev');
	 
h1 = axes('Parent',h0, ...
	'Units','pixels', ...
	'CameraUpVector',[0 1 0], ...
	'Color',[1 1 1], ...
	'ColorOrder',mat1, ...
	'Position',[210 144 454 436], ...
	'Tag','playaxes', ...
	'XColor',[0 0 0], ...
	'YColor',[0 0 0], ...
	'ZColor',[0 0 0]);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[0.5 0.5 10], ...
	'Tag','Axes1Text4', ...
	'VerticalAlignment','cap');
set(get(h2,'Parent'),'XLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[0.5 0.5 10], ...
	'Rotation',90, ...
	'Tag','Axes1Text3', ...
	'VerticalAlignment','baseline');
set(get(h2,'Parent'),'YLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','right', ...
	'Position',mat2, ...
	'Tag','Axes1Text2', ...
	'Visible','off');
set(get(h2,'Parent'),'ZLabel',h2);
h2 = text('Parent',h1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',mat3, ...
	'Tag','Axes1Text1', ...
	'VerticalAlignment','bottom');
set(get(h2,'Parent'),'Title',h2);
%  Title of Interface
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Units','points', ...
	'FontSize',24, ...
	'Position',[200 560 130 40], ...
	'String','Playtime', ...
	'Style','text');
%  Button to return to Epochs
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','playtime return', ...
	'Position',[400 570 150 30], ...
	'String','Return to Epochs', ...
	'Tag','Pushbutton2');
%  Lists the epochs selected in previous interface, removing from list as directed
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Max',10, ...
	'Min',1, ...
	'Position',[15 500 50 100], ...
	'String',str, ...
	'Style','listbox', ...
	'Tag','epochstr');
%  Button to remove epoch from current analysis
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','playtime remove', ...
	'Position',[15 450 150 30], ...
	'String','Deselect highlighted epoch(s)', ...
	'Tag','Pushbutton1');
%  Lists the original epochs selected in previous interface, does not change
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Position',[670 10 40 600], ...
	'String',oristr, ...
	'Style','text',...
	'Tag','oldlist');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Position',[715 10 40 600], ...
	'String',oristr2, ...
	'Style','text',...
	'Tag','oldlist2');
%  Displays epoch number currently on graph
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Units','points', ...
	'Position',[110 550 30 20], ...
	'String','', ...
	'Style','text', ...
	'Tag','curepoch');
%  Label for previous box
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Units','points', ...
	'Position',[80 570 100 20], ...
	'String','Current Epoch', ...
	'Style','text');	
%  Displays epoch segment currently on graph
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Units','points', ...
	'Position',[110 510 30 20], ...
	'String',segnum, ...
	'Style','text', ...
	'Tag','tdataseg');
%  Label for previous box
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Units','points', ...
	'Position',[80 530 100 20], ...
	'String','Current Segment', ...
	'Style','text');	
%  Radio buttons give user option of leaving analysis or previous plots on graph
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Position',[15 410 150 30], ...
	'String','Leave analysis on graph', ...
	'Style','radiobutton', ...
	'Tag','analysis',...
	'Userdata',0);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Position',[15 370 150 30], ...
	'String','Leave all plots on graph', ...
	'Style','radiobutton', ...
	'Tag','allplots');
%  Allows the user to change the axes by using cursors on the graph
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','playtime diffaxes', ...
	'Position',[15 330 150 30], ...
	'String','Set new axes', ...
	'Tag','Pushbutton1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','playtime resetaxes', ...
	'Position',[15 290 150 30], ...
	'String','Change axes', ...
	'Tag','Pushbutton2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','playtime defaultaxes', ...
	'Position',[15 250 150 30], ...
	'String','Default Axes', ...
	'Tag','Pushbutton3');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','playtime clear', ...
	'Position',[15 210 150 30], ...
	'String','Clear Graph', ...
	'Tag','Pushbutton4');
%  Buttons allow user to flip through epochs or go to a specific one
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','playtime jump', ...
	'Position',[210 60 100 50], ...
	'String','Previous', ...
	'Tag','Previous');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','playtime jump', ...
	'Position',[320 60 100 50], ...
	'String','Next', ...
	'Tag','Next');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','playtime jump', ...
	'Position',[430 60 100 50], ...
	'String','Or go to epoch #', ...
	'Tag','Choose');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','playtime jump', ...
	'Position',[550 70 70 30], ...
	'Style','edit', ...
	'Tag','newepoch');
%  Button saves any analysis done to file
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','playtime finalsave',...
	'Position',[495 10 154 39], ...
	'String','Save Current Analysis', ...
	'Tag','Pushbutton3');
%  Create the matrix of all epochs that will be used for calculations
%  groupepochs is a cell where the 1st matrix is the data and the 
%  1st cell of the 2nd matrix is list of epoch numbers
datasize = length(listepochs);
groupepochs=cell(2,datasize);
for ep = 1:datasize
	each = listepochs(ep);
	[dataepochs, ecode1] = ITCReadEpoch(each, segnum, fileptr);
	groupepochs(1,ep)={dataepochs};
end
groupepochs{2,1}=listepochs;
cat(1,groupepochs{2,1});
groupepochs;
epochnumber = listepochs(1);
epochnow = findobj(gcf,'tag','curepoch');
set(epochnow,'string',listepochs(1));
segnow = findobj(gcf,'tag','tsegnum');
set(segnow,'string',segnum);
[EpochData, ecode1] = ITCReadEpoch(epochnumber, segnum, fileptr);
if ecode1 ~= 0
	playtime error
end
[NumSegments, ecode1] = ITCGetEpochSegments(epochnumber,fileptr);
if ecode1 ~= 0
	playtime error
end
if NumSegments > 1 & segnum == 0 
	otherseg = 1;
elseif NumSegments > 1 & segnum == 1
	otherseg = 0;
else otherseg = [];
end
if ~isempty(otherseg)
	groupepochs2=cell(2,datasize);	
	for ep = 1:datasize
		each = listepochs(ep);
		[dataepochs2, ecode10] = ITCReadEpoch(each, otherseg, fileptr);
		if ecode10 ~= 0
			playtime error
		end
		groupepochs2(1,ep)={dataepochs2};
	end
end
%  Plotpref
%  1 means only one segment, no stim shown
%  2 Other segment shown, stim not shown
%  3 Other segment not shown, stim shown
%  4 Both segment and stim shown
%  Stim will be in a separate figure, since no analysis is actually done with it.  
if plotpref == 1 
	%disp('business as usual');
	plot(EpochData,'c-')
elseif plotpref == 3 
	playtime stimfigure
elseif plotpref == 2 | plotpref == 4
	playtime newfigure2
end
playfig = findobj('tag','Fig1');
figure(playfig);
%  If I were to do this again, I would have a separate command that would call up the epochdata and 
%  zero it, and than use this to call it for plotting.  May still do this.
case 'stimfigure'
	%  disp('At stimfigure');
	%  Adds a separate figure for the stim, only one segment shown
	[EpochData, ecode9] = ITCReadEpoch(epochnumber, segnum, fileptr);
	[StimData, ecode20] = ITCReadEpochStm(epochnumber, segnum, fileptr);
		if ecode20 ~= 0
			playtime error
		end
		if ecode9 ~= 0
			playtime error
		end	
		plot(EpochData,'-c')
		stimfig = figure('Position',[10 30 400 200]);
		plot(StimData,'-c')
case 'newfigure2'
	[EpochData, ecode9] = ITCReadEpoch(epochnumber, segnum, fileptr);
	if ecode9 ~= 0
		playtime error
	end
	if ~isempty(otherseg)
		[EpochData2, ecode10] = ITCReadEpoch(epochnumber, otherseg, fileptr);
		if ecode10 ~= 0
			playtime error
		end
		epochplot = subplot('Position',[.35 .60 .55 .25]);
		plot(EpochData,'-c','EraseMode','none')
		title('Segment 0');
		epochplot2 = subplot('Position',[.35 .30 .55 .25]);
		plot(EpochData2,'-y')
		title('Segment 1');
	elseif isempty(otherseg)
		plot(EpochData,'-c')
		disp('No other segment for this epoch');
	end
	if plotpref == 4	
		[StimData, ecode20] = ITCReadEpochStm(epochnumber, segnum, fileptr);
		if ecode20 ~= 0
			playtime error
		end
		if ~isempty(otherseg)	
			[StimData2, ecode21] = ITCReadEpochStm(epochnumber, otherseg, fileptr);
			if ecode21 ~= 0
				playtime error
			end		
			stimfig = figure('Position',[10 30 400 200]);
			stimsplot = subplot('Position',[.15 .60 .85 .30]);
			plot(StimData,'-c')
			title('Stimulus 0');
			stimplot2 = subplot('Position',[.15 .10 .85 .30]);
			plot(StimData2,'-y')
			title('Stimulus 1');
			%disp('suposedly worked');
		elseif isempty(otherseg) 
			stimfig = figure('Position',[10 30 400 200]);
			plot(StimData,'-c')
		end	
	end	
case 'newfigure3'
	%disp('At newfigure3');	
	%find out if there is analysis up that should be left up, and/or leave up 
	%all plots
	leavplot = findobj(gcf,'tag','analysis');
	analysis = get(leavplot,'value');
	set(leavplot,'Userdata',1);
	leavplots = findobj(gcf,'tag','allplots');
	allplots = get(leavplots,'value');
	% Leave up analysis only
	if analysis == 1 & allplots == 0
		%disp('newfigure3 averages')
		if isempty(otherseg)
			delete(gca)
			axes('Position',[.28 .30 .55 .55]);
			analine = plot(average,'r-','EraseMode','none');
			hold on
		elseif ~isempty(otherseg)	
			subplot(epochplot)
			cla
			epochplot = subplot('Position',[.35 .60 .55 .25]);
			analine = plot(average,'r-','EraseMode','none');
			hold on
			subplot(epochplot2)
			cla
			epochplot2 = subplot('Position',[.35 .30 .55 .25]);
			analine2 = plot(average2,'r-');	
			hold on
		end
	% Nothing left, erase graph	
	elseif analysis == 0 & allplots == 0
		delete(gca)
		delete(gca)
		axes('Position',[.28 .30 .55 .55]);	
	% leave up previous graph	
	elseif allplots == 1
		if isempty(otherseg)
			hold on
		else
		subplot(epochplot)
		hold on
		subplot(epochplot2)
		hold on
		end
	end
	%epochnumber;
	%segnum;
	%otherseg;
	[EpochData, ecode9] = ITCReadEpoch(epochnumber, segnum, fileptr);
	if ecode9 ~= 0
		playtime error
	end
	if ~isempty(otherseg)
		[EpochData2, ecode10] = ITCReadEpoch(epochnumber, otherseg, fileptr);
		if ecode10 ~= 0
			playtime error
		end
		epochplot = subplot('Position',[.35 .60 .55 .25]);
		plot(EpochData,'-c','EraseMode','none')
		title('Segment 0');
		epochplot2 = subplot('Position',[.35 .30 .55 .25]);
		plot(EpochData2,'-y')
		title('Segment 1');
	elseif isempty(otherseg)
		%axes('Position',[.28 .30 .55 .55]);
		plot(EpochData,'-c')
	end	
case 'jump'
	%disp('jumped');
	%find out if there is analysis up that should be left up, and/or leave up 
	%all plots
	leavplot = findobj(gcf,'tag','analysis');
	analysis = get(leavplot,'value');
	set(leavplot,'Userdata',1);
	leavplots = findobj(gcf,'tag','allplots');
	allplots = get(leavplots,'value');
	%find out if we are going forwards, backwards or jumping to a certain epoch
	%listepochs is the list of epochs that user is analyzing
	epochnow = findobj(gcf,'tag','curepoch');
	jumpdir = gcbo;
	movedir = get(jumpdir,'Tag');
	movedirf = strncmp(movedir,'Previous',1);
	now=epochnumber==listepochs;
	%  now is a matrix of zeros with a one marking the epoch showing before the user jumped
	%  old will tell us the index number of this epoch
	%  If the user is at an epoch that has been removed from the list, than any will not
	%  have any ones, so we go to the first epoch that has not been removed.
	if ~any(now)
		now(end) = 1;
	end
	old=find(now);
	max1=length(listepochs);
	if movedirf == 1
		if old - 1 < 1
			epochnumber = listepochs(max1);
			set(epochnow,'string',epochnumber);
			new = max1;
		else
			epochnumber = listepochs(old - 1);	
			set(epochnow,'string',epochnumber);
			new = old - 1;
		end
	elseif movedirf == 0
		movedirf2 = strncmp(movedir,'Next',1);
		if movedirf2 == 1
			if old + 1 > max1
				epochnumber = listepochs(1);
				set(epochnow,'string',epochnumber);
				new = 1;
			else
				epochnumber = listepochs(old + 1);
				set(epochnow,'string',epochnumber);
				new = old + 1;
			end
		elseif movedirf2 == 0
			changenum = findobj(gcf,'tag','newepoch');
			newepoch = get(changenum,'String');
			if isempty(newepoch)
				error('Must enter an epoch #');
				playtime error
			else
				newepoch = str2num(newepoch);
				existnum = find(listepochs == newepoch);
				if isempty(existnum)
					error('Epoch not available, choose an epoch from the list');
				else
					epochnumber = newepoch;
					new=existnum;
					set(epochnow,'string',epochnumber);
				end
			end
		end
	end
	[EpochData, ecode1] = ITCReadEpoch(epochnumber, segnum, fileptr);
	if ecode1 ~= 0
		playtime error
	end
	if plotpref == 1 | plotpref == 3
		if analysis == 1 & allplots == 0
			%disp('should leave up just averages')
			cla
			analine = plot(average,'r-');
		elseif analysis == 0 & allplots == 0
			clear analine
			hold off
		elseif allplots == 1 
			hold on
		else
			clear analine
			hold on
		end
		plot(EpochData,'c-')
	elseif plotpref == 2 | plotpref == 4
		playtime newfigure3
	end
	%  Check for standard deviation, and update window if present
	standev = findobj(gcf,'tag','standev');
	if ~isempty(standev)
		set(standev,'String',stdev(new));
	end
case 'remove'
	%  Remove is the actual epochnumber to be removed, epochdel is the index to 
	%  that number.
	%  Remove data cells
	epochstr = findobj(gcf,'tag','epochstr');
	epochdel = get(epochstr,'value');
	listrem = ones(size(listepochs));
	listrem(epochdel) = 0;
	listind = find(listrem);
	listnum = sum(listrem);
	tempgroup = groupepochs;	
	groupend = length(groupepochs);                      
	for i=1:listnum
		ind = listind(i);
		tempgroup2{1,i} = tempgroup{1,ind};
	end
	groupepochs = tempgroup2;
	% Must recreate cell with correct epoch #s
	groupepochs{2,1}=listepochs(listrem==1);
	cat(1,groupepochs{2,1});
	%  Now remove from gui
	listepochs = listepochs(listrem==1);
	set(epochstr,'Value',1);
	set(epochstr,'String',listepochs);
case 'diffaxes'
	newaxx=zeros(2,1);
	newaxy=zeros(2,1);
    n = 0;
    but = 1;
    for but=1:2
    	[newaxx(but),newaxy(but),but] = ginput(1);
      	if but~=1 & n<2
           but=1;
           disp('Pick two points please.');
       	end
      	if n==2
      	   but=2;
      	else
      	   n = n + 1;
       	  newaxx(:);
       	  newaxy(:);
      	end
   	end
case 'resetaxes'
	if newaxx(1)>newaxx(2)
      xmin=newaxx(2);
      xmax=newaxx(1);
   elseif newaxx(1)<newaxx(2)
      xmin=newaxx(1);
      xmax=newaxx(2);
   end
   if newaxy(1)>newaxy(2)
      ymin=newaxy(2);
      ymax=newaxy(1);
   elseif newaxy(1)<newaxy(2)
      ymin=newaxy(1);
      ymax=newaxy(2);
   end
   axis([xmin xmax ymin ymax])
case 'defaultaxes'
    set(gca,'XLimMode','auto','YLimMode','auto');
case 'clear'
	cla
case 'average'
	otherseg;
	groupepochs;
	size(groupepochs);
	want = size(groupepochs,2);
	condense = cat(1,groupepochs{1,1:want});
	size(condense);
	average = mean(condense);
	if isempty(otherseg)
		analine = plot(average,'r-','EraseMode','xor');
		hold on
	elseif ~isempty(otherseg)	
		want2 = size(groupepochs2,2);
		condense2 = cat(1,groupepochs2{1,1:want});
		average2 = mean(condense2);
		subplot(epochplot)
		analine = plot(average,'r-');
		subplot(epochplot2)
		analine2 = plot(average2,'r-');	
	end
case 'deviation'
	want = size(groupepochs,2);
	condense = cat(1,groupepochs{1,1:want});
	size(condense);
	condense = (condense)';
	stdev = std(condense);
	h0 = gcf;
	now=epochnumber==listepochs;
	ind=find(now);
	h1 = uicontrol('Parent',h0, ...
            'BackgroundColor',[0.8 0.8 0.8], ...
			'Units','points', ...
			'Position',[75 210 60 20], ...
			'String',stdev(ind), ...
			'Style','text', ...
			'Tag','standev', ...
			'Userdata','');
	h1 = uicontrol('Parent',h0, ...
            'BackgroundColor',[0.8 0.8 0.8], ...
			'Units','points', ...
			'Position',[15 210 60 20], ...
    		'String','Std Dev', ...
			'Style','text');	
case 'maximum'
	% Want data to be offset to zero			
	% want = How many epochs, condense is the data
	want = size(groupepochs,2);
	condense = cat(1,groupepochs{1,1:want});
	BaseLinePoints = zeros(want,1);
	each = 1;
	while each < want
		foo = listepochs(each);
		[PrePoints,ErrorFlag] = ITCGetStmPrePts(foo, 0, 0, fileptr);
		if ErrorFlag 
			fprintf(1, 'error1 reading from file \n');
			break
		end
		BaseLinePoints(each) = PrePoints;
		each = each + 1;
	end	
	% BaseLinePoints is the number of data points in the PrePts section
	% Get the data in the PrePts section
	trial = condense(:,1:BaseLinePoints);
	BaseLineAverage = mean(trial,2);
	for i = 1:want
		condense(i,:) = condense(i,:) - BaseLineAverage(i);
	end
	condens = abs(condense)';
	histoplot = max(condens);
	playtime histogram
case 'minimum'			
	% Want data to be offset to zero			
	% want = How many epochs, condense is the data
	want = size(groupepochs,2);
	condense = cat(1,groupepochs{1,1:want});
	BaseLinePoints = zeros(want,1);
	each = 1;
	while each < want
		foo = listepochs(each);
		[PrePoints,ErrorFlag] = ITCGetStmPrePts(foo, 0, 0, fileptr);
		if ErrorFlag 
			fprintf(1, 'error1 reading from file \n');
			break
		end
		BaseLinePoints(each) = PrePoints;
		each = each + 1;
	end	
	% BaseLinePoints is the number of data points in the PrePts section
	% Get the data in the PrePts section
	trial = condense(:,1:BaseLinePoints);
	BaseLineAverage = mean(trial,2);
	for i = 1:want
		condense(i,:) = condense(i,:) - BaseLineAverage(i);
	end
	condens = abs(condense)';
	histoplot = min(condens);
	playtime histogram
case 'threshold'
% 	makevisible = findobj(gcf,'tag','thresholdremove');
% 	set(makevisible,'visible','on');
 	wherefrom = gcbo;
 	mom = get(wherefrom,'userdata');
 	eval(mom);
	def = {'15'};
	info = inputdlg('What is the threshold?','Threshold',1,def);
	info = str2num(info{1});
	reject = find(histoplot>info);
	rejectplots = condense(reject,:);
	rejectlist = listepochs(reject);
	rejectnow = rejectlist(1);
	rejectind = 1;
	rejectnumlist = [];
gh0 = figure('Position',[256 334 600 400], ...
	'Tag','Fig2');
gh1 = axes('Parent',gh0, ...
	'Units','pixels', ...
	'CameraUpVector',[0 1 0], ...
	'Color',[1 1 1], ...
	'Position',[223 114 250 250], ...
	'Tag','Axes1', ...
	'XColor',[0 0 0], ...
	'YColor',[0 0 0], ...
	'ZColor',[0 0 0]);
gh2 = text('Parent',gh1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[0.4961240310077519 -0.06967213114754101 9.160254037844386], ...
	'Tag','Axes1Text4', ...
	'VerticalAlignment','cap');
set(get(gh2,'Parent'),'XLabel',gh2);
gh2 = text('Parent',gh1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[-0.08914728682170536 0.4959016393442622 9.160254037844386], ...
	'Rotation',90, ...
	'Tag','Axes1Text3', ...
	'VerticalAlignment','baseline');
set(get(gh2,'Parent'),'YLabel',gh2);
gh2 = text('Parent',gh1, ...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','right', ...
	'Position',[-0.8643410852713178 1.102459016393443 9.160254037844386], ...
	'Tag','Axes1Text2', ...
	'Visible','off');
set(get(gh2,'Parent'),'ZLabel',gh2);
gh2 = text('Parent',gh1, ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[0.4961240310077519 1.016393442622951 9.160254037844386], ...
	'Tag','Axes1Text1', ...
	'VerticalAlignment','bottom');
set(get(gh2,'Parent'),'Title',gh2);
gh1 = uicontrol('Parent',gh0, ...
	'BackGroundColor',[1 1 1], ...
	'Units','points', ...
	'Position',[10 10 30 350], ...
	'String',rejectlist, ...
	'Style','text', ...
	'Tag','rejectnums');
gh1 = uicontrol('Parent',gh0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...
	'Units','points', ...
	'Position',[60 160 100 30], ...
	'String','Current Epoch', ...
	'Style','text', ...
	'Tag','StaticText1');
gh1 = uicontrol('Parent',gh0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...
    'Units','points', ...
	'Position',[60 120 100 30], ...
	'String',rejectnow, ...
	'Style','text', ...
	'Tag','epochstr', ...
	'Userdata',rejectind);
gh1 = uicontrol('Parent',gh0, ...
	'Units','points', ...
	'Callback','playtime prevreject', ...
	'Position',[240 30 100 30], ...
	'String','Previous', ...
	'Tag','Pushbutton1');
gh1 = uicontrol('Parent',gh0, ...
	'Units','points', ...
	'Callback','playtime nextreject', ...
	'Position',[360 30 100 30], ...
	'String','Next', ...
	'Tag','Pushbutton2');
gh1 = uicontrol('Parent',gh0, ...
	'Units','points', ...
	'Callback','playtime removereject', ...
	'Position',[100 30 100 30], ...
	'String','Reject', ...
	'Tag','Pushbutton3');
gh1 = uicontrol('Parent',gh0, ...
	'BackGroundColor',[1 1 1], ...
	'Units','points', ...
	'Position',[485 10 30 350], ...
	'String',rejectnumlist, ...
	'Style','text', ...
	'Tag','newlist');
gh1 = uicontrol('Parent',gh0, ...
	'Units','points', ...
	'Callback','playtime histremove', ...
	'Position',[350 350 100 30], ...
	'String','Save', ...
	'Tag','thresholdremove');
plot(rejectplots(1,:));
case 'removereject'
	rejectnum = findobj(gcf,'tag','epochstr');
	rejectnum = get(rejectnum,'string');
	rejectnum = str2num(rejectnum);
	rejectnumlist = findobj(gcf,'tag','newlist');
	rejectnumlist = get(rejectnumlist,'string');
	rejectnumlist = str2num(rejectnumlist);
	if isempty(rejectnumlist)	
		rejectnumlist = rejectnum;
	else 
		rejind = length(rejectnumlist);
		rejectnumlist(rejind + 1) = rejectnum;
	end
	newlist = findobj(gcf,'tag','newlist');
	set(newlist,'string',rejectnumlist);
	rejeclength = length(rejectnumlist);
	newrejectlist(rejeclength) = find(listepochs==rejectnum);
	newrejectlist;
case 'histremove'
    newlist = findobj(gcf,'tag','newlist');
	epochstoremove = get(newlist,'string');
    epochstoremove = str2num(epochstoremove);
    close(gcf)
    epochdel = zeros(length(epochstoremove),1);
    for i = 1:length(epochstoremove)
        epochdel(i) = find(listepochs == epochstoremove(i));
    end
    playfig = findobj('tag','Fig1');
    figure(playfig);
    epochstr = findobj(gcf,'tag','epochstr');
	listrem = ones(size(listepochs));
	listrem(epochdel) = 0;
	listind = find(listrem);
	listnum = sum(listrem);
	tempgroup = groupepochs;	
	groupend = length(groupepochs);                      
	for i=1:listnum
		ind = listind(i);
		tempgroup2{1,i} = tempgroup{1,ind};
	end
	groupepochs = tempgroup2;
	% Must recreate cell with correct epoch #s
	groupepochs{2,1}=listepochs(listrem==1);
	cat(1,groupepochs{2,1});
	%  Now remove from gui
	listepochs = listepochs(listrem==1);
	set(epochstr,'Value',1);
	set(epochstr,'String',listepochs);
	set(epochstr,'value',newrejectlist);
case 'nextreject'
	rejectnow = findobj(gcf,'tag','rejectnums');
	rejectlist = get(rejectnow,'string');
	lastreject = findobj(gcf,'tag','epochstr');
	lastind = get(lastreject,'userdata');
	lastind = lastind + 1;
	rejectlength = size(rejectplots,1);
	if lastind > rejectlength
		lastind = 1;
	end
	cla
	lastnum = rejectlist(lastind,:);
	plot(rejectplots(lastind,:));
	set(lastreject,'userdata',lastind);
	set(lastreject,'string',lastnum);
case 'prevreject'
	rejectnow = findobj(gcf,'tag','rejectnums');
	rejectlist = get(rejectnow,'string');
	lastreject = findobj(gcf,'tag','epochstr');
	lastind = get(lastreject,'userdata');
	lastind = lastind - 1;
	rejectlength = size(rejectplots,1);
	if lastind < 1
		lastind = rejectlength;
	end
	cla
	lastnum = rejectlist(lastind,:);
	plot(rejectplots(lastind,:))
	set(lastreject,'userdata',lastind);
	set(lastreject,'string',lastnum);
case 'stdev'
	want = size(groupepochs,2);
	condense = cat(1,groupepochs{1,1:want});
	size(condense);
	condense = (condense)';
	stdev = std(condense);
	histoplot = stdev;
	playtime histogram	
case 'histogram'
	def = {'50'};
	bins = inputdlg('What is the bin size you prefer?','Threshold',1,def);
	bins = str2num(bins{1});
	figure;
	histoplot;
	hist(histoplot,bins);
	fortitle = gcbo;
	fortitle2 = get(fortitle,'label');
	title(fortitle2);
case 'return'
	delete(gcf)			
	delete(gcf)
	ITCFinishAnalysis(fileptr);
	% If this button was hit before the save button, then the user does not want this condition saved 
	% so we need to get rid of it.
	if saved == 0
		if structtype == 1
			CellInfo.EpochCondition(strindflag) = [];
			save(structname,'CellInfo')
		elseif structtype == 2
			CellInfo.FamilyCondition(strindflag) = [];
			save(structname,'CellInfo')
		end
	end
	step1 = strcmp(structname,'analyzethis');
	step2 = ~isfield(CellInfo,'comment');
	if strcmp(structname,'analyzethis') & ~isfield(CellInfo,'Comment')
		disp('new file')
		epochs('init',fp,plotpref,0);
	else	
		epochs('init',structname,plotpref,1);
	end
	%disp('Oops');
case 'finalsave'
	%disp('finalsave')
	% Create the Exclude matrix by comparing the two lists, the original and the list of 
	% only the epochs the user wants to keep.
	oldlis = findobj(gcf,'tag','oldlist');
	oldlis2 = findobj(gcf,'tag','oldlist2');
	ollist = get(oldlis,'string');
	ollist2 = get(oldlis2,'string');
	olist = str2num(ollist);
	olist2 = str2num(ollist2);
	oldlist = cat(1,olist,olist2)';
	lenlis = length(listepochs);
	lenold = length(oldlist);
	exclude = ones(size(oldlist));
	for a = 1:lenold
		for b = 1:lenlis
			if oldlist(a)==listepochs(b);
				exclude(a)=0;
			end
		end
	end
	if strindflag > 1
	  if structtype == 1
	    labDef = CellInfo.EpochCondition(strindflag - 1).Label;
	    decDef = CellInfo.EpochCondition(strindflag	- 1).DecimatePts;
	  elseif structtype == 2
	    labDef = CellInfo.FamilyCondition(strindflag - 1).Label;
	    decDef = CellInfo.FamilyCondition(strindflag - 1).DecimatePts;
	  end
	  decDef = {num2str(decDef)};
	else
	  labDef = {'Label'};
	  decDef = {'1'};
	end
    [epochLabel,Decimate] = getInput('init',labDef,decDef);
	Decimate = str2num(Decimate{1});
	%  This is stuff for CellInfo, if CellInfo does not already exist
    if newflag == 1
        allVars = getInput2;
        CellInfo.CellFile = fp;
        CellInfo.UserInfo = [];
        CellInfo.Comment = {allVars.comment};
        CellInfo.Label = {allVars.label}; 
        CellInfo.Rig = {allVars.rig};
        % This list must be the same as the list in getInput2
        cellTypes = {'rod' 'lCone' 'mCone' 'sCone' 'uvCone' 'ganglion'};
        CellInfo.CellType = cellTypes(allVars.type);
        numinput = length(allVars.allInputs);
        CellInfo.OutputScaleFactor = cell(numinput,1);
        CellInfo.OutputConfiguration = cell(numinput,1);
        CellInfo.NDFConfiguration = cell(numinput,1);
        info2 = allVars.allInputs;
        for i = 1:numinput
            CellInfo.OutputScaleFactor{i} = str2num(char(info2{i}{1}));
            CellInfo.OutputConfiguration{i} = info2{i}{2};
            CellInfo.NDFConfiguration{i} = str2num(char(info2{i}{3}));    
        end
    end
	if structtype == 1
		CellInfo.EpochCondition(strindflag).ExcludeEpochs = exclude;
		CellInfo.EpochCondition(strindflag).Label = epochLabel;
		CellInfo.EpochCondition(strindflag).DecimatePts = Decimate;
	elseif structtype == 2
		CellInfo.FamilyCondition(strindflag).ExcludeEpochs = exclude;
		CellInfo.FamilyCondition(strindflag).Label = epochLabel;
		CellInfo.FamilyCondition(strindflag).DecimatePts = Decimate;
	end	
	save(structname,'CellInfo')
	assignin('base','CellInfo',CellInfo)
	%familystructure
	exitnow = questdlg('Would you like to exit now?','Yes','No');
	if strcmp(exitnow,'Yes')
		playtime close
	else
		%  If we have saved this EpochCondition or FamilyCondition, then when we return to epochs 
		%  everything is fine, but if we decide to return to epochs without saving this 
		%  EpochCondition or FamilyCondition, we need to delete the last EpochCondition or 
		%  FamilyCondition, since these were already added by either epochs or familysearch.
		saved = 1;
	end
case 'close'
	filexist = which('analyzethis.mat');
	if ~isempty(filexist)
		delete(filexist)
	end
	if isempty(fileptr)
		delete(gcf)
		clear
	else
		ITCFinishAnalysis(fileptr);
		%disp('finished');
		delete(gcf)
		delete(gcf)
        close all
		clear all
	end
end

