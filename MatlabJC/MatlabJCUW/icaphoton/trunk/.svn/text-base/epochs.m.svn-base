function epochs(varargin)
%  
%  Summary:
%  This is the second of three interfaces, which together will select and
%  analyze a group of epochs in a selected file. Icaphoton selects the file,
%  and some preferences, epochs selects the epochs to analyze, and playtime
%  is for preliminary analysis and refined selection.  Epochs can be entered 
%  directly by typing epochs at the 
%  command line, with or without inputs.  It is recommended that this gui be 
%  entered by way of icaphoton, but this is not necessary, please see the m-file 
%  for details, particularly the comments at the very beginning of the file
%  (following the help comments) and after the persitent variables are declared.
%
%  Within epochs, begin by entering an epoch number or choosing a stimulus
%  type, and then pushing 'Push to search or clear'.  This button must be 
%  pushed whenever stim type is chosen to get back correct epochs.
%  The computer will give back the general stats on that epoch, then hit 'Push
%  to Plot' to get the full list of checkboxes, and to see the plot.
%
%  The computer uses the checkboxes marked with the values shown to search 
%  through the epochs.  The computer always uses, at the minimum, the
%  output type (ie. stimulus type) and the ampmode (voltage or current),
%  all other criteria are chosen by the user.  To find epochs that are
%  identical, push the button to set all checkboxes.  If both an epoch
%  number and a stim type are entered by the user, and the epoch number
%  does not actually correspond to the type of  stimuli, the computer will
%  return the actual stimuli type for that epoch in  the dialog box.  If
%  you flip thru the epochs using Next, Previous or Go to  Epoch #, epochs
%  will flip thru all epochs, unless you have a stim type chose,  in which
%  case it only flips thru that stim type.
%
%  Searches can be limited using the 'choose epochs to search from' button, this
%  list will remain in memory until changed, or until the epochs interface is 
%  exited.  Under 'epochs to search from' user can enter numbers, (ie. 0:99 or 
%  2,4,7:10), a file name containing numbers, or the word 'all' to change back to 
%  searching all epochs. Epoch numbers start with zero.
%

%
%  Special care must be used when your data file includes epochs with stimulus type 'sum'

%  epochs(y,fp,plotpref,cellstruct)
%  y is internal function (ie. y = 'init') to begin 
%		Inputs from icaphoton, or can be entered manually
%  fp			Variable received from Icaphoton that identifies data file,
%				default fp is given in file, 5th line of inner code
%	plotpref	Variable recieved from Icaphoton, represents whether the user
%				wants to look at the other segment or the stimulus or both
%				  Options are 
%				  1 Nothing specified
%				  2 Both segments shown, stim not shown
%				  3 Other segment not shown, stim shown
%				  4 Both segment and stim shown
%				default plotpref is 1
%	cellstruct is a boolean, if it is zero there is no cellinfo structure yet 
%           created for this data file, if it is one, user will be prompted for the name of
%           the file.

%  To add a new type of stimuli, first add the name of the stimuli to the 
%  string 'str' on line 16 and line  , and then add a new case with the same 
%  name.
%

%  Summary of persitent variables
%  	Variables received from Icaphoton (also sent to Playtime, except searnum)
%  	See above for summary.
%
%  	Variables used internally:
%  	fileptr  		Variable retrieved by ITCInitializeAnalusis, and used for all
%  						other ITC commands to indicate which data file active
%  	epochnumber		Variable that indicates the epoch that is currently active 
%		segnum  			Variable that indicates the segment that is currently active
%		NumEpochs		The total number of epochs in the active data file
%		NewOutputType	Number representing stimulus type, see case 'none_specified'
%							for more info
%		finalsearch		List of epochs used to search from, can be set either at the end 
%							of Icaphoton or at any time during epochs	
%
%  Rest of Variables sent to Playtime (through file analyzethis.mat)
%	datanum1	the values of the criteria used to select data ie. 
% 	datanum2	list of criteria used to select data, ie. 'tailpoints'
%	finalepochs	List of epochs chosen to analyze
%
%	Analyzethis.mat is a temporary file that is erased by Playtime when exiting.  If 
%   want to create a file to use in Playtime at a later date, use the save file button
%   in Playtime
%  Created January, 2001  MKMK
%
%	Changes in version 2.0:
%  Now there is a link to a function called familysearch, through the button
%  'Find Families', familysearch looks through the file, and pulls out all 
%  of the families, according to the information the user gives it.
%  There are two things one must remember when using this function, the epochs
%  interface must be loaded with an epoch that is in a family, and the epochs
%  must be on the plot, ie. Push to Plot button was pushed (this also gives
%  you the stimulus attributes). The user than precedes directly to
%  playtime, where epochs can be excluded.
%
%  2nd version completed June, 2001 MKMK
%
%  Minor bug fixed Feb, 2003 MKMK, error when the epoch number was missing or a 
%  number greater than # epochs given was fixed.  Still tells user what happened, but sets
%  the epoch # to zero, and proceeds normally.
%
%  Jul 2004  CHANGES MADE TO THIS VERSION OF EPOCHS:
%	The data structure is introduced earlier.  Now epochs initiates the structure
%  CellInfo before proceeding to playtime.  For more information about the structure
%  see help StructCellInfo.  
%    
%  Oct 2004 MKMK Changed ITCFinishAnalysis to use evalin in order to get around
%  wierd bug in Matlab version 7.  Also fixed bug when using icaphoton to
%  select epochs to choose from.  Did not ever open epochs before.
%   
%  Nov 2004 MKMK Fixed bugs with epochs window opening multiple times when
%  searching, and added amp to gaussiannoise as a search parameter.

persistent fp fileptr plotpref cellstruct newsearch epochnumber segnum NumEpochs NewOutputType datanum1 datanum2 finalsearch finalepochs structname CellInfo
varargin{:};
%  This first if statement was set as an easy way for me to enter epochs for 
%  troubleshooting and editing purposes, but means that I am always looking at 
%  the same data file, unless I change it here, in the file.
newsearch = 0;
if nargin<1
	y='init';
	fp = '/Users/maria/Data/original/testConverted';
	plotpref = 1;
	cellstruct = [];
	%  This is what is used internally to switch cases
elseif nargin<2
	y=varargin{1};
	%  plotpref default is one, so if nothing is specified, only one segment,
	%  and no stim is plotted.  This is only applicable if not entering through icaphoton
elseif nargin<3
	plotpref = 1;
elseif nargin<5
	[y,fp,plotpref,cellstruct] = deal(varargin{1:4});
	finalsearch = [];
else
    [y,fp,plotpref,cellstruct,newsearch] = deal(varargin{1:5});
end

switch y

case 'init'
	if isempty(fp)
		epochs close
		icaphoton	
	end
	if cellstruct == 1
		structname = fp;
		load(structname)
		fp = CellInfo.CellFile;
	else
		structname = [];
	end
	[fileptr, ecode] = ITCInitializeAnalysis(2e6, fp);
	if ecode~=0
		epochs error
	end
	%disp('initialized');
	[NumEpochs, ecode5] = ITCGetNumberEpochs(fileptr);
	if ecode5~=0
		epochs error
	end
str = {'none_specified' 'sinusoid' 'gaussiannoise' 'pulse' 'ramp',...
		'filewave' 'generic' 'sum' 'trigseries' 'segmentednoise' 'binarynoise' 'gclamp' 'circlepulse'};	
str2 = {fp 'Number of Epochs:' NumEpochs};		
load epochs
h0 = figure('Color',[0.8 0.8 0.8], ...
	'CloseRequestFcn','epochs close',...
	'Colormap',mat0, ...
    'MenuBar','none', ...
	'Position',[60 40 900 700], ...
	'Tag','Fig1');
h1 = axes('Parent',h0, ...
	'Units','pixels', ...
	'CameraUpVector',[0 1 0], ...
	'CameraUpVectorMode','manual', ...
	'Color',[1 1 1], ...
	'ColorOrder',mat1, ...
	'Nextplot','replace',...
	'Position',[500 250 380 380], ...
	'Tag','Axes1', ...
	'XColor',[0 0 0], ...
	'YColor',[0 0 0], ...
	'ZColor',[0 0 0]);
%  Title of Interface
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...
	'Units','points', ...
	'FontSize',24, ...
	'Position',[15 590 243 43], ...
	'String','Epochs', ...
	'Style','text', ...
	'Tag','Title');
%  This shows the user the current file
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
    'FontSize',12, ...
	'Position',[250 610 400 20], ...
	'String',str2(1), ...
	'Style','text', ...
	'Tag','filepath',...
	'Userdata',fp);
%  This is just a label
h1 = uicontrol('Parent',h0, ... 
	'Units','points', ...
    'BackgroundColor',[0.8 0.8 0.8], ...   
	'Position',[250 590 100 20], ...
	'String',str2(2), ...
	'Style','text', ...
	'Tag','filepath',...
	'Userdata',fp);
%  This shows the user the # epochs
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
    'BackgroundColor',[0.8 0.8 0.8], ...   
	'Position',[350 590 50 20], ...
	'String',str2(3), ...
	'Style','text', ...
	'Tag','filepath',...
	'Userdata',fp);
%  Instructional label
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Units','points', ...
	'FontSize',13, ...
	'Position',[25 560 140 35], ...
	'String','To begin, enter an epoch number', ...
	'Style','text', ...
	'Tag','StaticText3');
%  This is where the user enter an epoch #
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','epochs GetInfo',...
	'Position',[175 560 50 20], ...
	'String','',...
	'Style','edit', ...
	'Tag','epochnum');
%  To continue with the second set of buttons...
h1 = uicontrol('Parent',h0, ...
	'Callback','epochs GetInfo',...
	'Units','points', ...
	'FontSize',13, ...
	'Position',[25 510 200 35], ...
	'String','Push to search or clear', ...
	'Tag','getstats');
%  To enter numbers to search from...
h1 = uicontrol('Parent',h0, ...
	'Callback','epochs listnum',...
	'Units','points', ...
	'FontSize',10, ...
	'Position',[25 470 200 35], ...
	'String','Push to enter epoch numbers to search', ...
	'Tag','getstats');
%  To see numbers currently searching from...
h1 = uicontrol('Parent',h0, ...
	'Callback','epochs showlistnum',...
	'Units','points', ...
	'FontSize',10, ...
	'Position',[25 430 200 35], ...
	'String','Push to see epochs currently searching', ...
	'Tag','getstats');
%  This tells the user how many segments are in the epoch chosen
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Position',[350 570 50 20], ...
	'Style','text', ...
	'Tag','datasegments');
%  Label for previous text box
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Units','points', ...
	'Position',[250 570 100 20], ...
    'String','# of Data Segments:',...
	'Style','Text');
%  Now the user knows how many data segments, and can specify one
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','epochs segmentinfo',...
	'Position',[350 540 50 20], ...
	'Style','edit', ...
	'String','0',...
	'Tag','chooseseg');
%  Label for the previous user input box 
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Units','points', ...
	'Callback','',...
	'Position',[250 540 100 20], ...
	'Style','Text', ...
	'String','Which Segment #?');
%  Textbox where the computer returns the type of stim for chosen epoch
h1 = uicontrol('Parent',h0, ...
	'BackgroundColor',[0.9 0.9 0.9], ...
	'Units','points', ...
	'FontSize',13, ...
	'Position',[260 490 180 20], ...
	'Style','text', ...
	'Tag','Stimtype');
%  Label for previous text box, type of stim for epoch
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Units','points', ...
	'FontSize',13, ...
	'Position',[260 510 182 20], ...
	'String','Stim type for epoch shown', ...
	'Style','text', ...
	'Tag','StaticText1');
%  Popupmenu with choices of types of stimuli
h1 = uicontrol('Parent',h0, ...
	'FontSize',12, ...
	'Units','points', ...
	'Position',[260 420 182 25], ...
	'String',str, ...
	'Style','popupmenu', ...
	'Tag','Popupstim', ...
	'Value',1);
%  Label for popupmenu
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Units','points', ...
	'FontSize',13, ...
	'Position',[260 450 182 39], ...
	'String','To pick an epoch with a different type of stim', ...
	'Style','text', ...
	'Tag','StaticText1');
%  To continue with the third set of buttons...
h1 = uicontrol('Parent',h0, ...
	'Callback','epochs segmentinfo',...
	'Units','points', ...
	'FontSize',13, ...
	'Position',[260 390 182 20], ...
	'String','Push to Plot', ...
	'Tag','contseg');
%  This is an instruction label
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Units','points', ...
	'Position',[40 380 170 40], ...
	'String','Use checkboxes to determine search criteria', ...
	'Style','text', ...
	'Tag','StaticText9');
%  These are the text boxes that come up with values,
%  Starting with the ones common to all segments in a given epoch
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Position',[30 360 50 20], ...
	'Style','edit', ...
	'Tag','SampInterv');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Position',[30 330 50 20], ...
	'Style','edit', ...
	'Tag','EpochPts  ');
%  These are the corresponding checkboxes
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Position',[90 360 125 20], ...
	'String','Sampling Interval',...
    'Style','checkbox', ...
	'Tag','Checkbox1', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Position',[90 330 125 20], ...
    'String','# Data Pts in Segment',...
	'Style','checkbox', ...
	'Tag','epochsize1');
%  Textboxes that require segment #
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Position',[260 330 50 20], ...
	'Style','edit', ...
	'Tag','InputChan '); 
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Position',[260 300 50 20], ...
	'String','',...
	'Style','edit', ...
	'Tag','OutputChan');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Position',[260 270 50 20], ...
	'Style','edit', ...
	'Tag','AmpGain   ', ...
	'Visible','off');
%  And the corresponding checkboxes 
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Position',[320 330 125 20], ...
    'String','Input Channel',...
	'Style','checkbox', ...
	'Tag','Checkbox4');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Position',[320 300 125 20], ...
    'String','Output Channel',...
	'Style','checkbox', ...
	'Tag','nothing', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Position',[275 360 60 20], ...
	'String','voltage', ...
	'Style','radiobutton', ...
	'Tag','voltage');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Position',[365 360 60 20], ...
	'String','current', ...
	'Style','radiobutton', ...
	'Tag','current');
%  This is a label for the buttons created after stim is known
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
    'Units','points', ...
	'Position',[30 300 100 20], ...
	'String','Stimulus Attributes', ...
	'Style','Text', ...
	'Tag','stmlabel', ...
	'Visible','off');
%  To set all checkboxes at once
h1 = uicontrol('Parent',h0, ...
	'Callback','epochs setbuts',...
	'Units','points', ...
	'FontSize',13, ...
	'Position',[20 10 180 20], ...
	'String','Set all Checkboxes', ...
	'Tag','setboxes');
%  Clear all checkboxes
h1 = uicontrol('Parent',h0, ...
	'Callback','epochs clearboxes', ...
	'Units','points', ...
	'FontSize',13, ...
	'Position',[240 10 180 20], ...
	'String','Clear all Checkboxes', ...
	'Tag','Pushbutton3');
%  To flip through the epochs one at a time
%  Either back
h1 = uicontrol('Parent',h0, ...
	'Callback','epochs jump', ...
	'Units','points', ...
	'FontSize',13, ...
	'Position',[705 190 100 30], ...
	'String','Next', ...
	'Tag','Next');
%  or forward
h1 = uicontrol('Parent',h0, ...
	'Callback','epochs jump', ...
	'Units','points', ...
	'FontSize',13, ...
	'Position',[565 190 100 30], ...
	'String','Previous', ...
	'Tag','Previous');
%  This is where the comments will come up	
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.9 0.9 0.9], ...
	'Position',[460 65 380 80], ...
	'Style','text', ...
	'Tag','epochcom');
%  Label for previous textbox
h1 = uicontrol('Parent',h0, ...
    'BackgroundColor',[0.8 0.8 0.8], ...  
	'Units','points', ...
	'FontSize',13, ...
	'Position',[460 150 150 20], ...
	'String','Comment from Epoch', ...
	'Style','text', ...
	'Tag','StaticText11');
%  To separate the epochs into families
	h1 = uicontrol('Parent',h0, ...
	'Callback','epochs family', ...
	'Units','points', ...
	'FontSize',13, ...
	'Position',[460 10 120 40], ...
	'String','Find Families', ...
	'Tag','Pushbutton3');
%  To continue to next interface
h1 = uicontrol('Parent',h0, ...
	'Callback','epochs play', ...
	'Units','points', ...
	'FontSize',13, ...
	'Position',[590 10 120 40], ...
	'String','Analyze', ...
	'Tag','Pushbutton3');
%  To choose a different file, go back to icaphoton	
h1 = uicontrol('Parent',h0, ...
	'Callback','epochs return',...
	'Units','points', ...
	'FontSize',13, ...
	'Position',[720 10 120 40], ...
	'String','Icaphoton', ...
	'Tag','Pushbutton2');
	
%  This case sends you to familysearch
case 'family'
	delete(gcf)
	%finishAnalysis(fileptr);
    ITCFinishAnalysis(fileptr);
	%disp('finished - family')
    %structname
	if isempty(structname)
		save('analyzethis','fp','plotpref','segnum','structname','epochnumber');
	else
		if isfield(CellInfo,'FamilyCondition')
			famindex = length(CellInfo.FamilyCondition) + 1;
		else
			famindex = 1;
		end
		save('analyzethis','fp','plotpref','segnum','structname','epochnumber','famindex','CellInfo');
	end
	familysearch;
case 'listnum'
	%  This case allows user to specify which epochs to search
	def = {'[10:100,114:117]'};
	strangesearch = inputdlg('Enter epochs to search thru','Search Epochs',1,def);
	strangesearch = eval(strangesearch{1});
	if isempty(strangesearch) 
		finalsearch=[];
	else 
		finalsearch = strangesearch;
	end
    epochs showlistnum
case 'showlistnum'
	disp('Note: If finalsearch = [], then all the epochs have been selected')	
	disp('Epochs searching thru:')
    disp(finalsearch)
    newsearch 
    if newsearch == 1
        epochs('init')
    end
case 'GetInfo'
	%disp('First button working');
	%  Since this button may be used repeatedly, first we clear any buttons from
	%  last stim type or segnum
	%  First the Stim stuff
	genedit = findobj(gcf,'userdata','genedit');
	gencheck = findobj(gcf,'tag','gencheck');
	for y=1:length(genedit)
		delete(genedit(y));
	end
	for z=1:length(gencheck)
		delete(gencheck(z));
	end	
	% Defaults:
	SampInterv = findobj(gcf,'Tag','Checkbox1');
	OutputChan = findobj(gcf,'Tag','nothing');
	set(SampInterv,'value',1);
	set(OutputChan,'value',1);
	%  Clear all checkboxes
	%checkb = findobj(gcf,'Style','checkbox');
	%set(checkb,'value',0);
	%   
	%  
	voltage = findobj(gcf,'tag','voltage');
	current = findobj(gcf,'tag','current');
	set(voltage,'value',0);
	set(current,'value',0);
	intype = findobj(gcf,'tag','Inputtype ');
	inchan = findobj(gcf,'tag','InputChan ');
	outchan = findobj(gcf,'tag','OutputChan');
	set(intype,'String','');
	set(inchan,'String','');
	set(outchan,'String','');
	epochn = findobj(gcf,'tag','epochnum');
	epochnumber1 = get(epochn,'string');
	epochnumber = str2num(epochnumber1);
	%  Should I change this?  It automatically looks up epoch 0 if none is 
	%  specified.
	if isempty(epochnumber1) | epochnumber < 0
		epochnumber = 0;
		set(epochn,'string','0');
	elseif epochnumber > NumEpochs - 1
		disp('Epoch number is greater than number of epochs in file');
		epochnumber = 0;
		set(epochn,'string','0');
	% Try ending the if statement here.  What happens if it always go through the following 
	% code instead of only going through it if there is a "good" number in the epoch # field.
	% Not sure why I wasn't going through this code before for cases w/o a "good" #, but 
	% it seems to improve things, so I will leave it this way.
	end
	%else
	%  Get the comment for the epoch
	com = findobj(gcf,'tag','epochcom');
	set(com,'string','');
	%[EpochComment, ecode1] = ITCGetEpochComment(epochnumber,fileptr);
	%if ecode1 ~= 0
	%	epochs error
	%	else 
	%set(com,'string', EpochComment);
	%end
	%  Get the number of Segments in the epoch
	segm = findobj(gcf,'tag','datasegments');
	[NumSegments, ecode1] = ITCGetEpochSegments(epochnumber,fileptr);
	if ecode1 ~= 0
		epochs error
	else set(segm,'string', NumSegments);
	end
	%  Get the size of the epoch (# of points sampled in each segment)
	epochsi = findobj(gcf,'tag','EpochPts  ');
	[Epochsize, ecode2] = ITCGetEpochSize(epochnumber,fileptr);
	if ecode2 ~= 0
		epochs error
	else set(epochsi,'String', Epochsize);
	end
	%  Get the Sampling Interval (microsec)
	sampint = findobj(gcf,'tag','SampInterv');
	[SamplingInterval, ecode3] = ITCGetSamplingInterval(epochnumber,fileptr);
	if ecode3~= 0
		epochs error
	else set(sampint,'String',SamplingInterval);
	end
	
	epochs none_specified
	%end
	%  
	%  Now we need to get the info that is specific for the type of stimulus
	%  chosen.  First we go to Epochs none-specified to determine if the stim
	%  has been chosen or not.

	%  This is to find epochs with certain search parameters, specified by the 
	%  user with checkboxes.
case 'none_specified'
	%disp('2nd part of first button, none_specified, checking the stimtype');
	%  This is the textbox to input the real stim type to
	realtype = findobj(gcf,'tag','Stimtype');
	%  Get the epochnumber to get the actual stim type of the epoch # chosen
	%  Since there is an option of no stim chosen on the popupmenu, the lists 
	%  are not the same on the popupmenu and the ITCGetOutputType list, and we 
	%  need to add an index number to the value retrieved from OutputType in 
	%  order to compare them
	epochn = findobj(gcf,'tag','epochnum');
	epochnumber1 = get(epochn,'string');
	epochnumber = str2num(epochnumber1);
	chooseseg = findobj(gcf,'tag','chooseseg');
	segnum1 = get(chooseseg,'string');
	segnum = str2num(segnum1);
	[OutputType,ecode4] = ITCGetOutputType(epochnumber, segnum, 0, fileptr);
	if ecode4 ~= 0
		epochs error
	end
	%  Use the list from the popupmenu to get the right stim for the index #
	calling = findobj(gcf,'tag','Popupstim');
	stimval = get(calling, 'value');
	NewOutputType = OutputType + 1;
	%  Check to see if we got to this point from the user hitting next or
	%  previous.  If so, check to see what the popupmenu is showing.  If it is 
	%  none_specified, continue as normal on to the next epoch, if a stimtype 
	%  is specified, and there is a conflict with the next epoch, we will search
	%  for the next epoch of that specified stimtype.
	unk = gcbo;
	unk = get(unk,'tag');
	found = strcmp(unk,'Next');
	found2 = strcmp(unk,'Previous');
	if (found == 1 | found2 == 1) & stimval ~= 1 & stimval~=NewOutputType
		%disp('searching stim type');
		epochs simplesearch
	end
	%  Now we have to check if the stim type for the epoch chosen is the same
	%  stimtype that the user really wanted to look at.
	%  If the user has not chosen none_specified, and the value on the popupmenu
	%  and the value found by the epoch # listed are not the same than the 
	%  program searches for the first epoch of the stim type chosen in the
	%  popupmenu
	if stimval~=1 & stimval~=NewOutputType 
		%disp('Searching stim type');
		epochs simplesearch
	else
		%disp('No stimtype selected, use stimtype of epoch selected');
		str = get(calling, 'string');
		OutputName = str(NewOutputType);
		%		SINUSOID		1
	 	%		GAUSSIANNOISE	2
	 	%		PULSE			3
	 	%		RAMP			4	
	 	%		FILEWAVE		5	
	 	%		GENERIC			6	
		set(realtype,'String',OutputName);
		%com = findobj(gcf,'tag','chooseseg');
		other = findobj(gcf,'tag','contseg');
		OutputName2 = char(OutputName);
		OutputName3 = 'epochs ';
		OutputName4 = cat(2,OutputName3,OutputName2);
		%set(com,'callback',OutputName4); 	
		set(other,'callback',OutputName4);
		%disp('callback to 2nd button should be set');
	end
	
	
	if found == 1
		epochs segmentinfo
	end
	
	if found2 == 1 
		epochs segmentinfo
	end
	%  This case gets the info that is dependent on segment #, but is the same for
	%  all stim types.  This is the callback for the 2nd Press to Continue button
case 'segmentinfo'
	%  First position the buttons for the specific stim type
	%disp('2nd button pressed, segmentinfo')
	%segnum
	genedit = findobj(gcf,'Userdata','genedit');
	gencheck = findobj(gcf,'tag','gencheck');
	for y=1:length(genedit)
		set(genedit(y),'Position',[30 300-(y*30) 50 20]);
	end
	for z=1:length(gencheck)
		set(gencheck(z),'Position',[90 300-(z*30) 125 20]);
	end	
	% Label for above buttons
	stmlabel = findobj(gcf,'tag','stmlabel');
	set(stmlabel,'visible','on');
	%  The rest of the info
	epochn = findobj(gcf,'tag','epochnum');
	epochnumber1 = get(epochn,'string');
	epochnumber = str2num(epochnumber1);
	segm = findobj(gcf,'tag','datasegments');
	dataseg1 = get(segm,'string');
	dataseg = str2num(dataseg1);
	if segnum==dataseg
		disp('Segment #s start at zero');
		epochs error
	end
	voltage = findobj(gcf,'tag','voltage');
	current = findobj(gcf,'tag','current');
	[mode, ecode6] = ITCGetAmpMode(epochnumber, segnum, fileptr);
	if ecode6 ~=0
		epochs error
	end
	if mode==2
		set(voltage,'value',1);
        set(current,'value',0);
	elseif mode==6
		set(current,'value',1);
        set(voltage,'value',0);
	else
		set(voltage,'userdata',mode);
	end
	inchan = findobj(gcf,'tag','InputChan ');
	[InputChan, ecode7] = ITCGetInputChan(epochnumber, segnum, fileptr);
	if ecode7 ~= 0
		epochs error
	else
		set(inchan,'String',InputChan);
	end
	outchan = findobj(gcf,'tag','OutputChan');
	[OutputChan, ecode8] = ITCGetOutputChan(epochnumber, segnum, fileptr);
	if ecode8 ~= 0
		epochs error
	else
		set(outchan,'String',OutputChan);
	end
	[EpochData, ecode9] = ITCReadEpoch(epochnumber, segnum, fileptr);
	if ecode9 ~= 0
		epochs error
	end
	if plotpref == 1 
		%disp('business as usual');
		plot(EpochData)
	elseif plotpref == 2 | plotpref == 3
		epochs newfigure
	elseif plotpref == 4
		epochs newfigure2
	end
	
case 'newfigure'
	%disp('Making new graphs');
	segm = findobj(gcf,'tag','datasegments');
	dataseg1 = get(segm,'string');
	if ischar(dataseg1)
		dataseg = str2num(dataseg1);
	else
		dataseg = dataseg1;
	end
	if dataseg > 1 & segnum == 0 
		otherseg = 1;
	elseif dataseg > 1 & segnum == 1
		otherseg = 0;
	else otherseg = [];
	end
	[EpochData, ecode9] = ITCReadEpoch(epochnumber, segnum, fileptr);
		if ecode9 ~= 0
			epochs error
		end	
	if plotpref == 2
		if ~isempty(otherseg)
			[EpochData2, ecode10] = ITCReadEpoch(epochnumber, otherseg, fileptr);
			if ecode10 ~= 0
				epochs error
			end	
			epocplot = subplot('Position',[.53 .7 .45 .25]);
			plot(EpochData,'-k')
			stimplot = subplot('Position',[.53 .4 .45 .25]);
			plot(EpochData2,'-b')
			%disp('suposedly worked');
		elseif isempty(otherseg)
			delete(gca)
			delete(gca)
			h0 = gcf
			axes('Parent',h0, ...
			'Units','pixels', ...
			'CameraUpVector',[0 1 0], ...
			'CameraUpVectorMode','manual', ...
			'Color',[1 1 1], ...
			'Nextplot','replace',...
			'Position',[500 230 380 380], ...
			'Tag','Axes1', ...
			'XColor',[0 0 0], ...
			'YColor',[0 0 0], ...
			'ZColor',[0 0 0]);
			disp('No other segment for this epoch');
			plot(EpochData,'-k')
		end	
	elseif plotpref == 3 
		[StimData, ecode20] = ITCReadEpochStm(epochnumber, segnum, fileptr);
		if ecode20 ~= 0
			epochs error
		end
		epocplot = subplot('Position',[.53 .7 .45 .25]);
		plot(EpochData,'-k')
		stimplot = subplot('Position',[.53 .4 .45 .25]);
		plot(StimData,'-r')
		%disp('two graphs, stim and data');
	end
case 'newfigure2'
	segm = findobj(gcf,'tag','datasegments');
	dataseg1 = get(segm,'string');
	dataseg = str2num(dataseg1);
	if dataseg > 1 & segnum == 0 
		otherseg = 1;
	elseif dataseg > 1 & segnum == 1
		otherseg = 0;
	else otherseg = [];
	end
	epochnumber;
	segnum;
	otherseg;
	[EpochData, ecode9] = ITCReadEpoch(epochnumber, segnum, fileptr);
	if ecode9 ~= 0
		epochs error
	end
	[StimData, ecode20] = ITCReadEpochStm(epochnumber, segnum, fileptr);
	if ecode20 ~= 0
		epochs error
	end
	if ~isempty(otherseg)
		delete(gca)
		delete(gca)
		[EpochData2, ecode10] = ITCReadEpoch(epochnumber, otherseg, fileptr);
		if ecode10 ~= 0
			epochs error
		end			
		[StimData2, ecode21] = ITCReadEpochStm(epochnumber, otherseg, fileptr);
		if ecode21 ~= 0
			epochs error
		end
		%disp('should have 4 plots');
		epochplot = subplot('Position',[.54 .85 .45 .12]);
		plot(EpochData,'-k')
		title('Segment 0');
		stimsplot = subplot('Position',[.54 .68 .45 .12]);
		plot(StimData,'-r')
		title('Stimulus 0');
		epochplot2 = subplot('Position',[.54 .49 .45 .12]);
		plot(EpochData2,'-b')
		title('Segment 1');
		stimplot2 = subplot('Position',[.54 .33 .45 .12]);
		plot(StimData2,'-c')
		title('Stimulus 1');
	elseif isempty(otherseg)
		disp('No other segment for this epoch');
		epocplot = subplot('Position',[.55 .7 .45 .25]);
		plot(EpochData,'-k')
		stimplot = subplot('Position',[.55 .4 .45 .25]);
		plot(StimData,'-r')
	end		
case 'jump'
	jumpdir = gcbo;
	movedir = get(jumpdir,'Tag');
	movedirf = strncmp(movedir,'Previous',1);
	if movedirf == 1
		epochnumber = epochnumber - 1;
		epochn = findobj(gcf,'tag','epochnum');
		set(epochn,'string',epochnumber);
		epochs GetInfo
	elseif movedirf == 0
		epochnumber = epochnumber + 1;
		epochn = findobj(gcf,'tag','epochnum');
		set(epochn,'string',epochnumber);
		epochs GetInfo
	end
%  Stuff specific to certain stim types, get the info, set the button
case 'seed'
	[Seed, ecode10] =  ITCGetSeed(epochnumber, segnum, 0, fileptr);
	if ecode10 == 0	
		h0 = gcf;
		h1 = uicontrol('Parent',h0, ...
			'Units','points', ...
			'String',Seed, ...
			'Style','edit', ...	
			'Tag','Seed      ', ...
			'Userdata','genedit');
		h1 = uicontrol('Parent',h0, ...
			'Units','points', ...
	    	'String','Seed', ...
			'Style','checkbox', ...
			'Tag','gencheck');
	end	
case 'amp'
	[StimulusAmp, ecode11] = ITCGetStmAmp(epochnumber, segnum, 0, fileptr);
	if ecode11 == 0	
		h0 = gcf;
		h1 = uicontrol('Parent',h0, ...
			'Units','points', ...
			'String',StimulusAmp, ...
			'Style','edit', ...
			'Tag','StmAmp    ',...
			'Userdata','genedit');
		h1 = uicontrol('Parent',h0, ...
			'Units','points', ...
    		'String','Amp', ...
			'Style','checkbox', ...
			'Tag','gencheck', ...
			'Value',1);
	end	
case 'mean'
	[StmMean, ecode12] = ITCGetStmMean(epochnumber, segnum, 0, fileptr);
	if ecode12 == 0	
		h0 = gcf;
		h1 = uicontrol('Parent',h0, ...
			'Units','points', ...
			'String',StmMean, ...
			'Style','edit', ...
			'Tag','StmMean   ',...
			'Userdata','genedit');
		h1 = uicontrol('Parent',h0, ...
			'Units','points', ...
    		'String','Mean', ...
			'Style','checkbox', ...
			'Tag','gencheck');
	end	
case 'prepts'
	[PrePoints, ecode13] = ITCGetStmPrePts(epochnumber, segnum, 0, fileptr);
	if ecode13 == 0	
		h0 = gcf;
		h1 = uicontrol('Parent',h0, ...
			'Units','points', ...
			'String',PrePoints, ...
			'Style','edit', ...
			'Tag','PrePoints ', ...
			'Userdata','genedit');
		h1 = uicontrol('Parent',h0, ...
			'Units','points', ...
		    'String','Prepoints', ...
			'Style','checkbox', ...
			'Tag','gencheck');
	end	
case 'points'
	[Points, ecode14] = ITCGetStmPts(epochnumber, segnum, 0, fileptr);
	if ecode14 == 0	
		h0 = gcf;
		h1 = uicontrol('Parent',h0, ...
			'Units','points', ...
			'String',Points, ...
			'Style','edit', ...
			'Tag','StimDur   ', ...
			'Userdata','genedit');
		h1 = uicontrol('Parent',h0, ...
			'Units','points', ...
    		'String','Points', ...
			'Style','checkbox', ...
			'Tag','gencheck', ...
			'Value',1);
	end	
case 'tailpts'
	[TailPoints, ecode15] = ITCGetStmTailPts(epochnumber, segnum, 0, fileptr);
	if ecode15 == 0	
		h0 = gcf;
		h1 = uicontrol('Parent',h0, ...
			'Units','points', ...
			'String',TailPoints, ...
			'Style','edit', ...
			'Tag','TailPoints', ...
			'Userdata','genedit');
		h1 = uicontrol('Parent',h0, ...
			'Units','points', ...
    		'String','Tailpoints', ...
			'Style','checkbox', ...
			'Tag','gencheck');
	end	
case 'sinusoid'
	%  Clear any buttons from last stim type or segnum
	genedit = findobj(gcf,'userdata','genedit');
	gencheck = findobj(gcf,'tag','gencheck');
	for y=1:length(genedit)
		delete(genedit(y));
	end
	for z=1:length(gencheck)
		delete(gencheck(z));
	end	
	num = findobj(gcf,'tag','Stimtype');
	set(num,'String','sinusoid');
	epochs amp
	epochs mean
	epochs prepts
	epochs points
	epochs tailpts
	epochs segmentinfo
case 'gaussiannoise'
	%  Clear any buttons from last stim type or segnum
	genedit = findobj(gcf,'tag','genedit');
	gencheck = findobj(gcf,'tag','gencheck');
	for y=1:length(genedit)
		delete(genedit(y));
	end
	for z=1:length(gencheck)
		delete(gencheck(z));
	end	
	num = findobj(gcf,'tag','Stimtype');
	set(num,'String','gaussiannoise');
    epochs amp
	epochs seed
	epochs mean
	epochs prepts
	epochs points
	epochs tailpts
	epochs segmentinfo
case 'pulse'
	%  Clear any buttons from last stim type or segnum
	genedit = findobj(gcf,'userdata','genedit');
	gencheck = findobj(gcf,'tag','gencheck');
	for y=1:length(genedit)
		delete(genedit(y));
	end
	for z=1:length(gencheck)
		delete(gencheck(z));
	end	
	num = findobj(gcf,'tag','Stimtype');
	set(num,'String','pulse');
	epochs amp
	epochs mean
	epochs prepts
	epochs points
	epochs tailpts
	epochs segmentinfo
case 'ramp'
	%  Clear any buttons from last stim type or segnum
	genedit = findobj(gcf,'userdata','genedit');
	gencheck = findobj(gcf,'tag','gencheck');
	for y=1:length(genedit)
		delete(genedit(y));
	end
	for z=1:length(gencheck)
		delete(gencheck(z));
	end	
	num = findobj(gcf,'tag','Stimtype');
	set(num,'String','ramp');
	epochs amp
	epochs mean
	epochs prepts
	epochs points
	epochs tailpts
	epochs segmentinfo
case 'filewave'
	%  Clear any buttons from last stim type or segnum
	genedit = findobj(gcf,'userdata','genedit');
	gencheck = findobj(gcf,'tag','gencheck');
	for y=1:length(genedit)
		delete(genedit(y));
	end
	for z=1:length(gencheck)
		delete(gencheck(z));
	end	
	num = findobj(gcf,'tag','Stimtype');
	set(num,'String','filewave');
	epochs segmentinfo
case 'generic'
	%  Clear any buttons from last stim type or segnum
	genedit = findobj(gcf,'userdata','genedit');
	gencheck = findobj(gcf,'tag','gencheck');
	for y=1:length(genedit)
		delete(genedit(y));
	end
	for z=1:length(gencheck)
		delete(gencheck(z));
	end	
	num = findobj(gcf,'tag','Stimtype');
	set(num,'String','generic');
	epochs amp
	epochs mean
	epochs prepts
	epochs points
	epochs tailpts
	epochs segmentinfo	
	%  If the user does not specify a type of stim, then the push to cont. or clear
	%  button first sends them here.  The program then determines the stimulant 
	%  type and continues with the appropriate command
case 'simplesearch'
	%disp('Made it to simple')
	%  This simply searches for all epochs of the given output type (ie. stimulus)
	%  Used during jump to determine if output changes while flipping thru epochs 
	%  This is necessary, because if there is a stim type chosen, the interface 
	%  is designed to only show those epochs of that stim type
	%  Check to see how we got here
	unk = gcbo;
	unk = get(unk,'tag');
	found = strcmp(unk,'Next');
	found2 = strcmp(unk,'Previous');	
	if found == 1 
		startepoch = epochnumber;
	else startepoch = 0;
	end	
	newepochs = ones(1,NumEpochs);
	calling = findobj(gcf,'tag','Popupstim');
	stimval = get(calling, 'value');
	if stimval == 1
		%disp('Why did this happen?');
		%  If the user entered none_specified than search on the type of stim
		%  for the epoch # listed.	
	else Outputtype = stimval - 1;
	end	
	segnum;
	if found2 == 1
		for searchepoch = NumEpochs-1:-1:0
			[OutputType,ecode4] = ITCGetOutputType(searchepoch, segnum, 0, fileptr);
			if OutputType ~= Outputtype & newepochs(searchepoch+1) == 1
				newepochs(searchepoch+1) = 0;
			end
		end
	else	
		for searchepoch = startepoch:NumEpochs-1
			[OutputType,ecode4] = ITCGetOutputType(searchepoch, segnum, 0, fileptr);
			if OutputType ~= Outputtype & newepochs(searchepoch+1) == 1
				newepochs(searchepoch+1) = 0;
			end
		end
	end
		typepoch = find(newepochs==1);
		typeepoch = typepoch - 1;
		if isempty(typeepoch)
			disp('There are no (more) epochs in this segment with this type of stimulus');	
			epochs error
		end
		epochnum = findobj(gcf,'Tag','epochnum');
		set(epochnum,'String',typeepoch(1));
		str = get(calling, 'string');
		NewOutputType = Outputtype + 1;
		OutputName = str(NewOutputType);
			
		%		SINUSOID		1
 		%		GAUSSIANNOISE	2
 		%		PULSE			3
 		%		RAMP			4
 		%		FILEWAVE		5	
 		%		GENERIC			6	
		realtype = findobj(gcf,'tag','Stimtype');
		set(realtype,'String',OutputName);
		com = findobj(gcf,'tag','chooseseg');
		other = findobj(gcf,'tag','contseg');
		OutputName2 = char(OutputName);
		OutputName3 = 'epochs ';
		OutputName4 = cat(2,OutputName3,OutputName2);
		%set(com,'callback',OutputName4); 	
		set(other,'callback',OutputName4);
		%disp('Went through search for type w/o mode');
		epochs GetInfo
case 'search'
	%disp('searching');
	% For ampmode value 1 means current, value 0 means voltage, vise-versa for
	% ampmode3.  All searches look at ampmode and stim type
	%  LOOK AT STIM TYPE FOR REGULAR SEARCHES!!!
	ampmode = findobj(gcf,'tag','current');
	ampmode2 = get(ampmode,'value');
	ampmode3 = findobj(gcf,'tag','voltage');
	ampmode4 = get(ampmode3,'value');
	ampmode5 = get(ampmode3,'userdata');
	if ampmode2 == 1 & ampmode4 == 0
		ampmodef = 6;
	elseif ampmode2 == 0 & ampmode4 == 1
		ampmodef = 2;
	else 
		disp('Warning:  Ampmode is neither current nor voltage')
		ampmodef = ampmode5
	end
	ampmodef;
	%disp('Searching for like epochs, have mode already');
	%  Now for the rest of the stuff, first find which boxes are checked
	searchbox = findobj(gcf,'Style','checkbox','value',1);
	ldat = length(searchbox);
	%  Create an empty matrix, the length corresponds to the # of boxes checked,
	%  the width to a position vector
	searchpos = zeros(ldat,4);
	s1 = zeros(ldat,1);
	s2 = zeros(ldat,1);
	%  Now find the label and position for each search parameter.  The position
	%  of the label and the input value are always the same in relation to each 
	%  other, so we can use one to get the other.
	for x = 1:ldat;
		data = get(searchbox(x),'string');
		foo = get(searchbox(x),'Position');	
		searchpos(x,:) = foo;
	end
	%  The position columns that matter are 1 and 2
	s1 = searchpos(:,1);
	s2 = searchpos(:,2);
	sdat = length(s1);
	datanum = zeros(sdat,1);
	datanum1 = zeros(sdat,1);
	datanum2 = zeros(sdat,10);
	datanum3 = zeros(sdat,1);
	for y = 1:sdat
		foo = findobj(gcf,'Position',[s1(y)-60 s2(y) 50 20]);
		datanum(y,:) = foo;
		foo1 = str2num(get(datanum(y),'string'));
		datanum1(y,:) = foo1;
		foo2 = (get(datanum(y),'tag'));	
		datanum2(y,:) = foo2;
	end
	datanum2 = char(datanum2);
	datanum2 = cellstr(datanum2);
	datanum1;
	datasize = length(datanum2);
	newepochs = zeros(1,NumEpochs);
	% Newepochs is a matrix of zeros as big as the data set.  When we look at the mode, we
	% change all zeros to ones that have the same mode as the reference (remember that we 
	% use zero for ITC commands, but there is no zero index in a matlab matrix, so the 
	% corresponding index will be the searchepoch number + 1).  From then on, if 
	% a parameter value is the same as the reference, then we leave it alone.  If the 
	% parameter value is not the same, (and the corresponding Newepochs number is one, 
	% although this check is really unnecessary) then we change the Newepochs number 
	% corresponding with this epoch to zero.
	if isempty(finalsearch)
		FinalSearch = 0:NumEpochs-1;
	else FinalSearch = finalsearch; 
	end
	%for searchepoch = 0:NumEpochs-1
	for searchepoch = FinalSearch
		[mode, ecode6] = ITCGetAmpMode(searchepoch, segnum, fileptr);
		if mode == ampmodef
			newepochs(searchepoch+1) = 1;
		end
	end	
	Outputtype = NewOutputType - 1;
	%for searchepoch = 0:NumEpochs-1
	for searchepoch = FinalSearch	
		searchepoch;
	[OutputType,ecode4] = ITCGetOutputType(searchepoch, segnum, 0, fileptr);
		if OutputType ~= Outputtype & newepochs(searchepoch+1) == 1
			newepochs(searchepoch+1) = 0;
		end
	end
	for k = 1:datasize
		if isempty(datanum2)
            break
		end
		tail(k) = strcmp(datanum2{k},'TailPoints');
		if tail(k)==1
			Tailpoints = datanum1(k);
			%for searchepoch = 0:NumEpochs - 1
			for searchepoch = FinalSearch
				[TailPoints, ecode15] = ITCGetStmTailPts(searchepoch, segnum, 0, fileptr);	
				if ecode15~=0
					disp('Problem with TailPoints');
				break
				end
				if TailPoints ~= Tailpoints & newepochs(searchepoch+1) == 1
					newepochs(searchepoch+1) = 0;
				end	
			end
		end
		points(k) = strcmp(datanum2{k},'StimDur');
		if points(k)==1
			points = datanum1(k);
			%for searchepoch = 0:NumEpochs - 1;
			for searchepoch = FinalSearch
				[Points, ecode14] = ITCGetStmPts(searchepoch, segnum, 0, fileptr);
				if ecode14~=0
					disp('Problem with Points');
				break
				end
				if Points ~= points & newepochs(searchepoch+1) == 1
					newepochs(searchepoch+1) = 0;
				end
			end
		end
		pre(k) = strcmp(datanum2{k},'Prepoints');
		if pre(k)==1
			Prepoints = datanum1(k);
			%for searchepoch = 0:NumEpochs - 1
			for searchepoch = FinalSearch
				[PrePoints, ecode13] = ITCGetStmPrePts(searchepoch, segnum, 0, fileptr);	
				if ecode13~=0
					disp('Problem with PrePoints');
				break
				end
				if PrePoints ~= Prepoints & newepochs(searchepoch+1) == 1
					newepochs(searchepoch+1) = 0;		
				end
			end	
		end
		mean(k) = strcmp(datanum2{k},'StmMean');
		if mean(k)==1
			Stmmean = datanum1(k);
			%for searchepoch = 0:NumEpochs - 1;
			for searchepoch = FinalSearch
				[StmMean, ecode12] = ITCGetStmMean(searchepoch, segnum, 0, fileptr);
				if ecode12~=0
					disp('Problem with StmMean');
				break
				end
				if StmMean ~= Stmmean & newepochs(searchepoch+1) == 1
					newepochs(searchepoch+1) = 0;
				end
			end	
		end
		amp(k) = strcmp(datanum2{k},'StmAmp');
		if amp(k)==1
			Stmamp = datanum1(k);
			%for searchepoch = 0:NumEpochs - 1;
			for searchepoch = FinalSearch
				[StmAmp, ecode11] = ITCGetStmAmp(searchepoch, segnum, 0, fileptr);
				if ecode11~=0
					disp('Problem with amp');
				break
				end
				if StmAmp ~= Stmamp & newepochs(searchepoch+1) == 1
					newepochs(searchepoch+1) = 0;
				end
			end	
		end
		seed(k) = strcmp(datanum2{k},'Seed');
		if seed(k)==1
			Nseed = datanum1(k);
			%for searchepoch = 0:NumEpochs - 1;
			for searchepoch = FinalSearch
                [OutputType,ecode4] = ITCGetOutputType(searchepoch, segnum, 0, fileptr);
                if ecode4 ~= 0
                    epochs error
                end
                if OutputType == 2
    				[Seed, ecode10] =  ITCGetSeed(searchepoch, segnum, 0, fileptr);
        			if ecode10~=0
            			disp('Problem with seed');
                	break
                    end
                    if Seed ~= Nseed & newepochs(searchepoch+1) == 1
                        newepochs(searchepoch+1) = 0;
                    end
                end
			end		
		end
		outchan(k) = strcmp(datanum2{k},'OutputChan');
		if outchan(k)==1
			Outputchan = datanum1(k);
			%for searchepoch = 0:NumEpochs - 1;
			for searchepoch = FinalSearch	
				[OutputChan, ecode8] = ITCGetOutputChan(searchepoch, segnum, fileptr);
				if OutputChan ~= Outputchan & newepochs(searchepoch+1) == 1
					newepochs(searchepoch+1) = 0;
				end
			end		
		end
		inchan(k) = strcmp(datanum2{k},'InputChan');
		if inchan(k)==1
			Inputchan = datanum1(k);
			%for searchepoch = 0:NumEpochs - 1;
			for searchepoch = FinalSearch
				[InputChan, ecode7] = ITCGetInputChan(searchepoch, segnum, fileptr);
				if InputChan ~= Inputchan & newepochs(searchepoch+1) == 1
					newepochs(searchepoch+1) = 0;
				end
			end		
		end
		eposiz(k) = strcmp(datanum2{k},'EpochPts');
		if eposiz(k)==1
			Epochsize = datanum1(k);
			%for searchepoch = 0:NumEpochs - 1;
			for searchepoch = FinalSearch
				[EpochSize, ecode2] = ITCGetEpochSize(searchepoch,fileptr);
				if EpochSize ~= Epochsize & newepochs(searchepoch+1) == 1
					newepochs(searchepoch+1) = 0;
				end
			end	
		end
		samint(k) = strcmp(datanum2{k},'SampInterv');
		if samint(k)==1
			Sampinterv = datanum1(k);
			%for searchepoch = 0:NumEpochs - 1;
			for searchepoch = FinalSearch
				[SamplingInterval, ecode3] = ITCGetSamplingInterval(searchepoch,fileptr);
				if SamplingInterval ~= Sampinterv & newepochs(searchepoch+1) == 1
					newepochs(searchepoch+1) = 0;
				end
			end	
		end
	end
	%  finalepoch is the matlab matrix, which starts at 1, finalepochs is 
	%  the actual epoch #s, which start at zero
	FinalSearch;
	finalepoch = find(newepochs==1);
	finalepochs = finalepoch - 1;
	if isempty(finalepochs)
		error('Warning: No epochs chosen')
	end
case 'clearboxes'
	boxes = findobj(gcf,'style','checkbox');
	set(boxes,'value',0);
case 'play'
	epochs search
	delete(gcf)
	delete(gcf)
	finalength = length(finalepochs);
	SingleCondition.EpochTime = zeros(1,finalength);
	for i = 1: finalength
		ind = finalepochs(i);
		[EpochTime, ecode] = ITCGetEpochTime(ind, fileptr);
		SingleCondition.EpochTime(i) = EpochTime;
	end
	%finalepochs
	SingleCondition.EpochNumbers = finalepochs;
	SingleCondition.SearchCrit = datanum2;
	SingleCondition.SearchPara = datanum1;
	SingleCondition.PlotPref = plotpref;
	SingleCondition.SegNum = segnum;
	SingleCondition.ExcludeEpochs = [];
	SingleCondition.Label = [];
	SingleCondition.ScaleFactorIndex = [];
	SingleCondition.AverageResponse = [];
	SingleCondition.VarianceResponse = [];
	SingleCondition.DecimatePts = 1;
	SingleCondition.UserInfo = [];
	structname;
	if isempty(structname)
		CellInfo.CellFile = fp;
		CellInfo.EpochCondition = deal(SingleCondition);
		newflag = 0;
	else
		newflag = 1;
		%  If CellInfo already exists, we need to see if there are already EpochConditions 
		%  stored.  Could be that CellInfo already exists because there are FamilyConditions.
		%  If there are EpochConditions, we add to the list, if not, we add the field
		if isfield(CellInfo,'EpochCondition')
			strind = length(CellInfo.EpochCondition) + 1;
			CellInfo.EpochCondition(1,strind) = deal(SingleCondition);
			%assignin('base','CellInfo',CellInfo)
			save(structname,'CellInfo')
		else
			strind = 1;
			CellInfo.EpochCondition = deal(SingleCondition);
			save(structname,'CellInfo')
		end
	end	
    %finishAnalysis(fileptr)
    ITCFinishAnalysis(fileptr);
	%disp('finished');
	%savnum = input('Would you like to save the epoch numbers in a separate file? y or n ','s');
	%if strcmp(savnum,'y')
	%	savfil = input('What would you like to call the file? \n','s');
	%	save(savfil,'finalepochs')
	%end	
	save('analyzethis','CellInfo');
	if newflag == 1
		playtime('init',structname,1,strind)
	else
		playtime('init','analyzethis',1)
	end	
case 'setbuts'
	searchbox = findobj(gcf,'Style','checkbox');
	set(searchbox,'value',1);	
case 'return'
	delete(gcf)
	ITCFinishAnalysis(fileptr);
	disp('returned');
	icaphoton('init')
case 'error'
	com = findobj(gcf,'tag','epochcom');
	set(com,'String','There is a problem, please see command window for details'); 
case 'close'
	%disp('Leaving Epochs')
	filexist = which('analyzethis.mat');
	if ~isempty(filexist)
		delete(filexist)
		%disp('removed')
	end
	if isempty(fileptr)
		delete(gcf)
		clear
	else
		ITCFinishAnalysis(fileptr);
		%disp('epochs closed');
		delete(gcf)
		clear
	end
end

		


