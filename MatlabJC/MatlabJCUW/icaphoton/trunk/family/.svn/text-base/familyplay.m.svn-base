function familyplay(a)
% function familyplay(a)
% Load CellInfo into the workspace, then run familyplay with no inputs
% Need to have epoch data loaded into the structure. The function uses a
% structure that has been created by epochs, using the find families
% option, and playtime, and currently loaded into the workspace.  After
% each plotting, there are two variables sent to the workspace called
% xvalue and yvalue.  These are the x and y values that were plotted.
% They must be saved as a different variable or to a file before anymore
% plotting occurs, or they will be overwritten.


% 
% Let's use the terminology cue to identify the type to sort on, ie. led,
% temp, etc., and the term cuevalue to identify the number that
% identifies which value of that cue, ie. 2 is the cuevalue that
% corresponds to the red led, which is, for example, cue 1.

% Created June, 2001 MKMK

persistent CellInfo EpochNumbers EpochData famnumber FamilyFlag FamilyStep ...
	 somerealfamilies FamilyLabel yvalue yaxislabel xaxislabel flag ...
	 cueave stepave tempfamilies cuerow xtime timeplot
%somerealnumbs 
if nargin<1
  a = 'init';
end

switch a
  
 case 'init'
  % Get the structure, get various structures from it to save typing.
  CellInfo = evalin('base','CellInfo');
  famnum = input('Which family structure would you like to analyze (#): ');
  famnumber = famnum;
  famcell = CellInfo.FamilyCondition(famnum);
  EpochNumbers = famcell.EpochNumbers;
  EpochData = CellInfo.EpochData.Data;
  FamilyStep = famcell.FamilyStep;
  FamilyFlag = famcell.FamilyFlag;
  %FamilyIndex = famcell.FamilyIndex;
  FamilyCues = famcell.FamilyCues;
  FamilyLabel = famcell.Label;
  FamilyCueGuide = famcell.FamilyCueGuide;
  FamilyStepGuide = famcell.FamilyStepGuide;
  % Sampling interval is in microseconds
  xtime = FindSearchPara(famcell,'SampInterv');
  xtime = xtime*1e-6;
  if isempty(xtime)
	 disp('Sampling Interval could not be found, using default');
	 xtime = 0.01;
  end
  % If the data has been decimated, need to find the adjusted sampling interval
  decimate = famcell.DecimatePts;
  xtime = xtime*decimate;
  famcue = 0;
  % Set up figure to get the wants and desires of the user
  h0 = figure('CloseRequestFcn','familyplay close',...
				  'Position',[30 30 600 850], ...
				  'Tag','Fig1');	
  %  Title of Interface
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'FontSize',30, ...
					  'Position',[200 650 170 40], ...
					  'String','Familyplay', ...
					  'Style','text');
  % Allows user to limit the families by epoch numbers
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Position',[10 630 120 20], ...
					  'String','Graph epochs from:', ...
					  'Style','text');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Position',[130 630 40 20], ...
					  'String','0', ...
					  'Tag','startrange');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Position',[170 630 30 20], ...
					  'String','to:', ...
					  'Style','text');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Position',[200 630 40 20], ...
					  'String','end', ...
					  'Style','edit', ...
					  'Tag','endrange');
  % Label
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'FontSize',20, ...
					  'Position',[10 600 100 30], ...
					  'String','Sort by:', ...
					  'Style','text');
  % choice for selection
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton4', ...
					  'Position',[10 580 170 20], ...
					  'String','Select families including ONLY:', ...
					  'Style','radiobutton', ...
					  'Tag','only', ...
					  'Userdata','buttongroup4');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton4', ...
					  'Position',[10 560 210 20], ...
					  'String','Select families including a minimum of:', ...
					  'Style','radiobutton', ...
					  'Tag','minim', ...
					  'Userdata','buttongroup4', ...
					  'Value',1);
  stepcontrol = uicontrol('Parent',h0, ...
								  'Units','points', ...
								  'Position',[25 540 200 20], ...
								  'String','steps same as family incl. epoch #', ...
								  'Style','checkbox', ...
								  'Tag','steplistepoch');
  stepcontrol = uicontrol('Parent',h0, ...
								  'Units','points', ...
								  'Callback','familyplay setbutton', ...
								  'Position',[220 540 50 20], ...
								  'String','', ...
								  'Style','edit', ...
								  'Tag','stepepoch');
  stepcontrol = uicontrol('Parent',h0, ...
								  'Units','points', ...
								  'Position',[25 520 150 20], ...
								  'String','all', ...
								  'Style','checkbox', ...
								  'Tag','steplistall');
  stepcontrol = uicontrol('Parent',h0, ...
								  'Units','points', ...
								  'Position',[150 520 150 20], ...
								  'String','average by step', ...
								  'Style','checkbox', ...
								  'Tag','stepave');
  % This looks in CellInfo, and gets the steps that define the family,
  % and makes a checkbox for each one
  stepvtest = finite(FamilyStep);
  stepval = unique(FamilyStep(stepvtest));
  numstep = length(stepval);
  for i = 1:numstep
	 if stepval(i)==0
		i = i+1;
	 end
	 posvector = 520-(i*20);
	 stepcontrol(i) = uicontrol('Parent',h0, ...
										 'Units','points', ...
										 'Position',[25 posvector 150 20], ...
										 'String',stepval(i), ...
										 'Style','checkbox', ...
										 'Tag','steplist');
  end
  % Now do the same thing for cues, only there are not only several
  % values for a cue (cuevalue), but also can be more than one cue.  Also
  % more complicated because we make checkboxes allowing the user to
  % do the following: the option to average all families having this value;
  % assign one cue as a title (ie. all lines on plot have this value);
  % or one/several as the legend (different values for different lines
  % plotted).
  
  % i is kept in userdata, so we know which row the cuevalue is in, in
  % the case that there are more than one cues (numcue is longer than one)  
  posvector = posvector - 40;
  numcue = length(FamilyCueGuide);
  for i = 1:numcue
	 cuename = FamilyCueGuide{i};
	 cuevalue = unique(FamilyCues(i,:));
	 cuelength = length(cuevalue);
	 % Title of the cues
	 cuecontrol = uicontrol('Parent',h0, ...
									'Units','points', ...
									'FontSize',20, ...
									'Position',[5  posvector 200 30], ...
									'String',cuename, ...
									'Style','text', ...
									'Tag','cuelist');
	 valvector = posvector;
	 avestr = ['average by ' cuename];
	 % Average option
	 avecuecontrol = uicontrol('Parent',h0, ...
										'Units','points', ...
										'Position',[150 valvector-20 150 20], ...
										'String',avestr, ...
										'Style','checkbox', ...
										'Tag','cueave', ...
										'Userdata',i);
	 % Title option
	 titlecontrol = uicontrol('Parent',h0, ...
									  'Units','points', ...
									  'Position',[150 valvector-40 150 20], ...
									  'String','Title', ...
									  'Style','checkbox', ...
									  'Tag','plottitle', ...
									  'Userdata',i);
	 % Legend option
	 legendcontrol = uicontrol('Parent',h0, ...
										'Units','points', ...
										'Position',[150 valvector-60 150 20], ...
										'String','Legend', ...
										'Style','checkbox', ...
										'Tag','plotlegend', ...
										'Userdata',i);
	 % Setting cue checkboxes									  
	 for j = 1:cuelength
		cuevalcontrol = uicontrol('Parent',h0, ...
										  'Units','points', ...
										  'Position',[25 valvector-20*j 100 20], ...
										  'String',cuevalue(j), ...
										  'Style','checkbox', ...
										  'Tag','cuevallist', ...
										  'Userdata',i);
		  posvector = valvector-(20*(j+2));
	 end
  end
  % A bunch of choices for the actual plots
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'FontSize',20, ...
					  'Position',[300 610 200 30], ...
					  'String','Graphing options:', ...
					  'Style','text');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton', ...
					  'Position',[300 590 170 20], ...
					  'String','plot', ...
					  'Style','radiobutton', ...
					  'Tag','simplot', ...
					  'Userdata','buttongroup');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton', ...
					  'Position',[300 570 210 20], ...
					  'String','histogram (uses y-values only)', ...
					  'Style','radiobutton', ...
					  'Tag','histo', ...
					  'Userdata','buttongroup');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'FontSize',20, ...
					  'Position',[300 520 100 30], ...
					  'String','X Axis:', ...
					  'Style','text');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton2', ...
					  'Position',[300 500 200 20], ...
					  'String','Flash Strength ("I", Step parameter)', ...
					  'Style','checkbox', ...
					  'Tag','stepx', ...
					  'Userdata','buttongroup2');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton2', ...
					  'Position',[300 480 100 20], ...
					  'String','time', ...
					  'Style','checkbox', ...
					  'Tag','timex', ...
					  'Userdata','buttongroup2');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Position',[300 440 100 20], ...
					  'String','log scale', ...
					  'Style','checkbox', ...
					  'Tag','logx');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'FontSize',20, ...
					  'Position',[300 390 100 30], ...
					  'String','Y Axis:', ...
					  'Style','text');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton3', ...
					  'Position',[300 370 200 20], ...
					  'String','Max. Response,"R"', ...
					  'Style','checkbox', ...
					  'Tag','maxy', ...
					  'Userdata','buttongroup3');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton3', ...
					  'Position',[300 350 200 20], ...
					  'String','Area Response', ...
					  'Style','checkbox', ...
					  'Tag','areay', ...
					  'Userdata','buttongroup3');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton3', ...
					  'Position',[300 330 200 20], ...
					  'String','Normalized Max.,"R/Rmax"', ...
					  'Style','checkbox', ...
					  'Tag','normmaxy', ...
					  'Userdata','buttongroup3');
  h1 = uicontrol('Parent',h0, ... 
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton3', ...
					  'Position',[300 310 200 20], ...
					  'String','Normalized Area', ...
					  'Style','checkbox', ...
					  'Tag','normareay', ...
					  'Userdata','buttongroup3');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton3', ...
					  'Position',[300 290 100 20], ...
					  'String','Response', ...
					  'Style','checkbox', ...
					  'Tag','response', ...
					  'Userdata','buttongroup3');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Position',[300 250 100 20], ...
					  'String','log scale', ...
					  'Style','checkbox', ...
					  'Tag','logy');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'FontSize',20, ...
					  'Position',[300 190 200 30], ...
					  'String','Normalize Options:', ...
					  'Style','text');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton5', ...
					  'Position',[300 170 200 20], ...
					  'String','All epochs this type', ...
					  'Style','checkbox', ...
					  'Tag','allepoch', ...
					  'Userdata','buttongroup5');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton5', ...
					  'Position',[300 150 200 20], ...
					  'String','All epochs in Title', ...
					  'Style','checkbox', ...
					  'Tag','titlepoch', ...
					  'Userdata','buttongroup5');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay exclusivebutton5', ...
					  'Position',[300 130 200 20], ...
					  'String','Only epochs in family group', ...
					  'Style','checkbox', ...
					  'Tag','famgp', ...
					  'Userdata','buttongroup5', ...
					  'Value',1);
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay prepare', ...
					  'Position',[200 100 150 30], ...
					  'String','Set Options');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'Callback','familyplay comment', ...
					  'Position',[380 100 150 20], ...
					  'String','Get Label');
  h1 = uicontrol('Parent',h0, ...
					  'Units','points', ...
					  'FontSize',10, ...
					  'Position',[380 130 150 40], ...
					  'String','', ...
					  'Style','text', ...
					  'Tag','comments', ...
					  'Visible','off');  
 case 'setbutton'
  % To set the checkbox if something is written in the corresponding edit
  % box, ie. make life easier so you don't have to write in the box AND
  % check the corresponding box
  findbutton = findobj(gcf,'Tag','steplistepoch');
  set(findbutton,'Value',1);
  % These cases ensure that the user doesn't try to check more than one
  % box where he/she isn't suppose to, ie. only one thing on the x axis 
 case 'exclusivebutton'
  h = findobj(gcbf,'userdata','buttongroup');
  set(gcbo,'value',1);
  set(h(h~=gcbo),'value',0);
 case 'exclusivebutton2'
  h = findobj(gcbf,'userdata','buttongroup2');
  set(gcbo,'value',1);
  set(h(h~=gcbo),'value',0);
 case 'exclusivebutton3'
  h = findobj(gcbf,'userdata','buttongroup3');
  set(gcbo,'value',1);
  set(h(h~=gcbo),'value',0);
 case 'exclusivebutton4'
  h = findobj(gcbf,'userdata','buttongroup4');
  set(gcbo,'value',1);
  set(h(h~=gcbo),'value',0);
 case 'exclusivebutton5'
  h = findobj(gcbf,'userdata','buttongroup5');
  set(gcbo,'value',1);
  set(h(h~=gcbo),'value',0);
 case 'comment'
  %disp('comment')	
  x = findobj(gcf,'tag','comments');
  set(x,'visible','on','string',FamilyLabel);
 case 'prepare'
  %  Variables for sorting
  startrange = findobj(gcf,'Tag','startrange');
  startrange = get(startrange,'string');
  endrange = findobj(gcf,'Tag','endrange');
  endrange = get(endrange,'string');
  %  We need to send the real epoch numbers to GetFamilies function, but
  %  we have this thing where epoch numbers start with a zero, which matlab
  %  doesn't like, so if we want all the epochs, this is zero to the last
  %  epoch, which is the length of the data minus one.
  reallength = length(EpochData)-1;
  %  Which epochs to search through: famrange
  startrange = str2num(startrange);
  if strcmp(endrange,'end')
	 endrange = reallength;
  else
	 endrange = str2num(endrange);
  end
  %  Find out which steps we want included in the families we pull out
  steplistepoch = findobj(gcf,'Tag','steplistepoch');
  steplistepoch = get(steplistepoch,'value');
  steplistall = findobj(gcf,'Tag','steplistall');
  steplistall = get(steplistall,'Value');
  %  No matter how the user searches the vector finalsteps will be returned.  
  %  This is a vector of all the steps the user wants included in the
  %  search, but we haven't checked yet to see if this is a minimum or not.
  %  If we are comparing to a certain epoch the variable steplistepoch is
  %  one.  If steplistall is one, then we want only the families that
  %  have all the steps listed (at this point doesn't matter if only or
  %  if - the minimim is the maximum)  If neither of these are a one,
  %  then the user has checked specific values and we must get these
  if steplistepoch == 1
	 stepepoch = findobj(gcf,'Tag','stepepoch');
	 stepepoch = get(stepepoch,'String');
	 stepepoch = str2num(stepepoch);
	 % Index number for the epoch to get the familysteps from
	 stepindex = find(EpochNumbers==stepepoch);
	 %  stepepoch will be the search reference epoch, if there is one, find out
	 %  what family it is in (stepfamily), which epochs are in that family 
	 %  (stepvalfamily) and then return the values of all steps in that family 
	 %  (finalsteps)
	 stepfamily = FamilyFlag(stepindex);
	 stepvalfamily = find(FamilyFlag==stepfamily);	 
	 finalsteps = FamilyStep(stepvalfamily);
	 %  If we don't want to limit the step values
  elseif steplistall == 1
	 %disp('unique')
	 finalsteps = unique(FamilyStep);
	 %  If we want to look at epochs that include particular epochs
  elseif steplistall == 0 & steplistepoch == 0
	 steplist = findobj(gcf,'Tag','steplist');
	 stepfamlist = get(steplist,'value');
	 steplength = length(stepfamlist);
	 ind = 0;
	 for i = 1:steplength
		if stepfamlist{i} == 1
		  ind = ind + 1;
		  stepvalue = get(steplist(i),'string');
		  stepvalist = str2num(stepvalue);
		  finalsteps(ind) = stepvalist;
		end
	 end
  end
  finalsteps;
  if isempty(finalsteps)
	 error('Must select a sort choice')
  end
  %  Find out whether we are looking at all families that include these
  %  steps, or only the ones that include exactly these steps.  Minim is
  %  one if we want all the families, zero if we want just the ones that
  %  have exactly the minimum steps. Since this is an either/or, we just
  %  check minimum, and assume if it isn't checked than only is.
  minim = findobj(gcf,'Tag','minim');
  minim = get(minim,'value');
  %  Which cues are being plotted?
  % First check if there is a cue designated as the title.  If there is
  % than we first get the epochnumbers for this cuevalue, as all the
  % families on the plot will be selected from this set.
  cueplottitle = findobj(gcf,'Tag','plottitle','value',1);
  titlepoch = findobj(gcf,'Tag','titlepoch');
  titlepoch = get(titlepoch,'value');
  if ~isempty(cueplottitle)
	 cueplotkey = get(cueplottitle,'Userdata');
	 cuekey = findobj(gcf,'Tag','cuevallist','Userdata',cueplotkey,'Value',1);
	 cuerow = get(cuekey,'userdata');
	 cuekey = get(cuekey,'string')
	 cuekey = str2num(cuekey);
	 if length(cuerow) > 1
		error('Must only have one value for cue in title')
	 end
	 [sortnumbers,sortfamilies] = GetFamilies(CellInfo,famnumber, ...
															finalsteps,cuekey, ...
															startrange,endrange,minim, ...
															cuerow);
	 % Now we will use sortnumbers as the input to further sort the
	 % families. We have used the startrange and endrange specified by the
    % user to get this subset of epoch numbers, so we can now use this
    % list as the range.  GetFamilies allows us to send in a vector as
    % startrange, instead of giving an upper and lower bound.  If there
    % is a vector given for startrange, it ignores endrange.  If there
    % was no title specified, than the original startrange and endrange
    % are still valid.
	 %
	 % If the user wants the max or max. area from all of the title cues,
    % we get that here, which means we need to get the data here just so
    % we can get the max and the max. area. 
    % 
  	 startrange = sortnumbers;
	 if titlepoch == 1
		[tempparseddata,familynumbs] = ParseFamilies(CellInfo, ...
																	sortnumbers, ...
																	sortfamilies);
		[tempmaxes,tempareas] = datamaxorarea(tempparseddata,xtime);
		tempfamilies.maxes = tempmaxes;
		tempfamilies.areas = tempareas;
	 end
  end
  cuelegend = findobj(gcf,'Tag','plotlegend','value',1);
  % If legend is checked, use that, otherwise get any cuevalue that is
  % checked.  Notice this means that if title is checked, but not legend,
  % then the cuevalue for title will be retrieve again. 
  if ~isempty(cuelegend)
	 cuelegend = get(cuelegend,'Userdata');
	 cueplot = findobj(gcf,'Tag','cuevallist','Userdata',cuelegend, ...
							 'value',1);	  
  else
	 cueplot = findobj(gcf,'Tag','cuevallist','value',1);
  end
  if isempty(cueplot)
	 cuevalue = '0';
	 cuerowind = 1;
	 cuelength = 1;
	 
	  %error('Must check one or more cues')
  else
  cuevalue = get(cueplot,'string')
  cuerowind = get(cueplot,'userdata')
  cuelength = length(cueplot);
  end
  % Can use GetFamilies for only one cueval at a time, so if cuelength is
  % more than one, must do a loop.
  if cuelength > 1
	 somerealnumbers = cell(1,cuelength);
	 somerealfamilies = cell(1,cuelength);
	 for i = 1:cuelength
		cueval(i) = str2num(cuevalue{i});
		cuerow(i) = cuerowind{i};
		[realnumbers,realfamilies] = GetFamilies(CellInfo,famnumber, ...
															  finalsteps,cueval(i), ...
															  startrange,endrange, ...
															  minim,cuerow(i)); 
		somerealnumbers{i} = realnumbers;
		somerealfamilies{i} = realfamilies;
	 end
  elseif cuelength == 1
	 cueval = str2num(cuevalue);
	 cuerow = cuerowind;
	 [realnumbers,realfamilies] = GetFamilies(CellInfo,famnumber, ...
															finalsteps,cueval, ...
															startrange,endrange,minim, ...
															cuerow);
	 somerealnumbers = cell(1);
	 somerealnumbers{:} = realnumbers;
	 somerealfamilies = cell(1);
	 somerealfamilies{:} = realfamilies;
  end

  somerealfamilies
  cueval

  ind = 0;
  temp = CellInfo.FamilyCondition(famnumber).EpochNumbers;
  for i = 1:cuelength
	 [tempparseddata,familynumbs] = ParseFamilies(CellInfo, ...
																  somerealnumbers{i}, ...
																 somerealfamilies{i});
	 templength = length(tempparseddata);
	 % Need to keep families separate if they come from different cues,
	 % because they will be normalized to different maxes.
	 tempfamilies.cue(i)=cueval(i);
	 tempfamilies.data{i}=tempparseddata;
	 tempfamilies.numbs{i}=familynumbs;
  end
  
  % We have three different ways of getting the max and area.  One is to get all the 
  % epochs that have the cuevalue selected, regardless of whether all of these epochs 
  % were actually selected (ie. maybe the user just chose epochs in the first 1/3 of
  % the data, but wants the max to consider all the epochs with th red
  % led).  This doesn't work very well for times when more than one cue
  % value was selected (ie. red LED where the points were 900), because
  % it is ambiguous which epochs to include for evaluating the max (ie.
  % all red and all points 900, or just epochs that are both red and 900,
  % or even all red, but doesn't matter which points).  In this case, the
  % user can choose to use the title epochs (if one is designated) to
  % determine the max, ie. if led is the title it looks at all the red
  % led to determine the max, regardless of how many points were in those
  % epochs. Finally, you can determine the max from only the data that has
  % been selected for each plot.  This will be the default. 
  
  %  For the case that the user wants to look at all epochs of a certain
  %  cuevalue for the max, regardless of whether they are actually
  %  chosen, we just get the max for all the cues.  Once we
  %  have it for each cue, it does not change for this data, no matter
  %  what particular family we are looking at.
  
  %  Regardless of how the max to normalize is chosen by, 
  %  we save these values, along with the data, in structures within a
  %  structure called tempfamilies, either as vectors or matrices.  For
  %  the case the max is for all cues of that type, order in the vectors
  %  is in ascending order of the cues, otherwise it is the same order as
  %  the data.
  
  % DO I CAUSE MEMORY PROBLEMS IF I DON'T INITIATE TEMPFAMILIES?
  cuetest = unique(cuerow);
  maxflag = 'n';
  % WE DON'T REALLY NEED TO DO THIS IF WE AREN'T NORMALIZING!!!
  allepoch = findobj(gcf,'Tag','allepoch');
  allepoch = get(allepoch,'Value');
  famgp = findobj(gcf,'Tag','famgp');
  famgp = get(famgp,'Value');
  if allepoch == 1 & cuetest == 1
		disp('Used max of all families with this cue')
	  	[cuemaxes,cueareas,checkcues] = CueMax(CellInfo,famnumber,cuetest);
	 	tempfamilies.rows = cuetest;
     	tempfamilies.maxes = cuemaxes;
	 	tempfamilies.areas = cueareas;
		maxflag = 1;
		% If the user is looking at a couple of different cues, for example, leds
		% and points, must use the (max, area) value unique to that data set, the
		% max of the red led at 600 pts. is going to be that anyway (well, I
		% guess the user could say he/she wanted to only look at the first 100
		% epochs, in which case there could be a higher max that used a red led
		% at 600 pts, but I think the user will have to deal with this on their
		% own.
  elseif famgp == 1
	 disp('The max is the max of just the families selected')
	 [allmaxes,allmaxarea] = datamaxorarea(tempfamilies.data,xtime);
	 tempfamilies.maxes = allmaxes;
	 tempfamilies.areas = allmaxarea;
  elseif titlepoch == 1
	 disp('Used epochs from title cue to select max')
  else
	 disp('If you want to look at all epochs with a certain cue to get');
	 disp('the max, need to only use one cuetype');
  end
  assignin('base','tempfamilies',tempfamilies);
  % tempfamilies.numbs is the list of epoch numbers in the families we have
  % chosen.  We want the index numbers, so we check in EpochNumbers for
  % the corresponding index value.
  % tempfamilies.data is the data we have selected, broken down into
  % families, so that there are many cells, and each one contains all the
  % data for a family
  % tempfamilies.cue gives the cue number that corresponds to each row of 
  % tempfamilies.data
  
  % Get y axis
  maxy = findobj(gcf,'Tag','maxy');
  maxy = get(maxy,'value');
  areay = findobj(gcf,'Tag','areay');
  areay = get(areay,'value');
  normmaxy = findobj(gcf,'Tag','normmaxy');
  normmaxy = get(normmaxy,'value');	
  normareay = findobj(gcf,'Tag','normareay');
  normareay = get(normareay,'value');
  response = findobj(gcf,'Tag','response');
  response = get(response,'value');
  if maxy == 1 | normmaxy == 1
	 yaxislabel = 'Maximum Response';
	 yvalue = cell(1,length(cueval));
	 for i=1:length(cueval)
		% When there is only one cue, and the user requests all epochs to
      % be included in the normalization max, then 
		% the maxes are indexed in ascending order, so use unique to
      % determine where the max will be located.  Otherwise, it has the
      % same index number as the data
		if titlepoch == 1
		  testcuenums = unique(cueval)
		  indexcuevalue = find(testcuenums==cueval(i))
		else
		  indexcuevalue = i;
		end
		thismax = tempfamilies.maxes(indexcuevalue);
		for j = 1:length(tempfamilies.data{i})
		  tempdata = tempfamilies.data{i}{j};
			yvalue{i}{j} = max(tempdata,[],2);
			if normmaxy == 1
			  yvalue{i}{j} = yvalue{i}{j}/thismax;
			end
		end
	 end
	 if normmaxy == 1
		yaxislabel = 'Normalized Maximum';
	 end
  elseif areay == 1 | normareay == 1
	 yaxislabel = 'Area';
	 yvalue = cell(1,length(cueval));
	 for i = 1:length(cueval)
		 % Area finds out the area for the entire epoch, including prepoints 
		 % and tailpoints. Might want to change this.
		 % When there is only one cue, and the user requests all epochs to
		 % be included in the normalization, then the maxes (of the areas)
       % are indexed in ascending order, so use unique to determine where
       % the max will be located.  Otherwise, it has the same index
       % number as the data 
		if titlepoch == 1
		  testcuenums = unique(cueval)
		  indexcuevalue = find(testcuenums==cueval(i))
		else
		  indexcuevalue = i;
		end
		thisarea = tempfamilies.areas(indexcuevalue);
		for j = 1:length(tempfamilies.data{i})
		  tempdata = tempfamilies.data{i}{j};
		  tempnumfam = size(tempdata,1);
		  for k = 1:tempnumfam
			 newtempdata = tempdata(k,:);
			 newtempdata(1) = newtempdata(1)/2;
			 newtempdata(end) = newtempdata(end)/2;
			 I = xtime*sum(newtempdata);
			 I = abs(I);
			 yvalue{i}{j}(k,1) = I; 
			 if normareay == 1
				yvalue{i}{j}(k,1) = yvalue{i}{j}(k)/thisarea;	 
			 end
		  end
		end
	 end
	 if normareay == 1
		yaxislabel = 'Normalized Area';
	 end
  elseif response == 1
	 yvalue = tempfamilies.data;
  end
  cueave = findobj(gcf,'Tag','cueave','value',1)
  stepave = findobj(gcf,'Tag','stepave','value',1);
  %stepave = get(stepave,'value')
  % This gives us the average 
  % I think if the user wants the average, we will just replace the data
  % structure with these values.
  yvalue
  %  If timex is tagged, there is no x value, other than the index of y,
  %  stepx or timex better be tagged!
  %  Need to know this before taking average.
  timex = findobj('Tag','timex');
  timex = get(timex,'value');
  if stepave == 1
	 % This takes the average of all cues, just combines all the data together and takes 
	 % average at each flash strength.  Since we may have started out with more than one
	 % group, we need to set somerealfamilies to 1, since they are all now combined.
	 if timex == 1
		ind = 0;
		yexample = size(yvalue{1}{1},1);
		epochave = cell(yexample,1);
		for i = 1:length(yvalue)
			yave = yvalue{i};
			for j = 1:length(yave)
				ind = ind + 1;
				for k = 1:yexample
					epochave{k}(ind,:) = yave{j}(k,:);
				end
			end
		end
		somerealfamilies = cell(1,yexample);
		meanvalue = cell(yexample,1);
		for i = 1:yexample
			yvalue{i} = mean(epochave{i},1);
			somerealfamilies{i} = i;
		end
	 else
		% Average of all cues for the case where not looking at epochs (ie. flash strength vs.
		% Max)
	 end
	 % For cueave, we take the average and replace ydata with it.  Does not care which cueave it
	 % is, because it only makes sense to ave by the cue that has more than one.  Not well tested,
  	 % since we don't have any data that has more than one useful cue.  Too bad we don't have a 
  	 % temp. assignment
  	 assignin('base','yvalue',yvalue);
  elseif ~isempty(cueave)
	  if timex == 1
  		for i = 1:length(yvalue)
			yave = yvalue{i};
			size(yave);
			disp('If given the message, CAT arguments dimensions are not')
			disp('consistent, this means data is not the same length, and')
			disp('cannot be averaged')
			ycat = cat(2,yave{:});
			ymean = mean(ycat,1);
			yvalue{i} = ymean;
		 end  
       else 
		 for i = 1:length(yvalue)
			yave = yvalue{i};
			size(yave);
			disp('If given the message, CAT arguments dimensions are not')
			disp('consistent, this means data is not the same length, and')
			disp('cannot be averaged')
			ycat = cat(2,yave{:});
			ymean = mean(ycat,2);
			yvalue{i} = ymean;
	 	end
	end
  end
  assignin('base','yvalue',yvalue);
  oldfigure = findobj('tag','Fig1');
  set(oldfigure,'visible','off');
  

  if timex == 1 & isempty(cueave) & isempty(stepave)
	 timeplot = input('Plot all families on the same plot? y/n ','s');  
	 if strcmp(timeplot,'n')
		timeplot = 1;
	 else
		timeplot = 0;
	 end
  else
	 timeplot = 0;
  end
  %
  h2 = figure('CloseRequestFcn','familyplay miniclose',...
				  'Position',[50 50 400 400], ...
				  'Tag','Fig2');	
  h3 = uicontrol('Parent',h2, ...
					  'Units','points', ...
					  'Callback','familyplay graphall', ...
					  'Position',[100 270 150 30], ...
					  'String','Graph all families', ...
					  'Tag','cuenumber', ...
					  'Userdata',cueval);
  h3 = uicontrol('Parent',h2, ...
					  'Units','points', ...
					  'FontSize',14, ...
					  'Position',[100 230 180 25], ...
					  'String','One Family at a time...', ...
					  'Style','text', ...
					  'Tag','gograph', ...
					  'Userdata',0);
  h3 = uicontrol('Parent',h2, ...
					  'Units','points', ...
					  'Callback','familyplay firstone', ...
					  'Position',[100 180 70 20], ...
					  'String','Start', ...
					  'Userdata','options');
  h3 = uicontrol('Parent',h2, ...
					  'Units','points', ...
					  'Callback','familyplay nextone', ...
					  'Position',[100 180 70 20], ...
					  'String','Next', ...
					  'Tag','nextgraph', ...
					  'Visible','Off');	
  h3 = uicontrol('Parent',h2, ...
					  'Units','points', ...
					  'Callback','familyplay breaktime', ...
					  'Position',[180 180 70 20], ...
					  'String','Interupt', ...
					  'Userdata','options');
  h3 = uicontrol('Parent',h2, ...
					  'Units','points', ...
					  'Callback','familyplay resetoptions', ...
					  'Position',[180 155 70 20], ...
					  'String','New Options', ...
					  'Tag','resetops', ...
					  'Visible','Off');  
  h3 = uicontrol('Parent',h2, ...
					  'Units','points', ...
					  'Callback','familyplay closeall', ...
					  'Position',[180 130 70 20], ...
					  'String','Quit', ...
					  'Tag','closeall');
  h3 = uicontrol('Parent',h2, ...
					  'Units','points', ...
					  'Position',[300 290 100 20], ...
					  'String','# of Bins', ...
					  'Style','checkbox', ...
					  'Tag','histbins', ...
					  'Visible','off');
  h3 = uicontrol('Parent',h2, ...
					  'Units','points', ...
					  'Position',[100 10 50 20], ...
					  'String','family #', ...
					  'Style','text', ...
					  'Tag','famcount', ...
					  'Userdata','options', ...
					  'Visible','off');
  h3 = uicontrol('Parent',h2, ...
					  'Units','points', ...
					  'Position',[150 10 50 20], ...
					  'String','', ...
					  'Style','text', ...
					  'Tag','numbers', ...
					  'Userdata','options', ...
					  'Visible','off');
  h3 = uicontrol('Parent',h2, ...
					  'Units','points', ...
					  'Position',[100 40 50 40], ...
					  'String','last family', ...
					  'Style','text', ...
					  'Userdata','options', ...
					  'Visible','off');
  h3 = uicontrol('Parent',h2, ...
					  'Units','points', ...
					  'Position',[150 50 50 20], ...
					  'String','', ...
					  'Style','text', ...
					  'Tag','morenumbers', ...
					  'Userdata','options', ...
					  'Visible','off');
  
 case 'graphtime'
  controlfig = gcf;
  y = findobj('tag','resetops');
  set(y,'visible','off');
  yvalue;
  yvalue{1};
  %  Variables for plotting	
  cuespot = findobj('Tag','cuelist');
  cuename = get(cuespot,'string');
  if isempty(cuename)
	  cuename = 'FamilyCondition';
  end
  cuespot2 = findobj('Tag','cuenumber');
  cueval = get(cuespot2,'Userdata');
  simplot = findobj('Tag','simplot');
  simplot = get(simplot,'value');
  histo = findobj('Tag','histo');	
  histo = get(histo,'value');
 
  % WHERE DO I PUT THIS??????
  if histo == 1
	 histbins = findobj('Tag','histbins')
	 set(histbins,'Visible','On');
  end
 
  logx = findobj('Tag','logx');
  logx = get(logx,'value');
  if logx
	 set(0,'DefaultAxesXScale','log')
  else
	 set(0,'DefaultAxesXScale','linear')
  end
  logy = findobj('Tag','logy');
  logy = get(logy,'value');
  if logy
	 set(0,'DefaultAxesYScale','log')
  else
	 set(0,'DefaultAxesYScale','linear')
  end
  timex = findobj('Tag','timex');
  timex = get(timex,'value');  
  stepx = findobj('Tag','stepx');
  stepx = get(stepx,'value');
  colorlist = ['mo-';'co-';'ro-';'go-';'bo-';'ko-';'mx-';'cx-';'rx-'; ...
					'gx-';'bx-';'kx-';'m*-';'c*-';'r*-';'g*-';'b*-';'k*-'];
  % Don't want circles on the plot, when plotting the responses, rather
  % than the max or the area.
  colorlist2 = ['m-';'c-';'r-';'g-';'b-';'k-';'m-'];
  % Once you start plotting, can't go back to the previous screen or quit
  % unless you interrupt first.
  y = findobj('tag','resetops');
  set(y,'visible','off');
  z = findobj('tag','closeall');
  set(z,'visible','off');
  % flag = 0  plot all families at once
  % flag = 1  plot next family
  % flag = 2  break
  % Here is how the loop works, when we hit a button, the userdata for
  % the object with tag 'gograph' is changed, and this is accessed
  % regularly by the variable flag. If we want to plot everything at
  % once, userdata/flag is changed to 0.  If we want to go through
  % families one at a time, it is changed to 1, and for interrupt (break),
  % changed to 2.  When flag is 1 we are in family sequence mode, where
  % we flip through families one at a time, this means that at the end of
  % each plotting, we check to see what the userdata of 'gograph' is.
  % This could be a 0 to finish plotting all at once, a 1 to keep going
  % one at a time, or a 2 to stop.  
  
  % Set the userdata back to 0 inside the while loop 
  xvalue = cell(1);
  flag;
  cnt = 0;
  cnt2 = 0;
  famset = 1;
  plotcolor = 0;
  startx = 1;
  while flag == 1 | flag == 0
	 flag;
	 numfamset = length(yvalue);
	 cueave
	 stepave
	 if isempty(cueave) & isempty(stepave)
		numfams = length(yvalue{famset})
	 else
		numfams = 1
	 end
	 nextstep = findobj('Userdata','options');
	 set(nextstep,'visible','On');
	 if isempty(stepave)
	 	cuestring = num2str(cueval(famset));
	else
		cuestring = 'All';
	end
	 cuename = cellstr(cuename);
	 temp = cuename(cuerow);
	 if strcmp(temp,'OutputChan') & isempty(stepave) 
		famlegend = CellInfo.OutputConfiguration{cueval(famset) + 1};
	 else
		famlegend = cuestring;
	 end
	 for family = 1:numfams
		 % cnt is number of times through the loop
		 % cnt2 starts over at one when time plot is one, but cnt keeps going.
		cnt = cnt + 1;
		cnt2 = cnt2 + 1;
		family
		% Why do I set userdata 0 here?
		whichone = findobj('Tag','gograph');
	 	set(whichone,'userdata',0);
		realfam = unique(somerealfamilies{famset});
		newfam = realfam(family);
		realnum = tempfamilies.numbs{famset};
		familycount = findobj('tag','numbers');
		set(familycount,'string',family);
		familytotal = findobj('tag','morenumbers');
		lengthallfams = length(yvalue{famset});
		set(familytotal,'string',lengthallfams);
		if stepx == 1
		  xaxislabel = 'Flash Strength';
		  if isfield(CellInfo.FamilyCondition(famnumber).UserInfo,'FamilyStep')
		  	temp = CellInfo.FamilyCondition(famnumber).UserInfo.FamilyStep;
		  else
			 disp('Warning: Data has not been calibrated')
			 temp = CellInfo.FamilyCondition(famnumber).FamilyStep;
		  end
		  famnum = realnum{family};
		  indexnumbers = famindex(CellInfo.FamilyCondition(famnumber),famnum);
		  nowxvalue = temp(indexnumbers);
		  if isempty(cueave)
			 xvalue{famset}{family} = nowxvalue';
		  else
			 xvalue{famset} = nowxvalue';
		  end
		elseif timex == 1
		  xaxislabel = 'Time';
		  disp('time')
		end
		nowlegend = ['family ' num2str(newfam)];
		if ~isempty(cueave) | stepave == 1
		  nowyvalue = yvalue{famset};
		  nowlegend = ['Average family'];
	    else	
		  nowyvalue = yvalue{famset}{family};
	    end	  
		%if isempty(cueave) & isempty(stepave)
		%  nowyvalue = yvalue{famset}{family};
		%elseif timex == 1 & isempty(cueave)
		%  nowyvalue = yvalue{famset}{family};
		%else
		%  nowyvalue = yvalue{famset};
		%  nowlegend = ['Average family'];
		%end
		if family == 1
		  nowlegend;
		  famlegend;
		  nowlegend = [nowlegend '  ' famlegend];
		end
		if cnt == 1 | timeplot == 1
		  % First time through? or doing one family per axes? Set up the
        % axis labels and title. 
		  newfigure = figure;
		  ylabel(yaxislabel);
		  xlabel(xaxislabel);
		  % what is the title?  if outputchan can use led.  cuerow comes
        % from the title designation, or is the only cue if none is
        % designated. 
		  cuerow;
		  temp = cuename(cuerow)';
		  if strcmp(temp,'OutputChan')
			 newtitle = [CellInfo.OutputConfiguration{cueval + 1} ' LED'];
		  else
			 cuestring = cellstr(cuestring);
			 newtitle = [temp cuestring];
		  end
		  title(newtitle);
		  cnt2 = 1;
		else
		  figure(newfigure);
		end
		hold on;
		plotcolor = plotcolor + 1;
		if plotcolor > length(colorlist)
		  plotcolor = 1;
		end
		colorlist(plotcolor,:);
		if logx == 1 & logy == 0
		  %disp('logx == 1 & logy == 0')
		  newplot = semilogx(nowxvalue,nowyvalue,colorlist(plotcolor,:));
		elseif logx == 0 & logy == 1
		  %disp('logx == 0 & logy == 1')
		  if timex == 1
			 newplot = semilogy(nowyvalue,colorlist(plotcolor,:));
		  else
			 newplot = semilogy(nowxvalue,nowyvalue,colorlist(plotcolor,:));
		  end
		elseif logx == 1 & logy == 1
		  %disp('logx == 1 & logy == 1')
		  newplot = loglog(nowxvalue,nowyvalue,colorlist(plotcolor,:));
		else
		  %disp('else')
		  if timex == 1
			 for i = 1:size(nowyvalue,1)
				if plotcolor > length(colorlist2)
				  plotcolor = 1;
				end
				if histo == 1
				  newplot = hist(nowyvalue(i,:),bins)
				else
				  newplot = plot(nowyvalue(i,:),colorlist2(plotcolor,:));
				end
			 end
		  else
			 if histo == 1
				newplot = hist(nowyvalue,bins)
			 else
				newplot = plot(nowxvalue,nowyvalue,colorlist(plotcolor,:));
			 end
		  end
		end
		leg_cell{1,cnt2}=[nowlegend];
		% set string to give family number
		familycount(cnt) = newfam;
		%legend(leg_cell,-1);
		if flag == 1 & family~=numfams
		  %disp('waiting')
		  figure(controlfig)
		  whichone = findobj('Tag','gograph');
		  waitfor(whichone,'userdata');
		  flag = get(whichone,'Userdata');
		  if flag == 2
			 break
		  end
		end
	 end
	 famset = famset + 1;
	 % if we broke out early, either famset or family will not be at the end
	 if famset == numfamset + 1 & family == numfams
		disp('All families plotted')
		flag = 2;
	 end
  end
  if flag == 2
	 disp('interupt or done')
	 x = findobj(gcf,'Tag','gograph');
	 set(x,'UserData',0);
	 %disp('restart loop')
	 x = findobj('tag','nextgraph');
	 set(x,'visible','off');
	 set(y,'visible','on');
	 set(z,'visible','on');
	 figure(controlfig)	 
  end
  assignin('base','xvalue',xvalue)
 case 'graphall'
  flag = 0;
  % ok this is klugey!
  y = findobj('tag','nextgraph');
  checky = get(y,'visible');
  checky = cellstr(checky);
  if strcmp(checky,'off')
	 familyplay graphtime
  else
	 % WHY DOESN'T THIS BLOODY WORK?
	 %whichone = findobj('Tag','gograph');
	 %set(whichone,'Userdata',0)
  end
 case 'firstone'
  flag = 1;
  x = findobj(gcf,'tag','gograph');
  set(x,'UserData',1)
  y = findobj('tag','nextgraph');
  set(y,'visible','on');
  familyplay graphtime
 case 'nextone'
  x = findobj(gcf,'Tag','gograph');
  set(x,'UserData',1);
 case 'breaktime'
  x = findobj(gcf,'Tag','gograph');
  set(x,'UserData',2);
  y = findobj('tag','nextgraph');
  set(y,'visible','off');
  z = findobj('tag','closeall');
  set(z,'visible','on');
 case 'resetoptions'
  newfigure = findobj('Tag','Fig2');
  close(newfigure)
  clear tempfamilies yvalue xvalue
 case 'miniclose'
  newfigure = findobj('Tag','Fig2');
  delete(newfigure)
  oldfigure = findobj('Tag','Fig1');
  set(oldfigure,'visible','on');
 case 'closeall'
  close all
  clear all
 case 'close'
  delete(gcf)
  clear all  
end

