function familysearch(a)
% familysearch(a)
% 
% Newest version, the epochs to search through are accurate, type in the epoch
% numbers you want to search through, and it will include all of those epochs in the
% search
% I think it really keeps all families, and only families.
%
% This is the first version of this function.  It is accessed through 
% epochs.  To use it in epochs, there are two things one must remember:
% the epochs interface must be loaded with an epoch that is in a family, 
% and the epochs must be fully loaded, ie. push to plot for the plot and
% stimulus attributes, before you hit the button find families.  Familysearch
% will find all epochs in families, according to the criteria given by the 
% user (the 'step').  For instance, if the user says the 'amp' changes in a
% family, then the program, searches for increases in amperage in
% successive  epochs and keeps track of the epochs in these families.  As
% long as there  is at least three epochs in a row that increase in
% amperage, that is considered  a family.   
%
% The program first asks you to mark all of the parameters that change in
% the family.  This is not just the parameters that mark the families you
% are interested in, but all parameters that change in the families.  If
% you have a parameter that changes sort of randomly, not marking it could
% cause you to not have any families at all in the end. Once you have
% selected all of the parameters that change, you are asked to select the
% parameter that defines the family, and given a choice of all of the
% parameters you marked as changing.   The program also sorts by 'cues',
% which means properties that change from family to family, such as temp.
% or LED color, or properties that change randomly throughout the families.
% If a parameter changes between families, like LED color, and you want all
% of the families of each different color, for instance, then you must mark
% the cues, otherwise it searches only for families where the LED color is
% the same as the epoch searched on.  The cues matrix will be fairly
% useless if the 'cues' you selected change randomly within families,
% rather than defining families you want to compare. Matrices the size of
% the number of epochs in a file are created, and coded so  that the
% program familyplay can sort them for plotting based on these cues and 
% steps.
%
%

% Created June, 2001 MKMK
% update February, 2003 MKMK - re-did epochs to search
% update March, 2004 MKMK - added in EpochTime to CellInfo structure
% update June, 2004 MKMK - fixed a problem getting the last family in a
% file, and using seed as a family criteria (now it makes sure an epoch is
% gassiannoise before checking the seed).
%  Oct 2004 MKMK Changed ITCFinishAnalysis to use evalin in order to get around
%  wierd bug in Matlab version 7
%  May 2005 MKMK Turns out it wasn't matlab; Fred fixed the bug in
%  ITCFinishAnalysis, so we don't have to use evalin anymore.


persistent fp fileptr plotpref segnum epochnumber epochdata familyitc datanum2 searchdata searchitc searchitc2 searchname searchfamstr sameitc samedata samename cuexist famindex NumEpochs FamilyCondition epochnums CellInfo structname newflag epochstart epochend


if nargin<1
    a = 'init';
    load analyzethis
end

switch a
    
    case 'init'
        %disp('familysearch')
        % if structname is empty, then this is a new file, clear CellInfo to get rid of old
        % junk
        if isempty(structname)
            CellInfo.CellFile = fp;
            newflag = 0;
            CellInfo.FamilyCondition = [];
            famindex = [];	
        else 
            newflag = 1;
        end	
        epochs close
        %  See if we want to limit the search (although this doesn't actually limit the
        %  search, after all searching is over it eliminates any epochs that aren't 
        %  suppose to be included in the search, this was because this feature was an 
        %  afterthought, and it is easier to adjust the epochs in the end than to 
        %  actually limit the search)
        home
        disp('Leave blank to include all epochs')
        searchepochs = input('Enter epoch # you would like to start search with: \n','s');
        test = str2num(searchepochs);
        if isempty(test) 
            epochstart = 0;
            epochend = 'end';
        else 
            epochstart = str2num(searchepochs);
            searchepochs = input('Enter epoch # you would like to end search with: \n','s');
            if isempty(searchepochs)
                epochend = 'end';
            else		
                epochend = str2num(searchepochs);
            end
        end
        epochstart;
        epochend;
        [fileptr, ecode] = ITCInitializeAnalysis(500000, fp);
        if ecode~=0
            error('Problem with ITCInitializeAnalysis')
        end
        h0 = figure('Color',[0.8 0.8 0.8], ...
            'MenuBar','none', ...
            'Position',[450 320 512 384], ...
            'Tag','Fig2');
        h1 = uicontrol('Parent',h0, ...
            'BackgroundColor',[0.8 0.8 0.8], ...  
            'FontSize',24, ...
            'Units','points', ...
            'Position',[32 331 455 28], ...
            'String','Which parameters change in the family?', ...
            'Style','text', ...
            'Tag','StaticText1');
        h1 = uicontrol('Parent',h0, ...
            'Callback','familysearch defined',...
            'Units','points', ...
            'FontSize',13, ...
            'Position',[320 20 150 35], ...
            'String','Push to continue');	  
        
        z = 1;
        familyitc = cell(20,1);
        datanum2 = cell(20,1);
        [epochsize, ecode2] = ITCGetEpochSize(epochnumber,fileptr);
        if ecode2 == 0
            stringname(z) = {'# Data Pts'};
            epochdata(z) = epochsize;
            familyitc(z) = {'ITCGetEpochSize'};
            datanum2{z} = 'EpochPts  ';
            z = z + 1;
        end	
        sampint = findobj(gcf,'tag','sampinterv');
        [SamplingInterval, ecode3] = ITCGetSamplingInterval(epochnumber,fileptr);
        if ecode3 == 0
            stringname(z) = {'Sampling Interval'};
            epochdata(z) = SamplingInterval;
            familyitc(z) = {'ITCGetSamplingInterval'};
            datanum2{z} = 'SampInterv';
            z = z + 1;
        end
        [OutputType,ecode4] = ITCGetOutputType(epochnumber, segnum, 0, fileptr);
        if ecode4 == 0
            stringname(z) = {'Output Type'};
            epochdata(z) = OutputType;
            familyitc(z) = {'ITCGetOutputType'};
            datanum2{z} = 'OutputType';
            checkmode = z;
            z = z + 1;
        end
        [mode, ecode6] = ITCGetAmpMode(epochnumber, segnum, fileptr);
        if ecode6 == 0
            stringname(z) = {'Amp Mode'};
            epochdata(z) = mode;
            familyitc(z) = {'ITCGetAmpMode'};
            datanum2{z} = 'AmpMode   ';
            z = z + 1;
        end
        [InputChan, ecode7] = ITCGetInputChan(epochnumber, segnum, fileptr);
        if ecode7 == 0
            stringname(z) = {'Input Channel'};
            epochdata(z) = InputChan;
            familyitc(z) = {'ITCGetInputChan'};
            datanum2{z} = 'InputChan ';
            z = z + 1;
        end
        [OutputChan, ecode8] = ITCGetOutputChan(epochnumber, segnum, fileptr);
        if ecode8 == 0
            stringname(z) = {'Output Channel'};
            epochdata(z) = OutputChan;
            familyitc(z) = {'ITCGetOutputChan'};
            datanum2{z} = 'OutputChan';
            z = z + 1;
        end
        if epochdata(checkmode) == 2
            [Seed, ecode10] =  ITCGetSeed(epochnumber, segnum, 0, fileptr);
            if ecode10 == 0	
                stringname(z) = {'Seed'};
                epochdata(z) = Seed;
                familyitc(z) = {'ITCGetSeed'};
                datanum2{z} = 'Seed      ';
                z = z + 1;
            end
        end
        [StmAmp, ecode11] = ITCGetStmAmp(epochnumber, segnum, 0, fileptr);
        if ecode11 == 0
            stringname(z) = {'Amp '};
            epochdata(z) = StmAmp;
            familyitc(z) = {'ITCGetStmAmp'};
            datanum2{z} = 'StmAmp    ';
            z = z + 1;	
        end
        [StmMean, ecode12] = ITCGetStmMean(epochnumber, segnum, 0, fileptr);
        if ecode12 == 0	
            stringname(z) = {'Mean'};
            epochdata(z) = StmMean;
            familyitc(z) = {'ITCGetStmMean'};
            datanum2{z} = 'StmMean   ';
            z = z + 1;
        end
        [PrePoints, ecode13] = ITCGetStmPrePts(epochnumber, segnum, 0, fileptr);
        if ecode13 == 0	
            stringname(z) = {'Prepoints'};
            epochdata(z) = PrePoints;
            familyitc(z) = {'ITCGetStmPrePts'};
            datanum2{z} = 'PrePoints ';
            z = z + 1;
        end
        [Points, ecode14] = ITCGetStmPts(epochnumber, segnum, 0, fileptr);
        if ecode14 == 0	
            stringname(z) = {'Points'};
            epochdata(z) = Points;
            familyitc(z) = {'ITCGetStmPts'};
            datanum2{z} = 'StimDur   ';
            z = z + 1;
        end
        [TailPoints, ecode15] = ITCGetStmTailPts(epochnumber, segnum, 0, fileptr);
        if ecode15 == 0	
            stringname(z) = {'Tailpoints'};
            epochdata(z) = TailPoints;
            familyitc(z) = {'ITCGetStmTailPts'};
            datanum2{z} = 'TailPoints';
            z = z + 1;
        end
        if z > 10
            for i = 1:10
                gencheckbox = uicontrol('Parent',h0, ...
                    'Position',[20 330-(i*30) 125 20], ...
                    'String',stringname(i), ... 
                    'Style','checkbox', ...
                    'Userdata',i);	
            end
            for i = 11:z-1
                gencheckbox = uicontrol('Parent',h0, ...
                    'Position',[180 330-((i-10)*30) 125 20], ...
                    'String',stringname(i), ... 
                    'Style','checkbox', ...
                    'Userdata',i);	
            end
        else
            for i = 1:z-1
                gencheckbox = uicontrol('Parent',h0, ...
                    'Position',[20 330-(i*30) 125 20], ...
                    'String',stringname(i), ... 
                    'Style','checkbox', ...
                    'Userdata',i);	
            end
        end
%        datanum2;
    case 'defined'
        %  First get all the checkboxes that are not checked.  We will use these
        %  to compare all values that should be the same in the family.
        %  While searching searchtag defines the family, searchbox are the
        %  parameters that change (i.e. led color), and all other parameters
        %  are the same.
        sametag = findobj(gcf,'Style','checkbox','Value',0);
        indata = get(sametag,'userdata');
        samedata = zeros(1,length(indata));
        sameitc = cell(1,length(indata));
        %samename = cell(1,length(indata));
        samename = zeros(length(indata),10);
        %datanum2;
        %datanum2{1};
        for a = 1:length(indata)
            samedata(a) = epochdata(indata{a});
            sameitc{a} = familyitc(indata{a});
            samename(a,:) = datanum2{indata{a}};
        end
        samename = char(samename);
        samename = cellstr(samename);
        %  Now get the checkboxes that are checked.  These values will change
        %  either between the families or will define the family.
        searchbox = findobj(gcf,'Style','checkbox','Value',1);
        searchstring = get(searchbox,'string');
        ind = get(searchbox,'userdata');
        cuexist = iscell(ind);
        %  If cuexist is a cell, then more than one box was checked, and we have to 
        %  find out which one is the step (changes within a family) and which one(s)
        %  is/are the cues (changes from family to family).  If there is only one 
        %  checked, then we assume all values except the step were the same.
        if cuexist == 1
            searchdata = zeros(1,length(ind));
            searchitc = cell(1,length(ind));
            searchname = zeros(length(ind),10);
            for a = 1:length(ind)
                searchdata(a) = epochdata(ind{a});
                searchitc{a} = familyitc(ind{a});
                searchname(a,:) = datanum2{ind{a}};
            end
            searchname = char(searchname);
            searchname = cellstr(searchname);
            strlen = length(searchbox);
            delete(gcf)
            h0 = figure('Color',[0.8 0.8 0.8], ...
                'Position',[450 320 512 384], ...
                'Tag','Fig2');
            h1 = uicontrol('Parent',h0, ...
                'Units','points', ...
                'Position',[32 331 455 28], ...
                'String','Which parameter defines the family (ie. what changes within each family)?', ...
                'Style','text', ...
                'Tag','StaticText1');
            h1 = uicontrol('Parent',h0, ...
                'Callback','familysearch selected',...
                'Units','points', ...
                'FontSize',13, ...
                'Position',[320 20 150 35], ...
                'String','Push to continue');	 
            for i = 1:strlen
                gencheckbox = uicontrol('Parent',h0, ...
                    'Position',[20 330-(i*30) 125 20], ...
                    'String',searchstring{i}, ... 
                    'Style','checkbox', ...
                    'Userdata',i);
            end	
        else
            searchdata = epochdata(ind);
            searchitc = familyitc(ind);
            searchname = datanum2{ind};
            familysearch selected
        end
    case 'selected'
        %  First get step stuff (what is checked)
        searchtag = findobj(gcf,'Style','checkbox','Value',1);
        ind = get(searchtag,'userdata');
        if length(ind)>1
            error('Check only one box, please');
        end
        if cuexist == 1
            %  This is for cues
            searchfamdat = searchdata(ind);
            searchfamitc = searchitc{ind}{1};
            searchitc{ind}{1}=[];
            searchfamstr = searchname(ind,:);
            searchfamstr = char(searchfamstr);
            searchfamstr = cellstr(searchfamstr);		
            % To get the cues, we made the searchitc empty for the one that will
            % be the step, so any others will be cues
            j = 1;
            for i = length(searchitc):-1:1
                if ~isempty(searchitc{i}{1})
                    searchitc2{j}{1} = searchitc{i}{1};
                    j = j + 1;
                end
            end
        else
            %  if no cues than 
            searchfamdat = searchdata;
            searchfamitc = searchitc{1};
            searchfamstr = searchname;
            searchfamstr = char(searchfamstr);
            searchfamstr = cellstr(searchfamstr);		
        end
        [NumEpochs, ecode5] = ITCGetNumberEpochs(fileptr);
        if ecode5~=0
            error('Problem with the number of epochs in the file')
        end
        flag1 = 1;
        flag2 = 1;
        i = epochstart;
        n = 1;
        start = [];
        finish = [];
        FamilyCondition.FamilyStep = zeros(1,NumEpochs);
        %NumEpochs;
        endfile = 0;
        if strcmp(epochend,'end') | epochend == NumEpochs - 1
            disp('here')
            endfile = 1;
            %  NumEpochs is the length of the list of epoch numbers, since the list 
            %  starts with zero, the last epoch number is NumEpochs - 1.  Since we
            %  are comparing i to i + 1, we can only go to the last epoch number - 1,
            %  if the user left it blank, this is 
            epochend = NumEpochs - 2;
        end
        %epochend;
        %  The loop makes a comparison between i and the next highest, we have put in 
        %  a failsafe, so that epochend is actually one less than the end of the file,
        %  so that when i is at a max, i + 1 still exists. This means we won't get the 
        %  last data point in Familystep, so we add in blah1 (So as long as the next to last
        %  is part of a family, the last is to) to the last Familystep.  
        %  FamilyStep at this point, is all of the values (since we were checking the
        %  values here anyway, we just pulled them out, values not in a family will
        %  be changed to NaNs later, when we know if all of the families sorted out 
        %  here are really going to be in the final list).
        %searchfamitc;
        while i < epochend + 1
            if strcmp(searchfamitc,'ITCGetEpochSize') | ...
                    strcmp(searchfamitc,'ITCGetSamplingInterval')
                [blah,error1] = feval(searchfamitc, i, fileptr);
                [blah1,error1] = feval(searchfamitc, i+1, fileptr);
            elseif strcmp(searchfamitc,'ITCGetAmpMode') | ...
                    strcmp(searchfamitc,'ITCGetInputChan') | ...
                    strcmp(searchfamitc,'ITCGetOutputChan')
                [blah,error1] = feval(searchfamitc, i, segnum, fileptr);
                [blah1,error1] = feval(searchfamitc, i+1, segnum, fileptr);
            else
                [blah,error1] = feval(searchfamitc, i, segnum, 0, fileptr);
                [blah1,error1] = feval(searchfamitc, i+1, segnum, 0, fileptr);
            end
            %  blah is just the value of the comparison value (ie. amp) for the current 
            %  epochnumber and blah1 is the value for the next epochnumber,
            %  we compare these two, blah (i) and blah1 (i+1). flag1 and flag2 start at
            %  1, if the numbers are increasing, the flag1 changes to zero, the start 
            %  is marked, 
            %  and i increases.  If the numbers increase again, flag1 remains zero, and 
            %  flag2 also becomes zero (this is to prevent a run of 2 to be considered 
            %  a family), and i increases.  If the numbers decrease immediately after the 
            %  first initial increase, then it was an 'accident' and not a family, so both 
            %  flags return to 1.  If the numbers are increasing, and it is not an accident, 
            %  both flags are zero until we detect the numbers stop increasing.  When this 
            %  happens, the finish is marked, flag1 and flag2 revert back to 1, and n 
            %  (family #) increases.  When we get to the end, if we are in a run (flag2 = 0), 
            %  then we make this the last in the run.
            %blah;
            %blah1;
            %flag1;
            %flag2;
            %i;
            %n;
            % if this epoch is smaller than the next one, and we are not in a family already,
            % then we start a family
            if blah < blah1 & flag1 == 1 & flag2 == 1
                start(n) = i;
                %disp('start')
                i = i + 1;
                flag1 = 0;
                % check to see if this is end of file    
            elseif i == epochend & flag2 == 0
                if endfile == 1
                    %disp('end of file, last finish')
                    blah;
                    blah1;
                    finish(n) = i + 1;
                    i = i + 1;
                else
                    %disp('yes - finish')
                    finish(n) = i;
                    i = i + 1;
                    flag1 = 1;
                    flag2 = 1;
                    n = n + 1;
                end
                % If 3 in a row have grown consequtively larger, than this is a family
            elseif blah < blah1 & flag2 == 1
                %disp('not an accident')
                i = i + 1;
                flag2 = 0;
                % We are in a family - keep going!
            elseif blah < blah1 & flag2 == 0
                %disp('in family')
                i = i+ 1;
                % But, if it was just 2 epochs that happened to be one larger than the other,
                % we can start again
            elseif blah > blah1 & flag2 == 1
                %disp('accident')
                i = i + 1;
                flag2 = 1;
                flag1 = 1;
                % If we are in a family (flags are zero), and values stop growing, we are at the 
                % end of a family
            elseif blah >= blah1 & flag2 == 0
                %disp('finish')
                % Not sure why I put this check in - this should never happen, but it can't 
                % hurt to leave it in.
                % If there was not a beginning to this family, there shouldn't be an end!
                if length(start) ~= n | isempty(start(n))
                    %disp('no beginning')
                    i = i + 1;
                    flag1 = 1;
                    flag2 = 1;
                else
                    %disp('yes - finish')
                    finish(n) = i;
                    i = i + 1;
                    flag1 = 1;
                    flag2 = 1;
                    n = n + 1;
                end
            else
                % not in a family - make sure flags are reset
                %disp('else')
                flag1 = 1;
                flag2 = 1;
                i = i + 1;
            end
            FamilyCondition.FamilyStep(i) = blah;
            FamilyCondition.FamilyStep(i+1) = blah1;
        end	
        % It is possible a family was started and never finished.  Turns out this doesn't 
        % make any difference, but lets make it pretty.
        if length(start) > length(finish)
            start = start(1:length(finish));
        end
        %start;
        %finish;
        FamilyCondition.FamilyStep(:);
        % We now have all of the possible family candidates from the range given, now 
        % check to see all parameters that are suppose to be the same are,
        % in fact, the same.
        i = 1;
        finalfamily = length(finish);
        FamilyCondition.FamilyFlag = nans(1,NumEpochs);
        FamilyCondition.EpochNumbers = zeros(1,NumEpochs);
        % i is the family #, j is the epoch #, k is an index for the ITC stuff
        % Want it so it checks each itc, and as soon as one comes up not equal, 
        % that epoch is flagged to not include
        flag1 = 1;
        k = 1;
        ind = 1;
        flag2 = 1;
        %finish(finalfamily)
        %finalfamily
        while i < finalfamily + 1
            j = start(i);
            while j < (finish(i) + 1)
                % If the start is zero, and we aren't on the first epoch, this is not a 
                % good family.
                if i ~= 1 & j == 0
                    j = finish(i) + 1;
                end
                %  Searches through all parameters (sameitc) unless one is different 
                %  from the example epoch, in which case flag1 is set to zero, and the 
                %  loop goes to the next epoch.
                %  Start at flag1 = 1, flag1 affects while loop, flag2 does not.
                %  flag1 = 1 continue in family, changes if we get an epoch that is not
                %  	the same as the example epoch
                %  flag1 = 0 epoch is not in a family, stop loop, make NaN, INCREASE IND??
                %  Loop is stopped for 2 reasons:
                %		flag1 = 0
                %		have checked all parameters (k>length(sameitc))
                %	When loop is stopped, check why:
                %		k max reached, reset k to one, increase epochnumber, go again, 
                %     	flag2 is set to zero meaning we are in a family.
                %   	or flag1 = 0, epoch not in family, k reset to one, Nans put in 
                % 			FamilyCondition, increase epochnumber, flag1 reset to one, 
                %			flag2 = 1,  meaning not in family.
                %  flag2 = 0 have looked at all attributes, and everything is good,
                %  	start k over, go to next epochnumber.   If we make it through 
                %  	family, and flag2 = 0, then family # increases.  
                %	PROBLEM  If the step increased to a non-family member by chance,
                %		this step would be included in the family, the family would end 
                %		with a NaN, and the family # would not increase.
                while flag1 == 1 & k < length(sameitc) + 1
                    if strcmp(sameitc{k}{1},'ITCGetEpochSize') | ...
                            strcmp(sameitc{k}{1},'ITCGetSamplingInterval')	
                        [blah,error1] = feval(sameitc{k}{1}, j, fileptr);
                        %sameitc{k}{1}
                    elseif strcmp(sameitc{k}{1},'ITCGetAmpMode') | ...
                            strcmp(sameitc{k}{1},'ITCGetInputChan') | ...
                            strcmp(sameitc{k}{1},'ITCGetOutputChan')
                        [blah,error1] = feval(sameitc{k}{1}, j, segnum, fileptr);
                        %sameitc{k}{1}
                    else
                        [blah,error1] = feval(sameitc{k}{1}, j, segnum, 0, fileptr);	
                        %sameitc{k}{1}
                    end
                    if blah == samedata(k)
                        k = k + 1;
                        flag1 = 1;
                    else
                        flag1 = 0;
                    end
                end
                %  Loop ended, did we get through all of the checks?  
                if k == length(sameitc) + 1
                    %disp('ok')
                    %  Yes, this epoch is ok, go to next one.
                    FamilyCondition.FamilyFlag(j+1) = ind;
                    FamilyCondition.EpochNumbers(j+1) = j;
                    j = j + 1;
                    k = 1;
                    %flag2 = 1;
                elseif flag1 == 0
                    %disp('not ok')		
                    %  If we didn't make it through the checks, No, this epoch is not ok.
                    %  Use NaNs to make sure unambiguous which values should be kept
                    %  What if this is the first or the last epoch?  Should just get rid
                    %  of this epoch and keep the rest of the family intact.
                    FamilyCondition.FamilyFlag(j+1) = NaN;
                    FamilyCondition.EpochNumbers(j+1) = NaN;
                    %  No matter what happens, flag1 = 1 and k = 1, so we can look at 
                    %  all of the parameters of the next epoch
                    flag1 = 1;
                    k = 1;
                    %  If we are at the first epoch, check to see how long the family is, 
                    %  without this epoch.  Must be more than 2 left, or three total, to be 
                    %  considered a family
                    if j == start(i)
                        %disp('start')
                        % See how long the family is, could be an accident, not a family
                        if finish(i) - start(i) < 4 
                            %disp('end really early')
                            FamilyCondition.FamilyFlag(start(i)+1:finish(i)+1) = NaN;
                            FamilyCondition.EpochNumbers(start(i)+1:finish(i)+1) = NaN;
                            j = finish(i) + 1;
                        else
                            %disp('get rid of first epoch')
                            j = j + 1;
                            flag1 = 1;
                        end		
                        %  Do the same check for the 2nd epoch.  Now getting rid of 2, so in order
                        %  to have 3 or more left must be 5 altogether.
                    elseif j == start(i) + 1
                        %disp('early?')
                        % See how long the family is, could be an accident, not a family
                        if finish(i) - start(i) < 5 
                            %disp('end early')
                            FamilyCondition.FamilyFlag(start(i)+1:finish(i)+1) = NaN;
                            FamilyCondition.EpochNumbers(start(i)+1:finish(i)+1) = NaN;
                            j = finish(i) + 1;
                        else
                            %disp('get rid of 2nd epoch')
                            j = j + 1;
                            flag1 = 1;	
                        end
                        %  If we are at the last epoch, get rid of this epoch, and increase
                        %  the family.
                    elseif j == finish(i)
                        %disp('finished family')
                        j = j + 1;
                        flag2 = 1;
                        %  If we are in the middle, something is wrong with the whole 
                        %  family, and we should go to the end, but not increase the family
                    else
                        %disp('else')
                        while j < finish(i) + 1
                            FamilyCondition.FamilyFlag(j+1) = NaN;
                            FamilyCondition.EpochNumbers(j+1) = NaN;
                            j = j + 1;
                            flag2 = 0;
                        end
                    end	
                end
                %j
            end
            %  Increase family if flag2 = 1, reset flag2 regardless
            if flag2 == 1
                ind = ind + 1;
                flag2 = 1;
            else
                flag2 = 1;
            end
            i = i + 1;
        end		
        %FamilyCondition.EpochNumbers(:);
        %FamilyCondition.FamilyFlag(:);
        FamilyCondition.SearchCrit = samename;
        FamilyCondition.SearchPara = samedata;
        %  If there are cues, go on to cues, otherwise continue
        if cuexist == 1
            familysearch cueselect
        else
            FamilyCondition.FamilyCues = zeros(size(FamilyCondition.EpochNumbers));
            FamilyCondition.FamilyCueGuide = [];
            familysearch continue
        end
    case 'cueselect'
        cuetag = findobj(gcf,'Style','checkbox','Value',0);
        test = get(cuetag,'userdata');
        assignin('base','test',test);
        assignin('base','searchname',searchname);
        %searchcuestr = zeros(length(test),10);
        notcue = findobj(gcf,'Style','checkbox','Value',1);
        nottest = get(notcue,'userdata');
        searchcuestr = cell(length(test),1);
        for y = 1:length(test)
            if iscell(test)
                newtest = test{y};
            else
                newtest = test(y);
            end
            searchcuestr{y} = searchname{newtest};
        end
        
        % 	for y = 1:length(test)
        % 		if iscell(test)
        % 			newtest=test{y};
        % 		else 
        % 			newtest = test(y);
        % 		end
        % 		temp = zeros(1,10)
        % 		temp = searchname{newtest}
        % 		searchcuestr(y,:) = temp 
        % 	end
        %searchcuestr;
        %segnum;
        %	searchcuestr = char(searchcuestr);
        %	searchcuestr = cellstr(searchcuestr);
        FamilyCondition.FamilyCues = zeros(length(searchitc2),NumEpochs);
        i = epochstart;
        while i < epochend + 1
            %while i < NumEpochs
            if FamilyCondition.FamilyFlag(i+1) ~= NaN
                for k = 1:length(searchitc2)
                    if strcmp(searchitc2{k}{1},'ITCGetEpochSize') | ...
                            strcmp(searchitc2{k}{1},'ITCGetSamplingInterval')	
                        [blah,error1] = feval(searchitc2{k}{1}, i, fileptr);
                    elseif strcmp(searchitc2{k}{1},'ITCGetAmpMode') | ...
                            strcmp(searchitc2{k}{1},'ITCGetInputChan') | ...
                            strcmp(searchitc2{k}{1},'ITCGetOutputChan')
                        [blah,error1] = feval(searchitc2{k}{1}, i, segnum, fileptr);
                    else
                        % can't check the seed, unless the stimulus is
                        % gaussiannoise.
                        if strcmp(searchitc2{k}{1},'ITCGetSeed')
                            [OutputType,ecode4] = ITCGetOutputType(i, segnum, 0, fileptr);  
                            if OutputType == 2
                                [blah,error1] = feval(searchitc2{k}{1}, i, segnum, 0, fileptr);
                            else 
                                blah=nan;
                            end
                        else
                            [blah,error1] = feval(searchitc2{k}{1}, i, segnum, 0, fileptr);
                        end
                    end
                    FamilyCondition.FamilyCues(k,i+1) = blah;
                end
            end
            i = i + 1;
        end
        FamilyCondition.FamilyCueGuide = searchcuestr;
        familysearch continue
    case 'continue'
        repochs = isnan(FamilyCondition.FamilyFlag);
        rmepochs = find(repochs);
        FamilyCondition.EpochNumbers(rmepochs) = [];
        finalength = length(FamilyCondition.EpochNumbers);
        FamilyCondition.EpochTime = zeros(1,finalength);
        for i = 1: finalength
            ind = FamilyCondition.EpochNumbers(i);
            [EpochTime, ecode] = ITCGetEpochTime(ind, fileptr);
            FamilyCondition.EpochTime(i) = EpochTime;
        end
        FamilyCondition.FamilyFlag(rmepochs) = [];
        FamilyCondition.FamilyStep(rmepochs) = [];
        FamilyCondition.FamilyCues(:,rmepochs) = [];		
        FamilyCondition.FamilyStepGuide = searchfamstr;
        FamilyCondition.SegNum = segnum;
        FamilyCondition.PlotPref = plotpref;
        FamilyCondition.ExcludeEpochs = [];
        FamilyCondition.Label = [];
        FamilyCondition.ScaleFactorIndex = [];
        FamilyCondition.DecimatePts = 1;
        FamilyCondition.UserInfo = [];
        %FamilyCondition.FamilyStep(:)
        %FamilyCondition.EpochNumbers(:)
        %FamilyCondition;
        if isempty(famindex)
            CellInfo.FamilyCondition = deal(FamilyCondition);
        else
            CellInfo.FamilyCondition(famindex) = deal(FamilyCondition);
        end
        delete(gcf)  
        ITCFinishAnalysis(fileptr);
        if newflag == 0
            save('analyzethis','CellInfo')
            playtime('init','analyzethis',2)
        else
            save(structname,'CellInfo')
            playtime('init',structname,2,famindex)
        end
end






