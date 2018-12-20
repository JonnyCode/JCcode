function [cuemaxes,cueareamaxes,checkcues] = CueMax(CellInfo,famnum,cuenum)
% function [cuemaxes,cueareamaxes,checkcues] = cuemax(CellInfo,famnum,cuenum)
% Looks at all of the cues (meaning all for the cues in the selected row
% (cuenum) of the cue matrix, ie. all values for the OutputChan if that
% is the row selected), for the FamilyCondition (famnum is the index
% number of which family we are interested in) given, returns a
% vector that has the max for each cue, and a vector that has the maximum
% area for each cue.  The maxes will be in ascending order of the cues.
% Ie., if the cues are 1,3,4, then cuemax(1) corresponds to cue = 1,
% cuemax(2) to cue = 3, etc. Same for the max of the areas.  Checkcues is
% the cues in same order as cuemaxes and cueareamaxes.

% Created MKMK Oct. 2001

famcell = CellInfo.FamilyCondition(famnum);
EpochNumbers = famcell.EpochNumbers;
EpochData = CellInfo.EpochData.Data;
FamilyStep = famcell.FamilyStep;
FamilyFlag = famcell.FamilyFlag;
FamilyCues = famcell.FamilyCues(cuenum,:);

% Do we use min or max?  Use two epochs to be sure.
testepoch = max(EpochData{EpochNumbers(1)+1});
testepochb = min(EpochData{EpochNumbers(1)+1});
testepoch2 = max(EpochData{EpochNumbers(2)+1});
testepoch2b = min(EpochData{EpochNumbers(2)+1});
if abs(testepoch)>abs(testepochb) & abs(testepoch2)>abs(testepoch2b)
  % use max
  mdefault = 1;
else
  % use min
  mdefault = 0;
end
checkcues = unique(FamilyCues);
checknum = length(checkcues);
cuemaxes = zeros(1,checknum);
cueareamaes = zeros(1,checknum);
for i = 1:checknum
  checknums = find(FamilyCues==checkcues(i));
  cueepochnums = EpochNumbers(checknums) + 1;
  numslen = length(cueepochnums);
  allmax = zeros(1,numslen);
  if mdefault == 1
	 for j = 1:numslen
		allmax(j) = max(EpochData{cueepochnums(j)});
	 end
  elseif mdefault == 0
	 for j = 1:numslen
		allmax(j) = abs(min(EpochData{cueepochnums(j)}));
	 end
  end
  cuemaxes(i) = max(allmax);
  % Now the area
  for j = 1:numslen
	 temp = EpochData{cueepochnums(j)};
	 h = 0.1;
	 temp(1) = temp(1)/2;
	 temp(end) = temp(end)/2;
	 I = h*sum(temp);
	 I = abs(I);
	 cuearea(j) = I;
  end
  cueareamaxes(i) = max(cuearea);
end
