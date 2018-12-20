function [realnumbers,epochfamilies] = GetFamilies(CellInfo,condition, ...
							  steps,cue,lowrange,highrange,minsteps,cueind)
%  [realnumbers,epochfamilies] =
%  GetFamilies(CellInfo,condition,steps,cue,lowrange,highrange,minsteps)
%
%  Gives back a list of epochnumbers that correspond to the steps and cue
%  given. Right now returns the actual epochnumbers, not an index number to the
%  epochdata in CellInfo (ie. must add a one to the numbers to index the
%  data). Epochfamilies is a list of which families those epochnumbers
%  belong to. 
%  Include only one cue.  If you want more cues, run it more than once.
%  Default for lowrange/highrange is all epochs in FamilyCondition, this
%  function assumes the range is the actual epoch numbers
%  Minsteps is 1 if the user wants all families that include these steps, and
%  0 if the user wants only the families that include exactly these
%  steps and no more, default is 1

%  Created MKMK Oct. 2001

if nargin < 5
  lowrange = 0;
end
if nargin < 6
  highrange = [];
end
if nargin < 7
  minsteps = 1;
end
famgp = CellInfo.FamilyCondition(condition);
%  allepochs is a list of all epochs in the familycondition (real epoch #s)
allepochs = famgp.EpochNumbers;
allsteps = famgp.FamilyStep;
allflags = famgp.FamilyFlag;
steplength = length(steps);
cueind;

%  First we look at first step, and find all the index numbers that the
%  step appears. 
stepepochs = find(allsteps==steps(1));
%  Now get the family #s for all these epochs
stepfamilies = allflags(stepepochs);
stepfamlength = length(stepfamilies);
%  Do a loop to check that all steps requested are in the families
%  pulled.
t = 1;
flag = 1;
for i = 1:stepfamlength
  % checkfam is the list of families that have the 1st step in them
  % checkfamsteps is a vector of all of the steps included in the family
  checkfam = find(allflags==stepfamilies(i));
  checkfamsteps = allsteps(checkfam);
  %  check the steps included in the family and make sure steps requested
  %  are present in each family.
  j = 1;
  while j < steplength + 1 & flag == 1;
	 if find(checkfamsteps==steps(j))
		flag = 1;
		j = j + 1;
	 else
		flag = 0;
	 end
  end
  if flag == 1
	 goodfamily(t) = i;
	 t = t + 1;
  end
end
%  goodfamily should be a list of families that include all of the steps
%  we are looking for.
goodfamily;
%  Do we care about all these families or only the ones that have exactly
%  the steps we requested?  We know these families have all of the steps
%  requested, so if there are more epochs than the steplength, than there
%  are more steps.  Check to see if we care.
if minsteps == 0
  %disp('minsteps')
  for i = 1:length(goodfamily)
	 famlength = find(allflags==goodfamily(i));
	 if length(famlength) > steplength
		goodfamily(i) = 0;
	 end
  end
  % WHY IS FAMILY 12 SO LONG -- MUST BE TWO FAMILIES TOGETHER!!!  WHY!!?!
  % THIS IS A PROBLEM WITH FAMILYSEARCH I THINK, IT COMES FROM CELLINFO
  % JUST NEED TO GET EPOCHNUMBERS INTO THE CORRECT VARIABLE NAME.  DON'T
  % REMEMBER EXACTLY WHAT THE VARIABLE NAME SHOULD BE.
end
stepepochs = cell(1,length(goodfamily));
for i = 1:length(goodfamily)
  stepepochs{i} = find(allflags==goodfamily(i));
end
catstepepochs = cat(2,stepepochs{:});
%  It could be that there are more steps in the family than what was
%  given.  We get the family number for the epochs, and than the epoch #s
%  for all of the steps.  This means that we will have all the families
%  that include the steps listed, including families that have more than
%  these steps, and all steps for each family will be included.  We can
%  do a check to see if the user wants the epochs that include more than
%  these epochs or not. 
%  Now we get the index numbers of the cue specified. 

%  Cueepochs should be index numbers, not the actual epoch numbers
if isempty(cue)
  newepochnums = catstepepochs;
else
  cueind;
  allcues = famgp.FamilyCues(cueind,:);
  cue;
  cueepochs = find(allcues==cue);
  cuelength = length(cueepochs);
  ind = 0;
  %  Figure out which index numbers are the same.
  for i = 1:cuelength
	 newlistepochs = find(catstepepochs==cueepochs(i));
	 if ~isempty(newlistepochs)
		ind = ind + 1;
		newepochs(ind) = newlistepochs;
	 end
  end
  %  newepochs are the index numbers in catstepepochs that are the same 
  %  as newcueepochs so we need to get the numbers from catstepepochs to find 
  %  out what the actual index numbers are.
  newepochnums = catstepepochs(newepochs);
end
%  newepochnums is a list of the epochs we are using.  How do we
%  find out which family each epoch is in?  FamilyFlag has all of the
%  epochs divided into families, so we get the index number of the epochs
%  we want from newepochnums, and check that index number in the
%  FamilyFlag list, we need to also check that index number with allepochs 
%  to get the epoch numbers.  But first we need to get rid of
%  epochnumbers that are not in the range specified by the inputs
%  lowrange and highrange.  Since we use
%  the index numbers to move between allepochs, allsteps, and allflags,
%  we need to leave the length of allepochs the same, but somehow get rid
%  of the epochnumbers we don't want.  Do this by using NaN (zero is no
%  good, since epoch #'s start with zero). We also change these index
%  numbers in the allflags, so that when we get family #s in the end, we
%  only get the families that we want.
if lowrange~=0
  rmlowepochs = find(allepochs<lowrange);
  allepochs(rmlowepochs) = NaN;
  allflags(rmlowepochs) = NaN;
end
if ~isempty(highrange)
  rmhighepochs = find(allepochs>highrange);
  allepochs(rmhighepochs) = NaN;
  allflags(rmhighepochs) = NaN;
end
newepochslength = length(newepochnums);
epochfamilies = allflags(newepochnums);
epochfamilies = epochfamilies(find(~isnan(epochfamilies)));
realnumbers = allepochs(newepochnums);
realnumbers = realnumbers(find(~isnan(realnumbers)));
%  Ok, now we should have a list of family#s that corresponds to the same
%  epochnumbers.
%  Right now returns the actual epochnumbers, not an index number to the
%  epochdata in CellInfo, so we leave the user to deal with the pesky
%  zero epoch.
%  Have to deal with that pesky zero epoch
%  realnumbers = realnumbers + 1;
