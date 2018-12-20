function [fammaxes,famareamaxes] = FamilyMax(CellInfo,famnum)
% function [fammaxes,famareamaxes] = FamilyMax(CellInfo,famnum)
% Returns a max for all the maxes of epochs in this
% familycondition. Uesful for normalizing.

%Created MKMK Oct 2001

famcell = CellInfo.FamilyCondition(famnum);
EpochNumbers = famcell.EpochNumbers;
EpochData = CellInfo.EpochData.Data;
FamilyStep = famcell.FamilyStep;
FamilyFlag = famcell.FamilyFlag;

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
epochnums = EpochNumbers + 1;
numslen = length(epochnums);
allmax = zeros(1,numslen);
if mdefault == 1
  for j = 1:numslen
	 allmax(j) = max(EpochData{epochnums(j)});
  end
elseif mdefault == 0
  for j = 1:numslen
	 allmax(j) = abs(min(EpochData{epochnums(j)}));
  end
end
famaxes = max(allmax);
for j = 1:numslen
  temp = EpochData{epochnums(j)};
  h = 0.1;
  temp(1) = temp(1)/2;
  temp(end) = temp(end)/2;
  I = h*sum(temp);
  I = abs(I);
  allarea(j) = I;
end
famareamaxes = max(allarea);
