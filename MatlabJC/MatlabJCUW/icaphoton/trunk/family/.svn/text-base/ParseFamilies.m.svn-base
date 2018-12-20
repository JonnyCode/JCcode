function [familydata,familynumbers] = ParseFamilies(CellInfo,realnumbers, ...
																  epochfamilies) 
%  [familydata] = ParseFamilies(CellInfo,realnumbers,epochfamilies)
%  Gives back data and epoch numbers parsed into families, when given a list of
%  epochnumbers and their corresponding familynumbers.  Assumes the
%  epochnumbers are the actual epoch numbers, and returns actual epoch
%  numbers.

%  Created MKMK Oct 2001

%  Since the EpochData has epoch zero as the 1st index, must add 1 to
%  realnumbers.
realnumbers = realnumbers + 1;
famnumbers = unique(epochfamilies);
familylength = length(famnumbers);
familydata = cell(1,familylength);
familyflag = 1;
j = 1;
firstnumber = j;
% We want to get the data so the max is a pos. value.  Test to see if
% currently pos. or neg.  Use two epochs to be sure.
data = CellInfo.EpochData.Data;
testepoch = max(data{realnumbers(1)});
testepochb = min(data{realnumbers(1)});
testepoch2 = max(data{realnumbers(2)});
testepoch2b = min(data{realnumbers(2)});
if abs(testepoch)>abs(testepochb) & abs(testepoch2)>abs(testepoch2b)
  flag = 0;
  % data is already set up to have a positive max.  Do nothing.
elseif abs(testepoch)<abs(testepochb) & abs(testepoch2)<abs(testepoch2b)
  % data has a negative max.  Multiply by -1
  flag = 1;
else
  test2 = input('Inconclusive test. Is data max. pos. or neg.? p/n:', ...
					 's');
  if strcmp(test2,'p')
	 flag = 0;
  elseif strcmp(test2,'n')
	 flag = 1;
  else
	 error('inconclusive test')
  end
end
flag;
%  Goes through the list of epochfamilies, as long as the numbers are the
%  same, each index number is in the same family, if the next number is
%  not the same, than we start another family.  
for i = 1:familylength
	while familyflag == 1
		if j + 1 > length(epochfamilies)
			familyflag = 0;
			lastnumber = j;
		elseif epochfamilies(j+1)==epochfamilies(j)
			j = j + 1;
		elseif epochfamilies(j+1)~=epochfamilies(j)
			lastnumber = j;
			familyflag = 0;
		end
	end
	if lastnumber > firstnumber		
		familyrange = realnumbers(firstnumber):realnumbers(lastnumber);
		tempdata = cat(1,CellInfo.EpochData.Data{familyrange});
		if flag == 1
		  tempdata = tempdata*(-1);
		end
		familydata{i} = tempdata;
		familynumbers{i} = familyrange-1;
		j = j + 1;
		familyflag = 1;
		firstnumber = j;
	end
end
