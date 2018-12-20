%script created to convert an old-style CellInfo in which the fields
% OutputScaleFactor and NDFConfiguration contain strings to ones 
% that have numbers.  Will go through every file in the current 
% directory and check to see if the mat files are ones that 
% contain the old-style fields.  Should not do anything to either
% newer CellInfo files or files that are not mat or do not contain 
% CellInfodircont = dir('*.mat');
for i = 1:length(dircont)
  filename = dircont(i).name
  load(filename)
  s = who;
  s(1)
  if strcmp(s(1),'CellInfo')
	  temp = CellInfo.OutputScaleFactor;
	  if ~isnumeric(temp)
		 temp = char(temp);	
		 temp = str2num(temp);
		 CellInfo.OutputScaleFactor = temp;
	  end
	  temp = CellInfo.NDFConfiguration;
	  if ~isnumeric(temp)
		 temp = char(temp);
		 temp = str2num(temp);
		 CellInfo.NDFConfiguration = temp;
	  end
	  save(filename,'CellInfo')
	  clear CellInfo
  else
	  disp('Did not change file')
	  clear s
  end
end


