function ReturnedCondition = ExcludeGlitches(CellInfo, Condition)
%
%  ReturnedCondition = ExcludeGlitches(CellInfo, Condition)
%
%  generates an amplitude histogram of epochs given in the condition, and 
%  then alows you to exclude epochs that have a glitch that
%  exceeds some threshold.
% 
%  Created ? GDF

ReturnedCondition = Condition;
[EpochData, Offset] = GetEpochData(CellInfo, Condition);
[NumIterations, EpochSize] = size(EpochData);
RectifiedEpochs = abs(EpochData);
Maximums = max(RectifiedEpochs, [], 2);
x = 0:.5:20;
hist(Maximums,x)

% Pick a criteria such that every epoch that you want to look at has a 
% maximum value greater than the criteria.  It will plot these epochs
% and alow you to choose which ones you want to keep and which ones 
% you want to exclude from the analysis.  
% A flag = 1 will be placed in the ExcludeEpoch vector in the structure for
% those epochs that you want to exclude from the analysis

GoOn = input('Would you like to set a criteria to sift through the data ? y or n \n','s');
c = 1;
for index = 1:c
	if strcmp(GoOn,'y')
		Criteria = input('What Criteria Would you like? \n');
		for cnt = 1:NumIterations
			if Maximums(cnt) > Criteria
				ReturnedCondition.ExcludeEpochs(cnt) = 1;
			end
		end
		hold off;
		a = 1;
		for cnt = 1:NumIterations
			if ReturnedCondition.ExcludeEpochs(cnt)
				plot(EpochData(cnt,:));
				Flag = input('is this epoch OK ? y or n \n','s');
				if strcmp(Flag,'y')
					ReturnedCondition.ExcludeEpochs(cnt) = 0;
				end
			end	
		end			
	end
end	
