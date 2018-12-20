function CheckStability(CellInfo, Condition)
%
%	function CheckStability(CellInfo, Condition)
%
% CheckStability.m
%	Check cell stability over time by comparing average response
%	in chunks of experiment with overall average.  Divides data into 
%	4 sections and compares average of each section with overall average.
%
%	Created 10/00 FMR
%	Edited 11/00 GDF
%		Altered to make it compatible with new data format

%Label = Condition.Label{1};
%fprintf(1, 'The condition is  %s\n', Label);
[EpochData, Offset] = GetGoodEpochData(CellInfo, Condition);
[Number, EpochPts] = size(EpochData)
tme = 1:EpochPts;
Ave = 1:EpochPts;
for segment = 0:3
	GoodEpochs = 0;
	Ave = 0;
	TempEpochs = EpochData((segment*floor(Number/4)+1):((segment+1)*floor(Number/4)),1:EpochPts);
	Ave = mean(TempEpochs);
	plot(tme, Ave, tme, Condition.AverageResponse);
	ExcludeFlag = input('Continue?');
	clear TempEpochs
end




