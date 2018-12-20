function ReturnedCellInfo = AdjustGain(CellInfo)
fp = CellInfo.CellFile;
[fileptr, ecode] = ITCInitializeAnalysis(500000, fp);
if ecode~=0
	error
end
data = CellInfo.EpochData.Data;
for i = 1:length(data)
	if ~isempty(data(i))
		[gain0, ecode] = ITCGetAmpGain(i-1, 0, fileptr);
		if ecode~=0
			error
		end
		[gain1, ecode] = ITCGetAmpGain(i-1, 1, fileptr);
		if ecode~=0
			error
		end
		gain = gain0/gain1;
		data{i} = data{i}*gain;
	end
end
ITCFinishAnalysis(fileptr);
CellInfo.EpochData.Data = data;
ReturnedCellInfo = CellInfo;
