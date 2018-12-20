function [ResponseStats, NoiseStats] = CalculateResponseStats(FlashCondition, NoiseCondition, CellParameters, ResponseSize, WindowSize)
%
% [ResponseStats, NoiseStats] = CalculateResponseStats(FlashCondition, NoiseCondition, CellParameters, ResponseSize, WindowSize)
%
%	Calculate statistics of single photon responses and dark noise.  Returns peak amplitude, area and projection
%	for each response in condition structure and noise structure.
%
% Created 1/01 FMR
% This version has some recent edits by Fred
% Revised: GDF 12/07/01
%


NumEpochs = length(FlashCondition.ExcludeEpochs);	
MaxAmplitude(1:NumEpochs) = 0;	
Area(1:NumEpochs) = 0; 	
RegressAmp(1:NumEpochs) = 0; 	
EpochNumbers(1:NumEpochs) = 0;
[Max, MaxLoc] = min(FlashCondition.AverageResponse);
SampInterv = FindSearchPara(FlashCondition, 'SampInterv');
SamplingInterval = FlashCondition.DecimatePts/SampInterv;
EpochPts = length(FlashCondition.EpochData.Data(1,:));
%PrePts = FlashCondition.UserInfo.StimulusStartTime(1) / FlashCondition.SamplingInterval;
PrePts = 1;

% template for response projection calculations
Template(1:ResponseSize) = CellParameters.ElementaryResponse(1:ResponseSize);
Template = Template / sum(Template .* Template);

% loop through FlashCondition structure
GoodValues = 0;
for cnt = 1:NumEpochs
	if (FlashCondition.ExcludeEpochs(cnt) == 0)
		GoodValues = GoodValues + 1;
		resp(1:ResponseSize) = FlashCondition.EpochData.Data(cnt, 1:ResponseSize);
		MaxAmplitude(GoodValues) = -sum(resp(MaxLoc-WindowSize/2:MaxLoc+WindowSize/2)) / WindowSize;
		RegressAmp(GoodValues) = sum(resp .* Template);
		resp2(1:EpochPts-PrePts+1) = FlashCondition.EpochData.Data(cnt, PrePts:EpochPts);
		Area(GoodValues) = -sum(resp) * SamplingInterval;;
		EpochNumbers(GoodValues) = FlashCondition.EpochNumbers(cnt);
	end
end

% load into return structure
ResponseStats.MaxAmplitude(1:GoodValues) = MaxAmplitude(1:GoodValues);	
ResponseStats.Area(1:GoodValues) = Area(1:GoodValues); 	
ResponseStats.RegressAmp(1:GoodValues) = RegressAmp(1:GoodValues); 	
ResponseStats.EpochNumbers(1:GoodValues) = EpochNumbers(1:GoodValues);

% loop through NoiseCondition structure
GoodValues = 0;
NumNoiseEpochs = length(NoiseCondition.ExcludeEpochs);
for cnt = 1:NumNoiseEpochs
	if (NoiseCondition.ExcludeEpochs(cnt) == 0)
		GoodValues = GoodValues + 1;
		resp(1:ResponseSize) = NoiseCondition.EpochData.Data(cnt, 1:ResponseSize);
		Area(GoodValues) = -sum(resp);
		MaxAmplitude(GoodValues) = -sum(resp(MaxLoc-WindowSize/2:MaxLoc+WindowSize/2)) / WindowSize;
		RegressAmp(GoodValues) = sum(resp .* Template);
		resp2(1:EpochPts-PrePts+1) = NoiseCondition.EpochData.Data(cnt, PrePts:EpochPts);
		Area(GoodValues) = -sum(resp) * SamplingInterval;
		%EpochNumbers(GoodValues) = FlashCondition.ExcludeEpochs(cnt);
	end
end

NoiseStats.MaxAmplitude(1:GoodValues) = MaxAmplitude(1:GoodValues);	
NoiseStats.Area(1:GoodValues) = Area(1:GoodValues); 	
NoiseStats.RegressAmp(1:GoodValues) = RegressAmp(1:GoodValues); 	
%NoiseStats.EpochNumbers(1:GoodValues) = EpochNumbers(1:GoodValues);
