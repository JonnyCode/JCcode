function [ReturnCondition, OutputComp] = FindResponses(ResponseStats, FlashCondition, CellParameters, NumConds, LowerAmpLimit, UpperAmpLimit)
%
%  function [ReturnCondition, OutputComp] = FindResponses(ResponseStats, FlashCondition, CellParameters, NumConds, LowerAmpLimit, UpperAmpLimit)
%
%  Find responses in a given window set by the loweramplimit and the upperamplimit
%  Created:  FMR ?
%  Revised  GDF 12/8/01
%			GDF: Changed so that it also outputs (via OutputComp) the number of responses isolated under each condition


EpochLength = length(FlashCondition(1).EpochData.Data(1,:));
ReturnCondition.SearchCrit = FlashCondition(1).SearchCrit;
ReturnCondition.SearchPara = FlashCondition(1).SearchPara;
ReturnCondition.PlotPref = FlashCondition(1).PlotPref;
ReturnCondition.DecimatePts = FlashCondition(1).DecimatePts;
ReturnCondition.Label = [];
ReturnCondition.UserInfo = [];
ReturnCondition.ScaleFactorIndex = FlashCondition(1).ScaleFactorIndex;

NumResponses = 0;
PrevNumResponses = 0;

for cond = 1:NumConds
	NumEpochs = length(ResponseStats(cond).EpochNumbers);	
	for epoch = 1:NumEpochs
		% if this epoch in response stats has desired amplitude
		if ((ResponseStats(cond).RegressAmp(epoch) > LowerAmpLimit) & (ResponseStats(cond).RegressAmp(epoch) < UpperAmpLimit))
			% find appropriate data in FlashCondition structure
			for cnt = 1:length(FlashCondition(cond).EpochNumbers)
				if (ResponseStats(cond).EpochNumbers(epoch) == FlashCondition(cond).EpochNumbers(cnt))
					NumResponses = NumResponses + 1;
					ReturnCondition.EpochData.Data(NumResponses, 1:EpochLength) = FlashCondition(cond).EpochData.Data(cnt, :);
					ReturnCondition.EpochData.Offset(NumResponses) = FlashCondition(cond).EpochData.Offset(cnt);
					ReturnCondition.ExcludeEpochs(NumResponses) = 0;
					ReturnCondition.EpochNumbers(NumResponses) = FlashCondition(cond).EpochNumbers(cnt);
					ReturnCondition.EpochTime(NumResponses) = FlashCondition(cond).EpochTime(cnt);
					ReturnCondition.UserInfo.StimulusAmp(NumResponses) = FlashCondition(cond).UserInfo.StimulusAmp;
					break;
				end
			end
		end
	end
	NumIsoms = CellParameters.CollectingArea * FlashCondition(cond).UserInfo.StimulusAmp;
	ExpectedNumResponses = exp(-NumIsoms) * NumIsoms * NumEpochs;
	ExpectedNumFails = exp(-NumIsoms) * NumEpochs;
	fprintf(1, 'Found %d responses; Expect %d singles, %d failures\n', NumResponses - PrevNumResponses, ExpectedNumResponses, ExpectedNumFails);
	OutputComp(cond).Found = NumResponses - PrevNumResponses;
	OutputComp(cond).ExpectedHits = ExpectedNumResponses;
	OutputComp(cond).ExpectedMisses = ExpectedNumFails;
	PrevNumResponses = NumResponses;
end
