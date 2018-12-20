function [ReturnedFlashConditionA, ReturnedFlashConditionB] = AddThermals(FlashConditionA, FlashConditionB, CellParameters, ThermRate);
% AddThermals
%
%	Add synthetic single photon responses to data in two conditions (A and B).  Single photon response is estimated from 
%   average responses in Condition A, scaled by flash strength and collecting area.  Identical responses are then added
% 	at random times at the rate specified.

SamplingInterval = FindSearchPara(FlashConditionA(1), 'SampInterv') * FlashConditionA(1).DecimatePts * 1e-6; 
StimulusDuration = FindSearchPara(FlashConditionA(1), 'StimDur') * SamplingInterval;
EpochPts =  size(FlashConditionA(1).EpochData.Data, 2);
ThermProb = ThermRate * EpochPts * SamplingInterval;
CombParameters = CellParameters.CombParameters(10);

% get flash strengths in photoisomerizations and estimate mean single photon response
for cond = 1:length(FlashConditionA)
	FlashStrength(cond) = FlashConditionA(cond).UserInfo.StimulusAmp * StimulusDuration * CellParameters.CollectingArea;
	if (cond == 1)
		ElementaryResponse = FlashConditionA(cond).AverageResponse / FlashStrength(cond);
	else
		ElementaryResponse = ElementaryResponse + FlashConditionA(cond).AverageResponse / FlashStrength(cond);
	end
end
ElementaryResponse = ElementaryResponse / length(FlashConditionA);

for cond = 1:length(FlashConditionA)
	fprintf(1, 'condition = %d\n', cond);
	for epoch = 1:size(FlashConditionA(cond).EpochData.Data, 1)
		if (unifrnd(0, 1) < ThermProb)
			RotatedElementaryResponse = RotateVector(round(unifrnd(0, EpochPts)), ElementaryResponse);
			FlashConditionA(cond).EpochData.Data(epoch, :) = FlashConditionA(cond).EpochData.Data(epoch, :) + normrnd(1, abs(CombParameters.UnitarySDA) / CombParameters.UnitaryAmpA) * RotatedElementaryResponse;
		end
	end
	for epoch = 1:size(FlashConditionB(cond).EpochData.Data, 1)
		if (unifrnd(0, 1) < ThermProb)
			RotatedElementaryResponse = RotateVector(round(unifrnd(0, EpochPts)), ElementaryResponse);
			FlashConditionB(cond).EpochData.Data(epoch, :) = FlashConditionB(cond).EpochData.Data(epoch, :) + normrnd(1, abs(CombParameters.UnitarySDA) / CombParameters.UnitaryAmpA) * RotatedElementaryResponse;
		end
	end
	ReturnedFlashConditionA(cond) = FlashConditionA(cond);
	ReturnedFlashConditionB(cond) = FlashConditionB(cond);
end
