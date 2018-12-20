% SimulateBipolarResponse
%
%	Simulate response of bipolar cell receiving input
%	from multiple rods.  Rod parameters as described in CellParameters,
%	which should be as defined from rod analysis routines.
%
%	Created 7/01 FMR
%	Revisions:
%	- 9/01 FMR: revised to include flag to set whether nonlinearity applied to each rod signal (Location = 0)
%			    or to summed of all rod singals (Location = 1)

function ResponseStats = SimulateBipolarResponse(Condition, CellParameters, NumRods, NumTrials, Mean, SD, Location)

PoissRate = Condition.StimulusAmp * CellParameters.CollectingArea * CellParameters.SamplingInterval;

% generate response for pool of rods based on rod parameters 

CombParameters = CellParameters.CombParameters;
PrePts = Condition.StimulusStartTime(1) / CellParameters.SamplingInterval;
StmPts = Condition.StimulusDuration / CellParameters.SamplingInterval;
EpochPts = round(Condition.EpochPts * Condition.SamplingInterval / CellParameters.SamplingInterval);
Response(NumTrials, 1:EpochPts) = 0;

for resp = 1:NumTrials
	for rod = 1:NumRods
		clear dat;
		dat(1:EpochPts) = 1e-10;
		dat(PrePts:PrePts+StmPts) = PoissRate;
		dat = poissrnd(dat);
		if (sum(dat) > 0)
			dat = conv(dat, CellParameters.ElementaryResponse);
			temp(1:EpochPts) = -dat(1:EpochPts) * normrnd(1, CombParameters.UnitarySDA / CombParameters.UnitaryAmpA);
			if (Location == 0)
				temp(1:EpochPts) = -normcdf(temp, Mean, SD) .* temp;
			end
			Response(resp, :) = Response(resp, :) + temp(1:EpochPts);
		end
	end
end

AverageResponse(1:EpochPts) = 0;
if (Location == 1)
	for resp = 1:NumTrials
		Response(resp, :) = -normcdf(Response(resp, :), Mean, SD) .* Response(resp, :);
	end
end

for resp = 1:NumTrials
	AverageResponse = AverageResponse + Response(resp, :);
end
ResponseStats.AverageResponse = AverageResponse / NumTrials;

VarianceResponse(1:EpochPts) = 0;
for resp = 1:NumTrials
	VarianceResponse = VarianceResponse + (ResponseStats.AverageResponse - Response(resp, :)).^2;
end
ResponseStats.VarianceResponse = VarianceResponse / NumTrials;

[ResponseStats.PeakHist, ResponseStats.PeakHistx] = hist(min(Response'), 40);

ResponseStats.EpochData = Response;
