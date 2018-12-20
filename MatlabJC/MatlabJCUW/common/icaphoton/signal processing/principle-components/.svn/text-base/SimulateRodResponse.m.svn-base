function ReturnCondition = SimulateRodResponse(Condition, CellParameters, PCAProjections, NumResponses, FlashTime, FlashStrength, PoissonFlag, VarScaleFact)

[NumEpochs, EpochPts] = size(Condition.EpochData.Data);
dat(1:EpochPts) = 0;
ReturnCondition.EpochNumbers(1:NumResponses) = 0;
ReturnCondition.ExcludeEpochs(1:NumResponses) = 0;
ReturnCondition.EpochData.Offset(1:NumResponses) = 0;
ReturnCondition.SearchCrit = Condition.SearchCrit;
ReturnCondition.SearchPara = Condition.SearchPara;
ReturnCondition.DecimatePts = Condition.DecimatePts;
ReturnCondition.UserInfo.StimulusStartTime(1:NumResponses) = 0;

%Time = 1:EpochPts;
%SampInterv = FindSearchPara(Condition, 'SampInterv');
%SamplingInterval = Condition.DecimatePts / SampInterv;
%Time = Time * SamplingInterval;
Time = Condition.UserInfo.TimePoints;
SamplingInterval = Time(2);


EigVec = PCAProjections.EigVec;
ThermalProb = CellParameters.ThermalRate * SamplingInterval;
ContinuousNoiseFilter = (Time/CellParameters.ContinuousNoiseTimeConstant).^(CellParameters.ContinuousNoisePow-1) .* exp(-Time/CellParameters.ContinuousNoiseTimeConstant);
ContinuousNoiseFilter = fft(ContinuousNoiseFilter);
NumPhotons = 0;

for resp = 1:NumResponses
	dat(1:EpochPts) = 0;
	if (PoissonFlag == 1)
		if (FlashStrength > 0)
			NumPhotons = poissrnd(FlashStrength);
		end
	else
		NumPhotons = FlashStrength;
	end
	if (NumPhotons > 0)
		for photon=1:NumPhotons
			for comp = 1:size(PCAProjections.EigVec, 2)
				StDev = (abs((PCAProjections.SinglesFitCoef(comp, 3)^2 - PCAProjections.FailuresFitCoef(comp, 3)^2) * VarScaleFact))^0.5;
				MeanProj = PCAProjections.SinglesFitCoef(comp, 2) - PCAProjections.FailuresFitCoef(comp, 2);
				dat = dat + (normrnd(MeanProj, StDev) * EigVec(:, comp)');
			end
		end
		ReturnCondition.EpochData.Data(resp, 1:EpochPts) = dat(1:EpochPts);
	else
		ReturnCondition.EpochData.Data(resp, 1:EpochPts) = 0;
	end
	
	% add thermals
	Thermals = poissrnd(ThermalProb, EpochPts, 1);
	if (sum(Thermals) > 0)
		for pnt = 1:EpochPts
			if (Thermals(pnt) == 1)
				dat = 0;
				for comp = 1:size(PCAProjections.EigVec, 2)
					StDev = ((PCAProjections.SinglesFitCoef(comp, 3)^2 - PCAProjections.FailuresFitCoef(comp, 3)^2) * VarScaleFact)^0.5;
					MeanProj = PCAProjections.SinglesFitCoef(comp, 2) - PCAProjections.FailuresFitCoef(comp, 2);
					dat = dat + normrnd(MeanProj, StDev) * EigVec(:, comp);
				end
				dat = dat';
%				fprintf(1, 'Thermal at %d (%d)\n', pnt, resp);
				ReturnCondition.EpochData.Data(resp, pnt:EpochPts) = ReturnCondition.EpochData.Data(resp, pnt:EpochPts) + dat(1:EpochPts-pnt+1);
			end		
		end
	end
	
	% add continuous noise
	ContinuousNoise = NORMRND(0, 1, 1, EpochPts);
	ContinuousNoise = fft(ContinuousNoise);
	ContinuousNoise = ContinuousNoise .* ContinuousNoiseFilter;
	ContinuousNoise = real(ifft(ContinuousNoise));	
	ContinuousNoise = ContinuousNoise * CellParameters.ContinuousNoiseSD / STD(ContinuousNoise);
	ReturnCondition.EpochData.Data(resp, 1:EpochPts) = ReturnCondition.EpochData.Data(resp, 1:EpochPts) + ContinuousNoise(1:EpochPts);
	
	flag = 1;
	if flag
		%add more continuous noise
		ContinuousNoise = NORMRND(0, 1, 1, EpochPts);
		ContinuousNoise = fft(ContinuousNoise);
		ContinuousNoise = ContinuousNoise .* ContinuousNoiseFilter;
		ContinuousNoise = real(ifft(ContinuousNoise));	
		ContinuousNoise = ContinuousNoise * CellParameters.ContinuousNoiseSD / STD(ContinuousNoise);
		ExtraData(resp, 1:EpochPts) = ContinuousNoise(1:EpochPts);	
	end
end

if flag
	ReturnCondition.EpochData.Data = cat(2, ReturnCondition.EpochData.Data, ExtraData);
end
