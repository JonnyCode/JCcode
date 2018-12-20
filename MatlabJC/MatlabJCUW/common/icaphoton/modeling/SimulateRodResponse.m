function ReturnCondition = SimulateRodResponse(Condition, CellParameters, PCAProjections, NumResponses, FlashTime, FlashStrength, PoissonFlag, VarScaleFact)

dat(1:Condition.EpochPts) = 0;
ReturnCondition.SamplingInterval = Condition.SamplingInterval;
ReturnCondition.EpochPts = Condition.EpochPts;
ReturnCondition.EpochNumbers(1:NumResponses) = 0;
ReturnCondition.ExcludeEpochs(1:NumResponses) = 0;

Time = 1:Condition.EpochPts;
Time = Time * Condition.SamplingInterval;

EigVec = PCAProjections.EigVec;
ThermalProb = CellParameters.ThermalRate * Condition.SamplingInterval;
ContinuousNoiseFilter = (Time/CellParameters.ContinuousNoiseTimeConstant).^(CellParameters.ContinuousNoisePow-1) .* exp(-Time/CellParameters.ContinuousNoiseTimeConstant);
ContinuousNoiseFilter = fft(ContinuousNoiseFilter);
NumPhotons = 0;

% generate simulated response
for resp = 1:NumResponses
	dat = 0;
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
				StDev = ((PCAProjections.SinglesFitCoef(comp, 3)^2 - PCAProjections.FailuresFitCoef(comp, 3)^2) * VarScaleFact)^0.5;
				MeanProj = PCAProjections.SinglesFitCoef(comp, 2) - PCAProjections.FailuresFitCoef(comp, 2);
				dat = dat + normrnd(MeanProj, StDev) * EigVec(:, comp);
			end
		end
		ReturnCondition.EpochData(resp, FlashTime:Condition.EpochPts) = dat(1:Condition.EpochPts-FlashTime+1);
	else
		ReturnCondition.EpochData(resp, 1:Condition.EpochPts) = 0;
	end
	
	% add thermals
	Thermals = poissrnd(ThermalProb, Condition.EpochPts, 1);
	if (sum(Thermals) > 0)
		for pnt = 1:Condition.EpochPts
			if (Thermals(pnt) == 1)
				dat = 0;
				for comp = 1:size(PCAProjections.EigVec, 2)
					StDev = ((PCAProjections.SinglesFitCoef(comp, 3)^2 - PCAProjections.FailuresFitCoef(comp, 3)^2) * VarScaleFact)^0.5;
					MeanProj = PCAProjections.SinglesFitCoef(comp, 2) - PCAProjections.FailuresFitCoef(comp, 2);
					dat = dat + normrnd(MeanProj, StDev) * EigVec(:, comp);
				end
				dat = dat';
%				fprintf(1, 'Thermal at %d (%d)\n', pnt, resp);
				ReturnCondition.EpochData(resp, pnt:Condition.EpochPts) = ReturnCondition.EpochData(resp, pnt:Condition.EpochPts) + dat(1:Condition.EpochPts-pnt+1);
			end		
		end
	end
	
	% add continuous noise
	ContinuousNoise = NORMRND(0, 1, 1, Condition.EpochPts);
	ContinuousNoise = fft(ContinuousNoise);
	ContinuousNoise = ContinuousNoise .* ContinuousNoiseFilter;
	ContinuousNoise = real(ifft(ContinuousNoise));	
	ContinuousNoise = ContinuousNoise * CellParameters.ContinuousNoiseSD / STD(ContinuousNoise);
	ReturnCondition.EpochData(resp, 1:Condition.EpochPts) = ReturnCondition.EpochData(resp, 1:Condition.EpochPts) + ContinuousNoise(1:Condition.EpochPts);
end


