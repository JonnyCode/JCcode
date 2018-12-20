function ResponseStats = SimBipFromComb(CombParameters, NumTrials, StartFlash, NumFlashConds, PoolSize, MeanCumulativeGauss, SDCumulativeGauss)

MeanResp(1:NumFlashConds) = 0;
VarResp(1:NumFlashConds) = 0;
MeanSumResp(1:NumFlashConds) = 0;
VarSumResp(1:NumFlashConds) = 0;
FlashStrength(1:NumFlashConds) = 0;

for cond = 1:NumFlashConds
	CombParameters.RateA = StartFlash * 2^(cond-1);
	CombParameters.RateB = StartFlash * 2^(cond-1);
	for trial = 1:NumTrials
		[respA, respB] = SampleFromCombDistribution2(CombParameters, PoolSize);
		sumresp = normcdf(mean(respA), MeanCumulativeGauss, SDCumulativeGauss) * sum(respA);
		respA = normcdf(respA, MeanCumulativeGauss, SDCumulativeGauss) .* respA;
		MeanResp(cond) = MeanResp(cond) + sum(respA);
		MeanSumResp(cond) = MeanSumResp(cond) + sumresp;
	end
	MeanResp(cond) = MeanResp(cond) / NumTrials;
	MeanSumResp(cond) = MeanSumResp(cond) / NumTrials;
	for trial = 1:NumTrials
		[respA, respB] = SampleFromCombDistribution2(CombParameters, PoolSize);
		sumresp = normcdf(mean(respA), MeanCumulativeGauss, SDCumulativeGauss) * sum(respA);
		respA = normcdf(respA, MeanCumulativeGauss, SDCumulativeGauss) .* respA;
		VarResp(cond) = VarResp(cond) + (MeanResp(cond) - sum(respA))^2;
		VarSumResp(cond) = VarSumResp(cond) + (MeanSumResp(cond) - sumresp)^2;
	end
	VarResp(cond) = VarResp(cond) / NumTrials;
	VarSumResp(cond) = VarSumResp(cond) / NumTrials;
	FlashStrength(cond) = CombParameters.RateA;
	fprintf(1, 'Flash = %d, Mean = %d, Var = %d, MeanSum = %d, VarSum = %d\n', FlashStrength(cond), MeanResp(cond), VarResp(cond), MeanSumResp(cond), VarSumResp(cond));
end

ResponseStats.MeanResp = MeanResp;
ResponseStats.VarResp = VarResp;
ResponseStats.MeanSumResp = MeanSumResp;
ResponseStats.VarSumResp = VarSumResp;
ResponseStats.FlashStrength = FlashStrength;
