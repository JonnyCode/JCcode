function RespProjection = CalcRespProjsandHists2(CellInfo, FlashConditionA, FlashConditionB, NoiseCondition, FreqCutoff, RespLength, shuffles)
%
%	 RespProjection = CalcRespProjsandHists2(CellIfno, FlashConditionA, FlashConditionB, NoiseCondition, FreqCutoff, RespLength, shuffles)
%
%	Generate the histograms from the response projections
%	Created:  FMR ?
%	Revised:  GDF 10/08/01
%			- made it compatible with the new data format

NumFlashConds = length(FlashConditionA);

global NumHistBins;

% calculate projections of responses through discriminant
for cond=1:NumFlashConds
	Projections = Calculate2AFCRespProjection2(CellInfo, FlashConditionA, FlashConditionB, NoiseCondition, cond, FreqCutoff, RespLength, shuffles);
	NumResponses = length(Projections.RespAProjection);
	RespProjection(cond).RespAProjection = Projections.RespAProjection;
	RespProjection(cond).RespBProjection = Projections.RespBProjection;
	RespProjection(cond).NoiseProjection = Projections.NoiseProjection;
end

% calculate histograms of response projections
clear Respx;
Respx(1:NumHistBins) = 1:NumHistBins;
Respx = Respx - NumHistBins/2;
Respx = (floor(max(abs(RespProjection(NumFlashConds).RespAProjection))) + 1) * 2 * Respx / NumHistBins;
max(Respx)
Respx = Respx';
for cond = 1:NumFlashConds
	temp = RespProjection(cond).RespAProjection;
	[HistA, Respx] = hist(temp', Respx);
	temp = RespProjection(cond).RespBProjection;
	[HistB, Respx] = hist(temp', Respx);
	temp = RespProjection(cond).NoiseProjection;
	[NoiseHist, Respx] = hist(temp', Respx);
	RespProjection(cond).RespAHist = HistA;
	RespProjection(cond).RespBHist = HistB;
	RespProjection(cond).NoiseHist = NoiseHist;
	RespProjection(cond).Respx = Respx;
	RespProjection(cond).NumHistBins = NumHistBins;
end


