function RespProjection = CalcRespProjsandHistsSing(CellInfo, FlashConditionA, NoiseCondition, FreqCutoff, RespLength, shuffles)
%
%	 RespProjection = CalcRespProjsandHistsSing(CellIfno, FlashConditionA, NoiseCondition, FreqCutoff, RespLength, shuffles)
%
%	Generate the histograms from the response projections
%	Created:  FMR ?
%	Revised:  GDF 10/08/01
%			- made it compatible with the new data format
%				GDF 10/01
%			- This version is designed to handle only a one flash condition and a noise condition
%				useful for isolating singles

NumFlashConds = length(FlashConditionA);

global NumHistBins;

% calculate projections of responses through discriminant
for cond=1:NumFlashConds
	Projections = Calc2AFCRespProjSingles(CellInfo, FlashConditionA, NoiseCondition, cond, FreqCutoff, RespLength, shuffles);
	NumResponses = length(Projections.RespAProjection);
	RespProjection(cond).RespAProjection = Projections.RespAProjection;
	RespProjection(cond).NoiseProjection = Projections.NoiseProjection;
end

% calculate histograms of response projections
clear Respx;
Respx(1:NumHistBins) = 1:NumHistBins;
Respx = Respx - NumHistBins/2;
Respx = (floor(max(RespProjection(NumFlashConds).RespAProjection)) + 1) * 2 * Respx / NumHistBins;
Respx = Respx';
for cond = 1:NumFlashConds
	temp = RespProjection(cond).RespAProjection;
	[HistA, Respx] = hist(temp', Respx);
	temp = RespProjection(cond).NoiseProjection;
	[NoiseHist, Respx] = hist(temp', Respx);
	RespProjection(cond).RespAHist = HistA;
	RespProjection(cond).NoiseHist = NoiseHist;
	RespProjection(cond).Respx = Respx;
	RespProjection(cond).NumHistBins = NumHistBins;
end


