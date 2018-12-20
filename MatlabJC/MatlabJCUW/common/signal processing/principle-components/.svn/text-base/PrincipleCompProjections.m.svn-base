function PCAProjections = PrincipleCompProjections(SinglesCondition, FailuresCondition, Components)
%
%	function PCAProjections = PrincipleCompProjections(SinglesCondition, FailuresCondition, Components)
%
% CalculateResponseStats
%
NumHistBins = 30;
[NumEpochs, EpochPts] = size(SinglesCondition.EpochData.Data);	
NumComp = size(Components, 2);

for cnt = 1:NumEpochs
	resp(1:EpochPts) = SinglesCondition.EpochData.Data(cnt, 1:EpochPts);
	for comp = 1:NumComp
		pcomp(1:EpochPts) = Components(1:EpochPts, comp);
		PCAProjections.SinglesProj(cnt, comp) = sum(resp .* pcomp);
	end
end

NumFailEpochs = length(FailuresCondition.ExcludeEpochs);	
for cnt = 1:NumFailEpochs
	resp(1:EpochPts) = FailuresCondition.EpochData.Data(cnt, 1:EpochPts);
	for comp = 1:NumComp
		pcomp(1:EpochPts) = Components(1:EpochPts, comp);
		PCAProjections.FailuresProj(cnt, comp) = sum(resp .* pcomp);
	end
end

coef = {100 10 10};
coef = cat(1, coef{:});

for comp = 1:NumComp
	[Hist, Histx] = hist(PCAProjections.SinglesProj(:, comp), NumHistBins);
	PCAProjections.SinglesHist(comp, :) = Hist;
	PCAProjections.SinglesHistx(comp, :) = Histx;
	HistFit = coef;
	HistFit = nlinfit(Histx', Hist', 'gaussian', coef);
	PCAProjections.SinglesFit(comp, :) = gaussian(HistFit, Histx);
	PCAProjections.SinglesFitCoef(comp, :) = HistFit;
	[Hist, Histx] = hist(PCAProjections.FailuresProj(:, comp), NumHistBins);
	PCAProjections.FailuresHist(comp, :) = Hist;
	PCAProjections.FailuresHistx(comp, :) = Histx;
	HistFit = coef;
	HistFit = nlinfit(Histx', Hist', 'gaussian', coef);
	PCAProjections.FailuresFit(comp, :) = gaussian(HistFit, Histx);
	PCAProjections.FailuresFitCoef(comp, :) = HistFit;
end
