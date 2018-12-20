% -------------------------------------------------------------------------------------
% This script will extract single photon responses from cells
% -------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------
if exist('IsoSingles') 
	NumIsoSingles = length(IsoSingles);
else
	NumIsoSingles = 0;
end


RespDuration = 0.6;		% duration of each response vector
StartTimeA = 0.0;			% starting time for response A 
StartTimeB = 0.2;		% starting time for response B 
errflag = 0;			% error return
FreqCutoff = 10;		% maximum frequency to retain for analysos
% Set you default list of EpochCondition.Labels here.
ListNames{1} = 'fl_01';
ListNames{2} = 'fl_02';
ListNames{3} = 'fl_04';
%ListNames{4} = 'fl_08';
NumFlashConds = length(ListNames);
NumConds = length(CellInfo.EpochCondition);
clear CondPointer;

% Extract Conditions from CellInfo to be used in 2AFC starting with the dimmest flash strength and work to the brightest
for cond = 1:NumFlashConds
	CondPointer(NumFlashConds) = 0;
	for cnt = 1:NumConds
		if strcmp(ListNames(cond), CellInfo.EpochCondition(cnt).Label)
			CondPointer(cond) = CondPointer(cond) + 1;
			break
		else
			CondPointer(cond) = CondPointer(cond) +1 ;
		end
	end
end

% Pull out some useful info
BaselinePts = CellInfo.EpochCondition(CondPointer(1)).UserInfo.BaselinePts;
CurrentEpoch = CellInfo.EpochCondition(CondPointer(1)).EpochNumbers(1);
EpochPts = length(CellInfo.EpochData.Data{CurrentEpoch + 1});
StimulusDuration = (1/1000) * FindSearchPara(CellInfo.EpochCondition(CondPointer(1)), 'StimDur');
SamplingInterval = (1 / FindSearchPara(CellInfo.EpochCondition(CondPointer(1)), 'SampInterv')) * CellInfo.EpochCondition(CondPointer(1)).DecimatePts;

% convert times to point numbers
TrialLength = EpochPts - BaselinePts - (StimulusDuration / SamplingInterval);
RespLength = round(RespDuration / SamplingInterval);
StartA = round(StartTimeA / SamplingInterval);
StartB = round(StartTimeB / SamplingInterval);

% clear residual crap from condition structures
clear FlashConditionA;
clear FlashConditionB;
clear NoiseCondition;

% extract A and B vectors
for cond = 1:NumFlashConds
	fprintf(1, 'Condition #%d\n', cond);
	% get A and B vectors
	[ConditionA, ConditionB] = Extract2AFCResponses2(CellInfo, CellInfo.EpochCondition(CondPointer(cond)), StartA, StartB, RespLength, TrialLength, errflag, CellInfo.UserInfo.RandomizeFlag);
	% limit frequency content
	[CellInfo, ConditionA] = ZeroBaseLine(CellInfo, ConditionA, 100);
	[CellInfo, ConditionB] = ZeroBaseLine(CellInfo, ConditionB, 100);
	[CellInfo, ConditionA] = ApplyFrequencyCutoff(CellInfo, ConditionA, FreqCutoff);
	[CellInfo, ConditionB] = ApplyFrequencyCutoff(CellInfo, ConditionB, FreqCutoff);
	ConditionA = AverageAndVariance(CellInfo, ConditionA);
	ConditionB = AverageAndVariance(CellInfo, ConditionB);
	FlashConditionA(cond) = ConditionA;
	FlashConditionB(cond) = ConditionB;
end

% extract noise vectors by splitting up dark noise records
clear NoiseCondition
noisecond = 0;
for cond = 1:length(CellInfo.EpochCondition);
	if (strcmp(CellInfo.EpochCondition(cond).Label, 'dark_noise'))
		noisecond = cond;
	end
end
NoiseCondition = SplitNoiseIntoEpochs(CellInfo, CellInfo.EpochCondition(noisecond), RespLength + 100);
[CellInfo, NoiseCondition] = ApplyFrequencyCutoff(CellInfo, NoiseCondition, FreqCutoff);
NoiseCondition = ExcludeGlitches(CellInfo, NoiseCondition);
[CellInfo, NoiseCondition] = ZeroMean(CellInfo, NoiseCondition);
[NoiseCondition] = AverageAndVariance(CellInfo, NoiseCondition);


for flash = 1:NumFlashConds
	CheckStability(CellInfo, FlashConditionA(flash))
end	


%-------------------------------------------------------------------------------------------------------------------------
%	Get estimate of collecting area and the single photon response
%-------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------


% start and end points for computing ratio of variance and mean squared
StartPoint = 50;
EndPoint = 320;

% compute mean number of isomerizations for each flash strength from average and variance
NumFlashConds = length(CondPointer);
clear FlashStrength;
clear NumIsoms;
clear ave;
clear var;
NumIsoms(1:NumFlashConds) = 0;
FlashStrength(1:NumFlashConds) = 0;

% for each flash strength compute light-dependent variance.  Then compute ratio
% of square of average response to light-dependent variance -- this provides
% estimate of mean number of isomerizations at this flash strength.
for cond = 1:NumFlashConds
	ave = FlashConditionA(cond).AverageResponse - NoiseCondition.AverageResponse;
	var = FlashConditionA(cond).VarianceResponse - NoiseCondition.VarianceResponse;
	avesq = cumsum(ave.^2);
	var = cumsum(var);
	NumIsoms(cond) = (avesq(EndPoint) - avesq(StartPoint)) / (var(EndPoint) - var(StartPoint));
	FlashStrength(cond) = FlashConditionA(cond).UserInfo.StimulusAmp;
end

% fit plot of photoisomerizations versus flash strength to obtain collecting area
FlashStrength = FlashStrength';
slope = {1};
slope = cat(1, slope{:});
slopefit = slope;
slopefit = nlinfit(FlashStrength, NumIsoms, 'ZeroInterceptLine', slope);
fit = ZeroInterceptLine(slopefit, FlashStrength);
plot(FlashStrength, NumIsoms, FlashStrength, fit);
CellParameters.CollectingArea = slopefit;
fprintf(1, 'Estimated collecting area = %d square microns\n', CellParameters.CollectingArea);

% sanity check on scaling factors -- compare variance to scaled square of average
for cond = 1:NumFlashConds
	ave = FlashConditionA(cond).AverageResponse - NoiseCondition.AverageResponse;
	var = FlashConditionA(cond).VarianceResponse - NoiseCondition.VarianceResponse;
	ave = ave.^2 / NumIsoms(cond);
	tme = 1:length(ave);
	plot(tme, var, tme, ave);
	GlitchFlag = input('hit return to continue');
end

% estimate single photon response
for cond = 1:NumFlashConds
	if (cond == 1)
		CellParameters.ElementaryResponse = FlashConditionA(cond).AverageResponse / (FlashStrength(cond) * CellParameters.CollectingArea);
	else
		CellParameters.ElementaryResponse = CellParameters.ElementaryResponse + FlashConditionA(cond).AverageResponse / (FlashStrength(cond) * CellParameters.CollectingArea);
	end
end

CellParameters.ElementaryResponse = CellParameters.ElementaryResponse / NumFlashConds;
plot(CellParameters.ElementaryResponse);



% -------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------
clear FlashConditionB
FlashConditionB = NoiseCondition
CollectingArea = CellInfo.UserInfo.CellParameters.CollectingArea;
global DarkNoiseSD;
global NumResponses;
global NumberCombConditions;
global NumHistBins;
global FlashStrength;
global MeanPhotoisoms;
NumFlashConds = length(CondPointer);
NumberCombConditions = NumFlashConds;
WindowSize = 50;
RespSize = 350;

% calculate statistics of response amplitude, area and projections
clear RespProjection;
clear ResponseStats;
clear NoiseStats;
for cond=1:NumFlashConds
	[Stats, NStats] = CalculateResponseStats(FlashConditionA(cond), NoiseCondition, CellParameters, RespSize, WindowSize);
	ResponseStats(cond) = Stats;
	NoiseStats(cond) = NStats;
end

% calculate histograms of response statistics
clear Respx;
clear Areax;
NumHistBins = 150;
Respx(1:NumHistBins) = -NumHistBins/2:NumHistBins/2-1;
Respx = (floor(max(ResponseStats(NumFlashConds).MaxAmplitude)) + 1) * Respx * 2 / NumHistBins;
Respx = Respx';
Areax(1:NumHistBins) = -NumHistBins/2:NumHistBins/2-1;
Areax = (floor(max(ResponseStats(NumFlashConds).Area)) + 1) * Areax * 2 / NumHistBins;
RegressAmpx(1:NumHistBins) = -NumHistBins/2:NumHistBins/2-1;
RegressAmpx = (floor(max(ResponseStats(NumFlashConds).RegressAmp)) + 1) * RegressAmpx * 2 / NumHistBins;
for cond = 1:NumFlashConds
	% flash response amplitude
	temp = ResponseStats(cond).MaxAmplitude;
	[HistA, Respx] = hist(temp', Respx);
	ResponseStats(cond).AmpHist = HistA;
	ResponseStats(cond).AmpHistx = Respx;
	% dark noise amplitude
	temp = NoiseStats(cond).MaxAmplitude;
	[HistA, Respx] = hist(temp', Respx);
	NoiseStats(cond).AmpHist = HistA;
	NoiseStats(cond).AmpHistx = Respx;
	% flash response area
	temp = ResponseStats(cond).Area;
	[HistA, Areax] = hist(temp', Areax);
	ResponseStats(cond).AreaHist = HistA;
	ResponseStats(cond).AreaHistx = Areax;
	% dark noise area
	temp = NoiseStats(cond).Area;
	[HistA, Areax] = hist(temp', Areax);
	NoiseStats(cond).AreaHist = HistA;
	NoiseStats(cond).AreaHistx = Areax;
	% flash response projection
	temp = ResponseStats(cond).RegressAmp;
	[HistA, RegressAmpx] = hist(temp', RegressAmpx);
	ResponseStats(cond).RegressAmpHist = HistA;
	ResponseStats(cond).RegressAmpHistx = RegressAmpx;
	% dark noise projection
	temp = NoiseStats(cond).RegressAmp;
	[HistA, RegressAmpx] = hist(temp', RegressAmpx);
	NoiseStats(cond).RegressAmpHist = HistA;
	NoiseStats(cond).RegressAmpHistx = RegressAmpx;
end

% examine histograms
for cond = 1:NumFlashConds
	plot(ResponseStats(cond).AmpHistx, ResponseStats(cond).AmpHist, NoiseStats(cond).AmpHistx, NoiseStats(cond).AmpHist)
	GlitchFlag = input('hit return to continue');
	plot(ResponseStats(cond).AreaHistx, ResponseStats(cond).AreaHist, NoiseStats(cond).AreaHistx, NoiseStats(cond).AreaHist)
	GlitchFlag = input('hit return to continue');
	plot(ResponseStats(cond).RegressAmpHistx, ResponseStats(cond).RegressAmpHist, NoiseStats(cond).RegressAmpHistx, NoiseStats(cond).RegressAmpHist)
	GlitchFlag = input('hit return to continue');
end

% for response projections
clear Respx;
Respx = RegressAmpx;
for cond = 1:NumFlashConds
	RespProjection(cond).RespAProjection = ResponseStats(cond).RegressAmp;
	RespProjection(cond).NoiseProjection = NoiseStats(cond).RegressAmp;
	RespProjection(cond).RespAHist = ResponseStats(cond).RegressAmpHist;
	RespProjection(cond).NoiseHist = NoiseStats(cond).RegressAmpHist;
	FlashStrength(cond) = FlashConditionA(cond).UserInfo.StimulusAmp;
	MeanPhotoisoms(cond) = FlashConditionA(cond).UserInfo.StimulusAmp * CellParameters.CollectingArea;
end

% calculate projections of responses through discriminant
clear RespProjection;
shuffles = 10;
RespLength = 700;
RespProjection = CalcRespProjsandHistsSing(CellInfo, FlashConditionA, NoiseCondition, FreqCutoff, RespLength, shuffles);
NumHistBins = RespProjection(1).NumHistBins;

% fit yoked comb distribution to histograms -- hold collecting area fixed
[CombParameters, CombFit] = FitYokedCombDistSingles(RespProjection);
CellParameters.CombParameters = CombParameters;
CombParameters
SinglesCondition.UserInfo.CellParameters = CellParameters;


% -------------------------------------------------------------------------------------
%  Isolate the Singles
% -------------------------------------------------------------------------------------
% -------------------------------------------------------------------------------------
SingleAmp = CombParameters.UnitaryAmpA

% ISOLATE SINGLES
LowerFailLimit = -1;			% smallest failure
UpperFailLimit = SingleAmp - (0.5 * SingleAmp);				% largest failure
LowerSingleLimit = SingleAmp - (0.5 * SingleAmp);			% smallest single
UpperSingleLimit = SingleAmp + (SingleAmp * 0.5);			% largest single

% ISOLATE SINGLES
LowerFailLimit = -1;			% smallest failure
UpperFailLimit = 0.4;				% largest failure
LowerSingleLimit = 0.4;			% smallest single
UpperSingleLimit = 1.7;			% largest single

% identify failures and singles based on response projections
clear FailuresCondition;
clear SinglesCondition;
[SinglesCondition, SinglesExp] = FindResponses(ResponseStats, FlashConditionA, CellParameters, NumFlashConds, LowerSingleLimit, UpperSingleLimit);
[FailuresCondition, FailsExp] = FindResponses(ResponseStats, FlashConditionA, CellParameters, NumFlashConds, LowerFailLimit, UpperFailLimit);
EpochPts = length(SinglesCondition.EpochData.Data(1,:));
tme=1:EpochPts;
tme = tme /1000;
plot(0);
SinglesCondition = AverageAndVariance(CellInfo, SinglesCondition);
DisplayEpochs(SinglesCondition);
plot(0);
FailuresCondition = AverageAndVariance(CellInfo, FailuresCondition);
DisplayEpochs(FailuresCondition);
plot(tme, CellParameters.ElementaryResponse, tme, SinglesCondition.AverageResponse, tme, FailuresCondition.AverageResponse);
CellParameters.SinglesExp = SinglesExp;
CellParameters.FailsExp = FailsExp;
SinglesCondition.UserInfo.CellParameters = CellParameters;


% simple plotting
CheckStability(CellInfo, SinglesCondition);
CheckStability(CellInfo, FailuresCondition);
[SinglesCondition, BadEpochs] = ScanEpochs(SinglesCondition, BadEpochs);
[FailuresCondition, BadEpochs] = ScanEpochs(FailuresCondition, BadEpochs);

% time-dependent variance
var = SinglesCondition.VarianceResponse - FailuresCondition.VarianceResponse;
avesq = SinglesCondition.AverageResponse.^2;
avesq = avesq / 10;
plot(tme, SinglesCondition.VarianceResponse, tme, avesq, tme, FailuresCondition.VarianceResponse, tme, var);
plot(tme, avesq, tme, var);

% sanity check: sort singles by flash strength they originate from and compare averages
for cond = 1:NumFlashConds
	clear Condition
	clear Stats
	Condition = FlashConditionA(cond);
	Stats = ResponseStats(cond);
	TestCondition(1) = Condition;
	RStats(1) = Stats;
	clear TestSinglesCondition;
	TestSinglesCondition = FindResponses(RStats, TestCondition, CellParameters, 1, LowerSingleLimit, UpperSingleLimit);
	TestSinglesCondition = AverageAndVariance(CellInfo, TestSinglesCondition);
	AveSingle(cond, :) = TestSinglesCondition.AverageResponse;
	VarSingle(cond, :) = TestSinglesCondition.VarianceResponse;
end

plot(0);
hold on;
for cond = 1:NumFlashConds
	PlotColor = input('What color do you want? \n','s');
	plot(tme, AveSingle(cond, :),PlotColor);
end
hold off;

plot(0);
hold on;
for cond = 1:NumFlashConds
	PlotColor = input('What color do you want? \n','s');
	plot(tme, VarSingle(cond, :),PlotColor);
end
hold off;


%----------------------------------------------------------------------------------
% CHECK INDEPENDENCE OF DARK NOISE AND SINGLES FLUCTUATIONS
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------

% divide singles into those starting with large and small dark current and compare means
MeanDarkCurrent = mean(SinglesCondition.EpochData.Offset);
LargeDarkCnt = 0;
SmallDarkCnt = 0;
clear LargeDarkAve;
clear SmallDarkAve;
LargeAveOffset = 0;
SmallAveOffset = 0;
for cnt=1:length(SinglesCondition.EpochNumbers)
	if (SinglesCondition.EpochData.Offset(cnt) > MeanDarkCurrent)
		if (LargeDarkCnt == 0)
			LargeDarkAve = SinglesCondition.EpochData.Data(cnt, :);
		else
			LargeDarkAve = LargeDarkAve + SinglesCondition.EpochData.Data(cnt, :);
		end
		LargeAveOffset = LargeAveOffset + SinglesCondition.EpochData.Offset(cnt);
		LargeDarkCnt = LargeDarkCnt + 1;
	else
		if (SmallDarkCnt == 0)
			SmallDarkAve = SinglesCondition.EpochData.Data(cnt, :);
		else
			SmallDarkAve = SmallDarkAve + SinglesCondition.EpochData.Data(cnt, :);
		end
		SmallDarkCnt = SmallDarkCnt + 1;
		SmallAveOffset = SmallAveOffset + SinglesCondition.EpochData.Offset(cnt);
	end
end
LargeDarkAve = LargeDarkAve / LargeDarkCnt + LargeAveOffset / LargeDarkCnt;
SmallDarkAve = SmallDarkAve / SmallDarkCnt + SmallAveOffset / SmallDarkCnt;
tme = 1:length(LargeDarkAve);
plot(tme, LargeDarkAve, tme, SmallDarkAve);

MeanDarkCurrent = mean(FailuresCondition.EpochData.Offset);
LargeDarkCnt = 0;
SmallDarkCnt = 0;
LargeAveOffset = 0;
SmallAveOffset = 0;
clear FailLargeDarkAve;
clear FailSmallDarkAve;
clear tme;
for cnt=1:length(FailuresCondition.EpochNumbers)
	if (FailuresCondition.EpochData.Offset(cnt) > MeanDarkCurrent)
		if (LargeDarkCnt == 0)
			FailLargeDarkAve = FailuresCondition.EpochData.Data(cnt, :);
		else
			FailLargeDarkAve = FailLargeDarkAve + FailuresCondition.EpochData.Data(cnt, :);
		end
		LargeDarkCnt = LargeDarkCnt + 1;
		LargeAveOffset = LargeAveOffset + FailuresCondition.EpochData.Offset(cnt);
	else
		if (SmallDarkCnt == 0)
			FailSmallDarkAve = FailuresCondition.EpochData.Data(cnt, :);
		else
			FailSmallDarkAve = FailSmallDarkAve + FailuresCondition.EpochData.Data(cnt, :);
		end
		SmallDarkCnt = SmallDarkCnt + 1;
		SmallAveOffset = SmallAveOffset + FailuresCondition.EpochData.Offset(cnt);
	end
end
FailLargeDarkAve = FailLargeDarkAve / LargeDarkCnt + LargeAveOffset / LargeDarkCnt;
FailSmallDarkAve = FailSmallDarkAve / SmallDarkCnt + SmallAveOffset / SmallDarkCnt;
tme = 1:length(FailLargeDarkAve);
plot(tme, FailLargeDarkAve, tme, FailSmallDarkAve, tme, LargeDarkAve, tme, SmallDarkAve);

LargeDarkAve = LargeDarkAve - FailLargeDarkAve;
SmallDarkAve = SmallDarkAve - FailSmallDarkAve;
plot(tme, LargeDarkAve, tme, SmallDarkAve);
IsoSingles(NumIsoSingles + 1).LargeDarkAve = LargeDarkAve;
IsoSingles(NumIsoSingles + 1).SmallDarkAve = SmallDarkAve;


%----------------------------------------------------------------------------------
% STATISTICAL ANALYSIS OF SINGLES
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------

% determine principle components of singles
SinglesCondition = DecimateCondition(SinglesCondition, 10);
FailuresCondition = DecimateCondition(FailuresCondition, 10);
SinglesCondition = NormalizeSingles(SinglesCondition);
FailuresCondition = NormalizeFailures(FailuresCondition, -1.1895);
SinglesCondition = AverageAndVariance(CellInfo, SinglesCondition);


SinglesCovarMatrix = PrincipleComponents(CellInfo, SinglesCondition);
FailuresCovarMatrix = PrincipleComponents(CellInfo, FailuresCondition);
mesh(SinglesCovarMatrix);

%clear SinglesEigVal;
%clear FailuresEigVal;
%clear EigVec;
%SinglesEigVal = eig(SinglesCovarMatrix);
%FailuresEigVal = eig(FailuresCovarMatrix);
%cnt = 1:length(SinglesEigVal);
%semilogy(cnt, SinglesEigVal, '*', cnt, FailuresEigVal, 'o');
%axis([50 70 10 1e6])

[pcSingles, SingleEigVals, SingleVarCap] = pcacov(SinglesCovarMatrix);
[pcFails, FailsEigVals, FailsVarCap] =pcacov(FailuresCovarMatrix);
cnt = 1:length(SingleEigVals);
semilogy(cnt, SingleEigVals, '*', cnt, FailsEigVals, 'o')
axis([0 20 1 1e6])
%figure
%semilogy(cnt, SingleVarCap, '*', cnt, FailsVarCap, 'o')
%axis([0 20 0.1 100])

plot(SinglesCondition.UserInfo.TimePoints, pcFails(:,1),SinglesCondition.UserInfo.TimePoints,pcFails(:,2))
plot(SinglesCondition.UserInfo.TimePoints, pcSingles(:,1),SinglesCondition.UserInfo.TimePoints,pcSingles(:,2))

%[EigVec, EigVal] = eigs(FailuresCovarMatrix, 2);
%plot(SinglesCondition.UserInfo.TimePoints, EigVec);
%[EigVec, EigVal] = eigs(SinglesCovarMatrix, 2);
%plot(SinglesCondition.UserInfo.TimePoints, EigVec);
%CellParameters.EigVec = EigVec;

EigVec(:,1) = pcSingles(:,1);
EigVec(:,2) = pcSingles(:,2);
EigVal(1) = SingleEigVals(1);
EigVal(2) = SingleEigVals(2);
FailEigVec(:,1) = pcFails(:,1);
FailEigVec(:,2) = pcFails(:,2);
FailEigVal(1) = FailsEigVals(1);
FailEigVal(2) = FailsEigVals(2);

% look at distribution of projections along principle components
clear PCAProjections;
PCAProjections = PrincipleCompProjections(SinglesCondition, FailuresCondition, EigVec);
for comp = 1:size(EigVec, 2)
	plot(0);
	hold on
	plot(PCAProjections.SinglesHistx(comp, :), PCAProjections.SinglesHist(comp, :), PCAProjections.SinglesHistx(comp, :), PCAProjections.SinglesFit(comp, :));
	plot(PCAProjections.FailuresHistx(comp, :), PCAProjections.FailuresHist(comp, :), PCAProjections.FailuresHistx(comp, :), PCAProjections.FailuresFit(comp, :));
	hold off
	GlitchFlag = input('hit return to continue');
end	
CellParameters.PCAProjections = PCAProjections;
stairs(PCAProjections.SinglesHistx(1, :), PCAProjections.SinglesHist(1, :), 'b')
hold
plot(PCAProjections.SinglesHistx(1, :), PCAProjections.SinglesFit(1, :), 'b')
stairs(PCAProjections.SinglesHistx(2, :), PCAProjections.SinglesHist(2, :), 'g')
plot(PCAProjections.SinglesHistx(2, :), PCAProjections.SinglesFit(2, :), 'g')
axis([-10 10 0 20])

% examine particular response and PCA fit
resp = 70;
Fit = EigVec(:, 1) * PCAProjections.SinglesProj(resp, 1) + EigVec(:, 2) * PCAProjections.SinglesProj(resp, 2);
plot(SinglesCondition.UserInfo.TimePoints, SinglesCondition.EpochData.Data(resp, :), SinglesCondition.UserInfo.TimePoints, Fit);

% look for correlations in principle component projections
plot(PCAProjections.SinglesProj(:, 1), PCAProjections.SinglesProj(:, 2), '*');

% check for time-dependence of principle component projections (sanity check for drift in experiment)
for comp = 1:size(EigVec, 2)
	plot(SinglesCondition.EpochTime, PCAProjections.SinglesProj(:, comp), '*');
	GlitchFlag = input('hit return to continue');
end

% look for set of independent response components
Component1 = PCAProjections.SinglesProj(:, 1);
Component2 = PCAProjections.SinglesProj(:, 2);
coef = {1 1};
coef = cat(1, coef{:});
linefit = coef;
linefit = nlinfit(Component1, Component2, 'StraightLineFit', coef);
fit = StraightLineFit(linefit, Component1);
plot(Component1, Component2, '*', Component1, fit);

% orthogonalize principle components
EigVec = FindIndependentComponents(EigVec, linefit(1));
PCAProjections = PrincipleCompProjections(SinglesCondition, FailuresCondition, EigVec);
PCAProjections.EigVec = EigVec;

plot(EigVec);

plot(SinglesCondition.UserInfo.TimePoints, EigVec(:,1), SinglesCondition.UserInfo.TimePoints, EigVec(:,2))

% recheck distribution of projections
for comp = 1:size(EigVec, 2)
	plot(0);
	hold on
	plot(PCAProjections.SinglesHistx(comp, :), PCAProjections.SinglesHist(comp, :), PCAProjections.SinglesHistx(comp, :), PCAProjections.SinglesFit(comp, :));
	plot(PCAProjections.FailuresHistx(comp, :), PCAProjections.FailuresHist(comp, :), PCAProjections.FailuresHistx(comp, :), PCAProjections.FailuresFit(comp, :));
	hold off
	GlitchFlag = input('hit return to continue');
end	
stairs(PCAProjections.SinglesHistx(1, :), PCAProjections.SinglesHist(1, :), 'b')
hold
plot(PCAProjections.SinglesHistx(1, :), PCAProjections.SinglesFit(1, :), 'b')
stairs(PCAProjections.SinglesHistx(2, :), PCAProjections.SinglesHist(2, :), 'g')
plot(PCAProjections.SinglesHistx(2, :), PCAProjections.SinglesFit(2, :), 'g')


SinglesCondition.UserInfo.PCAProjections = PCAProjections;
SinglesCondition.UserInfo.PCVectors = EigVec;
SinglesCondition.UserInfo.PCVals = EigVal;
FailuresCondition.UserInfo.PCVectors = FailEigVec;
FailuresCondition.UserInfo.PCVals = FailEigVal;


% Store information about singles into this structure
IsoSingles(NumIsoSingles + 1).FileName = CellInfo.CellFile;
IsoSingles(NumIsoSingles + 1).SinglesCondition = SinglesCondition;
IsoSingles(NumIsoSingles + 1).FailuresCondition = FailuresCondition;

IsoSingles(NumIsoSingles + 1).IsolationParameters.LowerFailLimit = LowerFailLimit;			
IsoSingles(NumIsoSingles + 1).IsolationParameters.UpperFailLimit = UpperFailLimit;				
IsoSingles(NumIsoSingles + 1).IsolationParameters.LowerSingleLimit = LowerSingleLimit;			
IsoSingles(NumIsoSingles + 1).IsolationParameters.UpperSingleLimit = UpperSingleLimit;			


%----------------------------------------------------------------------------------
% SIMULATED SINGLE PHOTON RESPONSES BASED ON PCA
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------

% primate rod continuous noise: time constant 430 msec, SD 0.43 pA
% generate simulated responses
tic
clear SimulateCondition
NumResponses = 100;
CellParameters.ContinuousNoisePow = 4;
CellParameters.ContinuousNoiseTimeConstant = 0.46;
CellParameters.ContinuousNoiseSD = 0.0;
CellParameters.ThermalRate = 0.006;
Comp2SDScaleFact = 0.01;
CellParameters.PCAProjections.SinglesFitCoef(2,3) = CellParameters.PCAProjections.SinglesFitCoef(2,3) * Comp2SDScaleFact;
SimulateCondition = SimulateRodResponse(SinglesCondition, CellParameters, PCAProjections, NumResponses, 1, 1, 0, 1.0);
SimulateCondition = AverageAndVariance(CellInfo, SimulateCondition);
avesqsim = SimulateCondition.AverageResponse.^2;
avesqsim = avesqsim / 10;
avesq = SinglesCondition.AverageResponse.^2;
avesq = avesq / 10;
var = SinglesCondition.VarianceResponse - FailuresCondition.VarianceResponse;
tme = (1:length(SimulateCondition.AverageResponse)) * 0.01;
plot(tme, SimulateCondition.VarianceResponse, tme, avesqsim, tme, var, tme, avesq);
plot(2 * tme, SimulateCondition.VarianceResponse, 2 * tme, avesqsim);
toc
pause
plot(SinglesCondition.UserInfo.TimePoints, SimulateCondition.VarianceResponse, SinglesCondition.UserInfo.TimePoints, SinglesCondition.VarianceResponse)


CellParameters.ContinuousNoisePow = 4;
CellParameters.ContinuousNoiseTimeConstant = 0.43;
CellParameters.ContinuousNoiseSD = 0.0;
CellParameters.ThermalRate = 0.01;
NumIsoms = 1;
SimulateCondition = SimulateRodResponse(SinglesCondition, CellParameters, PCAProjections, NumResponses, 1, NumIsoms, 0, 1.0);
SimulateCondition = AverageAndVariance(SimulateCondition);
avesqsim = SimulateCondition.AverageResponse.^2 / 10;
plot(tme, SimulateCondition.VarianceResponse, tme, avesqsim);
NoiseCondition = SimulateRodResponse(SinglesCondition, CellParameters, PCAProjections, NumResponses, 4, 0, 0, 1.0);
NoiseCondition = AverageAndVariance(NoiseCondition);
plot(tme, SimulateCondition.VarianceResponse, tme, NoiseCondition.VarianceResponse, tme, avesqsim, tme, SimulateCondition.VarianceResponse - NoiseCondition.VarianceResponse);

% do two-alternative forced choice discrimination between response group A and B
shuffles = 5;
FreqCutoff = 10;
Verbose = 1;
PCorrect = DoDiscrimination(SimulateCondition, NoiseCondition, 1, shuffles, FreqCutoff, Verbose);


%----------------------------------------------------------------------------------
% Simulated pool of Rods surface with Principle Components
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------

NumTrials = 2000;
clear SimulateCondition
NumResponses = 200;
CellParameters.ContinuousNoisePow = 4;
CellParameters.ContinuousNoiseTimeConstant = 0.46;
CellParameters.ContinuousNoiseSD = 0.00;
CellParameters.ThermalRate = 0.006;
Comp2SDScaleFact = 1;
CellParameters.PCAProjections.SinglesFitCoef(2,3) = CellParameters.PCAProjections.SinglesFitCoef(2,3) * Comp2SDScaleFact;

ResponseLength = 1400
FlashStrengths = [0.25 0.5 0.75 1 1.5 2 3 4]
StartFlashStrength = 0.2;
NumFlashConds = length(FlashStrengths);
StartShift = 0.01;
NumShifts = 6;
StartTimeA = 0;
Verbose = 0;

errflag = 0;			% error return
FreqCutoff = 10;
%shuffles = 10;		% how many bootstrap resamplings
Verbose = 0;		% print lots of crap to command window?
SamplingInterval = 0.01

clear SimulateConditionA
for flash = 1:NumFlashConds
	clear Condition
	Condition = SimulateRodResponse(SinglesCondition, CellParameters, PCAProjections, NumResponses, 1, FlashStrengths(flash), 1, 1.0);
	SimulateConditionA(flash) = AverageAndVariance(CellInfo, Condition);
end

% Resample the data so that 1 pt is 1 ms apart.  
InterpFactor = 10;
for flash = 1:NumFlashConds
	clear Condition
	Condition = InterpolateCondition(SimulateConditionA(flash), InterpFactor);
	NewSimCondition(flash) = Condition;
end
clear SimulateConditionA;
SimulateConditionA = NewSimCondition;
clear NewSimCondition;

for flash = 1:NumFlashConds
	clear Condition
	Condition = DecimateCondition(SimulateConditionA(flash), 10);
	NewSimCondition(flash) = Condition;
end
clear SimulateConditionA
SimulateConditionA = NewSimCondition;
clear NewSimCondition

clear RotatedSimConditionB
for shift = 1:NumShifts
	clear PCorrect 
	% rotate B vectors
	ShiftB = floor((StartTimeA + 2^(shift-1) * StartShift) / SamplingInterval);
	for flash = 1:NumFlashConds
		if (Verbose)
			fprintf(1, 'Rotating condition %d\n', cond);
		end
		clear Condition
		Condition = RotateResponseVectors(SimulateConditionA(flash), ShiftB);
		RotatedSimConditionB(flash) = Condition;
	end

	fprintf(1, 'Starting shift %d\n', ShiftB * SamplingInterval);
		
	%----------------------------------------------------------------------------------
	% DISCRIMINATION VERSUS TEMPORAL SHIFT: 2D DISCRIMINATION SURFACES
	%----------------------------------------------------------------------------------

	% do two-alternative forced choice discrimination between response group A and B
	PCorrect(1:NumFlashConds) = 0;
	Verbose = 0;
	for cond = 1:NumFlashConds
		PCorrect(cond) = Do2AFCDiscrimination(CellInfo, SimulateConditionA, RotatedSimConditionB, cond, NumTrials, FreqCutoff, 1400, Verbose);
	end
	
	% store curve
	PCorrect_2D(shift, 1:length(PCorrect)) = PCorrect;
	TimeShift(shift) = floor((2^(shift-1) * StartShift) / SamplingInterval) * SamplingInterval;
	
	Verbose = 0;
	if (Verbose)
		plot(NumIsoms, PCorrect);
		pause(1)
		%GlitchFlag = input('hit return to continue');
	end
end	





