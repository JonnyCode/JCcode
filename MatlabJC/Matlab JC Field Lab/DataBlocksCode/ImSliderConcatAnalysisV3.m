function ForIgor = ImSliderConcatAnalysisV3(DataBlock, DB, Params)


% this function is adapted from 'ImSliderConcatAnalysisV2' and uses an OLE
% instead of vector sum

% JC 8/24/2017 

LoadOldMatFileFlag = false; % true skip calculations and just load
PlotFigsFlag = true ;
saveFigPath = ['/Users/jcafaro/Documents/AnalysisFigures/'] ;
numElectrodeLayers = 2 ; % number of electrode layers surrounding the max that the center is calculated over
eiDistBin = 100 ; % (um) bins 
eiDistMax = 1000 ; % (um) max distance ei can be separated
StimOnOffset = 1 ; % (sec) time after stim begins to start counting spikes
DgOffTime = 8 ; % (sec) time ds typing drifiting grating turned off
minEiRad = 100 ; % (um) minimum distance to be within a RGC set
minDenominator = 10^-5 ; % added to denominator to avoid Nans
psthBinTime = .001 ; % size of psth Bin step
psthSmoothTime = 1 ; % size of psth Bin (smooth sliding window)

ConcatPathNum = 3 ; %TEMP -SHOULD BE INPUT
ImPathNum = 3 ;
DsPathNum = 2 ;

% parameters
Color_list = {'k','r','b','g','c','y','k','r','b','g','c','y','k','r','b','g','c','y','k','r','b','g','c','y'} ; % order of colors for each 

saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/DsIdsDb',num2str(DB),'ConcatPathNum',num2str(ConcatPathNum)] ;
    
% load stimulus 
slashi = strfind(DataBlock(DB).ImSlide{ImPathNum},'/') ; % find the /
StimPath = [DataBlock(DB).ImSlide{ImPathNum}(1:slashi(end-1)),'stimuli/s',DataBlock(DB).ImSlide{ImPathNum}(end-1:end)] ;
load(StimPath) ;

% correct directions > 360
stimulus.directions(stimulus.directions>=360)=stimulus.directions(stimulus.directions>=360)-360 ;
stimulus.directions = sort(stimulus.directions) ;
stimulus.directionsShown(stimulus.directionsShown>=360)=stimulus.directionsShown(stimulus.directionsShown>=360)-360 ;

% load data
dataRun = load_data(DataBlock(DB).ImSlideConcat{ConcatPathNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun.piece.array_id = 1551 ; % CAUTION
dataRun = load_ei(dataRun, 'all') ;

% ImSlide stimulus time
dataRunIm = load_data(DataBlock(DB).ImSlide{ImPathNum}) ;
dataRunIm = load_neurons(dataRunIm) ;
ImStimTime = dataRunIm.duration ;
DgStimTime = dataRun.duration-ImStimTime ;
clear dataRunIm ;

% get ei center of mass 
for cells = 1:length(dataRun.spikes) ; % for each cell
    EiCtr(cells,:) = get_ei_com(dataRun, dataRun.cell_ids(cells), numElectrodeLayers) ;
end

% get triggers from 
triggs = dataRun.triggers(dataRun.triggers>DgStimTime) ; % image slide

% stim params
if iscell(stimulus.image) ;
    NumImages = length(stimulus.image) ; % number of images shown
else
    NumImages = 1 ;
end
NumDir = size(stimulus.directionsShown,2) ; % number of directions
NumTrialsAttempted = size(stimulus.directionsShown,1) ;
NumStimShort = NumImages * NumDir * NumTrialsAttempted -length(triggs) ; % number of stimuli short

% look for single missing triggers and estimate
TriggOutlierFactor = 1.75 ;
if NumStimShort>0 ; 
    triggsFixed = triggs ;
    triggDiffTime = median(diff(triggs)) ;
    for t=2:length(triggs) ;
        if (triggs(t)-triggs(t-1)) > triggDiffTime*TriggOutlierFactor ; % if the trigger is much later than other triggers
            NewTrigg = (triggs(t)-triggs(t-1))/2 
            triggsFixed = [triggsFixed(1:t-1); NewTrigg; triggs(t:end)] ; % add new trigger
        end
    end
    triggs = triggsFixed ;
end


% take out trial sets that were not complete
NumCompleteTrials = floor(length(triggs)/(NumDir*NumImages)) ; % number of trials with all directions and all images
directionsShown = stimulus.directionsShown(1:NumCompleteTrials,:) ; % take out the entire unfinished trials
triggs = triggs(1:NumCompleteTrials*NumImages*NumDir) ; % triggers for stim included in completed trials
 
% expand directionsShown for mutliple images
directionsShownFull = [] ;
for t = 1:NumCompleteTrials ; % for each trial
    directionsShownFull = [directionsShownFull;repmat(directionsShown(t,:),[NumImages 1])];
end

% select only one of the images (comment out if want all images)
% selectImage = 1 ;
% temp = nans(size(directionsShownFull)) ;
% temp([selectImage:NumImages:size(directionsShownFull,1)],:) = ...
%     directionsShownFull([selectImage:NumImages:size(directionsShownFull,1)],:) ;
% directionsShownFull = temp ;

% stim params
StimDuration = mean(diff(triggs)) ; % average stimulus duration
StimFrameRate = StimDuration/stimulus.num_frames ; % average frame rate sec/frame
dirShownFullTranspose = directionsShownFull' ;

% id ds cells in DsPath
Params.TimeBounds = [0,DgStimTime] ;
Params.DsPathNum = DsPathNum ;
try load(saveDsIdsPath) ; % if they are ds ids already saved
catch % if not find them

    DataBlock(DB).DsConcatPath = DataBlock(DB).ImSlideConcat{ConcatPathNum} ;
    ForIgor = DsCellFinder(DataBlock, DB, Params) ;

    save(saveDsIdsPath, 'ForIgor')
end

% On-Off cell from DsPath
cell_id = ForIgor.ds_id{2} ; 
DsTypeName = ForIgor.dsName{2} ; 
dsi = ForIgor.dsi{2} ;

for DsType=1:length(DsTypeName) ;
    DsCelli{DsType} = get_cell_indices(dataRun,cell_id{DsType}) ;
    DsCellTypeNum(DsType) = length(DsCelli{DsType}) ; 
end

% spike count 
for cells=1:length(dataRun.spikes) ; % for each cell
    for t=1:length(triggs) ; % for each trigger
        spk = dataRun.spikes{cells}-triggs(t) ;
        spk = spk(spk>=0+StimOnOffset & spk<=StimDuration) ;
        spkCount(cells,t) = length(spk) ; %#ok<AGROW>
    end
end

% spike count for each image direction
for d = 1:NumDir ; % for each diection
   di = find(dirShownFullTranspose==stimulus.directions(d)) ; % trials with that direction
   spkCount_direction_mean(:,d) = mean(spkCount(:,di),2) ;
   spkCount_direction_std(:,d) = std(spkCount(:,di),[],2) ;
end
 
% taylor/vaney dsi
for cells=1:length(dataRun.spikes) ; % for each cell
    PV = PolarVectorAddition([stimulus.directions',spkCount_direction_mean(cells,:)']) ;
    PreferedDirection(cells) = PV(1) ;
    DSindex(cells) = PV(2)/sum(spkCount_direction_mean(cells,:)) ;
end
    
% psth
%psthTime = [0:StimFrameRate:StimDuration] ; 
psthTime = [0:psthBinTime:StimDuration] ;
psthSmoothPnts = floor(psthSmoothTime/psthBinTime) ; % number of psthBins per smooth window 
for cells=1:length(dataRun.spikes) ; % for each cell
    for st = 1:NumDir ; % for each stimulus direction
        ti = find(dirShownFullTranspose(:)==stimulus.directions(st)) ; % index of triggers for that stim
        for t=1:NumCompleteTrials ; % for each trial
            spk = dataRun.spikes{cells}-triggs(ti(t)) ;
            spk = spk(spk>=0 & spk<=StimDuration) ;
            %psth{cells}{st}(t,:) = histc(spk,psthTime) ;
            
            spikeTrain = zeros(1,length(psthTime)) ; % empty spike train
            spkPnts = 1+floor(spk/psthBinTime) ; % spike points 
            spikeTrain(spkPnts) = 1 ; % spike train
            psth{cells}{st}(t,:) = smooth(spikeTrain,psthSmoothPnts) ; % psth
        end
        psth_mean{cells}(st,:)=mean(psth{cells}{st},1) ;
        psth_var{cells}(st,:)=var(psth{cells}{st},[],1) ;
    end
end

%% drifting grating data

% stim path
TrialTrigInterval = 10 ;
slashi = strfind(DataBlock(DB).DsPath{DsPathNum},'/') ; % find the /
dataRunDg.triggers = dataRun.triggers(dataRun.triggers<=DgStimTime) ;
dataRunDg.names.stimulus_path = [DataBlock(DB).DsPath{Params.DsPathNum}(1:slashi(end-1)),'stimuli/s',DataBlock(DB).DsPath{DsPathNum}(end-1:end),'.txt'] ;
dataRunDg = load_stim(dataRunDg,'user_defined_trigger_interval', TrialTrigInterval) ;

% psth
%psthTimeDg = [0:StimFrameRate:TrialTrigInterval] ; 
psthTimeDg = [0:psthBinTime:TrialTrigInterval] ;
for cells=1:length(dataRun.spikes) ; % for each cell
    for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus 
        ti = find(dataRunDg.stimulus.trial_list==st) ; % index of triggers for that stim
        for t=1:length(ti) ; % for each trial
            spk = dataRun.spikes{cells}-dataRunDg.stimulus.triggers(ti(t)) ;
            spk = spk(spk>=0 & spk<=TrialTrigInterval) ;
            %psthDg{cells}{st}(t,:) = histc(spk,psthTimeDg) ;
            
            spikeTrain = zeros(1,length(psthTimeDg)) ; % empty spike train
            spkPnts = 1+floor(spk/psthBinTime) ; % spike points 
            spikeTrain(spkPnts) = 1 ; % spike train
            psthDg{cells}{st}(t,:) = smooth(spikeTrain,psthSmoothPnts) ; % psth
        end
        psthDg_mean{cells}(st,:)=mean(psthDg{cells}{st},1) ;
        psthDg_var{cells}(st,:)=var(psthDg{cells}{st},[],1) ;
    end
end

% grating vector sum
for cells=1:length(dataRun.spikes) ; % for each cell
    VectTemp = nans(length(dataRunDg.stimulus.combinations),2) ; % prep mat
    for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus 
        VectTemp(st,:) = [dataRunDg.stimulus.combinations(st).DIRECTION,...
            sum(psthDg_mean{cells}(st,ceil(StimOnOffset/diff(psthTimeDg(1:2))):ceil(DgOffTime/diff(psthTimeDg(1:2)))))] ; % prep vector (dir, spike rate)
    end
    VectSumDg(cells,:) = PolarVectorAddition(VectTemp) ;
end      

% grating tuning curves divided across speeds
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus 
    di = find(dataRunDg.stimulus.params.DIRECTION == dataRunDg.stimulus.combinations(st).DIRECTION) ;
    ti = find(dataRunDg.stimulus.params.TEMPORAL_PERIOD == dataRunDg.stimulus.combinations(st).TEMPORAL_PERIOD) ;
    for cells=1:length(dataRun.spikes) ; % for each cell
        DgTuningCurve{ti}(cells,di) = sum(psthDg_mean{cells}(st,ceil(StimOnOffset/diff(psthTimeDg(1:2))):ceil(DgOffTime/diff(psthTimeDg(1:2))))) ; %  ; % NOT FINISHED
    end
end

%% OLE (wieghts for each cell) - entire population

% direction mat
NumTpnts = length(psthDg_mean{1}(1,:)) ; % number of time points
DirMat = nans(NumTpnts*length(dataRunDg.stimulus.combinations),2) ; % prep mat
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus 
    DirX(st) = cosd(dataRunDg.stimulus.combinations(st).DIRECTION) ; % x-values of direction
    DirY(st) = sind(dataRunDg.stimulus.combinations(st).DIRECTION) ; % y-values of direction
    
    StartPnt = (st-1)*NumTpnts+1 ;
    DirMat(StartPnt:StartPnt+NumTpnts-1,1) = DirX(st) ;
    DirMat(StartPnt:StartPnt+NumTpnts-1,2) = DirY(st) ;
end

% response mat
RspMat = nans(NumTpnts*length(dataRunDg.stimulus.combinations),sum(DsCellTypeNum)) ; % prep mat
CellCnt = 1 ;
for DsType=1:length(DsCellTypeNum) ;
    for cells = 1:DsCellTypeNum(DsType) ;
        Temp = psthDg_mean{DsCelli{DsType}(cells)}' ; % transpose
        RspMat(:,CellCnt) = Temp(:) ; % linearize so [st1time1,st1time2,....stNtimeQ]
        
        CellCnt = CellCnt + 1 ;
    end
end

% OLE
Ole = RspMat\DirMat ; 


%% stimulus estimate - grating prefered vector sum on *natural images*
     
% organize response mat
for st = 1:NumDir ; % for each stimulus direction
    CellCnt = 1 ; % ds cell count
    for DsType = 1:length(DsTypeName) ; % for each ds type
        for cells=1:length(DsCelli{DsType}) ; % for each cell
            RspMatIm{st}(:,CellCnt) = psth_mean{DsCelli{DsType}(cells)}(st,:)' ;
            CellCnt = CellCnt+1 ;
        end
    end
end

% estimate stim direction with OLE
for st = 1:NumDir ; % for each stimulus direction
    DirEstXY = RspMatIm{st}*Ole ; % calculate estimate
    DirEstimate(st,:) = atan2d(DirEstXY(:,2), DirEstXY(:,1)) ; % four quadrant inverse tangent
    DirEstimate(st,DirEstimate(st,:)<0) =  DirEstimate(st,DirEstimate(st,:)<0) +360 ; % no negatives
end

% correct for grating vs image 180 deg difference
Temp = DirEstimate ;
Temp(DirEstimate<180) = DirEstimate(DirEstimate<180) + 180 ;
Temp(DirEstimate>=180) = DirEstimate(DirEstimate>=180) - 180 ;
DirEstimate = Temp ;

% error 
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        AcuteAngleError(st,t) = acuteAngle(DirEstimate(st,t),stimulus.directions(st)) ;
    end
end

%% stimulus estimate - grating prefered vector sum on *ds typing drifting grating*

% estimate stim direction with OLE
DirEstXY = RspMat*Ole ; % calculate estimate
DirEstLin = atan2d(DirEstXY(:,2), DirEstXY(:,1)) ; % four quadrant inverse tangent
DirEstLin(DirEstLin<0) = DirEstLin(DirEstLin<0)+360 ; % no negatives
DirEstimateDg = reshape(DirEstLin,[size(psthDg_mean{1},2),size(psthDg_mean{1},1)])' ; %put back by into mat(st,t)

% error 
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus 
    for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
        AcuteAngleErrorDg(st,t) = acuteAngle(DirEstimateDg(st,t),dataRunDg.stimulus.combinations(st).DIRECTION) ;
    end
end


%% histograms of error

TempError = AcuteAngleErrorDg(:,ceil(StimOnOffset/diff(psthTimeDg(1:2))):ceil(DgOffTime/diff(psthTimeDg(1:2)))) ;
histDg = hist(TempError(:),[0:3:180]) ;

TempError = AcuteAngleError(:,ceil(StimOnOffset/diff(psthTime(1:2))):end) ;
histIm = hist(TempError(:),[0:3:180]) ;

TempError = AcuteAngleErrorQuadSum(:,ceil(StimOnOffset/diff(psthTime(1:2))):end) ;
histQSumIm = hist(TempError(:),[0:3:180]) ;

TempError = AcuteAngleErrorQuadUnitAv(:,ceil(StimOnOffset/diff(psthTime(1:2))):end) ;
histQUnitAvIm = hist(TempError(:),[0:3:180]) ;

TempError = AcuteAngleErrorQuadAvNorm(:,ceil(StimOnOffset/diff(psthTime(1:2))):end) ;
histQnormIm = hist(TempError(:),[0:3:180]) ;

TempError = AcuteAngleErrorQuadAvNorm2(:,ceil(StimOnOffset/diff(psthTime(1:2))):end) ;
histQnorm2Im = hist(TempError(:),[0:3:180]) ;

TempError = AcuteAngleErrorQuadBest(:,ceil(StimOnOffset/diff(psthTime(1:2))):end) ;
histQBestIm = hist(TempError(:),[0:3:180]) ;

TempError = AcuteAngleErrorQuadLongest(:,ceil(StimOnOffset/diff(psthTime(1:2))):end) ;
histQLongestIm = hist(TempError(:),[0:3:180]) ;

TempError = AcuteAngleErrorQuadStrongest(:,ceil(StimOnOffset/diff(psthTime(1:2))):end) ;
histQStrongestIm = hist(TempError(:),[0:3:180]) ;

TempError = AcuteAngleErrorQuadLongestNorm(:,ceil(StimOnOffset/diff(psthTime(1:2))):end) ;
histQLongestNormIm = hist(TempError(:),[0:3:180]) ;

TempError = AcuteAngleErrorQuadMostSpread(:,ceil(StimOnOffset/diff(psthTime(1:2))):end) ;
histQMostSpreadIm = hist(TempError(:),[0:3:180]) ;

%% figures


