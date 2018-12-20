function ForIgor = ImSliderConcatAnalysisV4(DataBlock, DB, Params)


% this function is adapted from 'ImSliderConcatAnalysisV3' to focus on OLE
% and OQE using functions.  Also, edited to separate individual images form NI and
% Temporal periods for DG 

% JC 8/24/2017 

LoadOldMatFileFlag = false; % true skip calculations and just load
PlotFigsFlag = true ;
saveFigPath = ['/Users/jcafaro/Documents/AnalysisFigures/'] ;
numElectrodeLayers = 2 ; % number of electrode layers surrounding the max that the center is calculated over

StimOnOffset = 1 ; % (sec) time after stim begins to start counting spikes
DgOffTime = 8 ; % (sec) time ds typing drifiting grating turned off
StimOffTime = .5 ; % (sec) time before the end of image response to ignore in ole test

psthBinTime = .001 ; % (sec) size of psthNi Bin step
psthSmoothTime = 2 ; % (sec) size of psthNi Bin (smooth sliding window)

ConcatPathNum = 1 ; %TEMP -SHOULD BE INPUT
ImPathNum = 1 ;
DsPathNum = 2 ;

% parameters
Color_list = {'k','r','b','g','c','y','k','r','b','g','c','y','k','r','b','g','c','y','k','r','b','g','c','y'} ; % order of colors for each 

saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/DsIdsDb',num2str(DB),'ConcatPathNum',num2str(ConcatPathNum),'DsPathNum',num2str(DsPathNum)] ;
    
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
NumCells = length(dataRun.spikes) ; % number of cells

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
NumDirNi = size(stimulus.directionsShown,2) ; % number of directions
NumTrialsAttempted = size(stimulus.directionsShown,1) ;
NumStimShort = NumImages * NumDirNi * NumTrialsAttempted -length(triggs) ; % number of stimuli short

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
NumCompleteTrials = floor(length(triggs)/(NumDirNi*NumImages)) ; % number of trials with all directions and all images
directionsShown = stimulus.directionsShown(1:NumCompleteTrials,:) ; % take out the entire unfinished trials
triggs = triggs(1:NumCompleteTrials*NumImages*NumDirNi) ; % triggers for stim included in completed trials
 
% expand directionsShown for mutliple images
directionsShownFull = [] ;
for t = 1:NumCompleteTrials ; % for each trial
    directionsShownFull = [directionsShownFull;repmat(directionsShown(t,:),[NumImages 1])];
end

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

% psthNi
psthTimeNi = [0:psthBinTime:StimDuration] ;
psthSmoothPnts = floor(psthSmoothTime/psthBinTime) ; % number of psthBins per smooth window 
for cells=1:NumCells ; % for each cell
    for st = 1:NumDirNi ; % for each stimulus direction
        sti = find(dirShownFullTranspose(:)==stimulus.directions(st)) ; % index of triggers for that direction
        for im = 1:NumImages ; % for each image
            ti = sti([im:NumImages:length(sti)]) ; % trigger index for that image and direction 
            
            for t=1:(NumCompleteTrials) ; % for each trial
                spk = dataRun.spikes{cells}-triggs(ti(t)) ;
                spk = spk(spk>=0 & spk<=StimDuration) ;

                spikeTrain = zeros(1,length(psthTimeNi)) ; % empty spike train
                spkPnts = 1+floor(spk/psthBinTime) ; % spike points 
                spikeTrain(spkPnts) = 1 ; % spike train
                psthNi{cells}{st}{im}(t,:) = smooth(spikeTrain,psthSmoothPnts) ; % psthNi
                spkCount{cells}{st}{im}(t) = sum(spikeTrain) ; % spike count
            end
            psthNi_mean{cells}{st}(im,:)=mean(psthNi{cells}{st}{im},1) ;
            psthNi_var{cells}{st}(im,:)=var(psthNi{cells}{st}{im},[],1) ;
        end
    end
end

% % normalize psth for each image
% for cells=1:NumCells ; % for each cell
%     for im = 1:NumImages ; % for each image
%         Temp = [] ;
%         for st = 1:NumDirNi ; % for each stimulus direction
%             Temp = [Temp,mean(psthNi{cells}{st}{im},1)] ;
%         end
%         for st = 1:NumDirNi ; % for each stimulus direction
%             psthNi{cells}{st}{im} = psthNi{cells}{st}{im}/median(Temp) ;
%         end
%     end
% end
% clear Temp

%% drifting grating data
% parameters
StimOnOffsetPnts = ceil(StimOnOffset/psthBinTime) ;
DgOffPnts = ceil(DgOffTime/psthBinTime) ;

% stim path
TrialTrigInterval = 10 ;
slashi = strfind(DataBlock(DB).DsPath{DsPathNum},'/') ; % find the /
dataRunDg.triggers = dataRun.triggers(dataRun.triggers<=DgStimTime) ;
dataRunDg.names.stimulus_path = [DataBlock(DB).DsPath{Params.DsPathNum}(1:slashi(end-1)),...
    'stimuli/s',DataBlock(DB).DsPath{DsPathNum}(end-1:end),'.txt'] ;
dataRunDg = load_stim(dataRunDg,'user_defined_trigger_interval', TrialTrigInterval) ;

NumDirDg = length(dataRunDg.stimulus.params.DIRECTION) ; % number of directions
NumTempPeriods = length(dataRunDg.stimulus.params.TEMPORAL_PERIOD) ; % number of directions

% stim trials
for t=1:length(dataRunDg.stimulus.trials) ; % for each DG trial
    DgParams(t,1) = dataRunDg.stimulus.trials(t).DIRECTION ; % direction
    DgParams(t,2) = dataRunDg.stimulus.trials(t).TEMPORAL_PERIOD; % temporal period
end

% psth Dg
psthTimeDg = [0:psthBinTime:TrialTrigInterval] ;
for cells=1:NumCells ; % for each cell
    for st = 1:NumDirDg ; % for each direction 
        for tmp = 1:NumTempPeriods ; % for each temp period
            prm = [dataRunDg.stimulus.params.DIRECTION(st),dataRunDg.stimulus.params.TEMPORAL_PERIOD(tmp)] ; % params set
            ti = find((DgParams(:,1)==prm(1)).*(DgParams(:,2)==prm(2))==1) ; % index of triggers for that stim
            for t=1:length(ti) ; % for each trial
                spk = dataRun.spikes{cells}-dataRunDg.stimulus.triggers(ti(t)) ;
                spk = spk(spk>=0 & spk<=TrialTrigInterval) ;

                spikeTrain = zeros(1,length(psthTimeDg)) ; % empty spike train
                spkPnts = 1+floor(spk/psthBinTime) ; % spike points 
                spikeTrain(spkPnts) = 1 ; % spike train
                psthDg{cells}{st}{tmp}(t,:) = smooth(spikeTrain,psthSmoothPnts) ; % psthNi
                spkCountDg{cells}{st}{tmp}(t) = sum(spikeTrain) ; % spike count
            end
            psthDg_mean{cells}{st}(tmp,:)=mean(psthDg{cells}{st}{tmp},1) ;
            psthDg_var{cells}{st}(tmp,:)=var(psthDg{cells}{st}{tmp},[],1) ;
        end
    end
end

% grating vector sum
for cells=1:length(dataRun.spikes) ; % for each cell
    VectTemp = nans(NumDirDg,2) ; % prep mat
    for tmp = 1:NumTempPeriods ; % for each temp period
        for st = 1:NumDirDg ; % for each stimulus 
            VectTemp(st,:) = [dataRunDg.stimulus.params.DIRECTION(st),...
                sum(psthDg_mean{cells}{st}(tmp,StimOnOffsetPnts:DgOffPnts))] ; % prep vector (dir, spike rate)
        end
        VectSumDg{tmp}(cells,:) = PolarVectorAddition(VectTemp) ;
    end 
end      

% grating tuning curves 
for st = 1:NumDirDg ; % for each stimulus 
    for tmp = 1:NumTempPeriods ; % for each temp period
        for cells=1:NumCells ; % for each cell
            DgTuningCurve{tmp}(cells,st) = sum(psthDg_mean{cells}{st}(tmp,StimOnOffsetPnts:DgOffPnts)) ; %  
        end
    end
end

%% location subset of all cells
LocalPopFlag = false ;
LocalPopRad = 125 ; % um

SeedSpot = mean(EiCtr,1) ; % cell at center of circle
for cells=1:length(dataRun.spikes) ; % for each cell
    EiDist(cells) = sqrt(sum((SeedSpot-EiCtr(cells,:)).^2)) ;
end
    
LocalCelli = find(EiDist<LocalPopRad) ;

if LocalPopFlag ; % if only the local population
    OleAllCellsi = LocalCelli ;
    OleDsCellsi = intersect(cell2mat(DsCelli),LocalCelli) ; %DS cells
    OleNonDsCellsi = setdiff(LocalCelli,OleDsCellsi) ; % Non DS cells
else
    OleAllCellsi = [1:NumCells] ;
    OleDsCellsi = cell2mat(DsCelli) ; %DS cells
    OleNonDsCellsi = setdiff([1:NumCells],OleDsCellsi) ; % Non DS cells
end

OleCellsi = OleDsCellsi ; % indecies of cells for Ole (can choose OleAllCellsi, OleDsCellsi, or OleNonDsCellsi)
NumOleCells = length(OleCellsi) ; % number of cells to train/test ole

%% OLE/OQE - train on DG and/or NI
StimOffPnts = StimOffTime/psthBinTime ;

OqeFlag = false ; % true = use OQE, false use OLE

DgTrainFlag = true ; % train on Dg
NiTrainFlag = false ; % train on Ni

TrainingTmp_Dg = [1] ; % temporal periods (eg. [1:NumTempPeriods])
TrainingTrials_Dg = [1:size(psthDg{1}{1}{1},1)] ;
TrainingPnts_Dg = [2000:4000] ; %[1:2500] ; % randperm(size(psthDg{1}{1},2),1000) ; % randperm(DgOffPnts,500) ; %[(StimOnOffsetPnts+1):DgOffPnts]

TrainingIm_Ni = [1:NumImages] ; % images
TrainingTrials_Ni = [1:NumCompleteTrials] ; % setdiff([1:size(psthNi{1}{1},1)],[6:6:size(psthNi{1}{1},1)])
TrainingPnts_Ni = [1:2500]; % randperm(size(psthNi{1}{1},2),2000) ; %randperm(size(psthNi{1}{1},2),250) ; %[(StimOnOffsetPnts+1):length(psthNi{1}{1}(1,:))-StimOffPnts] ;

% number of trials
NumTrainingTmp_Dg = length(TrainingTmp_Dg) ;
NumTrainingIm_Ni = length(TrainingIm_Ni) ;

NumTpnts_Dg = length(TrainingPnts_Dg) ;
NumTpnts_Ni = length(TrainingPnts_Ni) ;

NumTrainingTrials_Dg = length(TrainingTrials_Dg) ;
NumTrainingTrials_Ni = length(TrainingTrials_Ni) ;

if ~DgTrainFlag ;
    NumTpnts_Dg = 0 ;
    NumTrainingTrials_Dg = 0 ;
    NumTrainingTmp_Dg = 0 ;
end

if ~NiTrainFlag ;
    NumTpnts_Ni = 0 ;
    NumTrainingTrials_Ni = 0 ;
    NumTrainingIm_Ni = 0 ;
end

stGroupNumNi = NumTpnts_Ni*NumTrainingTrials_Ni*NumTrainingIm_Ni ; % number of training 'responses' for each direction
stGroupNumDg = NumTpnts_Dg*NumTrainingTrials_Dg*NumTrainingTmp_Dg ; % number of training 'responses' for each direction
NumAllPnts = (stGroupNumDg*NumDirDg)+(stGroupNumNi*NumDirNi) ;
    
if OqeFlag ;
    disp((NumAllPnts*(NumOleCells^2-NumOleCells)/2)/10^8) ; % times bigger than possible mat size?
else
    disp((NumAllPnts*NumOleCells/2)/10^8) ; %
end

% prep mats
DirVect = nans(NumAllPnts,1) ; % prep mat 
RspMat = nans(NumAllPnts,NumOleCells) ; % prep mat

if DgTrainFlag %organize dg data
     
    % direction vector [st1time1trial1,st1time1trial2,...st1time2trial1,....stNtimeQtrialT]
    for st = 1:NumDirDg ; % for each stimulus 
        StartPnt = (st-1)*stGroupNumDg+1 ;
        DirVect(StartPnt:StartPnt+stGroupNumDg-1) = dataRunDg.stimulus.params.DIRECTION(st) ;
    end

    % response mat
    for cells = 1:NumOleCells ; % for each cell in ole
        TempRspVect = [] ; % prep mat
        for st = 1:NumDirDg ; % for each stimulus 
            for tmp = 1:NumTrainingTmp_Dg ; % for each training temp period
                Temp = psthDg{OleCellsi(cells)}{st}{TrainingTmp_Dg(tmp)}(TrainingTrials_Dg,TrainingPnts_Dg)' ; % transpose
                TempRspVect= [TempRspVect;Temp(:)] ; % linearize so [st1time1trial1,st1time1trial2,...st1time2trial1,....stNtimeQtrialT]
            end
        end
        RspMat(1:length(TempRspVect),cells) = TempRspVect ;
    end
end

if NiTrainFlag % add Ni data
        
    % direction vector
    for st = 1:NumDirNi ; % for each stimulus 
        if 180<stimulus.directions(st) ; 
           TempDir = stimulus.directions(st)-180 ; % rotate to put direction like Dg
        else
           TempDir = stimulus.directions(st)+180 ;
        end
        StartPnt = (st-1)*stGroupNumNi+1 + stGroupNumDg*NumDirDg ;
        DirVect(StartPnt:StartPnt+stGroupNumNi-1) = TempDir ;
    end

    % response mat
    for cells = 1:NumOleCells ; % for each cell in ole
        TempRspVect = [] ; % prep mat
        for st = 1:NumDirNi ; % for each stimulus
            for im = 1:NumTrainingIm_Ni ; % for each image
                Temp = psthNi{OleCellsi(cells)}{st}{TrainingIm_Ni(im)}(TrainingTrials_Ni,TrainingPnts_Ni)' ; % transpose
                TempRspVect = [TempRspVect;Temp(:)] ; % linearize so [st1time1trial1,st1time1trial2,...st1time2trial1,....stNtimeQtrialT]
            end
        end
        RspMat(stGroupNumDg*NumDirDg+1:end,cells) = TempRspVect ;
    end
end

% OLE
if OqeFlag ;
    [OleWeights, OlePolarWeights] = OleFinder(DirVect,RspMat,'OQE') ;
else
    [OleWeights, OlePolarWeights] = OleFinder(DirVect,RspMat) ;
end

%% OLE - test on NI and DG
KillDsRspFlag = false ;

TestingTmp_Dg = 1 ; % [1:NumTempPeriods] ;
TestingTrials_Dg = [1:5] ;
TestingPnts_Dg = [4000:6000]  ; % [(StimOnOffsetPnts+1):DgOffPnts] ; % setdiff([1:DgOffPnts],TrainingPnts_Dg) ;% [2600:DgOffPnts] ; % randperm(DgOffPnts,500) ; %[(StimOnOffsetPnts+1):DgOffPnts]

TestingIm_Ni = [1:NumImages] ;
TestingTrials_Ni = [1:NumCompleteTrials] ;
TestingPnts_Ni = setdiff([1:size(psthNi{1}{1}{1},2)],TrainingPnts_Ni) ;% [2600:4000]; %randperm(size(psthNi{1}{1},2),250) ; %[(StimOnOffsetPnts+1):length(psthNi{1}{1}(1,:))-StimOffPnts] ;

OleWeights_alt = OleWeights*-1 ; % rotate 180 (difference between DG and NI)

% test on NI
for st = 1:NumDirNi ; % for each stimulus 
    for im = 1:NumImages ;
        % direction estimate
        for t=1:NumCompleteTrials ; % for each trial

            RspMatTest = nans(length(psthTimeNi),NumOleCells) ; % prep mat

            for cells = 1:NumOleCells ; % for each cell in ole

                RspMatTest(:,cells) = psthNi{OleCellsi(cells)}{st}{im}(t,:) ;

                if sum(ismember(OleCellsi(cells),OleDsCellsi))>0 && KillDsRspFlag ; % if you want to kill DS responses
                    RspMatTest(:,cells) = RspMatTest(:,cells)*0 ;
                end
            end

            if OqeFlag ;
                [Temp1,Temp2] = OleTester(OleWeights_alt,RspMatTest,'OQE') ;
                DirEstimate{st}{im}(t,:) = Temp1 ;
                DirEstimate_mag{st}{im}(t,:) = Temp2 ;
            else
                [Temp1,Temp2]= OleTester(OleWeights_alt,RspMatTest) ;
                DirEstimate{st}{im}(t,:) = Temp1 ;
                DirEstimate_mag{st}{im}(t,:) = Temp2 ;
            end
        end
    
        % error
        for t=1:NumCompleteTrials; % for each trial
            for tp=1:length(DirEstimate{st}{im}(1,:)) ; % for each time point
                DirEstimate_Error{st}{im}(t,tp) = acuteAngle(DirEstimate{st}{im}(t,tp),stimulus.directions(st)) ;
            end
        end   
    end
end

% average estimate across trials for each image and direction 
for st = 1:NumDirNi ; % for each stimulus
    for im = 1:NumImages ; % for each image
        for tp=1:length(DirEstimate{st}{im}(1,:)) ; % for each time point
            Temp = PolarVectorAddition([DirEstimate{st}{im}(:,tp),ones(NumCompleteTrials,1)]) ; % average across trials
            DirEstimate_Mean{st}(im,tp) = Temp(1) ;
        end
    end
end

% average error as a function of time for each image and direction 
for st = 1:NumDirNi ; % for each stimulus  
    for im = 1:NumImages ; % for each image
        DirEstimate_Error_Mean{st}(im,:) = mean(DirEstimate_Error{st}{im},1) ;
        DirEstimate_Error_SEM{st}(im,:) = std(DirEstimate_Error{st}{im},[],1)/sqrt(NumCompleteTrials) ;
    end
end

% error histogram - error across all tested images, times and directions 
ErrorHistX = [0:3:183] ;

Temp = [] ;
for st = 1:NumDirNi ; % for each stimulus 
    for im = 1:length(TestingIm_Ni) ; % for each image to test
        Temp = [Temp,DirEstimate_Error{st}{TestingIm_Ni(im)}(TestingTrials_Ni,TestingPnts_Ni)] ; % error in steady state
    end
end
ErrorHist = hist(Temp(:),ErrorHistX) ;
ErrorMean = mean(Temp(:)) ;
ErrorStd = std(Temp(:)) ;

% error histrogram - error at most acurate time point
% find most acurate time point
Temp = [] ; % prep
for st = 1:NumDirNi ; % for each stimulus
    for im = 1:length(TestingIm_Ni) ; % for each image to test
        [m,mxi] = min(DirEstimate_Error_Mean{st}(TestingIm_Ni(im),TestingPnts_Ni)) ; % test time point with min error
        Temp = [Temp,DirEstimate_Error{st}{TestingIm_Ni(im)}(:,TestingPnts_Ni(mxi))] ; % error at that time point for each trial
    end
end
MinErrorHist = hist(Temp(:),ErrorHistX) ;
MinErrorMean = mean(Temp(:)) ;
MinErrorStd = std(Temp(:)) ;

% error histrogram - error at least acurate time point
% find least acurate time point
Temp = [] ; % prep
for st = 1:NumDirNi ; % for each stimulus
    for im = 1:length(TestingIm_Ni) ; % for each image to test
        [m,mxi] = max(DirEstimate_Error_Mean{st}(TestingIm_Ni(im),TestingPnts_Ni)) ; % test time point with max error
        Temp = [Temp,DirEstimate_Error{st}{TestingIm_Ni(im)}(:,TestingPnts_Ni(mxi))] ; % error at that time point for each trial
    end
end
MaxErrorHist = hist(Temp(:),ErrorHistX) ;
MaxErrorMean = mean(Temp(:)) ;
MaxErrorStd = std(Temp(:)) ;

% test on DG 
for st = 1:NumDirDg ; % for each stimulus
    for tmp = 1:NumTempPeriods ; % for each temp period
        for t=1:size(psthDg{1}{1}{1},1) ; % for each trial

            RspMatTest = nans(length(psthTimeDg),NumOleCells) ; % prep mat

            for cells = 1:NumOleCells ; % for each cell in ole
                RspMatTest(:,cells) = psthDg{OleCellsi(cells)}{st}{tmp}(t,:) ;

                if sum(ismember(OleCellsi(cells),OleDsCellsi))>0 && KillDsRspFlag ; % if you want to kill DS responses
                    RspMatTest(:,cells) = RspMatTest(:,cells)*0 ;
                end    
            end
            if OqeFlag ;
                DirEstimateDg{st}{tmp}(t,:) = OleTester(OleWeights,RspMatTest,'OQE') ;
            else     
                DirEstimateDg{st}{tmp}(t,:) = OleTester(OleWeights,RspMatTest) ;
            end
        end

        % error
        for t=1:size(psthDg{1}{1}{1},1) ; % for each trial
            for tp=1:length(DirEstimateDg{st}{tmp}(t,:)) ; % for each time point
                DirEstimateDg_Error{st}{tmp}(t,tp) = acuteAngle(DirEstimateDg{st}{tmp}(t,tp),dataRunDg.stimulus.params.DIRECTION(st)) ;
            end
        end   
    end
end

% average estimate across trials for each image and direction 
for st = 1:NumDirDg ; % for each direction
    for tmp = 1:NumTempPeriods ; % for temp period
        for t=1:length(DirEstimateDg{st}{tmp}(1,:)) ; % for each time point
            Temp = PolarVectorAddition([DirEstimateDg{st}{tmp}(:,t),ones(size(DirEstimateDg{st}{tmp}(:,t)))]) ; % average across trials
            DirEstimateDg_Mean{st}(tmp,t) = Temp(1) ;
        end
    end
end

% average error as a function of time for each temp period and direction 
for st = 1:NumDirDg ; % for each stimulus 
    for tmp = 1:NumTempPeriods ;
        DirEstimateDg_Error_Mean{st}(tmp,:) = mean(DirEstimateDg_Error{st}{tmp},1) ;
        DirEstimateDg_Error_SEM{st}(tmp,:) = std(DirEstimateDg_Error{st}{tmp},[],1)/sqrt(size(psthDg{1}{1}{1},1)) ; % CHECK
    end
end

% error histogram for all tested temp periods, times and directions
Temp = [] ;
for st = 1:NumDirDg ; % for each stimulus 
    for tmp = 1:length(TestingTmp_Dg) ; % for each tested tmp
        Temp = [Temp,DirEstimateDg_Error{st}{TestingTmp_Dg(tmp)}(TestingTrials_Dg,TestingPnts_Dg)] ; % error in steady state
    end
end
ErrorHistDg = hist(Temp(:),ErrorHistX) ;
ErrorMeanDg = mean(Temp(:)) ;
ErrorStdDg = std(Temp(:)) ;

% error histrogram - error at most acurate time point
% find most acurate time point
Temp = [] ; % prep
for st = 1:NumDirNi ; % for each stimulus
    for tmp = 1:length(TestingTmp_Dg) ;
        [m,mxi] = min(DirEstimateDg_Error_Mean{st}(TestingTmp_Dg(tmp),TestingPnts_Dg)) ; % test time point with min error
        Temp = [Temp,DirEstimateDg_Error{st}{TestingTmp_Dg(tmp)}(:,TestingPnts_Dg(mxi))] ; % error at that time point for each trial
    end
end
MinErrorHistDg = hist(Temp(:),ErrorHistX) ;
MinErrorMeanDg = mean(Temp(:)) ;
MinErrorStdDg = std(Temp(:)) ;

% error histrogram - error at least acurate time point
% find least acurate time point
Temp = [] ; % prep
for st = 1:NumDirNi ; % for each stimulus
    for tmp = 1:length(TestingTmp_Dg) ;
        [m,mxi] = max(DirEstimateDg_Error_Mean{st}(TestingTmp_Dg(tmp),TestingPnts_Dg)) ; % test time point with max error
        Temp = [Temp,DirEstimateDg_Error{st}{TestingTmp_Dg(tmp)}(:,TestingPnts_Dg(mxi))] ; % error at that time point for each trial
    end
end
MaxErrorHistDg = hist(Temp(:),ErrorHistX) ;
MaxErrorMeanDg = mean(Temp(:)) ;
MaxErrorStdDg = std(Temp(:)) ;


%% crossterm analysis for OQE
if OqeFlag ; 
    
    % which crossterms belong to which cell
    CrossTermCelli = [] ;
    for cells = 2:NumOleCells ;
        CrossTermCelli = [CrossTermCelli,[repmat(cells,[1,cells-1]);[1:cells-1]]] ; % cellA;cellB in order of crossterms
    end

    % direction of cross terms 
    for cells = 1:NumOleCells ; % for each cell in the OQE
        [temp,cti] = find(CrossTermCelli==cells) ; % find cross terms with that cell 
        cti = cti + NumOleCells ; % adjust for first set of terms which are not cross terms - index of OleWeights

        CrsTermVecSum(cells,:) = PolarVectorAddition(OlePolarWeights(cti,:)) ; % vector sum of cross terms
        
        TermCrsTermAngDiff(cells) = acuteAngle(CrsTermVecSum(cells,1),OlePolarWeights(cells,1))
    end
    
    % strength of cross terms
     for cells = 1:NumOleCells ; % for each cell in the OQE
        [temp,cti] = find(CrossTermCelli==cells) ; % find cross terms with that cell 
        cti = cti + NumOleCells ; % adjust for first set of terms which are not cross terms - index of OleWeights

        CrsTermsMag{cells} = OlePolarWeights(cti,2) ; % magnitude
        for c = 1:length(cti) ; % for each cross cell
            CrsTermsDirDiff{cells}(c) = acuteAngle(OlePolarWeights(cells,1),OlePolarWeights(cti(c),1)) ; % direciton difference
            CrsTermsDot{cells}(c) = OleWeights(cells,:)*OleWeights(cti(c),:)' ; % dot product of cross term and single cell wieghts 
        end
        
        ci = setdiff([1:NumOleCells],cells) ; % order of cells 
        for c2 = 1:length(ci) ;
            CrsTermsPairDist{cells}(c2) = sqrt(sum((EiCtr(OleCellsi(cells),:) - EiCtr(OleCellsi(ci(c2)),:)).^2)) ;
        end  
     end
end

%% figures

figure % mean direction estimate - NI
for st = 1:NumDirNi ; % for each stimulus
    for im = 1:NumImages ;
        plot(psthTimeNi,DirEstimate_Mean{st}(im,:),Color_list{st})
        hold on
        plot(psthTimeNi,ones(size(psthTimeNi))*stimulus.directions(st),[Color_list{st},'--'])
         pause; hold off
    end
end

figure % mean direction estimate - DG
for st = 1:NumDirDg ; % for each stimulus
    for tmp = 1:NumTempPeriods ;
        plot(psthTimeDg,DirEstimateDg_Mean{st}(tmp,:),Color_list{st})
        hold on
        plot(psthTimeDg,ones(size(psthTimeDg))*dataRunDg.stimulus.params.DIRECTION(st),[Color_list{st},'--'])
        %pause; hold off
    end
end

figure % direction estimates on each trial - NI
for st = 1:NumDirNi ; % for each stimulus
    for im=1:NumImages ; % for each image
        for t=1:NumCompleteTrials ;
            subplot(2,1,1)
            plot(psthTimeNi,DirEstimate{st}{im}(t,:))
            hold on
            plot(psthTimeNi,ones(size(psthTimeNi))*stimulus.directions(st),'k--')
            hold off

            subplot(2,1,2)
            plot(psthTimeNi,DirEstimate_Error{st}{im}(t,:))
            pause
        end
    end
end

figure % direction estimates on each trial - DG
for st = 1:NumDirNi ; % for each stimulus
    for tmp = 1:NumTempPeriods ;
        for t = 1:size(psthDg{1}{st}{tmp},1) ; % for each trial
            subplot(2,1,1)
            plot(psthTimeDg,DirEstimateDg{st}{tmp}(t,:))
            hold on
            plot(psthTimeDg,ones(size(psthTimeDg))*dataRunDg.stimulus.params.DIRECTION(st),'k--')
            hold off

            subplot(2,1,2)
            plot(psthTimeDg,DirEstimateDg_Error{st}{tmp}(t,:))
            pause
        end
    end
end

figure % error - NI
subplot(2,1,1) % time dependant error 
for st = 1:NumDirNi ;
    for im = 1:NumImages ; 
        errorbar(psthTimeNi,DirEstimate_Error_Mean{st}(im,:),DirEstimate_Error_SEM{st}(im,:),'k')
        %pause
    end
end
xlabel('time (sec)')
ylabel('error (deg)')

subplot(2,1,2) % error histogram
plot(ErrorHistX,ErrorHist/sum(ErrorHist),'k')
hold on
plot(ErrorHistX,MinErrorHist/sum(MinErrorHist),'k--')
plot(ErrorHistX,MaxErrorHist/sum(MaxErrorHist),'k:')
xlabel('error (deg)')
ylabel('fraction')

figure % error - DG

subplot(2,1,1) % time dependant error 
for st = 1:NumDirNi ;
    for tmp = 1:NumTempPeriods ; 
        errorbar(psthTimeDg,DirEstimateDg_Error_Mean{st}(tmp,:),DirEstimateDg_Error_SEM{st}(tmp,:),'r')
        %pause
    end
end
xlabel('time (sec)')
ylabel('error (deg)')

subplot(2,1,2) % error histogram
plot(ErrorHistX,cumsum(ErrorHistDg/sum(ErrorHistDg)),'r')
hold on
plot(ErrorHistX,cumsum(MinErrorHistDg/sum(MinErrorHistDg)),'r--')
xlabel('error (deg)')
ylabel('fraction')

figure % wieghts
subplot(3,1,1) % weights for each cell
plot([1:NumOleCells],OlePolarWeights(1:NumOleCells,2),'bo')
hold on
plot(OleDsCellsi,OlePolarWeights(OleDsCellsi,2),'r*')
plot([NumOleCells+1:length(OlePolarWeights(:,2))],OlePolarWeights(NumOleCells+1:end,2),'ro')
xlabel('cell and pairs')
ylabel('weight magnitude')

subplot(3,1,2) % histogram of weight magnitudes
Temp = hist(OlePolarWeights(1:NumOleCells,2),10.^[-2:.1:2]) ;
Temp2 = hist(OlePolarWeights(NumOleCells+1:end,2),10.^[-2:.1:2]) ;
plot(10.^[-2:.1:2],Temp/sum(Temp),'b--')
hold on
plot(10.^[-2:.1:2],Temp2/sum(Temp2),'r--')
ylabel('fraction')
xlabel('weight magnitude (fraction)')
set(gca,'xscale','log')

subplot(3,1,3) % histogram of weight directions
Temp = hist(OlePolarWeights(1:NumOleCells,1),[0:10:360]) ;
Temp2 = hist(OlePolarWeights(NumOleCells+1:end,1),[0:10:360]) ;
plot([0:10:360],Temp/sum(Temp),'b--')
hold on
plot([0:10:360],Temp2/sum(Temp2),'r--')
xlabel('weight direction')
ylabel('fraction')
legend('cell weights','pair weights')

figure % polar weights
polar(0,max(OlePolarWeights(cells,2))) ;
for cells = 1:size(OlePolarWeights,1) ; % for each cell
    polar([0,OlePolarWeights(cells,1)*pi/180],[0,OlePolarWeights(cells,2)])
    hold on
end

figure % cross terms wieghts - all cells
for cells = 1:NumOleCells ;
    [temp,cti] = find(CrossTermCelli==cells) ; % find cross terms with that cell
    cti = cti + NumOleCells ; % adjust for first set of terms which are not cross terms
    
    ctVs = PolarVectorAddition(OlePolarWeights(cti,:)) ; % vector sum of cross terms
    
    clf ;
    polar(0,max(OlePolarWeights(cti,2)),'b') ; % adjust axis so longest vector can be seen
    hold on
    for cells2 = 1:length(cti) ; % plot cross terms
        polar([0,OlePolarWeights(cti(cells2),1)*pi/180],[0,OlePolarWeights(cti(cells2),2)],'b')
    end
    polar([0,OlePolarWeights(cells,1)*pi/180],[0,OlePolarWeights(cells,2)],'r') ; % plot cell
    polar([0,ctVs(1)*pi/180],[0,ctVs(2)],'g') ; % plot cell
    pause
end
    
figure % cross terms wieghts - DS cells only
for cells = 1:length(OleDsCellsi) ;
    ci = find(OleCellsi == OleDsCellsi(cells)) ; % index of OleCellsi
    [temp,cti] = find(CrossTermCelli==ci) ; % find cross terms with that cell
    cti = cti + NumOleCells ; % adjust for first set of terms which are not cross terms - index of OleWieghts
    
    ctVs = PolarVectorAddition(OlePolarWeights(cti,:)) ; % vector sum of cross terms
    
    clf ;
    polar(0,max(OlePolarWeights(cti,2)),'b') ; % adjust axis so longest vector can be seen
    hold on
    for cells2 = 1:length(cti) ; % plot cross terms
        polar([0,OlePolarWeights(cti(cells2),1)*pi/180],[0,OlePolarWeights(cti(cells2),2)],'b')
    end
    polar([0,OlePolarWeights(ci,1)*pi/180],[0,OlePolarWeights(ci,2)],'r') ; % plot cell
    polar([0,ctVs(1)*pi/180],[0,ctVs(2)],'g') ; % plot cell
    pause
end

figure % cross terms vs single cell term histograms
subplot(2,1,1) % all cells
hist(TermCrsTermAngDiff,[0:10:200])
xlabel('Angle Cross Term and single cell') 
ylabel('number of cells obs')

subplot(2,1,2) % ds cells only
hist(TermCrsTermAngDiff(ismember(OleCellsi,OleDsCellsi)),[0:10:200])
xlabel('Angle Cross Term and single cell') 
ylabel('number of DS cells obs')

figure % direction and mag of crossterm as function of distance - DS cells only
subplot(3,1,1) % magnitude
for cells = 1:length(OleDsCellsi) ; % for each DS cell
    ci = find(OleCellsi == OleDsCellsi(cells)) ; % index of OleCellsi
    plot(CrsTermsPairDist{ci},CrsTermsMag{ci},'k*')
    hold on
end
xlabel('distance between cells (um)')
ylabel('crossterm magnitude')

subplot(3,1,2) % direction difference
for cells = 1:length(OleDsCellsi) ; % for each DS cell
    ci = find(OleCellsi == OleDsCellsi(cells)) ; % index of OleCellsi
    plot(CrsTermsPairDist{ci},CrsTermsDirDiff{ci},'k*')
    hold on
end
xlabel('distance between cells (um)')
ylabel('direction difference between crossterm and single cell')

subplot(3,1,3)
for cells = 1:length(OleDsCellsi) ; % for each DS cell
    ci = find(OleCellsi == OleDsCellsi(cells)) ; % index of OleCellsi
     plot(CrsTermsPairDist{ci},CrsTermsDot{ci},'k*')
    hold on
end
xlabel('distance between cells (um)')
ylabel('dot product difference between crossterm and single cell')

figure % error vs. estimate magnitude
for t=1:NumCompleteTrials ; % for each trial
    for st = 1:NumDirNi ; % for each stimulus
        plot(DirEstimate_Error{st}(t,:),DirEstimate_mag{st}(t,:),'.')

        hold on
    end
end


figure % DS cell EI positions  
for DsType = 1:length(DsTypeName) ;
    %subplot(length(DsTypeName),1,DsType)
    %figure
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        plot(EiCtr(DsCelli{DsType}(cells),1),EiCtr(DsCelli{DsType}(cells),2),[Color_list{DsType},'+'])
        drawCircle(EiCtr(DsCelli{DsType}(cells),1),EiCtr(DsCelli{DsType}(cells),2),50,'color',Color_list{DsType})
        hold on
    end
    axis([-500 500 -500 500])
end
drawCircle(SeedSpot(1),SeedSpot(2),LocalPopRad,'color','c')

figure % Dg tuning curves of cells used in Ole estimate
tmp = 1 ; % temporal period
for cells = 1:length(OleDsCellsi) ; % for each DS cell
    %plot(dataRunDg.stimulus.params.DIRECTION,DgTuningCurve{1}(OleDsCellsi(cells),:))
    polar(VectSumDg{tmp}(OleDsCellsi(cells),1)*pi/180,1,'*')
    hold on
end

