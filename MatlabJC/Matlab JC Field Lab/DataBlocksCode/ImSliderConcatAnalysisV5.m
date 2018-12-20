function ForIgor = ImSliderConcatAnalysisV5(DataBlock, DB, Params)


% this function is adapted from 'ImSliderConcatAnalysisV4' to test and train on NI
% only and train always outside of test image.  Also added extras
% calculations for conditional histograms, quad response histograms, correlation by distance. 

% added ability to run as a script and call DataBlocks

% JC 11/14/2017 
RunAsScript = true ;
if RunAsScript ;
    DB = 15 ;
    [DataBlock,Params] = DataBlocks_NaturalStim ;
end

% parameters - CAUTION flags distributed throughout 
ConcatPathNum = 1 ; %
ImPathNum = 1 ;
DsPathNum = 1 ;

psthBinTime = .001 ; % (sec) size of psthNi Bin step
psthSmoothTime = .02 ; % (sec) size of psthNi Bin (smooth sliding window)

StimOnOffset = 1 ; % (sec) time after stim begins to start counting spikes
DgOffTime = 8 ; % (sec) time ds typing drifiting grating turned off
StimOffTime = .5 ; % (sec) time before the end of image response to ignore in ole test
numElectrodeLayers = 2 ; % number of electrode layers surrounding the max that the center is calculated over

Color_list = {'k','r','b','g','c','y','k','r','b','g','c','y','k','r','b','g','c','y','k','r','b','g','c','y'} ; % order of colors for each 
saveFigPath = ['/Users/jcafaro/Documents/AnalysisFigures/'] ;
saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/DsIdsDb',num2str(DB),'ConcatPathNum',num2str(ConcatPathNum),'DsPathNum',num2str(DsPathNum)] ;
    
% parameter time-->pnts
StimOnOffsetPnts = ceil(StimOnOffset/psthBinTime) ;
StimOffPnts = StimOffTime/psthBinTime ;
DgOffPnts = ceil(DgOffTime/psthBinTime) ;

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
            NewTrigg = (triggs(t)-triggs(t-1))/2 ;
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
SteadyStatePntsNi = [(StimOnOffsetPnts+1):length(psthTimeNi)-StimOffPnts] ;
SteadyStatePntsDg = [(StimOnOffsetPnts+1):DgOffPnts] ;

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
                psthNi{cells}{st}{im}(t,:) = smooth(spikeTrain,psthSmoothPnts)/psthBinTime ; % psthNi (hz)
                spkCount{cells}{st}{im}(t) = sum(spikeTrain) ; % spike count
            end
            psthNi_mean{cells}{st}(im,:)=mean(psthNi{cells}{st}{im},1) ; % average across trials
            psthNi_var{cells}{st}(im,:)=var(psthNi{cells}{st}{im},[],1) ;
        end
    end
end

% firing rate distribution of all DS cells - NI
psthNiHistX = [0:1:psthSmoothPnts]/(psthSmoothPnts*psthBinTime) ; % possible spike rates
psthNiHist = psthNiHistX*0 ; % prep mat
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(DsCelli{DsType}) ; % for each cell
        for st = 1:NumDirNi ; % for each stimulus direction
            for im = 1:NumImages ; % for each image
                psthNiHist = psthNiHist+ hist(psthNi{DsCelli{DsType}(cells)}{st}{im}(:),psthNiHistX) ;
            end
        end
    end
end
psthNiHist = psthNiHist/sum(psthNiHist) ;

% Ni tuning curves (averaged over images, trials, time)
for st = 1:NumDirNi ; % for each stimulus 
    for cells=1:NumCells ; % for each cell
        Temp = psthNi_mean{cells}{st}(:,SteadyStatePntsNi) ; %  
        NiTuningCurve(cells,st) = mean(Temp(:)) ;
    end
end
clear Temp

%% drifting grating data

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
                psthDg{cells}{st}{tmp}(t,:) = smooth(spikeTrain,psthSmoothPnts)/psthBinTime ; % psthDg (hz)
                spkCountDg{cells}{st}{tmp}(t) = sum(spikeTrain) ; % spike count
            end
            psthDg_mean{cells}{st}(tmp,:)=mean(psthDg{cells}{st}{tmp},1) ; % average across trials
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

%% find quads and calculate quad response histograms 
QuadEiCircRad = 200 ; % minimum pairwise distance to be in quad

% get ei center of mass and make TypeCellMat so can figure out which center
% belongs to which cell
CellCnt = 1 ;
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(DsCelli{DsType}) ; % for each cell
        ctr(CellCnt,:) = EiCtr(DsCelli{DsType}(cells),:) ;
        
        TypeCellMat(CellCnt,:) = [DsType,cells] ;

        CellCnt = CellCnt+1 ;
    end
end

% distance between cells
CellCnt = 1 ; % start cell counter
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(DsCelli{DsType}) ; % for each cell
        deltaCtr = ctr-repmat(ctr(CellCnt,:),size(ctr,1),1) ; % X,Y distance between that cell and all others
        distCtr{CellCnt} = sqrt(deltaCtr(:,1).^2+deltaCtr(:,2).^2) ; % euclidean distance
        [~,si(CellCnt,:)] = sort(distCtr{CellCnt}) ; % in order of nearest to farthest

        CellCnt = CellCnt+1 ;
    end
end

% find potential quadruplets
CellCnt = 1 ; % start cell counter
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(DsCelli{DsType}) ; % for each cell
        quad{CellCnt} =[DsType,cells] ; % all the cells in this quadruplet
         for NearCell = 2:size(si,2) ; % for each ds cell
             Ni = si(CellCnt,NearCell) ; % near cell index
             psth_i = TypeCellMat(Ni,:) ;
             if distCtr{CellCnt}(Ni)<QuadEiCircRad && ~ismember(psth_i(1),quad{CellCnt}(:,1)) ; % if the cell is very nearby and that type has not yet part of the quad
                quad{CellCnt}=[quad{CellCnt};psth_i] ; % add it to the quad     
             end
         end
         CellCnt = CellCnt+1 ;
    end
end
       
% identify real quadruplets (4 different type ds cells with max pairwise dist<QuadEiCircRad)
RealQuad_i = [] ; % index of quadruplets in quad

for q=1:length(quad) ; % for each putative quadruplet
    if size(quad{q},1)==4 ; % if it is four cells
        for cells=1:4 ;
            qloc(cells,:) = EiCtr(DsCelli{quad{q}(cells,1)}(quad{q}(cells,2)),:) ; % location of each cell in quad
        end
        
        if max(pdist(qloc))<QuadEiCircRad ; % if all pairwise distances are below threshold
            urq=true ; % assume it is a unique real quad 
            if ~isempty(RealQuad_i) ; % if there are any real quadruplets
                for rq=RealQuad_i ; % for each real quad 
                    [~,rqOrder] = sort(quad{rq}(:,1)) ;
                    [~,qOrder] = sort(quad{q}(:,1)) ;
                    if quad{rq}(rqOrder,:) == quad{q}(qOrder,:); % if its the same as this real quad
                        urq=false ;
                    end
                end
            end
      
            if urq ; % if its a unique real quad
                RealQuad_i = [RealQuad_i,q] ;
            end
        end
    end
end

% calc quadruplet stats
quadNumCellsHistX = [0:4] ; % number of cells bins
quadRateRatioHistX = [.1:.1:1] ; % rate ratio bins
if ~isempty(RealQuad_i) ; % if there are any quadruplets
    quadNiNumCellsHist = quadNumCellsHistX*0 ; % prep mat
    quadNiRateRatioHist = quadRateRatioHistX*0 ; % prepMat
    quadDgNumCellsHist = quadNumCellsHistX*0 ; % prep mat
    quadDgRateRatioHist = quadRateRatioHistX*0 ; % prepMat

    for q=RealQuad_i ; % for each real quad
        
        psthBlock = [] ;
        for st = 1:NumDirNi ; % for each stimulus direction NI
            for im=1:NumImages ; % for each image
                for cells=1:4 ; % for each cell in quad
                    Temp = psthNi{DsCelli{quad{q}(cells,1)}(quad{q}(cells,2))}{st}{im}(:,SteadyStatePntsNi) ; % concat
                    psthBlock(cells,:) = Temp(:) ; % linearize
                end
                quadNiNumCellsHist = quadNiNumCellsHist + hist(sum(psthBlock>0,1),quadNumCellsHistX) ; % how many cells are firing at the same time (0-4)
                quadNiRateRatioHist = quadNiRateRatioHist + hist([max(psthBlock,[],1)./sum(psthBlock,1)],quadRateRatioHistX) ; % ratio
            end
        end

        psthBlock = [] ;
        for st = 1:NumDirDg ; % for each stimulus direction 
            for tmp=1:NumTempPeriods ; % for each temp period
                for cells=1:4 ; % for each cell in quad
                    Temp = psthDg{DsCelli{quad{q}(cells,1)}(quad{q}(cells,2))}{st}{tmp}(:,SteadyStatePntsDg) ; % concat
                    psthBlock(cells,:) = Temp(:) ;
                end
                quadDgNumCellsHist = quadDgNumCellsHist + hist(sum(psthBlock>0,1),quadNumCellsHistX) ; % how many cells are firing at the same time (0-4)
                quadDgRateRatioHist = quadDgRateRatioHist + hist([max(psthBlock,[],1)./sum(psthBlock,1)],quadRateRatioHistX) ; % ratio
            end
        end
    end
    
    %change hist to pdfs
    quadNiNumCellsHist = quadNiNumCellsHist/sum(quadNiNumCellsHist) ; % 
    quadNiRateRatioHist = quadNiRateRatioHist/sum(quadNiRateRatioHist); %
    quadDgNumCellsHist = quadDgNumCellsHist/sum(quadDgNumCellsHist) ; %
    quadDgRateRatioHist = quadDgRateRatioHist/sum(quadDgRateRatioHist) ; 
end


for TEMP_LOOP = [1:7] ;
%% select cells for OLE
LocalPopFlag = false ; % use local population
ExcludeOnDSFlag = true ; 
UseQuadFlag = false ; % use single quad only
RealQuadToUse = 1 ; % real quad i to use (could be 'TEMP_LOOP')
UseDsCellsOnlyFlag = false ; 
UseNonDsCellsOnlyFlag = false ;
LocalPopRad = 70 ; % um
SeedSpot = [0,0] ; % mean(EiCtr,1) ; % center of circle

clear OleCellsi ; % in case its been run before
if UseQuadFlag % use a quad
    RandomQuadFlag = false ; % construct a random quadruplet

    if ~RandomQuadFlag ; % if you want to use a real quad
        for cells=1:4 ; % for each cell
            OleCellsi(cells) = DsCelli{quad{RealQuad_i(RealQuadToUse)}(cells,1)}...
                (quad{RealQuad_i(RealQuadToUse)}(cells,2)) ;
        end
    else
        for cells=1:4 ; % for each cell
            OleCellsi(cells) = DsCelli{cells}(randi(length(DsCelli{cells}))) ; % pick randomly for each type
        end
    end
else % if not using a quad
    if LocalPopFlag ; % if only the local population
        % local group
        for cells=1:length(dataRun.spikes) ; % for each cell
            EiDist(cells) = sqrt(sum((SeedSpot-EiCtr(cells,:)).^2)) ;
        end
        LocalCelli = find(EiDist<LocalPopRad) ;
        
        OleAllCellsi = LocalCelli ;
        OleDsCellsi = intersect(cell2mat(DsCelli),LocalCelli) ; %DS cells
        OleNonDsCellsi = setdiff(LocalCelli,OleDsCellsi) ; % Non DS cells
    else
        OleAllCellsi = [1:NumCells] ;
        OleDsCellsi = cell2mat(DsCelli) ; %DS cells
        OleNonDsCellsi = setdiff([1:NumCells],OleDsCellsi) ; % Non DS cells
    end

    if UseDsCellsOnlyFlag ; % use DS cells only 
        OleCellsi = OleDsCellsi ; 
    elseif UseNonDsCellsOnlyFlag ; % use non-DS cells only
        OleCellsi = OleNonDsCellsi ; 
    else % use both DS and NonDS
        OleCellsi = OleAllCellsi ; 
    end

    if ExcludeOnDSFlag ; % if you want to exclude On-DS cells
        OleCellsi = setdiff(OleCellsi,get_cell_indices(dataRun,cell2mat(ForIgor.ds_id{1}))) ;
    end
end

NumOleCells = length(OleCellsi) ; % number of cells to train/test ole

%% OLE/OQE - train and test on NI only
OqeFlag = true ; % true = use OQE, false use OLE
RegularizeOqeFlag = false ; % true = regularize OQE
RegLambda = 1 ; % value used for regularization
ShuffTrainFlag = false ; % shuffle the responses for the training data
ShuffTestFlag = false ; % shuffle the responses for the test set
ShuffTrialTestFlag = false ; % shuffle trials
BestCellTestFlag = false ; % use only the strongest responding cell at any time point in estimate

KillNonDsSingleTermFlag  = false ;
KillDsSingleTermFlag = false ;
KillNonDsNonDsCrossTermFlag = false ;
KillDsNonDsCrossTermFlag = false ;
KillDsDsCrossTermFlag = false ;

% crossterms indicies
CrossTermCelli = [] ;
for cells = 2:NumOleCells ;
    CrossTermCelli = [CrossTermCelli,[repmat(cells,[1,cells-1]);[1:cells-1]]] ; % cellA;cellB in order of crossterms
end
CrossTermCellOlei = OleCellsi(CrossTermCelli) ;

CtKilli = [] ;
if KillNonDsSingleTermFlag ;
    CtKilli = [CtKilli,find(ismember(OleCellsi,OleNonDsCellsi))] ;
end

if KillDsSingleTermFlag ;
    CtKilli = [CtKilli,find(ismember(OleCellsi,OleDsCellsi))] ;
end

if KillNonDsNonDsCrossTermFlag ;
    CtKilli = [CtKilli,find(sum(ismember(CrossTermCellOlei,OleNonDsCellsi),1)==2)+NumOleCells];
end

if KillDsNonDsCrossTermFlag ;
    CtKilli = [CtKilli,find(sum(ismember(CrossTermCellOlei,OleNonDsCellsi),1)==1)+NumOleCells] ;
end
 
if KillDsDsCrossTermFlag ;
    CtKilli = [CtKilli,find(sum(ismember(CrossTermCellOlei,OleNonDsCellsi),1)==0)+NumOleCells] ;
end

for im = 1:NumImages ; % for each Ni image 

    TrainingIm_Ni = im ; % images: use setdiff([1:NumImages],im) for ole ; use 1 or im for OQE
    TrainingTrials_Ni = [5] ;% [1:NumCompleteTrials] ; % setdiff([1:size(psthNi{1}{1},1)],[6:6:size(psthNi{1}{1},1)])
    TrainingPnts_Ni = SteadyStatePntsNi(randperm(length(SteadyStatePntsNi),1000)) ; % use SteadyStatePntsNi for ole; % randperm(size(psthNi{1}{1}{1},2)-StimOnOffsetPnts,1000)+StimOnOffsetPnts; % randperm(size(psthNi{1}{1},2),2000) ; %[(StimOnOffsetPnts+1):length(psthNi{1}{1}(1,:))-StimOffPnts] ;

    % number of trials
    NumTrainingIm_Ni = length(TrainingIm_Ni) ;
    NumTpnts_Ni = length(TrainingPnts_Ni) ;
    NumTrainingTrials_Ni = length(TrainingTrials_Ni) ;

    stGroupNumNi = NumTpnts_Ni*NumTrainingTrials_Ni*NumTrainingIm_Ni ; % number of training 'responses' for each direction
    NumAllPnts = stGroupNumNi*NumDirNi ;

    if OqeFlag ;
        disp((NumAllPnts*(NumOleCells^2-NumOleCells)/2)/10^8) ; % times bigger than 'possible' mat size?
    else
        %disp((NumAllPnts*NumOleCells/2)/10^8) ; %
    end

    % prep mats
    DirVect = nans(NumAllPnts,1) ; % prep mat 
    RspMat = nans(NumAllPnts,NumOleCells) ; % prep mat

    % direction vector
    for st = 1:NumDirNi ; % for each stimulus 
        if 180<stimulus.directions(st) ; 
           TempDir = stimulus.directions(st)-180 ; % rotate to put direction like Dg
        else
           TempDir = stimulus.directions(st)+180 ;
        end
        StartPnt = (st-1)*stGroupNumNi+1 ;
        DirVect(StartPnt:StartPnt+stGroupNumNi-1) = TempDir ;
    end

    % response mat
    for cells = 1:NumOleCells ; % for each cell in ole
        TempRspVect = [] ; % prep mat
        for st = 1:NumDirNi ; % for each stimulus
            for imTemp = 1:NumTrainingIm_Ni ; % for each image
                Temp = psthNi{OleCellsi(cells)}{st}{TrainingIm_Ni(imTemp)}(TrainingTrials_Ni,TrainingPnts_Ni)' ; % transpose
                TempRspVect = [TempRspVect;Temp(:)] ; % linearize so [st1time1trial1,st1time1trial2,...st1time2trial1,....stNtimeQtrialT]
            end
        end
        RspMat(1:end,cells) = TempRspVect ;
    end

    if ShuffTrainFlag ; % if Shuffle training flag - circshift shuffle within image 
        for cells = 1:NumOleCells ; % for each cell in ole
            TempRspVect = [] ; % prep mat
            for st = 1:NumDirNi ; % for each stimulus
                for imTemp = 1:NumTrainingIm_Ni ; % for each image
                    ShuffPntsTemp = circshift(TrainingPnts_Ni,[0,randi(length(TrainingPnts_Ni))]) ; % same shuffle for all trials
                    Temp = psthNi{OleCellsi(cells)}{st}{TrainingIm_Ni(imTemp)}(TrainingTrials_Ni,ShuffPntsTemp)' ; % transpose
                    TempRspVect = [TempRspVect;Temp(:)] ; % linearize so [st1time1trial1,st1time1trial2,...st1time2trial1,....stNtimeQtrialT]
                end
            end
            RspMat(1:end,cells) = TempRspVect ;
        end 
        
        % orig code - shuffles across images
%         DirVectUnique = unique(DirVect) ;
%         for st = 1:NumDirNi ; % for each direction
%             for cells = 1:NumOleCells ; % for each cell in ole
%                 Tempi = find(DirVect==DirVectUnique(st)) ;
%                 RspMat(Tempi,cells) = RspMat(Tempi(randperm(length(Tempi))),cells) ; % shuffle keeping within direction and cell
%             end
%         end
    end
    
    % OLE
    if OqeFlag ;
        if ~RegularizeOqeFlag ;
            [OleWeights, OlePolarWeights] = OleFinder(DirVect,RspMat,'OQE') ;
        else
            [OleWeights, OlePolarWeights] = OleFinder(DirVect,RspMat,'OQE','Regularize',RegLambda) ;
        end
    else
        [OleWeights, OlePolarWeights] = OleFinder(DirVect,RspMat) ;
    end

    % OLE - test on NI and DG
    KillDsRspFlag = false ;

    TestingIm_Ni = im ;
    TestingTrials_Ni = [1:NumCompleteTrials] ;
    TestingPnts_Ni = setdiff(SteadyStatePntsNi,TrainingPnts_Ni); % SteadyStatePntsNi ; % setdiff([1:size(psthNi{1}{1}{1},2)],TrainingPnts_Ni) ;% [2600:4000]; %randperm(size(psthNi{1}{1},2),250) ; %[(StimOnOffsetPnts+1):length(psthNi{1}{1}(1,:))-StimOffPnts] ;

    OleWeights_alt = OleWeights*-1 ; % rotate 180 (difference between DG and NI)
    
    % set flagged terms to zero
    OleWeights_altTest = OleWeights_alt ;
    OleWeights_altTest(CtKilli,:) = 0 ;
    
    % test on NI
    for st = 1:NumDirNi ; % for each direction 

        % direction estimate
        for t=1:NumCompleteTrials ; % for each trial

            RspMatTest = nans(length(psthTimeNi),NumOleCells) ; % prep mat

            for cells = 1:NumOleCells ; % for each cell in ole

                RspMatTest(:,cells) = psthNi{OleCellsi(cells)}{st}{im}(t,:) ;

                if ShuffTrialTestFlag ; % if you want to shuffle the trials
                    RspMatTest(:,cells) = psthNi{OleCellsi(cells)}{st}{im}(randi(NumCompleteTrials),:) ;
                end

                if sum(ismember(OleCellsi(cells),OleDsCellsi))>0 && KillDsRspFlag ; % if you want to kill DS responses
                    RspMatTest(:,cells) = RspMatTest(:,cells)*0 ;
                end
                
                if ShuffTestFlag ; % if you want to shuffle the test set (only shuffle within testing pnts and trial)
                    ShuffledPnts = circshift(TestingPnts_Ni,[0,randi(length(TestingPnts_Ni))]) ; % randperm(TestingPnts_Ni) ;
                    RspMatTest(TestingPnts_Ni,cells) = RspMatTest(ShuffledPnts,cells) ;
                end
                    
            end
            
            if BestCellTestFlag ; % if using only the strongest response at each time point
                for tp = 1:size(RspMatTest,1) ; % for each time point
                    [mx,mpnt] = max(RspMatTest(tp,:)) ; % find max response
                    RspMatTest(tp,:) = RspMatTest(tp,:)*0 ; % zero all responses
                    RspMatTest(tp,mpnt) = mx ; % keep max response only
                end
            end

            if OqeFlag ;
                [Temp1,Temp2] = OleTester(OleWeights_altTest,RspMatTest,'OQE') ;
                DirEstimate{st}{im}(t,:) = Temp1 ;
                DirEstimate_mag{st}{im}(t,:) = Temp2 ;
            else
                [Temp1,Temp2]= OleTester(OleWeights_altTest,RspMatTest) ;
                DirEstimate{st}{im}(t,:) = Temp1 ;
                DirEstimate_mag{st}{im}(t,:) = Temp2 ;
            end
            
            % average response 
            RspMatAverage{st}{im}(t,:) = mean(RspMatTest,2) ;
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
ErrorHistX = [0:3:180] ;

Temp = [] ;
for st = 1:NumDirNi ; % for each stimulus 
    for im = 1:length(TestingIm_Ni) ; % for each image to test
        Temp = [Temp,DirEstimate_Error{st}{TestingIm_Ni(im)}(TestingTrials_Ni,TestingPnts_Ni)] ; % error where tested
    end
end
ErrorAll = Temp(:) ;
ErrorHist = hist(ErrorAll,ErrorHistX) ;
ErrorMean = mean(ErrorAll) ;
ErrorStd = std(ErrorAll) ;
ErrorMed = median(ErrorAll) ; 

% response mag histogram
MagHistX = [0:.1:3] ;

Temp = [] ;
for st = 1:NumDirNi ; % for each stimulus 
    for im = 1:length(TestingIm_Ni) ; % for each image to test
        Temp = [Temp,DirEstimate_mag{st}{TestingIm_Ni(im)}(TestingTrials_Ni,TestingPnts_Ni)] ; % mag where tested
    end
end
MagAll = Temp(:) ;
MagHist = hist(MagAll,MagHistX) ;
MagMean = mean(MagAll) ;
MagStd = std(MagAll) ;

%% make pdf of estimate error vs estimate magnitude (or mean acitivity)

ErrorMagPdf = nan(length(ErrorHistX),length(MagHistX)) ;
ErrorBinDelta = (ErrorHistX(2)-ErrorHistX(1))/2 ; % assumes a uniform error X
MagBinDelta = (MagHistX(2)-MagHistX(1))/2 ; % assumes a uniform error X

for ErrorB = 1:length(ErrorHistX) ; % for each error bin
    for MagB = 1:length(MagHistX) ; % for each mag bin
        ErrorMagPdf(ErrorB,MagB) = sum( ...
            [ErrorAll>=(ErrorHistX(ErrorB)-ErrorBinDelta) & ErrorAll<(ErrorHistX(ErrorB)+ErrorBinDelta)].* ...
            [MagAll>=(MagHistX(MagB)-MagBinDelta) & MagAll<(MagHistX(MagB)+MagBinDelta)]) ; % number of observations
    end
end
ErrorMagPdf = ErrorMagPdf/sum(ErrorMagPdf(:)) ; % make into pdf from hist

% info in jpdf
ErrorMagPdf_MI = MutualInfoInJointPdf(ErrorMagPdf) ; % mutual information

% find median error as function of mag estimate
for MagB = 1:length(MagHistX) ; % for each mag bin
    ErrorMedian_eachBin(MagB) = median(ErrorAll([MagAll>=(MagHistX(MagB)-MagBinDelta)...
         & MagAll<(MagHistX(MagB)+MagBinDelta)])) ; % interpolate to median for each bin
    
    ErrorMedianCurve(MagB) = median(ErrorAll(MagAll>=MagHistX(MagB))) ; % median above threshold
end

% find median for top 10 % mag estimates
MagPercentThreshold = 0.9 ; %(fraction below cut out)
MagAllSort = sort(MagAll) ; % mag sorted
magThresh = MagAllSort(round(length(MagAllSort)*0.9)) ; % threshold
ErrorBestMedian = median(ErrorAll(MagAll>magThresh)) ;

ErrorMed_TEMPSET(TEMP_LOOP) = ErrorMed ;
ErrorBestMedian_TEMPSET(TEMP_LOOP) = ErrorBestMedian ;
end % TEMP LOOP
%% crossterm analysis for OQE
if OqeFlag ; 
    
    % direction of cross terms 
    for cells = 1:NumOleCells ; % for each cell in the OQE
        [temp,cti] = find(CrossTermCelli==cells) ; % find cross terms with that cell 
        cti = cti + NumOleCells ; % adjust for first set of terms which are not cross terms - index of OleWeights

        CrsTermVecSum(cells,:) = PolarVectorAddition(OlePolarWeights(cti,:)) ; % vector sum of cross terms
        
        TermCrsTermAngDiff(cells) = acuteAngle(CrsTermVecSum(cells,1),OlePolarWeights(cells,1)) ;
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

%% calculate correlation by distance/type 
CalcSynchronyFlag = false; % if you want fraction of time responding together instead of corr coef

NumAllDsCells = size(TypeCellMat,1) ; % number of all ds cells

% correlation mat for Natural Image and Dg
corrMatNi = nan(NumAllDsCells-1,NumAllDsCells) ;
for cell1 = 1:NumAllDsCells-1 ; % for each ds cell
    for cell2 = cell1+1:NumAllDsCells ; % for each ds cell that hasn't been looked at 
        corrMatNi(cell1,cell2) = 0 ; % get rid of nan 
        Counter = 0 ;
        for st = 1:NumDirNi ; % for each stimulus
            for im = 1:NumImages ; % for each image
                for t=1:NumCompleteTrials ; % for each trial
                    TempPsth1Ni = psthNi{DsCelli{TypeCellMat(cell1,1)}(TypeCellMat(cell1,2))}{st}{im}(t,SteadyStatePntsNi) ;
                    TempPsth2Ni = psthNi{DsCelli{TypeCellMat(cell2,1)}(TypeCellMat(cell2,2))}{st}{im}(t,SteadyStatePntsNi) ;

                    if sum(TempPsth1Ni)>0 && sum(TempPsth2Ni)>0 ; % if both cells respond on that trial
                        if ~CalcSynchronyFlag ;
                            corrMatNi(cell1,cell2) = corrMatNi(cell1,cell2)+...
                                corr(TempPsth1Ni(:),TempPsth2Ni(:)) ; % corr coef between cells
                        else % if you want fraction of synchronous response
                            corrMatNi(cell1,cell2) = corrMatNi(cell1,cell2)+ ...
                                sum((TempPsth1Ni(:).*TempPsth2Ni(:))>0)/sum((TempPsth1Ni(:)+TempPsth2Ni(:))>0) ;
                        end
                        Counter = Counter+1 ;
                    end
                end
            end
        end
        corrMatNi(cell1,cell2) = corrMatNi(cell1,cell2)/Counter ; % average from sum 
    end
end
  
% distance between cells and pair type

for cell1 = 1:NumAllDsCells-1 ; % for each ds cell
    for cell2 = cell1+1:NumAllDsCells ; % for each ds cell that hasn't been looked at 
        PairDist(cell1,cell2) = sqrt(sum((EiCtr(DsCelli{TypeCellMat(cell1,1)}(TypeCellMat(cell1,2)),:)-...
            EiCtr(DsCelli{TypeCellMat(cell2,1)}(TypeCellMat(cell2,2)),:)).^2)) ; % distance between pair
            
        DsType1 = TypeCellMat(cell1,1) ;
        DsType2 = TypeCellMat(cell2,1) ;
        
        if DsType1==DsType2 % if they are the same cell type
            Ctemp = 1 ;
        elseif abs(DsType1-DsType2)==1 | abs(DsType1-DsType2)==3 ; % if they are 90
            Ctemp = 2 ;
        else % if they are 180
            Ctemp = 3 ;
        end
        
        PairType(cell1,cell2) =  Ctemp ; % angle relationship between pair (0,90, or 180)
    end
end
      
% correlation averaged within PairTypes as function of distance
eiDistBin = 100 ;
eiDistMax = 1000 ;
DistBins = [.0001:eiDistBin:eiDistMax] ;
for pt = 1:3 ; % for each possible pair type
    pti = PairType==pt ; % indicies
    distTemp = PairDist(pti) ;
    corrNiTemp = corrMatNi(pti) ;
    for b=1:length(DistBins)-1 ;
        corrNiByDist{pt}(b) = mean(corrNiTemp(distTemp>=DistBins(b) & distTemp<DistBins(b+1))) ;
    end
end
    
%% DS response hisotgrams conditional on NonDS activity
ThreshSpike = 0 ; %
ThreshRad = 100 ; % (um) Distance for conditional

for DsType = 1:length(DsCelli) ; % for each DS cell type
    for cells = 1:length(DsCelli{DsType}) ; % for each DS cell
        
        DirDistSolo{DsType}{cells} = zeros(1,NumDirNi) ; % for each direction number of times cell fires alone 
        DirDistWithNonDs{DsType}{cells} = zeros(1,NumDirNi) ; % for each direction number of times cell is active at the same time as a local NonDS cell
        DirDistWithOtherDs{DsType}{cells} = zeros(1,NumDirNi) ; % for each direction number of times cell fires is active at the same time as another local DS cell 
        
        for cells2=1:length(dataRun.spikes) ; % for each cell
            EiDistTemp(cells2) = sqrt(sum((EiCtr(cells,:)-EiCtr(cells2,:)).^2)) ;
        end
    
        LocalCelliTemp = find(EiDistTemp<ThreshRad) ;
        LocalNonDsCellsiTemp = setdiff(LocalCelliTemp,cell2mat(DsCelli)) ; % Non DS cells
        LocalDsCellsiTemp = intersect(LocalCelliTemp,cell2mat(DsCelli)) ; %DS cells
    
        for st = 1:NumDirNi ; % for each stimulus
            for im = 1:NumImages ; % for each image
                for t=1:NumCompleteTrials ; % for each trial
                    si = find(psthNi{DsCelli{DsType}(cells)}{st}{im}(t,:)>ThreshSpike) ; % times that cell spiked above threshold
                    for spk = 1:length(si) ; % for each spike
                        NonDsSpkCheck = 0 ; % initiate
                        DsSpkCheck = 0 ;
                        
                        % check for simultaneous activity in NonDS cells
                        for cells2=1:length(LocalNonDsCellsiTemp) ; % for each NonDS cell within local area
                            if psthNi{LocalNonDsCellsiTemp(cells2)}{st}{im}(t,si(spk))>ThreshSpike ; % does it spike above threshold
                                NonDsSpkCheck = 1 ;
                            end
                        end
                        
                        if NonDsSpkCheck == 0 ;
                            % check for simultaneous activity in NonDS cells
                            for cells2=1:length(LocalDsCellsiTemp) ; % for each DS cell within local area
                                if LocalCelliTemp(cells2)~=cells && ...
                                        psthNi{LocalNonDsCellsiTemp(cells2)}{st}{im}(t,si(spk))>ThreshSpike ; % does it spike above threshold
                                    DsSpkCheck = 1 ;
                                end
                            end
                        end
                        
                        DirDistWithNonDs{DsType}{cells}(st) = DirDistWithNonDs{DsType}{cells}(st)+NonDsSpkCheck ; % add an observation
                        DirDistWithOtherDs{DsType}{cells}(st) = DirDistWithOtherDs{DsType}{cells}(st)+DsSpkCheck ; % add an observation           
                        DirDistSolo{DsType}{cells}(st) = DirDistSolo{DsType}{cells}(st)+((NonDsSpkCheck+DsSpkCheck)==0) ; % add an observation
                    end     
                end
            end
        end
    end
end

% combine distributions across all cells
r=1 ;
for DsType = 1:length(DsCelli) ; % for each DS cell type
    for cells = 1:length(DsCelli{DsType}) ; % for each DS cell
        histFolded_WithNonDs(r,:) = AverageHistAroundPeak(DirDistWithNonDs{DsType}{cells}) ;
        histFolded_WithDs(r,:) = AverageHistAroundPeak(DirDistWithOtherDs{DsType}{cells}) ;
        histFolded_Solo(r,:) = AverageHistAroundPeak(DirDistSolo{DsType}{cells}) ;
        
        histFolded_WithNonDs_pdf(r,:) = histFolded_WithNonDs(r,:)/sum(histFolded_WithNonDs(r,:));
        histFolded_WithDs_pdf(r,:) = histFolded_WithDs(r,:)/sum(histFolded_WithDs(r,:));
        histFolded_Solo_pdf(r,:) = histFolded_Solo(r,:)/sum(histFolded_Solo(r,:));
        r=r+1 ;
    end
end
histFolded_WithNonDs_pdf_mean = mean(histFolded_WithNonDs_pdf) ;
histFolded_WithDs_pdf_mean = mean(histFolded_WithDs_pdf) ;
histFolded_Solo_pdf_mean = mean(histFolded_Solo_pdf) ;

histFolded_WithNonDs_pdf_std = std(histFolded_WithNonDs_pdf) ;
histFolded_WithDs_pdf_std = std(histFolded_WithDs_pdf) ;
histFolded_Solo_pdf_std = std(histFolded_Solo_pdf) ;
    
%% calculate response histograms given DS and population response
RhistX = [0:1:psthSmoothPnts]/(psthSmoothPnts*psthBinTime) ; % possible spike rates
RhistXDelta = RhistX(2)-RhistX(1) ;
% calc P(R|stim)
Rhist_AllCells = nan(length(cell2mat(DsCelli)),NumDirNi,length(RhistX)) ; % prep empty mat
RallHist_AllCells = nan(length(cell2mat(DsCelli)),NumDirNi,length(RhistX)) ; % prep empty mat
RGivenHist_AllCells = nan(length(cell2mat(DsCelli)),length(RhistX),NumDirNi,length(RhistX)) ; % prep empty mat

centerPnt = 4 ; % where to center the tuning curves for across cell data
c = 0 ;
for DsType = 1:length(DsCelli) ; % for each DS cell type
    for cells = 1:length(DsCelli{DsType}) ; % for each DS cell
        c=c+1 ; % update counter
        
        Rhist = zeros(NumDirNi,length(RhistX)) ; % prep empty mat
        RallHist = zeros(NumDirNi,length(RhistX)) ; % prep empty mat
        RGivenHist = zeros(length(RhistX),NumDirNi,length(RhistX)) ; % prep empty mat
        
        % get covariance with all other cells
        for cells2 = 1:length(dataRun.spikes) ; % for each cell  
            TempA = cell2mat(psthNi_mean{DsCelli{DsType}(cells)}) ; % colapse across directions and images
            TempB = cell2mat(psthNi_mean{cells2}) ;
            Cv(cells2) = corr(TempA(:),TempB(:)) ; % corr coef
        end
        
        for st = 1:NumDirNi ; % for each direction
            for im = 1:NumImages ; % for each image
                for t=1:NumCompleteTrials ; % for each trial
                    R = psthNi{DsCelli{DsType}(cells)}{st}{im}(t,:) ; % response of cell 
                    Rhist(st,:) = Rhist(st,:)+ hist(R,RhistX) ; % hist of response
                    
                    Rall = zeros(1,length(psthNi{DsCelli{DsType}(cells)}{st}{im}(t,:))) ; % prep vector
                    for cells2 = 1:length(dataRun.spikes) ; % for each cell
                        if ~ismember(cells2,cell2mat(DsCelli)); % if its not a DS cell  
                            Rall = Rall + psthNi{cells2}{st}{im}(t,:)*Cv(cells2)/sum(Cv) ; % wieghted sum
                        end
                    end
                    Rall = Rall/max(Rall) ;
                    
                    RallHist(st,:) = RallHist(st,:) + hist(Rall,RhistX) ; % histogram
                    
                    for RallBin = 1:length(RhistX) ; % for each bin of Rall
                        Rselect = R(Rall<RhistX(RallBin)+RhistXDelta & Rall>=RhistX(RallBin)-RhistXDelta) ; % givens

                        RGivenHist(RallBin,st,:) = squeeze(RGivenHist(RallBin,st,:))' + hist(Rselect,RhistX) ;
                    end
                end
            end
        end
        [temp,tempi] = max(mean(Rhist(:,2:end),2)) ; % direction with max response
        shiftTemp = centerPnt - tempi ;
        
        Rhist_AllCells(c,:,:) = circshift(Rhist,[shiftTemp,0]) ; % 
        RallHist_AllCells(c,:,:) = circshift(RallHist,[shiftTemp,0]) ; % 
        RGivenHist_AllCells(c,:,:,:) = circshift(RGivenHist,[0,shiftTemp,0]) ;
    end
end
   
% average across cells
Rhist_AllCells_mean = squeeze(mean(Rhist_AllCells,1)) ; % 
RallHist_AllCells_mean = squeeze(mean(RallHist_AllCells,1)) ; % 
RGivenHist_AllCells_mean = squeeze(mean(RGivenHist_AllCells,1)) ; %       
        
%% estimate direction based on estimate magnitude threshold

EstimateMagThreshold = 1 ;

for st = 1:NumDirNi ; % for each stimulus
    for im = 1:NumImages ; % for each image to test
        for t=1:NumCompleteTrials ; % for each trial
            DirEstimate_HighMag{st}{im}(t,:) = zeros(1,length(DirEstimate{st}{im}(t,:))) ; % set default  
            tpi = find(DirEstimate_mag{st}{im}(t,:)>EstimateMagThreshold) ; % find time points above mag threshold
            for tp = 1:length(tpi) ; % for each threshold crossing
                DirEstimate_HighMag{st}{im}(t,tpi(tp):end) = DirEstimate{st}{im}(t,tpi(tp)) ; % set all future points to the estimate at that time
            end
            
            for tp=1:length(DirEstimate{st}{im}(1,:)) ; % for each time point
                DirEstimate_HighMag_Error{st}{im}(t,tp) = acuteAngle(DirEstimate_HighMag{st}{im}(t,tp),stimulus.directions(st)) ;
            end
        end
    end
end

% error histogram - error across all tested images, times and directions 
Temp = [] ;
for st = 1:NumDirNi ; % for each stimulus 
    for im = 1:length(TestingIm_Ni) ; % for each image to test
        Temp = [Temp,DirEstimate_HighMag_Error{st}{TestingIm_Ni(im)}(TestingTrials_Ni,TestingPnts_Ni)] ; % error in steady state
    end
end
HighMagErrorAll = Temp(:) ;
HighMagErrorHist = hist(HighMagErrorAll,ErrorHistX) ;

%% figures

figure % mean direction estimate - NI
for st = 1:NumDirNi ; % for each stimulus
    %for im = 1:NumImages ;
        plot(psthTimeNi,DirEstimate_Mean{st}(im,:),Color_list{st})
        hold on
        plot(psthTimeNi,ones(size(psthTimeNi))*stimulus.directions(st),[Color_list{st},'--'])
         %pause; hold off
    %end
end

figure % mean direction estimate - DG
for st = 1:NumDirDg ; % for each stimulus
    for tmp = 1:NumTempPeriods ;
        plot(psthTimeDg,DirEstimateDg_Mean{st}(tmp,:),Color_list{st})
        hold on
        plot(psthTimeDg,ones(size(psthTimeDg))*dataRunDg.stimulus.params.DIRECTION(st),[Color_list{st},'--'])
        pause; hold off
    end
end

figure % direction estimates on each trial - NI
for st = 1:NumDirNi ; % for each stimulus
    for im=1:NumImages ; % for each image
        for t=1:NumCompleteTrials ;
            %subplot(3,1,1)
            plot(psthTimeNi,DirEstimate{st}{im}(t,:),'r')
            plot(psthTimeNi,DirEstimate_mag{st}{im}(t,:),'r')
            %plot(psthTimeNi,DirEstimate_HighMag{st}{im}(t,:))
            hold on
            plot(psthTimeNi,ones(size(psthTimeNi))*stimulus.directions(st),'r--')
            hold off

            %subplot(3,1,2)
            %plot(psthTimeNi,DirEstimate_Error{st}{im}(t,:))
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
subplot(3,1,1) % error histogram
plot(ErrorHistX,cumsum(ErrorHist)/sum(ErrorHist),'r')
hold on
%plot(ErrorHistX,cumsum(HighMagErrorHist)/sum(HighMagErrorHist),'r*')
xlabel('error (deg)')
ylabel('fraction')

subplot(3,1,2) % distribution of confidence
plot(MagHistX,cumsum(MagHist)/sum(MagHist),'r')
hold on
xlabel('Mag (au)')
ylabel('fraction')

subplot(3,1,3) % error median as function of confidence
plot(MagHistX,ErrorMedianCurve,'k')
hold on
plot(MagHistX,ErrorMedian_eachBin,'r:')
xlabel('estimate mag')
ylabel('median error (deg)')

figure % Ni confidence error 
[ax,h1,h2] = plotyy(MagHistX,ErrorMedianCurve,MagHistX,cumsum(MagHist)/sum(MagHist)) ;
hold on
set(ax(1),'YLim',[0 90])
set(ax(2),'YLim',[0 1])
set(ax(1),'XLim',[0 3])
set(ax(2),'XLim',[0 3])
set(ax(1),'YTick',[0:15:90])
set(ax(2),'YTick',[0:.1:1])

figure % time dependant error -NI
for st = 1:NumDirNi ;
    for im = 1:NumImages ; 
        errorbar(psthTimeNi,DirEstimate_Error_Mean{st}(im,:),DirEstimate_Error_SEM{st}(im,:),'g')
        %pause
    end
end
xlabel('time (sec)')
ylabel('error (deg)')


figure ; % error box plot
boxplot(ErrorAll)

figure % error - DG

subplot(2,1,1) % time dependant error 
for st = 1:NumDirNi ;
    for tmp = 1:NumTempPeriods ; 
        errorbar(psthTimeDg,DirEstimateDg_Error_Mean{st}(tmp,:),DirEstimateDg_Error_SEM{st}(tmp,:),'r')
        pause
    end
end
xlabel('time (sec)')
ylabel('error (deg)')

subplot(2,1,2) % error histogram
plot(ErrorHistX,ErrorHistDg/sum(ErrorHistDg),'r')
hold on
plot(ErrorHistX,MinErrorHistDg/sum(MinErrorHistDg),'r--')
xlabel('error (deg)')
ylabel('fraction')

figure % wieghts
subplot(3,1,1) % weights for each cell
plot([1:NumOleCells],OlePolarWeights(1:NumOleCells,2),'bo')
hold on
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
subplot(1,3,1) % error vs. estimate magnitude each point
for t=1:NumCompleteTrials ; % for each trial
    for im = 1:NumImages ;
        for st = 1:NumDirNi ; % for each stimulus
            plot(DirEstimate_Error{st}{im}(t,:),DirEstimate_mag{st}{im}(t,:),'.')
            hold on
            % pause
        end
    end
end

subplot(1,3,2) % pdf
imagesc(flipud(ErrorMagPdf'/sum(ErrorMagPdf(:))))
colorbar

subplot(1,3,3) % slices of pdf by magnitude
for a=1:size(ErrorMagPdf,2) ; % for each magnitude
    plot(DirEstimate_Error_bins(1:end-1),cumsum(ErrorMagPdf(:,a)/sum(ErrorMagPdf(:,a))))
    hold on
    %pause
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


figure % plot EI positions for quads
for ri = 1:length(RealQuad_i) ; % for each quad
    for cells = 1:4 ; % for each cell in the quad
        ct = quad{RealQuad_i(ri)}(cells,1) ; % cell type
        cn = quad{RealQuad_i(ri)}(cells,2) ; % cell number
        plot(EiCtr(DsCelli{ct}(cn),1),EiCtr(DsCelli{ct}(cn),2),[Color_list{ct},'+'])
        drawCircle(EiCtr(DsCelli{ct}(cn),1),EiCtr(DsCelli{ct}(cn),2),150,'color',Color_list{ct})
        hold on 
    end
    pause
    hold off
end

figure % Dg tuning curves of cells used in Ole estimate
tmp = 1 ; % temporal period
for cells = 1:length(OleDsCellsi) ; % for each DS cell
    %plot(dataRunDg.stimulus.params.DIRECTION,DgTuningCurve{1}(OleDsCellsi(cells),:))
    polar(VectSumDg{tmp}(OleDsCellsi(cells),1)*pi/180,1,'*')
    hold on
end

figure % All Dg tuning curves of DS cells 
tmp = 1 ; % temporal period
for DsType = 1:4 ;
    for cells = 1:length(DsCelli{DsType}) ; % for each DS cell
        plot(dataRunDg.stimulus.params.DIRECTION,DgTuningCurve{1}(DsCelli{DsType}(cells),:),'color',Color_list{DsType})
        %polar(VectSumDg{tmp}(DsCelli{DsType}(cells),1)*pi/180,1,'*')
        hold on
    end
end

figure ; % Ni tuning curves
for DsType = 1:4 ;
    for cells = 1:length(DsCelli{DsType}) ; % for each DS cell
        plot(stimulus.directions,NiTuningCurve(DsCelli{DsType}(cells),:),'color',Color_list{DsType})
        hold on
    end
end

figure % plot raster of all cells for each trial 
t=1 ;
Temp = 1 ;
for DsType = 1:4 ; % for each ds
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        
        spk = dataRun.spikes{DsCelli{DsType}(cells)}-triggs(t) ;
        spk = spk(spk>=0 & spk<=StimDuration) ;
        plot(spk,ones(1,length(spk))*Temp,'.','color',Color_list{DsType})
        axis([0,4,0,length(dataRun.spikes)]) 
        hold on
        Temp = Temp + 1 ;
    end
end
NonDsi = setdiff([1:NumCells],cell2mat(DsCelli)) ;
for cells=1:NumCells ;
    if ismember(cells,NonDsi) ;
        spk = dataRun.spikes{cells}-triggs(t) ;
        spk = spk(spk>=0 & spk<=StimDuration) ;
        if ~isempty(spk) ;
            plot(spk,ones(1,length(spk))*Temp,'k.')
            axis([0,4,0,length(dataRun.spikes)]) 
        end
        Temp = Temp + 1 ;
    end
end

figure % plot raster of example cell for each trial 
cells =  DsCelli{1}(1) ;
for st = 1:NumDirNi ; % for each stimulus direction
    sti = find(dirShownFullTranspose(:)==stimulus.directions(st)) ; % index of triggers for that direction
    for im = 1:NumImages ; % for each image
        ti = sti([im:NumImages:length(sti)]) ; % trigger index for that image and direction 

        Temp = 1 ;
        subplot(2,1,1)
        for t=1:(NumCompleteTrials) ; % for each trial
            spk = dataRun.spikes{cells}-triggs(ti(t)) ;
            spk = spk(spk>=0 & spk<=StimDuration) ;

            for pnt=1:length(spk) ;
                plot([spk(pnt),spk(pnt)],[Temp-.5,Temp+.5],'k-')
                axis([0,4,-.5,NumCompleteTrials+.5]) 
                hold on
            end
            Temp = Temp + 1 ;
        end
        hold off
        subplot(2,1,2)
        plot(psthTimeNi, psthNi_mean{cells}{st}(im,:)) ; % average across trials
        xlim([0,4])
        pause
    end
end



figure % Dg raster
t=1 ; % trial
Temp = 1 ;
for DsType = 1:4 ; % for each ds
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        spk = dataRun.spikes{DsCelli{DsType}(cells)}-dataRunDg.stimulus.triggers(t) ;
        spk = spk(spk>=0 & spk<=TrialTrigInterval) ;
        plot(spk,ones(1,length(spk))*Temp,'.','color',Color_list{DsType})
        axis([0,4,0,length(dataRun.spikes)]) 
        hold on
        Temp = Temp + 1 ;
    end
end
NonDsi = setdiff([1:NumCells],cell2mat(DsCelli)) ;
for cells=1:NumCells ;
    if ismember(cells,NonDsi) ;
        spk = dataRun.spikes{cells}-dataRunDg.stimulus.triggers(t) ;
        spk = spk(spk>=0 & spk<=TrialTrigInterval) ;
        if ~isempty(spk) ;
            plot(spk,ones(1,length(spk))*Temp,'k.')
            axis([0,4,0,length(dataRun.spikes)]) 
        end
        Temp = Temp + 1 ;
    end
end
        
figure % plot raster for DG and NI for a example quad
EgQuadi = 2 ; % index of quad example
t = 3 ; % trial
tmp = 1 ; % DG temporal period

subplot(2,1,1) % dg
for st = 1:NumDirDg ; % for each direction 
    prm = [dataRunDg.stimulus.params.DIRECTION(st),dataRunDg.stimulus.params.TEMPORAL_PERIOD(tmp)] ; % params set
    ti = find((DgParams(:,1)==prm(1)).*(DgParams(:,2)==prm(2))==1) ; % index of triggers for that stim
    for cells=1:4 ; % for each cell
        ct = quad{RealQuad_i(EgQuadi)}(cells,1) ; % cell type
        cn = quad{RealQuad_i(EgQuadi)}(cells,2) ; % cell number
        c = DsCelli{ct}(cn) ; % cell i

        spk = dataRun.spikes{c}-dataRunDg.stimulus.triggers(ti(t)) ;
        spk = spk(spk>=0 & spk<=TrialTrigInterval) ;

        for pnt=1:length(spk) ;
            plot([spk(pnt),spk(pnt)],[cells-.5,cells+.5],'k-')
            axis([0,4,-.5,4.5]) 
            hold on
        end
    end

    pause
    hold off
end
        
subplot(2,1,2) % Ni
for st = 1:NumDirNi ; % for each stimulus direction
    sti = find(dirShownFullTranspose(:)==stimulus.directions(st)) ; % index of triggers for that direction
    for im = 1:NumImages ; % for each image
        for cells=1:4 ; % for each cell in the quad
            ct = quad{RealQuad_i(EgQuadi)}(cells,1) ; % cell type
            cn = quad{RealQuad_i(EgQuadi)}(cells,2) ; % cell number
            c = DsCelli{ct}(cn) ;
            
            ti = sti([im:NumImages:length(sti)]) ; % trigger index for that image and direction 
            spk = dataRun.spikes{c}-triggs(ti(t)) ;
            spk = spk(spk>=0 & spk<=StimDuration) ;

            for pnt=1:length(spk) ;
                plot([spk(pnt),spk(pnt)],[cells-.5,cells+.5],'k-')
                axis([0,4,-.5,4.5]) 
                hold on
            end
        end
        pause
        hold off
    end
end
        

    

figure % conditional tuning curves
errorbar([0:45:180],histFolded_WithNonDs_pdf_mean, histFolded_WithNonDs_pdf_std,'k')
hold on
errorbar([0:45:180],histFolded_Solo_pdf_mean, histFolded_Solo_pdf_std,'r')
xlabel('angle off prefered (deg)')
ylabel('fraction responses')
legend('spike with NonDS','spike alone')

figure % correlation by distance and pair type - NI (Fig 1E)
plot(PairDist,corrMatNi,'k.')
hold on
plot(PairDist(PairType==1),corrMatNi(PairType==1),'k.','MarkerSize',10)
plot(PairDist(PairType==2),corrMatNi(PairType==2),'r.','MarkerSize',10)
plot(PairDist(PairType==3),corrMatNi(PairType==3),'c.','MarkerSize',10)
plot(DistBins(2:end)-eiDistBin/2,corrNiByDist{1},'k-','LineWidth',2)
plot(DistBins(2:end)-eiDistBin/2,corrNiByDist{2},'r-','LineWidth',2)
plot(DistBins(2:end)-eiDistBin/2,corrNiByDist{3},'c-','LineWidth',2)
xlabel('distance (um)')
ylabel('corr coef')
hgsave('CorrByDist200msFig')
print(gcf, '-dpdf','CorrByDist200msFig')

figure % correlation by distance and pair type - DG
plot(PairDist,corrMatDg,'k*')
hold on
plot(PairDist(PairType==1),corrMatDg(PairType==1),'k*')
plot(PairDist(PairType==2),corrMatDg(PairType==2),'r*')
plot(PairDist(PairType==3),corrMatDg(PairType==3),'b*')
plot(DistBins(2:end),corrDgByDist{1},'k-')
plot(DistBins(2:end),corrDgByDist{2},'r-')
plot(DistBins(2:end),corrDgByDist{3},'b-')

figure % firing rate histograms - NI only 
plot(psthNiHistX,psthNiHist)
hold on
xlabel('firing rate (hz)')
ylabel('probability')

figure % number of cells firing in a quad (Fig 1D)
plot(quadNumCellsHistX,quadNiNumCellsHist,'k') % NI
hold on
plot(quadNumCellsHistX,quadDgNumCellsHist) % DG
xlabel('number active cells in quad')
ylabel('probability')

figure % fraction of strongest cell to all others firing in a quad 
plot(quadRateRatioHistX,quadNiRateRatioHist) % NI
hold on
%plot(RateRatioX,quadDgRateRatio_dist) % DG
xlabel('most active / all cells')
ylabel('probability')

figure % conditional histograms
subplot(3,1,1) % DS cells 
for a=1:6 ; % for each spike rate bin (not the last for noise)
    plot(Rhist_AllCells_mean(:,a)/sum(Rhist_AllCells_mean(:,a)),'LineWidth',a) ;
    hold on
end
xlabel('direction')
ylabel('probability')

subplot(3,1,2) % DS cells given DS low 
for a=1:length(RhistX)-5 ; % for each spike rate bin (not the last for noise)
    plot(RGivenHist_AllCells_mean(a,:,1)/sum(RGivenHist_AllCells_mean(a,:,1)),'LineWidth',a) ;
    hold on
end
xlabel('direction')
ylabel('probability')

subplot(3,1,3) % DS cells given DS high 
for a=1:length(RhistX)-5 ; % for each spike rate bin (not the last for noise)
    plot(RGivenHist_AllCells_mean(a,:,4)/sum(RGivenHist_AllCells_mean(a,:,4)),'LineWidth',a) ;
    hold on
end
xlabel('direction')
ylabel('probability')




