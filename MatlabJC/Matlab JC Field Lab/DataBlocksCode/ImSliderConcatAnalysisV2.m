function ForIgor = ImSliderConcatAnalysisV2(DataBlock, DB, Params)


% this function is adapted from 'ImSliderConcatAnalysis' to deal with
% multiple images

% JC 6/22/2017 

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
psthSmoothTime = .02 ; % size of psth Bin (smooth sliding window)
EstimateFromTrialFlag = true ;

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

% Make all DS cells have the same prefered direction - OPTIONAL
% for DsType = 1:length(DsTypeName) ; % for each ds type
%     Temp = nans(length(DsCelli{DsType}),2) ;
%     for cells=1:length(DsCelli{DsType}) ; % for each cell
%         Temp(cells,:) = VectSumDg(DsCelli{DsType}(cells),:) ;
%     end
%     VectSumDg_byDsType(DsType,:) = PolarVectorAddition(Temp) ;
%     
%     for cells=1:length(DsCelli{DsType}) ; % for each cell
%         VectSumDg(DsCelli{DsType}(cells),:) = VectSumDg_byDsType(DsType,:) ; % replace original with average
%     end
% end

%% select stimulus estimate input

if EstimateFromTrialFlag ;
    TestTrial = 2 ;
    for cells=1:length(dataRun.spikes) ; % for each cell
        for st = 1:NumDir ; % for each stimulus direction
            psth_EstIn{cells}(st,:) = psth{cells}{st}(TestTrial,:) ;
        end
       
        for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus (speed and direction)
            psthDg_EstIn{cells}(st,:) = psthDg{cells}{st}(TestTrial,:) ;
        end
    end
else
    psth_EstIn = psth_mean ;
    psthDg_EstIn = psthDg_mean ;
end


%% stimulus estimate - grating prefered vector sum on *natural images*

% correct VectSumDg for 180 degree difference between dg and image
for cells=1:length(dataRun.spikes) ; % for each cell 
    if VectSumDg(cells,1)+180<360 ;
        VectSumDg_alt(cells,1) = VectSumDg(cells,1)+180 ;
    else
        VectSumDg_alt(cells,1) = VectSumDg(cells,1)-180 ;
    end
end
        
% estimate stim direction
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        clear CellCnt VectTemp
        
        CellCnt = 1 ; % ds cell count
        for DsType = 1:length(DsTypeName) ; % for each ds type
            for cells=1:length(DsCelli{DsType}) ; % for each cell
            %for cells=1:min(DsCellTypeNum) ; % for each cell (max number in a particular direction so all directions are equally represented)
                VectTemp(CellCnt,:) = [VectSumDg_alt(DsCelli{DsType}(cells),1),...
                    psth_EstIn{DsCelli{DsType}(cells)}(st,t)] ;
                CellCnt = CellCnt+1 ;
            end
        end
        
        VectSumImTemp = PolarVectorAddition(VectTemp) ;
        DirEstimate(st,t) = VectSumImTemp(1) ;
        DirEstimate_vecLength(st,t) = VectSumImTemp(2) ;
    end
end

% error 
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        AcuteAngleError(st,t) = acuteAngle(DirEstimate(st,t),stimulus.directions(st)) ;
    end
end

%% stimulus estimate - grating prefered vector sum on *ds typing drifting grating*

% estimate stim direction
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus (speed and direction)
    for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
        clear CellCnt VectTemp
        
        CellCnt = 1 ; % ds cell count
        for DsType = 1:length(DsTypeName) ; % for each ds type
            %for cells=1:length(DsCelli{DsType}) ; % for each cell
            for cells=1:min(DsCellTypeNum) ; % for each cell (max number in a particular direction so all directions are equally represented)
                VectTemp(CellCnt,:) = [VectSumDg(DsCelli{DsType}(cells),1),...
                    psthDg_EstIn{DsCelli{DsType}(cells)}(st,t)] ;
                CellCnt = CellCnt+1 ;
            end
        end
        
        VectSumImTemp = PolarVectorAddition(VectTemp) ;
        DirEstimateDg(st,t) = VectSumImTemp(1) ;
        DirEstimateDg_vecLength(st,t) = VectSumImTemp(2) ;
    end
end

% error 
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus 
    for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
        AcuteAngleErrorDg(st,t) = acuteAngle(DirEstimateDg(st,t),dataRunDg.stimulus.combinations(st).DIRECTION) ;
    end
end

%% local population vector estimate - images

% correct VectSumDg for 180 degree difference between dg and image
for cells=1:length(dataRun.spikes) ; % for each cell 
    if VectSumDg(cells,1)+180<360 ;
        VectSumDg_alt(cells,1) = VectSumDg(cells,1)+180 ;
    else
        VectSumDg_alt(cells,1) = VectSumDg(cells,1)-180 ;
    end
end

% find distance between all cells
for cells=1:length(dataRun.spikes) ; % for each cell
    for cells2=1:length(dataRun.spikes) ; % for each other cell
        EiDist(cells,cells2) = sqrt(sum((EiCtr(cells,:)-EiCtr(cells2,:)).^2)) ;
    end
end

% find quadruplets of DS cells
QuadSet = [] ; % prep quads
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(DsCelli{DsType}) ; % for each cell
        TempQuad = DsCelli{DsType}(cells) ; 
        for DsType2 = 1:length(DsTypeName) ; % for each ds type
            if DsType2~=DsType ; % if its not the same type
                [mdist,mi] = min(EiDist(DsCelli{DsType}(cells),DsCelli{DsType2})) ; % nearest cell
                if mdist<minEiRad ; % if the disance is small enough
                    TempQuad = [TempQuad,DsCelli{DsType2}(mi)] ; 
                end
            
                if length(TempQuad)==4 ;
                    QuadSet = [QuadSet;TempQuad] ; 
                end
            end
        end
    end
end
QuadSetUnique = unique(sort(QuadSet,2),'rows') ;        
       
% calculate direction estimate for each quadruplet
for q = 1:size(QuadSetUnique,1) ; % for each quad
    for st = 1:NumDir ; % for each stimulus direction
        for t=1:length(psth_mean{1}) ; % for each time point in the psth
            clear VectTemp
            for DsType=1:length(DsTypeName) ; % for each ds cell in the quad
                VectTemp(DsType,:) = [VectSumDg_alt(QuadSetUnique(q,DsType),1),...
                        psth_EstIn{QuadSetUnique(q,DsType)}(st,t)] ; % quad (prefered direction,spike rate)
            end
            
            DsQuadSum{q}(st,t) = sum(VectTemp(:,2)) ; % sum of firing rate (not vector sum)
            DsQuadSpreadIndex{q}(st,t) = (sum(VectTemp(:,2))-max(VectTemp(:,2)))/(sum(VectTemp(:,2))+minDenominator) ; % fraction of rate that not max firing cell
            DsQuadSpreadIndex2{q}(st,t) = sum(VectTemp(:,2)>max(VectTemp(:,2))/20) ; % number of cells more than 5% the firing rate of the max
            
            VectSumTemp = PolarVectorAddition(VectTemp) ; % vector sum for quad
            VectSumImQuad_dirEstimate{q}(st,t) = VectSumTemp(1) ; % direction of vector sum
            VectSumImQuad_vecLength{q}(st,t) = VectSumTemp(2) ; % length of vector sum
        end
    end
end

% vector sum of quads
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        VectTemp = nans(size(QuadSetUnique,1),2) ; % prep mat
        for q = 1:size(QuadSetUnique,1) ; % for each quad
            VectTemp(q,1) = VectSumImQuad_dirEstimate{q}(st,t) ; % direction of quad vector
            VectTemp(q,2) = VectSumImQuad_vecLength{q}(st,t) ; % length of quad vector
        end
        Temp = PolarVectorAddition(VectTemp) ;
        VectSumImQuad_dirEstimate_QuadSum(st,t) = Temp(1) ;
    end
end

% average across unit normalized quads
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        VectTemp = ones(size(QuadSetUnique,1),2) ; % set quad vector length to 1
        for q = 1:size(QuadSetUnique,1) ; % for each quad
            VectTemp(q,1) = VectSumImQuad_dirEstimate{q}(st,t) ; % direction of quad vector
        end
        Temp = PolarVectorAddition(VectTemp) ;
        VectSumImQuad_dirEstimate_QuadUnitAv(st,t) = Temp(1) ;
    end
end
        
% error of each quad
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        for q = 1:size(QuadSetUnique,1) ; % for each quad    
            AcuteAngleErrorQuads{q}(st,t) = acuteAngle(VectSumImQuad_dirEstimate{q}(st,t),stimulus.directions(st)) ;
        end
    end
end

% error of most accurate quads
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        Temp = nans(size(QuadSetUnique,1),1) ; % prep mat
        for q = 1:size(QuadSetUnique,1) ; % for each quad    
            Temp(q) = acuteAngle(VectSumImQuad_dirEstimate{q}(st,t),stimulus.directions(st)) ; % error of that quad
        end
        [AcuteAngleErrorQuadBest(st,t),BestQuadi(st,t)] = min(Temp) ; % error and index of most accurate quad
    end
end

% error of vector sum
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        AcuteAngleErrorQuadSum(st,t) = acuteAngle(VectSumImQuad_dirEstimate_QuadSum(st,t),stimulus.directions(st)) ;
    end
end

% error of unit vector sum
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        AcuteAngleErrorQuadUnitAv(st,t) = acuteAngle(VectSumImQuad_dirEstimate_QuadUnitAv(st,t),stimulus.directions(st)) ;
    end
end

% error of quad with longest vector sum (winner take all)
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        clear Temp
        for q = 1:size(QuadSetUnique,1) ; % for each quad    
            Temp(q) = VectSumImQuad_vecLength{q}(st,t) ;
        end
        AcuteAngleErrorQuadLongest(st,t) = acuteAngle(VectSumImQuad_dirEstimate{find(Temp==max(Temp),1,'first')}(st,t),stimulus.directions(st)) ; 
    end
end

% Estimate of sum of quad estimates wieghted by ds set firing rate
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        clear VectTemp
        for q = 1:size(QuadSetUnique,1) ; % for each quad 
            VectTemp(q,1) = VectSumImQuad_dirEstimate{q}(st,t); % sum of DS cell activity
            VectTemp(q,2) = VectSumImQuad_vecLength{q}(st,t)*DsQuadSum{q}(st,t) ; % sum of DS cell activity
        end
        Temp = PolarVectorAddition(VectTemp) ;
        AcuteAngleErrorQuadStrongest(st,t) = acuteAngle(Temp(1),stimulus.directions(st)) ; % weighted sum
        %AcuteAngleErrorQuadStrongest(st,t) = acuteAngle(VectSumImQuad_dirEstimate{find(VectTemp(:,2)==max(VectTemp(:,2)),1,'first')}(st,t),stimulus.directions(st)) ; % winner take all 
    end
end

% Estimate of sum of quad estimates wieghted by populations with best spread
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        clear VectTemp
        for q = 1:size(QuadSetUnique,1) ; % for each quad 
            VectTemp(q,1) = VectSumImQuad_dirEstimate{q}(st,t); % sum of DS cell activity
            VectTemp(q,2) = VectSumImQuad_vecLength{q}(st,t)*DsQuadSpreadIndex{q}(st,t) ;
        end
        Temp = PolarVectorAddition(VectTemp) ;
        AcuteAngleErrorQuadMostSpread(st,t) = acuteAngle(Temp(1),stimulus.directions(st)) ; 
        %AcuteAngleErrorQuadMostSpread(st,t) = acuteAngle(VectSumImQuad_dirEstimate{find(Temp==max(Temp),1,'first')}(st,t),stimulus.directions(st)) ; 
    end
end

%% local population vector estimate - Drifting grating

% calculate direction estimate for each quadruplet
for q = 1:size(QuadSetUnique,1) ; % for each quad
    for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus (speed and direction)
        for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
            clear VectTemp
            for DsType=1:length(DsTypeName) ; % for each ds cell in the quad
                VectTemp(DsType,:) = [VectSumDg(QuadSetUnique(q,DsType),1),...
                        psthDg_EstIn{QuadSetUnique(q,DsType)}(st,t)] ; % quad (prefered direction,spike rate)
            end
            
            DsQuadSumDg{q}(st,t) = sum(VectTemp(:,2)) ; % sum of firing rate (not vector sum)
            %DsQuadSpreadIndexDg{q}(st,t) = (sum(VectTemp(:,2))-max(VectTemp(:,2)))/(sum(VectTemp(:,2))+minDenominator) ; % fraction of rate that not max firing cell
            DsQuadSpreadIndexDg{q}(st,t) = sum(VectTemp(:,2)>max(VectTemp(:,2))/20) ; % number of cells more than 5% the firing rate of the max
            
            VectSumTemp = PolarVectorAddition(VectTemp) ; % vector sum for quad
            VectSumDgQuad_dirEstimate{q}(st,t) = VectSumTemp(1) ; % direction of vector sum
            VectSumDgQuad_vecLength{q}(st,t) = VectSumTemp(2) ; % length of vector sum
        end
    end
end

% vector sum of quads
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus (speed and direction)
    for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
        VectTemp = nans(size(QuadSetUnique,1),2) ; % prep mat
        for q = 1:size(QuadSetUnique,1) ; % for each quad
            VectTemp(q,1) = VectSumDgQuad_dirEstimate{q}(st,t) ; % direction of quad vector
            VectTemp(q,2) = VectSumDgQuad_vecLength{q}(st,t) ; % length of quad vector
        end
        Temp = PolarVectorAddition(VectTemp) ;
        VectSumDgQuad_dirEstimate_QuadSum(st,t) = Temp(1) ;
    end
end

% average across unit normalized quads
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus (speed and direction)
    for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
        VectTemp = ones(size(QuadSetUnique,1),2) ; % set quad vector length to 1
        for q = 1:size(QuadSetUnique,1) ; % for each quad
            VectTemp(q,1) = VectSumDgQuad_dirEstimate{q}(st,t) ; % direction of quad vector
        end
        Temp = PolarVectorAddition(VectTemp) ;
        VectSumDgQuad_dirEstimate_QuadUnitAv(st,t) = Temp(1) ;
    end
end
        
% error of each quad
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus (speed and direction)
    for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
        for q = 1:size(QuadSetUnique,1) ; % for each quad    
            AcuteAngleErrorQuadsDg{q}(st,t) = acuteAngle(VectSumDgQuad_dirEstimate{q}(st,t),dataRunDg.stimulus.combinations(st).DIRECTION) ;
        end
    end
end

% error of most accurate quads
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus (speed and direction)
    for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
        Temp = nans(size(QuadSetUnique,1),1) ; % prep mat
        for q = 1:size(QuadSetUnique,1) ; % for each quad    
            Temp(q) = acuteAngle(VectSumDgQuad_dirEstimate{q}(st,t),dataRunDg.stimulus.combinations(st).DIRECTION) ; % error of that quad
        end
        [AcuteAngleErrorQuadBestDg(st,t),BestQuadiDg(st,t)] = min(Temp) ; % error and index of most accurate quad
    end
end

% error of vector sum
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus (speed and direction)
    for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
        AcuteAngleErrorQuadSumDg(st,t) = acuteAngle(VectSumDgQuad_dirEstimate_QuadSum(st,t),dataRunDg.stimulus.combinations(st).DIRECTION) ;
    end
end

% error of unit vector sum
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus (speed and direction)
    for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
        AcuteAngleErrorQuadUnitAvDg(st,t) = acuteAngle(VectSumDgQuad_dirEstimate_QuadUnitAv(st,t),dataRunDg.stimulus.combinations(st).DIRECTION) ;
    end
end

% error of quad with longest vector sum (winner take all)
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus (speed and direction)
    for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
        clear Temp
        for q = 1:size(QuadSetUnique,1) ; % for each quad    
            Temp(q) = VectSumDgQuad_vecLength{q}(st,t) ;
        end
        AcuteAngleErrorQuadLongestDg(st,t) = acuteAngle(VectSumDgQuad_dirEstimate{find(Temp==max(Temp),1,'first')}(st,t),dataRunDg.stimulus.combinations(st).DIRECTION) ; 
    end
end

% Quads estimate from wieghted sum by quad firing rate
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus (speed and direction)
    for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
        clear VectTemp
        for q = 1:size(QuadSetUnique,1) ; % for each quad 
            VectTemp(q,1) = VectSumDgQuad_dirEstimate{q}(st,t); % sum of DS cell activity
            VectTemp(q,2) = DsQuadSumDg{q}(st,t) ; % sum of DS cell activity
        end
        Temp = PolarVectorAddition(VectTemp) ;
        AcuteAngleErrorQuadStrongestDg(st,t) = acuteAngle(Temp(1),dataRunDg.stimulus.combinations(st).DIRECTION) ; % weighted sum
        %AcuteAngleErrorQuadStrongestDg(st,t) = acuteAngle(VectSumDgQuad_dirEstimate{find(VectTemp(:,2)==max(VectTemp(:,2)),1,'first')}(st,t),dataRunDg.stimulus.combinations(st).DIRECTION) ; % winner take all 
    end
end

% error of quad with best firing rate spread (quads least dominated by 1 ds cells firing rate)
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus (speed and direction)
    for t=1:length(psthDg_mean{1}) ; % for each time point in the psth
        clear Temp
        for q = 1:size(QuadSetUnique,1) ; % for each quad    
            Temp(q) = DsQuadSpreadIndexDg{q}(st,t) ;
        end
        AcuteAngleErrorQuadMostSpreadDg(st,t) = acuteAngle(VectSumDgQuad_dirEstimate{find(Temp==max(Temp),1,'first')}(st,t),dataRunDg.stimulus.combinations(st).DIRECTION) ; 
    end
end

%% local vector estimates normalized by local contrast estimates

% find distance between quad center and all cells
for q = 1:size(QuadSetUnique,1) ; % for each quad
    QuadSetAll{q} = [] ; % prep vect
    QuadCenter = mean(EiCtr(QuadSetUnique(q,:),:),1) ;
    
    for cells=1:length(dataRun.spikes) ; % for each cell
        QuadEiDist = sqrt(sum((QuadCenter-EiCtr(cells,:)).^2)) ; % distance between quad center and other cell

        if QuadEiDist<minEiRad ; % if the cell is close enough
            QuadSetAll{q} = [QuadSetAll{q},cells] ; % add index of cell
        end
    end
end
    
% firing rate of all cells near the quads
for q = 1:size(QuadSetUnique,1) ; % for each quad
    Qpsth{q} = zeros(size(psth_EstIn{1})) ; % prep mat
    for cells=1:length(QuadSetAll{q}) ; % for all cells in the within distance of the DS set
        Qpsth{q} = Qpsth{q}+psth_EstIn{QuadSetAll{q}(cells)}/length(QuadSetAll{q}) ; % average across all cell rates
    end
end

% divide vector length by average firing rate of nearby cells
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        clear VectTemp
        for q = 1:size(QuadSetUnique,1) ; % for each quad
            VectTemp(q,1) = VectSumImQuad_dirEstimate{q}(st,t) ;
            VectTemp(q,2) = VectSumImQuad_vecLength{q}(st,t)/(Qpsth{q}(st,t)+minDenominator) ;
        end

        Temp = PolarVectorAddition(VectTemp) ;
        VectSumImQuad_dirEstimate_QuadAvNorm(st,t) = Temp(1) ;
    end
end

% mutliple normalized vector length by average firing rate of nearby cells
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        clear VectTemp
        for q = 1:size(QuadSetUnique,1) ; % for each quad
            VectTemp(q,1) = VectSumImQuad_dirEstimate{q}(st,t) ;
            VectTemp(q,2) = VectSumImQuad_vecLength{q}(st,t)*Qpsth{q}(st,t) ;
        end

        Temp = PolarVectorAddition(VectTemp) ;
        VectSumImQuad_dirEstimate_QuadAvNorm2(st,t) = Temp(1) ;
    end
end

% error quad av
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        AcuteAngleErrorQuadAvNorm(st,t) = acuteAngle(VectSumImQuad_dirEstimate_QuadAvNorm(st,t),stimulus.directions(st)) ;
    end
end

for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        AcuteAngleErrorQuadAvNorm2(st,t) = acuteAngle(VectSumImQuad_dirEstimate_QuadAvNorm2(st,t),stimulus.directions(st)) ;
    end
end

% error of quad with longest normalized vector (winner take all)
for st = 1:NumDir ; % for each stimulus direction
    for t=1:length(psth_mean{1}) ; % for each time point in the psth
        clear Temp
        for q = 1:size(QuadSetUnique,1) ; % for each quad    
            %Temp(q) = VectSumImQuad_vecLength{q}(st,t)*Qpsth{q}(st,t) ;
            Temp(q) = VectSumImQuad_vecLength{q}(st,t)/(Qpsth{q}(st,t)+.00001) ;
        end
        AcuteAngleErrorQuadLongestNorm(st,t) = acuteAngle(VectSumImQuad_dirEstimate{find(Temp==max(Temp),1,'first')}(st,t),stimulus.directions(st)) ; 
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

%% noise correlations
% psth residuals
for cells=1:length(dataRun.spikes) ; % for each cell
    for st = 1:length(stimulus.directions) ; % for each stimulus direction
        for t=1:NumTrialsPerDir ; % for each trial
            psth_residual{cells}{st}(t,:) = psth{cells}{st}(t,:)-psth_mean{cells}(st,:) ;
        end
    end
end

% noise correlations as function of time (CAUTION- not local residuals)
for st = 1:length(stimulus.directions) ; % for each stimulus direction
    corr{st} = nan(length(stimulus.directions),length(stimulus.directions),length(psthTime)) ;
    for cells=1:length(dataRun.spikes) ; % for each cell
        for cells2=1:length(dataRun.spikes) ; % for each cell
            for pt = 1:length(psthTime) ; % for each time point of psth
                temp = corrcoef(psth_residual{cells}{st}(:,pt)',psth_residual{cells2}{st}(:,pt)') ;
                corr{st}(cells,cells2,pt)= temp(1,2) ;
            end
        end
    end
end
                
% corr histograms 
for st = 1:length(stimulus.directions) ; % for each stimulus direction                
    corrhist(st,:) = hist(corr{st}(:),[-1:.1:1]) ; 
end


%% figures

figure % plot tuning curves (all cells)
for cells=1:length(dataRun.spikes) ; % for each cell
    errorbar(stimulus.directions,spkCount_direction_mean(cells,:),spkCount_direction_std(cells,:))
    hold on
end
   
figure % plot tuning curves (ds cells only)
for DsType = 1:length(DsTypeName) ;
    subplot(length(DsTypeName),1,DsType)
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        %errorbar(stimulus.directions,spkCount_direction_mean(DsCelli{DsType}(cells),:),spkCount_direction_std(DsCelli{DsType}(cells),:))
        plot(stimulus.directions,spkCount_direction_mean(DsCelli{DsType}(cells),:),'-*')
        hold on
    end
end
   
figure % dsi
plot(DSindex,'k*')
hold on
for DsType = 1:length(DsTypeName) ;
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        plot(DsCelli{DsType}(cells),DSindex(DsCelli{DsType}(cells)),'ro')
    end
end

figure % dsi and directions (polar)
polar([0,0],[0,1])
hold on
for cells=1:length(dataRun.spikes) ; % for each cell
    polar([PreferedDirection(cells),PreferedDirection(cells)]*pi/180,[DSindex(cells),DSindex(cells)],'*')
    hold on
end

figure % dsi and direction (linear)
plot(DSindex,PreferedDirection,'*')

figure % plot raster (ordered by stimulus direction)
NumTrialsPerDir = size(stimulus.directionsShown,1) ; % number of trial for each direction
for cells=1:length(dataRun.spikes) ; % for each cell
    [temp,ti] = sort(dirShownFullTranspose(:)) ; % get direction
    for t=1:length(ti) ;
        spk = dataRun.spikes{cells}-triggs(ti(t)) ;
        spk = spk(spk>=0 & spk<=StimDuration) ;
        plot(spk,ones(1,length(spk))*t,'k.')
        axis([0,4,0,length(triggs)]) 
        hold on
    end
    title(num2str(cells))
    hold off
    pause
end

figure % plot raster (ordered by stimulus direction) - DS cells only
for DsType = 1:length(DsTypeName) ;
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        [temp,ti] = sort(dirShownFullTranspose(:)) ; % get direction in order they were presented
        for t=1:length(ti) ;
            spk = dataRun.spikes{DsCelli{DsType}(cells)}-triggs(ti(t)) ;
            spk = spk(spk>=0 & spk<=StimDuration) ;
            plot(spk,ones(1,length(spk))*t,'k.')
            axis([0,4,0,length(triggs)]) 
            hold on
        end
        title(num2str(DsCelli{DsType}(cells)))
        
        pause
        clf
    end
end


figure % plot raster of all cells for each trial - DS cells only
for t=1:length(triggs) ;
    clf
    c=0 ;
    for DsType = 1:length(DsTypeName) ;
        for cells=1:length(DsCelli{DsType}) ; % for each cell
            c=c+1 ;
            spk = dataRun.spikes{DsCelli{DsType}(cells)}-triggs(t) ;
            spk = spk(spk>=0 & spk<=StimDuration) ;
            plot(spk,ones(1,length(spk))*c,[Color_list{DsType},'.']) 
            hold on
        end
    end
    axis([0,4,0,c])
    title(num2str(t))
    
    pause
end

figure % plot raster for one stim direction - DS cells only
for DsType = 1:length(DsTypeName) ;
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        temp = dirShownFullTranspose ;
        [temp,ti] = sort(temp(:)) ; % get direction
        for t=1:length(ti) ;
            spk = dataRun.spikes{DsCelli{DsType}(cells)}-triggs(ti(t)) ;
            spk = spk(spk>=0 & spk<=StimDuration) ;
            plot(spk,ones(1,length(spk))*t,'r.')
            axis([0,4,0,length(triggs)]) 
            hold on
        end
        title(num2str(DsCelli{DsType}(cells)))
        
        pause
        clf
    end
end

figure % plot raster of all cells for each trial
for t=1:length(triggs) ;
    for cells=1:length(dataRun.spikes) ; % for each cell
        spk = dataRun.spikes{cells}-triggs(t) ;
        spk = spk(spk>=0 & spk<=StimDuration) ;
        plot(spk,ones(1,length(spk))*cells,color)
        axis([0,4,0,length(triggs)]) 
        hold on
    end
    title(num2str(cells))
    hold off
    pause
end

figure % grating ds index vs image slide dsi
for cells=1:length(dataRun.spikes) ; % for each cell
    plot(mean([ForIgor.dsi_all{1}(cells),ForIgor.dsi_all{2}(cells)]), DSindex(cells),'*')
    hold on
end
xlabel('grading dsi')
ylabel('image dsi')

figure % noise correlations
for st = 1:length(stimulus.directions) ;
    for cells=1:length(dataRun.spikes) ; % for each cell
        for cells2=1:length(dataRun.spikes) ; % for each cell
            plot(sqrt(psth_mean{cells}(st,:).*psth_mean{cells2}(st,:)),squeeze(corr{st}(cells,cells2,:)),'*')
            pause
        end
    end
end
   
figure % noise correlations
for st = 1:length(stimulus.directions) ;
    for cells=1:length(dataRun.spikes) ; % for each cell
        for cells2=1:length(dataRun.spikes) ; % for each cell
            plot(psthTime,psth_mean{cells}(st,:),'b')
            hold on
            plot(psthTime,psth_mean{cells2}(st,:),'r')
            plot(psthTime,squeeze(corr{st}(cells,cells2,:)),'k')
            pause
            clf
        end
    end
end

figure % DS cell EI positions  
for DsType = 1:length(DsTypeName) ;
    %subplot(length(DsTypeName),1,DsType)
    figure
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        plot(EiCtr(DsCelli{DsType}(cells),1),EiCtr(DsCelli{DsType}(cells),2),[Color_list{DsType},'+'])
        drawCircle(EiCtr(DsCelli{DsType}(cells),1),EiCtr(DsCelli{DsType}(cells),2),50,'color',Color_list{DsType})
        hold on
    end
    axis([-500 500 -500 500])
end

figure % DS EI positions firing rate movie
TypeAngle = [270,0,90,180] ;
EiSpikeRateFactor = 50 ; % (um/spike)
for st = 1:length(stimulus.directions) ; 
    title(num2str(st))
    for PsthPnt=1:length(psthTime) ; % for each point in psth
        for DsType = 1:length(DsTypeName) ;
            for cells=1:length(DsCelli{DsType}) ; % for each cell
                plot(EiCtr(DsCelli{DsType}(cells),1),EiCtr(DsCelli{DsType}(cells),2),[Color_list{DsType},'o'])
                hold on
                
                spikeNorm = EiSpikeRateFactor*psth_mean{DsCelli{DsType}(cells)}(st,PsthPnt) ; % normalized spike rate

                deltaX = spikeNorm*cosd(TypeAngle(DsType)) ;
                deltaY = spikeNorm*sind(TypeAngle(DsType)) ;

                X=[EiCtr(DsCelli{DsType}(cells),1),EiCtr(DsCelli{DsType}(cells),1)+deltaX] ;
                Y=[EiCtr(DsCelli{DsType}(cells),2),EiCtr(DsCelli{DsType}(cells),2)+deltaY] ;

                plot(X,Y,'Color',Color_list{DsType})
            end
        end
        pause
        clf
    end 
end

figure % DG tuning curves - divided by tp
for DsType = 1:length(DsTypeName) ; % ds type
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        plot(dataRunDg.stimulus.params.DIRECTION, DgTuningCurve{1}(DsCelli{DsType}(cells),:),'k*-')
        hold on
        plot(dataRunDg.stimulus.params.DIRECTION, DgTuningCurve{2}(DsCelli{DsType}(cells),:),'r*-')
        title(['DsType: ',num2str(DsType), ' cell: ',num2str(cells)])
        hold off
        pause
    end
end

figure % DG vector sums
for DsType = 1:length(DsTypeName) ; % ds type
    for cells=1:length(DsCelli{DsType}) ; % for each cell
        polar(VectSumDg(DsCelli{DsType}(cells),1)*pi/180,VectSumDg(DsCelli{DsType}(cells),2),['*',Color_list{DsType}])
        hold on
    end
end

figure % direction estimate from vector sum
for st = 1:length(stimulus.directions) ; 
    plot([psthTime(1),psthTime(end)],[1,1]*[stimulus.directions(st)],[':',Color_list{st}])
    hold on
    plot(psthTime,DirEstimate(st,:),['-',Color_list{st}])
    %plot(psthTime,VectSumImQuad_dirEstimate_QuadAv(st,:),['',Color_list{st}])
    %plot(psthTime,VectSumImQuad_dirEstimate_QuadAvNorm(st,:),['--',Color_list{st}])
    %plot(psthTime,VectSumImQuad_dirEstimate_QuadAvNorm2(st,:),[':',Color_list{st}])
    
%     pause
%     hold off
end
xlabel('time (sec)')
ylabel('direction (deg)')

figure % looking at quad data
for q=1:length(QuadSetUnique) ; % for each ds quad
    for st = 1:length(stimulus.directions) ;
        plot(VectSumImQuad_vecLength{q}(st,:))
        hold on
        plot(Qpsth{q}(st,:))
        hold off
        pause
    end
end

figure % direction estimate from vector sum - for ds typing dg 
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus
    plot(psthTimeDg,DirEstimateDg(st,:),['-',Color_list{st}])
    hold on
    plot([psthTimeDg(1),psthTimeDg(end)],[1,1]*[dataRunDg.stimulus.combinations(st).DIRECTION],[':',Color_list{st}])
      %pause
      %hold off
end
xlabel('time (sec)')
ylabel('direction (deg)')

figure % direction estimate from 1 select quad - for ds typing dg
q=1 ; % select quad
for st = 1:length(dataRunDg.stimulus.combinations) ; % for each stimulus
%for st = [1:4,6,9,15,16] ; % for each stimulus
%for st = [5,7:8,10:14] ; % for each stimulus
    plot(psthTimeDg,VectSumDgQuad_dirEstimate{q}(st,:),['-',Color_list{st}])
    hold on
    plot([psthTimeDg(1),psthTimeDg(end)],[1,1]*[dataRunDg.stimulus.combinations(st).DIRECTION],[':',Color_list{st}])
      pause
      hold off
end
xlabel('time (sec)')
ylabel('direction (deg)')

figure % direction estimate average
plot(stimulus.directions,mean(DirEstimate,2),'*')
xlabel('direction')
ylabel('mean direction estimate')

figure % direction estimate wieghted average by vector sum length
for st = 1:length(stimulus.directions) ; 
    plot(stimulus.directions(st),DirEstimate(st,:)*DirEstimate_vecLength(st,:)'/sum(DirEstimate_vecLength(st,:)),'ro')
    hold on
end
xlabel('direction')
ylabel('mean direction estimate')

figure % error in dir estimate
TempError = AcuteAngleError(:,ceil(StimOnOffset/diff(psthTime(1:2))):end) ;
%hist(AcuteAngleError(:),[0:3:180])
hist(TempError(:),[0:3:180])
ylabel('number of observations (time bins)')
xlabel('direction error (deg)')

figure % error in dir estimate - ds typing dg
TempError = AcuteAngleErrorDg(:,ceil(StimOnOffset/diff(psthTimeDg(1:2))):ceil(DgOffTime/diff(psthTimeDg(1:2)))) ;
hist(TempError(:),[0:3:180])
ylabel('number of observations (time bins)')
xlabel('direction error (deg)')

figure % error of each quad - image
for q=1:length(QuadSetUnique) ;
    TempError = AcuteAngleErrorQuads{q}(:,ceil(StimOnOffset/diff(psthTime(1:2))):end) ;
    TempHist = hist(TempError(:),[0:3:180]) ;
    
    subplot(2,1,1)
    plot([0:3:180],TempHist/sum(TempHist))
    hold on
    
    subplot(2,1,2)
    plot(q,mean(TempError(:)),'*')
    hold on
end
 
figure % error of each quad - dg
for q=1:length(QuadSetUnique) ;
    TempError = AcuteAngleErrorQuadsDg{q}(:,ceil(StimOnOffset/diff(psthTimeDg(1:2))):ceil(DgOffTime/diff(psthTimeDg(1:2)))) ;
    TempHist = hist(TempError(:),[0:3:180]) ;
    
    subplot(2,1,1)
    plot([0:3:180],TempHist/sum(TempHist))
    hold on
    
    subplot(2,1,2)
    plot(q,mean(TempError(:)),'*')
    hold on
end

figure % error in dir estimate - comparing

plot([0:3:180],histIm/sum(histIm),'k') % all cells
hold on
plot([0:3:180],histQSumIm/sum(histQSumIm),'g') % quad sum
plot([0:3:180],histQBestIm/sum(histQBestIm),'r') % best of the quads

plot([0:3:180],histQUnitAvIm/sum(histQUnitAvIm),'g-.') % 
plot([0:3:180],histQnormIm/sum(histQnormIm),'g--')
plot([0:3:180],histQnorm2Im/sum(histQnorm2Im),'g:')

plot([0:3:180],histQLongestIm/sum(histQLongestIm),'r--')
plot([0:3:180],histQStrongestIm/sum(histQStrongestIm),'r:')
plot([0:3:180],histQLongestNormIm/sum(histQLongestNormIm),'r-.')
plot([0:3:180],histQMostSpreadIm/sum(histQMostSpreadIm),'y')

plot([0:3:180],histDg/sum(histDg),'b')

ylabel('number of observations (time bins)')
xlabel('direction error (deg)')

figure % best quad index
subplot(2,1,1)
plot(psthTime,BestQuadi(1,:),'k') ;
hold on
plot(psthTime,BestQuadi(2,:),'b') ;
plot(psthTime,BestQuadi(5,:),'r') ;
xlabel('time (sec)')
ylabel('best quad index')

subplot(2,1,2)
hist(BestQuadi(:),[1:length(QuadSetUnique)])
xlabel('best quad index')
ylabel('number observations')

figure % comparing stats for each quad as function of error
statstr = {'VectSumImQuad_vecLength{q}(st,:)','DsQuadSum{q}(st,:)','DsQuadSpreadIndex{q}(st,:)'...
    'DsQuadSpreadIndex2{q}(st,:)','Qpsth{q}(st,:)','VectSumImQuad_vecLength{q}(st,:)./Qpsth{q}(st,:)'...
    'VectSumImQuad_vecLength{q}(st,:).*Qpsth{q}(st,:)','(DsQuadSum{q}(st,:)/4)./Qpsth{q}(st,:)',...
    '(DsQuadSum{q}(st,:)/4).*Qpsth{q}(st,:)','VectSumImQuad_vecLength{q}(st,:).*DsQuadSpreadIndex2{q}(st,:)',...
    'VectSumImQuad_vecLength{q}(st,:).*DsQuadSpreadIndex{q}(st,:)'} ;
for statNum = 1:length(statstr) ;
    for q=1:length(QuadSetUnique) ; % for each ds quad
        for st = 1:length(stimulus.directions) ;
            eval(['plot(AcuteAngleErrorQuads{q}(st,:),',statstr{statNum},',''k.'')'])
            hold on
        end
    end
    xlabel('quad error (deg)')
    ylabel(statstr{statNum}) 
    pause
    print(gcf, '-djpeg',[saveFigPath,'DsStatsFigsDB',num2str(DB),'20msPsth','stat',num2str(statNum)])
    hold off
end
        

% for joel

ForJoel(DB,ImPathNum,DsPathNum).stimulus = struct(stimulus) ;
ForJoel(DB,ImPathNum,DsPathNum).stimulus.directionsShownFull = directionsShownFull ;
ForJoel(DB,ImPathNum,DsPathNum).psth = psth ; % {cell}{stimulus}(trial,:)
ForJoel(DB,ImPathNum,DsPathNum).psth_mean = psth_mean ; % {cell}(stimulus,:)
ForJoel(DB,ImPathNum,DsPathNum).psth_timeBins = psthTime ; % 
ForJoel(DB,ImPathNum,DsPathNum).EiCenters = EiCtr ; % (cellX,cellY)
ForJoel(DB,ImPathNum,DsPathNum).DsCellsIndices = DsCelli ; % {DsType}(cell)
ForJoel(DB,ImPathNum,DsPathNum).GratingDsi = ForIgor.dsi_all ; % {grating speed}(cell)
ForJoel(DB,ImPathNum,DsPathNum).psthDg_timeBins = psthTimeDg ; % 
ForJoel(DB,ImPathNum,DsPathNum).psthDg = psthDg ; % {cell}{stimulus}(trial,:)
ForJoel(DB,ImPathNum,DsPathNum).psthDg_mean = psthDg_mean ; % {cell}(stimulus,:)
ForJoel(DB,ImPathNum,DsPathNum).stimulusDg = dataRunDg.stimulus.combinations ; %


