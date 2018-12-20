function ForIgor = ImJitterAnalysisV3(DataBlock, DB, Params)

% edited from "ImJitterAnalysisV2" to allow for psth time scale and focus
% on paper topics

% JC 5/23/2018

RunAsScript = true ;
if RunAsScript ;
    DB = 19 ; % 20,19,8
    [DataBlock,Params] = DataBlocks_NaturalStim ;
end

ImPathNum = 1 ; 
DsPathNum = 1 ;

% parameters
PsthBinNum = 100 ; % number of frame bins in each psth bin

saveFigPath = ['/Users/jcafaro/Documents/AnalysisFigures/'] ;
saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/DsIdsDb',num2str(DB),'ImJitterPathNum',num2str(ImPathNum)] ;

Color_list = {'c','r','b','g','k','y','c','r','b','g','k','y','c','r','b','g','k','y',...
    'c','r','b','g','k','y','c','r','b','g','k','y','c','r','b','g','k','y'} ; % order of colors for each 
FrameRate = 60.35 ; % rate of monitor display
FramesPerTrig = 100 ; % number of frames between each trigger
StaTime = [-3:1/FrameRate:3] ;% (sec) time before spike:time after spike

% load data
dataRun = load_data(DataBlock(DB).ImJitterConcat{ImPathNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun.piece.array_id = 1551 ; % CAUTION
dataRun = load_ei(dataRun, 'all') ;

% identify DS cells
dataRunTemp = load_data(DataBlock(DB).DsPath{DsPathNum}) ;
dataRunTemp = load_neurons(dataRunTemp) ;

Params.TimeBounds = [0,dataRunTemp.duration] ; % times of Dg data
Params.SkipBwMap = true ;
Params.DsPathNum = DsPathNum ;
DataBlock(DB).DsConcatPath = DataBlock(DB).ImJitterConcat{ImPathNum} ;
clear dataRunTemp

try load(saveDsIdsPath) ; % if they are ds ids already saved
catch % if not find them

    ForIgor = DsCellFinder(DataBlock, DB, Params) ;

    cell_id = ForIgor.ds_id{2} ; % On-Off cell 
    DsTypeName = ForIgor.dsName{2} ; 
    dsi = ForIgor.dsi{2} ;

    for DsType=1:length(DsTypeName) ;
        cell_i{DsType} = get_cell_indices(dataRun, cell_id{DsType}) ;
    end

    save(saveDsIdsPath, 'cell_id','DsTypeName','cell_i')
end

% load stimulus
slashi = strfind(DataBlock(DB).ImJitterConcat{ImPathNum},'/') ; % find the /
dashi = strfind(DataBlock(DB).ImJitterConcat{ImPathNum},'-') ; % find the -
StimPath = [DataBlock(DB).ImJitterConcat{ImPathNum}(1:slashi(end-1)),'stimuli/s0',DataBlock(DB).ImJitterConcat{ImPathNum}(dashi(end)+1)] ;
load(StimPath) ;

StaPnts = floor(StaTime*(FrameRate/stimulus.frame_interval)) ; 
StaTime = StaPnts/(FrameRate/stimulus.frame_interval) ; % this is the precise time (entered above an estimate)

Trigs_orig = dataRun.triggers(dataRun.triggers>Params.TimeBounds(2)) ; % Triggers in image jitter block
Trigs = Trigs_orig - Trigs_orig(1) ; % Triggers adjusted for first frame

% sta
if ~ischar(stimulus.XImageJitterVector) && stimulus.image_jitter_std>0 ; % if there was a jittered image
    XJitterVectorDiff = [0,diff(stimulus.XImageJitterVector)'] ;
    YJitterVectorDiff = [0,diff(stimulus.YImageJitterVector)'] ;
end

if ~ischar(stimulus.XSquareJitterVector) && stimulus.square_jitter_std>0 ; % if there was a jittered square
    XJitterVectorDiff = [0,diff(stimulus.XSquareJitterVector)'] ;
    YJitterVectorDiff = [0,diff(stimulus.YSquareJitterVector)'] ;
end
    
for DsType=1:length(DsTypeName) ; % for each DS cell type
    for c = 1:length(cell_i{DsType}) ; % for each DS cell

        spikeTimes_orig = dataRun.spikes{cell_i{DsType}(c)}(dataRun.spikes{cell_i{DsType}(c)}>=Trigs_orig(1)) ; % spikes during stimulus
        spikeTimes = spikeTimes_orig - Trigs_orig(1) ; % spike times relative to trial start

        numSpikes{DsType}(c)= 0 ;
        StaXJitterVector{DsType}(c,:) = zeros(1,length(StaPnts)) ;
        StaYJitterVector{DsType}(c,:) = zeros(1,length(StaPnts)) ;

        for s=1:length(spikeTimes) ; % for each spike
            if spikeTimes(s)<Trigs(end) ; % if the spike time is not past the last trigger
                T1i = find(Trigs<spikeTimes(s),1,'last') ; % index of trigger preceding spike        
                FrameRateEstimate = (Trigs(T1i+1) - Trigs(T1i))/FramesPerTrig ; % estimated time of each image frame between triggers
                FramesAfterTrig = floor((spikeTimes(s) - Trigs(T1i))/FrameRateEstimate) ; % number of frames after the trigger
                FramesBeforeTrig = (T1i-1)*FramesPerTrig ;
                spikePnt = FramesBeforeTrig + FramesAfterTrig ;
                VectorPnts = spikePnt+StaPnts ;
                if min(VectorPnts)>0 ; % if its late enough
                    StaXJitterVector{DsType}(c,:) = StaXJitterVector{DsType}(c,:) + XJitterVectorDiff(VectorPnts) ;
                    StaYJitterVector{DsType}(c,:) = StaYJitterVector{DsType}(c,:) + YJitterVectorDiff(VectorPnts) ;
                    numSpikes{DsType}(c) = numSpikes{DsType}(c) + 1 ;
                end
            end
        end
        StaXJitterVector{DsType}(c,:) = StaXJitterVector{DsType}(c,:)/numSpikes{DsType}(c) ;
        StaYJitterVector{DsType}(c,:) = StaYJitterVector{DsType}(c,:)/numSpikes{DsType}(c) ;

        StaMagJitterVector{DsType}(c,:) = sqrt(StaXJitterVector{DsType}(c,:).^2+...
            StaYJitterVector{DsType}(c,:).^2) ; % magnitude of movement

        StaDirJitterVector{DsType}(c,:) = cart2pol(StaXJitterVector{DsType}(c,:),...
            StaYJitterVector{DsType}(c,:)) ; % direction of movement (radians)              

        [m,mi] = max(StaMagJitterVector{DsType}(c,:)) ;
        StaMagPeakJitterVector{DsType}(c) = m ;
        StaDirAtPeakJitterVector{DsType}(c) = StaDirJitterVector{DsType}(c,mi) ; % direction at peak of movement mag           
    end
end

% distribution of directions vs average for all cells
AngleDiffDistX = [5:10:180] ; % 
AngleDiffDist = AngleDiffDistX*0 ; % prep dist
for DsType=1:length(DsTypeName) ; % for each DS cell type
    for c = 1:length(cell_i{DsType}) ; % for each DS cell

        spikeTimes_orig = dataRun.spikes{cell_i{DsType}(c)}(dataRun.spikes{cell_i{DsType}(c)}>=Trigs_orig(1)) ; % spikes during stimulus
        spikeTimes = spikeTimes_orig - Trigs_orig(1) ; % spike times relative to trial start

        [m,mi] = max(StaMagJitterVector{DsType}(c,:)) ;
        MeanDirTemp = StaDirJitterVector{DsType}(c,mi)*180/pi ; % direction at peak of movement mag  
        
        for s=1:length(spikeTimes) ; % for each spike
            if spikeTimes(s)<Trigs(end) ; % if the spike time is not past the last trigger
                T1i = find(Trigs<spikeTimes(s),1,'last') ; % index of trigger preceding spike        
                FrameRateEstimate = (Trigs(T1i+1) - Trigs(T1i))/FramesPerTrig ; % estimated time of each image frame between triggers
                FramesAfterTrig = floor((spikeTimes(s) - Trigs(T1i))/FrameRateEstimate) ; % number of frames after the trigger
                FramesBeforeTrig = (T1i-1)*FramesPerTrig ;
                spikePnt = FramesBeforeTrig + FramesAfterTrig ;
                VectorPnts = spikePnt+StaPnts ;
                if min(VectorPnts)>0 ; % if its late enough
                    DirVectTemp = cart2pol(XJitterVectorDiff(VectorPnts),...
                        YJitterVectorDiff(VectorPnts))*180/pi ;
                    AngleDiffDistTemp = hist(acuteAngle(MeanDirTemp,DirVectTemp(mi)),AngleDiffDistX) ;
                    AngleDiffDist = AngleDiffDist+AngleDiffDistTemp ;
                end
            end
        end
    end
end

% find frame time 
frame_timeBins = 0 ;
for t = 1:length(Trigs)-1 ; % for each trigger
    FrameRateEstimate = (Trigs(t+1) - Trigs(t))/FramesPerTrig ; % estimated time of each image frame between triggers
    frame_timeBins = [frame_timeBins,[Trigs(t)+FrameRateEstimate:FrameRateEstimate:Trigs(t+1)]] ; % bins do not go beyond last trigger (though data may)
end
NumStimPnts = length(frame_timeBins) ; % number of stim frames actually delivered

% psth - at each frame
psthTimeBin = mean(diff(frame_timeBins))*PsthBinNum ; 

for DsType=1:length(DsTypeName) ; % for each DS cell type      
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        spk = dataRun.spikes{cell_i{DsType}(c)}(dataRun.spikes{cell_i{DsType}(c)}>=Trigs_orig(1)) ; % spikes during stimulus
        spk = spk - Trigs_orig(1) ; % spike times relative to trial start

        for tb = 1:NumStimPnts ; % for each frame bin
            strtPnt = max([1,tb-PsthBinNum]) ; 
            psth{DsType}(c,tb) = sum((spk>frame_timeBins(strtPnt)).*(spk<=frame_timeBins(tb)))/psthTimeBin ; % spikes per sec within window
        end    
    end
end

% spike rate of each cell
for DsType=1:length(DsTypeName) ; % for each DS cell type      
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        SpikeRateAvOverStim{DsType}(c) = numSpikes{DsType}(c)/Trigs(end) ; % spikes/sec
    end
end
SpikeRateAvOverStim_max = max(cell2mat(SpikeRateAvOverStim)) ;
SpikeRateAvOverStim_m = min(cell2mat(SpikeRateAvOverStim)) ;

% fraction of psth bins with response>0
Temp = cell2mat(psth') ;
FracAboveZero = sum(Temp(:)>0)/length(Temp(:)) ;


%% peak to trough ratio and half width of Sta

% peak to trough ratio 
for DsType=1:length(DsTypeName) ; % for each DS cell type      
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        StaXMax = max(StaXJitterVector{DsType}(c,:)) ;
        StaXMin = min(StaXJitterVector{DsType}(c,:)) ;
        StaXMaxOverMin{DsType}(c) = max(abs(StaXMax),abs(StaXMin))/min(abs(StaXMax),abs(StaXMin)) ; % dominatant peak over minor peak
        
        StaYMax = max(StaYJitterVector{DsType}(c,:)) ;
        StaYMin = min(StaYJitterVector{DsType}(c,:)) ;
        StaYMaxOverMin{DsType}(c) = max(abs(StaYMax),abs(StaYMin))/min(abs(StaYMax),abs(StaYMin)) ; % dominatant peak over minor peak  
    end
end

MaxOverMinMat = [cell2mat(StaXMaxOverMin),cell2mat(StaYMaxOverMin)] ;
StaMaxOverMin_mean = mean(MaxOverMinMat) ;
StaMaxOverMin_sem = std(MaxOverMinMat)/sqrt(length(MaxOverMinMat)) ;

% half width
for DsType=1:length(DsTypeName) ; % for each DS cell type      
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        StaXAbs = abs(StaXJitterVector{DsType}(c,:)) ;
        [StaXAbsMax,StaXAbsMaxi] = max(StaXAbs) ;
        p1 = find(StaXAbs(1:StaXAbsMaxi)>StaXAbsMax/2,1,'first') ; % first point above half max
        p2 = find(StaXAbs(StaXAbsMaxi:end)<StaXAbsMax/2,1,'first') + StaXAbsMaxi - 1; % first point below half max
        StaXhalfwidth{DsType}(c) = p2-p1 ;
        
        StaYAbs = abs(StaYJitterVector{DsType}(c,:)) ;
        [StaYAbsMax,StaYAbsMaxi] = max(StaYAbs) ;
        p1 = find(StaYAbs(1:StaYAbsMaxi)>StaYAbsMax/2,1,'first') ; % first point above half max
        p2 = find(StaYAbs(StaYAbsMaxi:end)<StaYAbsMax/2,1,'first') + StaYAbsMaxi - 1; % first point below half max
        StaYhalfwidth{DsType}(c) = p2-p1 ;
    end
end

HalfwidthMat = [cell2mat(StaXhalfwidth),cell2mat(StaYhalfwidth)] ;
StaHalfwidth_mean = mean(HalfwidthMat) ;
StaHalfwidth_sem = std(HalfwidthMat)/sqrt(length(HalfwidthMat)) ;

%% find quads and calculate quad response histograms 
QuadEiCircRad = 200 ; % minimum pairwise distance to be in quad
numElectrodeLayers = 2 ;

% get ei center of mass 
for cells = 1:length(dataRun.spikes) ; % for each cell
    EiCtr(cells,:) = get_ei_com(dataRun, dataRun.cell_ids(cells), numElectrodeLayers) ;
end

% get ei center of mass and make TypeCellMat so can figure out which center
% belongs to which cell
CellCnt = 1 ;
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
        ctr(CellCnt,:) = EiCtr(cell_i{DsType}(cells),:) ;
        
        TypeCellMat(CellCnt,:) = [DsType,cells] ;

        CellCnt = CellCnt+1 ;
    end
end

% distance between cells
CellCnt = 1 ; % start cell counter
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
        deltaCtr = ctr-repmat(ctr(CellCnt,:),size(ctr,1),1) ; % X,Y distance between that cell and all others
        distCtr{CellCnt} = sqrt(deltaCtr(:,1).^2+deltaCtr(:,2).^2) ; % euclidean distance
        [~,si(CellCnt,:)] = sort(distCtr{CellCnt}) ; % in order of nearest to farthest

        CellCnt = CellCnt+1 ;
    end
end

% find potential quadruplets
CellCnt = 1 ; % start cell counter
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
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
            qloc(cells,:) = EiCtr(cell_i{quad{q}(cells,1)}(quad{q}(cells,2)),:) ; % location of each cell in quad
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

for TEMP_LOOP = 1:19 ;
%% OLE (optimal linear estimator)
UseQuadFlag = false ; % use single quad only
RealQuadToUse = TEMP_LOOP ;
BestCellTestFlag = false ; % use only the strongest responding cell at any time point in estimate

clear OleCellsi ; % in case its been run before
if UseQuadFlag % use a quad
    RandomQuadFlag = false ; % construct a random quadruplet

    if ~RandomQuadFlag  % if you want to use a real quad
        for cells=1:4  % for each cell
            OleCellsi(cells) = find(TypeCellMat(:,1)==quad{RealQuad_i(RealQuadToUse)}(cells,1),1)-1+...
                quad{RealQuad_i(RealQuadToUse)}(cells,2) ;
      
        end
    else
        for cells=1:4  % for each cell
            OleCellsi(cells) = cell_i{cells}(randi(length(cell_i{cells}))) ; % pick randomly for each type
        end
    end
else % if not using a quad
    OleCellsi = [1:length(cell2mat(cell_i))] ; % = all ds cells
end

TrainPnts = [1:floor(NumStimPnts*3/4)] ; % points to train on
TestPnts = [1:floor(NumStimPnts*1/4)]+floor(NumStimPnts*3/4); % points to test on
OptLinDecoderFlag = false ; % run optimal linear decoder code instead of OLE code
FilterLength = 20 ; % number of points to be used in optimal linear decoding 

psthMat = cell2mat(psth')' ; % mat of responses cell=column, time=row
psthMatTrain = psthMat(TrainPnts,OleCellsi) ;

MagJitterVector = sqrt(XJitterVectorDiff.^2 + YJitterVectorDiff.^2) ; % magnitude

DirJitterVector = cart2pol(XJitterVectorDiff,YJitterVectorDiff) ; % direction of movement (radians)
DirJitterVectorDeg = DirJitterVector*180/pi+180 ;

if OptLinDecoderFlag 
    [LinFiltersX,ConstantX] = MultiCellLinFilterFinder(XJitterVectorDiff(1:size(psthMatTrain,1)),psthMatTrain,FilterLength) ;
    [LinFiltersY,ConstantY] = MultiCellLinFilterFinder(YJitterVectorDiff(1:size(psthMatTrain,1)),psthMatTrain,FilterLength) ;
    Xest = MultiCellLinFilterTester(LinFiltersX,ConstantX,psthMat);
    Yest = MultiCellLinFilterTester(LinFiltersY,ConstantY,psthMat);
    Est= cart2pol(Xest,Yest)*180/pi+180 ;
else
    psthMatTest = psthMat ; % added so can manipulate psthMat for single cell
    
    for td = 0:FilterLength  % for each time delay
        [OleWeights, OlePolarWeights] = OleFinder(DirJitterVectorDeg(1:size(psthMatTrain,1)-td),psthMatTrain((td+1):end,:)) ;
        
        if BestCellTestFlag ; % if using only the strongest response at each time point
            for tp = 1:size(psthMatTest,1) ; % for each time point
                [mx,mpnt] = max(psthMatTest(tp,:)) ; % find max response
                psthMatTest(tp,:) = psthMatTest(tp,:)*0 ; % zero all responses
                psthMatTest(tp,mpnt) = mx ; % keep max response only
            end
        end
        
        [DirectionEstimate,DirectionEstimateMag] = OleTester(OleWeights,psthMatTest((td+1):end,OleCellsi));
        Est{td+1} = DirectionEstimate ;
    end
end

% error
if OptLinDecoderFlag 
    for s=1:length(TestPnts)  % for each stim bin
        EstError(s) = acuteAngle(Est(TestPnts(s)),DirJitterVectorDeg(TestPnts(s))) ;
    end
else
    for td = 0:FilterLength  % for each time delay
        for s=1:(length(TestPnts)-FilterLength)  % for each stim bin
            EstError{td+1}(s) = acuteAngle(Est{td+1}(TestPnts(s)),DirJitterVectorDeg(TestPnts(s))) ;
        end
    end
end

% median error
if OptLinDecoderFlag ;
    EstError_med = median(EstError) ;
else   
    for td = 0:FilterLength ; % for each time delay
        EstError_med(td+1) = median(EstError{td+1}) ;
    end
end

% transfer function of estimate and actual x,y
[~,Besti] = min(EstError_med) ;
EstX = -cosd(Est{Besti}([TestPnts(1)-Besti+1:(TestPnts(end)-FilterLength)])) ; % best estimate X
EstY = -sind(Est{Besti}([TestPnts(1)-Besti+1:(TestPnts(end)-FilterLength)])) ; % best estimate Y
XestTransfer = xcov(EstX,XJitterVectorDiff([TestPnts(1):(TestPnts(end)-FilterLength+Besti-1)]),'coef') ;
YestTransfer = xcov(EstY,YJitterVectorDiff([TestPnts(1):(TestPnts(end)-FilterLength+Besti-1)]),'coef') ;
  
EstError_med_LOOP(TEMP_LOOP,:) = EstError_med ;

% number of cells responding at each time point
NumCellsResp = [sum(psthMat(Besti:end,OleCellsi)>0,2);zeros(1,Besti-1)'] ;
NumCellsResp_DistX = [0:4] ;

% distribution of number of cells
NumCellsResp_Dist = hist(NumCellsResp,NumCellsResp_DistX)/length(NumCellsResp) ;
NumCellsResp_Dist_LOOP(TEMP_LOOP,:) = NumCellsResp_Dist ;

% median error for each number
for nc = 1:length(NumCellsResp_DistX) ; 
    Tempi = NumCellsResp([TestPnts(1):(TestPnts(end)-FilterLength)])==NumCellsResp_DistX(nc) ; 
    EstErrorForCellNumber(nc) = median(EstError{Besti}(Tempi)) ;
end

EstErrorForCellNumber_LOOP(TEMP_LOOP,:) = EstErrorForCellNumber ;

end ; %TEMP LOOP




%% figures

figure % spike rasters 
RasterDuration = 10 ; % (sec)
rnd = 1 ;
for DsType=1:length(DsTypeName) ; % for each DS cell type
        for c = 1:length(cell_i{DsType}) ; % for each DS cell
            for s=1:length(dataRun.spikes{cell_i{DsType}(c)}) ; % spikes
                if dataRun.spikes{cell_i{DsType}(c)}(s)>Trigs_orig(1) && dataRun.spikes{cell_i{DsType}(c)}(s)<Trigs_orig(1)+RasterDuration ;
                    line([dataRun.spikes{cell_i{DsType}(c)}(s),...
                        dataRun.spikes{cell_i{DsType}(c)}(s)],[rnd-.5,rnd+.5],...
                    'color',Color_list{DsType})
                    hold on
                end
            end
            rnd = rnd+1 ;
        end
end
xlabel('time (sec)')
ylabel('cell')
   
figure % sta 
for DsType=1:length(DsTypeName) ; % for each DS cell type
    for c = 1:length(cell_i{DsType}) ; % for each DS cell 
        subplot(3,1,1)
        plot(StaTime,StaXJitterVector{DsType}(c,:),'r')
        hold on
        plot(StaTime,StaYJitterVector{DsType}(c,:),'b')
        hold off
        title(num2str(numSpikes{DsType}(c)))
        ylabel('delta X, Y')

        subplot(3,1,2)
        plot(StaTime,StaMagJitterVector{DsType}(c,:))
        ylabel('delta magnitude')

        subplot(3,1,3)
        plot(StaTime,StaDirJitterVector{DsType}(c,:))
        ylabel('delta direction')
        xlabel('time to spike (sec)')

        pause
    end
end
  
figure % all cells sta for image
subplot(1,3,3)
polar(cell2mat(StaDirAtPeakJitterVector),cell2mat(StaMagPeakJitterVector),'*') ; % to scale plot correctly
hold on 

for DsType=1:length(DsTypeName) ; % for each DS cell type
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        subplot(1,3,1)
        plot(StaTime,StaMagJitterVector{DsType}(c,:),Color_list{DsType})
        hold on
        xlabel('time (s)')
        ylabel('magnitude')

        subplot(1,3,2)
        plot(StaTime,StaDirJitterVector{DsType}(c,:)*180/pi,Color_list{DsType})
        hold on
        xlabel('time (s)')
        ylabel('direction')

        subplot(1,3,3)
        polar(StaDirAtPeakJitterVector{DsType}(c),StaMagPeakJitterVector{DsType}(c),[Color_list{DsType},'*'])
        hold on
    end
end
 
figure % spike rasters for repeats 
RasterDuration = 10 ; % (sec)
triggs = [401:36:length(dataRun.triggers)] ;
rnd = 1 ;
for DsType=1:length(DsTypeName) ; % for each DS cell type
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        for a = 1:length(triggs) ;
            t= dataRun.triggers(triggs(a)) ;
            for s=1:length(dataRun.spikes{cell_i{DsType}(c)}) ; % spikes
                if dataRun.spikes{cell_i{DsType}(c)}(s)>t && dataRun.spikes{cell_i{DsType}(c)}(s)<t+RasterDuration ;
%                     line([dataRun.spikes{cell_i{DsType}(c)}(s),...
%                         dataRun.spikes{cell_i{DsType}(c)}(s)],[rnd-.5,rnd+.5],...
%                     'color',Color_list{DsType})
                    plot(dataRun.spikes{cell_i{DsType}(c)}(s)-t,rnd,[Color_list{DsType},'.'])
                    hold on
                end
            end
            rnd = rnd+1 ;
        end
    end
end

% OLE estimates
figure
for td=1:20 ;
    plot(DirJitterVectorDeg(1:length(Est{td})),'k')
    hold on
    plot(Est{td},'r')
    hold off
    pause
end

figure
for td=1:20 ;
    plot(td,mean(EstError{td}),'b*')
    %plot(frame_timeBins(td),mean(EstError{td}),'b*')
    hold on
end
ylabel('average error (deg)')
xlabel('psth time bins (post direction)') 
saveas(gcf,'JitterErrorLagFig')
print(gcf, '-dpdf','JitterErrorLagFig')

figure % OLE best estimate and error and spikes
displayX = [3000,3002] ; % time
subplot(3,1,1) % spike raster
rnd = 1 ;
for cells = 1:length(OleCellsi) % for each cell used in ole
    DsType = TypeCellMat(OleCellsi(cells),1) ; % type of ds cell
    c = TypeCellMat(OleCellsi(cells),2) ; % cell 
    
    spk = dataRun.spikes{cell_i{DsType}(c)}(dataRun.spikes{cell_i{DsType}(c)}>=Trigs_orig(1)) ; % spikes during stimulus
    spk = spk - Trigs_orig(1) ; % spike times relative to trial start
    spk = spk  ;
    
    for s=1:length(spk) % spikes
        if spk(s)>displayX(1) && spk(s)<displayX(2) ;
            line([spk(s),spk(s)],[rnd-.5,rnd+.5],...
            'color','k')
            hold on
        end
    end
    rnd = rnd+1 ;
end
xlim([displayX])
xlabel('time (sec)')
ylabel('cell')

% subplot(4,1,2) % psth
% plot(frame_timeBins,psthMat(:,OleCellsi))
% xlim([displayX])

subplot(3,1,2) % estimate
plot(frame_timeBins,DirJitterVectorDeg(1:NumStimPnts),'k')
hold on
plot(frame_timeBins(Besti:end),Est{Besti},'g')
xlim([displayX])
ylim([0 360])

subplot(3,1,3) % error
plot(frame_timeBins(TestPnts(1:end-100)),EstError{Besti}(1:length(TestPnts)-100),'g')
xlim([displayX])
ylim([0 180])

hgsave('EstErrorSpikesDb19_pop_1and10bin')
print(gcf, '-dpdf', 'EstErrorSpikesDb19_pop_1and10bin')


figure % OLE median error as function of lag
plot(EstError_med)
plot(EstError_med_LOOP','b--')
hold on
ylabel('median error (deg)')
xlabel('stim frames') 

figure % cross correlation
plot([1:length(XestTransfer)]-ceil(length(XestTransfer)/2),XestTransfer)
hold on
plot([1:length(YestTransfer)]-ceil(length(YestTransfer)/2),YestTransfer)

figure % DS cell EI positions  
for DsType = 1:length(DsTypeName) ;
    %subplot(length(DsTypeName),1,DsType)
    %figure
    for cells=1:length(cell_i{DsType}) ; % for each cell
        plot(EiCtr(cell_i{DsType}(cells),1),EiCtr(cell_i{DsType}(cells),2),[Color_list{DsType},'+'])
        drawCircle(EiCtr(cell_i{DsType}(cells),1),EiCtr(cell_i{DsType}(cells),2),150,'color',Color_list{DsType})
        hold on
    end
    axis([-500 500 -500 500])
end



figure % plot EI positions for quads
for ri = 1:length(RealQuad_i) ; % for each quad
    for cells = 1:4 ; % for each cell in the quad
        ct = quad{RealQuad_i(ri)}(cells,1) ; % cell type
        cn = quad{RealQuad_i(ri)}(cells,2) ; % cell number
        plot(EiCtr(cell_i{ct}(cn),1),EiCtr(cell_i{ct}(cn),2),[Color_list{ct},'+'])
        drawCircle(EiCtr(cell_i{ct}(cn),1),EiCtr(cell_i{ct}(cn),2),150,'color',Color_list{ct})
        hold on 
    end
    axis([-600 600 -600 600])
    pause
    hold off
end

figure % direction distribution
bar(AngleDiffDistX,AngleDiffDist/sum(AngleDiffDist),'k')
hold on
plot([1,1]*(AngleDiffDistX*AngleDiffDist')/sum(AngleDiffDist),[0,.1],'k:')
xlabel('Angle off mean (deg)')
ylabel('fraction')

hgsave('DirectionDistDb19')
print(gcf, '-dpdf', 'DirectionDistDb19')

figure % cell number distribition
errorbar(NumCellsResp_DistX,mean(NumCellsResp_Dist_LOOP),std(NumCellsResp_Dist_LOOP))
xlabel('Number of active cells')
ylabel('fraction')

hgsave('NumCellsActiveDistDb19')
print(gcf, '-dpdf', 'NumCellsActiveDistDb19')

figure % error vs cell number 
errorbar(NumCellsResp_DistX,mean(EstErrorForCellNumber_LOOP),std(EstErrorForCellNumber_LOOP))
xlabel('Number of active cells')
ylabel('median error (deg)')

hgsave('ErrorVsNumCellsActiveDistDb19')
print(gcf, '-dpdf', 'ErrorVsNumCellsActiveDistDb19')


