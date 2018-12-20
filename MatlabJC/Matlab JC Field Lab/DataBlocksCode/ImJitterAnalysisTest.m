
% data set select
RunAsScript = true ;
if RunAsScript ;
    DB = 19 ; % 20,19
    [DataBlock,Params] = DataBlocks_NaturalStim ;
end

ImPathNum = 1 ; 

% flags and parameters 
UseImportDsIdsPath = false ; % use ImportDsIdsPath to id DS cells
OnDsFlag = false ; % if true, analayze ON not ON-OFF DS cells

PsthBinNum = 1 ; % number of frame bins in each psth bin
StimBinNum = 1 ; % number of frame bins in Stim to estimate

QuadEiCircRad = 200 ; % minimum pairwise distance to be in quad
numElectrodeLayers = 2 ; % for EI com calculation

Color_list = {'c','r','b','g','k','y','c','r','b','g','k','y','c','r','b','g','k','y',...
    'c','r','b','g','k','y','c','r','b','g','k','y','c','r','b','g','k','y'} ; % order of colors for each 
FramesPerTrig = 100 ; % number of frames between each trigger

% save and loads paths
saveFigPath = ['/Users/jcafaro/Documents/AnalysisFigures/'] ;
saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/ImJitterAnalysisV4DsSelection/DsIdsDb',num2str(DB),'ImJitterPathNum',num2str(ImPathNum)] ;
ImportDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/ImJitterAnalysisV4DsSelection/ImportDsIdsDb',num2str(DB),'ImJitterPathNum',num2str(ImPathNum)] ;

% load data
dataRun = load_data(DataBlock(DB).ImJitter{ImPathNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun.piece.array_id = 1551 ; % CAUTION - assumes dense array
dataRun = load_ei(dataRun, 'all') ;

% load stimulus
slashi = strfind(DataBlock(DB).ImJitter{ImPathNum},'/') ; % find the /
dashi = strfind(DataBlock(DB).ImJitter{ImPathNum},'-') ; % find the -
StimPath = [DataBlock(DB).ImJitter{ImPathNum}(1:slashi(end-1)),'stimuli/s0',DataBlock(DB).ImJitter{ImPathNum}(dashi(end)+1)] ;
load(StimPath) ;

Trigs= dataRun.triggers ; % Triggers in image jitter block

% stimulus velocity X,Y
if ~ischar(stimulus.XImageJitterVector) && stimulus.image_jitter_std>0 ; % if there was a jittered image
    XJitterVectorDiff = [0,diff(stimulus.XImageJitterVector)'] ;
    YJitterVectorDiff = [0,diff(stimulus.YImageJitterVector)'] ;
end

if ~ischar(stimulus.XSquareJitterVector) && stimulus.square_jitter_std>0 ; % if there was a jittered square
    XJitterVectorDiff = [0,diff(stimulus.XSquareJitterVector)'] ;
    YJitterVectorDiff = [0,diff(stimulus.YSquareJitterVector)'] ;
end
   
% bin stimulus (same as psth - bin is average of X bins preceding)
if StimBinNum>1 ; % if there is binning
    frame_timeBins = 0 ;
    for tb = 1:length(XJitterVectorDiff) ; % for each stim point
        strtPnt = max([1,tb-StimBinNum+1]) ; %
        XJitterVectorDiff(tb) = mean(XJitterVectorDiff(strtPnt:tb)) ;
        YJitterVectorDiff(tb) = mean(YJitterVectorDiff(strtPnt:tb)) ;
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

psth= nan(length(dataRun.spikes),NumStimPnts) ; % prep matrix for speed
for c = 1:length(dataRun.spikes) ; % for each DS cell
    spk = dataRun.spikes{c} ; % spikes during stimulus
    spk = spk - Trigs(1) ; % spike times relative to trial start

    psth(c,:) = hist(spk,frame_timeBins) ;  
end
psth = psth/psthTimeBin ;

% spike rate of each cell
for DsType=1:length(DsTypeName) ; % for each DS cell type      
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        SpikeRateAvOverStim{DsType}(c) = mean(psth{DsType}(c,:)) ; % spikes/sec
    end
end
SpikeRateAvOverStim_max = max(cell2mat(SpikeRateAvOverStim)) ;
SpikeRateAvOverStim_min = min(cell2mat(SpikeRateAvOverStim)) ;
SpikeRateAvOverStim_mean = mean(cell2mat(SpikeRateAvOverStim)) ;

% fraction of psth bins with response>0
for DsType=1:length(DsTypeName) ; % for each DS cell type      
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        FracAboveZero{DsType}(c) = sum(psth{DsType}(c,:)>0)/length(psth{DsType}(c,:)) ;
    end
end
FracAboveZero_mean = mean(cell2mat(FracAboveZero)) ;

% sta
for c = 1:length(dataRun.spikes) ; % for each DS cell
    tempCorrX = xcov(psth(c,:),XJitterVectorDiff) ; % 
    tempCorrY = xcov(psth(c,:),YJitterVectorDiff) ; %

    StaXJitterVector(c,:) = tempCorrX()/sum(psth(c,:)) ;
    StaYJitterVector(c,:) = tempCorrY()/sum(psth(c,:)) ;
end

for c = 1:length(dataRun.spikes) ; % for each DS cell

    StaXJitterVector(c,:) = StaXJitterVector(c,:)-mean(StaXJitterVector) ;
    StaYJitterVector(c,:) = StaYJitterVector(c,:)-mean(StaYJitterVector) ;
    
    StaMagJitterVector(c,:) = sqrt(StaXJitterVector(c,:).^2+...
        StaYJitterVector(c,:).^2) ; % magnitude of movement

    StaDirJitterVector(c,:) = cart2pol(StaXJitterVector(c,:),...
        StaYJitterVector(c,:)) ; % direction of movement (radians)              

    [m,mi] = max(StaMagJitterVector(c,:)) ;
    StaMagPeakJitterVector(c) = m ;
    StaDirAtPeakJitterVector(c) = StaDirJitterVector(c,mi) ; % direction at peak of movement mag 
end

StaPnts = [1:length(tempCorrX)]-ceil(length(tempCorrX)/2) ; % X-axis of sta in points
StaWindowPnts = find(StaPnts>=0 & StaPnts<50) ; % window around sta peak for analysis
StaWindowPnts2 = find(StaPnts>=-50 & StaPnts<0) ; 

for c = 1:length(dataRun.spikes) ; % for each DS cell
    NormMax(c) = max(StaMagJitterVector(c,StaWindowPnts))/std(StaMagJitterVector(c,StaWindowPnts2)) ;
    [m,mi] = max(StaMagJitterVector(c,StaWindowPnts)) ;
    p1 = find(StaMagJitterVector(c,StaWindowPnts)>m/2,1,'first') ;
    p2 = find(StaMagJitterVector(c,StaWindowPnts(mi:end))<m/2,1,'first') ;
    if ~isempty(p1) && ~isempty(p2) ;
        halfmax(c) = mi+p2-p1 ; 
    end
    
end
    
    




% distribution of directions vs average for all cells (assumes StaWindows starts at 0)
DirVect = cart2pol(XJitterVectorDiff,YJitterVectorDiff)*180/pi ; % direction of stimulus (deg)
AngleDiffDistX = [0:10:180] ; % 
AngleDiffDist = AngleDiffDistX*0 ; % prep dist
for DsType=1:length(DsTypeName) ; % for each DS cell type
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        
        [m,mi] = max(StaMagJitterVector{DsType}(c,StaWindowPnts)) ; % peak time of mag vector
        MeanDirTemp = StaDirJitterVector{DsType}(c,StaWindowPnts(mi))*180/pi ; % direction at peak of movement mag  
        deltaDirVect = acuteAngle(MeanDirTemp*ones(1,length(DirVect)),DirVect) ; % acute angle between mean and all directions
        
        for ad = 1:length(AngleDiffDistX)-1 ; % for each bin
            temp = deltaDirVect>=AngleDiffDistX(ad) & deltaDirVect<AngleDiffDistX(ad+1) ; % identify times 
            AngleDiffDist(ad) = AngleDiffDist(ad)+temp(1:length(psth{DsType}(c,mi:end)))*psth{DsType}(c,mi:end)' ; % wieghted sum by spike rate
        end
    end
end
AngleDiffDist = AngleDiffDist/sum(AngleDiffDist) ;


%% peak to trough ratio and half width of Sta

% peak to trough ratio 
for DsType=1:length(DsTypeName) ; % for each DS cell type      
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        [~,PkPnt] = max(abs(StaXJitterVector{DsType}(c,StaWindowPnts))) ;
        StaXMax = max(StaXJitterVector{DsType}(c,StaWindowPnts(PkPnt))) ; % peak
        StaXMin = min(StaXJitterVector{DsType}(c,StaWindowPnts(PkPnt):StaWindowPnts(end))) ; % trough
        StaXMaxOverMin{DsType}(c) = abs(StaXMax/StaXMin) ; % peak/trough
        
        [~,PkPnt] = max(abs(StaYJitterVector{DsType}(c,StaWindowPnts))) ;
        StaYMax = max(StaYJitterVector{DsType}(c,StaWindowPnts(PkPnt))) ; % peak
        StaYMin = min(StaYJitterVector{DsType}(c,StaWindowPnts(PkPnt):StaWindowPnts(end))) ; % trough
        StaYMaxOverMin{DsType}(c) = abs(StaYMax/StaYMin) ; % peak/trough
    end
end

MaxOverMinMat = [cell2mat(StaXMaxOverMin),cell2mat(StaYMaxOverMin)] ;
StaMaxOverMin_mean = mean(MaxOverMinMat) ;
StaMaxOverMin_sem = std(MaxOverMinMat)/sqrt(length(MaxOverMinMat)) ;
