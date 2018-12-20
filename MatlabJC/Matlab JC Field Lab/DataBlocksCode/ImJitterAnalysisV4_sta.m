function ForIgor = ImJitterAnalysisV4_sta(DataBlock, DB, Params)

% edited from "ImJitterAnalysisV4" (NOTE- V4 was never fully tested) 
% moved decoding code to ImJitterAnalysisV4_decode and added more exporting
% code

% JC 10/20/2018

% data set select
RunAsScript = true ;
if RunAsScript ;
    DB = 19 ; % 19,20,21
    [DataBlock,Params] = DataBlocks_NaturalStim ;
end

ImPathNum = 1 ; 
DsPathNum = 1 ;

% flags and parameters 
UseImportDsIdsPath = false ; % use ImportDsIdsPath to id DS cells
OnDsFlag = false ; % if true, analayze ON not ON-OFF DS cells
PsthBinNum = 1 ; % number of frame bins in each psth bin

QuadEiCircRad = 200 ; % minimum pairwise distance to be in quad
numElectrodeLayers = 2 ; % for EI com calculation

Color_list = {'c','r','b','g','k','y','c','r','b','g','k','y','c','r','b','g','k','y',...
    'c','r','b','g','k','y','c','r','b','g','k','y','c','r','b','g','k','y'} ; % order of colors for each 
FramesPerTrig = 100 ; % number of frames between each trigger

% save and loads paths
saveFigPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/2018Paper/ImJitterAnalysisV4_sta/Figs/'] ;
saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/ImJitterAnalysisV4DsSelection/DsIdsDb',num2str(DB),'ImJitterPathNum',num2str(ImPathNum)] ;
ImportDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/ImJitterAnalysisV4DsSelection/ImportDsIdsDb',num2str(DB),'ImJitterPathNum',num2str(ImPathNum)] ;
RunId = datestr(now,'yyyymmddhhMM') ; % RunId is date and time of code run

% load data
dataRun = load_data(DataBlock(DB).ImJitterConcat{ImPathNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun.piece.array_id = 1551 ; % CAUTION - assumes dense array
dataRun = load_ei(dataRun, 'all') ;

% identify DS cells
dataRunTemp = load_data(DataBlock(DB).DsPath{DsPathNum}) ;
dataRunTemp = load_neurons(dataRunTemp) ;

Params.TimeBounds = [0,dataRunTemp.duration] ; % times of Dg data
Params.OutlierFlag = true ; 
Params.DsPathNum = DsPathNum ;
DataBlock(DB).DsConcatPath = DataBlock(DB).ImJitterConcat{ImPathNum} ;
clear dataRunTemp

if UseImportDsIdsPath ; % if you want the DS ids selected elsewhere
    load(ImportDsIdsPath) ;
else
    try load(saveDsIdsPath) ; % if they are ds ids already saved
    catch % if not find them
        ForIgor = DsCellFinder(DataBlock, DB, Params) ;
        save(saveDsIdsPath, 'ForIgor')
    end
end

if OnDsFlag ; % if you want 'ON-DS' not 'ON-OFF DS' cells
    ci = 1 ; % On cells
else
    ci = 2 ; % On-Off cell 
end

cell_id = ForIgor.ds_id{ci} ; 
DsTypeName = ForIgor.dsName{ci} ; 

for DsType=1:length(DsTypeName) ; % cell indicies of DsTypes
    cell_i{DsType} = get_cell_indices(dataRun, cell_id{DsType}) ;
end

% load stimulus
slashi = strfind(DataBlock(DB).ImJitterConcat{ImPathNum},'/') ; % find the /
dashi = strfind(DataBlock(DB).ImJitterConcat{ImPathNum},'-') ; % find the -
StimPath = [DataBlock(DB).ImJitterConcat{ImPathNum}(1:slashi(end-1)),'stimuli/s0',DataBlock(DB).ImJitterConcat{ImPathNum}(dashi(end)+1)] ;
load(StimPath) ;

Trigs_orig = dataRun.triggers(dataRun.triggers>Params.TimeBounds(2)) ; % Triggers in image jitter block
Trigs = Trigs_orig - Trigs_orig(1) ; % Triggers adjusted for first frame

% stimulus velocity X,Y
if ~ischar(stimulus.XImageJitterVector) && stimulus.image_jitter_std>0 ; % if there was a jittered image
    XJitterVectorDiff = [0,diff(stimulus.XImageJitterVector)'] ;
    YJitterVectorDiff = [0,diff(stimulus.YImageJitterVector)'] ;
end

if ~ischar(stimulus.XSquareJitterVector) && stimulus.square_jitter_std>0 ; % if there was a jittered square
    XJitterVectorDiff = [0,diff(stimulus.XSquareJitterVector)'] ;
    YJitterVectorDiff = [0,diff(stimulus.YSquareJitterVector)'] ;
end
  
% find frame time 
frame_timeBins = 0 ;
for t = 1:length(Trigs)-1 ; % for each trigger
    FrameRateEstimate = (Trigs(t+1) - Trigs(t))/FramesPerTrig ; % estimated time of each image frame between triggers
    frame_timeBins = [frame_timeBins,[Trigs(t)+FrameRateEstimate:FrameRateEstimate:Trigs(t+1)]] ; % bins do not go beyond last trigger (though data may)
end
NumStimPnts = length(frame_timeBins) ; % number of stim frames actually delivered
FrameRateMean = mean(diff(frame_timeBins)) ; % avearge frame rate

% psth - at each frame
psthTimeBin = mean(diff(frame_timeBins))*PsthBinNum ; 

for DsType=1:length(DsTypeName) ; % for each DS cell type 
    psth{DsType} = nan(length(cell_i{DsType}),NumStimPnts) ; % prep matrix for speed
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        spk = dataRun.spikes{cell_i{DsType}(c)}(dataRun.spikes{cell_i{DsType}(c)}>=Trigs_orig(1)) ; % spikes during stimulus
        spk = spk - Trigs_orig(1) ; % spike times relative to trial start

        for tb = 1:NumStimPnts ; % for each frame bin
            strtPnt = max([1,tb-PsthBinNum]) ; 
            psth{DsType}(c,tb) = sum((spk>frame_timeBins(strtPnt)).*(spk<=frame_timeBins(tb)))/psthTimeBin ; % spikes per sec within window
        end    
    end
end
%% get ei center of mass 
for cells = 1:length(dataRun.spikes) ; % for each cell
    EiCtr(cells,:) = get_ei_com(dataRun, dataRun.cell_ids(cells), numElectrodeLayers) ;
end

%% sta and spike stats

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
for DsType=1:length(DsTypeName) ; % for each DS cell type
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        
        tempCorrX = xcov(psth{DsType}(c,:),XJitterVectorDiff) ; % 
        tempCorrY = xcov(psth{DsType}(c,:),YJitterVectorDiff) ; %
        
        StaXJitterVector{DsType}(c,:) = tempCorrX()/sum(psth{DsType}(c,:)) ;
        StaYJitterVector{DsType}(c,:) = tempCorrY()/sum(psth{DsType}(c,:)) ;

        StaMagJitterVector{DsType}(c,:) = sqrt(StaXJitterVector{DsType}(c,:).^2+...
            StaYJitterVector{DsType}(c,:).^2) ; % magnitude of movement

        StaDirJitterVector{DsType}(c,:) = cart2pol(StaXJitterVector{DsType}(c,:),...
            StaYJitterVector{DsType}(c,:)) ; % direction of movement (radians)              

        [m,mi] = max(StaMagJitterVector{DsType}(c,:)) ;
        StaMagPeakJitterVector{DsType}(c) = m ;
        StaDirAtPeakJitterVector{DsType}(c) = StaDirJitterVector{DsType}(c,mi) ; % direction at peak of movement mag 
    end
end

StaPnts = [1:length(tempCorrX)]-ceil(length(tempCorrX)/2) ; % X-axis of sta in points
StaWindowPnts = find(StaPnts>=0 & StaPnts<50) ; % window around sta peak for analysis
StaTime = StaPnts*FrameRateMean ; % change to frame rate

% average+/-Sem magnitude filter
StaMagJitterVector_mean = mean(cell2mat(StaMagJitterVector')) ; % mean across all cells
StaMagJitterVector_sem = std(cell2mat(StaMagJitterVector'))/sqrt(length(cell2mat(cell_id))) ; % mean across all cells

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

% half width
for DsType=1:length(DsTypeName) ; % for each DS cell type      
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        
        StaXAbs = abs(StaXJitterVector{DsType}(c,StaWindowPnts)) ;
        [StaXAbsMax,StaXAbsMaxi] = max(StaXAbs) ;
        p1 = find(StaXAbs(1:StaXAbsMaxi)>StaXAbsMax/2,1,'first') ; % first point above half max
        p2 = find(StaXAbs(StaXAbsMaxi:end)<StaXAbsMax/2,1,'first') + StaXAbsMaxi - 1; % first point below half max
        StaXhalfwidth{DsType}(c) = p2-p1 ;
        
        StaYAbs = abs(StaYJitterVector{DsType}(c,StaWindowPnts)) ;
        [StaYAbsMax,StaYAbsMaxi] = max(StaYAbs) ;
        p1 = find(StaYAbs(1:StaYAbsMaxi)>StaYAbsMax/2,1,'first') ; % first point above half max
        p2 = find(StaYAbs(StaYAbsMaxi:end)<StaYAbsMax/2,1,'first') + StaYAbsMaxi - 1; % first point below half max
        StaYhalfwidth{DsType}(c) = p2-p1 ;
    end
end

HalfwidthMat = [cell2mat(StaXhalfwidth),cell2mat(StaYhalfwidth)] ;
StaHalfwidth_mean = mean(HalfwidthMat) ;
StaHalfwidth_sem = std(HalfwidthMat)/sqrt(length(HalfwidthMat)) ;

%% response distributions (all cells using fixed bin psth for speed)
cell_i1Mat = get_cell_indices(dataRun, cell2mat(ForIgor.ds_id{1})) ; % selected DS cells 1
cell_i2Mat = get_cell_indices(dataRun, cell2mat(ForIgor.ds_id{2})) ; % selected DS cells 2
NonDsi = setdiff([1:length(dataRun.spikes)],[cell_i1Mat,cell_i2Mat]) ; % non-DS cells
RspDistX = [0:10] ; % number of spikes in a bin

for c=1:length(dataRun.spikes) ;
    spk = dataRun.spikes{c}(dataRun.spikes{c}>=Trigs_orig(1)) ; % spikes during stimulus
    spk = spk - Trigs_orig(1) ; % spike times relative to trial start
    
    psthTemp = hist(spk,frame_timeBins) ; % fixed bin psth
    RspDist(c,:) = hist(psthTemp,RspDistX) ; 
end

RspDist_meanNonDs = mean(RspDist(NonDsi,:)) ;
Ntemp = sum(RspDist_meanNonDs) ;
RspDist_meanNonDs = RspDist_meanNonDs/Ntemp ;
RspDist_semNonDs = std(RspDist(NonDsi,:))/sqrt(length(NonDsi)) ;
RspDist_semNonDs = RspDist_semNonDs/Ntemp ;

RspDist_meanDs1 = mean(RspDist(cell_i1Mat,:)) ;
Ntemp = sum(RspDist_meanDs1) ;
RspDist_meanDs1 = RspDist_meanDs1/Ntemp ;
RspDist_semDs1 = std(RspDist(cell_i1Mat,:))/sqrt(length(cell_i1Mat)) ;
RspDist_semDs1 = RspDist_semDs1/Ntemp ;

RspDist_meanDs2 = mean(RspDist(cell_i2Mat,:)) ;
Ntemp = sum(RspDist_meanDs2) ; 
RspDist_meanDs2 = RspDist_meanDs2/Ntemp ;
RspDist_semDs2 = std(RspDist(cell_i2Mat,:))/sqrt(length(cell_i2Mat)) ;
RspDist_semDs2 = RspDist_semDs2/Ntemp ;

%% figures

figure % spike rasters 
RasterDuration = 1 ; % (sec)
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
        plot(StaTime(StaWindowPnts),StaXJitterVector{DsType}(c,StaWindowPnts),'r')
        hold on
        plot(StaTime(StaWindowPnts),StaYJitterVector{DsType}(c,StaWindowPnts),'b')
        hold off
        %title(num2str(numSpikes{DsType}(c)))
        ylabel('delta X, Y')

        subplot(3,1,2)
        plot(StaTime(StaWindowPnts),StaMagJitterVector{DsType}(c,StaWindowPnts))
        ylabel('delta magnitude')

        subplot(3,1,3)
        plot(StaTime(StaWindowPnts),StaDirJitterVector{DsType}(c,StaWindowPnts))
        ylabel('delta direction')
        xlabel('time to spike (sec)')

        pause
    end
end
  
figure % all cells sta for image
% subplot(1,3,3)
% polar(cell2mat(StaDirAtPeakJitterVector),cell2mat(StaMagPeakJitterVector),'*') ; % to scale plot correctly
% hold on 

for DsType=1:length(DsTypeName) ; % for each DS cell type
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        subplot(1,2,1)
        plot(StaTime(StaWindowPnts),StaMagJitterVector{DsType}(c,StaWindowPnts),Color_list{DsType})
        hold on
        xlabel('time (s)')
        ylabel('magnitude')
        xlim([0 1])
        box off

        subplot(1,2,2)
        plot(StaTime(StaWindowPnts),StaDirJitterVector{DsType}(c,StaWindowPnts)*180/pi,Color_list{DsType})
        hold on
        xlabel('time (s)')
        ylabel('direction')
        xlim([0 1])
        ylim([-180 180])
        set(gca,'YTick',[-180:45:180])
        box off

%         subplot(1,3,3)
%         polar(StaDirAtPeakJitterVector{DsType}(c),StaMagPeakJitterVector{DsType}(c),[Color_list{DsType},'*'])
%         hold on
    end
end

saveas(gcf,[saveFigPath,'StaDb',num2str(DB),'Fig_',RunId])
print(gcf, '-dpdf',[saveFigPath,'StaDb',num2str(DB),'Fig_',RunId])
 
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

figure % direction distribution
AngleDiffDistXcenters = interp1([1:length(AngleDiffDistX)],AngleDiffDistX,[1:length(AngleDiffDistX)]+.5,'linear','extrap') ;

bar(AngleDiffDistXcenters,AngleDiffDist/sum(AngleDiffDist),'k')
hold on
plot([1,1]*(AngleDiffDistXcenters*AngleDiffDist')/sum(AngleDiffDist),[0,.1],'k:')
xlabel('Angle off mean (deg)')
ylabel('fraction')

hgsave([saveFigPath,'DirectionDistDb',num2str(DB),'_',RunId])
print(gcf, '-dpdf', [saveFigPath,'DirectionDistDb',num2str(DB),'_',RunId])

figure % mean sta mag vector
plot(StaTime(StaWindowPnts),StaMagJitterVector_mean(StaWindowPnts),'k')
hold on
plot(StaTime(StaWindowPnts),StaMagJitterVector_mean(StaWindowPnts)+StaMagJitterVector_sem(StaWindowPnts),'k:')
plot(StaTime(StaWindowPnts),StaMagJitterVector_mean(StaWindowPnts)-StaMagJitterVector_sem(StaWindowPnts),'k:')
xlim([0 1])

hgsave([saveFigPath,'StaMagJitterVectorMean',num2str(DB),'_',RunId])
print(gcf, '-dpdf', [saveFigPath,'StaMagJitterVectorMean',num2str(DB),'_',RunId])

figure ; % response distribtution
errorbar(RspDistX,RspDist_meanNonDs,RspDist_semNonDs,'g')
hold on
errorbar(RspDistX,RspDist_meanDs1,RspDist_semDs1,'b')
errorbar(RspDistX,RspDist_meanDs2,RspDist_semDs2,'r')
set(gca,'Yscale','log')
xlabel('spikes per bin')
ylabel('fraction')
legend('nonDs','OnDs','ooDS')

hgsave([saveFigPath,'RspDistDb',num2str(DB),'_',RunId])
print(gcf, '-dpdf', [saveFigPath,'RspDistDb',num2str(DB),'_',RunId])

%% save mat file analysis
save([saveFigPath,'ImJitterAnalysisV4_sta_DB',num2str(DB),'_',RunId],'-v7.3') ; % save mat file

