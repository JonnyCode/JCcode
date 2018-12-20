function ForIgor = ImJitterAnalysisV4_decode(DataBlock, DB, Params)

% edited from "ImJitterAnalysisV4" (NOTE- V4 was never fully tested) 
% moved sta code to ImJitterAnalysisV4_sta, added BinLoop, and export
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
StimBinFlag = false ; % if true, bin the stim same as psthbin number

BinLoopSet = [1,5,10,20] ; % number of frame bins in each bin

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

%% find quads and calculate quad response histograms 

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

%% decode
for BinLoop = 1:length(BinLoopSet) ;
    
    PsthBinNum = BinLoopSet(BinLoop) ; % number of bins for psth
    StimBinNum = PsthBinNum ; % number of frame bins in Stim to estimate

    % bin stimulus (same as psth - bin is average of X bins preceding)
    if StimBinFlag ; % if there is binning
        frame_timeBins = 0 ;
        for tb = 1:length(XJitterVectorDiff) ; % for each stim point
            strtPnt = max([1,tb-StimBinNum+1]) ; %
            XJitterVectorDiff(tb) = mean(XJitterVectorDiff(strtPnt:tb)) ;
            YJitterVectorDiff(tb) = mean(YJitterVectorDiff(strtPnt:tb)) ;
        end
    end

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

    %% OLE (optimal linear estimator)
    UseQuadFlag = false ; % true = use single quad only
    BestCellTestFlag = false ; % use only the strongest responding cell at any time point in estimate
    RandomQuadFlag = false ; % construct a random quadruplet

    OptLinDecoderFlag = false ; % run optimal linear decoder code instead of OLE code
    FilterLength = 20 ; % number of points to be used in optimal linear decoding

    if UseQuadFlag ; % if you want to use quads
        QuadLoopSet = [1:length(RealQuad_i)] ;
    else
        QuadLoopSet = 1 ; % use the entire population
    end

    for QuadLoop = QuadLoopSet ;

        clear OleCellsi ; % in case its been run before
        if UseQuadFlag % use a quad

            if ~RandomQuadFlag  % if you want to use a real quad
                for cells=1:4  % for each cell
                    OleCellsi(cells) = find(TypeCellMat(:,1)==quad{RealQuad_i(QuadLoop)}(cells,1),1)-1+...
                        quad{RealQuad_i(QuadLoop)}(cells,2) ;

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
            EstError_med(QuadLoop,:) = median(EstError) ;
        else   
            for td = 0:FilterLength ; % for each time delay
                EstError_med(QuadLoop,td+1) = median(EstError{td+1}) ;
            end
        end

        % transfer function of estimate and actual x,y
        [~,Besti] = min(EstError_med(QuadLoop,td+1)) ;
        EstX = -cosd(Est{Besti}([TestPnts(1)-Besti+1:(TestPnts(end)-FilterLength)])) ; % best estimate X
        EstY = -sind(Est{Besti}([TestPnts(1)-Besti+1:(TestPnts(end)-FilterLength)])) ; % best estimate Y
        XestTransfer = xcov(EstX,XJitterVectorDiff([TestPnts(1):(TestPnts(end)-FilterLength+Besti-1)]),'coef') ;
        YestTransfer = xcov(EstY,YJitterVectorDiff([TestPnts(1):(TestPnts(end)-FilterLength+Besti-1)]),'coef') ;

        % number of cells responding at each time point
        NumCellsResp = [sum(psthMat(Besti:end,OleCellsi)>0,2);zeros(1,Besti-1)'] ;
        NumCellsResp_DistX = [0:4] ;

        % distribution of number of cells
        NumCellsResp_Dist(QuadLoop,:) = hist(NumCellsResp,NumCellsResp_DistX)/length(NumCellsResp) ;

        % median error for each number
        for nc = 1:length(NumCellsResp_DistX) ; 
            Tempi = NumCellsResp([TestPnts(1):(TestPnts(end)-FilterLength)])==NumCellsResp_DistX(nc) ; 
            EstErrorForCellNumber(QuadLoop,nc) = median(EstError{Besti}(Tempi)) ;
        end

    end ; %QuadLoop
end ; % BinLoop
    % 


%% figures
StaTime = StaPnts*.025 ; % change to frame rate

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
subplot(1,3,3)
polar(cell2mat(StaDirAtPeakJitterVector),cell2mat(StaMagPeakJitterVector),'*') ; % to scale plot correctly
hold on 

for DsType=1:length(DsTypeName) ; % for each DS cell type
    for c = 1:length(cell_i{DsType}) ; % for each DS cell
        subplot(1,3,1)
        plot(StaTime(StaWindowPnts),StaMagJitterVector{DsType}(c,StaWindowPnts),Color_list{DsType})
        hold on
        xlabel('time (s)')
        ylabel('magnitude')

        subplot(1,3,2)
        plot(StaTime(StaWindowPnts),StaDirJitterVector{DsType}(c,StaWindowPnts)*180/pi,Color_list{DsType})
        hold on
        xlabel('time (s)')
        ylabel('direction')

        subplot(1,3,3)
        polar(StaDirAtPeakJitterVector{DsType}(c),StaMagPeakJitterVector{DsType}(c),[Color_list{DsType},'*'])
        hold on
    end
end

%saveas(gcf,'StaDb19Fig')
%print(gcf, '-dpdf','StaDb19Fig')
 
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
AngleDiffDistXcenters = interp1([1:length(AngleDiffDistX)],AngleDiffDistX,[1:length(AngleDiffDistX)]+.5,'linear','extrap') ;

bar(AngleDiffDistXcenters,AngleDiffDist/sum(AngleDiffDist),'k')
hold on
plot([1,1]*(AngleDiffDistXcenters*AngleDiffDist')/sum(AngleDiffDist),[0,.1],'k:')
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


