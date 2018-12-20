function ForIgor = ImJitterAnalysis(DataBlock, DB, Params)


% this function will analyze DS cells during a jittered image.

% JC 1/3/2017 

LoadOldMatFileFlag = false; % true skip calculations and just load
PlotFigsFlag = true ;
saveFigPath = ['/Users/jcafaro/Documents/AnalysisFigures/'] ;
numElectrodeLayers = 2 ; % number of electrode layers surrounding the max that the center is calculated over
eiDistBin = 100 ; % (um) bins 
eiDistMax = 1000 ; % (um) max distance ei can be separated

% load vaiables if saved
% if LoadOldMatFileFlag ;
%     load([saveFigPath,'entire'])
%     PlotFigsFlag = true ;
% else

ImNum = 4 ; %TEMP -SHOULD LOOP THIS EVENTUALLY!

% parameters
Color_list = {'c','r','b','g','k','y','c','r','b','g','k','y','c','r','b','g','k','y',...
    'c','r','b','g','k','y','c','r','b','g','k','y','c','r','b','g','k','y'} ; % order of colors for each 
FrameRate = 60.35 ; % rate of monitor display
FramesPerTrig = 100 ; % number of frames between each trigger
StaTime = [-3:1/FrameRate:3] ;% (sec) time before spike:time after spike
saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/DsIdsDb',num2str(DB),'ImNum',num2str(ImNum)] ;


% load data
dataRun = load_data(DataBlock(DB).ImJitterConcat{ImNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ;

% identify DS cells
dataRunTemp = load_data(DataBlock(DB).DsPath) ;
dataRunTemp = load_neurons(dataRunTemp) ;

Params.TimeBounds = [0,dataRunTemp.duration] ; % times of Dg data
Params.SkipBwMap = true ;
DataBlock(DB).DsConcatPath = DataBlock(DB).ImJitterConcat{ImNum} ;
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
slashi = strfind(DataBlock(DB).ImJitterConcat{ImNum},'/') ; % find the /
dashi = strfind(DataBlock(DB).ImJitterConcat{ImNum},'-') ; % find the -
StimPath = [DataBlock(DB).ImJitterConcat{ImNum}(1:slashi(end-1)),'stimuli/s0',DataBlock(DB).ImJitterConcat{ImNum}(dashi(end)+1:end)] ;
load(StimPath) ;

StaPnts = floor(StaTime*(FrameRate/stimulus.frame_interval)) ; 
StaTime = StaPnts/(FrameRate/stimulus.frame_interval) ; % this is the precise time (entered above an estimate)

Trigs_orig = dataRun.triggers(dataRun.triggers>Params.TimeBounds(2)) ; % Triggers in image jitter block
Trigs = Trigs_orig - Trigs_orig(1) ; % Triggers adjusted for first frame

% sta
if ~ischar(stimulus.XImageJitterVector) && stimulus.image_jitter_std>0 ; % if there was a jittered image
    
    XImageJitterVectorDiff = [0,diff(stimulus.XImageJitterVector)'] ;
    YImageJitterVectorDiff = [0,diff(stimulus.YImageJitterVector)'] ;
    
    for DsType=1:length(DsTypeName) ; % for each DS cell type
        for c = 1:length(cell_i{DsType}) ; % for each DS cell

            spikeTimes_orig = dataRun.spikes{cell_i{DsType}(c)}(dataRun.spikes{cell_i{DsType}(c)}>=Trigs_orig(1)) ; % spikes during stimulus
            spikeTimes = spikeTimes_orig - Trigs_orig(1) ; % spike times relative to trial start

            numSpikes{DsType}(c)= 0 ;
            StaXImageJitterVector{DsType}(c,:) = zeros(1,length(StaPnts)) ;
            StaYImageJitterVector{DsType}(c,:) = zeros(1,length(StaPnts)) ;

            for s=1:length(spikeTimes) ; % for each spike
                if spikeTimes(s)<Trigs(end) ; % if the spike time is not past the last trigger
                    T1i = find(Trigs<spikeTimes(s),1,'last') ; % index of trigger preceding spike        
                    FrameRateEstimate = (Trigs(T1i+1) - Trigs(T1i))/FramesPerTrig ; % estimated time of each image frame between triggers
                    FramesAfterTrig = floor((spikeTimes(s) - Trigs(T1i))/FrameRateEstimate) ; % number of frames after the trigger
                    FramesBeforeTrig = (T1i-1)*FramesPerTrig ;
                    spikePnt = FramesBeforeTrig + FramesAfterTrig ;
                    VectorPnts = spikePnt+StaPnts ;
                    if min(VectorPnts)>0 ; % if its late enough
                        StaXImageJitterVector{DsType}(c,:) = StaXImageJitterVector{DsType}(c,:) + XImageJitterVectorDiff(VectorPnts) ;
                        StaYImageJitterVector{DsType}(c,:) = StaYImageJitterVector{DsType}(c,:) + YImageJitterVectorDiff(VectorPnts) ;
                        numSpikes{DsType}(c) = numSpikes{DsType}(c) + 1 ;
                    end
                end
            end
            StaXImageJitterVector{DsType}(c,:) = StaXImageJitterVector{DsType}(c,:)/numSpikes{DsType}(c) ;
            StaYImageJitterVector{DsType}(c,:) = StaYImageJitterVector{DsType}(c,:)/numSpikes{DsType}(c) ;
            
            StaMagImageJitterVector{DsType}(c,:) = sqrt(StaXImageJitterVector{DsType}(c,:).^2+...
                StaYImageJitterVector{DsType}(c,:).^2) ; % magnitude of movement
            
            StaDirImageJitterVector{DsType}(c,:) = cart2pol(StaXImageJitterVector{DsType}(c,:),...
                StaYImageJitterVector{DsType}(c,:)) ; % direction of movement (radians)              
            
            [m,mi] = max(StaMagImageJitterVector{DsType}(c,:)) ;
            StaMagPeakImageJitterVector{DsType}(c) = m ;
            StaDirAtPeakImageJitterVector{DsType}(c) = StaDirImageJitterVector{DsType}(c,mi) ; % direction at peak of movement mag           
        end
    end
end
                

if ~ischar(stimulus.XSquareJitterVector) && stimulus.square_jitter_std>0 ; % if there was a jittered square
    
    XSquareJitterVectorDiff = [0,diff(stimulus.XSquareJitterVector)'] ;
    YSquareJitterVectorDiff = [0,diff(stimulus.YSquareJitterVector)'] ;
    
    for DsType=1:length(DsTypeName) ; % for each DS cell type
        for c = 1:length(cell_i{DsType}) ; % for each DS cell

            spikeTimes_orig = dataRun.spikes{cell_i{DsType}(c)}(dataRun.spikes{cell_i{DsType}(c)}>=Trigs_orig(1)) ; % spikes during stimulus
            spikeTimes = spikeTimes_orig - Trigs_orig(1) ; % spike times relative to trial start

            numSpikes{DsType}(c)= 0 ;
            StaXSquareJitterVector{DsType}(c,:) = zeros(1,length(StaPnts)) ;
            StaYSquareJitterVector{DsType}(c,:) = zeros(1,length(StaPnts)) ;

            for s=1:length(spikeTimes) ; % for each spike
                if spikeTimes(s)<Trigs(end) ; % if the spike time is not past the last trigger
                    T1i = find(Trigs<spikeTimes(s),1,'last') ; % index of trigger preceding spike         
                    FrameRateEstimate = (Trigs(T1i+1) - Trigs(T1i))/FramesPerTrig ; % estimated time of each image frame between triggers
                    FramesAfterTrig = floor((spikeTimes(s) - Trigs(T1i))/FrameRateEstimate) ; % number of frames after the trigger
                    FramesBeforeTrig = (T1i-1)*FramesPerTrig ;
                    spikePnt = FramesBeforeTrig + FramesAfterTrig ;
                    VectorPnts = spikePnt+StaPnts ;
                    if min(VectorPnts)>0 ; % if its late enough
                        StaXSquareJitterVector{DsType}(c,:) = StaXSquareJitterVector{DsType}(c,:) + XSquareJitterVectorDiff(VectorPnts) ;
                        StaYSquareJitterVector{DsType}(c,:) = StaYSquareJitterVector{DsType}(c,:) + YSquareJitterVectorDiff(VectorPnts) ;
                        numSpikes{DsType}(c) = numSpikes{DsType}(c) + 1 ;
                    end
                end
            end
            StaXSquareJitterVector{DsType}(c,:) = StaXSquareJitterVector{DsType}(c,:)/numSpikes{DsType}(c) ;
            StaYSquareJitterVector{DsType}(c,:) = StaYSquareJitterVector{DsType}(c,:)/numSpikes{DsType}(c) ;
            
            StaMagSquareJitterVector{DsType}(c,:) = sqrt(StaXSquareJitterVector{DsType}(c,:).^2+...
                StaYSquareJitterVector{DsType}(c,:).^2) ; % magnitude of movement
            
            StaDirSquareJitterVector{DsType}(c,:) = cart2pol(StaXSquareJitterVector{DsType}(c,:),...
                StaYSquareJitterVector{DsType}(c,:)) ; % direction of movement (radians)
            
            [m,mi] = max(StaMagSquareJitterVector{DsType}(c,:)) ;
            StaMagPeakSquareJitterVector{DsType}(c) = m ;
            StaDirAtPeakSquareJitterVector{DsType}(c) = StaDirSquareJitterVector{DsType}(c,mi) ; % direction at peak of movement mag
        end
    end
end                   

% weighted directions by the STA
DirectionWieghts = StaMagImageJitterVector{1}(1,(StaTime>-.3 & StaTime<=0)) ;

% linear prediction
% find psth time bins
psth_timeBins = 0 ;
for t = 1:length(Trigs)-1 ; % for each trigger
    FrameRateEstimate = (Trigs(t+1) - Trigs(t))/FramesPerTrig ; % estimated time of each image frame between triggers
    psth_timeBins = [psth_timeBins,[Trigs(t)+FrameRateEstimate:FrameRateEstimate:Trigs(t+1)]] ; % bins do not go beyond last trigger (though data may)
end

NumStimPnts = length(psth_timeBins) ; % number of stim frames actually delivered

if ~ischar(stimulus.XImageJitterVector) && stimulus.image_jitter_std>0 ; % if there was a jittered image
    
    XImageJitterVectorDiff = [0,diff(stimulus.XImageJitterVector)'] ;
    YImageJitterVectorDiff = [0,diff(stimulus.YImageJitterVector)'] ;
    
    MagImageJitterVector = sqrt(XImageJitterVectorDiff.^2 + YImageJitterVectorDiff.^2) ; % magnitude
    DirImageJitterVector = cart2pol(XImageJitterVectorDiff,YImageJitterVectorDiff) ; % direction of movement (radians)
    
    for DsType=1:length(DsTypeName) ; % for each DS cell type
        
        xF = mean(StaXImageJitterVector{DsType}(:,(StaTime>-.3 & StaTime<=0)),1) ; % x filter
        yF = mean(StaYImageJitterVector{DsType}(:,(StaTime>-.3 & StaTime<=0)),1) ; % y filter
        magF = sqrt(xF.^2+yF.^2) ; % filter mag
        
        for s=1:NumStimPnts ; % for each stimulus bin
            startPnt = max(1,s-length(DirectionWieghts)+1) ;
            endOffset = min(s-1,length(DirectionWieghts)-1) ;
            
            lp(DsType,s) = sum(XImageJitterVectorDiff(startPnt:s).*xF(end-endOffset:end) +...
                YImageJitterVectorDiff(startPnt:s).*xF(end-endOffset:end)) ;
        end
        
%         xF = fliplr(mean(StaXImageJitterVector{DsType},1)) ; % x filter
%         yF = fliplr(mean(StaYImageJitterVector{DsType},1)) ; % y filter
%         magF = sqrt(xF.^2+yF.^2) ; % filter mag
%         
%         XImageJitter_lp = conv(XImageJitterVector,xF,'same') ; % linear prediction of X
%         YImageJitter_lp = conv(YImageJitterVector,yF,'same') ; % linear prediction of Y
%         lp(DsType,:) = XImageJitter_lp + YImageJitter_lp  ; % linear prediction (dot product)
        
        psth_acrossCells(DsType,:) = zeros(1,length(psth_timeBins)) ; % prep vector
        for c = 1:length(cell_i{DsType}) ; % for each DS cell

            spikeTimes_orig = dataRun.spikes{cell_i{DsType}(c)}(dataRun.spikes{cell_i{DsType}(c)}>=Trigs_orig(1)) ; % spikes during stimulus
            spikeTimes = spikeTimes_orig - Trigs_orig(1) ; % spike times relative to trial start

            psth{DsType}(c,:) = histc(spikeTimes,psth_timeBins)' ;
            psth_acrossCells(DsType,:) = psth_acrossCells(DsType,:)+psth{DsType}(c,:) ;
        end
    end
end
    

% Correlations
% get ei center of mass and make TypeCellMat so can figure out which center belongs to which cell
CellCnt = 1 ;
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
        ctr(CellCnt,:) = get_ei_com(dataRun, cell_id{DsType}(cells), numElectrodeLayers) ;
        
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
        [s,si(CellCnt,:)] = sort(distCtr{CellCnt}) ; % in order of nearest to farthest

        CellCnt = CellCnt+1 ;
    end
end

% correlation mat
numCells = size(TypeCellMat,1) ; % number of all ds cells
corrMat_max = diag(ones(1,numCells)) ; % diagonal is 1 otherwise 0
corrMat_zero = diag(ones(1,numCells)) ; % diagonal is 1 otherwise 0
for cell1 = 1:numCells-1 ; % for each ds cell
    for cell2 = cell1+1:numCells ; % for each ds cell that hasn't been looked at
        TempCorr = xcov(psth{TypeCellMat(cell1,1)}(TypeCellMat(cell1,2),:),...
            psth{TypeCellMat(cell2,1)}(TypeCellMat(cell2,2),:),'coef') ;
        
        TempCorr_max = max(TempCorr) ; 
        corrMat_max(cell1,cell2) = TempCorr_max ;
        corrMat_max(cell2,cell1) = TempCorr_max ;
        
        TempCorr_zero = TempCorr(ceil(length(TempCorr)/2)) ; % correlation at time shift=0
        corrMat_zero(cell1,cell2) = TempCorr_zero ;
        corrMat_zero(cell2,cell1) = TempCorr_zero ;
    end
end
      

% correlation by distance
DistBins = [0,1:eiDistBin:eiDistMax] ;
DistMat = cell2mat(distCtr) ;

CorrByDistMean = cell(1,length(DsTypeName)) ; % prep array
for DsType=1:length(DsTypeName) ; % for each type
    CorrByDistMean{DsType} = nan(length(DsTypeName),length(DistBins)-1) ; % prep mats
end
CorrByDistStd = CorrByDistMean ;

for DsType1=1:length(DsTypeName) ; % for each type
    i1 = find(TypeCellMat(:,1)==DsType1) ; % cells of type A
    for DsType2=1:length(DsTypeName) ; % for each type
        corrVect = [] ; % make empty vector for exponential fit
        i2 = find(TypeCellMat(:,1)==DsType2) ; % cells of type B 
        for cells1=1:length(i1) ; % for each cell of type A
            for cells2=1:length(i2) ; % for each cell of type A
                corrVect = [corrVect;[DistMat(i1(cells1),i2(cells2)),corrMat_zero(i1(cells1),i2(cells2))]] ;
            end
        end
        % fit an exponential
        nonani = ~isnan(corrVect(:,2)) ;
        StatSetOptions = statset('FunValCheck','off') ; % creates a structure that prevents nlinfit from crashing with Nans
        fitCoef=nlinfit(corrVect(nonani,1),corrVect(nonani,2),@exponential,[.2,.5,100],StatSetOptions) ; % fit an exponential
        CorrExpFit{DsType1}(DsType2,:) = exponential(fitCoef,DistBins) ;
        
        % average within bins
        for b = 1:length(DistBins)-1 ;
            DistBinCenters(b) = mean([DistBins(b),DistBins(b+1)]) ;
            bi = find(DistBins(b)<=corrVect(:,1) & corrVect(:,1)<DistBins(b+1)) ;
            CorrByDistMean{DsType1}(DsType2,b) = nanmean(corrVect(bi,2)) ; 
            CorrByDistStd{DsType1}(DsType2,b) = nanstd(corrVect(bi,2)) ; 
        end
    end
end

% synchrony psth
for DsType = 1:length(DsTypeName) ; % for each ds type
    psth_sync(DsType,:) = sum(psth{DsType}>0,1)/size(psth{DsType},1) ; % fraction of cell firing at the same time
end

% stim correlation with high/low synchoronous events
for DsType = 1:length(DsTypeName) ; % for each ds type
    stimCorr_meanPsth(DsType,:) = xcov(mean(psth{DsType},1),MagImageJitterVector(1:length(psth{DsType})),'coef') ;
    
    stimCorr_SyncFull(DsType,:) = xcov(psth_sync(DsType,:),MagImageJitterVector(1:length(psth{DsType})),'coef') ;
    syncValues{DsType} = unique(psth_sync(DsType,:)) ; % all increments of syncrony
    for sv = 1:length(syncValues{DsType}) ; % for each value of synchrony
        stimCorr_SyncSelect{DsType}(sv,:) = xcov(psth_sync(DsType,:)==syncValues{DsType}(sv),MagImageJitterVector(1:length(psth{DsType})),'coef') ;
    end   
end

% OLE (optimal linear estimator)
% optimal weights = pinv(responses)*directions
% precision = 1/sum(optimal weights*responses - directions)^2
psthMat = cell2mat(psth')' ; % mat of responses cell=column, time=row
for s = 1:size(psthMat,1) ; 
    if sum(psthMat(s,:))~=0 ;
        psthMat(s,:) = psthMat(s,:)/sum(psthMat(s,:)) ; 
    end
end

DirectionWieghts = StaMagImageJitterVector{1}(1,(StaTime>-.3 & StaTime<=0)) ;
for s=1:NumStimPnts ; % for each stim bin
    startPnt = max(1,s-length(DirectionWieghts)+1) ;
    DirImageJitterVector_wieghted(s) = circ_mean(DirImageJitterVector(startPnt:s)',...
        DirectionWieghts(end-min(s-1,length(DirectionWieghts)-1):end)') ; % % direction as filtered through sta
end
Wopt = pinv(psthMat)*DirImageJitterVector_wieghted(1:length(psthMat))' ; % optimal weights
Ole = psthMat*Wopt ; % optimal linear estimate

% MLE (max likelihood estimate)
% response histograms
StimDir = DirImageJitterVector_wieghted(1:length(psthMat)) ;
StimDirBlocks = [-pi:.3:pi];
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
        for s=1:length(StimDirBlocks)-1 ; % for each stimulus direcition block 
            Stim_i= find(StimDir>StimDirBlocks(s) & StimDir<=StimDirBlocks(s+1)) ;
            NumStim(s) = length(Stim_i) ;
            DirResponse_mean{DsType}(cells,s) = mean(psth{DsType}(cells,Stim_i)) ;
            DirResponse_std{DsType}(cells,s) = std(psth{DsType}(cells,Stim_i)) ;
        end
    end
end
            

% figures

figure % spike rasters 
RasterDuration = 100 ; % (sec)
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
   
figure % sta for image
for DsType=1:length(DsTypeName) ; % for each DS cell type
        for c = 1:length(cell_i{DsType}) ; % for each DS cell 
            subplot(3,1,1)
            plot(StaTime,StaXImageJitterVector{DsType}(c,:),'r')
            hold on
            plot(StaTime,StaYImageJitterVector{DsType}(c,:),'b')
            hold off
            title(num2str(numSpikes{DsType}(c)))
            ylabel('delta X, Y')
            
            subplot(3,1,2)
            plot(StaTime,StaMagImageJitterVector{DsType}(c,:))
            ylabel('delta magnitude')
            
            subplot(3,1,3)
            plot(StaTime,StaDirImageJitterVector{DsType}(c,:))
            ylabel('delta direction')
            xlabel('time to spike (sec)')
            
            pause
        end
end
     

figure % sta for square
for DsType=1:length(DsTypeName) ; % for each DS cell type
    for c = 1:length(cell_i{DsType}) ; % for each DS cell 
        subplot(3,1,1)
        plot(StaTime,StaXSquareJitterVector{DsType}(c,:),'r')
        hold on
        plot(StaTime,StaYSquareJitterVector{DsType}(c,:),'b')
        hold off
        title(num2str(numSpikes{DsType}(c)))

        subplot(3,1,2)
        plot(StaTime,StaMagSquareJitterVector{DsType}(c,:))

        subplot(3,1,3)
        plot(StaTime,StaDirSquareJitterVector{DsType}(c,:))

        pause
    end
end
                
     
figure % all cells sta for image
subplot(1,3,3)
polar(cell2mat(StaDirAtPeakImageJitterVector),cell2mat(StaMagPeakImageJitterVector),'*') ; % to scale plot correctly
hold on 

for DsType=1:length(DsTypeName) ; % for each DS cell type
        for c = 1:length(cell_i{DsType}) ; % for each DS cell
            subplot(1,3,1)
            plot(StaTime,StaMagImageJitterVector{DsType}(c,:),Color_list{DsType})
            hold on
            xlabel('time (s)')
            ylabel('magnitude')
            
            subplot(1,3,2)
            plot(StaTime,StaDirImageJitterVector{DsType}(c,:),Color_list{DsType})
            hold on
            xlabel('time (s)')
            ylabel('direction')
            
            subplot(1,3,3)
            polar(StaDirAtPeakImageJitterVector{DsType}(c),StaMagPeakImageJitterVector{DsType}(c),[Color_list{DsType},'*'])
            hold on
        end
end
    
figure % all cells sta for square
subplot(1,3,3)
polar(cell2mat(StaDirAtPeakSquareJitterVector),cell2mat(StaMagPeakSquareJitterVector),'*') ; % to scale plot correctly
hold on 

for DsType=1:length(DsTypeName) ; % for each DS cell type
        for c = 1:length(cell_i{DsType}) ; % for each DS cell
            subplot(1,3,1)
            plot(StaTime,StaMagSquareJitterVector{DsType}(c,:),Color_list{DsType})
            hold on
            
            subplot(1,3,2)
            plot(StaTime,StaDirSquareJitterVector{DsType}(c,:),Color_list{DsType})
            hold on
            
            subplot(1,3,3)
            polar(StaDirAtPeakSquareJitterVector{DsType}(c),StaMagPeakSquareJitterVector{DsType}(c),[Color_list{DsType},'*'])
            hold on  
        end
end
       
% linear prediction
figure % linear prediction and mean psth
for DsType=1:length(DsTypeName) ; % for each DS cell type
    subplot(length(DsTypeName),1,DsType)
    plot(psth_acrossCells(DsType,:)/std(psth_acrossCells(DsType,:)),'b')
    hold on
    plot((lp(DsType,:)-mean(lp(DsType,:)))/std(lp(DsType,:)-mean(lp(DsType,:))),'r')
end
    
figure % cross correlation of lp and psth
for DsType=1:length(DsTypeName) ; % for each DS cell type
    subplot(length(DsTypeName),1,DsType)
    plot(xcov(psth_acrossCells(DsType,:),lp(DsType,1:length(psth_acrossCells)),'coef'))
end
    
figure % lp vs psth
for DsType=1:length(DsTypeName) ; % for each DS cell type
    subplot(length(DsTypeName),1,DsType)
    plot(lp(DsType,1:length(psth_acrossCells)),psth_acrossCells(DsType,:),'*')
    hold on
end

    
% correlations
figure
plot(DistMat(:),corrMat_zero(:),'b*')
hold on
plot(DistMat(:),corrMat_max(:),'ro')
xlabel('dist (um)')
ylabel('corr coef')

for DsType1=1:length(DsTypeName) ; % for each type
    i1 = find(TypeCellMat(:,1)==DsType1) ;
    figure
    for DsType2=1:length(DsTypeName) ; % for each type
        i2 = find(TypeCellMat(:,1)==DsType2) ;
        for cells1=1:length(i1) ;
            for cells2=1:length(i2) ;
                plot(DistMat(i1(cells1),i2(cells2)),corrMat_zero(i1(cells1),i2(cells2)),'*','Color',Color_list{DsType2})
                hold on
            end
        end
        %plot(DistBins,CorrExpFit{DsType1}(DsType2,:),'-','Color',Color_list{DsType2})
        errorbar(DistBinCenters,CorrByDistMean{DsType1}(DsType2,:),CorrByDistStd{DsType1}(DsType2,:),'-','Color',Color_list{DsType2})
    end
    title(num2str(DsType1))
end

% synchorony and stim correlations
figure
for DsType=1:length(DsTypeName) ; % for each type
     subplot(length(DsTypeName),1,DsType)
%     plot(stimCorr_meanPsth(DsType,length(stimCorr_SyncFull)/2:end),Color_list{1}) ;
%     hold on
%     plot(stimCorr_SyncFull(DsType,length(stimCorr_SyncFull)/2:end),Color_list{2}) ;
    for sv = 1:length(syncValues{DsType}) ; % for each value of synchrony
        plot(stimCorr_SyncSelect{DsType}(sv,length(stimCorr_SyncFull)/2:end),Color_list{2+sv})
        hold on
    end 
    xlabel('frame')
    ylabel('corr coef')
end

    
% tuning curves
figure
%polar(1,1)
%hold on
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
        %plot(StimDirBlocks(1:end-1),DirResponse_mean{DsType}(cells,:),Color_list{DsType})
        %polar(StimDirBlocks(1:end-1),DirResponse_mean{DsType}(cells,:),Color_list{DsType})
        errorbar(StimDirBlocks(1:end-1),DirResponse_mean{DsType}(cells,:),DirResponse_std{DsType}(cells,:),Color_list{DsType})
        hold on
        %pause
    end
end
xlabel('stim direction')
ylabel('response')

% look at prespike X,T trajectories
figure
plotPnts = find(StaTime>-0.14 & StaTime<-0.07) ;
DsType = 1 ; % for example ds type
c = 1 ; % for example cell 
spikeTimes_orig = dataRun.spikes{cell_i{DsType}(c)}(dataRun.spikes{cell_i{DsType}(c)}>=Trigs_orig(1)) ; % spikes during stimulus
spikeTimes = spikeTimes_orig - Trigs_orig(1) ; % spike times relative to trial start

for s=1:length(spikeTimes) ; % for each spike
    if spikeTimes(s)<Trigs(end) ; % if the spike time is not past the last trigger
        T1i = find(Trigs<spikeTimes(s),1,'last') ; % index of trigger preceding spike         
        FrameRateEstimate = (Trigs(T1i+1) - Trigs(T1i))/FramesPerTrig ; % estimated time of each image frame between triggers
        FramesAfterTrig = floor((spikeTimes(s) - Trigs(T1i))/FrameRateEstimate) ; % number of frames after the trigger
        FramesBeforeTrig = (T1i-1)*FramesPerTrig ;
        spikePnt = FramesBeforeTrig + FramesAfterTrig ;
        VectorPnts = spikePnt+StaPnts ;
        if min(VectorPnts)>0 ; % if its late enough
            plot(stimulus.XImageJitterVector(VectorPnts(plotPnts))-stimulus.XImageJitterVector(VectorPnts(plotPnts(1))),...
                stimulus.YImageJitterVector(VectorPnts(plotPnts))-stimulus.YImageJitterVector(VectorPnts(plotPnts(1)))) ;
            hold on
        end
    end
end

% compare dsi and sta mag amp
figure
for DsType=1:length(DsTypeName) ; % for each type
    plot(dsi{DsType}(:,1),StaMagPeakImageJitterVector{DsType},[Color_list{DsType},'*'])
    hold on
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


% for Joel

ForJoel(DB).stimulus = struct(stimulus) ;
ForJoel(DB).psth = psth ; % {DsType}(cells,:)
ForJoel(DB).psth_timeBins = psth_timeBins ; % {DsType}(cells,:)
ForJoel(DB).StaTime = StaTime ;
ForJoel(DB).StaMag = StaMagImageJitterVector ; % {DsType}(cells,:)
ForJoel(DB).StaDir = StaDirImageJitterVector ; % {DsType}(cells,:)
ForJoel(DB).WieghtedDirectionVector = DirImageJitterVector_wieghted ;     
for c = 1:size(TypeCellMat,1) ;
    ForJoel(DB).EiCenters{TypeCellMat(c,1)}(TypeCellMat(c,2),:) = ctr(c,:) ;
end


