function ForIgor = RepStimNecSufAnalysisV2(DataBlock, DB, Params) 

% adapted from 'RepStimNecSufAnalyis' to use PopDistFinder2 instead of
% PopDistFinder

% JC 7/18/2017


mapEiFlag = false ;
SquareGridSizeUpperBound = 300 ;% (um) ~size of boxes on one side of a square grid dividing up array space
numElectrodeLayers = 2 ; % number of electrode layers used to calc ei center
RepStimNum = 1 ; % select repstim set
BwPathNum = 1 ;
numTrigPerRep = 3 ; % (integer) number of triggs per stimulus repeat
PsthBinTime = 0.01 ; % (sec) 
BinSearchNumber = 2 ; % (integer) number of psth bins that can be compared 
StimFrameNumber = 300 ; % (frames) number of frames in repeated stimulus
StimFrameStimNumber = 10 ; % (frames) number of frames that RGCs can integrate over
StimPath = '/Volumes/lab/Documents/Movies/CatCam/cat_mean117_sd62_0to255' ; 
eiDistBin = 100 ;
eiDistMax = 1000 ;

% load data
dataRun = load_data(DataBlock(DB).MovieRepPath{RepStimNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun,'all') ;

% load BW data with identified neurons
dataRunBw = load_data(DataBlock(DB).BwPath{BwPathNum}) ;
dataRunBw = load_neurons(dataRunBw) ;
dataRunBw = load_params(dataRunBw,'cell_type_depth', 5) ;
dataRunBw = load_sta(dataRunBw) ;

marks_params.thresh = 4.5;
dataRunBw= get_sta_summaries(dataRunBw, 'all','marks_params', marks_params);


% calculate ei centers
for c = 1:length(dataRun.spikes) ; % for each cell
    eiCnt(c,:) = get_ei_com(dataRun, dataRun.cell_ids(c), numElectrodeLayers) ;
end

% find groups of cells in each grid bin

SquareGridNum_X = ceil(range(dataRun.ei.array_bounds_x)/SquareGridSizeUpperBound) ; % bin number
SquareGridNum_Y = ceil(range(dataRun.ei.array_bounds_y)/SquareGridSizeUpperBound) ;

SquareGridSize_X = range(dataRun.ei.array_bounds_x)/SquareGridNum_X ; % bin size
SquareGridSize_Y = range(dataRun.ei.array_bounds_y)/SquareGridNum_Y ;

r=1 ;
for x=1:SquareGridNum_X ;
    for y = 1:SquareGridNum_Y ;
        xbounds = [(x-1)*SquareGridSize_X+dataRun.ei.array_bounds_x(1),...
            x*SquareGridSize_X+dataRun.ei.array_bounds_x(1)] 
        ybounds = [(y-1)*SquareGridSize_Y+dataRun.ei.array_bounds_y(1),...
            y*SquareGridSize_Y+dataRun.ei.array_bounds_y(1)] 
        
        xi = find(eiCnt(:,1)>=xbounds(1) & eiCnt(:,1)<xbounds(2)) ;
        yi = find(eiCnt(:,2)>=ybounds(1) & eiCnt(:,2)<ybounds(2)) ;
        Group_i{r} = intersect(xi,yi) ; % indicies of cells within the box
        r=r+1 ;
    end
end
        
% identified cell types

for a = 1:length(dataRunBw.cell_types) ;
    celltypes{a} = dataRunBw.cell_types{a}.name ;
end

UniqueCellTypes = unique(celltypes) ;
if isempty(UniqueCellTypes{1}) ;
    UniqueCellTypes = UniqueCellTypes(2:end) ;
end

% map dim flash onto bw data
if mapEiFlag ; % if using map ei cells

    % map using electrical images
    cell_list_map = map_ei(dataRunBw, dataRun) ;

    % cells ids in slave for each UniqueCellType set in master data
    for uc = 1:length(UniqueCellTypes) ;
        Masteri{uc} = get_cell_indices(dataRunBw, UniqueCellTypes{uc}) ;
        cell_ids{uc} = cell2mat(cell_list_map(Masteri{uc})) ;
    end
else % if not using map ei
    for uc = 1:length(UniqueCellTypes) ;
        cell_ids{uc} = intersect(dataRun.cell_ids, get_cell_ids(dataRunBw,UniqueCellTypes{uc})) ;
        cell_i{uc} = get_cell_indices(dataRun,cell_ids{uc}) ;
    end
end

% relationship between distance and spike rate
for c=1:length(dataRun.spikes) ; % for each cell
    PopDist = PopDistFinderV2(dataRun.spikes(c), dataRun.triggers([1:numTrigPerRep:length(dataRun.triggers)]),...
            'PsthBinTime', PsthBinTime) ;
    psth(c,:) = PopDist.r_mean ;
    Rdist(c,:) = PopDist.AcrossStim ;   
    
    tempCorr = corrcoef(psth(c,:),Rdist(c,:)) ;
    RateDistCorr(c) = tempCorr(2,1) ;
end

% relationship between correlations in firing rate vs correlations in
% distance
for c=1:length(dataRun.spikes) ; % for each cell
    for c2=1:length(dataRun.spikes) ; % for each cell
        tempCorr = corrcoef(psth(c,:),psth(c2,:)) ;
        psthCorr(c,c2) = tempCorr(2,1) ;
        tempCorr = corrcoef(Rdist(c,:),Rdist(c2,:)) ;
        RdistCorr(c,c2) = tempCorr(2,1) ;
    end
end

% cell type spike rate and rdist correlations

% correlation by distance
DistBins = [0,1:eiDistBin:eiDistMax] ;
NearbyThreshRangei = [2:4] ;

for uc = 1:length(UniqueCellTypes) ; % for each cell type
    for uc2 = 1:length(UniqueCellTypes) ; % for each cell type
        RdistCorr_ByDist_inbins{uc}{uc2} = cell(1,length(DistBins)) ;
        PsthCorr_ByDist_inbins{uc}{uc2} = cell(1,length(DistBins)) ;
        
        for c=1:length(cell_i{uc}) ; % for each cell of type 1  
            for c2=1:length(cell_i{uc2}) ; % for each cell of type 2
                distEiCnt = sqrt(sum((eiCnt(cell_i{uc}(c),:) - eiCnt(cell_i{uc2}(c2),:)).^2)) ;
                bin = find(DistBins<=distEiCnt,1,'last') ; % appropriate distance bin
                
                RdistCorr_ByDist_inbins{uc}{uc2}{bin} = [RdistCorr_ByDist_inbins{uc}{uc2}{bin}, RdistCorr(cell_i{uc}(c),cell_i{uc2}(c2))] ;
                PsthCorr_ByDist_inbins{uc}{uc2}{bin} = [PsthCorr_ByDist_inbins{uc}{uc2}{bin}, psthCorr(cell_i{uc}(c),cell_i{uc2}(c2))] ;
            end
        end
    end
end
     
% find range and average cells within same bins
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    for uc2 = 1:length(UniqueCellTypes) ; % for each cell type
        for bin=1:length(DistBins) ; % for each bin
              RdistCorr_ByDist{uc}(uc2,bin) = nanmean(RdistCorr_ByDist_inbins{uc}{uc2}{bin}) ;
              PsthCorr_ByDist{uc}(uc2,bin) = nanmean(PsthCorr_ByDist_inbins{uc}{uc2}{bin}) ;
              
              if length(RdistCorr_ByDist_inbins{uc}{uc2}{bin})>1 ;
                RdistCorr_ByDist_max{uc}(uc2,bin) = max(RdistCorr_ByDist_inbins{uc}{uc2}{bin}) ;
                PsthCorr_ByDist_max{uc}(uc2,bin) = max(PsthCorr_ByDist_inbins{uc}{uc2}{bin}) ;  
                  
                RdistCorr_ByDist_range{uc}(uc2,bin) = range(RdistCorr_ByDist_inbins{uc}{uc2}{bin}) ;
                PsthCorr_ByDist_range{uc}(uc2,bin) = range(PsthCorr_ByDist_inbins{uc}{uc2}{bin}) ;
              else
                RdistCorr_ByDist_max{uc}(uc2,bin) = nan ;
                PsthCorr_ByDist_max{uc}(uc2,bin) = nan ;
                  
                RdistCorr_ByDist_range{uc}(uc2,bin) = nan ;
                PsthCorr_ByDist_range{uc}(uc2,bin) = nan ;
              end          
        end
        RdistCorr_NearbyMat(uc,uc2) = nanmean(RdistCorr_ByDist{uc}(uc2, NearbyThreshRangei)) ;
        PsthCorr_NearbyMat(uc,uc2) = nanmean(PsthCorr_ByDist{uc}(uc2, NearbyThreshRangei)) ;
        
        RdistCorr_NearbyMat_range(uc,uc2) = nanmean(RdistCorr_ByDist_range{uc}(uc2, NearbyThreshRangei)) ;
        PsthCorr_NearbyMat_range(uc,uc2) = nanmean(PsthCorr_ByDist_range{uc}(uc2, NearbyThreshRangei)) ;
        
        RdistCorr_NearbyMat_max(uc,uc2) = nanmean(RdistCorr_ByDist_max{uc}(uc2, NearbyThreshRangei)) ;
        PsthCorr_NearbyMat_max(uc,uc2) = nanmean(PsthCorr_ByDist_max{uc}(uc2, NearbyThreshRangei)) ;
    end
end
            
            
    
            
        
        
        
    
    

% figures
PsthTime = [1:size(psth,2)]*PsthBinTime ;

figure % psth vs Rdist
for c=1:length(dataRun.spikes) ; % for each cell
    subplot(2,1,1)
    plotyy(PsthTime, psth(c,:),PsthTime,Rdist(c,:))
    xlabel('time')
    ylabel('spike rate , discriminability')
    
    
    subplot(2,1,2)
    plot(psth(c,:),Rdist(c,:),'*')
    xlabel('spike rate')
    ylabel('discriminability')
    pause
end



figure % correlation of psth compared to correlation of rdist
for c=1:length(dataRun.spikes) ; % for each cell
    for c2=1:length(dataRun.spikes) ; % for each cell
        plot(xcov(psth(c,:),psth(c2,:),'coef'))
        hold on
        plot(xcov(Rdist(c,:),Rdist(c2,:),'coef'))
        pause
        hold off
    end
end

figure % correlation of psth vs correlation of Rdist - all cell pairs
plot(psthCorr(:),RdistCorr(:),'*')
hold on
plot(psthCorr([1:length(dataRunBw.spikes)],[1:length(dataRunBw.spikes)]),RdistCorr([1:length(dataRunBw.spikes)],[1:length(dataRunBw.spikes)]),'ro')

figure % correlation of psth vs correlation of Rdist - averaged across types for nearby cells only
plot(PsthCorr_NearbyMat,RdistCorr_NearbyMat,'k*')
hold on
plot(diag(PsthCorr_NearbyMat),diag(RdistCorr_NearbyMat),'ro')
xlabel('spike rate correlation')
ylabel('discriminability correlation')

figure
imagesc(RdistCorr_NearbyMat)
colorbar

figure
imagesc(PsthCorr_NearbyMat)
colorbar

figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    for uc2 = 1:length(UniqueCellTypes) ; % for each cell type
        plot(PsthCorr_NearbyMat(uc,uc2),RdistCorr_NearbyMat(uc,uc2),'*')
        axis([-1,1,-1,1])
        title(['uc =', UniqueCellTypes{uc},' uc2=', UniqueCellTypes{uc2}])
        pause
    end
end

figure % range vs mean of correlations
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    for uc2 = 1:length(UniqueCellTypes) ; % for each cell type
        plot(RdistCorr_ByDist_max{uc}(uc2,:), RdistCorr_ByDist_range{uc}(uc2,:),'*')
        hold on
    end
end

figure
plot(RdistCorr_NearbyMat_max, RdistCorr_NearbyMat_max-RdistCorr_NearbyMat_range,'*')
xlabel('max correlation')
ylabel('min correlation')

figure
plot(PsthCorr_NearbyMat_max, PsthCorr_NearbyMat_range,'*')

figure % rdist by distance
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    for uc2 = 1:length(UniqueCellTypes) ; % for each cell type
        if ~isempty(RdistCorr_ByDist{uc}) ;
            subplot(1,2,1)
            plot(DistBins, RdistCorr_ByDist{uc}(uc2,:),'*') ;
            hold on
            plot(DistBins, PsthCorr_ByDist{uc}(uc2,:),'o') ;
            hold off
            title(['uc =', num2str(uc),' uc2=', num2str(uc2)])
            
            subplot(1,2,2)
            plot(RdistCorr_ByDist{uc}(uc2,:),PsthCorr_ByDist{uc}(uc2,:),'*') ;
            hold on
            
            pause
        end 
    end
end
 

figure % rdist by distance
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    for uc2 = 1:length(UniqueCellTypes) ; % for each cell type
        if ~isempty(RdistCorr_ByDist{uc}) ;
            subplot(length(UniqueCellTypes),length(UniqueCellTypes),(uc-1)*length(UniqueCellTypes)+uc2)
            plot(DistBins, RdistCorr_ByDist{uc}(uc2,:),'*') ;
            hold on
            plot(DistBins, PsthCorr_ByDist{uc}(uc2,:),'o') ;
        end 
    end
end
 
figure % Srf
norm_rf = norm_image(dataRunBw.stas.rfs{304});
imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
colormap(brewermap([],'RdBu'))
caxis([0,1]) 

figure % Trf
plot(dataRunBw.stas.time_courses{304})

