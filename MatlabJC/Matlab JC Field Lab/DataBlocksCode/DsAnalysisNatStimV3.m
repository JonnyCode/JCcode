function ForIgor = DsAnalysisNatStimV3(DataBlock, DB, Params) 

% this function will analyze DS cells during natural stimulus,
% adapted from 'DsAnalysisNatStimV2', to move figures logic to bottom of
% code, deal with multi-movies, save DS cell_ids, calculate spike rate statistics,
% calculate synchrony within a cell type, and vector sum in large groups

% JC 11/17/2016 

LoadOldMatFileFlag = true ; % true skip calculations and just load
PlotFigsFlag = false ;
saveFigPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/Db',num2str(DB),'Movie',num2str(Params.MovieNum)] ;

% load vaiables if saved
if LoadOldMatFileFlag ;
    load([saveFigPath,'entire'])
else

% parameters
Color_list = {'c','r','b','g','k','y','c','r','b','g','k','y','c','r','b','g','k','y',...
    'c','r','b','g','k','y','c','r','b','g','k','y','c','r','b','g','k','y'} ; % order of colors for each 
EiCircRad = 100 ; % (um) radius of possible dendrites around ei center of mass
eiDistBin = 50 ; % (um) bins 
eiDistMax = 1000 ; % (um) max distance ei can be separated
FramesPerTrig = 100 ;
FramesPerSecond = 60.35 ;
eiAxisLims = [-450, 450, -450, 450] ; % Axis limits figures ([-425, 425, -425, 425] covers entire array) ;
numElectrodeLayers = 2 ; % number of electrode layers surrounding the max that the center is calculated over
MaxNumCells = 10 ; % number of psths of the nearest cells that you want to plot
bin_size =  0.1 ; % (s) psth bin size (0.1 is good default)
EiSpikeRateFactor = 100 ; % (um/spike)
VectorSumDistMax = 300 ; % (um) distance at which center firing rates between different cell types works
synchronyDistMax = 300 ; % (um) distance at which synchrony between same cell type matters
SpikeToLineWidthFact = 0.1 ; % (norm spikes per change in line width)
TypeAngle = [120,30,300,210] ; % prefered angle for each DsType
TimeStep = .016 ; % (s) time steps of composite movie/psth composite (0.016 is good default)
psthMaxOverSum_distX = [0:.1:1] ; % distribution X
psth_histX = [0,1:20:500] ; % (hz) distribution X
saveDsIdsPath = ['/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/DsIdsDb',num2str(DB),'Movie',num2str(Params.MovieNum)] ;
DsPathNum = 1 ;

BinaryFlag = false ; % use binary psth to assess correlation ; 

MovieStixelWidth = DataBlock(DB).MovieStixWidth(Params.MovieNum) ; 
FrameInterval = DataBlock(DB).MovieFrameInterval(Params.MovieNum) ; % frame is switched every...

% load data
dataRun = load_data(DataBlock(DB).MovieDsConcat{Params.MovieNum}) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ;

% identify DS cells
dataRunTemp = load_data(DataBlock(DB).DsPath{DsPathNum}) ;
dataRunTemp = load_neurons(dataRunTemp) ;

Params.TimeBounds = [0,dataRunTemp.duration] ; % times of Dg data
Params.SkipBwMap = true ;
Params.DsPathNum = DsPathNum ;
DataBlock(DB).DsConcatPath = DataBlock(DB).MovieDsConcat{Params.MovieNum} ;
clear dataRunTemp

try load(saveDsIdsPath) ; % if they are ds ids already saved
catch % if not find them

    ForIgor = DsCellFinder(DataBlock, DB, Params) ;

    cell_id = ForIgor.ds_id{2} ; % On-Off cell 
    DsTypeName = ForIgor.dsName{2} ; 

    for DsType=1:length(DsTypeName) ;
        cell_i{DsType} = get_cell_indices(dataRun, cell_id{DsType}) ;
    end

    save(saveDsIdsPath, 'cell_id','DsTypeName','cell_i')
end

% calculate movie triggers
trigNum = ceil(DataBlock(DB).MovieRepFrameNum(Params.MovieNum)/FramesPerTrig) ; % number of triggers per repeat
FirstTrig = find(dataRun.triggers>=Params.TimeBounds(2),1,'first') ; % first trigger of the movie
trigs = dataRun.triggers(FirstTrig:trigNum:end) ; % begining of each movie 
trigEnds = dataRun.triggers(FirstTrig+trigNum-1:trigNum:end) ; % last trigger of each movie 
repeatNum= floor(length(dataRun.triggers(FirstTrig:end))/trigNum) ; % show how many repeats this is
stim_time =[0, DataBlock(DB).MovieRepFrameNum(Params.MovieNum)/FramesPerSecond] ;
                
% get psths for all ds cells
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
        if ~isnan(cell_i{DsType}(cells)) ;
            spikes = dataRun.spikes{cell_i{DsType}(cells)} ;
            [psth{DsType}(cells,:), bins] = get_psth(spikes(spikes>Params.TimeBounds(2)), trigs(trigs>Params.TimeBounds(2)),...
                'stop',stim_time(2),'bin_size',bin_size) ;
            %[psth{DsType}(cells,:),psthTime] = get_smooth_psth(spikes(spikes>Params.TimeBounds(2)), trigs(trigs>Params.TimeBounds(2)), 'stop', stim_time(2)) ;
        else
            psth{DsType}(cells,:) = nans(1,ceil(stim_time(2)/bin_size)) ;
        end
        psth_hist{DsType}(cells,:) = hist(psth{DsType}(cells,:),psth_histX) ;
        psth_norm{DsType}(cells,:) = psth{DsType}(cells,:)/max(psth{DsType}(cells,:)) ; % normalized by max
        psth_binary{DsType}(cells,:) = psth{DsType}(cells,:)>0 ; % make psth 1 if firing 0 if not
    end
end

% get ei center of mass and make TypeCellMat so can figure out which center
% belongs to which cell
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
        
        if BinaryFlag ;
            TempCorr = xcorr(double(psth_binary{TypeCellMat(cell1,1)}(TypeCellMat(cell1,2),:)),...
            double(psth_binary{TypeCellMat(cell2,1)}(TypeCellMat(cell2,2),:)),'unbiased') ; % CHECK NORMALIZATION - fraction of time they fire together
        end
        
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
%         % fit an exponential
%         nonani = ~isnan(corrVect(:,2)) ;
%         StatSetOptions = statset('FunValCheck','off') ; % creates a structure that prevents nlinfit from crashing with Nans
%         fitCoef=nlinfit(corrVect(nonani,1),corrVect(nonani,2),@exponential,[.2,.5,100],StatSetOptions) ; % fit an exponential
%         CorrExpFit{DsType1}(DsType2,:) = exponential(fitCoef,DistBins) ;
        
        % average within bins
        for b = 1:length(DistBins)-1 ;
            DistBinCenters(b) = mean([DistBins(b),DistBins(b+1)]) ;
            bi = find(DistBins(b)<=corrVect(:,1) & corrVect(:,1)<DistBins(b+1)) ;
            CorrByDistMean{DsType1}(DsType2,b) = nanmean(corrVect(bi,2)) ; 
            CorrByDistStd{DsType1}(DsType2,b) = nanstd(corrVect(bi,2)) ; 
        end
    end
end

% correlation by type combination
CorrHistX = [-1:.1:1] ;
TypeGroupCorr = cell(1,3) ;
TypeGroupDist = cell(1,3) ;

for DsType1=1:length(DsTypeName) ; % for each type
    i1 = find(TypeCellMat(:,1)==DsType1) ;
    for DsType2=1:length(DsTypeName) ; % for each type
        i2 = find(TypeCellMat(:,1)==DsType2) ;

        if abs(DsType1-DsType2)==1 | abs(DsType1-DsType2)==3 ; % if they are 90
            Ctemp = 1 ;
        elseif DsType1==DsType2 ; % if they are same prefered
            Ctemp = 2 ;
        else % if they are 180
            Ctemp = 3 ;
        end

        for cells1=1:length(i1) ;
            for cells2=1:length(i2) ;
                TypeGroupCorr{Ctemp} = [TypeGroupCorr{Ctemp},corrMat_zero(i1(cells1),i2(cells2))] ;
                TypeGroupDist{Ctemp} = [TypeGroupDist{Ctemp},DistMat(i1(cells1),i2(cells2))] ;
            end
        end
        
        TypeGroupHist{Ctemp} = hist(TypeGroupCorr{Ctemp}(TypeGroupDist{Ctemp}(:)<=EiCircRad...
            & TypeGroupDist{Ctemp}(:)>0),CorrHistX) ; % histogram of correlation for nearby cells 
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
             if distCtr{CellCnt}(Ni)<EiCircRad && ~ismember(psth_i(1),quad{CellCnt}(:,1)) ; % if the cell is very nearby and that type has not yet part of the quad
                quad{CellCnt}=[quad{CellCnt};psth_i] ; % add it to the quad     
             end
         end
         CellCnt = CellCnt+1 ;
    end
end
              
% calc vectors in quadruplets
psthNumCells = nan(length(quad),length(bins)) ;
psthMax = nan(length(quad),length(bins)) ;
psthMin = nan(length(quad),length(bins)) ;
psthSum = nan(length(quad),length(bins)) ;
psthMaxVect = nan(length(quad),length(bins)) ;
RealQuad_i = [] ; % index of quadruplets in quad

for q=1:length(quad) ; % for each putative quadruplet
    if size(quad{q},1)==4 ; % if it is a true quadruplet
        RealQuad_i = [RealQuad_i,q] ;
        for t=1:length(bins) ; % for each time point in the psth
            for cells=1:4 ;
                psthBlock(cells,:) = psth_norm{quad{q}(cells,1)}(quad{q}(cells,2),t) ;
            end
            
            psthNumCells(q,t) = sum(psthBlock>0) ; % how many cells are firing at the same time (0-4)
            
            psthMax(q,t) = max(psthBlock) ; % max relative firing rate
            psthMin(q,t) = min(psthBlock) ; % minimum relative firing rate
            psthSum(q,t) = sum(psthBlock) ; % sum of firing rates
            

            PV = PolarVectorAddition([TypeAngle(quad{q}(:,1))',psthBlock]) ; %        
            quad_vectAngle{q}(t) = PV(1) ;
            quad_vectMag{q}(t) = PV(2) ;
            quad_vectAngle_sum(t) = quad_vectAngle_sum(t) + PV(1) ;
            
% 
%             OppPairs_i = [find(quad{q}(:,1)==1),find(quad{q}(:,1)==3),...
%                 find(quad{q}(:,1)==2),find(quad{q}(:,1)==4)] ; % index of opposing cell pairs
%             psthMaxVect(q,t) = sqrt(max(psthBlock(OppPairs_i(1:2)))^2 + max(psthBlock(OppPairs_i(3:4)))^2) ; % max vector if opposite directions were not counteracting
        end
    end
end
  
if ~isempty(RealQuad_i) ; % if there are any quadruplets
    psthMaxOverSum_dist = hist(psthMax(:)./psthSum(:),psthMaxOverSum_distX) ;
    psthNumCells_dist = hist(psthNumCells(:),[0:4]) ;
end

% sum all quad vectors
for t=1:length(bins) ; % for each time point in the psth
    for rqi = RealQuad_i ; 
        TempBlock(rqi,:) = [quad_vectAngle{rqi}(t),quad_vectMag{rqi}(t)] ; %
    end
    quad_vect_sum(t,:) = PolarVectorAddition(TempBlock) ;
end

% calculate vector sum over larger groups (average within a cell type)
Group_vectAngle = nans(numCells,length(bins)) ;
Group_vectMag = nans(numCells,length(bins)) ;
for t=1:length(bins) ; % for each time point in the psth
    for cells1=1:numCells ; % for each cell
        psthBlock = nan(1,4) ; % prep a space for each type
        psthBlockNum = zeros(1,4) ; % number of cells of each type
        for cells2=1:numCells ; % for each cell
            if distCtr{cells1}(cells2)< VectorSumDistMax ; % if its within the distance
                psthBlock(TypeCellMat(cells2,1)) = nansum([psth_norm{TypeCellMat(cells2,1)}(TypeCellMat(cells2,2),t),...
                    psthBlock(TypeCellMat(cells2,1))]) ; 
                psthBlockNum(TypeCellMat(cells2,1)) = psthBlockNum(TypeCellMat(cells2,1)) + 1 ;
            end
        end
        if sum(isnan(psthBlock))==0 ; % if there is at least one of each type
            PV = PolarVectorAddition([TypeAngle',psthBlock'./psthBlockNum']) ; %        
            Group_vectAngle(cells1,t) = PV(1) ;
            Group_vectMag(cells1,t) = PV(2) ;
        end      
        Group_vectMaxDivSum(cells1,t) = max(psthBlock'./psthBlockNum')/sum((psthBlock'./psthBlockNum')) ;
    end
end
Group_MaxOverSum_dist = hist(Group_vectMaxDivSum(:),psthMaxOverSum_distX) ;

% movie and transforms
Temp = load(DataBlock(DB).MoviePath{Params.MovieNum}) ;% load movie
for t=1:size(Temp.mov,3) ;
    NatStimMovieMat(:,:,t) = Temp.mov(:,:,t)' ; % the movie is loaded as transpose of displayed movie
end
clear Temp ;
NumFrames = size(NatStimMovieMat,3) ; % number of movie frames

load(DataBlock(DB).TformEiPath) ;% load EI-->monitor transform

CompMovieTime = [0:TimeStep:stim_time(2)] ; % time points of movie/psth composite
MovieTime = [0:NumFrames-1]*(FrameInterval/FramesPerSecond) ; % time points of movie

Temp = tformfwd(Tform, reshape(eiAxisLims,2,2)) ;
MonAxisLims = [sort(Temp(:,1))',sort(Temp(:,2))'] ;

MovieX = [1+DataBlock(DB).xstart(Params.MovieNum) ,...
    (DataBlock(DB).xstart(Params.MovieNum)+size(NatStimMovieMat,2))*MovieStixelWidth] ;
MovieY = [1+DataBlock(DB).ystart(Params.MovieNum),...
    (DataBlock(DB).ystart(Params.MovieNum)+size(NatStimMovieMat,1))*MovieStixelWidth] ;

%StdGaussBlur = 4 ;
%NatStimMovieMat_Blurred = MovBlurring(NatStimMovieMat,StdGaussBlur) ;

% save entire workspace
save([saveFigPath,'entire'],'-v7.3')

end ; % END LoadOldMatFilesFlag
% ForIgor

% psth histogram
DbId = ['Db',num2str(DB),'Movie',num2str(Params.MovieNum)] ;
Identifier = [DbId,'psthHistX'] ;
ForIgor.(Identifier) = psth_histX ;

for DsType = 1:length(DsTypeName) ;
    Identifier = [DbId,'psthHistType',num2str(DsType)] ;
    ForIgor.(Identifier) = mean(psth_hist{DsType},1) ;
end

% pairwise correlations as function of distance
for DsType1=1:length(DsTypeName) ; % for each type
    for DsType2=1:length(DsTypeName) ; % for each type
        Identifier = [DbId,'CorrX','Type',num2str(DsType1),'v',num2str(DsType2)] ;
        ForIgor.(Identifier) = DistBinCenters ;

        Identifier = [DbId,'Corr','Type',num2str(DsType1),'v',num2str(DsType2)] ;
        ForIgor.(Identifier) = CorrByDistMean{DsType1}(DsType2,:) ;  
        
        Identifier = [DbId,'CorrStd','Type',num2str(DsType1),'v',num2str(DsType2)] ;
        ForIgor.(Identifier) = CorrByDistStd{DsType1}(DsType2,:) ; 
    end
end

% quadurplet histograms
Identifier = [DbId,'QuadMaxFracHistX'] ;
ForIgor.(Identifier) = psthMaxOverSum_distX ;

Identifier = [DbId,'QuadMaxFracHist'] ;
ForIgor.(Identifier) = psthMaxOverSum_dist/sum(psthMaxOverSum_dist) ;

Identifier = [DbId,'QuadNumSynchHistX'] ;
ForIgor.(Identifier) = [0:4] ;

Identifier = [DbId,'QuadNumSynchHist'] ;
ForIgor.(Identifier) = psthNumCells_dist/sum(psthNumCells_dist) ;

% larger vector group histogram
Identifier = [DbId,'VectGroupMaxFracHistX'] ;
ForIgor.(Identifier) = psthMaxOverSum_distX ;

Identifier = [DbId,'VectGroupMaxFracHist'] ;
ForIgor.(Identifier) = Group_MaxOverSum_dist/sum(Group_MaxOverSum_dist) ;

% for matlab
ForMatlab = struct ;

% pairwise correlation histograms by perfered direction relationship
ForMatlab.(['DB',num2str(DB)]).(['Movie',num2str(Params.MovieNum)]).CorrHistX = CorrHistX ;
ForMatlab.(['DB',num2str(DB)]).(['Movie',num2str(Params.MovieNum)]).TypeGroupHist = TypeGroupHist ;


% figures
if PlotFigsFlag ;

    % DS spike rates during movie
    figure
    plotQuadsFlag = 0;
    for t=1:length(CompMovieTime) ; % for each time point in composite image
        PsthPnt = interp1(bins,[1:length(bins)],CompMovieTime(t),'nearest') ;
        MoviePnt = interp1(MovieTime,[1:NumFrames], CompMovieTime(t),'nearest') ;

        colormap gray
        imagesc(MovieX,MovieY,NatStimMovieMat(:,:,MoviePnt),[0,255]) ; % plot movie frame
        hold on

        CellCnt = 1 ;
        for DsType=1:length(DsTypeName) ; 
            for cells = 1:length(cell_id{DsType}) ; 

                ctrMon = tformfwd(Tform,ctr(CellCnt,1),ctr(CellCnt,2)) ;

                plot(ctrMon(1),ctrMon(2),'o','Color',Color_list{DsType})
                hold on

                if ~isnan(psth_norm{DsType}(cells,PsthPnt)) ;

                    spikeNorm = EiSpikeRateFactor*psth_norm{DsType}(cells,PsthPnt) ; % normalized spike rate

                    deltaX = spikeNorm*cosd(TypeAngle(DsType)) ;
                    deltaY = spikeNorm*sind(TypeAngle(DsType)) ;

                    X=[ctr(CellCnt,1),ctr(CellCnt,1)+deltaX] ;
                    Y=[ctr(CellCnt,2),ctr(CellCnt,2)+deltaY] ;

                    [xMon,yMon] = tformfwd(Tform, X',Y') ;

                    plot(xMon,yMon,'Color',Color_list{DsType})

                    % plot quadruplet vectors
                    if plotQuadsFlag ;
                        for q=1:length(quad) ; % for each putative quadruplet
                            if size(quad{q},1)==4 ; % if it is a true quadruplet
                                if sum(ismember(quad{q},[DsType,cells],'rows'))>0 ; % if this cell is a member of a quadruplet
                                    spikeNorm = EiSpikeRateFactor*quad_vectMag{q}(PsthPnt) ; % normalized spike rate

                                    deltaX = spikeNorm*cosd(quad_vectAngle{q}(PsthPnt)) ;
                                    deltaY = spikeNorm*sind(quad_vectAngle{q}(PsthPnt)) ;

                                    X = [ctr(CellCnt,1),ctr(CellCnt,1)+deltaX] ;
                                    Y = [ctr(CellCnt,2),ctr(CellCnt,2)+deltaY] ;

                                    [xMon,yMon] = tformfwd(Tform, X',Y') ;

                                    plot(xMon,yMon,'y')
                                end
                            end
                        end 
                    end
                end

                CellCnt = CellCnt + 1 ; 
            end
        end
        
        
        axis(MonAxisLims)
        M(t)=getframe ;
        hold off
    end

    save([saveFigPath,'movieZoomQuads'],'M','-v7.3') ; % save movie

    % Quad sum vector during movie
    figure
    quadFactor = 10 ;
    quadCenterXY = [0,0] ;
    for t=1:length(CompMovieTime) ; % for each time point in composite image
        PsthPnt = interp1(bins,[1:length(bins)],CompMovieTime(t),'nearest') ;
        MoviePnt = interp1(MovieTime,[1:NumFrames], CompMovieTime(t),'nearest') ;

        colormap gray
        imagesc(MovieX,MovieY,NatStimMovieMat(:,:,MoviePnt),[0,255]) ; % plot movie frame
        hold on

        deltaX = quadFactor*quad_vect_sum(PsthPnt,2)*cosd(quad_vect_sum(PsthPnt,1)) ;
        deltaY = quadFactor*quad_vect_sum(PsthPnt,2)*sind(quad_vect_sum(PsthPnt,1)) ;

        X = [quadCenterXY(1),quadCenterXY(1)+deltaX] ;
        Y = [quadCenterXY(2),quadCenterXY(2)+deltaY] ;

        [xMon,yMon] = tformfwd(Tform, X',Y') ;

        plot(xMon,yMon,'r')
        
        axis(MonAxisLims)
        M(t)=getframe ;
        hold off
    end
        

    
    % Ei mosaics
    CellCnt = 1 ;
    for DsType=1:length(DsTypeName) ; 
        figure
        for cells = 1:length(cell_id{DsType}) ;
            plot_ei(dataRun,cell_id{DsType}(cells),'pos_color',Color_list{cells},'neg_color',Color_list{cells},'cutoff',.25,'scale',2, 'coordinates', 'array')
            %plot_ei(dataRun,cell_id{DsType}(cells),'pos_color',Color_list{cells},'neg_color',Color_list{cells})
            hold on
        end

        figure
        for cells = 1:length(cell_id{DsType}) ;
            plot(ctr(CellCnt,1),ctr(CellCnt,2),'+','Color',Color_list{DsType})
            hold on
            drawCircle(ctr(CellCnt,1),ctr(CellCnt,2),100,'Color',Color_list{DsType})
            CellCnt = CellCnt + 1 ;
        end
    end

    % Spike Rasters for each cell
    for DsType = 1:length(DsTypeName) ; % for each ds type
        for cells = 1:length(cell_i{DsType}) ; % for each cell

            spikes = dataRun.spikes{cell_i{DsType}(cells)} ;

            figure(1); clf
            get_raster(spikes(spikes>trigs(1)), trigs(1:repeatNum),'tic_color','k',...
                'foa',-1,'axis_range', [stim_time 0 repeatNum]) ;
            drawnow

            title([DsTypeName{DsType},' Cell id: ',num2str(cell_id{DsType}(cells))])
            %print(gcf, '-djpeg', [saveFigPath,'rasterCell',num2str(cells)])
            pause
        end
    end

    % all ei centers
    figure 
    CellCnt = 1 ;
    for DsType = 1:length(DsTypeName) ; % for each ds type
        for cells = 1:length(cell_i{DsType}) ; % for each cell       
            plot(ctr(CellCnt,1),ctr(CellCnt,2),'+','Color',Color_list{DsType})
            hold on
            drawCircle(ctr(CellCnt,1),ctr(CellCnt,2),EiCircRad,'Color',Color_list{DsType})

            CellCnt = CellCnt+1 ;
        end
    end

    % plot ei and psths
    psthPlots = figure ;
    rfPlots = figure ;

    CellCnt = 1 ; % start cell counter
    for DsType = 1:length(DsTypeName) ; % for each ds type
        for cells = 1:length(cell_i{DsType}) ; % for each cell
            figure(psthPlots); clf
            ax(1) = subplot(MaxNumCells,1,1)
            plot(psth{DsType}(cells,:),'Color',Color_list{DsType}) ;

            figure(rfPlots) ; clf
            plot(ctr(CellCnt,1),ctr(CellCnt,2),'+','Color',Color_list{DsType})
            hold on
            drawCircle(ctr(CellCnt,1),ctr(CellCnt,2),EiCircRad,'Color',Color_list{DsType},'lineWidth',2,'lineStyle','--')

            numCellsPloted = 2 ;
            for NearCell = 1:size(si,2) ; % for each ds cell
                if numCellsPloted<=MaxNumCells ; % if you havn't plotted the max number of cells yet
                    Ni = si(CellCnt,NearCell+1) ; % near cell index
                    psth_i = TypeCellMat(Ni,:) ; % psth index
                    if ~isnan(psth{psth_i(1)}(psth_i(2),1)) ; % if the psth exists
                        figure(psthPlots)
                        ax(numCellsPloted)=subplot(MaxNumCells,1,numCellsPloted) ;
                        plot(psth{psth_i(1)}(psth_i(2),:),'Color',Color_list{psth_i(1)})

                        figure(rfPlots)
                        plot(ctr(Ni,1),ctr(Ni,2),'+','Color',Color_list{psth_i(1)})
                        hold on
                        drawCircle(ctr(Ni,1),ctr(Ni,2),EiCircRad,'Color',Color_list{psth_i(1)})

                        numCellsPloted = numCellsPloted + 1 ;
                    end
                end    
            end
            linkaxes(ax(:),'x')
            CellCnt = CellCnt+1 ;
            pause
        end
    end

    % look at putative quadruplets
    for q=1:length(quad) ; % for each putative quadruplet
        if size(quad{q},1)==4 ; % if it is a true quadruplet
            figure
            for cells=1:4 ; % for each member of the quad
                ci = find(ismember(TypeCellMat,quad{q}(cells,:),'rows')==1) ; % cell index
                subplot(6,3,1)
                plot(ctr(ci,1),ctr(ci,2),'+','Color',Color_list{quad{q}(cells,1)})
                hold on
                drawCircle(ctr(ci,1),ctr(ci,2),EiCircRad,'Color',Color_list{quad{q}(cells,1)})

                ax(cells) = subplot(6,3,(cells*3+1):(cells*3+3))
                plot(bins,psth{quad{q}(cells,1)}(quad{q}(cells,2),:),'Color',Color_list{quad{q}(cells,1)})
            end
            ax(5) = subplot(6,3,16:18)
            plotyy(bins,quad_vectMag{q},bins,quad_vectAngle{q})
            plot(bins,quad_vectMag{q},'k')
            hold on 
            plot(bins, psthMin(q,:),'r:')
            plot(bins, psthMax(q,:),'c--')
            plot(bins, psthMaxVect(q,:),'r:')
            linkaxes(ax(:),'x')
        end
    end

    % firing of DS cells without stimulus
    figure
    for t=1:length(bins) ; % for each PSTH time point
        CellCnt = 1 ;
        for DsType=1:length(DsTypeName) ; 
            for cells = 1:length(cell_id{DsType}) ;

                plot(ctr(CellCnt,1),ctr(CellCnt,2),'o','Color',Color_list{DsType})
                hold on

                spikeNorm = EiSpikeRateFactor*psth_norm{DsType}(cells,t) ; % normalized spike rate

                deltaX = spikeNorm*cosd(TypeAngle(DsType)) ;
                deltaY = spikeNorm*sind(TypeAngle(DsType)) ;

                plot([ctr(CellCnt,1),ctr(CellCnt,1)+deltaX],[ctr(CellCnt,2),ctr(CellCnt,2)+deltaY],'Color',Color_list{DsType})
                %quiver(ctr(1),ctr(2),deltaX,deltaY,'Color',Color_list{DsType})
                %hold on

                % plot quadruplet vectors
                for q=1:length(quad) ; % for each putative quadruplet
                    if size(quad{q},1)==4 ; % if it is a true quadruplet
                        if sum(ismember(quad{q},[DsType,cells],'rows'))>0 ; % if this cell is a member of a quadruplet
                            spikeNorm = EiSpikeRateFactor*quad_vectMag{q}(t) ; % normalized spike rate

                            deltaX = spikeNorm*cosd(quad_vectAngle{q}(t)) ;
                            deltaY = spikeNorm*sind(quad_vectAngle{q}(t)) ;

                            plot([ctr(CellCnt,1),ctr(CellCnt,1)+deltaX],[ctr(CellCnt,2),ctr(CellCnt,2)+deltaY],'k')
                        end
                    end
                end 
                CellCnt = CellCnt + 1 ;
            end
        end
        %axis(eiAxisLims)
        M(t)=getframe ;
        hold off
    end

    % psth histograms
    figure
    for DsType = 1:length(DsTypeName) ;
        plot(psth_histX,mean(psth_hist{DsType},1)/sum(mean(psth_hist{DsType},1)),'Color',Color_list{DsType})
        hold on
    end

    % quadurplet histograms
    figure
    plot(psthMaxOverSum_distX,psthMaxOverSum_dist/sum(psthMaxOverSum_dist)) ;

    figure
    plot([0:4],psthNumCells_dist/sum(psthNumCells_dist)) ;

    % correlations
    figure
    plot(DistMat(:),corrMat_zero(:),'b*')
    hold on
    plot(DistMat(:),corrMat_max(:),'ro')

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
            %errorbar(DistBinCenters,CorrByDistMean{DsType1}(DsType2,:),CorrByDistStd{DsType1}(DsType2,:),'-','Color',Color_list{DsType2})
        end
        title(num2str(DsType1))
    end

    figure % correlation of all pairs by prefered direction relationship
    for DsType1=1:length(DsTypeName) ; % for each type
        i1 = find(TypeCellMat(:,1)==DsType1) ;
        for DsType2=1:length(DsTypeName) ; % for each type
            i2 = find(TypeCellMat(:,1)==DsType2) ;
            
            if abs(DsType1-DsType2)==1 | abs(DsType1-DsType2)==3 ; % if they are 90
                Ctemp = 'k' ; Ntemp = 1 ;
            elseif DsType1==DsType2 ; % if they are same prefered
                Ctemp = 'c' ; Ntemp = 2 ;
            else % if they are 180
                Ctemp = 'r' ; Ntemp = 3 ;
            end
            
            for cells1=1:length(i1) ;
                for cells2=1:length(i2) ;
                    subplot(2,1,1)
                    plot(DistMat(i1(cells1),i2(cells2)),corrMat_zero(i1(cells1),i2(cells2)),'*','Color',Ctemp)
                    hold on
                end
            end
            xlabel('dist (um)')
            ylabel('corr coef')
            
            subplot(2,1,2)
            plot(CorrHistX,TypeGroupHist{Ntemp},'Color',Ctemp)
            hold on
            xlabel('corr coef')
            ylabel('number of observations')
        end
    end

    % DS spike synchrony during movie
    figure
    for t=1:length(CompMovieTime) ; % for each time point in composite image
        PsthPnt = interp1(bins,[1:length(bins)],CompMovieTime(t),'nearest') ;
        MoviePnt = interp1(MovieTime,[1:NumFrames], CompMovieTime(t),'nearest') ;

        colormap gray
        imagesc(MovieX,MovieY,NatStimMovieMat(:,:,MoviePnt),[0,255]) ; % plot movie frame
        hold on

        CellCnt = 1 ;
        for DsType=1:length(DsTypeName) ; 
            for cells = 1:length(cell_id{DsType}) ; 

                ctrMon = tformfwd(Tform,ctr(CellCnt,1),ctr(CellCnt,2)) ;

                plot(ctrMon(1),ctrMon(2),'o','Color',Color_list{DsType})
                hold on

                for cells2 = 1:length(cell_id{DsType}) ; % for each other cell of that type
                    if cells~=cells2 ; % if its not the same cell
                        ci2 = find(ismember(TypeCellMat,[DsType,cells2],'rows')==1) ;
                        if distCtr{CellCnt}(ci2)<synchronyDistMax ; % if its within the max distance
                            geoMean = sqrt(psth_norm{DsType}(cells,PsthPnt)*...
                                psth_norm{TypeCellMat(ci2,1)}(TypeCellMat(ci2,2),PsthPnt)) ; % geometric mean
                            if geoMean>0 ;
                                ctrMon2 = tformfwd(Tform,ctr(ci2,1),ctr(ci2,2)) ;
                                plot([ctrMon(1),ctrMon2(1)],[ctrMon(2),ctrMon2(2)],'Color',Color_list{DsType},...
                                    'LineWidth',geoMean/SpikeToLineWidthFact)
                            end
                        end
                    end
                end
                CellCnt = CellCnt + 1 ;
            end
        end
        axis(MonAxisLims)
        M(t)=getframe ;
        hold off
    end

    save([saveFigPath,'SynchronyMovie'],'M','-v7.3') ; % save movie

    % vector sum in larger groups
    figure
    for t=1:length(CompMovieTime) ; % for each PSTH time point
        PsthPnt = interp1(bins,[1:length(bins)],CompMovieTime(t),'nearest') ;
        MoviePnt = interp1(MovieTime,[1:NumFrames], CompMovieTime(t),'nearest') ;

        colormap gray
        imagesc(MovieX,MovieY,NatStimMovieMat(:,:,MoviePnt),[0,255]) ; % plot movie frame
        hold on

        CellCnt = 1 ;
        for DsType=1:length(DsTypeName) ; 
            for cells = 1:length(cell_id{DsType}) ;

                ctrMon = tformfwd(Tform,ctr(CellCnt,1),ctr(CellCnt,2)) ;

                plot(ctrMon(1),ctrMon(2),'o','Color',Color_list{DsType})

                % group vectors
                if ~isnan(Group_vectMag(CellCnt,PsthPnt)) ; % if this cell has a group vector
                    spikeNorm = EiSpikeRateFactor*Group_vectMag(CellCnt,PsthPnt) ; % normalized spike rate

                    deltaX = spikeNorm*cosd(Group_vectAngle(CellCnt,PsthPnt)) ;
                    deltaY = spikeNorm*sind(Group_vectAngle(CellCnt,PsthPnt)) ;

                    X=[ctr(CellCnt,1),ctr(CellCnt,1)+deltaX] ;
                    Y=[ctr(CellCnt,2),ctr(CellCnt,2)+deltaY] ;

                    [xMon,yMon] = tformfwd(Tform, X',Y') ;

                    plot(xMon,yMon,'r')
                end

                CellCnt = CellCnt + 1 ;
            end
        end
        axis(MonAxisLims)
        M(t)=getframe ;
        hold off
    end    

    save([saveFigPath,'VectorMovieZoom'],'M','-v7.3') ; % save movie

    figure % entire movie
    colormap gray
    for t=1:NumFrames ; % for each movie point
        imagesc(MovieX,MovieY,NatStimMovieMat(:,:,t),[0,255]) ;
        pause
    end
    
    % Spike Rasters for each trial - all cells
    for t=0:bins(end) ; % for each second
    figure(1); clf
    Temp = 0 ;

        for DsType = 1:length(DsTypeName) ; % for each ds type
            for cells = 1:length(cell_i{DsType}) ; % for each cell
                Temp = Temp+1 ;
                spikes = dataRun.spikes{cell_i{DsType}(cells)} ;
                spkT = spikes(spikes>Params.TimeBounds(2)+t & spikes<(Params.TimeBounds(2)+t+1)) ;
                if ~isempty(spkT) ;
                    %plot(spkT, Temp,'.','Color',Color_list{DsType}) ;
                    for s=1:length(spkT) ;
                        plot([spkT(s),spkT(s)], [-.5,.5]+Temp,'-','Color',Color_list{DsType}) ;
                        hold on
                    end
                    
                    drawnow
                end
                hold on
            end
        end
        hold off
        pause
    end
 
    
    % quadruplets ei cntrs and rasters
    for q=1:length(quad) ; % for each putative quadruplet
        if size(quad{q},1)==4 ; % if it is a true quadruplet
            for cells=1:4 ; % for each member of the quad
                ci = find(ismember(TypeCellMat,quad{q}(cells,:),'rows')==1) ; % cell index
                figure(1); clf
                plot(ctr(ci,1),ctr(ci,2),'+','Color',Color_list{quad{q}(cells,1)})
                hold on
                drawCircle(ctr(ci,1),ctr(ci,2),EiCircRad,'Color',Color_list{quad{q}(cells,1)})
            end
        end
    end

   
   for q=1:length(quad) ; % for each putative quadruplet
        if size(quad{q},1)==4 ; % if it is a true quadruplet
            for t=0:bins(end) ; % for each second
                figure(2); clf
                Temp = 0 ;
                
                for cells=1:4 ; % for each member of the quad
                    Temp = Temp+1 ;
                    spikes = dataRun.spikes{cell_i{quad{q}(cells,1)}(quad{q}(cells,2))} ;
                    spkT = spikes(spikes>Params.TimeBounds(2)+t & spikes<(Params.TimeBounds(2)+t+1)) ;
                    if ~isempty(spkT) ;
                        for s=1:length(spkT) ;
                            plot([spkT(s),spkT(s)], [-.5,.5]+Temp,'-','Color',Color_list{cells}) ;
                            hold on
                        end
                        drawnow
                    end
                end
                hold off
                pause
            end
        end
    end
            
            
    % poplution figures

%     % prep population stats
%     if ~isfield(ForIgor,'psth_hist') ;
%         ForIgor.N = 0 ;
%         ForIgor.psth_hist = 0 ;
%         ForIgor.psthMaxOverSum_dist = 0 ;
%     end
%     ForIgor.N = ForIgor.N + 1 ;  
% 
%     ForIgor.psth_hist = ForIgor.psth_hist + mean(psth_hist{DsType},1) ;
%     ForIgor.psthMaxOverSum_dist = ForIgor.psthMaxOverSum_dist + psthMaxOverSum_dist/sum(psthMaxOverSum_dist) ;
% 
%     % psth histograms
%     figure
%     for DsType = 1:length(DsTypeName) ;
%         plot(psth_histX,mean(psth_hist{DsType},1),'Color',Color_list{DsType})
%         hold on
%     end
% 
%     % quadurplet histograms
%     figure
%     plot(psthMaxOverSum_distX,psthMaxOverSum_dist/sum(psthMaxOverSum_dist)) ;
%     hold on
% 
%     figure
%     plot([0:4],psthNumCells_dist/sum(psthNumCells_dist)) ;
%     hold on

figure % make a movie with all cells
cellIdsArray{1} = dataRun.cell_ids ; % make all 1 type
StimPsthMovie = StimPsthMovieMaker(dataRun, cellIdsArray, Tform, ...
    DataBlock(DB).MoviePath{Params.MovieNum}, trigs,...
    'MovieStixelWidth',MovieStixelWidth,'MovieXstart',DataBlock(DB).xstart(Params.MovieNum),...
    'MovieYstart',DataBlock(DB).ystart(Params.MovieNum),'MonitorBounds',MonAxisLims,...
    'MovieFrameInterval',FrameInterval,'SpikeDisplayScale', .5)
end

