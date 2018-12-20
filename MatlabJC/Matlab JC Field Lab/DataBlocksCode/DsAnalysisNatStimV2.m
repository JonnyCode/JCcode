function ForIgor = DsAnalysisNatStimV2(DataBlock, DB, Params) 

% this function will analyze DS cells during natural stimulus,
% adapted from 'DsAnalysisNatStim', to analyze data from concatinated data
% [drifting gratings, movie]

% JC 10/31/2016 

% parameters
Color_list = {'k','r','b','g','y','c','k','r','b','g','y','c','k','r','b','g','y','c',...
    'k','r','b','g','y','c','k','r','b','g','y','c','k','r','b','g','y','c'} ; % order of colors for each 
TrialTrigInterval = 10 ; % sec
EiCircRad = 100 ; % (um) radius of possible dendrites around ei center of mass
FramesPerTrig = 100 ;
FramesPerSecond = 60.35 ;
FrameInterval = 2 ; % frame is switched every...
saveFigPath = '/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/Db2' ;
eiAxisLims = [-200, 200, -200, 200] ; % Axis limits figures ([-425, 425, -425, 425] covers entire array) ;

% load data
dataRun = load_data(DataBlock(DB).MovieDsConcat) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ;

% identify DS cells
dataRunTemp = load_data(DataBlock(DB).DsPath) ;
dataRunTemp = load_neurons(dataRunTemp) ;

Params.TimeBounds = [0,dataRunTemp.duration] ; % times of Dg data
Params.SkipBwMap = true ;
DataBlock(DB).DsConcatPath = DataBlock(DB).MovieDsConcat ;
clear dataRunTemp

ForIgor = DsCellFinder(DataBlock, DB, Params) ;

cell_id = ForIgor.ds_id{2} ; % On-Off cell 
DsTypeName = ForIgor.dsName{2} ; 
    
for DsType=1:length(DsTypeName) ;
    cell_i{DsType} = get_cell_indices(dataRun, cell_id{DsType}) ;
end

% calculate movie triggers
trigNum = ceil(DataBlock(DB).MovieRepFrameNum(1)/FramesPerTrig) ; % number of triggers per repeat
FirstTrig = find(dataRun.triggers>=Params.TimeBounds(2),1,'first') ; % first trigger of the movie
trigs = dataRun.triggers(FirstTrig:trigNum:end) ; % begining of each movie 
trigEnds = dataRun.triggers(FirstTrig+trigNum-1:trigNum:end) ; % last trigger of each movie 
repeatNum= floor(length(dataRun.triggers(FirstTrig:end))/trigNum) ; % show how many repeats this is
stim_time =[0, DataBlock(DB).MovieRepFrameNum(1)/FramesPerSecond] ;

% Question # how good are the ds mosaics?
for DsType=1:length(DsTypeName) ; 
    figure
    for cells = 1:length(cell_id{DsType}) ;
        plot_ei(dataRun,cell_id{DsType}(cells),'pos_color',Color_list{cells},'neg_color',Color_list{cells},'cutoff',.25,'scale',2, 'coordinates', 'array')
        hold on
    end
    
    figure
    for cells = 1:length(cell_id{DsType}) ;
        ctr = get_ei_com(dataRun, cell_id{DsType}(cells), 2) ;
        
        plot(ctr(1),ctr(2),'+','Color',Color_list{DsType})
        hold on
        drawCircle(ctr(1),ctr(2),100,'Color',Color_list{DsType})
    end
end

% Question #1 do the DS cells respond during the natural movie?
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell

        spikes = dataRun.spikes{cell_i{DsType}(cells)} ;

        figure(1); clf
        get_raster(spikes(spikes>trigs(1)), trigs(1:repeatNum),'tic_color','k',...
            'foa',-1,'axis_range', [stim_time 0 repeatNum]) ;
        drawnow

        title([DsTypeName{DsType},' Cell id: ',num2str(cell_id{DsType}(cells))])
        %print(gcf, '-djpeg', [saveFigPath,'Cell',num2str(cells)])
        pause

    end
end
                
% Question #2 - what are the relative firing rates of similar / disimilar cell types
MaxNumCells = 10 ; % number of psths of the nearest cells that you want to plot
bin_size =  0.1 ; % (s) psth bin size

% get psths for all ds cells
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
        if ~isnan(cell_i{DsType}(cells)) ;
            spikes = dataRun.spikes{cell_i{DsType}(cells)} ;
            [psth{DsType}(cells,:), bins] = get_psth(spikes(spikes>Params.TimeBounds(2)), trigs(trigs>Params.TimeBounds(2)),...
                'stop',stim_time(2),'bin_size',bin_size) ;
            %[psth,psthTime] = get_smooth_psth(spikes(spikes>Params.TimeBounds(2)), trigs(trigs>Params.TimeBounds(2)), 'stop', stim_time(2)) ;
        else
            psth{DsType}(cells,:) = nans(1,ceil(stim_time(2)/bin_size)) ;
        end
    end
end

% plot eis
numElectrodeLayers = 2 ;

figure 
CellCnt = 1 ;
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
        ctr(CellCnt,:) = get_ei_com(dataRun, cell_id{DsType}(cells), numElectrodeLayers) ;
        
        TypeCellMat(CellCnt,:) = [DsType,cells] ;
        
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
        deltaCtr = ctr-repmat(ctr(CellCnt,:),size(ctr,1),1) ; % X,Y distance between that cell and all others
        distCtr{CellCnt} = sqrt(deltaCtr(:,1).^2+deltaCtr(:,2).^2) ; % euclidean distance
        [s,si(CellCnt,:)] = sort(distCtr{CellCnt}) ; % in order of nearest to farthest

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
        %pause
    end
end

% Question #3 - how often do quaduplets have minimal vector sum depsite high individiual firing rates?

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
              
% look at putative quadruplets
for q=1:length(quad) ; % for each putative quadruplet
    if size(quad{q},1)==4 ; % if it is a true quadruplet
        for t=1:length(bins) ; % for each time point in the psth
            for cells=1:4 ;
                psthBlock(cells,:) = psth{quad{q}(cells,1)}(quad{q}(cells,2),t)/max(psth{quad{q}(cells,1)}(quad{q}(cells,2),:)) ;
            end
            PV = PolarVectorAddition([quad{q}(:,1)*90-90,psthBlock]) ; % assumes [up, right, down, left]
            OppPairs_i = [find(quad{q}(:,1)==1),find(quad{q}(:,1)==3),...
                find(quad{q}(:,1)==2),find(quad{q}(:,1)==4)] ; % index of opposing cell pairs
            psthMax(q,t) = sqrt(max(psthBlock(OppPairs_i(1:2)))^2 + max(psthBlock(OppPairs_i(3:4)))^2) ; % max vector if opposite directions were not counteracting
            psthMin(q,t) = min(psthBlock) ; % minimum relative firing rate
            
            quad_vectAngle{q}(t) = PV(1) ;
            quad_vectMag{q}(t) = PV(2) ;
        end
            
%         figure
%         for cells=1:4 ; % for each member of the quad
%             ci = find(ismember(TypeCellMat,quad{q}(cells,:),'rows')==1) ; % cell index
%             subplot(6,3,1)
%             plot(ctr(ci,1),ctr(ci,2),'+','Color',Color_list{quad{q}(cells,1)})
%             hold on
%             drawCircle(ctr(ci,1),ctr(ci,2),EiCircRad,'Color',Color_list{quad{q}(cells,1)})
%             
%             ax(cells) = subplot(6,3,(cells*3+1):(cells*3+3))
%             plot(bins,psth{quad{q}(cells,1)}(quad{q}(cells,2),:)/max(psth{quad{q}(cells,1)}(quad{q}(cells,2),:)),'Color',Color_list{quad{q}(cells,1)})
%         end
%         ax(5) = subplot(6,3,16:18)
        %plotyy(bins,quad_vectMag{q},bins,quad_vectAngle{q})
%         plot(bins,quad_vectMag{q},'k')
%         hold on 
%         plot(bins, psthMin(q,:),'r:')
%         plot(bins, psthMax(q,:),'c--')
%         linkaxes(ax(:),'x')
    end
end
            
            
% Question #4 - how strongly correlated are quadruplet vectors
EiSpikeRateFactor = 100 ; % (um/spike)
TypeAngle = [120,30,300,210] ; % prefered angle for each DsType

figure
for t=1:length(bins) ; % for each PSTH time point
    for DsType=1:length(DsTypeName) ; 
        for cells = 1:length(cell_id{DsType}) ;
            ctr = get_ei_com(dataRun, cell_id{DsType}(cells), 2) ;

            plot(ctr(1),ctr(2),'o','Color',Color_list{DsType})
            hold on
            
            spikeNorm = EiSpikeRateFactor*(psth{DsType}(cells,t)/max(psth{DsType}(cells,:))) ; % normalized spike rate
            
            deltaX = spikeNorm*cosd(TypeAngle(DsType)) ;
            deltaY = spikeNorm*sind(TypeAngle(DsType)) ;

            plot([ctr(1),ctr(1)+deltaX],[ctr(2),ctr(2)+deltaY],'Color',Color_list{DsType})
            %quiver(ctr(1),ctr(2),deltaX,deltaY,'Color',Color_list{DsType})
            %hold on
            
            % plot quadruplet vectors
            for q=1:length(quad) ; % for each putative quadruplet
                if size(quad{q},1)==4 ; % if it is a true quadruplet
                    if sum(ismember(quad{q},[DsType,cells],'rows'))>0 ; % if this cell is a member of a quadruplet
                        spikeNorm = EiSpikeRateFactor*quad_vectMag{q}(t) ; % normalized spike rate

                        deltaX = spikeNorm*cosd(quad_vectAngle{q}(t)) ;
                        deltaY = spikeNorm*sind(quad_vectAngle{q}(t)) ;

                        plot([ctr(1),ctr(1)+deltaX],[ctr(2),ctr(2)+deltaY],'c')
                    end
                end
            end        
        end
    end
    axis(eiAxisLims)
    M(t)=getframe ;
    hold off
end


% movie with ds cell responses

Temp = load(DataBlock(DB).MoviePath{1}) ;% load movie
for t=1:size(Temp.mov,3) ;
    NatStimMovieMat(:,:,t) = Temp.mov(:,:,t)' ; % the movie is loaded as transpose of displayed movie
end
clear Temp ;
NumFrames = size(NatStimMovieMat,3) ; % number of movie frames

load(DataBlock(2).TformEiPath) ;% load EI-->monitor transform
%tforminv(Tform,[0;1],[0;0]) ;

TimeStep = .016 ; % (s) time steps of composite movie/psth composite
CompMovieTime = [0:TimeStep:mean(trigEnds-trigs(1:repeatNum))] ; % time points of movie/psth composite
MovieTime = [0:NumFrames-1]*(FrameInterval/FramesPerSecond) ; % time points of movie

Temp = tformfwd(Tform, reshape(eiAxisLims,2,2)) ;
MonAxisLims = [sort(Temp(:,1))',sort(Temp(:,2))'] ;

figure
for t=1:length(CompMovieTime) ; % for each time point in composite image
    PsthPnt = interp1(bins,[1:length(bins)],CompMovieTime(t),'nearest') ;
    MoviePnt = interp1(MovieTime,[1:NumFrames], CompMovieTime(t),'nearest') ;
    
    colormap gray
    imagesc(NatStimMovieMat(:,:,MoviePnt),[0,255]) ; % plot movie frame
    hold on
    
    for DsType=1:length(DsTypeName) ; 
        for cells = 1:length(cell_id{DsType}) ; 
            ctr = get_ei_com(dataRun, cell_id{DsType}(cells), 2) ;
            
            ctrMon = tformfwd(Tform,ctr(1),ctr(2)) ;
            
            plot(ctrMon(1),ctrMon(2),'o','Color',Color_list{DsType})
            hold on
            
            spikeNorm = EiSpikeRateFactor*(psth{DsType}(cells,PsthPnt)/max(psth{DsType}(cells,:))) ; % normalized spike rate
            
            deltaX = spikeNorm*cosd(TypeAngle(DsType)) ;
            deltaY = spikeNorm*sind(TypeAngle(DsType)) ;

            X=[ctr(1),ctr(1)+deltaX] ;
            Y=[ctr(2),ctr(2)+deltaY] ;
            
            [xMon,yMon] = tformfwd(Tform, X',Y') ;
            
            plot(xMon,yMon,'Color',Color_list{DsType})

        end
    end
    axis(MonAxisLims)
    M(t)=getframe ;
    hold off
end


end
    
    



