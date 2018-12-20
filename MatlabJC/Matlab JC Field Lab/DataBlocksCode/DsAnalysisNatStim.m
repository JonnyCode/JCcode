function ForIgor = DsAnalysisNatStim(DataBlock, DB, Params) 

% this function will analyze DS cells during natural stimulus

% JC 10/31/2016 

% parameters
Color_list = {'k','r','b','g','y','v'} ; % order of colors for each 
TrialTrigInterval = 10 ; % sec
saveFigPath = '/Users/jcafaro/Documents/AnalysisFigures/NatStimDs/Db2' ;

% binary white data run
% load Master data
dataRunMaster = load_data(DataBlock(DB).BwPath{1}) ;
dataRunMaster = load_neurons(dataRunMaster) ;
dataRunMaster = load_ei(dataRunMaster, 'all') ;
dataRunMaster = load_params(dataRunMaster,'cell_type_depth', 5) ;
dataRunMaster = load_sta(dataRunMaster) ; % only necessary to get trf_time

% load receptive fields and filter
marks_params.thresh = 4.5;
dataRunMaster = get_sta_summaries(dataRunMaster, 'all','marks_params', marks_params);
 
filt_params.radius = 0.75;
dataRunMaster = get_rfs_filtered(dataRunMaster, 'all', 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');

% load Natural Stim response
dataRun = load_data(DataBlock(DB).MovieRepPath{1}) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ;

% map and find ds cells
cell_list_map = map_ei(dataRunMaster, dataRun) ;

DsTypeCnt = 1 ;
for cellType = 1:length(dataRunMaster.cell_types) ; % for each cell type
    if strncmp(dataRunMaster.cell_types{cellType}.name,'Ds',2) ; % if its a 'Ds' cell
        Master_i{DsTypeCnt} = get_cell_indices(dataRunMaster,dataRunMaster.cell_types{cellType}.name) ;
        DsTypeName{DsTypeCnt} = dataRunMaster.cell_types{cellType}.name ; % name of ds type
        for cells=1:length(Master_i{DsTypeCnt}) ; % for each ds master cell
            if ~isempty(cell_list_map{Master_i{DsTypeCnt}(cells)}) ; % if there is a mapped cell
                cell_id = cell_list_map{Master_i{DsTypeCnt}(cells)} ;
                cell_i{DsTypeCnt}(cells) = get_cell_indices(dataRun,cell_id) ; % cell indicy in dataRun
            else
                cell_i{DsTypeCnt}(cells) = nan ;
            end
        end
        DsTypeCnt = DsTypeCnt + 1 ;
    end
end

% Question #1 do the DS cells respond during the natural movie?
trigNum = ceil(DataBlock(DB).MovieRepFrameNum(1)/100) ; % number of triggers per repeat
trigs = dataRun.triggers(1:trigNum:end) ; % 
repeatNum= floor(length(dataRun.triggers)/trigNum) ; % show how many repeats this is
stim_time =[0, DataBlock(DB).MovieRepFrameNum(1)/60] ;

for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
        if ~isnan(cell_i{DsType}(cells)) ;
            
            figure(1); clf
            get_raster(dataRun.spikes{cell_i{DsType}(cells)}, trigs,'tic_color','k',...
                'foa',-1,'axis_range', [stim_time 0 repeatNum]) ;
            drawnow

            title([DsTypeName{DsType},' Cell id: ',num2str(dataRunMaster.cell_ids(Master_i{DsType}(cells)))])
            print(gcf, '-djpeg', [saveFigPath,'Cell',num2str(cells)])
        end
    end
end
                
% Question #2 - what are the relative firing rates of similar / disimilar cell types
MaxNumCells = 10 ; % number of psths of the nearest cells that you want to plot
bin_size =  0.1 ;
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
        if ~isnan(cell_i{DsType}(cells)) ;
            [psth{DsType}(cells,:), bins] = get_psth(dataRun.spikes{cell_i{DsType}(cells)}, trigs, 'stop',stim_time(2),'bin_size',bin_size) ;
            %[psth,psthTime] = get_smooth_psth(dataRun.spikes{cell_i{DsType}(cells)}, trigs, 'stop', stim_time(2)) ;
        else
            psth{DsType}(cells,:) = nans(1,ceil(stim_time(2)/bin_size)) ;
        end
    end
end

figure 
coord_tform = coordinate_transform(dataRunMaster,'sta');
CellCnt = 1 ;
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
        Mi = Master_i{DsType}(cells) ;
        ctr(CellCnt,:) = dataRunMaster.stas.fits{Mi}.mean ;
        rad(CellCnt,:) = dataRunMaster.stas.fits{Mi}.sd ;
        angle(CellCnt,:) = dataRunMaster.stas.fits{Mi}.angle ;
        
        TypeCellMat(CellCnt,:) = [DsType,cells] ;
        
        [X,Y] = drawEllipse([ctr(CellCnt,:) rad(CellCnt,:) angle(CellCnt,:)]) ;
        if ~any(isnan([X,Y])) ;
            [X,Y] = tformfwd(coord_tform, X, Y) ;
            plot(X,Y,'Color',Color_list{DsType})
            hold on
        end
        
        CellCnt = CellCnt+1 ;
    end
end
  
CellCnt = 1 ;
for DsType = 1:length(DsTypeName) ; % for each ds type
    for cells = 1:length(cell_i{DsType}) ; % for each cell
        deltaCtr = ctr-repmat(ctr(CellCnt,:),size(ctr,1),1) ;
        distCtr = sqrt(deltaCtr(:,1).^2+deltaCtr(:,2).^2) ;
        [s,si(CellCnt,:)] = sort(distCtr) ;
        
        if ~isnan(cell_i{DsType}(cells)) ; % if this cell is mapped
            psthPlots = figure ;
            rfPlots = figure ;
            
            figure(psthPlots)
            ax(1) = subplot(MaxNumCells,1,1)
            plot(psth{DsType}(cells,:),'Color',Color_list{DsType}) ;
            
            figure(rfPlots)
            [X,Y] = drawEllipse([ctr(CellCnt,:) rad(CellCnt,:) angle(CellCnt,:)]) ;
            if ~any(isnan([X,Y])) ;
                [X,Y] = tformfwd(coord_tform, X, Y) ;
                plot(X,Y,'Color',Color_list{DsType},'lineWidth',2,'lineStyle','--')
                hold on
            end
            
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
                        [X,Y] = drawEllipse([ctr(Ni,:) rad(Ni,:) angle(Ni,:)]) ;
                        if ~any(isnan([X,Y])) ;
                            [X,Y] = tformfwd(coord_tform, X, Y) ;
                            plot(X,Y,'Color',Color_list{psth_i(1)})
                            hold on
                        end    
                        numCellsPloted = numCellsPloted + 1 ;
                    end
                end
            end    
        end
        linkaxes(ax(:),'x')
        CellCnt = CellCnt+1 ;
    end
end
        
        
        


% Question #3 - how often do quaduplets have minimal vector sum depsite high individiual firing rates?

% Question #4 - how strongly correlated are quadruplet vectors

l_cells = length(dataRun.spikes) ; % number of cells


    
end
    
    



