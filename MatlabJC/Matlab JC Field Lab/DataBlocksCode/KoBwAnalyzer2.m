function ForIgor = KoBwAnalyzer2(DataBlock, DB, Params)

% this function will analyze binary white stimuli response from array data
% before and after cells are knocked out (KO)

% This function assumes that cells are only typed in first data set.  Cells
% that are not mapped onto this are thrown out.

% JC 11/11/14
%DB=3; % TEMP

% parameters
subplot_X = 10 ;
mapEiFlag = DataBlock(DB).mapEiFlag ; % if datablocks should be mapped using map_ei
Color_list = {'k','r','g','y','b'} ; % order of colors for each 
numPCs = 3 ; % number of principle components
numDs = length(DataBlock(DB).DataPath) ; % number of data sets

% compare spatial and temporal receptive fields across data defined sets 
for ds = 1:numDs % for each data set

    % load data
    dataRun = load_data(DataBlock(DB).DataPath{ds}) ;
    dataRun = load_neurons(dataRun) ;

    % how different are the preko and postko response fits
    % spatial rf?
    marks_params.thresh = 3 ;
    dataRun = load_sta(dataRun) ;
    dataRun = load_params(dataRun) ;
    dataRun = get_sta_summaries(dataRun, 'all','marks_params', marks_params) ;
    dataRun = get_sta_fits_from_vision(dataRun) ;

    if ds == 1 ; % if Master Run
        % get list of cell types from Master dataRun
        for a = 1:length(dataRun.cell_types) ;
            celltypes{a} = dataRun.cell_types{a}.name ;
        end

        UniqueCellTypes = unique(celltypes) ;
        if isempty(UniqueCellTypes{1}) ;
            UniqueCellTypes = UniqueCellTypes(2:end) ;
        end
        
        for uc = 1:length(UniqueCellTypes) ;
            cell_ids{uc}{ds} = get_cell_ids(dataRun, UniqueCellTypes{uc}) ;
        end
        
    else % if slave run
        
        % load Master data
        dataRunMaster = load_data(DataBlock(DB).DataPath{1}) ;
        dataRunMaster = load_neurons(dataRunMaster) ;
        dataRunMaster = load_ei(dataRunMaster, 'all') ;
        dataRunMaster = load_params(dataRunMaster) ;
        
        if mapEiFlag ; % if using map ei cells
            % load electrical images

            dataRun = load_ei(dataRun, 'all') ;

            % map using electrical images
            cell_list_map{ds} = map_ei(dataRunMaster, dataRun) ;

            % cells ids in slave for each UniqueCellType set in master data
            for uc = 1:length(UniqueCellTypes) ;
                Tempi = get_cell_indices(dataRunMaster, UniqueCellTypes{uc}) ;
                cell_ids{uc}{ds} = cell2mat(cell_list_map{ds}(Tempi)) ;
            end
        else % if not using map ei
            for uc = 1:length(UniqueCellTypes) ;
                cell_ids{uc}{ds} = intersect(dataRun.cell_ids, get_cell_ids(dataRunMaster,UniqueCellTypes{uc})) ;
            end
        end
    end

    % plot spatial rfs
    for uc = 1:length(UniqueCellTypes) ;
        subplot_Y(uc) = ceil((length(cell_ids{uc}{1})+numPCs+3)/subplot_X) ;
        
        figure(uc) 
        subplot(subplot_Y(uc),subplot_X,1)
        title(UniqueCellTypes{uc})
        hold on
      
        if ~isempty(cell_ids{uc}) ; % if there are cells 
            %plot_rf_summaries(dataRun, cell_ids{uc}{ds}, 'plot_fits', true, 'foa',uc,'fit_color',Color_list{ds}) ; % plot mosiac
            
            coord_tform = coordinate_transform(dataRun,'sta');
            
            for c = 1:length(cell_ids{uc}{ds}) ; % for each cell
                ci = get_cell_indices(dataRun,cell_ids{uc}{ds}(c)) ; % cell index
                ctr = dataRun.stas.fits{ci}.mean ;
                rad = dataRun.stas.fits{ci}.sd ;
                angle = dataRun.stas.fits{ci}.angle ;
                [X,Y] = drawEllipse([ctr rad angle]) ;
                if ~any(isnan([X,Y])) ;
                    [X,Y] = tformfwd(coord_tform, X, Y) ;
                    plot(X,Y,'Color',Color_list{ds})
                    hold on
                end
            end   
            
        end

    end

    % plot and concatinate temporal rf
    for uc = 1:length(UniqueCellTypes) ;        
        if ~isempty(cell_ids{uc}{ds}) ;
            cell_i{uc}{ds} = get_cell_indices(dataRun, cell_ids{uc}{ds}) ;
            for b=1:length(cell_i{uc}{ds}) ;
                Trf{uc}{ds}(b,:) = flipud(dataRun.stas.time_courses{cell_i{uc}{ds}(b)}) ;
                Trf_norm{uc}{ds}(b,:) = Trf{uc}{ds}(b,:)/norm(Trf{uc}{ds}(b,:)) ;
            end
            figure(uc) 
            subplot(subplot_Y(uc),subplot_X,2)
            plot(mean(Trf_norm{uc}{ds}),Color_list{ds})
            hold on
            plot(mean(Trf_norm{uc}{ds})+std(Trf_norm{uc}{ds}),[Color_list{ds},':'])
            plot(mean(Trf_norm{uc}{ds})-std(Trf_norm{uc}{ds}),[Color_list{ds},':'])
            axis tight
            %title(UniqueCellTypes{uc})
        end
    end
end
    
% run PCA on termporal receptive fields
for uc = 1:length(UniqueCellTypes) ; % for each cell type in the master
    % find the empty groups
    ds_range = [] ;
    for ds=1:numDs ; % for each data set
        if ~isempty(cell_ids{uc}{ds}) ;
            ds_range = [ds_range,ds] ;
        end
    end
       
    % make stas all same length
    clear sta_pnts ;
    for ds=ds_range ; % for each data set
        sta_pnts(ds) = size(Trf_norm{uc}{ds},2) ; % sta lenth of 
    end
    sta_pntsMin = min(sta_pnts) ; % shortest sta

    for ds=ds_range ; % for each data set
        Trf_short{uc}{ds} = Trf_norm{uc}{ds}(:,1:sta_pntsMin) ; % cut it to the same length
    end
           
    % subtract off mean and keep data in grouped by celltype only
    Trf_mat{uc} = cell2mat(Trf_short{uc}') ; % make a matrix from all conditions 
    for b=1:size(Trf_mat{uc},1) ; % for each cell in all data sets
        Trf_mat_dev{uc}(b,:) = Trf_mat{uc}(b,:) - mean(Trf_mat{uc},1) ; % subtract off mean
    end

    % subtract off mean and keep data grouped by celltype and data set
    for ds=ds_range ; % for each data set defined
        for c=1:size(Trf_short{uc}{ds},1) ; % for each cell 
            Trf_dev{uc}{ds}(c,:) = Trf_short{uc}{ds}(c,:) - mean(Trf_mat{uc},1) ; % subtract off mean
        end
    end

    % find eigs of deviation matrix      
    Trf_cov = cov(Trf_mat_dev{uc}) ; % take covariance matrix
    [eigVec,eigVal] = eig(Trf_cov) ; % eigenvectors and values

    PCs{uc} = eigVec(:,end-numPCs+1:end) ;

    % project top PCs onto deviations
    temp_mat = nan(length(cell_ids{uc}{1}),numDs) ; % make matrix
    for c=1:numPCs ; % for each PC
        projValues{uc}{c} = temp_mat ; 
    end
    
    for ds=ds_range ; % for each data set defined
        if ds==1 ;
            i = [1:length(cell_ids{uc}{ds})] ;
        else 
            clear i
            if mapEiFlag ;
                Tempi = get_cell_indices(dataRunMaster, UniqueCellTypes{uc}) ;
                r=1 ;
                for c=1:length(Tempi) ; % for each cell
                    if ~isempty(cell_list_map{ds}{Tempi(c)}) ;
                        i(r)=c ;
                        r=r+1 ;
                    end
                end
                        
            else
                [temp,i] = intersect(cell_ids{uc}{1},cell_ids{uc}{ds}) ;
            end
        end
        
        for p=1:numPCs ; % for each PC
            projValues{uc}{p}(i,ds) = Trf_dev{uc}{ds} * PCs{uc}(:,end+1-p) ; % dot product
        end
        
        % plot each cells trf in the right spot
        for c=1:size(Trf_short{uc}{ds},1) ; % for each cell ;
            figure(uc)
            subplot(subplot_Y(uc),subplot_X,3+numPCs+i(c))
            plot(Trf_short{uc}{ds}(c,:),Color_list{ds})
            hold on
            axis tight
            if ds == 1 ;
                title(num2str(cell_ids{uc}{ds}(c)))
            end
        end
            
    end
    
    % plot PCA results
    figure(uc)

    % principle component
    subplot(subplot_Y(uc),subplot_X,3)
    plot(PCs{uc}(:,[end-numPCs+1:end]))
    axis tight

    for p=1:numPCs ; % for each PC
        subplot(subplot_Y(uc),subplot_X,(3+p))
        plot([1:numDs],projValues{uc}{p},'*-')
    end
    
end

check = [] ;


% how well are the models at predicting the response?

    


