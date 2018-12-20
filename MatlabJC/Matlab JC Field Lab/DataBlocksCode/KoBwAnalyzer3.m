function ForIgor = KoBwAnalyzer3(DataBlock, DB, Params)

% this function is modified from 'KoBwAnalyzer2'.  Got rid of PCA, added srf and trf parameter and nl comparisons.  

% JC 4/1/16

% parameters
subplot_X = 10 ;
mapEiFlag = DataBlock(DB).mapEiFlag ; % if datablocks should be mapped using map_ei
Color_list = {'k','r','g','y','b'} ; % order of colors for each 
numPs = 2 ; % number of prameters to compare
numDs = length(DataBlock(DB).DataPath) ; % number of data sets

NlFlag = true ; % if you want to calc nonlinearities
srfParamCellFlag = false ;
srfParamAllCellFlag = false ;

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
    filt_params.radius = 0.75;
    dataRun = get_rfs_filtered(dataRun, 'all', 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');
    if NlFlag ;
        dataRun = load_java_movie(dataRun, DataBlock(DB).BwMoviePath, dataRun.triggers) ;
        dataRun = get_snls(dataRun, 'all') ;
    end

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
        subplot_Y(uc) = ceil((length(cell_ids{uc}{1})+numPs+3)/subplot_X) ;
        
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
        if ~isempty(cell_ids{uc}{ds}) ; % if there are any cells of this type
            cell_i{uc}{ds} = get_cell_indices(dataRun, cell_ids{uc}{ds}) ;
            for b=1:length(cell_i{uc}{ds}) ;
                Trf{uc}{ds}(b,:) = flipud(dataRun.stas.time_courses{cell_i{uc}{ds}(b)}) ;
                Trf_norm{uc}{ds}(b,:) = Trf{uc}{ds}(b,:)/norm(Trf{uc}{ds}(b,:)) ;
            
                if NlFlag ; 
                    slope{uc}{ds}(b) = dataRun.stas.snls{cell_i{uc}{ds}(b)}.fit_params.a ;
                    offset{uc}{ds}(b) = dataRun.stas.snls{cell_i{uc}{ds}(b)}.fit_params.b ;
                    Nl_x{uc}{ds}{b} = unique(dataRun.stas.snls{cell_i{uc}{ds}(b)}.gen_signal) ;
                    Nl_y{uc}{ds}{b} = exp(slope{uc}{ds}(b)*Nl_x{uc}{ds}{b} + offset{uc}{ds}(b)) ; 
                end
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
    
    % find srf params
    for uc = 1:length(UniqueCellTypes) ; % for each cell type in the master
         if ~isempty(cell_ids{uc}) ; % if there are cells 
            for c = 1:length(cell_ids{uc}{ds}) ; % for each cell
                ci = get_cell_indices(dataRun,cell_ids{uc}{ds}(c)) ; % cell index
                
                [srfProfileXTemp, srfProfileTemp] = get_rf_profiles(dataRun,cell_ids{uc}{ds}(c),'normalize','none') ;
                [srfProfileXTemp,ui] = unique(srfProfileXTemp) ;
                srfProfileX = [0:.1:20] ;
                srfProfile{uc}{ds}(c,:) = interp1(srfProfileXTemp,srfProfileTemp(ui),srfProfileX,'nearest','extrap') ; 
                
                if srfParamCellFlag ; % if you want to get the parameters for each cell's srf
                    sta = dataRun.stas.stas{ci} ; % get the sta
                    temp_marks_sta = significant_stixels(sta, 'thresh', 5, 'time', 'max');

                    if nnz(temp_marks_sta)>0 ; %if there are any sig stixels
                        fit_ins = struct('sig_stixels', temp_marks_sta);

                        fit_temp = fit_sta_sequence(sta,'fit_instructions', fit_ins); % fit the sta of an individual cell
                        param{1}{uc}(ds,c) = fit_temp.surround_sd_scale ;
                        param_list{1} = 'surround std' ;

                        param{2}{uc}(ds,c) = fit_temp.surround_amp_scale ;
                        param_list{2} = 'surround amp' ;
                    else
                        param{1}{uc}(ds,c) = nan ;
                        param{2}{uc}(ds,c) = nan ;
                    end
                end
                if srfParamAllCellFlag ; % if you want to average across all cells and then get params

    %                 sta_std = robust_std(sta(:),3) ;
    %                 sta_snr = sta/sta_std ;
    %                 sta_ucWm = sta_ucWm = 
                end          
            end
        end
    end
end


% find rf params
for uc = 1:length(UniqueCellTypes) ; % for each cell type in the master
    figure(uc)
    for p=1:numPs-1 ; % for each Parameter
        if srfParamCellFlag
            subplot(subplot_Y(uc),subplot_X,(2+p))
            plot([1:numDs],param{p}{uc},'*-')
            title(param_list{p})
        end
    end
    
    subplot(subplot_Y(uc),subplot_X,(2+numPs))
    plot(srfProfileX,mean(srfProfile{uc}{1}),Color_list{1})
    hold on
    plot(srfProfileX,mean(srfProfile{uc}{2}),Color_list{2})
    plot(srfProfileX,mean(srfProfile{uc}{3}),Color_list{3})
    title('srf profile')
end

figure
for uc = 1:length(UniqueCellTypes) ;
    subplot(ceil(length(UniqueCellTypes)/4),4,uc)
    plot(srfProfileX,mean(srfProfile{uc}{1}),Color_list{1})
    hold on
    plot(srfProfileX,mean(srfProfile{uc}{2}),Color_list{2})
    plot(srfProfileX,mean(srfProfile{uc}{3}),Color_list{3})
    title(UniqueCellTypes{uc})
end

check = [] ;


% how well are the models at predicting the response?

    


