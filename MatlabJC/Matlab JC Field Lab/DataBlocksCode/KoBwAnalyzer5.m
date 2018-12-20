function ForIgor = KoBwAnalyzer5(DataBlock, DB, Params)

% this function is modified from 'KoBwAnalyzer4'. 

% JC 11/29/16

% data set select
RunAsScript = true ;
if RunAsScript ;
    DB = 56 ; % 
    [DataBlock,Params] = DataBlocks_KO ;
end


% parameters
subplot_X = 10 ;
mapEiFlag = 1 ; % if datablocks should be mapped using map_ei
Color_list = {'k','r','g','b','y'} ; % order of colors for each 
numPs = 2 ; % number of prameters to compare
numDs = length(DataBlock(DB).BwPath) ; % number of data sets

srfProfileX = [0:.1:20] ;

NlFlag = false; % if you want to calc nonlinearities
srfParamCellFlag = false ;
srfParamAllCellFlag = false ;

% compare spatial and temporal receptive fields across data defined sets 
for ds = 1:numDs % for each data set

    % load data
    dataRun = load_data(DataBlock(DB).BwPath{ds}) ;
    dataRun = load_neurons(dataRun) ;

    % how different are the preko and postko response fits
    % spatial rf?
    marks_params.thresh = 3 ;
    dataRun = load_sta(dataRun) ;
    dataRun = load_params(dataRun) ;
    dataRun = load_ei(dataRun, 'all') ;
    dataRun = get_sta_summaries(dataRun, 'all','marks_params', marks_params) ;
    dataRun = get_sta_fits_from_vision(dataRun) ;
    filt_params.radius = 0.75;
    dataRun = get_rfs_filtered(dataRun, 'all', 'filt_params', filt_params,'save_filt_params', 'filt_rf_params', 'save_name','filt_rfs');
    if NlFlag ;
        dataRun = load_java_movie(dataRun, DataBlock(DB).BwMoviePath, dataRun.triggers) ;
        dataRun = get_snls(dataRun, 'all') ;
    end

    if ds == 1 ; % if Master Run
        dataRunMaster = dataRun ;
        
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
            
            Master_i{uc}{ds} = get_cell_indices(dataRun, cell_ids{uc}{ds}) ;
        end
        
    else % if slave run

        if mapEiFlag ; % if using map ei cells
            % map using electrical images
            cell_list_map{ds} = map_ei(dataRunMaster, dataRun) ;
            % cell list struct to mat
            for a=1:length(cell_list_map{ds}) ;
                if isempty(cell_list_map{ds}{a}) ;
                    cell_list_map_mat{ds}(a) = nan ;
                else
                    cell_list_map_mat{ds}(a) = cell_list_map{ds}{a} ;
                end
            end
            
            % cells ids in slave for each UniqueCellType set in master data
            for uc = 1:length(UniqueCellTypes) ;
                Tempi = get_cell_indices(dataRunMaster, UniqueCellTypes{uc}) ;
                cell_ids{uc}{ds} = cell2mat(cell_list_map{ds}(Tempi)) ;
                
                [temp, Master_i{uc}{ds}] = intersect(cell_list_map_mat{ds},cell_ids{uc}{ds}) ; % master cell index
            end
        else % if not using map ei
            for uc = 1:length(UniqueCellTypes) ;
                cell_ids{uc}{ds} = intersect(dataRun.cell_ids, get_cell_ids(dataRunMaster,UniqueCellTypes{uc})) ;
            end
        end
    end
    
    % plot spatial rfs fits
    for uc = 1:length(UniqueCellTypes) ;
        subplot_Y(uc) = ceil((length(cell_ids{uc}{1})+numPs+3)/subplot_X) ;
        
        figure(uc) 
        subplot(2,4,1)
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

    % get trf, srf, nl, and srf profiles
    for uc = 1:length(UniqueCellTypes) ;      
        cell_i{uc}{ds} = get_cell_indices(dataRun, cell_ids{uc}{ds}) ;
        
        Trf{uc}{ds} = nans(length(dataRunMaster.spikes),length(dataRun.stas.time_courses{1})) ;
        Trf_norm{uc}{ds} = nans(length(dataRunMaster.spikes),length(dataRun.stas.time_courses{1})) ;
        srf{uc}{ds} = cell(1,length(dataRunMaster.spikes)) ;
        slope{uc}{ds} = nans(1,length(dataRunMaster.spikes)) ;
        offset{uc}{ds} = nans(1,length(dataRunMaster.spikes)) ;
        srfProfile{uc}{ds} = nans(length(dataRunMaster.spikes),length(srfProfileX)) ;
        
        RfParam{1}{uc}(ds,:) = nans(1,length(dataRunMaster.spikes)) ;
        RfParam{2}{uc}(ds,:) = nans(1,length(dataRunMaster.spikes)) ;
        
        if ~isempty(cell_ids{uc}{ds}) ; % if there are any cells of this type
            for b=1:length(cell_i{uc}{ds}) ; % for each cell
                mi = Master_i{uc}{ds}(b) ;

                Trf{uc}{ds}(mi,:) = flipud(dataRun.stas.time_courses{cell_i{uc}{ds}(b)}) ;
                Trf_norm{uc}{ds}(mi,:) = Trf{uc}{ds}(Master_i{uc}{ds}(b),:)/norm(Trf{uc}{ds}(mi,:)) ;
            
                % srf
                srf{uc}{ds}{mi} = dataRun.stas.filt_rfs{cell_i{uc}{ds}(b)}*...
                    dataRunMaster.stas.polarities{cell_i{uc}{ds}(b)} ; % give the correct pix polarity 

                
                if NlFlag ; 
                    slope{uc}{ds}(mi) = dataRun.stas.snls{cell_i{uc}{ds}(b)}.fit_params.a ;
                    offset{uc}{ds}(mi) = dataRun.stas.snls{cell_i{uc}{ds}(b)}.fit_params.b ;
                    Nl_x{uc}{ds}{mi} = unique(dataRun.stas.snls{cell_i{uc}{ds}(b)}.gen_signal) ;
                    Nl_y{uc}{ds}{mi} = exp(slope{uc}{ds}(mi)*Nl_x{uc}{ds}{mi} + offset{uc}{ds}(mi)) ; 
                end
                
                [srfProfileXTemp, srfProfileTemp] = get_rf_profiles(dataRun,cell_ids{uc}{ds}(b),'normalize','none') ;
                [srfProfileXTemp,ui] = unique(srfProfileXTemp) ;
                srfProfile{uc}{ds}(mi,:) = interp1(srfProfileXTemp,srfProfileTemp(ui),srfProfileX,'nearest','extrap') ; 
            
                if srfParamCellFlag ; % if you want to get the parameters for each cell's srf
                    sta = dataRun.stas.stas{cell_i{uc}{ds}(b)} ; % get the sta
                    temp_marks_sta = significant_stixels(sta, 'thresh', 5, 'time', 'max');

                    if nnz(temp_marks_sta)>0 ; %if there are any sig stixels
                        fit_ins = struct('sig_stixels', temp_marks_sta);

                        fit_temp = fit_sta_sequence(sta,'fit_instructions', fit_ins); % fit the sta of an individual cell
                        RfParam{1}{uc}(ds,mi) = fit_temp.surround_sd_scale ;
                        RfParam_list{1} = 'surround std' ;

                        RfParam{2}{uc}(ds,mi) = fit_temp.surround_amp_scale ;
                        RfParam_list{2} = 'surround amp' ;
                    else
                        RfParam{1}{uc}(ds,mi) = nan ;
                        RfParam{2}{uc}(ds,mi) = nan ;
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
 
for ds = 1:numDs % for each data set
    for uc = 1:length(UniqueCellTypes) ;
        figure(uc)
        subplot(2,4,2)
        plot(srfProfileX,nanmean(srfProfile{uc}{ds}),Color_list{ds})
        hold on
        axis tight
        title('srf profile')

        subplot(2,4,3)
        plot(nanmean(Trf_norm{uc}{ds}),Color_list{ds})
        hold on
        plot(nanmean(Trf_norm{uc}{ds})+nanstd(Trf_norm{uc}{ds}),[Color_list{ds},':'])
        plot(nanmean(Trf_norm{uc}{ds})-nanstd(Trf_norm{uc}{ds}),[Color_list{ds},':'])
        axis tight
        title('temporal rf')

        if NlFlag ;
            subplot(2,4,4)
            plot(Nl_x{1}{1}{Master_i{1}{1}(1)},exp(nanmean(slope{uc}{ds})*Nl_x{1}{1}{Master_i{1}{1}(1)} + nanmean(offset{uc}{ds})),Color_list{ds})
            hold on
        end 
    end
end

 
% find rf params
for uc = 1:length(UniqueCellTypes) ; % for each cell type in the master
    figure(uc)
    for p=1:numPs ; % for each Parameter
        if srfParamCellFlag
            subplot(2,4,4+p)
            plot([1:numDs],RfParam{p}{uc},'*-')
            title(RfParam_list{p})
        end
    end
end

% plot all cells
for uc = 1:length(UniqueCellTypes) ;
    figure
    NumCellsMaster = length(Master_i{uc}{1}) ; % number cells
    for ds = 1:numDs % for each data set
        NumCells = length(Master_i{uc}{ds}) ; % number cells
        for b=1:NumCells ; % for each cell
            mi = Master_i{uc}{ds}(b) ; % master index
            mpi = find(Master_i{uc}{1}==mi) ; % plot inex
            
            subplot(numDs+2,NumCellsMaster,mpi+NumCellsMaster*(ds-1)) ; % srf
            norm_rf = norm_image(srf{uc}{ds}{mi});
            imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
            colormap(brewermap([],'RdBu'))
            caxis([0,1]) 
            set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
            title(num2str(mi))
            
            subplot(numDs+2,NumCellsMaster,mpi+NumCellsMaster*numDs) ; % trf
            plot(Trf_norm{uc}{ds}(mi,:),Color_list{ds})
            hold on
            
            subplot(numDs+2,NumCellsMaster,mpi+NumCellsMaster*(numDs+1)) ; % nl
        end
    end
end
        


clear dataRun*
save(['/Users/jcafaro/Desktop/Temp/Matfiles/KoBwAnalyzer4Db',num2str(DB)],'-v7.3') ;


% how well are the models at predicting the response?

    


