function ForIgor = BwAnalyzerV2(DataBlock, DB, Params)

% this function will analyze binary white stimuli response from array data

% This function assumes that cells are only typed in first data set.  Cells
% that are not mapped onto this are thrown out.

% JC 5/4/15
%DB=3; % TEMP

% parameters
NlFlag = false ; % get NLs
mapEiFlag = DataBlock(DB).mapEiFlag ; % if datablocks should be mapped using map_ei
Color_list = {'k','r','g','y','b'} ; % order of colors for each 
numDs = length(DataBlock(DB).BwPath) ; % number of data sets

% compare spatial and temporal receptive fields across data defined sets 
for ds = 1:numDs % for each data set

    % load data
    dataRun = load_data(DataBlock(DB).BwPath{ds}) ;
    dataRun = load_neurons(dataRun) ;

    marks_params.thresh = 3 ;
    dataRun = load_sta(dataRun) ;
    dataRun = load_params(dataRun,'cell_type_depth', 5) ;
    dataRun = get_sta_summaries(dataRun, 'all','marks_params', marks_params) ;
    dataRun = get_sta_fits_from_vision(dataRun) ;
    if NlFlag ;
        dataRun = load_java_movie(dataRun, DataBlock(DB).BwMoviePath{ds}, dataRun.triggers) ;
        dataRun = get_snls(dataRun, 'all') ;
    end
    
    dataRun.stimulus.monitor_refresh = 60.35 ; % the monitor_refresh rate is not accurate
    trf_time(ds,:) = -[0:dataRun.stas.depth-1]*dataRun.stimulus.interval/dataRun.stimulus.monitor_refresh ; % time axis of temporal receptive field
    PixPerStix = dataRun.stimulus.stixel_height ;

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
        dataRunMaster = load_data(DataBlock(DB).BwPath{1}) ;
        dataRunMaster = load_neurons(dataRunMaster) ;
        dataRunMaster = load_ei(dataRunMaster, 'all') ;
        dataRunMaster = load_params(dataRunMaster,'cell_type_depth', 5) ;
        
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

    % spatial rfs
    for uc = 1:length(UniqueCellTypes) ;
        if ~isempty(cell_ids{uc}) ; % if there are cells           
            coord_tform = coordinate_transform(dataRun,'sta') ;
            
            for c = 1:length(cell_ids{uc}{ds}) ; % for each cell
                ci = get_cell_indices(dataRun,cell_ids{uc}{ds}(c)) ; % cell index
                ctr = dataRun.stas.fits{ci}.mean * PixPerStix ;
                rad = dataRun.stas.fits{ci}.sd * PixPerStix ;
                angle = dataRun.stas.fits{ci}.angle ;
                [x,y] = drawEllipse([ctr rad angle]) ;
                if ~any(isnan([x,y])) ;
                    [X{uc}{ds}(c,:),Y{uc}{ds}(c,:)] = tformfwd(coord_tform, x, y) ;
                end                
            end         
        end
    end

    % temporal rf and NL
    for uc = 1:length(UniqueCellTypes) ;        
        if ~isempty(cell_ids{uc}{ds}) ;
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
            
            Trf_norm_mean{uc}{ds} = mean(Trf_norm{uc}{ds}) ;
            Trf_norm_std{uc}{ds} = std(Trf_norm{uc}{ds}) ;
            
            Trf_mean{uc}{ds} = mean(Trf{uc}{ds}) ;
            Trf_std{uc}{ds} = std(Trf{uc}{ds}) ;
            
            if NlFlag ;
                Nl_x_mean{uc}{ds} = unique(cell2mat(Nl_x{uc}{ds}')) ;
                for b=1:length(cell_i{uc}{ds}) ;
                     y_temp(b,:) = exp(slope{uc}{ds}(b)*Nl_x_mean{uc}{ds} + offset{uc}{ds}(b)) ;   
                end
                Nl_y_mean{uc}{ds} = mean(y_temp,1) ;
                clear y_temp
            end
        end
    end
end

% get parameters for STA time course
for ds = 1:numDs % for each data set
    for uc = 1:length(UniqueCellTypes) ;  % for each cell type
        for i = 1:length(cell_i{uc}{ds}) ;
            [tcFit{uc}{ds}(i,:), final_params(i,:)] = fit_time_course(fliplr(Trf{uc}{ds}(i,:))', 'verbose', false) ;
            tcFit_norm{uc}{ds}(i,:) = tcFit{uc}{ds}(i,:)/norm(tcFit{uc}{ds}(i,:)) ;
        end
            
        tcParams{uc}{ds} = tc_params_finder(fliplr(tcFit{uc}{ds}),trf_time(ds,:)) ;
        
        [firstPeakTimeHist{uc}{ds}, firstPeakTimeHistX{uc}{ds}] = hist(tcParams{uc}{ds}.firstPeakt,[-1:.01:0]) ;
        firstPeakTime_mean(uc,ds) = nanmean(tcParams{uc}{ds}.firstPeakt) ;
        firstPeakTime_sem(uc,ds) = nanstd(tcParams{uc}{ds}.firstPeakt)/sqrt(sum(~isnan(tcParams{uc}{ds}.firstPeakt))) ;
        
        [firstPeakHist{uc}{ds}, firstPeakHistX{uc}{ds}] = hist(tcParams{uc}{ds}.firstPeak,[-1:.01:1]) ;
        firstPeak_mean(uc,ds) = nanmean(tcParams{uc}{ds}.firstPeak) ;
        firstPeak_sem(uc,ds) = nanstd(tcParams{uc}{ds}.firstPeak)/sqrt(sum(~isnan(tcParams{uc}{ds}.firstPeak))) ;
        
        [DoTHist{uc}{ds}, DoTHistX{uc}{ds}] = hist(tcParams{uc}{ds}.DoT,[0:.25:10]) ;
        DoT_mean(uc,ds) = nanmean(tcParams{uc}{ds}.DoT) ;
        DoT_sem(uc,ds) = nanstd(tcParams{uc}{ds}.DoT)/sqrt(sum(~isnan(tcParams{uc}{ds}.DoT))) ;  
        
    end
end

% figures
for ds = 1:numDs % for each data set
    for uc = 1:length(UniqueCellTypes) ;
        figure(uc) 
        
        if ~isempty(cell_ids{uc}) ; % if there are cells 
            % spatial receptive field
            subplot(3,3,1:3)
            title(UniqueCellTypes{uc})
            hold on
            for c = 1:length(cell_ids{uc}{ds}) ; % for each cell
                if ~any(isnan([X{uc}{ds}(c,:),Y{uc}{ds}(c,:)])) ;
                    plot(X{uc}{ds}(c,:),Y{uc}{ds}(c,:),'Color','r')
                    hold on
                    plot_array(dataRun,'k',coord_tform,dataRun.piece.array_id) % array outline
                end
            end         
        
            % temporal receptive field and nl
            figure(uc) 
    
            % for each cell
            for c = 1:length(cell_ids{uc}{ds}) ; % for each cell
                if ~any(isnan([X{uc}{ds}(c,:),Y{uc}{ds}(c,:)])) ;
                    subplot(3,3,4:5)
                    plot(trf_time(ds,:), Trf_norm{uc}{ds}(c,:),[Color_list{ds},':'],'LineWidth', .05)
                    hold on
                    
                    if NlFlag ;
                        subplot(3,3,6)
                        plot(Nl_x{uc}{ds}{c}, Nl_y{uc}{ds}{c},[Color_list{ds},':'],'LineWidth', .05)
                        hold on
                        axis tight
                    end
                end
            end
            
            % averages
            subplot(3,3,4:5)
            plot(trf_time(ds,:), Trf_norm_mean{uc}{ds},Color_list{ds}, 'LineWidth', 2)
            hold on
            plot(trf_time(ds,:), Trf_norm_mean{uc}{ds}+Trf_norm_std{uc}{ds},[Color_list{ds},'--'], 'LineWidth', 2)
            plot(trf_time(ds,:), Trf_norm_mean{uc}{ds}-Trf_norm_std{uc}{ds},[Color_list{ds},'--'], 'LineWidth', 2)
            axis tight
            
            if NlFlag ;
                subplot(3,3,6)
                plot(Nl_x_mean{uc}{ds}', Nl_y_mean{uc}{ds}, 'LineWidth', 2)
            end

            % histograms
            subplot(3,3,7)
            plot(firstPeakTimeHistX{uc}{ds},firstPeakTimeHist{uc}{ds},Color_list{ds})
            xlabel('peak1 t')
            hold on

            subplot(3,3,8)
            plot(firstPeakHistX{uc}{ds},firstPeakHist{uc}{ds},Color_list{ds})
            xlabel('peak1 amp')
            hold on

            subplot(3,3,9)
            plot(DoTHistX{uc}{ds},DoTHist{uc}{ds},Color_list{ds})
            xlabel('DoT')
            hold on
            
            title(UniqueCellTypes{uc})      
        end  
    end
end




% For Igor

ForIgor = struct() ;

Example_Srf_uc = 3 ;
Example_Srf_ds = 1 ;


for ds = 1:numDs % for each data set
    for uc = 1:length(UniqueCellTypes) ;  % for each cell type
        
        VecName = ['Trf','Uc',num2str(uc),'Ds',num2str(ds),'Db',num2str(DB)] ;
        ForIgor = setfield(ForIgor,VecName,mean(Trf_norm{uc}{ds})) ; % time course

        VecName = ['TrfX','Ds',num2str(ds),'Db',num2str(DB)] ;
        ForIgor = setfield(ForIgor,VecName,trf_time(ds,:)) ; % time course
    end
end

for c = 1:length(cell_ids{Example_Srf_uc}{Example_Srf_ds}) ; % for each cell
    if ~any(isnan([X{Example_Srf_uc}{Example_Srf_ds}(c,:),Y{Example_Srf_uc}{Example_Srf_ds}(c,:)])) ;
        TempX(:,c) = X{Example_Srf_uc}{Example_Srf_ds}(c,:) ;
        TempY(:,c) = Y{Example_Srf_uc}{Example_Srf_ds}(c,:) ;
    end
end   
VecName = ['SrfX','Uc',num2str(Example_Srf_uc),'Ds',num2str(Example_Srf_ds),'Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,TempX) ; % spatial X

VecName = ['SrfY','Uc',num2str(Example_Srf_uc),'Ds',num2str(Example_Srf_ds),'Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,TempY) ; % spatial Y


check = [] ;








        
    
    
    


