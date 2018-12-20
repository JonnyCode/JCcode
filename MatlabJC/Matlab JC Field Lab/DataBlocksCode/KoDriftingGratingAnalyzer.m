function ForIgor = KoDriftingGratingAnalyzer(DataBlock, DB, Params)

% this function will analyze drfiting grating responses from array data
% before and after cells are knocked out (KO)

% This function assumes that cells are only typed in first data set from white noise stim.  Cells
% that are not mapped onto this are thrown out.

% JC 2/25/15
%DB=9; % TEMP

% parameters
mapEiFlag = DataBlock(DB).mapEiFlag ; % if datablocks should be mapped using map_ei
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
numPCs = 3 ; % number of principle components
numDs = length(DataBlock(DB).DataPath) ; % number of data sets
subplot_X = 5 ;
TrialTrigInterval = 10 ; % sec
InterTrialTime = 2 ; % sec
StimTransTime = 1 ; % sec
spatial_plot = 2 ; % mat index of stim_group
temporal_plot = 2 ; % mat index of stim_group
ds_thresh = 0.5 ; % threshold of ds index where cells will be noted as ds

% compare spatial and temporal receptive fields across data defined sets 
for ds = 1:numDs % for each data set

    % load data
    dataRun = load_data(DataBlock(DB).DataPath{ds}) ;
    dataRun = load_neurons(dataRun) ;

    if ds == 1 ; % if Master Run
        
        % load data
        marks_params.thresh = 3 ;
        dataRun = load_sta(dataRun) ;
        dataRun = load_params(dataRun) ;
        dataRun = get_sta_summaries(dataRun, 'all','marks_params', marks_params) ;
        dataRun = get_sta_fits_from_vision(dataRun) ;

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
        
        % load stimulus
        dataRun.names.stimulus_path = DataBlock(DB).StimPath{ds} ;
        dataRun = load_stim(dataRun,'user_defined_trigger_interval', TrialTrigInterval) ;
        
        % stim trials with the same spatial and temporal periods
        for sp=1:length(dataRun.stimulus.params.SPATIAL_PERIOD) ; % for each spatial period
            for tp=1:length(dataRun.stimulus.params.TEMPORAL_PERIOD) ; % for each temporal period
                stim_groups_trials{sp,tp} = [] ; % empty group
                for t=1:length(dataRun.stimulus.trials) ; % for each trial
                    if dataRun.stimulus.trials(t).SPATIAL_PERIOD == dataRun.stimulus.params.SPATIAL_PERIOD(sp) ...
                            & dataRun.stimulus.trials(t).TEMPORAL_PERIOD == dataRun.stimulus.params.TEMPORAL_PERIOD(tp) ;
                        stim_groups_trials{sp,tp} = [stim_groups_trials{sp,tp}(:)',t] ; % this trial included
                    end
                end
            end
        end
        
        % trigger times for each
        for sp=1:length(dataRun.stimulus.params.SPATIAL_PERIOD) ; % for each spatial period
            for tp=1:length(dataRun.stimulus.params.TEMPORAL_PERIOD) ; % for each temporal period
                stim_group_triggers{sp,tp} = [] ;
                tempTrigs = dataRun.stimulus.triggers(stim_groups_trials{sp,tp}) ;
                for a = 1:length(tempTrigs) ; % for each trigger
                    stim_group_triggers{sp,tp} = [stim_group_triggers{sp,tp}(:); dataRun.triggers(dataRun.triggers>=tempTrigs(a)+StimTransTime & ...
                        dataRun.triggers < tempTrigs(a)+TrialTrigInterval-InterTrialTime)] ; 
                end
            end
        end
        
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

    if ds~=1 ; % if its not the first data set and has a grating response
        
        % check direction selectivity
        for uc = 1:length(UniqueCellTypes) ;
            [NumSpikesCell_temp, StimComb_temp] = get_spikescellstim(dataRun, cell_ids{uc}{ds}, 0);
            [mag  dsindex_temp  magmax  magave  angle  rho  theta  num  U  V ] = dscellanalysis(NumSpikesCell_temp, StimComb_temp);
            
            dsi{uc}{ds}= dsindex_temp{temporal_plot,spatial_plot} ;
        end

 
        % plot step responses
        for uc = 1:length(UniqueCellTypes) ;
            subplot_Y = ceil((2 + length(cell_ids{uc}{1}))/subplot_X) ;
            figure(uc) 
            title(UniqueCellTypes{uc})
            
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

            if ~isempty(cell_ids{uc}) ;
                cell_i{uc}{ds} = get_cell_indices(dataRun, cell_ids{uc}{ds}) ;
                for b=1:length(cell_i{uc}{ds}) ;
                    % [psthTemp, binsTemp] =
                    % get_psth(dataRun.spikes{cell_i{uc}{ds}(b)}, dataRun.triggers(:)) ; 
                    % [psthTemp, binsTemp] = get_psth(dataRun.spikes{cell_i{uc}{ds}(b)}, dataRun.stimulus.triggers(stim_groups_trials{spatial_plot,temporal_plot})) ;
                    %[psthTemp, binsTemp] = get_psth(dataRun.spikes{cell_i{uc}{ds}(b)}, dataRun.stimulus.triggers) ;
                    [psthTemp, binsTemp] = get_psth(dataRun.spikes{cell_i{uc}{ds}(b)}, stim_group_triggers{spatial_plot,temporal_plot}(:)','stop',diff(stim_group_triggers{spatial_plot,temporal_plot}(1:2)),...
                        'bin_size', .05) ;
                    psth{uc}{ds}(b,:) = psthTemp ;
                    bins{uc}{ds}(b,:) = binsTemp ;
                
                    subplot(subplot_Y,subplot_X,i(b)+2)
                    plot(binsTemp,psthTemp,Color_list{ds})
                    axis tight
                    if dsi{uc}{ds}(b)>ds_thresh ;
                        text(.9,1-ds/10,'*','units','normalized','color',Color_list{ds},'fontsize', 15)
                    end
                    hold on
                end
                
                subplot(subplot_Y,subplot_X,2)
                plot(bins{uc}{ds}(1,:),mean(psth{uc}{ds}),Color_list{ds})
                hold on
                plot(bins{uc}{ds}(1,:),mean(psth{uc}{ds})+std(psth{uc}{ds}),[Color_list{ds},':'])
                plot(bins{uc}{ds}(1,:),mean(psth{uc}{ds})-std(psth{uc}{ds}),[Color_list{ds},':'])
                axis tight
                
            end
        end
    end
end

dsi_mat = [] ;
for uc = 1:length(UniqueCellTypes) ;
    for ds = 1:numDs % for each data set
        dsi_mat = [dsi_mat;dsi{uc}{ds}(:)] ;
    end
end
    

check = [] ;

