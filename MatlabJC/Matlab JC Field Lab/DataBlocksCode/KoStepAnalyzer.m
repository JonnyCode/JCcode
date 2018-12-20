function ForIgor = KoStepAnalyzer(DataBlock, DB, Params)

% this function will analyze full field light step response from array data
% before and after cells are knocked out (KO)

% This function assumes that cells are only typed in first data set from white noise stim.  Cells
% that are not mapped onto this are thrown out.

% JC 2/13/15
%DB=9; % TEMP

% parameters
mapEiFlag = DataBlock(DB).mapEiFlag ; % if datablocks should be mapped using map_ei
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
numPCs = 3 ; % number of principle components
numDs = length(DataBlock(DB).DataPath) ; % number of data sets
subplot_X = 5 ;

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

    if ds~=1 ; % if its not the first data set and has a step response        
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
                    [psthTemp, binsTemp] = get_psth(dataRun.spikes{cell_i{uc}{ds}(b)}, dataRun.triggers([4:4:end])) ;
                    psth{uc}{ds}(b,:) = psthTemp ;
                    bins{uc}{ds}(b,:) = binsTemp ;
                
                    subplot(subplot_Y,subplot_X,i(b)+2)
                    plot(binsTemp,psthTemp,Color_list{ds})
                    axis tight
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

check = [] ;

