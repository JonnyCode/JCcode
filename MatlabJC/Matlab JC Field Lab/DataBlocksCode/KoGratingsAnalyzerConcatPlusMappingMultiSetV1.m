function ForIgor = KoGratingsAnalyzerConcatPlusMappingMultiSetV1(DataBlock, DB, Params)

% this function is adapted from
% 'KoReverseGratingsAnalyzerConcatPlusMappingMultiSetV1'  

% JC 5/16/2016 

% data set select
RunAsScript = true ;
if RunAsScript ;
    DB = 57 ; % 56,57
    [DataBlock,Params] = DataBlocks_KO ;
end


% parameters
Color_list = {'k','r','b','g','y','v'} ; % order of colors for each 
TrialTrigInterval = 10 ; % sec
InterTrialTime = 2 ; % sec
StimTransTime = .1 ; % (sec) skip cycles before this time after starting new stimulus
MicronPerPix = 4 ; % (um/pix)
EsemFlag = false; % calc psth esemble
FigureFlag = false ; % if you want to run figures
MseRelWashThreshold = 1 ; % (unitless) below this point response is "stable"
MseRelDrugDivWashThreshold = 1 ; % (unitless) above this point response has big drug effect 
LfSpThresh = 300 ; % (pix) spatial periods equal and above are termed "low frequency"

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

% save spatial and temporal rf so can delete dataRunMaster later
for a=1:length(dataRunMaster.stas.filt_rfs) ;
    srf{a} = dataRunMaster.stas.filt_rfs{a}*dataRunMaster.stas.polarities{a} ; % give the correct pix polarity 
end
trf = dataRunMaster.stas.time_courses ; 

dataRunMaster.stimulus.monitor_refresh = 60.35 ; % the monitor_refresh rate is not accurate
trf_time = -[0:dataRunMaster.stas.depth-1]*dataRunMaster.stimulus.interval/dataRunMaster.stimulus.monitor_refresh ; 

% get cell types 
for a = 1:length(dataRunMaster.cell_types) ;
    celltypes{a} = dataRunMaster.cell_types{a}.name ;
end

UniqueCellTypes = unique(celltypes) ;
if isempty(UniqueCellTypes{1}) ;
    UniqueCellTypes = UniqueCellTypes(2:end) ;
end
        
% load stimulus and get triggers for each part of concat data set
for dset=1:length(DataBlock(DB).DgPath) ;
    % load data 
    dataRun = load_data(DataBlock(DB).DgPath{dset}) ;
    dataRun = load_neurons(dataRun) ;

    % load stim
    dataRun.names.stimulus_path = [DataBlock(DB).DgPath{dset}(1:end-15),'stimuli/s',DataBlock(DB).DgPath{dset}(end-1:end),'.txt'] ;
    dataRun = load_stim(dataRun,'user_defined_trigger_interval', TrialTrigInterval) ;

    setTrialNum(dset) = length(dataRun.stimulus.trials) ; % number of trials for this stim set
    setTrialNumShift(dset) = sum(setTrialNum(1:dset))-setTrialNum(dset) ; % number of trials in preceding sets
    setTrigNum(dset) = length(dataRun.triggers) ; % number of triggers in this data set
    setTrigNumShift(dset) = sum(setTrigNum(1:dset))-setTrigNum(dset) ; % number of triggers in preceding sets
    
    dataRun.stimulus.params.RGB = unique(dataRun.stimulus.params.RGB) ; % not sure why this is not already done as for TP and SP

    l_sp = length(dataRun.stimulus.params.SPATIAL_PERIOD) ;
    l_tp = length(dataRun.stimulus.params.TEMPORAL_PERIOD) ;
    l_cntrst = length(dataRun.stimulus.params.RGB) ;
    
    % stim trials with the same spatial and temporal periods
    for sp=1:l_sp ; % for each spatial period
        for tp=1:l_tp ; % for each temporal period
            for cntrst=1:l_cntrst ; % for each contrast
                stim_groups_trials{sp,tp,cntrst} = [] ; % empty group
                for t=1:length(dataRun.stimulus.trials) ; % for each trial
                    if dataRun.stimulus.trials(t).SPATIAL_PERIOD == dataRun.stimulus.params.SPATIAL_PERIOD(sp) ...
                            & dataRun.stimulus.trials(t).TEMPORAL_PERIOD == dataRun.stimulus.params.TEMPORAL_PERIOD(tp)...
                            & dataRun.stimulus.trials(t).RGB(1) == dataRun.stimulus.params.RGB(cntrst) ;
                        stim_groups_trials{sp,tp,cntrst} = [stim_groups_trials{sp,tp,cntrst}(:)',t] ; % this trial included
                    end
                end
            end
        end
    end

    % trigger times for each (THIS ASSUMES the sp, tp,cntrst are all the
    % same between data sets, and only the trigg times change)
    for sp=1:l_sp ; % for each spatial period
        for tp=1:l_tp ; % for each temporal period
            for cntrst=1:l_cntrst ; % for each contrast
                stim_group_triggers_i{dset}{sp,tp, cntrst} = [] ;
                tempTrigs = dataRun.stimulus.triggers(stim_groups_trials{sp,tp,cntrst}) ; % the first trigger for each stimulus 
                for a = 1:length(tempTrigs) ; % for each first trigger
                    stim_group_triggers_i{dset}{sp,tp,cntrst} = [stim_group_triggers_i{dset}{sp,tp,cntrst}(:); find(dataRun.triggers>=tempTrigs(a)+StimTransTime & ...
                        dataRun.triggers <= tempTrigs(a)+TrialTrigInterval-InterTrialTime)+setTrigNumShift(dset)] ; % find the triggers within that set
                end
            end
        end
    end
    % stim params
    SpatialPeriod = dataRun.stimulus.params.SPATIAL_PERIOD ;
    TemporalPeriod = dataRun.stimulus.params.TEMPORAL_PERIOD ;
    Contrast = dataRun.stimulus.params.RGB ;
    
    clear dataRun stim_groups_trials
end

Spatial_frequency = 1./(SpatialPeriod*MicronPerPix) ; % cycles/um

% load concatinated data
dataRun = load_data(DataBlock(DB).DgConcat{1}) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_ei(dataRun, 'all') ; % load electrical images

l_cells = length(dataRun.spikes) ; % number of cells
% map using electrical images
cell_list_map = map_ei(dataRunMaster, dataRun) ;

% cells ids in slave for each UniqueCellType set in master data
for uc = 1:length(UniqueCellTypes) ;
    Tempi = get_cell_indices(dataRunMaster, UniqueCellTypes{uc}) ;
    cell_ids{uc} = cell2mat(cell_list_map(Tempi)) ;
    if ~isempty(cell_ids{uc}) ;
        cell_i{uc} = get_cell_indices(dataRun, cell_ids{uc}) ;
    else
        cell_i{uc} = [] ;
    end
end

% prep cells and get trigger times
for dset=1:length(DataBlock(DB).DgPath) ;
    for sp=1:l_sp ; % for each spatial period
        for cntrst=1:l_cntrst ; % for each contrast
            for tp=1:l_tp ; % for each temporal period          
                stim_group_triggers{dset}{sp,tp,cntrst}= dataRun.triggers(stim_group_triggers_i{dset}{sp,tp,cntrst}) ; % trigger times
                stpPnt{dset}(sp,tp) = ceil(min(diff(stim_group_triggers{dset}{sp,tp,1}'))*dataRun.sampling_rate) ; % length of psth (prevents single index errors)
                psth{dset}{tp}{cntrst}{sp} = nans(l_cells,stpPnt{dset}(sp,tp)) ;
            end
        end
    end

    for tp=1:l_tp ; % for each temporal period
        for sp=1:l_sp ; % for each spatial period
            for cells = 1:l_cells ; % for each cell
                for cntrst=1:l_cntrst ; % for each phase
                    if ~isempty(stim_group_triggers{dset}{sp,tp,cntrst}) ;
                       
                        [psthTemp,PsthTimeTemp] = get_smooth_psth(dataRun.spikes{cells},stim_group_triggers{dset}{sp,tp,cntrst}','stop',min(diff(stim_group_triggers{dset}{sp,tp,cntrst}')),'tailTime',.01) ;
                        psth{dset}{tp}{cntrst}{sp}(cells,:) = psthTemp(1:stpPnt{dset}(sp,tp)) ;
                        PsthTime{tp} = PsthTimeTemp(1:stpPnt{dset}(sp,tp)) ;

                        [powerspec_xvalues, mean_powerspec] = PowerSpectrumFinder(psth{dset}{tp}{cntrst}{sp}(cells,:),dataRun.sampling_rate) ;
                        psth_varF1{dset}{tp}{sp}(cells,cntrst) = mean_powerspec(2) ;
                        psth_varF2{dset}{tp}{sp}(cells,cntrst) = mean_powerspec(3) ;
                        
                        psth_var{dset}{tp}{sp}(cells,cntrst) = var(psth{dset}{tp}{cntrst}{sp}(cells,:)) ; 
                        psth_mean{dset}{tp}{sp}(cells,cntrst) = mean(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                        psth_range{dset}{tp}{sp}(cells,cntrst) = range(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;

                        psth_peak{dset}{tp}{sp}(cells,cntrst) = max(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                        psth_peakDivMean{dset}{tp}{sp}(cells,cntrst) = psth_peak{dset}{tp}{sp}(cells,cntrst)/psth_mean{dset}{tp}{sp}(cells,cntrst) ;
                        psth_duty{dset}{tp}{sp}(cells,cntrst) = sum(psth{dset}{tp}{cntrst}{sp}(cells,:)>(psth_peak{dset}{tp}{sp}(cells,cntrst)/2))/length(psthTemp) ; % fraction of time above 50%  
                    end
                end
            end        
        end
    end
end

% averaged over contrast
for dset=1:length(DataBlock(DB).DgPath) ;
    for tp=1:l_tp ; % for each temporal period
        for sp=1:l_sp ; % for each spatial period
            for cells = 1:l_cells ; % for each cell
                 psth_varF1_cntrstMean{dset}{tp}(cells,sp) = mean(psth_varF1{dset}{tp}{sp}(cells,:)) ;
                 psth_varF2_cntrstMean{dset}{tp}(cells,sp) = mean(psth_varF2{dset}{tp}{sp}(cells,:)) ;
                 psth_var_cntrstMean{dset}{tp}(cells,sp) = mean(psth_var{dset}{tp}{sp}(cells,:)) ;
            end
        end
    end
end

% change comared to 1st dset
for dset=2:length(DataBlock(DB).DgPath) ;
    for tp=1:l_tp ; % for each temporal period
         psth_varF1_cntrstMean_DivCntrl{dset}{tp} = psth_varF1_cntrstMean{dset}{tp}./psth_varF1_cntrstMean{1}{tp} ;
         psth_varF2_cntrstMean_DivCntrl{dset}{tp} = psth_varF2_cntrstMean{dset}{tp}./psth_varF2_cntrstMean{1}{tp} ;
         psth_var_cntrstMean_DivCntrl{dset}{tp} = psth_var_cntrstMean{dset}{tp}./psth_var_cntrstMean{1}{tp} ;
    end
end

% cell list struct to mat
for a=1:length(cell_list_map) ;
    if isempty(cell_list_map{a}) ;
        cell_list_map_mat(a) = nan ;
    else
        cell_list_map_mat(a) = cell_list_map{a} ;
    end
end

% temporal receptive fields by cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    if length(cell_ids{uc})>0 ; 
        for a=1:length(cell_ids{uc}) ; % for each cell of this type
            Mi = find(cell_list_map_mat==cell_ids{uc}(a)) ; % master cell index
            Mi = Mi(1) ; % some slaves seem to have two masters?
            if isempty(trf{Mi})
                disp(['empty trf: ',num2str(Mi)])
                trf_uc{uc}(a,:) = nan(size(trf_time)) ;
            else
                trf_uc{uc}(a,:) = trf{Mi}' ;
            end
        end
        trf_uc_mean(uc,:) = nanmean(trf_uc{uc},1) ;
    else
        trf_uc_mean(uc,:) = nan(size(trf_time)) ;
    end
end

% average tuning curves across identified unique cell types
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    for tp=1:l_tp ; % for each temporal period
        for dset=1:length(DataBlock(DB).DgPath) ;
            clear tempMat*
            
            lc = length(cell_i{uc}) ; % number of cell of that type
            
            if lc>0 ; % if there are identified cells of this type
                for cells=1:lc ; % for each cell of this type  
                    tempMat(cells,:) = psth_var_cntrstMean{dset}{tp}(cell_i{uc}(cells),:) ;
                    tempMatF1(cells,:) = psth_varF1_cntrstMean{dset}{tp}(cell_i{uc}(cells),:) ;
                    tempMatF2(cells,:) = psth_varF2_cntrstMean{dset}{tp}(cell_i{uc}(cells),:) ;            
                end
                psth_var_cntrstMean_CellMean{dset}{tp}(uc,:) = mean(tempMat,1) ;
                psth_var_cntrstMean_CellSem{dset}{tp}(uc,:) = std(tempMat,[],1)/sqrt(lc) ; 

                psth_varF1_cntrstMean_CellMean{dset}{tp}(uc,:) = mean(tempMatF1,1) ;
                psth_varF1_cntrstMean_CellSem{dset}{tp}(uc,:) = std(tempMatF1,[],1)/sqrt(lc) ;

                psth_varF2_cntrstMean_CellMean{dset}{tp}(uc,:) = mean(tempMatF2,1) ;
                psth_varF2_cntrstMean_CellSem{dset}{tp}(uc,:) = std(tempMatF2,[],1)/sqrt(lc) ;
            end      
        end
    end
end

% average tuning curve change across identified unique cell types
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    for tp=1:l_tp ; % for each temporal period
        for dset=2:length(DataBlock(DB).DgPath) ;
            clear tempMat*
            
            lc = length(cell_i{uc}) ; % number of cell of that type
            
            if lc>0 ; % if there are identified cells of this type
                for cells=1:lc ; % for each cell of this type  
                    tempMat(cells,:) = psth_var_cntrstMean_DivCntrl{dset}{tp}(cell_i{uc}(cells),:) ;
                    tempMatF1(cells,:) = psth_varF1_cntrstMean_DivCntrl{dset}{tp}(cell_i{uc}(cells),:) ;
                    tempMatF2(cells,:) = psth_varF2_cntrstMean_DivCntrl{dset}{tp}(cell_i{uc}(cells),:) ;            
                end
                psth_var_cntrstMean_DivCntrl_CellMean{dset}{tp}(uc,:) = mean(tempMat,1) ;
                psth_var_cntrstMean_DivCntrl_CellSem{dset}{tp}(uc,:) = std(tempMat,[],1)/sqrt(lc) ; 

                psth_varF1_cntrstMean_DivCntrl_CellMean{dset}{tp}(uc,:) = mean(tempMatF1,1) ;
                psth_varF1_cntrstMean_DivCntrl_CellSem{dset}{tp}(uc,:) = std(tempMatF1,[],1)/sqrt(lc) ;

                psth_varF2_cntrstMean_DivCntrl_CellMean{dset}{tp}(uc,:) = mean(tempMatF2,1) ;
                psth_varF2_cntrstMean_DivCntrl_CellSem{dset}{tp}(uc,:) = std(tempMatF2,[],1)/sqrt(lc) ;
            end      
        end
    end
end

% rename some datafields so can delete dataRun struct
dataRun_cell_ids = dataRun.cell_ids ; % so can delete dataRun
dataRunMaster_cell_ids = dataRunMaster.cell_ids ;
cell_type_list = get_cell_type_list(dataRunMaster) ;
coord_tform = coordinate_transform(dataRunMaster,'sta');
for c = 1:length(dataRunMaster.cell_ids) ; % for each with a BW mapped ref
    ctr{c} = dataRunMaster.stas.fits{c}.mean ;
    rad{c} = dataRunMaster.stas.fits{c}.sd ;
    angle{c} = dataRunMaster.stas.fits{c}.angle ;
    if ~isempty(cell_list_map{c}) ; % if there is a mapped cell
        slave_c(c) = get_cell_indices(dataRun,cell_list_map{c}) ;
    end
end

for uc = 1:length(UniqueCellTypes) ; % for each cell type
    Master_i{uc} = get_cell_indices(dataRunMaster,UniqueCellTypes{uc}) ;
end

clear dataRun dataRunMaster
save(['/Users/jcafaro/Desktop/Temp/Matfiles/KoGratMultiSetV1',num2str(DB)],'-v7.3')

% FIGURES
if FigureFlag ;
    
% psths of each cell in control, +psem, wash
cells=1 ;
cntrst = 1 ;
while cells <=l_cells ; % for each cell
    figure(1)
    set(gcf,'name',[num2str(cells)])
    clf
    for dset=1:length(DataBlock(DB).DgPath) ;
        for sp=1:l_sp ; % for each spatial period 
            for tp=1:l_tp ; % for each temporal period
                subplot(l_sp,l_tp,l_tp*(sp-1)+tp) ;
                plot(psth{dset}{tp}{cntrst}{sp}(cells,:),'color',Color_list{dset})
                axis tight
                hold on
                title(['']) ;
            end
        end
    end
    
    figure(2)
    clf
    % master cell indicy
    masteri=[] ;
    for clm=1:length(cell_list_map) ;
        if cell_list_map{clm}==dataRun_cell_ids(cells) ;
            masteri=clm ;
        end
    end
    if ~isempty(masteri) && ~isempty(srf{masteri}) ; % if there is a mapped master cell that has a rf
        subplot(2,3,1) % spatial receptive field
        temp_rf = srf{masteri} ;
        norm_rf = norm_image(temp_rf);
        imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
        colormap(brewermap([],'RdBu'))
        caxis([0,1]) 
        set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
        title(num2str(dataRunMaster_cell_ids(masteri)))
        
         subplot(2,3,2) % temporal receptive field
        plot(fliplr(trf_time),trf{masteri}) ;
        title(cell_type_list{masteri})
    end
    for dset=1:length(DataBlock(DB).DgPath) ;
        subplot(2,3,3) % psth variance (spatial frequency)
        plot(log(Spatial_frequency),psth_varF1_cntrstMean{dset}{1}(cells,:),'color',Color_list{dset})
        hold on
        plot(log(Spatial_frequency),psth_varF2_cntrstMean{dset}{1}(cells,:),'color',Color_list{dset},'LineStyle','--')
        set(gca,'Yscale','log')
        xlabel('log 2 spatial f')
        ylabel('psth variance')
        title(num2str(dataRun_cell_ids(cells)))

        subplot(2,3,4) % psth variance (spatial frequency)
        plot(log(Spatial_frequency),psth_var_cntrstMean{dset}{1}(cells,:),'color',Color_list{dset})
        xlabel('log 2 spatial f')
        ylabel('psth variance')
        title(num2str(dataRun_cell_ids(cells)))
        hold on

        subplot(2,3,5) % psth variance (spatial frequency)
        plot(log(Spatial_frequency),psth_varF1_cntrstMean{dset}{1}(cells,:),'color',Color_list{dset})
        xlabel('log 2 spatial f')
        ylabel('psth variance (F1)')
        title(num2str(dataRun_cell_ids(cells)))
        hold on

        subplot(2,3,6) % psth variance (spatial frequency)
        plot(log(Spatial_frequency),psth_varF2_cntrstMean{dset}{1}(cells,:),'color',Color_list{dset})
        xlabel('log 2 spatial f')
        ylabel('psth variance (F2)')
        title(num2str(dataRun_cell_ids(cells)))
        hold on
    end
    
    nxtval = input('next (cell num, 0=back, return=forward)') ;
    if isempty(nxtval) ;
        cells=cells+1 ;
    elseif nxtval == 0 ;
        cells=cells-1 ;
    elseif nxtval>0 ;
        cells=nxtval ;
    end
end

% save some figs
saveFigPath = '/Users/jcafaro/Documents/AnalysisFigures/PsemDb38/' ;
for cells =1:l_cells ; % for each cell
    figure(1)
    clf
    set(gcf,'name',[num2str(cells)])
    for dset=1:length(DataBlock(DB).DgPath) ;
        for sp=1:l_sp ; % for each spatial period 
            for tp=1:l_tp ; % for each temporal period
                subplot(l_sp,l_tp,l_tp*(sp-1)+tp) ;
                plot(psth{dset}{tp}{1}{sp}(cells,:),'color',Color_list{dset})
                axis tight
                hold on
                title(['']) ;
            end
        end
    end
    print(gcf, '-djpeg', [saveFigPath,'Cell',num2str(cells)])
end

% psth [control, +PSEM, wash] organized by unique cell types for a single spatial freq and contrast
sp = 2 ; % spatial frequency default
cntrst = 1 ; % contrast default
tp = 1 ;
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=[1,2] ; 
            plot(PsthTime{tp},psth{dset}{tp}{cntrst}{sp}(cell_i{uc}(cells),:),'color',Color_list{dset})
            hold on
        end
        title(num2str(cell_i{uc}(cells)))
    end
end


% contrast gain functions for each cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=1:3 ; 
            plot(Spatial_frequency,psth_var_cgain{dset}{1}(cell_i{uc}(cells),:),'color',Color_list{dset})
            hold on
        end
        set(gca,'Xscale','log')
        title(num2str(cell_i{uc}(cells)))
    end
end
    
% spatial tuning functions averaged across contrasts for each cell type
tp = 1 ;
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=1:length(DataBlock(DB).DgPath) ;
            plot(Spatial_frequency,psth_var_cntrstMean{dset}{tp}(cell_i{uc}(cells),:),'color',Color_list{dset})
            %plot(Spatial_frequency,psth_varF1_cntrstMean{dset}{tp}(cell_i{uc}(cells),:),'color',Color_list{dset},'LineStyle','--')
            hold on 
            %plot(Spatial_frequency,psth_varF2_cntrstMean{dset}{tp}(cell_i{uc}(cells),:),'color',Color_list{dset})
            set(gca,'yscale','log')
        end
        set(gca,'Xscale','log')
        title(num2str(cell_i{uc}(cells)))
    end
end

% response as function of sp averaged across contrast - averaged within cell type
NumFigRows = 4 ;
for tp=1:l_tp ;
    figure
    set(gcf,'name',num2str(tp))
    for uc = 1:length(UniqueCellTypes) ; % for each cell type
        subplot(NumFigRows, ceil(length(UniqueCellTypes)/NumFigRows),uc) ;
        for dset=1:length(DataBlock(DB).DgPath) ;       
            errorbar(Spatial_frequency,psth_var_cntrstMean_CellMean{dset}{tp}(uc,:),...
                psth_var_cntrstMean_CellSem{dset}{tp}(uc,:),'color',Color_list{dset}) ;
            hold on
        end
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        axis tight
        title(UniqueCellTypes{uc})
    end
end

% response change as function of sp averaged across contrast - averaged within cell type
NumFigRows = 4 ;
figure
for tp=1:l_tp ;
    for uc = 1:length(UniqueCellTypes) ; % for each cell type
        subplot(NumFigRows, ceil(length(UniqueCellTypes)/NumFigRows),uc) ;
        for dset=2:length(DataBlock(DB).DgPath) ;       
            errorbar(Spatial_frequency,psth_var_cntrstMean_DivCntrl_CellMean{dset}{tp}(uc,:),...
                psth_var_cntrstMean_DivCntrl_CellSem{dset}{tp}(uc,:),'color',Color_list{tp}) ;
            hold on
        end
        set(gca,'xscale','log')
        plot(Spatial_frequency,ones(1,length(Spatial_frequency)),'k:')
        axis tight
        title(UniqueCellTypes{uc})
    end
end

% response change as function of sp avearged across contrast - grouped within cell type
NumFigRows = 4 ;
for tp=1:l_tp ;
    figure
    set(gcf,'name',num2str(tp))
    for uc = 1:length(UniqueCellTypes) ; % for each cell type
        subplot(NumFigRows, ceil(length(UniqueCellTypes)/NumFigRows),uc) ;
        plot(Spatial_frequency,ones(1,length(Spatial_frequency)),'k')
        hold on 
        set(gca,'xscale','log')
        lc = length(cell_i{uc}) ; % number of cell of that type
        for cells=1:lc ; % for each cell of this type
            for dset=2:length(DataBlock(DB).DgPath) ;
                plot(Spatial_frequency,psth_var_cntrstMean_DivCntrl{dset}{tp}(cell_i{uc}(cells),:),'color',Color_list{dset})
            end
           
            axis tight
            set(gca,'ylim',[0,2])
            
            title(UniqueCellTypes{uc})
        end
    end
end

% response change vs. mean response - is change just drop in weak responses?
figure
for tp=1:l_tp ;
    
    set(gcf,'name',num2str(tp))
    for uc = 1:length(UniqueCellTypes) ; % for each cell type
        subplot(NumFigRows, ceil(length(UniqueCellTypes)/NumFigRows),uc) ;
        lc = length(cell_i{uc}) ; % number of cell of that type
        
        for dset=2:length(DataBlock(DB).DgPath) ;
            plot(psth_var_cntrstMean_CellMean{dset}{tp}(uc,:),...
                psth_var_cntrstMean_DivCntrl_CellMean{dset}{tp}(uc,:),'*','color',Color_list{dset}) ;
            hold on
        end
        set(gca,'xscale','log')
        axis tight
        
        title(UniqueCellTypes{uc})
    end
end
 
%% map location of cells on epiflourescence image
% ei centers
numElectrodeLayers = 2 ;
for cells = 1:length(dataRun.spikes)  ; % for each cell
    EiCtr(cells,:) = get_ei_com(dataRun, dataRun.cell_ids(cells), numElectrodeLayers) ;
end


Tform = map_Ei_to_camera(DataBlock(DB).IrImagePath, DataBlock(DB).BwRepeat{1}) ; % get transform to camera

%figure
CheckOn = get_cell_indices(dataRun,[947,1877,1878,2328,2912,3124,4816,6633,6635,6961,7156]) ;

im_array = imread(DataBlock(DB).EpiImagePath);

figure
imshow(im_array);
hold on
cnd = 2 ;
for cells=1:length(dataRun.spikes) ;
    if ismember(cells,CheckOn)
        col = 'r' ;
    else
        col = 'c' ;
    end
    
    ctrImage = tformfwd(Tform,EiCtr(cells,1),EiCtr(cells,2)) ;
    plot(ctrImage(1),ctrImage(2),'+','Color',col)
    %plot(ctrImage(1),ctrImage(2),'o','MarkerSize',100*abs(psth_mean_dprime{cells}(cnd))+.001,'Color',col)
end



%%

% For Igor export
ForIgor = struct ;

VecName = ['MseHistX','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,MseRelHistX) ; % mse bins

VecName = ['MseDrugOverWashCdf','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,MseRelDrug_DivMseRelWash_Hist_cumsum) ; % mse drug:wash ratio cumulative pdf

VecName = ['MseLfmDrugOverWashCdf','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,MseRelDrug_DivMseRelWashLfMean_Hist_cumsum) ; % mse drug:wash ratio cumulative pdf

VecName = ['MseLfmWashCdf','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,MseRelWashLfMean_Hist_cumsum) ; % mse drug:wash ratio cumulative pdf

VecName = ['MseLfmDrugCdf','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,MseRelDrugLfMean_Hist_cumsum) ; % mse drug:wash ratio cumulative pdf

exportStructToHDF5(ForIgor,['/Users/jcafaro/Desktop/Temp/Matfiles/KoSpatialV3Db',num2str(DB),'.h5'],'/');

else
ForIgor.Empty = [] ;
end



