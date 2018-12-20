function ForIgor = KoSpatialTuningAnalyzerConcatPlusMapping(DataBlock, DB, Params)

% this function is adapted from 'KoSpatialTuningAnalyzerConcat' to include
% cell mapping from a binary white noise stimulus

% JC 12/17/2015 

% parameters
Color_list = {'k','r','b','g','c','v'} ; % order of colors for each 
TrialTrigInterval = 10 ; % sec
InterTrialTime = 2 ; % sec
StimTransTime = .1 ; % (sec) skip cycles before this time after starting new stimulus
MicronPerPix = 4 ; % (um/pix)
EsemFlag = false; % calc psth esemble

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
for dset=1:3 ;
    % load data 
    dataRun = load_data(DataBlock(DB).DgPath{dset}) ;
    dataRun = load_neurons(dataRun) ;

    % load stim
    dataRun.names.stimulus_path = [DataBlock(DB).DgPath{dset}(1:end-15),'stimuli/s',DataBlock(DB).DgPath{dset}(end-1:end)] ;
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

Spatial_frequency = 1./(SpatialPeriod*MicronPerPix) ;

% load concatinated data
dataRun = load_data(DataBlock(DB).DgConcat) ;
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
    end
end
clear dataRunMaster

% prep cells and get trigger times
for dset=1:3 ;
    for sp=1:l_sp ; % for each spatial period
        for cntrst=1:l_cntrst ; % for each contrast
            for tp=1:l_tp ; % for each temporal period
                stim_group_triggers{dset}{sp,tp, cntrst}= dataRun.triggers(stim_group_triggers_i{dset}{sp,tp,cntrst}) ; % trigger times
                psth{dset}{tp}{cntrst}{sp} = nans(l_cells,ceil(min(diff(stim_group_triggers{dset}{sp,tp,cntrst}'))*dataRun.sampling_rate)) ;
            end
        end
    end

    for tp=1:l_tp ; % for each temporal period
        for cntrst=1:l_cntrst ; % for each contrast
            for cells = 1:l_cells ; % for each cell
                for sp=1:l_sp ; % for each spatial period
                    [psthTemp,PsthTime] = get_smooth_psth(dataRun.spikes{cells},stim_group_triggers{dset}{sp,tp,cntrst}','stop',min(diff(stim_group_triggers{dset}{sp,tp,cntrst}'))) ;
                    psth{dset}{tp}{cntrst}{sp}(cells,:) = psthTemp ;
                    
                    if EsemFlag ; % if you want psth error bars
                        psthEsem = nans(length(stim_group_triggers{dset}{sp,tp,cntrst}'),ceil(min(diff(stim_group_triggers{dset}{sp,tp,cntrst}'))/.1)) ; % prep for speed
                        for trl=1:length(stim_group_triggers{dset}{sp,tp,cntrst}') ; % for each trial 
                            psthEsem(trl,:) = get_psth(dataRun.spikes{cells},stim_group_triggers{dset}{sp,tp,cntrst}(trl)','stop',min(diff(stim_group_triggers{dset}{sp,tp,cntrst}'))) ; % notice that this is a histogram not a sliding window
                            psthEsem_var(trl) = var(psthEsem(trl,:)) ;
                        end
                        psthEsem_var_mean{dset}{tp}{cntrst}(cells,sp) = mean(psthEsem_var) ; % average variance of this cell and condition
                        psthEsem_var_sem{dset}{tp}{cntrst}(cells,sp) = std(psthEsem_var)/sqrt(trl) ; % sem
                        clear psthEsem psthEsem_var
                    end    
                        
                    psth_mean{dset}{tp}{cntrst}(cells,sp) = mean(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                    
                    mn =  min(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                    psth_min{dset}{tp}{cntrst}(cells,sp) = mn ;
                    
                    [mx,mi] =  max(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                    psth_peak{dset}{tp}{cntrst}(cells,sp) = mx ;
                    psth_peak_time{dset}{tp}{cntrst}(cells,sp) = PsthTime(mi) ;
                    
                    psth_peakDivMean{dset}{tp}{cntrst}(cells,sp) = mx/psth_mean{dset}{tp}{cntrst}(cells,sp) ;
                      
                    [powerspec_xvalues, mean_powerspec, mean_phase] = PowerPhaseFinder(psth{dset}{tp}{cntrst}{sp}(cells,:),dataRun.sampling_rate) ;
                    if sum(mean_powerspec) ~=0 ; % if there is any power 
                        psth_var_f1{dset}{tp}{cntrst}(cells,sp) = mean_powerspec(2)/sum(mean_powerspec) ; % f1 power
                        psth_var_f2{dset}{tp}{cntrst}(cells,sp) = mean_powerspec(3)/sum(mean_powerspec) ; % f2 power
                    
                        psth_phase_f1{dset}{tp}{cntrst}(cells,sp) = mean_phase(2) ; % (deg) phase of f1 
                        psth_phase_f2{dset}{tp}{cntrst}(cells,sp) = mean_phase(3) ; % (deg) phase of f2
                    else
                        psth_var_f1{dset}{tp}{cntrst}(cells,sp) = 0 ; % f1 power
                        psth_var_f2{dset}{tp}{cntrst}(cells,sp) = 0 ; % f2 power
                        
                        psth_phase_f1{dset}{tp}{cntrst}(cells,sp) = nan ; % (deg) phase of f1 
                        psth_phase_f2{dset}{tp}{cntrst}(cells,sp) = nan ;
                    end
                    psth_var{dset}{tp}{cntrst}(cells,sp) = nanvar(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                end
                
                [mx,mi] = max(psth_var{dset}{tp}{cntrst}(cells,:)) ;
                psth_var_norm{dset}{tp}{cntrst}(cells,:) = psth_var{dset}{tp}{cntrst}(cells,:)/mx ;
                psth_var_peak{tp}{cntrst}(cells,dset) = SpatialPeriod(mi) ;

                [mx,mi] = max(psth_var_f1{dset}{tp}{cntrst}(cells,:)) ;
                psth_var_f1_norm{dset}{tp}{cntrst}(cells,:) = psth_var_f1{dset}{tp}{cntrst}(cells,:)/mx ;
                psth_var_f1_peak{tp}{cntrst}(cells,dset) = SpatialPeriod(mi) ;

                [mx,mi] = max(psth_var_f2{dset}{tp}{cntrst}(cells,:)) ;
                psth_var_f2_norm{dset}{tp}{cntrst}(cells,:) = psth_var_f2{dset}{tp}{cntrst}(cells,:)/mx ;
                psth_var_f2_peak{tp}{cntrst}(cells,dset) = SpatialPeriod(mi) ;

                [mx,mi] = max(psth_peak{dset}{tp}{cntrst}(cells,:)) ;
                psth_peak_norm{dset}{tp}{cntrst}(cells,:) = psth_peak{dset}{tp}{cntrst}(cells,:)/mx ;
            end

            psth_var_mean{dset}{tp}{cntrst} = nanmean(psth_var{dset}{tp}{cntrst},1) ; % nans can happen when no spikes are detected
            psth_var_norm_mean{dset}{tp}{cntrst} = nanmean(psth_var_norm{dset}{tp}{cntrst},1) ;

            psth_var_f1_mean{dset}{tp}{cntrst} = nanmean(psth_var_f1{dset}{tp}{cntrst},1) ;
            psth_var_f1_norm_mean{dset}{tp}{cntrst} = nanmean(psth_var_f1_norm{dset}{tp}{cntrst},1) ; 

            psth_var_f2_mean{dset}{tp}{cntrst} = nanmean(psth_var_f2{dset}{tp}{cntrst},1) ;
            psth_var_f2_norm_mean{dset}{tp}{cntrst} = nanmean(psth_var_f2_norm{dset}{tp}{cntrst},1) ;
            
            psth_peak_mean{dset}{tp}{cntrst} = nanmean(psth_peak{dset}{tp}{cntrst},1) ; % nans can happen when no spikes are detected
            psth_peak_norm_mean{dset}{tp}{cntrst} = nanmean(psth_peak_norm{dset}{tp}{cntrst},1) ;
            
            psth_mean_mean{dset}{tp}{cntrst} = nanmean(psth_mean{dset}{tp}{cntrst},1) ;
            
            psth_peakDivMean_mean{dset}{tp}{cntrst} = nanmean(psth_peakDivMean{dset}{tp}{cntrst},1) ; 
        end
    end
end

% compare control, drug and wash conditions

% find psth var drug effect index
for cntrst=1:l_cntrst ; % for each contrast
    for cells = 1:l_cells ; % for each cell
%         TempNm = psth_var{2}{1}{cntrst}(cells,:) - (psth_var{1}{1}{cntrst}(cells,:)+psth_var{3}{1}{cntrst}(cells,:))/2 ; % drug-(cnt+wash)/2
%         TempDnm = abs(psth_var{1}{1}{cntrst}(cells,:)-psth_var{3}{1}{cntrst}(cells,:)) ; % abs(cntrl-wash)
%         psth_var_deltaDrug{cntrst}(cells,:)  = TempNm./(TempDnm+10) ;

        TempNm = psth_var{2}{1}{cntrst}(cells,:) - psth_var{1}{1}{cntrst}(cells,:) ; % drug-cnt
        TempDnm = 1-abs((psth_var{1}{1}{cntrst}(cells,:)-psth_var{3}{1}{cntrst}(cells,:))./(psth_var{1}{1}{cntrst}(cells,:)+psth_var{3}{1}{cntrst}(cells,:))) ; % 1-abs(cntrl-wash/cnt+wash)
        psth_var_deltaDrug{cntrst}(cells,:)  = TempNm.*TempDnm ; % units are same as variance scaled by certainty (+ drug is bigger, - drug is smaller)
    end
end

% find psth drug effect index
for cntrst=1:l_cntrst ; % for each contrast
    for sp=1:l_sp ; % for each spatial period
        for cells = 1:l_cells ; % for each cell
            % norm vectors (v-mean(v))/norm(v-mean(v))
            psth1n = (psth{1}{1}{cntrst}{sp}(cells,:)-mean(psth{1}{1}{cntrst}{sp}(cells,:)))/norm(psth{1}{1}{cntrst}{sp}(cells,:)-mean(psth{1}{1}{cntrst}{sp}(cells,:))) ; 
            psth2n = (psth{2}{1}{cntrst}{sp}(cells,:)-mean(psth{2}{1}{cntrst}{sp}(cells,:)))/norm(psth{2}{1}{cntrst}{sp}(cells,:)-mean(psth{2}{1}{cntrst}{sp}(cells,:))) ; 
            psth3n = (psth{3}{1}{cntrst}{sp}(cells,:)-mean(psth{3}{1}{cntrst}{sp}(cells,:)))/norm(psth{3}{1}{cntrst}{sp}(cells,:)-mean(psth{3}{1}{cntrst}{sp}(cells,:))) ; 
            
            TempNm = (psth1n*psth3n') ;% (cntrl*wash)
            if TempNm<0 ;
                TempNm = 0 ; % rectify (neg dot products are as bad as unrelated)
            end
            
            TempDnm = (psth1n*psth2n') ; %(cntrl*drug)
            if TempDnm<0 ;
                TempDnm = 0 ; % rectify (neg dot products are as bad as unrelated)
            end
            
            MseWash = sum((psth{1}{1}{cntrst}{sp}(cells,:)-psth{3}{1}{cntrst}{sp}(cells,:)).^2) ; % mse of cntrl-wash
            Var1 = psth{1}{1}{cntrst}{sp}(cells,:)*psth{1}{1}{cntrst}{sp}(cells,:)' ;
            Var3 = psth{3}{1}{cntrst}{sp}(cells,:)*psth{3}{1}{cntrst}{sp}(cells,:)' ;
            MseWashIndex = (Var1+Var3-MseWash)/(Var1+Var3) ; % index (0:1 1=same exact wave)
            
            psth_deltaDrug{cntrst}{sp}(cells)  = TempNm*(TempNm-TempDnm) ; % (index -1:1, 1=big drug effect)
            psth_deltaWash{cntrst}{sp}(cells)  = TempNm*MseWashIndex ; % (index 0:1, punished for both mse and norm dot differences)
        
            if isnan(psth_deltaWash{cntrst}{sp}(cells)) ; % nan if no spikes detected
                psth_deltaWash{cntrst}{sp}(cells) = 0 ;
            end
                
        end
    end
end
    
% phase, time, and var change
for dset=1:3 ;
    for tp=1:l_tp ; % for each temporal period
        for sp=1:l_sp ; % for each spatial period
            for cntrst=1:l_cntrst ; % for each contrast
                for cells = 1:l_cells ; % for each cell
                    psth_phase_f1_deltaCntrl{dset}{tp}{cntrst}(cells,sp) = psth_phase_f1{dset}{tp}{cntrst}(cells,sp)- psth_phase_f1{1}{tp}{cntrst}(cells,sp) ; % condition-control
                    psth_phase_f2_deltaCntrl{dset}{tp}{cntrst}(cells,sp) = psth_phase_f2{dset}{tp}{cntrst}(cells,sp)- psth_phase_f2{1}{tp}{cntrst}(cells,sp) ; % condition-control              
                
                    psth_peak_time_deltaCntrl{dset}{tp}{cntrst}(cells,sp) = psth_peak_time{dset}{tp}{cntrst}(cells,sp)- psth_peak_time{1}{tp}{cntrst}(cells,sp) ; % condition-control
                
                    psth_var_divCntrl{dset}{tp}{cntrst}(cells,sp) = psth_var{dset}{tp}{cntrst}(cells,sp)/psth_var{1}{tp}{cntrst}(cells,sp) ; % condition/control
                    psth_var_minCntrl{dset}{tp}{cntrst}(cells,sp) = psth_var{dset}{tp}{cntrst}(cells,sp)-psth_var{1}{tp}{cntrst}(cells,sp) ; % condition-control
                    
                    tempcc = xcov(psth{1}{tp}{cntrst}{sp}(cells,:),psth{dset}{tp}{cntrst}{sp}(cells,:),length(psthTemp)/2) ; % keeping the max lag within half a cycle to prevent correlation with earlier phase
                    timecc = ([1:length(tempcc)]-ceil(length(tempcc)/2))/dataRun.sampling_rate ;
                    [mx,mi] = max(tempcc) ;
                    psth_cct_deltaCntrl{dset}{tp}{cntrst}(cells,sp) = timecc(mi) ; % negative if dset psth is trailing
                end
            end
        end
    end 
end

dataRunCell_ids = dataRun.cell_ids ; % so can delete dataRun

clear dataRun

% FIGURES

lineWidthRange = 3 ;
% psths of each cell in control, +psem, wash
cells=1 ;
while cells <=l_cells ; % for each cell
    figure(1)
    clf
    for dset=1:3 ;
        for sp=1:l_sp ; % for each spatial period 
            for cntrst=1:l_cntrst ; % for each contrast
                subplot(l_sp,l_cntrst,l_cntrst*(sp-1)+cntrst) ;
                %plot(PsthTime,psth{dset}{1}{cntrst}{sp}(cells,:),'color',Color_list{dset},'linewidth',lineWidthRange*psth_deltaWash{cntrst}{sp}(cells)+.1)
                plot(PsthTime,psth{dset}{1}{cntrst}{sp}(cells,:),'color',Color_list{dset})
                axis tight
                hold on
            end
        end
    end
    title(num2str(cells))
    
    figure(2)
    clf
    % master cell indicy
    masteri=[] ;
    for clm=1:length(cell_list_map) ;
        if cell_list_map{clm}==dataRunCell_ids(cells) ;
            masteri=clm ;
        end
    end
    if ~isempty(masteri) ; % if there is a mapped master cell
        subplot(6,l_cntrst,1) % spatial receptive field
        temp_rf = srf{masteri} ;
        norm_rf = norm_image(temp_rf);
        imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
        colormap(brewermap([],'RdBu'))
        caxis([0,1]) 
        set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])

        subplot(6,l_cntrst,2) % temporal receptive field
        plot(fliplr(trf_time),trf{masteri}) ;
    end
    for dset=1:3 ;
        for cntrst=1:l_cntrst ; % for each contrast

            subplot(6,l_cntrst,l_cntrst+cntrst) % psth variance (spatial frequency)
            plot(log(Spatial_frequency),psth_var{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            %errorbar(log(Spatial_frequency),psthEsem_var_mean{dset}{1}{cntrst}(cells,:),psthEsem_var_sem{dset}{1}{cntrst}(cells,:),psthEsem_var_sem{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth variance')
            hold on

            subplot(6,l_cntrst,l_cntrst*2+cntrst) % psth mean (spatial frequency)
            plot(log(Spatial_frequency),psth_mean{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth mean')
            hold on

            subplot(6,l_cntrst,l_cntrst*3+cntrst) % psth peak/mean  (spatial frquency)
            plot(log(Spatial_frequency),psth_peakDivMean{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth peak/mean')
            hold on
            
            subplot(6,l_cntrst,l_cntrst*4+cntrst) % psth peak/mean  (spatial frquency)
            plot(log(Spatial_frequency),psth_cct_deltaCntrl{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth cc peak time')
            hold on
            
            subplot(6,l_cntrst,l_cntrst*5+cntrst) % psth variance (spatial frequency)
            plot(log(Spatial_frequency),psth_var_norm{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('norm psth variance')
            set(gca,'yscale','log')
            hold on
            
        end
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

figure % spatial/temporal receptive fields + variance, peak and peak_time (spatial freq) 
sp = [2,10,11] ; % spatial frequency default
cntrst = [2,2,2] ; % contrast default
cells=1 ;
while cells <=l_cells ; % for each cell
    clf
    % master cell indicy
    masteri=[] ;
    for clm=1:length(cell_list_map) ;
        if cell_list_map{clm}==dataRunCell_ids(cells) ;
            masteri=clm ;
        end
    end
    if ~isempty(masteri) ; % if there is a mapped master cell
        subplot(3,3,1) % spatial receptive field
        temp_rf = srf{masteri} ;
        norm_rf = norm_image(temp_rf);
        imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
        colormap(brewermap([],'RdBu'))
        caxis([0,1]) 
        set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])

        subplot(3,3,2) % temporal receptive field
        plot(fliplr(trf_time),trf{masteri}) ;
    end
    for dset=1:3 ;

        subplot(3,3,4) % psth variance (spatial frequency)
        plot(log(Spatial_frequency),psth_var{dset}{1}{cntrst(1)}(cells,:),'color',Color_list{dset})
        xlabel('log 2 spatial f')
        ylabel('psth variance')
        hold on
        
        subplot(3,3,5) % psth peak (spatial frequency)
        plot(log(Spatial_frequency),psth_peak{dset}{1}{cntrst(1)}(cells,:),'color',Color_list{dset})
        xlabel('log 2 spatial f')
        ylabel('psth peak')
        hold on
        
        subplot(3,3,6) % psth peak time (spatial frquency)
        plot(log(Spatial_frequency),psth_peak_time{dset}{1}{cntrst(1)}(cells,:),'color',Color_list{dset})
        xlabel('log 2 spatial f')
        ylabel('psth peak time')
        hold on
        
        subplot(3,3,7) % psth
        plot(PsthTime,psth{dset}{1}{cntrst(1)}{sp(1)}(cells,:),'color',Color_list{dset})
        xlabel('time (s)')
        ylabel('spike rate (hz)')
        hold on
        
        subplot(3,3,8) % psth
        plot(PsthTime,psth{dset}{1}{cntrst(2)}{sp(2)}(cells,:),'color',Color_list{dset})
        xlabel('time (s)')
        ylabel('spike rate (hz)')
        hold on
        
        subplot(3,3,9) % psth
        plot(PsthTime,psth{dset}{1}{cntrst(3)}{sp(3)}(cells,:),'color',Color_list{dset})
        xlabel('time (s)')
        ylabel('spike rate (hz)')
        hold on
    end
        
    title(num2str(cells))
    nxtval = input('next (cell num, 0=back, return=forward)') ;
    if isempty(nxtval) ;
        cells=cells+1 ;
    elseif nxtval == 0 ;
        cells=cells-1 ;
    elseif nxtval>0 ;
        cells=nxtval ;
    end
end


% psth [control, +PSEM, wash] organized by unique cell types for a single spatial freq and contrast
sp = 7 ; % spatial frequency default
cntrst = 1 ; % contrast default
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=1:3 ; 
            plot(PsthTime,psth{dset}{1}{cntrst}{sp}(cell_i{uc}(cells),:),'color',Color_list{dset})
            hold on
        end
        title(num2str(cell_i{uc}(cells)))
    end
end
    
figure % spatial tuning plots normalized with peak comparisons
for cntrst=1:l_cntrst ; % for each contrast
    for dset=1:3 ;
        subplot(2,l_cntrst,cntrst)
        plot(log(Spatial_frequency), psth_var_norm_mean{dset}{1}{cntrst},Color_list{dset},'linewidth',4)
        %set(gca,'xscale','log')
        hold on
    end
    subplot(2,l_cntrst,l_cntrst+cntrst)
    plot([1:3],psth_var_peak{1}{cntrst})
    hold on
    plot([1:3],mean(psth_var_peak{1}{cntrst}),'k','linewidth',4)
end

figure % spatial tuning plots of f1 normalized with peak comparisons
for cntrst=1:l_cntrst ; % for each contrast
    for dset=1:3 ;
        subplot(2,l_cntrst,cntrst)
        plot(SpatialPeriod, psth_var_f1_norm_mean{dset}{1}{cntrst},Color_list{dset},'linewidth',4) 
        hold on
    end
    subplot(2,l_cntrst,l_cntrst+cntrst)
    plot([1:3],psth_var_f1_peak{1}{cntrst})
    hold on
    plot([1:3],mean(psth_var_f1_peak{1}{cntrst}),'k','linewidth',4)
end

figure % spatial tuning plots of f2 normalized with peak comparisons
for cntrst=1:l_cntrst ; % for each contrast
    for dset=1:3 ;
        subplot(2,l_cntrst,cntrst)
        plot(SpatialPeriod, psth_var_f2_norm_mean{dset}{1}{cntrst},Color_list{dset},'linewidth',4) 
        hold on
    end
    subplot(2,l_cntrst,l_cntrst+cntrst)
    plot([1:3],psth_var_f2_peak{1}{cntrst})
    hold on
    plot([1:3],mean(psth_var_f2_peak{1}{cntrst}),'k','linewidth',4)
end

figure % spatial tuning plots non normalized with drug comparisons
for cntrst=1:l_cntrst ; % for each contrast
    for dset=1:3 ;
        subplot(2,l_cntrst,cntrst)
        plot(SpatialPeriod, psth_var_mean{dset}{1}{cntrst},Color_list{dset},'linewidth',4) 
        hold on
        
    end
    subplot(2,l_cntrst,l_cntrst+cntrst)
    plot(SpatialPeriod, psth_var_deltaDrug{cntrst}) 
    hold on
    errorbar(SpatialPeriod, mean(psth_var_deltaDrug{cntrst}),std(psth_var_deltaDrug{cntrst}),std(psth_var_deltaDrug{cntrst}),'k','linewidth',4)  
end

figure % spatial tuning peak plots normalized with peak comparisons
for cntrst=1:l_cntrst ; % for each contrast
    for dset=1:3 ;
        subplot(2,l_cntrst,cntrst)
        plot(SpatialPeriod, psth_peak_norm_mean{dset}{1}{cntrst},Color_list{dset}) 
        hold on
    end
    for dset=1:3 ;
        subplot(2,l_cntrst,cntrst+l_cntrst)
        plot(SpatialPeriod, psth_peak_mean{dset}{1}{cntrst},Color_list{dset}) 
        hold on
    end
end

figure 



% phase changes
figure
for cntrst=1:l_cntrst ;
    for sp=1:l_sp ;
        subplot(l_sp,l_cntrst,l_cntrst*(sp-1)+cntrst)
        plot([1:3],nanmean(psth_phase_f1_drugDelta{1}{cntrst}{sp},1))
        hold on
        subplot(l_sp,l_cntrst,l_cntrst*(sp-1)+cntrst)
        plot([1:3],nanmean(psth_phase_f2_drugDelta{1}{cntrst}{sp},1))
        hold on
    end
end

        
%set(gca,'XScale','log')    



