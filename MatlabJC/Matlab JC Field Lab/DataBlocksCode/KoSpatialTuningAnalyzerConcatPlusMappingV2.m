function ForIgor = KoSpatialTuningAnalyzerConcatPlusMappingV2(DataBlock, DB, Params)

% this function is adapted from 'KoSpatialTuningAnalyzerConcatPlusMapping'.
%  Changed parameter analysis (eliminated lots that didn't seem useful)

% JC 1/21/2016 

% parameters
Color_list = {'k','r','b','g','c','v'} ; % order of colors for each 
TrialTrigInterval = 10 ; % sec
InterTrialTime = 2 ; % sec
StimTransTime = .1 ; % (sec) skip cycles before this time after starting new stimulus
MicronPerPix = 4 ; % (um/pix)
EsemFlag = false; % calc psth esemble
MseRelThreshold = 0.5 ; % (unitless) below this point response is "stable"
MseRelDrugDivWashThreshold = 1.5 ; % (unitless) above this point response has big drug effect 
LfSpThresh = 150 ; % (pix) spatial periods equal and above are termed "low frequency"

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
%clear dataRunMaster

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
                    
                    psth_var{dset}{tp}{cntrst}(cells,sp) = var(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                    psth_mean{dset}{tp}{cntrst}(cells,sp) = mean(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                    psth_range{dset}{tp}{cntrst}(cells,sp) = range(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                    
                    psth_peak{dset}{tp}{cntrst}(cells,sp) = max(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                    psth_peakDivMean{dset}{tp}{cntrst}(cells,sp) = psth_peak{dset}{tp}{cntrst}(cells,sp)/psth_mean{dset}{tp}{cntrst}(cells,sp) ;
                    psth_duty{dset}{tp}{cntrst}(cells,sp) = sum(psth{dset}{tp}{cntrst}{sp}(cells,:)>(psth_peak{dset}{tp}{cntrst}(cells,sp)/2))/length(psthTemp) ; % fraction of time above 50%
                    
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
                end

                % normalize psth_var for each dataset, contrast, and cell 
                [mx,mi] = max(psth_var{dset}{tp}{cntrst}(cells,:)) ;
                psth_var_norm{dset}{tp}{cntrst}(cells,:) = psth_var{dset}{tp}{cntrst}(cells,:)/mx ; 
                psth_var_peak{tp}{cntrst}(cells,dset) = SpatialPeriod(mi) ; % peak of spatial frequency
            end

            % average across cells
            psth_var_mean{dset}{tp}{cntrst} = nanmean(psth_var{dset}{tp}{cntrst},1) ; % nans can happen when no spikes are detected
            psth_var_norm_mean{dset}{tp}{cntrst} = nanmean(psth_var_norm{dset}{tp}{cntrst},1) ;
 
            psth_mean_mean{dset}{tp}{cntrst} = nanmean(psth_mean{dset}{tp}{cntrst},1) ;            
        end
    end
end

% compare control, drug and wash conditions

% find psth drug mse and mse/mean
for cells = 1:l_cells ; % for each cell 
    for cntrst=1:l_cntrst ; % for each contrast
        for sp=1:l_sp ; % for each spatial period   
            MseWash{cells}(cntrst,sp) = mean((psth{1}{1}{cntrst}{sp}(cells,:)-psth{3}{1}{cntrst}{sp}(cells,:)).^2) ; % mse of cntrl-wash
            MseDrug{cells}(cntrst,sp) = mean((psth{1}{1}{cntrst}{sp}(cells,:)-psth{2}{1}{cntrst}{sp}(cells,:)).^2) ; % mse of cntrl-drug 
            
            MseRelWash{cells}(cntrst,sp) = sqrt(MseWash{cells}(cntrst,sp))/mean((psth{1}{1}{cntrst}{sp}(cells,:)+psth{3}{1}{cntrst}{sp}(cells,:))/2) ; % mse(Wash)/mean(cntrl+wash) 
            MseRelDrug{cells}(cntrst,sp) = sqrt(MseDrug{cells}(cntrst,sp))/mean((psth{1}{1}{cntrst}{sp}(cells,:)+psth{2}{1}{cntrst}{sp}(cells,:))/2) ; % mse(Wash)/mean(cntrl)
    
            MseRelDrug_DivMseRelWash{cells}(cntrst,sp) = MseRelDrug{cells}(cntrst,sp)/MseRelWash{cells}(cntrst,sp) ;
            
            % put responses into groups accordin to stability of wash and relative drug effect
            if MseRelDrug{cells}(cntrst,sp)>=MseRelWash{cells}(cntrst,sp) ; % if drug change >= wash change
                if MseRelWash{cells}(cntrst,sp)<=MseRelThreshold ; % if wash change <= threshold
                    MseGroup{cells}(cntrst,sp) = 1 ; % group 1 (wash is stable, drug has bigger effect)
                elseif MseRelWash{cells}(cntrst,sp)>MseRelThreshold ; % if wash change > threshold
                    MseGroup{cells}(cntrst,sp) = 2 ; % group 1 (wash is not stable, drug has bigger effect)
                end
            elseif MseRelDrug{cells}(cntrst,sp)<MseRelWash{cells}(cntrst,sp) ; % if drug change < wash change
                if MseRelWash{cells}(cntrst,sp)<=MseRelThreshold ; % if wash change <= threshold
                    MseGroup{cells}(cntrst,sp) = 3 ; % group 3 (wash is stable, drug has no effect)
                elseif MseRelWash{cells}(cntrst,sp)>MseRelThreshold ; % if wash change > threshold
                    MseGroup{cells}(cntrst,sp) = 4 ; % group 1 (wash is not stable, drug has small effect)
                end
            end
            if mean(psth{1}{1}{cntrst}{sp}(cells,:))==0 ; % if the control mean is zero
                MseGroup{cells}(cntrst,sp) = 0 ;
            end
        end
    end
    
    % average across contrasts and all low frequency stimuli
    lfSet = find(SpatialPeriod>=LfSpThresh) ;
    hfSet = find(SpatialPeriod<LfSpThresh) ;
    MseRelWashLfMean(cells) = nanmean(mean(MseRelWash{cells}(:,lfSet))) ; 
    MseRelDrugLfMean(cells) = nanmean(mean(MseRelDrug{cells}(:,lfSet))) ; 

    MseRelDrug_DivMseRelWashLfMean(cells) = MseRelDrugLfMean(cells)/MseRelWashLfMean(cells) ;
    
    % MseGroup analyses
    num_r = length(MseGroup{cells}(:)) ;
    tempGroup = MseGroup{cells}(:) ;

    MseGroup_sameFract(cells) = sum(tempGroup==mode(tempGroup))/num_r ; % fraction of responses...in the same group
    MseGroup_drugChangeFract(cells) = max(sum(tempGroup==1 | tempGroup==2),sum(tempGroup==3 | tempGroup==4))/num_r ; % in either 1/2 or 3/4
    MseGroup_stableFract(cells) = max(sum(tempGroup==1 | tempGroup==3),sum(tempGroup==2 | tempGroup==4))/num_r ; % in either 1/3 or 2/4
    
end

MseRelHistX = [0:.1:10] ;
tempMat = cell2mat(MseRelDrug_DivMseRelWash) ;
MseRelDrug_DivMseRelWash_Hist = hist(tempMat(:),MseRelHistX) ;
MseRelDrug_DivMseRelWash_FractAboveThresh = sum(tempMat(:)>MseRelDrugDivWashThreshold)/length(tempMat(:)) ;

MseRelDrug_DivMseRelWashLfMean_Hist = hist(MseRelDrug_DivMseRelWashLfMean,MseRelHistX) ;
MseRelDrug_DivMseRelWashLfMean_FractAboveThresh = sum(MseRelDrug_DivMseRelWashLfMean>MseRelDrugDivWashThreshold)/length(MseRelDrug_DivMseRelWashLfMean) ;
   
numResponsesAbove = sum(tempMat(:)>MseRelDrugDivWashThreshold) ;
MseRelDrug_DivMseRelWash_FractAboveThreshHist = zeros(l_cntrst, l_sp) ;
for cells=1:l_cells ;
    MseRelDrug_DivMseRelWash_FractAboveThreshHist = MseRelDrug_DivMseRelWash_FractAboveThreshHist + (MseRelDrug_DivMseRelWash{cells}>MseRelDrugDivWashThreshold)/numResponsesAbove ;
end



% phase, amplitude and transients change
for dset=1:3 ;
    for tp=1:l_tp ; % for each temporal period
        for sp=1:l_sp ; % for each spatial period
            for cntrst=1:l_cntrst ; % for each contrast
                for cells = 1:l_cells ; % for each cell
                    psth_var_divCntrl{dset}{tp}{cells}(cntrst,sp) = psth_var{dset}{tp}{cntrst}(cells,sp)/psth_var{1}{tp}{cntrst}(cells,sp) ; % condition/control  
                    psth_range_divCntrl{dset}{tp}{cells}(cntrst,sp) = psth_range{dset}{tp}{cntrst}(cells,sp)/psth_range{1}{tp}{cntrst}(cells,sp) ; % condition/control               
                    psth_mean_divCntrl{dset}{tp}{cells}(cntrst,sp) = psth_mean{dset}{tp}{cntrst}(cells,sp)/psth_mean{1}{tp}{cntrst}(cells,sp) ; % condition/control            
                    psth_duty_divCntrl{dset}{tp}{cells}(cntrst,sp) = psth_duty{dset}{tp}{cntrst}(cells,sp)/psth_duty{1}{tp}{cntrst}(cells,sp) ; % condition/control          
                    psth_peakDivMean_divCntrl{dset}{tp}{cells}(cntrst,sp) = psth_peakDivMean{dset}{tp}{cntrst}(cells,sp)/psth_peakDivMean{1}{tp}{cntrst}(cells,sp) ; % condition/control
                    
                    tempcc = xcov(psth{1}{tp}{cntrst}{sp}(cells,:),psth{dset}{tp}{cntrst}{sp}(cells,:),round(length(psthTemp)/2)) ; % keeping the max lag within half a cycle to prevent correlation with earlier phase
                    timecc = ([1:length(tempcc)]-ceil(length(tempcc)/2))/dataRun.sampling_rate ;
                    [mx,mi] = max(tempcc) ;
                    psth_cct_deltaCntrl{dset}{tp}{cells}(cntrst,sp) = timecc(mi) ; % negative if dset psth is trailing  
                end
            end
        end
    end 
end

% parameters of responses that have relatively large drug effects compared to wash
% histograms of parameter changes for responses that have relatively large drug effects compared to wash


% average parameters changes within cell types for responses that have relatively large drug effects compared to wash  


% parameters averaged across all cntrst and low frequencies stim
% average across contrasts and spatial frequencies within low frequency stimuli for drug change
for cells = 1:l_cells ; % for each cell
    psth_var_divCntrl_LfMean(cells) =  nanmean(nanmean(psth_var_divCntrl{2}{1}{cells}(:,lfSet))); 
    psth_range_divCntrl_LfMean(cells) = nanmean(nanmean(psth_range_divCntrl{2}{1}{cells}(:,lfSet))) ;
    psth_mean_divCntrl_LfMean(cells) =  nanmean(nanmean(psth_mean_divCntrl{2}{1}{cells}(:,lfSet))); 
    psth_duty_divCntrl_LfMean(cells) = nanmean(nanmean(psth_duty_divCntrl{2}{1}{cells}(:,lfSet))) ;
    psth_peakDivMean_divCntrl_LfMean(cells) = nanmean(nanmean(psth_peakDivMean_divCntrl{2}{1}{cells}(:,lfSet))) ;
    psth_cct_deltaCntrl_LfMean(cells) = nanmean(nanmean(psth_cct_deltaCntrl{2}{1}{cells}(:,lfSet))) ; 
end

% histogram for cells that had strong drug effects
[psth_var_divCntrl_LfMean_hist,psth_var_divCntrl_LfMean_histX] =  hist(psth_var_divCntrl_LfMean(MseRelDrug_DivMseRelWashLfMean>MseRelDrugDivWashThreshold),[0:.1:10]) ; 
[psth_range_divCntrl_LfMean_hist,psth_range_divCntrl_LfMean_histX] = hist(psth_range_divCntrl_LfMean(MseRelDrug_DivMseRelWashLfMean>MseRelDrugDivWashThreshold),[0:.1:10]) ;
[psth_mean_divCntrl_LfMean_hist,psth_mean_divCntrl_LfMean_histX] =  hist(psth_mean_divCntrl_LfMean(MseRelDrug_DivMseRelWashLfMean>MseRelDrugDivWashThreshold),[0:.1:10]) ;
[psth_duty_divCntrl_LfMean_hist,psth_duty_divCntrl_LfMean_histX] = hist(psth_duty_divCntrl_LfMean(MseRelDrug_DivMseRelWashLfMean>MseRelDrugDivWashThreshold),[0:.1:10]) ;
[psth_peakDivMean_divCntrl_LfMean_hist,psth_peakDivMean_divCntrl_LfMean_histX] = hist(psth_peakDivMean_divCntrl_LfMean(MseRelDrug_DivMseRelWashLfMean>MseRelDrugDivWashThreshold),[0:.1:10]) ;
[psth_cct_deltaCntrl_LfMean_hist,psth_cct_deltaCntrl_LfMean_histX] = hist(psth_cct_deltaCntrl_LfMean(MseRelDrug_DivMseRelWashLfMean>MseRelDrugDivWashThreshold),[-.5:.01:.5]) ; 

% change in spatial 
for cells = 1:l_cells ; % for each cell
    psth_var_divCntrl_HfMean(cells) =  nanmean(nanmean(psth_var_divCntrl{2}{1}{cells}(:,hfSet))); 
    psth_range_divCntrl_HfMean(cells) = nanmean(nanmean(psth_range_divCntrl{2}{1}{cells}(:,hfSet))) ;
end

% contrast gain as function of spatial frequency (finding cgain with linear fit)
for dset=1:3 ;
    for tp=1:l_tp ; % for each temporal period
        for sp=1:l_sp ; % for each spatial period
            for cells = 1:l_cells ; % for each cell
                for cntrst=1:l_cntrst ; % for each contrast
                    v1(cntrst) = psth_var{dset}{tp}{cntrst}(cells,sp) ;
                    r1(cntrst) = psth_range{dset}{tp}{cntrst}(cells,sp) ;
                end
                temp = polyfit(Contrast,v1,1) ;
                psth_var_cgain{dset}{tp}(cells,sp) = temp(1) ;
                
                temp = polyfit(Contrast,r1,1) ;
                psth_range_cgain{dset}{tp}(cells,sp) = temp(1) ;
            end
        end
    end
end

% contrast gain as function of spatial frequency (finding cgain with
% response/contrast)
for dset=1:3 ;
    for tp=1:l_tp ; % for each temporal period
        for cells = 1:l_cells ; % for each cell
            for sp=1:l_sp ; % for each spatial period
                for cntrst=1:l_cntrst ; % for each contrast
                    v1(cntrst,sp) = psth_var{dset}{tp}{cntrst}(cells,sp) ;
                    r1(cntrst,sp) = psth_range{dset}{tp}{cntrst}(cells,sp) ;
                end
            end
            vdm = abs(v1-max(v1(:))/2) ; % difference from half max of all responses
            rdm = abs(r1-max(r1(:))/2) ;
            
            for sp=1:l_sp ; % for each spatial period    
                [m,i] = min(vdm(:,sp)) ;
                psth_var_cgain{dset}{tp}(cells,sp) = v1(i,sp)/Contrast(i) ;
                
                [m,i] = min(rdm(:,sp)) ;
                psth_range_cgain{dset}{tp}(cells,sp) = r1(i,sp)/Contrast(i) ;
            end
        end
    end
end

for dset=1:3 ;
    for cells = 1:l_cells ; % for each cell
        psth_var_cgain_DivCntrl{dset}(cells,:) = psth_var_cgain{dset}{1}(cells,:)./psth_var_cgain{1}{1}(cells,:) ;
        psth_range_cgain_DivCntrl{dset}(cells,:) = psth_range_cgain{dset}{1}(cells,:)./psth_range_cgain{1}{1}(cells,:) ;
    end
end
    
% average across contrasts 
for dset=1:3 ;
    for tp=1:l_tp ; % for each temporal period
        for cells = 1:l_cells ; % for each cell
            va = zeros(1,l_sp) ;
            ra = zeros(1,l_sp) ;
            for cntrst=1:l_cntrst ; % for each contrast
                va = va + psth_var{dset}{tp}{cntrst}(cells,:) ;
                ra = ra + psth_range{dset}{tp}{cntrst}(cells,:) ;
            end
            psth_var_cMean{dset}{tp}(cells,:) = va/l_cntrst ;
            psth_range_cMean{dset}{tp}(cells,:) = ra/l_cntrst ;
        end
    end
end

% rename some datafields so can delete dataRun struct
dataRunCell_ids = dataRun.cell_ids ; % so can delete dataRun

%clear dataRun

% FIGURES

lineWidthRange = 3 ;
% psths of each cell in control, +psem, wash
cells=1 ;
while cells <=l_cells ; % for each cell
    figure(1)
    set(gcf,'name',[num2str(cells),' LfMse W ',num2str(ceil(MseRelWashLfMean(cells)*100)/100),...
                    ' D ',num2str(ceil(MseRelDrugLfMean(cells)*100)/100),...
                    ' D/W ',num2str(ceil(MseRelDrug_DivMseRelWashLfMean(cells)*100)/100)])
    clf
    for dset=1:3 ;
        for sp=1:l_sp ; % for each spatial period 
            for cntrst=1:l_cntrst ; % for each contrast
                subplot(l_sp,l_cntrst,l_cntrst*(sp-1)+cntrst) ;
                %plot(PsthTime,psth{dset}{1}{cntrst}{sp}(cells,:),'color',Color_list{dset},'linewidth',lineWidthRange*psth_deltaWash{cntrst}{sp}(cells)+.1)
                plot(PsthTime,psth{dset}{1}{cntrst}{sp}(cells,:),'color',Color_list{dset})
                axis tight
                hold on
                title(['Mse W ',num2str(ceil(MseRelWash{cells}(cntrst,sp)*100)/100),...
                    ' D ',num2str(ceil(MseRelDrug{cells}(cntrst,sp)*100)/100),...
                    'D/W ',num2str(ceil(MseRelDrug_DivMseRelWash{cells}(cntrst,sp)*100)/100)]) ;
            end
        end
    end
    
    figure(2)
    clf
    % master cell indicy
    masteri=[] ;
    for clm=1:length(cell_list_map) ;
        if cell_list_map{clm}==dataRunCell_ids(cells) ;
            masteri=clm ;
        end
    end
    if ~isempty(masteri) && ~isempty(srf{masteri}) ; % if there is a mapped master cell that has a rf
        subplot(8,l_cntrst,1) % spatial receptive field
        temp_rf = srf{masteri} ;
        norm_rf = norm_image(temp_rf);
        imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
        colormap(brewermap([],'RdBu'))
        caxis([0,1]) 
        set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])

        subplot(8,l_cntrst,2) % temporal receptive field
        plot(fliplr(trf_time),trf{masteri}) ;
    end
    for dset=1:3 ;
        for cntrst=1:l_cntrst ; % for each contrast

            subplot(8,l_cntrst,l_cntrst+cntrst) % psth variance (spatial frequency)
            plot(log(Spatial_frequency),psth_var{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            %errorbar(log(Spatial_frequency),psthEsem_var_mean{dset}{1}{cntrst}(cells,:),psthEsem_var_sem{dset}{1}{cntrst}(cells,:),psthEsem_var_sem{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth variance')
            hold on
            
            subplot(8,l_cntrst,l_cntrst*2+cntrst) % psth norm variance (spatial frequency)
            plot(log(Spatial_frequency),psth_var_norm{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('norm psth variance')
            %set(gca,'yscale','log')
            hold on
            
            subplot(8,l_cntrst,l_cntrst*3+cntrst) % psth range (spatial frequency)
            plot(log(Spatial_frequency),psth_range{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth range')
            hold on

            subplot(8,l_cntrst,l_cntrst*4+cntrst) % psth mean (spatial frequency)
            plot(log(Spatial_frequency),psth_mean{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth mean')
            hold on

            subplot(8,l_cntrst,l_cntrst*5+cntrst) % psth peak/mean  (spatial frquency)
            plot(log(Spatial_frequency),psth_peakDivMean{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth peak/mean')
            hold on
            
            subplot(8,l_cntrst,l_cntrst*6+cntrst) % psth duty  (spatial frquency)
            plot(log(Spatial_frequency),psth_duty{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth duty')
            hold on
            
            subplot(8,l_cntrst,l_cntrst*7+cntrst) % psth cc shift  (spatial frquency)
            plot(log(Spatial_frequency),psth_cct_deltaCntrl{dset}{1}{cells}(cntrst,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth cc peak time')
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

% mse 
figure
subplot(2,2,1)
for cells=1:l_cells ; % for each cell of this type 
    plot(MseRelWash{cells}(:),MseRelDrug{cells}(:),'k*') ;
    hold on
    plot([min(MseRelWash{cells}(:)),max(MseRelWash{cells}(:))],[min(MseRelWash{cells}(:)),max(MseRelWash{cells}(:))],'g')
    plot([MseRelThreshold,MseRelThreshold],[min(MseRelDrug{cells}(:)),max(MseRelDrug{cells}(:))],'r')
    plot([min(MseRelWash{cells}(:)),max(MseRelWash{cells}(:))],MseRelDrugDivWashThreshold*[min(MseRelWash{cells}(:)),max(MseRelWash{cells}(:))],'r')
    set(gca,'xscale','log')
    set(gca,'yscale','log')
end
axis tight
xlabel('Rel Mse Wash')
ylabel('Rel Mse Drug')
title('all responses')

subplot(2,2,2)
plot(MseRelWashLfMean,MseRelDrugLfMean,'k*') ;
hold on
plot([min(MseRelWashLfMean),max(MseRelWashLfMean)],[min(MseRelWashLfMean),max(MseRelWashLfMean)],'g')
plot([MseRelThreshold,MseRelThreshold],[min(MseRelDrugLfMean),max(MseRelDrugLfMean)],'r')
plot([min(MseRelWashLfMean),max(MseRelWashLfMean)],MseRelDrugDivWashThreshold*[min(MseRelWashLfMean),max(MseRelWashLfMean)],'r')
set(gca,'xscale','log')
set(gca,'yscale','log')
axis tight
xlabel('Rel Mse Wash')
ylabel('Rel Mse Drug')
title('cell mean lf')
 
subplot(2,2,3)
for cells=1:l_cells ; % for each cell of this type 
    plot(MseRelHistX,MseRelDrug_DivMseRelWash_Hist) ;
    text(.5,.9,['%aboveThresh=',num2str(100*MseRelDrug_DivMseRelWash_FractAboveThresh)],'units','norm')
end
xlabel('Rel Mse drug/wash')
ylabel('num obs')
title('all responses')

subplot(2,2,4)
for cells=1:l_cells ; % for each cell of this type 
    plot(MseRelHistX,MseRelDrug_DivMseRelWashLfMean_Hist) ;
    text(.5,.9,['%aboveThresh=',num2str(100*MseRelDrug_DivMseRelWashLfMean_FractAboveThresh)],'units','norm')
end
xlabel('Rel Mse drug/wash')
ylabel('num obs')
title('cell mean lf')

% does magintude of relative drug change depend on array location or cell type?
figure

subplot(2,1,1)
coord_tform = coordinate_transform(dataRunMaster,'sta');
            
for c = 1:length(dataRunMaster.cell_ids) ; % for each with a BW mapped ref
    if ~isempty(cell_list_map{c}) ; % if there is a mapped cell
        slave_c = get_cell_indices(dataRun,cell_list_map{c}) ;
        lw = MseRelDrug_DivMseRelWashLfMean(slave_c)/max(MseRelDrug_DivMseRelWashLfMean) ; 
        
        if ~isnan(lw) ;
            ctr = dataRunMaster.stas.fits{c}.mean ;
            rad = dataRunMaster.stas.fits{c}.sd ;
            angle = dataRunMaster.stas.fits{c}.angle ;
            [X,Y] = drawEllipse([ctr rad angle]) ;
            if ~any(isnan([X,Y])) ;
                [X,Y] = tformfwd(coord_tform, X, Y) ;
                plot(X,Y,'k','color',[1-lw,1-lw,1-lw],'linewidth',lw*3)
                hold on
            end
        end
    end
end   

subplot(2,1,2)
for uc = 1:length(UniqueCellTypes) ;
    tempMean = 0 ;
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        plot(uc,MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells)),'ko')
        if ~isnan(MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells))) 
            tempMean = tempMean+MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells))/lc ;
        end
        hold on
    end
    plot(uc,tempMean,'r+','MarkerSize',30)
end

% parameter comarisons for cells with relatively large drug changes
figure
for uc = 1:length(UniqueCellTypes) ;
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        if MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells))>MseRelDrugDivWashThreshold ; % if this cell had a drug effect    
            if strcmp('ON',UniqueCellTypes{uc}(1:2))
                SignPnt = 'o' ;
            elseif strcmp('OFF',UniqueCellTypes{uc}(1:3))
                SignPnt = '*' ;
            else
                SignPnt = '+' ;
            end

            subplot(6,2,1)
            plot(uc,psth_var_divCntrl_LfMean(cell_i{uc}(cells)),['k',SignPnt])
            hold on
            set(gca,'Yscale','log')
            xlabel('cell type')
            ylabel('Variance drug/cntrl')
            
            subplot(6,2,3)
            plot(uc,psth_range_divCntrl_LfMean(cell_i{uc}(cells)),['k',SignPnt])
            hold on
            set(gca,'Yscale','log')
            xlabel('cell type')
            ylabel('Range drug/cntrl')
            
            subplot(6,2,5)
            plot(uc,psth_mean_divCntrl_LfMean(cell_i{uc}(cells)),['k',SignPnt])
            hold on
            set(gca,'Yscale','log')
            xlabel('cell type')
            ylabel('Mean drug/cntrl')
            
            subplot(6,2,7)
            plot(uc,psth_duty_divCntrl_LfMean(cell_i{uc}(cells)),['k',SignPnt])
            hold on
            set(gca,'Yscale','log')
            xlabel('cell type')
            ylabel('Duty drug/cntrl')
            
            subplot(6,2,9)
            plot(uc,psth_peakDivMean_divCntrl_LfMean(cell_i{uc}(cells)),['k',SignPnt])
            hold on
            set(gca,'Yscale','log')
            xlabel('cell type')
            ylabel('peak/mean drug/cntrl')
            
            subplot(6,2,11)
            plot(uc,psth_cct_deltaCntrl_LfMean(cell_i{uc}(cells)),['k',SignPnt])
            hold on
            xlabel('cell type')
            ylabel('cc peak t drug,cntrl')
        end
    end
end
subplot(6,2,1); plot([1,length(UniqueCellTypes)],[1,1])
subplot(6,2,3); plot([1,length(UniqueCellTypes)],[1,1])
subplot(6,2,5); plot([1,length(UniqueCellTypes)],[1,1])
subplot(6,2,7); plot([1,length(UniqueCellTypes)],[1,1])
subplot(6,2,9); plot([1,length(UniqueCellTypes)],[1,1])
subplot(6,2,11); plot([1,length(UniqueCellTypes)],[0,0])

subplot(6,2,2)
plot(psth_var_divCntrl_LfMean_histX,psth_var_divCntrl_LfMean_hist)
set(gca,'Xscale','log')
xlabel('var drug/cntrl')
ylabel('# obs')

subplot(6,2,4)
plot(psth_range_divCntrl_LfMean_histX,psth_range_divCntrl_LfMean_hist)
set(gca,'Xscale','log')
xlabel('range drug/cntrl')
ylabel('# obs')

subplot(6,2,6)
plot(psth_mean_divCntrl_LfMean_histX,psth_mean_divCntrl_LfMean_hist)
set(gca,'Xscale','log')
xlabel('mean drug/cntrl')
ylabel('# obs')

subplot(6,2,8)
plot(psth_duty_divCntrl_LfMean_histX,psth_duty_divCntrl_LfMean_hist)
set(gca,'Xscale','log')
xlabel('duty drug/cntrl')
ylabel('# obs')

subplot(6,2,10)
plot(psth_peakDivMean_divCntrl_LfMean_histX,psth_peakDivMean_divCntrl_LfMean_hist)
set(gca,'Xscale','log')
xlabel('peak/mean drug/cntrl')
ylabel('# obs')

subplot(6,2,12)
plot(psth_cct_deltaCntrl_LfMean_histX,psth_cct_deltaCntrl_LfMean_hist)
xlabel('ccPeakTime drug/cntrl')
ylabel('# obs')


% mseGroups - consistency of drug and wash for all stim conditions 
subplot(3,3,2) ; % hist of groups
tempmat = cell2mat(MseGroup) ;
temphist = hist(tempmat(:),[0:4]) ;
plot([0:4],temphist/sum(temphist))

subplot(3,3,3) ; % fraction of responses in each cell in the same group
hist(MseGroup_sameFract,[0:.1:1]) % is each cell's response similarly changed in drug and wash?

subplot(3,3,4) ; % fraction of responses in each cell in either 1/2 or 3/4
hist(MseGroup_drugChangeFract,[0:.1:1]) % is each cell's response similarly changed in drug?

subplot(3,3,5) ; % fraction of responses in each cell in either 1/2 or 3/4
hist(MseGroup_stableFract,[0:.1:1]) % is each cell's response similarly changed in drug?

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

% spatial tuning for cells with strong lf Drug effects
figure 
for cells=1:l_cells ; % for each cell of this type
    subplot(2,1,1)
    if MseRelDrug_DivMseRelWashLfMean(cells)>MseRelDrugDivWashThreshold ;
        plot(psth_var_divCntrl_LfMean(cells),psth_var_divCntrl_HfMean(cells),'k*')
        hold on
    end
    set(gca,'Xscale','log')
    set(gca,'Yscale','log')
    xlabel('lf mean var drug/cntl')
    ylabel('hf mean var drug/cntl')
    
    subplot(2,1,2)
    if MseRelDrug_DivMseRelWashLfMean(cells)>MseRelDrugDivWashThreshold ;
        plot(psth_range_divCntrl_LfMean(cells),psth_range_divCntrl_HfMean(cells),'k*')
        hold on
    end
    set(gca,'Xscale','log')
    set(gca,'Yscale','log')
    xlabel('lf mean range drug/cntl')
    ylabel('hf mean range drug/cntl')      
end

subplot(2,1,1)
axis tight
ar = [get(gca,'Xlim'),get(gca,'YLim')] ;
plot([min(ar),max(ar)],[min(ar),max(ar)],'r')

subplot(2,1,2)
axis tight
ar = [get(gca,'Xlim'),get(gca,'YLim')] ;
plot([min(ar),max(ar)],[min(ar),max(ar)],'r')

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
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=1:3 ; 
            plot(Spatial_frequency,psth_var_cMean{dset}{1}(cell_i{uc}(cells),:),'color',Color_list{dset})
            hold on
        end
        set(gca,'Xscale','log')
        title(num2str(cell_i{uc}(cells)))
    end
end
    
        
%set(gca,'XScale','log')    



