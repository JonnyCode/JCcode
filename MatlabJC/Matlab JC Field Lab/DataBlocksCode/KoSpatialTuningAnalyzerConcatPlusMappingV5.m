function ForIgor = KoSpatialTuningAnalyzerConcatPlusMappingV5(DataBlock, DB, Params)

% this function is adapted from
% 'KoSpatialTuningAnalyzerConcatPlusMapping4'.  Edited to look at changes 
% in spatial tuning curves across contrasts 


% JC 9/26/2016 

% parameters
Color_list = {'k','r','b','g','c','v'} ; % order of colors for each 
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

Spatial_frequency = 1./(SpatialPeriod*MicronPerPix) ; % cycles/um

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
    else
        cell_i{uc} = [] ;
    end
end

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
        end
    end
    
    % average across contrasts and all low frequency stimuli
    lfSet = find(SpatialPeriod>=LfSpThresh) ;
    hfSet = find(SpatialPeriod<LfSpThresh) ;
    MseRelWashLfMean(cells) = nanmean(mean(MseRelWash{cells}(:,lfSet))) ; 
    MseRelDrugLfMean(cells) = nanmean(mean(MseRelDrug{cells}(:,lfSet))) ; 

    MseRelDrug_DivMseRelWashLfMean(cells) = MseRelDrugLfMean(cells)/MseRelWashLfMean(cells) ;
end

% all responses
MseRelHistX = [0:.1:10] ;
tempMat = cell2mat(MseRelDrug_DivMseRelWash) ;
MseRelDrug_DivMseRelWash_Hist = hist(tempMat(:),MseRelHistX) ;
MseRelDrug_DivMseRelWash_FractAboveThresh = sum(tempMat(:)>MseRelDrugDivWashThreshold)/length(tempMat(:)) ;
MseRelDrug_DivMseRelWash_Hist_cumsum = cumsum(MseRelDrug_DivMseRelWash_Hist)/sum(MseRelDrug_DivMseRelWash_Hist) ; 

% cell responses averaged at low sf
MseRelDrug_DivMseRelWashLfMean_Hist = hist(MseRelDrug_DivMseRelWashLfMean,MseRelHistX) ;
MseRelDrug_DivMseRelWashLfMean_Hist_cumsum = cumsum(MseRelDrug_DivMseRelWashLfMean_Hist)/sum(MseRelDrug_DivMseRelWashLfMean_Hist) ;   

MseRelWashLfMean_Hist = hist(MseRelWashLfMean,MseRelHistX) ;
MseRelWashLfMean_Hist_cumsum = cumsum(MseRelWashLfMean_Hist)/sum(MseRelWashLfMean_Hist) ;

MseRelDrugLfMean_Hist = hist(MseRelDrugLfMean,MseRelHistX) ;
MseRelDrugLfMean_Hist_cumsum = cumsum(MseRelDrugLfMean_Hist)/sum(MseRelDrugLfMean_Hist) ;

% fraction above and below thresholds
MseRel_FractAboveDrugDivWashThresh = sum(MseRelDrug_DivMseRelWashLfMean>MseRelDrugDivWashThreshold)/length(MseRelDrug_DivMseRelWashLfMean) ;
MseRel_FractBelowWashThresh = sum(MseRelWashLfMean<MseRelWashThreshold)/length(MseRelWashLfMean) ;
MseRel_FractBelowWashAndAboveDrugDivWashThresh = sum(MseRelWashLfMean<MseRelWashThreshold.*MseRelDrug_DivMseRelWashLfMean>MseRelDrugDivWashThreshold)/length(MseRelWashLfMean) ;

% phase, amplitude and transients change (a-b/abs(a)+abs(b))
for dset=1:3 ;
    for tp=1:l_tp ; % for each temporal period
        for sp=1:l_sp ; % for each spatial period
            for cntrst=1:l_cntrst ; % for each contrast
                for cells = 1:l_cells ; % for each cell
                    psth_var_divCntrl{dset}{tp}{cells}(cntrst,sp) = DiffOverSum(psth_var{dset}{tp}{cntrst}(cells,sp),psth_var{1}{tp}{cntrst}(cells,sp)) ; % condition,control  
                    psth_range_divCntrl{dset}{tp}{cells}(cntrst,sp) = DiffOverSum(psth_range{dset}{tp}{cntrst}(cells,sp),psth_range{1}{tp}{cntrst}(cells,sp)) ; % condition,control               
                    psth_mean_divCntrl{dset}{tp}{cells}(cntrst,sp) = DiffOverSum(psth_mean{dset}{tp}{cntrst}(cells,sp),psth_mean{1}{tp}{cntrst}(cells,sp)) ; % condition,control            
                    psth_duty_divCntrl{dset}{tp}{cells}(cntrst,sp) = DiffOverSum(psth_duty{dset}{tp}{cntrst}(cells,sp),psth_duty{1}{tp}{cntrst}(cells,sp)) ; % condition,control          
                    psth_peakDivMean_divCntrl{dset}{tp}{cells}(cntrst,sp) = DiffOverSum(psth_peakDivMean{dset}{tp}{cntrst}(cells,sp),psth_peakDivMean{1}{tp}{cntrst}(cells,sp)) ; % condition,control
                end
            end
        end
    end 
end

% drug and wash changes averaged across contrasts
for cells = 1:l_cells ; % for each cell
    psth_var_DrugdivCntrl_cMean(cells,:) =  nanmean(psth_var_divCntrl{2}{1}{cells}); 
    psth_range_DrugdivCntrl_cMean(cells,:) = nanmean(psth_range_divCntrl{2}{1}{cells}) ;
    psth_mean_DrugdivCntrl_cMean(cells,:) =  nanmean(psth_mean_divCntrl{2}{1}{cells}); 
    psth_duty_DrugdivCntrl_cMean(cells,:) = nanmean(psth_duty_divCntrl{2}{1}{cells}) ;
    psth_peakDivMean_DrugdivCntrl_cMean(cells,:) = nanmean(psth_peakDivMean_divCntrl{2}{1}{cells}) ;

    psth_var_WashdivCntrl_cMean(cells,:) =  nanmean(psth_var_divCntrl{3}{1}{cells}); 
    psth_range_WashdivCntrl_cMean(cells,:) = nanmean(psth_range_divCntrl{3}{1}{cells}) ;
    psth_mean_WashdivCntrl_cMean(cells,:) =  nanmean(psth_mean_divCntrl{3}{1}{cells}); 
    psth_duty_WashdivCntrl_cMean(cells,:) = nanmean(psth_duty_divCntrl{3}{1}{cells}) ;
    psth_peakDivMean_WashdivCntrl_cMean(cells,:) = nanmean(psth_peakDivMean_divCntrl{3}{1}{cells}) ;
end


% drug and wash changes averaged across cntrst and low frequencies stim
for cells = 1:l_cells ; % for each cell
    psth_var_DrugdivCntrl_LfMean(cells) =  nanmean(psth_var_DrugdivCntrl_cMean(cells,:)); 
    psth_range_DrugdivCntrl_LfMean(cells) = nanmean(psth_range_DrugdivCntrl_cMean(cells,:)) ;
    psth_mean_DrugdivCntrl_LfMean(cells) =  nanmean(psth_mean_DrugdivCntrl_cMean(cells,:)); 
    psth_duty_DrugdivCntrl_LfMean(cells) = nanmean(psth_duty_DrugdivCntrl_cMean(cells,:)) ;
    psth_peakDivMean_DrugdivCntrl_LfMean(cells) = nanmean(psth_peakDivMean_DrugdivCntrl_cMean(cells,:)) ;
    
    psth_var_WashdivCntrl_LfMean(cells) =  nanmean(psth_var_WashdivCntrl_cMean(cells,:)); 
    psth_range_WashdivCntrl_LfMean(cells) = nanmean(psth_range_WashdivCntrl_cMean(cells,:)) ;
    psth_mean_WashdivCntrl_LfMean(cells) =  nanmean(psth_mean_WashdivCntrl_cMean(cells,:)); 
    psth_duty_WashdivCntrl_LfMean(cells) = nanmean(psth_duty_WashdivCntrl_cMean(cells,:)) ;
    psth_peakDivMean_WashdivCntrl_LfMean(cells) = nanmean(psth_peakDivMean_WashdivCntrl_cMean(cells,:)) ; 
end

% histogram for cells that were stable
[psth_var_DrugdivCntrl_LfMean_hist,psth_var_divCntrl_LfMean_histX] =  hist(psth_var_DrugdivCntrl_LfMean(MseRelWashLfMean<MseRelWashThreshold),[-1.1:.1:1.1]) ; 
[psth_range_DrugdivCntrl_LfMean_hist,psth_range_divCntrl_LfMean_histX] = hist(psth_range_DrugdivCntrl_LfMean(MseRelWashLfMean<MseRelWashThreshold),[-1.1:.1:1.1]) ;
[psth_mean_DrugdivCntrl_LfMean_hist,psth_mean_divCntrl_LfMean_histX] =  hist(psth_mean_DrugdivCntrl_LfMean(MseRelWashLfMean<MseRelWashThreshold),[-1.1:.1:1.1]) ;
[psth_duty_DrugdivCntrl_LfMean_hist,psth_duty_divCntrl_LfMean_histX] = hist(psth_duty_DrugdivCntrl_LfMean(MseRelWashLfMean<MseRelWashThreshold),[-1.1:.1:1.1]) ;
[psth_peakDivMean_DrugdivCntrl_LfMean_hist,psth_peakDivMean_divCntrl_LfMean_histX] = hist(psth_peakDivMean_DrugdivCntrl_LfMean(MseRelWashLfMean<MseRelWashThreshold),[-1.1:.1:1.1]) ;

[psth_var_WashdivCntrl_LfMean_hist,psth_var_divCntrl_LfMean_histX] =  hist(psth_var_WashdivCntrl_LfMean(MseRelWashLfMean<MseRelWashThreshold),[-1.1:.1:1.1]) ; 
[psth_range_WashdivCntrl_LfMean_hist,psth_range_divCntrl_LfMean_histX] = hist(psth_range_WashdivCntrl_LfMean(MseRelWashLfMean<MseRelWashThreshold),[-1.1:.1:1.1]) ;
[psth_mean_WashdivCntrl_LfMean_hist,psth_mean_divCntrl_LfMean_histX] =  hist(psth_mean_WashdivCntrl_LfMean(MseRelWashLfMean<MseRelWashThreshold),[-1.1:.1:1.1]) ;
[psth_duty_WashdivCntrl_LfMean_hist,psth_duty_divCntrl_LfMean_histX] = hist(psth_duty_WashdivCntrl_LfMean(MseRelWashLfMean<MseRelWashThreshold),[-1.1:.1:1.1]) ;
[psth_peakDivMean_WashdivCntrl_LfMean_hist,psth_peakDivMean_divCntrl_LfMean_histX] = hist(psth_peakDivMean_WashdivCntrl_LfMean(MseRelWashLfMean<MseRelWashThreshold),[-1.1:.1:1.1]) ;

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
        psth_var_cgain_DivCntrl{dset}(cells,:) = DiffOverSum(psth_var_cgain{dset}{1}(cells,:),psth_var_cgain{1}{1}(cells,:)) ;
        psth_range_cgain_DivCntrl{dset}(cells,:) = DiffOverSum(psth_range_cgain{dset}{1}(cells,:),psth_range_cgain{1}{1}(cells,:)) ;
    end
end
    
% average parameters across contrasts 
for dset=1:3 ;
    for tp=1:l_tp ; % for each temporal period
        for cells = 1:l_cells ; % for each cell
            va = zeros(1,l_sp) ;
            ra = zeros(1,l_sp) ;
            mn = zeros(1,l_sp) ;
            dty = zeros(1,l_sp) ;
            pkdmn = zeros(1,l_sp) ;
            for cntrst=1:l_cntrst ; % for each contrast
                va = va + psth_var{dset}{tp}{cntrst}(cells,:) ;
                ra = ra + psth_range{dset}{tp}{cntrst}(cells,:) ;
                mn = mn + psth_mean{dset}{tp}{cntrst}(cells,:) ;
                dty = dty + psth_duty{dset}{tp}{cntrst}(cells,:) ;
                pkdmn = pkdmn + psth_peakDivMean{dset}{tp}{cntrst}(cells,:) ;
            end
            psth_var_cMean{dset}{tp}(cells,:) = va/l_cntrst ;
            psth_range_cMean{dset}{tp}(cells,:) = ra/l_cntrst ;
            psth_mean_cMean{dset}{tp}(cells,:) = mn/l_cntrst ;
            psth_dty_cMean{dset}{tp}(cells,:) = dty/l_cntrst ;
            psth_pkdmn_cMean{dset}{tp}(cells,:) = pkdmn/l_cntrst ;
        end
    end
end

% average within cell types
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    cellsStable{uc} = find(MseRelWashLfMean(cell_i{uc})<MseRelWashThreshold) ; % cells that were x stable
    cellsEffect{uc} = find(MseRelDrug_DivMseRelWashLfMean(cell_i{uc})<MseRelDrugDivWashThreshold) ; % cells that had a x times bigger drug change
    cellsStableEffect{uc} = find(MseRelWashLfMean(cell_i{uc})<MseRelWashThreshold & MseRelDrug_DivMseRelWashLfMean(cell_i{uc})<MseRelDrugDivWashThreshold) ; % cells
    
    % drug 
    psth_var_DrugdivCntrl_cMean_AllCells(uc,:) = nanmean(psth_var_DrugdivCntrl_cMean(cell_i{uc},:),1) ; % all cells of type uc
    psth_var_DrugdivCntrl_cMean_StableCells(uc,:) = nanmean(psth_var_DrugdivCntrl_cMean(cell_i{uc}(cellsStable{uc}),:),1) ; % cells of type uc that pass mse threshold
    psth_var_DrugdivCntrl_cMean_EffectCells(uc,:) = nanmean(psth_var_DrugdivCntrl_cMean(cell_i{uc}(cellsEffect{uc}),:),1) ;
    psth_var_DrugdivCntrl_cMean_StableEffectCells(uc,:) = nanmean(psth_var_DrugdivCntrl_cMean(cell_i{uc}(cellsStableEffect{uc}),:),1) ;
    
    psth_mean_DrugdivCntrl_cMean_AllCells(uc,:) = nanmean(psth_mean_DrugdivCntrl_cMean(cell_i{uc},:),1) ; % all cells of type uc
    psth_mean_DrugdivCntrl_cMean_StableCells(uc,:) = nanmean(psth_mean_DrugdivCntrl_cMean(cell_i{uc}(cellsStable{uc}),:),1) ;
    psth_mean_DrugdivCntrl_cMean_EffectCells(uc,:) = nanmean(psth_mean_DrugdivCntrl_cMean(cell_i{uc}(cellsEffect{uc}),:),1) ;
    psth_mean_DrugdivCntrl_cMean_StableEffectCells(uc,:) = nanmean(psth_mean_DrugdivCntrl_cMean(cell_i{uc}(cellsStableEffect{uc}),:),1) ;
    
    psth_duty_DrugdivCntrl_cMean_AllCells(uc,:) = nanmean(psth_duty_DrugdivCntrl_cMean(cell_i{uc},:),1) ; % all cells of type uc
    psth_duty_DrugdivCntrl_cMean_StableCells(uc,:) = nanmean(psth_duty_DrugdivCntrl_cMean(cell_i{uc}(cellsStable{uc}),:),1) ;
    psth_duty_DrugdivCntrl_cMean_EffectCells(uc,:) = nanmean(psth_duty_DrugdivCntrl_cMean(cell_i{uc}(cellsEffect{uc}),:),1) ;
    psth_duty_DrugdivCntrl_cMean_StableEffectCells(uc,:) = nanmean(psth_duty_DrugdivCntrl_cMean(cell_i{uc}(cellsStableEffect{uc}),:),1) ;
    
    % wash
    psth_var_WashdivCntrl_cMean_AllCells(uc,:) = nanmean(psth_var_WashdivCntrl_cMean(cell_i{uc},:),1) ; % all cells of type uc
    psth_var_WashdivCntrl_cMean_StableCells(uc,:) = nanmean(psth_var_WashdivCntrl_cMean(cell_i{uc}(cellsStable{uc}),:),1) ; % cells of type uc that pass mse threshold
    psth_var_WashdivCntrl_cMean_EffectCells(uc,:) = nanmean(psth_var_WashdivCntrl_cMean(cell_i{uc}(cellsEffect{uc}),:),1) ;
    psth_var_WashdivCntrl_cMean_StableEffectCells(uc,:) = nanmean(psth_var_WashdivCntrl_cMean(cell_i{uc}(cellsStableEffect{uc}),:),1) ;
    
    psth_mean_WashdivCntrl_cMean_AllCells(uc,:) = nanmean(psth_mean_WashdivCntrl_cMean(cell_i{uc},:),1) ; % all cells of type uc
    psth_mean_WashdivCntrl_cMean_StableCells(uc,:) = nanmean(psth_mean_WashdivCntrl_cMean(cell_i{uc}(cellsStable{uc}),:),1) ;
    psth_mean_WashdivCntrl_cMean_EffectCells(uc,:) = nanmean(psth_mean_WashdivCntrl_cMean(cell_i{uc}(cellsEffect{uc}),:),1) ;
    psth_mean_WashdivCntrl_cMean_StableEffectCells(uc,:) = nanmean(psth_mean_WashdivCntrl_cMean(cell_i{uc}(cellsStableEffect{uc}),:),1) ;
    
    psth_duty_WashdivCntrl_cMean_AllCells(uc,:) = nanmean(psth_duty_WashdivCntrl_cMean(cell_i{uc},:),1) ; % all cells of type uc
    psth_duty_WashdivCntrl_cMean_StableCells(uc,:) = nanmean(psth_duty_WashdivCntrl_cMean(cell_i{uc}(cellsStable{uc}),:),1) ;  
    psth_duty_WashdivCntrl_cMean_EffectCells(uc,:) = nanmean(psth_duty_WashdivCntrl_cMean(cell_i{uc}(cellsEffect{uc}),:),1) ;
    psth_duty_WashdivCntrl_cMean_StableEffectCells(uc,:) = nanmean(psth_duty_WashdivCntrl_cMean(cell_i{uc}(cellsStableEffect{uc}),:),1) ;
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

% variance changes from lowest to highest contrast used
for dset=1:3 ;
    for cells = 1:l_cells ; % for each cell
        for sp=1:l_sp ; % for each spatial period 
            psth_var_cntrstChangeIndex{dset}{1}{cntrst}(cells,sp) = DiffOverSum(psth_var{dset}{1}{end}(cells,sp),psth_var{dset}{1}{1}(cells,sp)) ;
        end
    end
end

clear dataRun dataRunMaster
save(['/Users/jcafaro/Desktop/Temp/Matfiles/KoSpatialV4',num2str(DB)],'-v7.3')

% FIGURES
if FigureFlag ;
    
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
                plot(PsthTime,psth{dset}{1}{cntrst}{sp}(cells,:),'color',Color_list{dset})
                axis tight
                hold on
                title(['Mse W ',num2str(ceil(MseRelWash{cells}(cntrst,sp)*100)/100),...
                    ' D ',num2str(ceil(MseRelDrug{cells}(cntrst,sp)*100)/100),...
                    ' D/W ',num2str(ceil(MseRelDrug_DivMseRelWash{cells}(cntrst,sp)*100)/100)]) ;
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
        subplot(8,l_cntrst,1) % spatial receptive field
        temp_rf = srf{masteri} ;
        norm_rf = norm_image(temp_rf);
        imagesc(matrix_scaled_up(norm_rf(:,:,1),8))
        colormap(brewermap([],'RdBu'))
        caxis([0,1]) 
        set(gca,'XTickLabel','','YTickLabel','','xtick', [], 'ytick', [])
        title(num2str(dataRunMaster_cell_ids(masteri)))
        
        subplot(8,l_cntrst,2) % temporal receptive field
        plot(fliplr(trf_time),trf{masteri}) ;
        title(cell_type_list{masteri})
    end
    for dset=1:3 ;
        for cntrst=1:l_cntrst ; % for each contrast

            subplot(8,l_cntrst,l_cntrst+cntrst) % psth variance (spatial frequency)
            plot(log(Spatial_frequency),psth_var{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            %errorbar(log(Spatial_frequency),psthEsem_var_mean{dset}{1}{cntrst}(cells,:),psthEsem_var_sem{dset}{1}{cntrst}(cells,:),psthEsem_var_sem{dset}{1}{cntrst}(cells,:),'color',Color_list{dset})
            xlabel('log 2 spatial f')
            ylabel('psth variance')
            title(num2str(dataRun_cell_ids(cells)))
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

% psth [control, +PSEM, wash] organized by unique cell types for a single spatial freq and contrast
sp = 7 ; % spatial frequency default
cntrst = 2 ; % contrast default
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

% mse (how does drug change of psth compare to wash change of psth?)
figure
subplot(2,2,1)
for cells=1:l_cells ; % for each cell of this type 
    plot(MseRelWash{cells}(:),MseRelDrug{cells}(:),'k*') ;
    hold on
    plot([min(MseRelWash{cells}(:)),max(MseRelWash{cells}(:))],[min(MseRelWash{cells}(:)),max(MseRelWash{cells}(:))],'g')
    plot([MseRelWashThreshold,MseRelWashThreshold],[min(MseRelDrug{cells}(:)),max(MseRelDrug{cells}(:))],'r')
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
plot([MseRelWashThreshold,MseRelWashThreshold],[min(MseRelDrugLfMean),max(MseRelDrugLfMean)],'r')
plot([min(MseRelWashLfMean),max(MseRelWashLfMean)],MseRelDrugDivWashThreshold*[min(MseRelWashLfMean),max(MseRelWashLfMean)],'r')
set(gca,'xscale','log')
set(gca,'yscale','log')
axis tight
xlabel('Rel Mse Wash')
ylabel('Rel Mse Drug')
title('cell mean lf')
 
subplot(2,2,3)
for cells=1:l_cells ; % for each cell of this type 
    plotyy(MseRelHistX,MseRelDrug_DivMseRelWash_Hist,MseRelHistX,MseRelDrug_DivMseRelWash_Hist_cumsum) ;
    text(.5,.9,['%aboveThresh=',num2str(100*MseRelDrug_DivMseRelWash_FractAboveThresh)],'units','norm')
end
xlabel('Rel Mse drug/wash')
ylabel('num obs')
title('all responses')

subplot(2,2,4)
for cells=1:l_cells ; % for each cell of this type 
    plotyy(MseRelHistX,MseRelDrug_DivMseRelWashLfMean_Hist,MseRelHistX,MseRelDrug_DivMseRelWashLfMean_Hist_cumsum) ;
    text(.5,.9,['%aboveThresh=',num2str(100*MseRel_FractAboveDrugDivWashThresh)],'units','norm')
end
xlabel('Rel Mse drug/wash')
ylabel('num obs')
title('cell mean lf')

% mse (does magintude of relative drug change depend on array location or cell type?)
figure

subplot(2,1,1)
for c = 1:length(dataRunMaster_cell_ids) ; % for each with a BW mapped ref
    if ~isempty(cell_list_map{c}) ; % if there is a mapped cell
        lw = MseRelDrug_DivMseRelWashLfMean(slave_c(c))/max(MseRelDrug_DivMseRelWashLfMean) ; 
        
        if ~isnan(lw) ;
            [X,Y] = drawEllipse([ctr{c} rad{c} angle{c}]) ;
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

% LfMean parameter comparisons for cells with (black) and without (cyan) relatively stable cells
figure
for uc = 1:length(UniqueCellTypes) ;
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        if MseRelWashLfMean(cell_i{uc}(cells))<MseRelWashThreshold ; % if this cell had a drug effect
            PlotColor = 'k' ;
        else
            PlotColor = 'c' ;
        end
        
        if strcmp('ON',UniqueCellTypes{uc}(1:2))
            SignPnt = 'o' ;
        elseif strcmp('OF',UniqueCellTypes{uc}(1:2))
            SignPnt = '*' ;
        else
            SignPnt = '+' ;
        end
        
        subplot(6,2,1)
        plot(uc,psth_var_DrugdivCntrl_LfMean(cell_i{uc}(cells)),[PlotColor,SignPnt])
        hold on
        xlabel('cell type')
        ylabel('Variance drug/cntrl')

        subplot(6,2,3)
        plot(uc,psth_range_DrugdivCntrl_LfMean(cell_i{uc}(cells)),[PlotColor,SignPnt])
        hold on
        xlabel('cell type')
        ylabel('Range drug/cntrl')

        subplot(6,2,5)
        plot(uc,psth_mean_DrugdivCntrl_LfMean(cell_i{uc}(cells)),[PlotColor,SignPnt])
        hold on
        xlabel('cell type')
        ylabel('Mean drug/cntrl')

        subplot(6,2,7)
        plot(uc,psth_duty_DrugdivCntrl_LfMean(cell_i{uc}(cells)),[PlotColor,SignPnt])
        hold on
        xlabel('cell type')
        ylabel('Duty drug/cntrl')

        subplot(6,2,9)
        plot(uc,psth_peakDivMean_DrugdivCntrl_LfMean(cell_i{uc}(cells)),[PlotColor,SignPnt])
        hold on
        xlabel('cell type')
        ylabel('peak/mean drug/cntrl')
    end
end
subplot(6,2,1); plot([1,length(UniqueCellTypes)],[0,0])
subplot(6,2,3); plot([1,length(UniqueCellTypes)],[0,0])
subplot(6,2,5); plot([1,length(UniqueCellTypes)],[0,0])
subplot(6,2,7); plot([1,length(UniqueCellTypes)],[0,0])
subplot(6,2,9); plot([1,length(UniqueCellTypes)],[0,0])
subplot(6,2,11); plot([1,length(UniqueCellTypes)],[0,0])

subplot(6,2,2)
plot(psth_var_divCntrl_LfMean_histX,psth_var_DrugdivCntrl_LfMean_hist)
xlabel('var drug/cntrl')
ylabel('# obs')

subplot(6,2,4)
plot(psth_range_divCntrl_LfMean_histX,psth_range_DrugdivCntrl_LfMean_hist)
xlabel('range drug/cntrl')
ylabel('# obs')

subplot(6,2,6)
plot(psth_mean_divCntrl_LfMean_histX,psth_mean_DrugdivCntrl_LfMean_hist)
xlabel('mean drug/cntrl')
ylabel('# obs')

subplot(6,2,8)
plot(psth_duty_divCntrl_LfMean_histX,psth_duty_DrugdivCntrl_LfMean_hist)
xlabel('duty drug/cntrl')
ylabel('# obs')

subplot(6,2,10)
plot(psth_peakDivMean_divCntrl_LfMean_histX,psth_peakDivMean_DrugdivCntrl_LfMean_hist)
xlabel('peak/mean drug/cntrl')
ylabel('# obs')

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

% parameter changes as function of sp averaged across contrast
% organized by cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        PlotColor = 'c' ;
        if MseRelWashLfMean(cell_i{uc}(cells))<MseRelWashThreshold ; % if this cell was stable
            PlotColor = 'b' ;
        end
        if MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells))>MseRelDrugDivWashThreshold ; % if this cell had a drug change > wash change
            PlotColor = 'r' ;
        end
        if MseRelWashLfMean(cell_i{uc}(cells))<MseRelWashThreshold & MseRelDrug_DivMseRelWashLfMean(cell_i{uc}(cells))>MseRelDrugDivWashThreshold ; % if this cell was stable and had a bigger drug effect
            PlotColor = 'k' ;
        end
        
        subplot(5,1,1) 
        plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean(cell_i{uc}(cells),:), PlotColor)
        hold on
        
        subplot(5,1,2) 
        plot(Spatial_frequency,psth_var_cgain_DivCntrl{2}(cell_i{uc}(cells),:), PlotColor)
        hold on
        
        subplot(5,1,3) 
        plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean(cell_i{uc}(cells),:), PlotColor)
        hold on

        subplot(5,1,4) 
        plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean(cell_i{uc}(cells),:), PlotColor)
        hold on
    end
    
    subplot(5,1,1) 
    %plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean_AllCells(uc,:),'b','linewidth',2)
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'r:')
    set(gca,'Xscale','log')
    axis tight
    ylabel('var drug/cntrl')
    xlabel('Spatial Frequency (c/um)')

    subplot(5,1,2) 
    %plot(Spatial_frequency,psth_var_cgain_DrugDivCntrl_AllCells(uc,:),'b','linewidth',2)
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'r:')
    set(gca,'Xscale','log')
    axis tight
    ylabel('var gain drug/cntrl')
    xlabel('Spatial Frequency (c/um)')

    subplot(5,1,3) 
    %plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean_AllCells(uc,:),'b','linewidth',2)
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'r:')
    set(gca,'Xscale','log')
    axis tight
    ylabel('mean drug/cntrl')
    xlabel('Spatial Frequency (c/um)')

    subplot(5,1,4) 
    %plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean_AllCells(uc,:),'b','linewidth',2)
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'r:')
    set(gca,'Xscale','log')
    axis tight
    ylabel('duty drug/cntrl')
    xlabel('Spatial Frequency (c/um)')
end

% parameter changes as function of sp averaged across contrast
% comparing averages across cell types
figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    subplot(6,1,1) 
    plot(Spatial_frequency,nanmean(psth_var_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight

    subplot(6,1,2) 
    plot(Spatial_frequency,nanmean(psth_range_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight

    subplot(6,1,3) 
    plot(Spatial_frequency,nanmean(psth_mean_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight

    subplot(6,1,4) 
    plot(Spatial_frequency,nanmean(psth_duty_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight

    subplot(6,1,5) 
    plot(Spatial_frequency,nanmean(psth_peakDivMean_DrugdivCntrl_cMean(cell_i{uc},:)),'linewidth',2)
    hold on
    plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g:')
    set(gca,'Xscale','log')
    axis tight
end

% average stats for all cell types for 'stable cells' 
NumFigRows = 6 ;
figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    
    subplot(NumFigRows,length(UniqueCellTypes),uc) ; % trfs 
    for c = 1:length(Master_i{uc}) ; % for each cell
        [X,Y] = drawEllipse([ctr{Master_i{uc}(c)} rad{Master_i{uc}(c)} angle{Master_i{uc}(c)}]) ;
        if ~any(isnan([X,Y])) ;
            [X,Y] = tformfwd(coord_tform, X, Y) ;
            if isempty(cell_list_map{Master_i{uc}(c)})    
                plot(X,Y,'k')
            else
                plot(X,Y,'k')
            end
            hold on
        end
    end   
    title(['N=',num2str(length(cell_ids{uc})),'/',num2str(length(Master_i{uc}))])
    
    if cell_ids{uc}>0 ;
        Nsqrt = sqrt(length(cell_ids{uc})) ;

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*1+uc) ; % trfs
        plot(trf_time,trf_uc_mean(uc,:))
        axis tight
        title(UniqueCellTypes{uc}) 
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*2+uc) ; % variance (spatial frequecy)
        errorbar(Spatial_frequency,mean(psth_var_cMean{1}{1}(cell_i{uc}(cellsStable{uc}),:),1),std(psth_var_cMean{1}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{1}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,'color',Color_list{1})
        hold on
        errorbar(Spatial_frequency,mean(psth_var_cMean{2}{1}(cell_i{uc}(cellsStable{uc}),:),1),std(psth_var_cMean{2}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{2}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,'color',Color_list{2})
        errorbar(Spatial_frequency,mean(psth_var_cMean{3}{1}(cell_i{uc}(cellsStable{uc}),:),1),std(psth_var_cMean{3}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{3}{1}(cell_i{uc}(cellsStable{uc}),:),[],1)/Nsqrt,'color',Color_list{3})
        set(gca,'Xscale','log')
        axis tight
        ylabel('var (sp/s)^2')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*3+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_var_WashdivCntrl_cMean_StableCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean_StableCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('var delta')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*4+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_mean_WashdivCntrl_cMean_StableCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean_StableCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('mean delta')
        xlabel('Sf (c/um)')
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*5+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_duty_WashdivCntrl_cMean_StableCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean_StableCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('duty delta')
        xlabel('Sf (c/um)')
    end
end
        
% cells with large Drug mse /wash mse 
NumFigRows = 6 ;
figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    
    subplot(NumFigRows,length(UniqueCellTypes),uc) ; % trfs
    for c = 1:length(Master_i{uc}) ; % for each cell
        [X,Y] = drawEllipse([ctr{Master_i{uc}(c)} rad{Master_i{uc}(c)} angle{Master_i{uc}(c)}]) ;
        if ~any(isnan([X,Y])) ;
            [X,Y] = tformfwd(coord_tform, X, Y) ;
            if isempty(cell_list_map{Master_i{uc}(c)})    
                plot(X,Y,'k')
            else
                plot(X,Y,'r')
            end
            hold on
        end
    end   
    title(['N=',num2str(length(cell_ids{uc})),'/',num2str(length(Master_i{uc}))])
    
    if cell_ids{uc}>0 ;
        Nsqrt = sqrt(length(cell_ids{uc})) ;

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*1+uc) ; % trfs
        plot(trf_time,trf_uc_mean(uc,:))
        axis tight
        title(UniqueCellTypes{uc}) 
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*2+uc) ; % variance (spatial frequecy)
        errorbar(Spatial_frequency,mean(psth_var_cMean{1}{1}(cell_i{uc}(cellsEffect{uc}),:),1),std(psth_var_cMean{1}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{1}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{1})
        hold on
        errorbar(Spatial_frequency,mean(psth_var_cMean{2}{1}(cell_i{uc}(cellsEffect{uc}),:),1),std(psth_var_cMean{2}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{2}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{2})
        errorbar(Spatial_frequency,mean(psth_var_cMean{3}{1}(cell_i{uc}(cellsEffect{uc}),:),1),std(psth_var_cMean{3}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{3}{1}(cell_i{uc}(cellsEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{3})
        set(gca,'Xscale','log')
        axis tight
        ylabel('var (sp/s)^2')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*3+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_var_WashdivCntrl_cMean_EffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean_EffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('var delta')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*4+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_mean_WashdivCntrl_cMean_EffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean_EffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('mean delta')
        xlabel('Sf (c/um)')
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*5+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_duty_WashdivCntrl_cMean_EffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean_EffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('duty delta')
        xlabel('Sf (c/um)')
    end
end

% cells that are stable and have large drug/wash mse
NumFigRows = 6 ;
figure
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    
    subplot(NumFigRows,length(UniqueCellTypes),uc) ; % trfs
    for c = 1:length(Master_i{uc}) ; % for each cell
        [X,Y] = drawEllipse([ctr(Master_i{uc}(c)) rad(Master_i{uc}(c)) angle(Master_i{uc}(c))]) ;
        if ~any(isnan([X,Y])) ;
            [X,Y] = tformfwd(coord_tform, X, Y) ;
            if isempty(cell_list_map{Master_i(c)})    
                plot(X,Y,'k')
            else
                plot(X,Y,'r')
            end
            hold on
        end
    end   
    title(['N=',num2str(length(cell_ids{uc})),'/',num2str(length(Master_i{uc}))])
    
    if cell_ids{uc}>0 ;
        Nsqrt = sqrt(length(cell_ids{uc})) ;

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*1+uc) ; % trfs
        plot(trf_time,trf_uc_mean(uc,:))
        axis tight
        title(UniqueCellTypes{uc}) 
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*2+uc) ; % variance (spatial frequecy)
        errorbar(Spatial_frequency,mean(psth_var_cMean{1}{1}(cell_i{uc}(cellsStableEffect{uc}),:),1),std(psth_var_cMean{1}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{1}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{1})
        hold on
        errorbar(Spatial_frequency,mean(psth_var_cMean{2}{1}(cell_i{uc}(cellsStableEffect{uc}),:),1),std(psth_var_cMean{2}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{2}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{2})
        errorbar(Spatial_frequency,mean(psth_var_cMean{3}{1}(cell_i{uc}(cellsStableEffect{uc}),:),1),std(psth_var_cMean{3}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,std(psth_var_cMean{3}{1}(cell_i{uc}(cellsStableEffect{uc}),:),[],1)/Nsqrt,'color',Color_list{3})
        set(gca,'Xscale','log')
        axis tight
        ylabel('var (sp/s)^2')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*3+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_var_WashdivCntrl_cMean_StableEffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_var_DrugdivCntrl_cMean_StableEffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('var delta')
        xlabel('Sf (c/um)')

        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*4+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_mean_WashdivCntrl_cMean_StableEffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_mean_DrugdivCntrl_cMean_StableEffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('mean delta')
        xlabel('Sf (c/um)')
        
        subplot(NumFigRows,length(UniqueCellTypes),length(UniqueCellTypes)*5+uc) ; % variance change (spatial frequecy)
        plot(Spatial_frequency,psth_duty_WashdivCntrl_cMean_StableEffectCells(uc,:),'b')
        hold on
        plot(Spatial_frequency,psth_duty_DrugdivCntrl_cMean_StableEffectCells(uc,:),'r')
        plot([Spatial_frequency(1),Spatial_frequency(end)],[0,0],'g--')
        set(gca,'Xscale','log')
        axis tight
        set(gca,'ylim',[-1,1])
        ylabel('duty delta')
        xlabel('Sf (c/um)')
    end
end
    
% plot mosaics and trfs for each cell type
for uc = 1:length(UniqueCellTypes) ; % for each cell type
     figure
     for c = 1:length(Master_i{uc}) ; % for each cell
        [X,Y] = drawEllipse([ctr{Master_i{uc}(c)} rad{Master_i{uc}(c)} angle{Master_i{uc}(c)}]) ;
        if ~any(isnan([X,Y])) ;
            subplot(2,1,1)
            [X,Y] = tformfwd(coord_tform, X, Y) ;
            if isempty(cell_list_map{Master_i{uc}(c)})    
                plot(X,Y,'k')
            else
                plot(X,Y,'k')
            end
            hold on

            subplot(2,1,2)
            plot(fliplr(trf_time),trf{Master_i{uc}(c)},'k')
            hold on
        end
     end  
end
    
% comparing spatial tuning curves for all contrasts (normalize by sum of spatial tuning curve)
for uc = 1:length(UniqueCellTypes) ; % for each cell type
     figure
     set(gcf,'name',UniqueCellTypes{uc})
     lc = length(cell_i{uc}) ; % number of cell of that type
     for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=1:3 ;
            for cntrst=1:l_cntrst ; % for each contrast 
                TempNormPlot = psth_var{dset}{1}{cntrst}(cell_i{uc}(cells),:)/sum(psth_var{dset}{1}{cntrst}(cell_i{uc}(cells),:)) ;
                plot(log(Spatial_frequency),TempNormPlot,'color',Color_list{dset},'linewidth',cntrst)
                hold on
            end
        end
     end
end

% comparing spatial tuning curves for all contrasts (normalize by sum of spatial tuning curve)
for uc = 1:length(UniqueCellTypes) ; % for each cell type
     figure
     set(gcf,'name',UniqueCellTypes{uc})
     lc = length(cell_i{uc}) ; % number of cell of that type
     for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=1:3 ;
            plot(log(Spatial_frequency),psth_var_cntrstChangeIndex{dset}{1}{cntrst}(cell_i{uc}(cells),:),'color',Color_list{dset})
            hold on
        end
     end
end


figure
for dset=1:3 ;
    for cntrst=1:l_cntrst ; % for each contrast 
%         for uc = 1:length(UniqueCellTypes) ; % for each cell type
%             subplot(ceil((length(UniqueCellTypes)+1)/3),3,uc) 
%             TempNormPlot = mean(psth_var{dset}{1}{cntrst}(cell_i{uc},:))/sum(mean(psth_var{dset}{1}{cntrst}(cell_i{uc},:))) ;
%             plot(log(Spatial_frequency),TempNormPlot,'color',Color_list{dset},'linewidth',cntrst)
%             hold on
%             title(UniqueCellTypes{uc})
%             axis tight
%         end
%             
%         subplot(ceil((length(UniqueCellTypes)+1)/3),3,length(UniqueCellTypes)+1)  
        TempNormPlot = mean(psth_var{dset}{1}{cntrst})/sum(mean(psth_var{dset}{1}{cntrst})) ;
        plot(log(Spatial_frequency),TempNormPlot,'color',Color_list{dset},'linewidth',cntrst)
        hold on
        title('all')
        axis tight
    end
end


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



