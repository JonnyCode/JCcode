function ForIgor = KoSpatialTuningAnalyzerConcat(DataBlock, DB, Params)

% this function will analyze drfiting grating responses from array data
% before and after cells are knocked out (KO)

% This function assumes that the concatinated data set is a concatinated set of
% [control,drug,wash].  And that those sets are also present as individual sets in DfGrating.

% JC 12/9/2015 

% parameters
Color_list = {'k','r','b','g','c','v'} ; % order of colors for each 
TrialTrigInterval = 10 ; % sec
InterTrialTime = 2 ; % sec
StimTransTime = .1 ; % (sec) skip cycles before this time after starting new stimulus
bufferTime = .01 ; % (sec) time after min threshold before threshold time is identified
threshDenominator = 2 ; % threshold = peak/threshDenominator

% load stimulus and get triggers for each part of concat data set (I don't think this
% cannot be done in concatinated set)

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


% load concatinated data
dataRun = load_data(DataBlock(DB).DgConcat) ;
dataRun = load_neurons(dataRun) ;

% prep cells and get trigger times
for dset=1:3 ;
    for sp=1:l_sp ; % for each spatial period
        for cntrst=1:l_cntrst ; % for each contrast
            for tp=1:l_tp ; % for each temporal period
                stim_group_triggers{dset}{sp,tp, cntrst}= dataRun.triggers(stim_group_triggers_i{dset}{sp,tp,cntrst}) ; % trigger times
                psth{dset}{tp}{cntrst}{sp} = nans(length(dataRun.spikes),ceil(min(diff(stim_group_triggers{dset}{sp,tp,cntrst}'))*dataRun.sampling_rate)) ;
            end
        end
    end

    for tp=1:l_tp ; % for each temporal period
        for cntrst=1:l_cntrst ; % for each contrast
            for cells = 1:length(dataRun.spikes) ; % for each cell
                for sp=1:l_sp ; % for each spatial period
                    [psthTemp,PsthTime] = get_smooth_psth(dataRun.spikes{cells},stim_group_triggers{dset}{sp,tp,cntrst}','stop',min(diff(stim_group_triggers{dset}{sp,tp,cntrst}'))) ;
                    psth{dset}{tp}{cntrst}{sp}(cells,:) = psthTemp ;
                    
                    mn =  min(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                    psth_min{dset}{tp}{cntrst}(cells,sp) = mn ;
                    
                    [mx,mi] =  max(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                    psth_peak{dset}{tp}{cntrst}(cells,sp) = mx ;
                    psth_peak_time{dset}{tp}{cntrst}(cells,sp) = PsthTime(mi) ;
%                     thresh = (mx-mn)/threshDenominator+mn ;
%                     i1 = find(psthTemp<thresh,1,'first') ; % first point below threshold
%                     i2 = find(psthTemp(i1+bufferTime*dataRun.sampling_rate:end)>thresh,1,'first') ; % first point above thresh after below thresh 
%                     psth_threshold_time{dset}{tp}{cntrst}(cells,sp) = PsthTime(i2+i1+bufferTime*dataRun.sampling_rate) ; % time of this point
%                      
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
                end
                
                psth_var{dset}{tp}{cntrst}(cells,sp) = nanvar(psth{dset}{tp}{cntrst}{sp}(cells,:)) ;
                [mx,mi] = max(psth_var{dset}{tp}{cntrst}(cells,:)) ;
                psth_var_norm{dset}{tp}{cntrst}(cells,:) = psth_var{dset}{tp}{cntrst}(cells,:)/mx ;
                psth_var_peak{tp}{cntrst}(cells,dset) = SpatialPeriod(mi) ;

                [mx,mi] = max(psth_var_f1{dset}{tp}{cntrst}(cells,:)) ;
                psth_var_f1_norm{dset}{tp}{cntrst}(cells,:) = psth_var_f1{dset}{tp}{cntrst}(cells,:)/mx ;
                psth_var_f1_peak{tp}{cntrst}(cells,dset) = SpatialPeriod(mi) ;

                [mx,mi] = max(psth_var_f2{dset}{tp}{cntrst}(cells,:)) ;
                psth_var_f2_norm{dset}{tp}{cntrst}(cells,:) = psth_var_f2{dset}{tp}{cntrst}(cells,:)/mx ;
                psth_var_f2_peak{tp}{cntrst}(cells,dset) = SpatialPeriod(mi) ;

            end

            psth_var_mean{dset}{tp}{cntrst} = nanmean(psth_var{dset}{tp}{cntrst},1) ; % nans can happen when no spikes are detected
            psth_var_norm_mean{dset}{tp}{cntrst} = nanmean(psth_var_norm{dset}{tp}{cntrst},1) ;

            psth_var_f1_mean{dset}{tp}{cntrst} = nanmean(psth_var_f1{dset}{tp}{cntrst},1) ;
            psth_var_f1_norm_mean{dset}{tp}{cntrst} = nanmean(psth_var_f1_norm{dset}{tp}{cntrst},1) ; 

            psth_var_f2_mean{dset}{tp}{cntrst} = nanmean(psth_var_f2{dset}{tp}{cntrst},1) ;
            psth_var_f2_norm_mean{dset}{tp}{cntrst} = nanmean(psth_var_f2_norm{dset}{tp}{cntrst},1) ;
        end
    end
end

% compare control, drug and wash conditions
% find psth var drug effect index
for cntrst=1:l_cntrst ; % for each contrast
    for cells = 1:length(dataRun.spikes) ; % for each cell
        TempNm = psth_var{2}{1}{cntrst}(cells,:) - (psth_var{1}{1}{cntrst}(cells,:)+psth_var{3}{1}{cntrst}(cells,:))/2 ; % drug-(cnt+wash)/2
        TempDnm = abs(psth_var{1}{1}{cntrst}(cells,:)-psth_var{3}{1}{cntrst}(cells,:)) ; % abs(cntrl-wash)
        psth_var_deltaDrug{cntrst}(cells,:)  = TempNm./(TempDnm+10) ;
    end
end

% find psth drug effect index
for cntrst=1:l_cntrst ; % for each contrast
    for sp=1:l_sp ; % for each spatial period
        for cells = 1:length(dataRun.spikes) ; % for each cell
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
    
% phase and amp change
for dset=1:3 ;
    for tp=1:l_tp ; % for each temporal period
        for sp=1:l_sp ; % for each spatial period
            for cntrst=1:l_cntrst ; % for each contrast
                for cells = 1:length(dataRun.spikes) ; % for each cell
                    psth_phase_f1_drugDelta{tp}{cntrst}{sp}(cells,dset) = psth_phase_f1{dset}{tp}{cntrst}(cells,sp)- psth_phase_f1{1}{tp}{cntrst}(cells,sp) ;
                    psth_phase_f2_drugDelta{tp}{cntrst}{sp}(cells,dset) = psth_phase_f2{dset}{tp}{cntrst}(cells,sp)- psth_phase_f2{1}{tp}{cntrst}(cells,sp) ;
                
                end
            end
        end
    end
end

lineWidthRange = 3 ;
figure % psths of each cell in control, +psem, wash
for cells = 1:length(dataRun.spikes) ; % for each cell
    for dset=1:3 ;
        for sp=1:l_sp ; % for each spatial period 
            for cntrst=1:l_cntrst ; % for each contrast
                subplot(l_sp,l_cntrst,l_cntrst*(sp-1)+cntrst) ;
                plot(PsthTime,psth{dset}{1}{cntrst}{sp}(cells,:),'color',Color_list{dset},'linewidth',lineWidthRange*psth_deltaWash{cntrst}{sp}(cells)+.1)
                axis tight
                hold on
            end
        end
    end
    title(num2str(cells))
    pause
    clf
end

figure % spatial tuning plots normalized with peak comparisons
for cntrst=1:l_cntrst ; % for each contrast
    for dset=1:3 ;
        subplot(2,l_cntrst,cntrst)
        plot(SpatialPeriod, psth_var_norm_mean{dset}{1}{cntrst},Color_list{dset},'linewidth',4) 
        hold on
    end
    subplot(2,l_cntrst,l_cntrst+cntrst)
    plot([1:3],psth_var_peak{1}{cntrst})
    hold on
    plot([1:3],mean(psth_var_peak{1}{cntrst}),'k','linewidth',4)
end

figure % spatial tuning plots normalized with peak comparisons
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



