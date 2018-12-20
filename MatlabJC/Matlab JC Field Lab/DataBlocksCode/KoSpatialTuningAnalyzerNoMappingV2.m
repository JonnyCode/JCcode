function ForIgor = KoSpatialTuningAnalyzerNoMappingV2(DataBlock, DB, Params)

% Adapted from "KoSpatialTuningAnalyzerNoMapping" to include parsing of
% data with diffferent contrast values

% JC 12/9/15
%DB=9; % TEMP

% parameters
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
TrialTrigInterval = 10 ; % sec
InterTrialTime = 2 ; % sec
StimTransTime = .1 ; % (sec) skip cycles before this time after starting new stimulus

numDs = length(DataBlock(DB).DgPath) ;
 
for ds = 1:numDs % for each data set

    % load data
    dataRun = load_data(DataBlock(DB).DgPath{ds}) ;
    dataRun = load_neurons(dataRun) ;
        
    % load stimulus
    dataRun.names.stimulus_path = [DataBlock(DB).DgPath{ds}(1:end-15),'stimuli/s',DataBlock(DB).DgPath{ds}(end-1:end)] ;
    dataRun = load_stim(dataRun,'user_defined_trigger_interval', TrialTrigInterval) ;

    dataRun.stimulus.params.RGB = unique(dataRun.stimulus.params.RGB) ; % not sure why this is not already done as for TP and SP
    
    % stim trials with the same spatial and temporal periods
    for sp=1:length(dataRun.stimulus.params.SPATIAL_PERIOD) ; % for each spatial period
        for tp=1:length(dataRun.stimulus.params.TEMPORAL_PERIOD) ; % for each temporal period
            for cntrst=1:length(dataRun.stimulus.params.RGB) ; % for each contrast
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

    % trigger times for each
    for sp=1:length(dataRun.stimulus.params.SPATIAL_PERIOD) ; % for each spatial period
        for tp=1:length(dataRun.stimulus.params.TEMPORAL_PERIOD) ; % for each temporal period
            for cntrst=1:length(dataRun.stimulus.params.RGB) ; % for each contrast
                stim_group_triggers{sp,tp, cntrst} = [] ;
                tempTrigs = dataRun.stimulus.triggers(stim_groups_trials{sp,tp,cntrst}) ; % the first trigger for each stimulus 
                for a = 1:length(tempTrigs) ; % for each first trigger
                    stim_group_triggers{sp,tp,cntrst} = [stim_group_triggers{sp,tp,cntrst}(:); dataRun.triggers(dataRun.triggers>=tempTrigs(a)+StimTransTime & ...
                        dataRun.triggers <= tempTrigs(a)+TrialTrigInterval-InterTrialTime)] ; % find the triggers within that set
                end
            end
        end
    end
    
    % prep cells
    for sp=1:length(dataRun.stimulus.params.SPATIAL_PERIOD) ; % for each spatial period
        for cntrst=1:length(dataRun.stimulus.params.RGB) ; % for each contrast
            for tp=1:length(dataRun.stimulus.params.TEMPORAL_PERIOD) ; % for each temporal period
                psth{tp}{cntrst}{sp} = nans(length(dataRun.spikes),ceil(min(diff(stim_group_triggers{sp,tp,cntrst}'))*dataRun.sampling_rate)) ;
            end
        end
    end
    
    for tp=1:length(dataRun.stimulus.params.TEMPORAL_PERIOD) ; % for each temporal period
        for cntrst=1:length(dataRun.stimulus.params.RGB) ; % for each contrast
            for cells = 1:length(dataRun.spikes) ; % for each cell
                for sp=1:length(dataRun.stimulus.params.SPATIAL_PERIOD) ; % for each spatial period
                    psth{tp}{cntrst}{sp}(cells,:) = get_smooth_psth(dataRun.spikes{cells},stim_group_triggers{sp,tp,cntrst}','stop',min(diff(stim_group_triggers{sp,tp,cntrst}'))) ;
                    psth_var{tp}{cntrst}(cells,sp) = nanvar(psth{tp}{cntrst}{sp}(cells,:)) ;

                    [powerspec_xvalues, mean_powerspec] = PowerSpectrumFinder(psth{tp}{cntrst}{sp}(cells,:),dataRun.sampling_rate) ;
                    if sum(mean_powerspec(2:3)) ~=0 ; % if there is any f1 or f2 power 
                        psth_var_f1{tp}{cntrst}(cells,sp) = mean_powerspec(2)/sum(mean_powerspec(2:3)) ; % f1/(f1+f2) power
                        psth_var_f2{tp}{cntrst}(cells,sp) = mean_powerspec(3)/sum(mean_powerspec(2:3)) ; % f2/(f1+f2) power
                    else
                        psth_var_f1{tp}{cntrst}(cells,sp) = 0 ; % f1 power
                        psth_var_f2{tp}{cntrst}(cells,sp) = 0 ; % f2 power
                    end
                end
                [mx,mi] = max(psth_var{tp}{cntrst}(cells,:)) ;
                psth_var_norm{tp}{cntrst}(cells,:) = psth_var{tp}{cntrst}(cells,:)/mx ;
                psth_var_peak{tp}{cntrst}(cells) = mi ;

                [mx,mi] = max(psth_var_f1{tp}{cntrst}(cells,:)) ;
                psth_var_f1_norm{tp}{cntrst}(cells,:) = psth_var_f1{tp}{cntrst}(cells,:)/mx ;
                psth_var_f1_peak{tp}{cntrst}(cells) = mi ;

                [mx,mi] = max(psth_var_f2{tp}{cntrst}(cells,:)) ;
                psth_var_f2_norm{tp}{cntrst}(cells,:) = psth_var_f2{tp}{cntrst}(cells,:)/mx ;
                psth_var_f2_peak{tp}{cntrst}(cells) = mi ;

            end
            
            psth_var_mean{tp}{cntrst} = mean(psth_var{tp}{cntrst},1) ;
            psth_var_norm_mean{tp}{cntrst} = mean(psth_var_norm{tp}{cntrst},1) ;

            psth_var_f1_mean{tp}{cntrst} = mean(psth_var_f1{tp}{cntrst},1) ;
            psth_var_f1_norm_mean{tp}{cntrst} = mean(psth_var_f1_norm{tp}{cntrst},1) ; 

            psth_var_f2_mean{tp}{cntrst} = mean(psth_var_f2{tp}{cntrst},1) ;
            psth_var_f2_norm_mean{tp}{cntrst} = mean(psth_var_f2_norm{tp}{cntrst},1) ;
        end
    end

    figure
    for cells = 1:length(dataRun.spikes) ; % for each cell
        for sp=1:length(dataRun.stimulus.params.SPATIAL_PERIOD) ; % for each spatial period 
            for cntrst=1:length(dataRun.stimulus.params.RGB) ; % for each contrast
                subplot(length(dataRun.stimulus.params.SPATIAL_PERIOD),1,sp) ;
                plot([1:length(psth{1}{1}{1}(1,:))]/20000,psth{1}{cntrst}{sp}(cells,:))
                axis tight
                hold on
            end
        end
        pause
        clf
    end
    
    figure
    for cntrst=1:length(dataRun.stimulus.params.RGB) ; % for each contrast
        subplot(length(dataRun.stimulus.params.RGB),1,cntrst)
        plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var{1}{cntrst})
        hold on
        plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_mean{1}{cntrst},'k*-','linewidth',4)
    end
    
    figure
    for cntrst=1:length(dataRun.stimulus.params.RGB) ; % for each contrast
        subplot(length(dataRun.stimulus.params.RGB),1,cntrst)
        plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_norm{1}{cntrst})
        hold on
        plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_norm_mean{1}{cntrst},'k*-','linewidth',4)
    end
  
    figure
    for cntrst=1:length(dataRun.stimulus.params.RGB) ; % for each contrast
        subplot(length(dataRun.stimulus.params.RGB),1,cntrst)
        plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_f1_norm{1}{cntrst})
        hold on
        plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_f1_norm_mean{1}{cntrst},'k*-','linewidth',4)
    end
    
    figure
    for cntrst=1:length(dataRun.stimulus.params.RGB) ; % for each contrast
        subplot(length(dataRun.stimulus.params.RGB),1,cntrst)
        plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_f2_norm{1}{cntrst})
        hold on
        plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_f2_norm_mean{1}{cntrst},'k*-','linewidth',4)
    end
    
    figure
    for sp=1:length(dataRun.stimulus.params.SPATIAL_PERIOD) ; % for each spatial period 
        for cntrst=1:length(dataRun.stimulus.params.RGB) ; % for each contrast
            subplot(length(dataRun.stimulus.params.SPATIAL_PERIOD),length(dataRun.stimulus.params.RGB),(sp-1)*length(dataRun.stimulus.params.RGB)+cntrst) ;
            i = find(psth_var_peak{tp}{cntrst} == sp) ;
            if length(i)>0 ;
                plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_norm{1}{cntrst}(i,:)) ;
                hold on
                plot(dataRun.stimulus.params.SPATIAL_PERIOD, mean(psth_var_norm{1}{cntrst}(i,:),1),'k*-','linewidth', 4) ;
                title(num2str(length(i)))
            end    
        end
    end

    %set(gca,'XScale','log')    
end


        

