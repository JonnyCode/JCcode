function ForIgor = KoSpatialTuningAnalyzerNoMapping(DataBlock, DB, Params)

% this function will analyze drfiting grating responses from array data
% before and after cells are knocked out (KO)

% JC 11/13/15 % 
%DB=9; % TEMP

% parameters
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
TrialTrigInterval = 10 ; % sec
InterTrialTime = 2 ; % sec
StimTransTime = 1 ; % sec

numDs = length(DataBlock(DB).DgPath) ;
 
for ds = 1:numDs % for each data set

    % load data
    dataRun = load_data(DataBlock(DB).DgPath{ds}) ;
    dataRun = load_neurons(dataRun) ;
        
    % load stimulus
    dataRun.names.stimulus_path = [DataBlock(DB).DgPath{1}(1:end-15),'stimuli/s',DataBlock(DB).DgPath{1}(end-1:end)] ;
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
    
    for sp=1:length(dataRun.stimulus.params.SPATIAL_PERIOD) ; % for each spatial period
        for tp=1:length(dataRun.stimulus.params.TEMPORAL_PERIOD) ; % for each temporal period
            psth{tp}{sp} = nans(length(dataRun.spikes),ceil(min(diff(stim_group_triggers{sp,tp}'))*dataRun.sampling_rate)) ; 
        end
    end
    
    for tp=1:length(dataRun.stimulus.params.TEMPORAL_PERIOD) ; % for each temporal period
        for cells = 1:length(dataRun.spikes) ; % for each cell
            for sp=1:length(dataRun.stimulus.params.SPATIAL_PERIOD) ; % for each spatial period
                psth{tp}{sp}(cells,:) = get_smooth_psth(dataRun.spikes{cells},stim_group_triggers{sp,tp}','stop',min(diff(stim_group_triggers{sp,tp}'))) ;
                psth_var{tp}(cells,sp) = nanvar(psth{tp}{sp}(cells,:)) ;
               
                [powerspec_xvalues, mean_powerspec] = PowerSpectrumFinder(psth{tp}{sp}(cells,:),dataRun.sampling_rate) ;
                if sum(mean_powerspec(2:3)) ~=0 ; % if there is any f1 or f2 power 
                    psth_var_f1{tp}(cells,sp) = mean_powerspec(2)/sum(mean_powerspec(2:3)) ; % f1/(f1+f2) power
                    psth_var_f2{tp}(cells,sp) = mean_powerspec(3)/sum(mean_powerspec(2:3)) ; % f2/(f1+f2) power
                else
                    psth_var_f1{tp}(cells,sp) = 0 ; % f1 power
                    psth_var_f2{tp}(cells,sp) = 0 ; % f2 power
                end
            end
            [mx,mi] = max(psth_var{tp}(cells,:)) ;
            psth_var_norm{tp}(cells,:) = psth_var{tp}(cells,:)/mx ;
            psth_var_peak{tp}(cells) = mi ;
            
            [mx,mi] = max(psth_var_f1{tp}(cells,:)) ;
            psth_var_f1_norm{tp}(cells,:) = psth_var_f1{tp}(cells,:)/mx ;
            psth_var_f1_peak{tp}(cells) = mi ;
            
            [mx,mi] = max(psth_var_f2{tp}(cells,:)) ;
            psth_var_f2_norm{tp}(cells,:) = psth_var_f2{tp}(cells,:)/mx ;
            psth_var_f2_peak{tp}(cells) = mi ;
            
        end
        psth_var_mean{tp} = mean(psth_var{tp},1) ;
        psth_var_norm_mean{tp} = mean(psth_var_norm{tp},1) ;
        
        psth_var_f1_mean{tp} = mean(psth_var_f1{tp},1) ;
        psth_var_f1_norm_mean{tp} = mean(psth_var_f1_norm{tp},1) ; 
        
        psth_var_f2_mean{tp} = mean(psth_var_f2{tp},1) ;
        psth_var_f2_norm_mean{tp} = mean(psth_var_f2_norm{tp},1) ;
    end

    figure
    plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_norm{1})
    hold on
    plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_norm_mean{1},'k*-','linewidth',4)
  
    figure
    plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_f1_norm{1})
    hold on
    plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_f1_norm_mean{1},'k*-','linewidth',4)
    
    figure
    plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_f2_norm{1})
    hold on
    plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_f2_norm_mean{1},'k*-','linewidth',4)
    
    figure
    for sp=1:length(dataRun.stimulus.params.SPATIAL_PERIOD) ; % for each spatial period 
        subplot(length(dataRun.stimulus.params.SPATIAL_PERIOD),1,sp) ;
        i = find(psth_var_peak{tp} == sp) ;
        if length(i)>0 ;
            plot(dataRun.stimulus.params.SPATIAL_PERIOD, psth_var_norm{1}(i,:)) ;
            hold on
            plot(dataRun.stimulus.params.SPATIAL_PERIOD, mean(psth_var_norm{1}(i,:),1),'k*-','linewidth', 4) ;
            title(num2str(length(i)))
        end    
    end
        
    
    
    set(gca,'XScale','log')
    
end


        

