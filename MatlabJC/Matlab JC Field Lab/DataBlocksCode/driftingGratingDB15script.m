% script to get MSE from DataBlock(15) which is missing triggers 
% this script assumes dataRun.triggers are the same in DataBlock(14).

% copied from KoSpatialTuningAnalyzerConcatPlusMappingV3 with some added
% parts as noted.

% JC 3/3/16

Color_list = {'k','r','b','g','c','v'} ; % order of colors for each 
TrialTrigInterval = 10 ; % sec
InterTrialTime = 2 ; % sec
StimTransTime = .1 ; % (sec) skip cycles before this time after starting new stimulus
MicronPerPix = 4 ; % (um/pix)
EsemFlag = false; % calc psth esemble
MseRelThreshold = 0.5 ; % (unitless) below this point response is "stable"
MseRelDrugDivWashThreshold = 1 ; % (unitless) above this point response has big drug effect 
LfSpThresh = 300 ; % (pix) spatial periods equal and above are termed "low frequency"

% ADDED- get triggers from previous data block (DB 14)
DB = 14; 
dset = 1 ;
dataRun = load_data(DataBlock(DB).DgPath{dset}) ;
dataRun = load_neurons(dataRun) ;
Triggers = dataRun.triggers ;
clear dataRun
dataRun = load_data(DataBlock(DB).DgConcat) ;
dataRun = load_neurons(dataRun) ;
TriggersConcat = dataRun.triggers ;
clear dataRun

DB = 15 ;

% load stimulus and get triggers for each part of concat data set
for dset=1:3 ;
    % load data 
    dataRun = load_data(DataBlock(DB).DgPath{dset}) ;
    dataRun = load_neurons(dataRun) ;
    
    dataRun.triggers = Triggers ; % ADDED

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

dataRun.triggers = TriggersConcat ; % ADDED

l_cells = length(dataRun.spikes) ; % number of cells

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
                    [powerspec_xvalues, mean_powerspec] = PowerSpectrumFinder(psth{dset}{tp}{cntrst}{sp}(cells,:),dataRun.sampling_rate) ;
                    psth_varF1{dset}{tp}{cntrst}(cells,sp) = mean_powerspec(2) ;
                    psth_varF2{dset}{tp}{cntrst}(cells,sp) = mean_powerspec(3) ;
                    
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
                
                % mean across sp
                psth_varF1_mean{dset}{tp}{cntrst}(cells) = mean(psth_varF1{dset}{tp}{cntrst}(cells,:)) ;
                psth_varF2_mean{dset}{tp}{cntrst}(cells) = mean(psth_varF2{dset}{tp}{cntrst}(cells,:)) ;
                psth_varF2F1ratio{dset}{tp}{cntrst}(cells) = psth_varF2_mean{dset}{tp}{cntrst}(cells)/psth_varF1_mean{dset}{tp}{cntrst}(cells) ;
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
MseRelDrug_DivMseRelWash_Hist_cumsum = cumsum(MseRelDrug_DivMseRelWash_Hist)/sum(MseRelDrug_DivMseRelWash_Hist) ; 

MseRelDrug_DivMseRelWashLfMean_Hist = hist(MseRelDrug_DivMseRelWashLfMean,MseRelHistX) ;
MseRelDrug_DivMseRelWashLfMean_FractAboveThresh = sum(MseRelDrug_DivMseRelWashLfMean>MseRelDrugDivWashThreshold)/length(MseRelDrug_DivMseRelWashLfMean) ;
MseRelDrug_DivMseRelWashLfMean_Hist_cumsum = cumsum(MseRelDrug_DivMseRelWashLfMean_Hist)/sum(MseRelDrug_DivMseRelWashLfMean_Hist) ;   

   
numResponsesAbove = sum(tempMat(:)>MseRelDrugDivWashThreshold) ;
MseRelDrug_DivMseRelWash_FractAboveThreshHist = zeros(l_cntrst, l_sp) ;
for cells=1:l_cells ;
    MseRelDrug_DivMseRelWash_FractAboveThreshHist = MseRelDrug_DivMseRelWash_FractAboveThreshHist + (MseRelDrug_DivMseRelWash{cells}>MseRelDrugDivWashThreshold)/numResponsesAbove ;
end

MseRelWashLfMean_Hist = hist(MseRelWashLfMean,MseRelHistX) ;
MseRelWashLfMean_Hist_cumsum = cumsum(MseRelWashLfMean_Hist)/sum(MseRelWashLfMean_Hist) ;

MseRelDrugLfMean_Hist = hist(MseRelDrugLfMean,MseRelHistX) ;
MseRelDrugLfMean_Hist_cumsum = cumsum(MseRelDrugLfMean_Hist)/sum(MseRelDrugLfMean_Hist) ;

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
   
    nxtval = input('next (cell num, 0=back, return=forward)') ;
    if isempty(nxtval) ;
        cells=cells+1 ;
    elseif nxtval == 0 ;
        cells=cells-1 ;
    elseif nxtval>0 ;
        cells=nxtval ;
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

clear dataRunMaster dataRun
save(['/Users/jcafaro/Desktop/Temp/Matfiles/KoSpatialV3',num2str(DB)],'-v7.3')

exportStructToHDF5(ForIgor,['/Users/jcafaro/Desktop/Temp/Matfiles/KoSpatialV3Db',num2str(DB),'.h5'],'/');



