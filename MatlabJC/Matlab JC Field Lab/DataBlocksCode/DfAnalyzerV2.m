function ForIgor = DfAnalyzerV2(DataBlock, DB, Params)

% modified from 'DfAnalyzer' to include data using old LED stim and
% different analyses

% JC 7/8/15

% CONDITION ORGANIZING BUG!!!!!!!!!!!!!!

% parameters
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
interFlashIntVar = .005 ; % (sec) expected precision of trigger intervals
minNumFlashes = 5 ; % min number of flashes in a block 
LED_settings_calibrations_path = '/Volumes/lab/Experiments/Calibration/LED_settings_calibrations' ; % location of LED calibration mat file
Snr_threshold = 3 ;

% compare spatial and temporal receptive fields across data defined sets 

% load data
dataRun = load_data(DataBlock(DB).DfPath) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_params(dataRun) ;

% params
NumCells = length(dataRun.spikes) ;

% estimate the relative stim intensity assuming flashes are within linear response range
if isfield(DataBlock(DB).DfParams,'Setting') ;
    load(LED_settings_calibrations_path) 
    for a=1:length(DataBlock(DB).DfParams.Ftime) ;
        ci = find(LED_settings_calibrations.FlashDuration==DataBlock(DB).DfParams.Ftime(a)) ;
        ri = find(LED_settings_calibrations.setting==DataBlock(DB).DfParams.Setting(a)) ;
        Irel(a) = LED_settings_calibrations.NormPower(ri,ci)*DataBlock(DB).DfParams.Ftime(a)./10.^DataBlock(DB).DfParams.NDF(a) ;
    end
else
    Irel = DataBlock(DB).DfParams.Ftime./10.^DataBlock(DB).DfParams.NDF ;
end
Irel_unique = unique(Irel) ;

% find sets of flash stimuli
st=1 ; 
trigger_set_i{1} = [] ;
for a=1:length(dataRun.triggers) ; % for each trigger
    trigger_set_i{st} = [trigger_set_i{st},a] ; % put it in a set   
    if a<length(dataRun.triggers) ; % if its not the last trigger
        if sum(abs(dataRun.triggers(a+1)-dataRun.triggers(a)-DataBlock(DB).DfParams.interFlashInt)<interFlashIntVar)>0 ; % next trigger is one of accepted intervals
            st = st ; % keep it in the set
        else
            st = st+1 ; % put it in a new set
            trigger_set_i{st} = [] ;
        end
    end
end

% arrange and concatinate trigger sets that have min number of flashes
Full_trigger_set_i = cell(1,length(Irel_unique)) ; % prep array
ds = 1 ; % data set
for ts=1:length(trigger_set_i) ; % for each trigger set  
    if length(trigger_set_i{ts})>minNumFlashes ; % if more than X flashes     
        [v,i] = intersect(Irel_unique,Irel(ds)) ; % index of unique intensity
        Full_trigger_set_i{i} = [Full_trigger_set_i{i},trigger_set_i{ds}] ;              
        
        NDF(i) = DataBlock(DB).DfParams.NDF(ds) ;
        Ftime(i) = DataBlock(DB).DfParams.Ftime(ds) ;
        if isfield(DataBlock(DB).DfParams,'Setting') ;
            Setting(i) = DataBlock(DB).DfParams.Setting(ds) ;
        end
        ds = ds+1 ;
    end   
end

if ds-1 ~= length(DataBlock(DB).DfParams.Ftime) ; % if the number of identified flash blocks is not as expected
    error('Params or data length mismatch')
end

% get psth in each cell for each set of stimuli 
for ds=1:length(Full_trigger_set_i) ; % for each trigger set  
    for c = 1:NumCells ; % for each cell
        [psthTemp,psthTime] = get_smooth_psth(dataRun.spikes{c}, dataRun.triggers(Full_trigger_set_i{ds}), 'start', -2, 'stop', DataBlock(DB).DfParams.interFlashInt) ; 
        psth{ds}(c,:) = psthTemp ;

        preStimPnts = [1:find(psthTime==0)] ; % points before flash
        preStimMean = mean(psthTemp(preStimPnts)) ; % mean before flash
        preStimStd = std(psthTemp(preStimPnts)) ; % std before flash

        psth_delta{ds}(c,:) = psthTemp - preStimMean ; % response change
        [v,i] = max(abs(psth_delta{ds}(c,:))) ; % largest response change
        psth_delta_peak(ds,c) = psth_delta{ds}(c,i) ;

        psth_delta_peak_abs(ds,c) = abs(psth_delta_peak(ds,c)) ; % abs value of change
    end
    psth_mean(ds,:) = mean(psth{ds}) ; % average across all cells
    psth_std(ds,:) = std(psth{ds}) ; 
end

% fit hockey stick - find threshold and gain for each cell
for c = 1:NumCells ; % for each cell
    Temp = nlinfit(Irel_unique,psth_delta_peak(:,c)',@HockeyStickFree2,[10^-4,1]) ;
    Theshold(c) = Temp(1) ;
    Gain(c) = Temp(2) ;
    Hs_fit(c,:) = HockeyStickFree2([Theshold(c),Gain(c)],Irel_unique) ;
end

% histograms
HistBins = [-200:5:200] ;
for ds=1:length(Full_trigger_set_i) ;
    HistValues= hist(psth_delta_peak(ds,:),HistBins) ; 
    pdf(ds,:) = HistValues/sum(HistValues) ;
end
    
% figures
Ftimes = unique(DataBlock(DB).DfParams.Ftime) ;
NDFs = unique(DataBlock(DB).DfParams.NDF) ;

for ds=1:length(Full_trigger_set_i) ; % for each trigger set 
    plot_i(ds) = (find(NDFs == NDF(ds))-1)*length(Ftimes) + find(Ftimes == Ftime(ds)) ; % subplot i
end

figure
for ds = 1:size(psth_mean,1) ;
    subplot(length(NDFs),length(Ftimes),plot_i(ds))
    plot(psthTime,psth_mean(ds,:),'k','linewidth',3) ;
    hold on
    %plot(psthTime,psth{ds},'c') ; % all cells
    plot(psthTime,psth_mean(ds,:)+psth_std(ds,:),'b:') ;
    plot(psthTime,psth_mean(ds,:)-psth_std(ds,:),'b:') ;
    text(.9,.9,['NDF',num2str(NDF(ds)),',fdur ',num2str(Ftime(ds))],'units', 'norm') 
end

figure
for ds = 1:size(psth_mean,1) ;
    subplot(length(NDFs),length(Ftimes),plot_i(ds))
    plot(HistBins,pdf(ds,:)) ;
    text(.9,.9,['NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),',fdur ',num2str(DataBlock(DB).DfParams.Ftime(ds))],'units', 'norm') 
end

figure
plot(log(Irel_unique),psth_delta_peak)
hold on
plot(log(Irel_unique),psth_delta_peak_Irel_mean,'k','linewidth',2)
plot(log(Irel_unique),psth_delta_peak_Irel_meanOn,'k:','linewidth',2)
plot(log(Irel_unique),psth_delta_peak_Irel_meanOff,'k:','linewidth',2)

figure
plot(log(Irel_unique),psth_absRspSum_Irel)
hold on
plot(log(Irel_unique),psth_absRspSum_Irel_mean,'k','linewidth',2)

% for igor
ForIgor = struct() ;

VecName = ['psthTime','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psthTime) ; % psth time vector

VecName = ['histBins','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,HistBins) ; % pdf bins

for ds = 1:size(psth_mean,1) ;
    VecName = ['psth','NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),'Fdur',num2str(DataBlock(DB).DfParams.Ftime(ds)),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor, VecName, psth_mean(ds,:)) ; % psth vector
    
    VecName = ['pdf','NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),'Fdur',num2str(DataBlock(DB).DfParams.Ftime(ds)),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor, VecName, pdf(ds,:)) ; % pdf vector
end

% avearge intensity respones curves
VecName = ['Irel','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,Irel_unique) ; % relative intensity vector

VecName = ['peak_mean','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psth_delta_peak_Irel_mean) ; % relative intensity vector

VecName = ['peak_mean_On','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psth_delta_peak_Irel_meanOn) ; % relative intensity vector

VecName = ['peak_mean_Off','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psth_delta_peak_Irel_meanOff) ; % relative intensity vector

    




