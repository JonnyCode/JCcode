function ForIgor = DfAnalyzerV4(DataBlock, DB, Params)

% adapted from DfAnalyzer3; this code selects ON cells only 

% JC 9/23/15

% parameters
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
interFlashIntVar = .005 ; % (sec) expected precision of trigger intervals
minNumFlashes = 5 ; % min number of flashes in a block 
psth_bin_size = 0.1 ; % (sec) 
StdAbove = 2 ; % number of std above mean
LED_settings_calibrations_path = '/Volumes/lab/Experiments/Calibration/LED_settings_calibrations' ; % location of LED calibration mat file
LoadCalibrationFlag = true ;
Off_frac = .10 ; % fraction of stim that show a "off" response needed to call a cell an OFF
PdfHistFlag = false ; % make histograms in pdfs if true

% compare spatial and temporal receptive fields across data defined sets 

% load data
dataRun = load_data(DataBlock(DB).DfPath) ;
dataRun = load_neurons(dataRun) ;
dataRun = load_params(dataRun) ;

% params
NumCells = length(dataRun.spikes) ;

% find sets of flash stimuli
st=1 ; 
trigger_set_i{1} = [] ;
for a=1:length(dataRun.triggers) ; % for each trigger
    trigger_set_i{st} = [trigger_set_i{st},a] ; % put it in a set   
    if a<length(dataRun.triggers) ; % if its not the last trigger
        if sum(abs(dataRun.triggers(a+1)-dataRun.triggers(a)-DataBlock(DB).DfParams.interFlashInt)<interFlashIntVar)>0 ; % next trigger is has the right interval
            st = st ; % keep it in the set
        else
            st = st+1 ; % put it in a new set
            trigger_set_i{st} = [] ;
        end
    end
end

% id OFF and ON cells 
Offp_i = [] ;
ds = 0 ; % data set
for ts=1:length(trigger_set_i) ; % for each trigger set  
    if length(trigger_set_i{ts})>minNumFlashes ; % if more than X flashes
        ds = ds+1 ;
        for c = 1:NumCells ; % for each cell
            [psthTemp,psthTime] = get_smooth_psth(dataRun.spikes{c}, dataRun.triggers(trigger_set_i{ts}), 'start', -2, 'stop', DataBlock(DB).DfParams.interFlashInt, 'bin_size', psth_bin_size) ; 
            
            StimPnt = find(psthTime==0) ; % flash point
            preStimPnts = [1:StimPnt] ; % points before flash
            postStimPnts = [StimPnt+1:length(psthTime)] ; % points after flash
            
            preStimMean = mean(psthTemp(preStimPnts)) ; % mean before flash
            preStimStd = std(psthTemp(preStimPnts)) ; % std before flash
            preStimRange = range(psthTemp(preStimPnts)) ; % max before flash
        
            psth_deltaTemp = psthTemp - preStimMean ; % response change
            
            [mxv,mxi] = max(psth_deltaTemp(postStimPnts)) ; % max response after stim
            [mnv,mni] = min(psth_deltaTemp(postStimPnts)) ; % min response after stim
            
            if mnv<min(psth_deltaTemp(preStimPnts))  && mni<mxi ; % if the min value is less than the min prestim value and precedes the time of the max value
                Offp_i = [Offp_i,c] ;
            end
        end
    end
end
uOffp_i = unique(Offp_i) ;
Off_i = [] ;
for c=1:length(uOffp_i) ; % for each putative off cell
    if sum(uOffp_i(c)==Offp_i) > ds*Off_frac ; % if a large fraction of stim show the off response  
        Off_i = [Off_i,uOffp_i(c)] ;
    end
end

On_i = setxor([1:NumCells],unique(Off_i)) ; % on cell i
NumOnCells = length(On_i) ;

% get psth in each cell for each set of stimuli 
ds = 1 ; % data set
for ts=1:length(trigger_set_i) ; % for each trigger set  
    if length(trigger_set_i{ts})>minNumFlashes ; % if more than X flashes
        for c = 1:NumOnCells ; % for each cell
            [psthTemp,psthTime] = get_smooth_psth(dataRun.spikes{On_i(c)}, dataRun.triggers(trigger_set_i{ts}), 'start', -2, 'stop', DataBlock(DB).DfParams.interFlashInt, 'bin_size', psth_bin_size) ; 
            psth{ds}(c,:) = psthTemp ;
            
            StimPnt = find(psthTime==0) ; % flash point
            preStimPnts = [1:StimPnt] ; % points before flash
            postStimPnts = [StimPnt+1:length(psthTime)] ; % points after flash
            
            preStimMean = mean(psthTemp(preStimPnts)) ; % mean before flash
            preStimStd = std(psthTemp(preStimPnts)) ; % std before flash
            preStimRange = range(psthTemp(preStimPnts)) ; % max before flash
        
            psth_delta{ds}(c,:) = psthTemp - preStimMean ; % response change
            [v,i] = max(abs(psth_delta{ds}(c,postStimPnts))) ; % largest response change after flash
            psth_delta_peak(ds,c) = psth_delta{ds}(c,postStimPnts(i)) ;
            
            psth_delta_peak_abs(ds,c) = abs(psth_delta_peak(ds,c)) ;
            psth_delta_peak_abs_Snr(ds,c) = psth_delta_peak_abs(ds,c)/preStimStd ;
            
            psth_delta_peak_time(ds,c) = psthTime(postStimPnts(i)) ; % time of peak after stim onset
             
            i1 = find(abs(psth_delta{ds}(c,1:postStimPnts(i)))<(v/2),1,'last') ; % first point
            i2 = find(abs(psth_delta{ds}(c,postStimPnts(i):end))<v/2,1,'first')+postStimPnts(i) ; % last point 
            if ~isempty(i1 + i2) ;
                psth_delta_halfWidth(ds,c) = psthTime(i2) - psthTime(i1) ; % half width
            else
                psth_delta_halfWidth(ds,c) = nan ; % no half width can be found
            end
            
            for trls = 1:length(trigger_set_i{ts}) ;
                trigTime = dataRun.triggers(trigger_set_i{ts}(trls))+psthTime(postStimPnts(i)) ;
                
                psth_peak_trls{ds}(c,trls) = sum(((trigTime-psth_bin_size/2)<dataRun.spikes{On_i(c)}).*(dataRun.spikes{On_i(c)}<(trigTime+psth_bin_size/2)))/psth_bin_size ;    
                Bg_rate_trls{ds}(c,trls) = sum(((dataRun.triggers(trigger_set_i{ts}(trls))-1)<dataRun.spikes{On_i(c)}).*(dataRun.spikes{On_i(c)}<dataRun.triggers(trigger_set_i{ts}(trls)))) ;
                psth_peak_delta_trls{ds}(c,trls) = psth_peak_trls{ds}(c,trls)- Bg_rate_trls{ds}(c,trls) ;
            end
                     
            psth_peak_trls_mean(ds,c) = mean(psth_peak_trls{ds}(c,:),2) ;
            Bg_rate_trls_mean(ds,c) = mean(Bg_rate_trls{ds}(c,:),2) ;
            
            psth_peak_delta_trls_mean(ds,c) = mean(psth_peak_delta_trls{ds}(c,:),2) ;
            psth_peak_delta_trls_std(ds,c) = std(psth_peak_delta_trls{ds}(c,:),[],2) ;
           
        end
        psth_mean(ds,:) = mean(psth{ds}) ; % average across all cells
        psth_std(ds,:) = std(psth{ds}) ;
        
        psth_delta_mean(ds,:) = mean(psth_delta{ds}) ; 
        
        ds = ds+1 ;
    end   
end

if size(psth_mean,1) ~= length(DataBlock(DB).DfParams.Ftime) ; % if the number of identified flash blocks is not as expected
    error('Params or data length mismatch')
end

% estimate the relative stim intensity assuming flashes are within linear response range
if isfield(DataBlock(DB).DfParams,'Setting') ;
    load(LED_settings_calibrations_path) 
    for a=1:length(DataBlock(DB).DfParams.Ftime) ;
        ci = find(LED_settings_calibrations.FlashDuration==DataBlock(DB).DfParams.Ftime(a)) ;
        ri = find(round(LED_settings_calibrations.setting*10)==round(DataBlock(DB).DfParams.Setting(a)*10)) ;
        Irel(a) = LED_settings_calibrations.NormPower(ri,ci)*DataBlock(DB).DfParams.Ftime(a)./10.^DataBlock(DB).DfParams.NDF(a) ;
    end
elseif LoadCalibrationFlag ;
    load /Volumes/lab/Experiments/Calibration/NdfCalibration
    Irel = (DataBlock(DB).DfParams.Ftime/1000).*NdfCalibration(2,DataBlock(DB).DfParams.NDF+1) ;
else
    Irel = DataBlock(DB).DfParams.Ftime./10.^DataBlock(DB).DfParams.NDF ;
end
Irel_unique = unique(Irel) ;


% Organise by and average across cells with the same intensity levels 
psth_delta_peak_Irel = zeros(length(Irel_unique), NumOnCells) ;
psth_delta_peak_abs_Irel = zeros(length(Irel_unique), NumOnCells) ;
psth_delta_peak_abs_Snr_Irel = zeros(length(Irel_unique), NumOnCells) ;
psth_peak_delta_trls_mean_Irel = zeros(length(Irel_unique), NumOnCells) ;
psth_peak_delta_trls_std_Irel = zeros(length(Irel_unique), NumOnCells) ;

psth_delta_peak_time_Irel = zeros(length(Irel_unique), NumOnCells) ;
psth_delta_halfWidth_Irel = zeros(length(Irel_unique), NumOnCells) ;

for ds = 1:size(psth_mean,1) ;
    [v,i] = intersect(Irel_unique,Irel(ds)) ; % index of unique intensity
    NumSame = sum(Irel_unique(i)==Irel) ; % number of same stim
    psth_delta_peak_Irel(i,:) = psth_delta_peak_Irel(i,:)+ psth_delta_peak(ds,:)/NumSame ;
    psth_delta_peak_abs_Irel(i,:) = psth_delta_peak_abs_Irel(i,:)+ psth_delta_peak_abs(ds,:)/NumSame ;
    psth_delta_peak_abs_Snr_Irel(i,:) = psth_delta_peak_abs_Snr_Irel(i,:)+ psth_delta_peak_abs_Snr(ds,:)/NumSame ;
    
    psth_peak_delta_trls_mean_Irel(i,:) = psth_peak_delta_trls_mean_Irel(i,:)+ psth_peak_delta_trls_mean(ds,:)/NumSame ;
    psth_peak_delta_trls_std_Irel(i,:) = psth_peak_delta_trls_std_Irel(i,:)+ psth_peak_delta_trls_std(ds,:)/NumSame ;
    
    psth_delta_peak_time_Irel(i,:) = psth_delta_peak_time_Irel(i,:)+ psth_delta_peak_time(ds,:)/NumSame ;
    psth_delta_halfWidth_Irel(i,:) = psth_delta_halfWidth_Irel(i,:)+ psth_delta_halfWidth(ds,:)/NumSame ;
end
psth_delta_peak_Irel_mean = mean(psth_delta_peak_Irel,2) ;
psth_delta_peak_Irel_meanOn = mean(psth_delta_peak_Irel(:,sum(psth_delta_peak_Irel,1)>0),2) ;
psth_delta_peak_Irel_meanOff = mean(psth_delta_peak_Irel(:,sum(psth_delta_peak_Irel,1)<0),2) ;

% Thresholds
for c = 1:NumOnCells ; % for each cell
    ti = find(abs(psth_peak_delta_trls_mean_Irel(:,c))>(psth_peak_delta_trls_std_Irel(:,c)*StdAbove),1,'first') ;
    if ~isempty(ti) ;
        Threshold(c) = log10(Irel_unique(ti)) ;
    else
        Threshold(c) = 100 ;
    end
end

% peak rates
last_peak = psth_peak_delta_trls_mean_Irel(end,:) ; % max of tuning curve 

% half width and peak time of the response to the max stim
last_peak_time = psth_delta_peak_time_Irel(end,:) ;

last_halfWidth = psth_delta_halfWidth_Irel(end,:) ;

    
% %fit hockey stick - find threshold and gain for each cell
% for c = 1:NumOnCells ; % for each cell
%     TempCoef = nlinfit(log10(Irel_unique),psth_delta_peak_abs_Irel(:,c)',@HockeyStick,[-8,10,10], 'Weights', 10.^fliplr([1:length(Irel_unique)])) ;
%     Threshold(c) = TempCoef(1) ;
%     Hs_fit(c,:) = HockeyStick(TempCoef,log10(Irel_unique)) ;
% end

% histograms
peak_HistBins = [-200:5:200] ;
for ds = 1:size(psth_mean,1) ;
    peak_HistValues(ds,:) = hist(psth_delta_peak(ds,:),peak_HistBins) ; 
end

Threshold_HistBins = [log10(Irel_unique),100] ;
Threshold_HistValues= hist(Threshold,Threshold_HistBins) ; 

last_peak_HistBins = peak_HistBins ;
last_peak_HistValues = hist(last_peak,last_peak_HistBins) ; 

last_peak_time_HistBins = 0:.01:1 ;  
last_peak_time_HistValues = hist(last_peak_time,last_peak_time_HistBins) ;

last_halfWidth_HistBins = 0:.01:1 ;  
last_halfWidth_HistValues = hist(last_halfWidth,last_halfWidth_HistBins) ;

if PdfHistFlag ; % if you want probability not number
    Threshold_HistValues = Threshold_HistValues/sum(Threshold_HistValues) ; 
end 


% figures
Ftimes = unique(DataBlock(DB).DfParams.Ftime) ;
NDFs = unique(DataBlock(DB).DfParams.NDF) ;

for ds = 1:size(psth_mean,1) ;
    plot_i(ds) = (find(NDFs==DataBlock(DB).DfParams.NDF(ds))-1)*length(Ftimes) + find(Ftimes==DataBlock(DB).DfParams.Ftime(ds)) ; % subplot i
end

figure
for ds = 1:size(psth_mean,1) ;
    subplot(length(NDFs),length(Ftimes),plot_i(ds))
    plot(psthTime,psth_mean(ds,:),'k','linewidth',3) ;
    hold on
    %plot(psthTime,psth{ds},'c') ;
    plot(psthTime,psth_mean(ds,:)+psth_std(ds,:),'b:') ;
    plot(psthTime,psth_mean(ds,:)-psth_std(ds,:),'b:') ;
    text(.9,.9,['NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),',fdur ',num2str(DataBlock(DB).DfParams.Ftime(ds))],'units', 'norm') 
end

figure
for ds = 1:size(psth_mean,1) ;
    subplot(length(NDFs),length(Ftimes),plot_i(ds))
    plot(psthTime,psth_delta_mean(ds,:)) ;
    hold on
    text(.9,.9,['NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),',fdur ',num2str(DataBlock(DB).DfParams.Ftime(ds))],'units', 'norm') 
end

figure
for ds = 1:size(psth_mean,1) ;
    subplot(length(NDFs),length(Ftimes),plot_i(ds))
    plot(peak_HistBins,peak_HistValues(ds,:)) ;
    text(.9,.9,['NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),',fdur ',num2str(DataBlock(DB).DfParams.Ftime(ds))],'units', 'norm') 
end

figure
subplot(4,1,1)
plot(Threshold_HistBins,Threshold_HistValues) ;

subplot(4,1,2)
plot(last_halfWidth_HistBins,last_halfWidth_HistValues) ;

subplot(4,1,3)
plot(Threshold_HistBins,Threshold_HistValues) ;

subplot(4,1,4)
plot(Threshold_HistBins,Threshold_HistValues) ;

figure
plot(log10(Irel_unique),psth_delta_peak_Irel)
hold on
plot(log10(Irel_unique),psth_delta_peak_Irel_mean,'k','linewidth',2)
plot(log10(Irel_unique),psth_delta_peak_Irel_meanOn,'k:','linewidth',2)
plot(log10(Irel_unique),psth_delta_peak_Irel_meanOff,'k:','linewidth',2)

figure
plot(log10(Irel_unique),psth_delta_peak_abs_Irel)
pause
for c=1:NumOnCells ;
    plot(log10(Irel_unique),psth_delta_peak_abs_Irel(:,c))
    hold on
    plot(log10(Irel_unique),Hs_fit(c,:),'r')
    hold off
    pause
end

figure
plot(log10(Irel_unique),psth_delta_peak_abs_Snr_Irel)

figure
for c=1:NumOnCells ;
    
    subplot(1,2,1)
    for ds = 1:size(psth_mean,1) ;
        %subplot(length(NDFs),length(Ftimes),plot_i(ds))
        plot(psthTime,psth_delta{ds}(c,:)) ;
        hold on
        xlim([-.2, .75])
        %text(.9,.9,['NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),',fdur ',num2str(DataBlock(DB).DfParams.Ftime(ds))],'units', 'norm') 
    end
    hold off
    
    subplot(1,2,2)
    errorbar(log10(Irel_unique),psth_peak_delta_trls_mean_Irel(:,c),psth_peak_delta_trls_std_Irel(:,c),psth_peak_delta_trls_std_Irel(:,c))
    hold on 
    if Threshold(c)~=100;
        plot(Threshold(c),psth_peak_delta_trls_mean_Irel(Threshold(c)==log10(Irel_unique),c),'ro')
    end
    plot(log10(Irel_unique),psth_delta_peak_Irel(:,c),'g')
    hold off
    
    title(num2str(c))
    pause
end


% for igor
ExampleCell = 183 ;
ForIgor = struct() ;

VecName = ['logIrel','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,log10(Irel_unique)) ; % relative intensity vector

VecName = ['RspSpikeRateMean','Cell',num2str(ExampleCell),'Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psth_peak_delta_trls_mean_Irel(:,ExampleCell)) ; % mean of trials

VecName = ['RspSpikeRateStd','Cell',num2str(ExampleCell),'Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psth_peak_delta_trls_std_Irel(:,ExampleCell)) ; % std of trials


VecName = ['psthTime','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psthTime) ; % psth time vector

for ds = 1:size(psth_mean,1) ;
    VecName = ['psth','NDF',num2str(DataBlock(DB).DfParams.NDF(ds)),'Fdur',num2str(DataBlock(DB).DfParams.Ftime(ds)),'Cell',num2str(ExampleCell),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor, VecName, psth_delta{ds}(ExampleCell,:)) ; % psth vector
end


VecName = ['ThresholdHistBins','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,Threshold_HistBins) ; % threshold bins

VecName = ['ThresholdHistValues','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,Threshold_HistValues) ; % threshold hist values

VecName = ['lastPeakHistBins','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,last_peak_HistBins) ; % bins

VecName = ['lastPeakHistValues','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,last_peak_HistValues) ; % hist values

VecName = ['lastPeakTimeHistBins','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,last_peak_time_HistBins) ; % threshold bins

VecName = ['lastPeakTimeHistValues','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,last_peak_time_HistValues) ; % threshold hist values

VecName = ['lastHalfWidthHistBins','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,last_halfWidth_HistBins) ; % threshold bins

VecName = ['lastHalfWidthHistValues','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,last_halfWidth_HistValues) ; % threshold hist values

% avearge intensity respones curves
VecName = ['peak_mean','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psth_delta_peak_Irel_mean) ; % relative intensity vector

VecName = ['peak_mean_On','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psth_delta_peak_Irel_meanOn) ; % relative intensity vector

VecName = ['peak_mean_Off','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psth_delta_peak_Irel_meanOff) ; % relative intensity vector

    




