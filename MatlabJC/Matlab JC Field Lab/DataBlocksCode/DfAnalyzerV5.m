function ForIgor = DfAnalyzerV5(DataBlock, DB, Params)

% adapted from DfAnalyzerV3 (not V4) to do 2AFC instead of variance

% JC 7/14/15

% parameters
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
interFlashIntVar = .005 ; % (sec) expected precision of trigger intervals
minNumFlashes = 5 ; % min number of flashes in a block 
psth_bin_size = 0.1 ; % (sec) 
LED_settings_calibrations_path = '/Volumes/lab/Experiments/Calibration/LED_settings_calibrations' ; % location of LED calibration mat file
LoadCalibrationFlag = true ;
PdfHistFlag = true ; % make histograms in pdfs if true
TestTime = 1.5 ; % (sec) before trigger (dark) and after (flash)

% load stim data
if LoadCalibrationFlag ;
    load /Volumes/lab/Experiments/Calibration/NdfCalibration
    Irel = (DataBlock(DB).DfParams.Ftime/1000).*NdfCalibration(2,DataBlock(DB).DfParams.NDF+1) ;
else
    Irel = DataBlock(DB).DfParams.Ftime./10.^DataBlock(DB).DfParams.NDF ;
end
Irel_unique = unique(Irel) ;

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

% select trigger sets that have enough trials to be counted
ds= 0 ;
for ts=1:length(trigger_set_i) ; % for each trigger set 
    if length(trigger_set_i{ts})>minNumFlashes ; % if more than X flashes
        ds=ds+1 ;
        select_trigger_set_i{ds} = trigger_set_i{ts} ;
    end
end

if length(select_trigger_set_i) ~= length(DataBlock(DB).DfParams.Ftime) ; % if the number of identified flash blocks is not as expected
    disp('Params or data length mismatch')
end

% group trigger sets with the same stimuli and order by stim intensity

for ds = 1:length(Irel_unique) ; % for each stimulus
    full_trigger_set_i{ds} = [] ;
    Ireli = find(Irel == Irel_unique(ds)) ; % find the 
    for ts=1:length(Ireli) ;
        full_trigger_set_i{ds} = [full_trigger_set_i{ds},select_trigger_set_i{Ireli(ts)}] ;
    end
end

UniqueStimNum = length(Irel_unique) ; % number of different flashes

% get psth in each cell for each set of stimuli 
psthTime = [-TestTime:psth_bin_size:TestTime] ;
for ds=1:UniqueStimNum ; % for each stimulus 
    for c = 1:NumCells ; % for each cell
        for trls = 1:length(full_trigger_set_i{ds}) ; % for each trial

            spikes = dataRun.spikes{c}-dataRun.triggers(full_trigger_set_i{ds}(trls));
            psthTemp = hist(spikes(spikes>=-TestTime & spikes<=TestTime),psthTime) ;
            
            %[psthTemp,psthTime] = get_psth(dataRun.spikes{c}, dataRun.triggers(full_trigger_set_i{ds}(trls)), 'start', -TestTime, 'stop', TestTime, 'bin_size', psth_bin_size) ; 
            psth{ds}{c}(trls,:) = psthTemp ;

            StimPnt = find(psthTime==0) ; % flash point
            preStimPnts = [1:StimPnt] ; % points before flash
            postStimPnts = [StimPnt:length(psthTime)] ; % points after flash

            preStimMean = mean(psthTemp(preStimPnts)) ; % mean before flash

            psthPostFlash{ds}{c}(trls,:) = psthTemp(postStimPnts) - preStimMean ; % post flash response change (flash response)
            psthPreFlash{ds}{c}(trls,:) = psthTemp(preStimPnts) - preStimMean ; % pre (dark) flash
        end
        psth_mean{ds}(c,:) = mean(psth{ds}{c}) ;
        psthPostFlash_sum{ds}(c,:) = sum(psthPostFlash{ds}{c}) ;
        psthPreFlash_sum{ds}(c,:) = sum(psthPreFlash{ds}{c}) ;
    end
end   

% 2 AFC test
for ds = 1:UniqueStimNum ; % for each stim
    NumTrials(ds) = length(psthPostFlash{ds}{c}) ;
    for c = 1:NumCells ; % for each cell
        for trls = 1:NumTrials(ds) ; % for each trial
            flashExempt = psthPostFlash_sum{ds}(c,:) - psthPostFlash{ds}{c}(trls,:) ; 
            flashExempt = flashExempt/norm(flashExempt) ;
            
            %darkExempt = psthPreFlash_sum{ds}(c,:) - psthPreFlash{ds}{c}(trls,:) ; 
            %darkExempt = darkExempt/norm(darkExempt) ;
            
            %Discriminant = flashExempt - darkExempt ;
            
            flashDot{ds}(c,trls) = psthPostFlash{ds}{c}(trls,:)* flashExempt' ;
            darkDot{ds}(c,trls) = psthPreFlash{ds}{c}(trls,:)* flashExempt' ;
        end
        Pc(c,ds) = (sum(flashDot{ds}(c,:) > 0) + sum(darkDot{ds}(c,:) <= 0))/(NumTrials(ds)*2); % percent correct
    end
end
    
% Thresholds
for c = 1:NumCells ; % for each cell
    ti = find(Pc(c,:)>0.84,1,'first') ;
    if ~isempty(ti) ;
        Threshold(c) = log10(Irel_unique(ti)) ;
    else
        Threshold(c) = 100 ;
    end
end

% histograms
Threshold_HistBins = [Irel_unique,100] ;  
Threshold_HistValues = hist(Threshold,Threshold_HistBins) ;

Pc_at_04_i = find(Irel_unique>.04,1,'first') ; 
Pc_at_04_HistBins = [0:.05:1] ; 
Pc_at_04_HistValues = hist(Pc(:,Pc_at_04_i),Pc_at_04_HistBins) ;


if PdfHistFlag ;
    Threshold_HistValues = Threshold_HistValues/sum(Threshold_HistValues) ;
    Pc_at_04_HistValues = Pc_at_04_HistValues/sum(Pc_at_04_HistValues) ;
end

% figures
figure % all cell flash mean
for ds = 1:UniqueStimNum ; % for each stim 
    subplot(1,UniqueStimNum,ds)
    plot(psthTime,mean(psth_mean{ds})) ;
    title(['Trial#: ',num2str(NumTrials(ds))])
end

figure % single cell curves
for c = 1:NumCells ; % for each cell
    subplot(1,2,1)
    for ds = 1:UniqueStimNum ; % for each stim  
        plot(psthTime,psth_mean{ds}(c,:)) ;
        xlabel('time (s)')
        ylabel('spike number')
        hold on
    end
    hold off
    
    subplot(1,2,2)
    plot(log10(Irel_unique),Pc(c,:))
    xlabel('log 10 Rh*/rod')
    ylabel('percent correct')
    
    pause
end
    
figure % histograms
subplot(1,2,1)
plot(log10(Threshold_HistBins),Threshold_HistValues)

subplot(1,2,2)
plot(Pc_at_04_HistBins,Pc_at_04_HistValues)

% for igor
ExampleCell = 78 ;
ForIgor = struct() ;

VecName = ['logIrel','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,log10(Irel_unique)) ; % relative intensity vector

VecName = ['Pc','Cell',num2str(ExampleCell),'Db',num2str(DB)] ;
ForIgor = setfield(ForIgor, VecName, Pc(ExampleCell,:)) ; % Pc as function of flash for example cell


% psth example cell
VecName = ['psthTime','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,psthTime) ; % psth time vector

for ds = 1:UniqueStimNum ; % for each stim
    VecName = ['psthIrel',num2str(Irel_unique(ds)),'Cell',num2str(ExampleCell),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor, VecName, psth_mean{ds}(ExampleCell,:)) ; % psth vector
end


VecName = ['ThresholdHistBins','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,Threshold_HistBins) ; % threshold bins

VecName = ['ThresholdHistValues','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,Threshold_HistValues) ; % threshold hist values


VecName = ['PcAtOneFlashHistBins','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,Pc_at_04_HistBins) ; % Pc hist bins

VecName = ['PcAtOneFlashHistValues','Db',num2str(DB)] ;
ForIgor = setfield(ForIgor,VecName,Pc_at_04_HistValues) ; % Pc hist values   




