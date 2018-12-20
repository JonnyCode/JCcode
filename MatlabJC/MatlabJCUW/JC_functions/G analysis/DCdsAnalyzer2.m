function ForIgor = DCdsAnalyzer2(Input,Parameters,id,id2,A) ;

% adapted from earlier version, DCdsAnalyzer, but this one will focus on
% shuffle corrected cross correlations (bad normalization!) and normalizes the std in order to
% compare across many cells

% JC 1/5/10


epochs = str2num(Input(A).(id)) ; %sim and shuffled g
epochsSimRepeat = str2num(Input(A).(id2)) ; % sim repeated

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);

for a = 1:length(epochs) ;
    [Data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ; 
    [excg(a,:), inhg(a,:), error(a,:)] = ITCReadEpochStmGClamp(epochs(a), 0, fp);
end

for a = 1:length(epochsSimRepeat) ;
    [dataSimRepeat(a,:), error] = ITCReadEpoch(epochsSimRepeat(a), 0, fp) ; 
    [excgSimRepeat(a,:), inhgSimRepeat(a,:), error(a,:)] = ITCReadEpochStmGClamp(epochsSimRepeat(a), 0, fp);
end

[SI, error] = ITCGetSamplingInterval(epochs(1), fp); % get sampling interval    
SI = SI * 1e-6; % Sampling interval in sec
if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end
time = [SI:SI:SI*length(Data)] ;

SpikePnts = SpikeDetection_WC(Data,-30,10000) ;
SpikePntsSimRepeat = SpikeDetection_WC(dataSimRepeat,-30,10000) ;

spikeTrainSimRepeat = zeros(length(SpikePntsSimRepeat),length(Data)) ;
spikeTrainShuff = spikeTrainSimRepeat ;
spikeTrainSim = spikeTrainSimRepeat ;


round = 0 ;
for a = 1:2:length(epochs) ; % for each spike epoch
    round = round +1 ;

    dataSim(round,:) = Data(a,:) ;       
    dataShuff(round,:) = Data(a+1,:) ;       

    spikeTrainSim(round,SpikePnts{a}) = 1 ;
    spikeTrainShuff(round,SpikePnts{a+1}) = 1 ;

    spikeTimes_Sim{round} = time(SpikePnts{a}) ;
    spikeTimes_Shuff{round} = time(SpikePnts{a+1}) ;

    excgSim(round,:) = excg(a,:) ;
    inhgSim(round,:) = inhg(a,:) ;    
  
    excgShuff(round,:) = excg(a+1,:) ;
    inhgShuff(round,:) = inhg(a+1,:) ;
    
end

for a = 1:length(epochsSimRepeat)
    spikeTrainSimRepeat(a,SpikePntsSimRepeat{a}) = 1 ; 
    spikeTimes_SimRepeat{a} = time(SpikePntsSimRepeat{a}) ;
end

excgSim_Mean = mean(excgSim) ;
inhgSim_Mean = mean(inhgSim) ;

excgShuff_Mean = mean(excgShuff) ;
inhgShuff_Mean = mean(inhgShuff) ;

excgSimRepeat_Mean = mean(excgSimRepeat) ; %mean conductances
inhgSimRepeat_Mean = mean(inhgSimRepeat) ;



% tuning curve with error bars
groupTime = 3.0667 ; % (sec) time of each bar 
OnTime = .9 ; % (sec) time of on response

prePnts = 2000 ; % pnts added to g before g stim
groupPnts = floor(groupTime/SI(1));
OnPnts = floor(OnTime/SI(1)) ;

for a=1:8 ; %for each bar
    for trial = 1:length(epochs)/2 ; % for each trial
        SpikeRateSim(trial,a) = sum(spikeTrainSim(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a))/groupTime ; % hz
        SpikeRateShuff(trial,a) = sum(spikeTrainShuff(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a))/groupTime ; % hz
      
        SimExcG{a}(trial,:) = excgSim(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a) ;
        SimInhG{a}(trial,:) = inhgSim(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a) ;
        
        ShuffExcG{a}(trial,:) = excgShuff(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a) ;
        ShuffInhG{a}(trial,:) = inhgShuff(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a) ;
        
        acSimExc{a}(trial,:) = xcov(SimExcG{a}(trial,:)) ;
        acSimInh{a}(trial,:) = xcov(SimInhG{a}(trial,:)) ;
        
        acShuffExc{a}(trial,:) = xcov(ShuffExcG{a}(trial,:)) ;
        acShuffInh{a}(trial,:) = xcov(ShuffInhG{a}(trial,:)) ;       
        
        ccSim{a}(trial,:) = xcov(SimExcG{a}(trial,:),SimInhG{a}(trial,:)) ;
        
        ccShuff{a}(trial,:) = xcov(ShuffExcG{a}(trial,:),ShuffInhG{a}(trial,:)) ;
        
    end    
    SimExcG_mean(a,:) = mean(SimExcG{a}) ;
    SimInhG_mean(a,:) = mean(SimInhG{a}) ;
    
    ShuffExcG_mean(a,:) = mean(ShuffExcG{a}) ;
    ShuffInhG_mean(a,:) = mean(ShuffInhG{a}) ;
    
    acSimMean_exc(a,:) = mean(acSimExc{a}) ;
    acSimMean_inh(a,:) = mean(acSimInh{a}) ;
    
    acShuffMean_exc(a,:) = mean(acShuffExc{a}) ;
    acShuffMean_inh(a,:) = mean(acShuffInh{a}) ;
    
    ccMeanSim(a,:) = xcov(SimExcG_mean(a,:), SimInhG_mean(a,:)) ;
    ccMeanShuff(a,:) = xcov(ShuffExcG_mean(a,:), ShuffInhG_mean(a,:)) ;
    
    ccSim_mean(a,:) = mean(ccSim{a}) ;
    ccShuff_mean(a,:) = mean(ccShuff{a}) ;
    
    ccSim_shuffleCorrected(a,:) = (ccSim_mean(a,:) - ccMeanSim(a,:))./(max(sqrt(acSimMean_exc(a,:).*acSimMean_inh(a,:))) - ccMeanSim(a,:)) ;
    ccShuff_shuffleCorrected(a,:) = (ccShuff_mean(a,:) - ccMeanShuff(a,:))./(max(sqrt(acShuffMean_exc(a,:).*acShuffMean_inh(a,:))) - ccMeanShuff(a,:)) ;
    
    CrossCorrSim(a) = ccSim_shuffleCorrected(a,ceil(length(ccSim_shuffleCorrected)/2)) ;
    CrossCorrShuff(a) = ccShuff_shuffleCorrected(a,ceil(length(ccShuff_shuffleCorrected)/2)) ;
end

% mean and std of firing rate for each bar angle
meanSpikeRateSim = mean(SpikeRateSim) ;
meanSpikeRateShuff = mean(SpikeRateShuff) ;

stdSpikeRateSim = std(SpikeRateSim) ;
stdSpikeRateShuff = std(SpikeRateShuff) ;

% normalize std
stdSpikeRateSim_norm = stdSpikeRateSim/max(max(stdSpikeRateSim),max(stdSpikeRateShuff)) ;
stdSpikeRateShuff_norm = stdSpikeRateShuff/max(max(stdSpikeRateSim),max(stdSpikeRateShuff)) ;

% figure
figure
plot(CrossCorrSim,stdSpikeRateSim_norm,'k*') ;
hold on
plot(CrossCorrShuff,stdSpikeRateShuff_norm,'g*') ;



% for Igor


% std and cross corr coefs for comparsion plot

identifier = ['stdFRSimNorm','cell',num2str(A)] ;
ForIgor.(identifier) = stdSpikeRateSim_norm ; 

identifier = ['stdFRShuffNorm','cell',num2str(A)] ;
ForIgor.(identifier) = stdSpikeRateShuff_norm ;

identifier = ['corrCoefSim','cell',num2str(A)] ;
ForIgor.(identifier) = CrossCorrSim ; 

identifier = ['corrCoefShuff','cell',num2str(A)] ;
ForIgor.(identifier) = CrossCorrShuff ;    
    

