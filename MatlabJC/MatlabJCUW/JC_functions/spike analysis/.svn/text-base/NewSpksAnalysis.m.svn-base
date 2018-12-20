function ForIgor = NewSpksAnalysis(Input,Parameters,id1,id2,A) ;

% this function will compare the intantaneous firing rates from different
% traces and assess the STA when one cell is firing stronger or weaker than
% the cell under a different condition

% get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

% get spike data and sampling interval from first spike set
epochs = str2num(Input(A).(id1)) ;
for a = 1:length(epochs) ; % for each epoch  
    [predata{a}, error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get cell attached data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp);
    SI(a) = SI(a) * 1e-6; % Sampling interval in sec
end

% get spike data and sampling interval from second spike set
epochs2 = str2num(Input(A).(id2)) ;
for a = 1:length(epochs2) ; % for each epoch  
    [predata2{a}, error] = ITCReadEpoch(epochs2(a), 0, fp) ;    % get cell attached data
    [SI2(a), error] = ITCGetSamplingInterval(epochs2(a), fp);
    SI2(a) = SI2(a) * 1e-6; % Sampling interval in sec
end

% get light stimuli used to record g
cd ~/Data_analysis/Index/
load DCGstm % stimuli used when recorded conductances used in dc
prestm = DCGstm ; % stm

% rectifier and zero stim
negstmPnts = find(prestm<0) ;  % find indicies which would be getting a negative stim voltage
prestm(negstmPnts) = 0 ;       % make these points zero
prestmGeneric = prestm - repmat(mean(prestm(:,Parameters.PrePnts:Parameters.StmPnts),2),1,size(prestm,2)) ; % subtract off mean of stimulus during time varying stimulus

prestm = repmat(prestmGeneric,length(predata)/5,1) ; %repeat in a matrix as many times as needed
prestm2 = repmat(prestmGeneric,length(predata2)/5,1) ; %repeat in a matrix as many times as needed

% interpolate all signals so they are pseudo sampled at the same rate
MSI = min(min(SI),min(SI2)) ; % find the signal with the highest sampling rate
SLsec = ((Parameters.PrePnts+Parameters.StmPnts+Parameters.PostPnts)/10000) ; % entire stimulus length in sec
time = [MSI:MSI:SLsec] ; % time vector in sec 

for a = 1:length(epochs) ; % for each epoch
    data(a,:) = interp1([SI(a):SI(a):SLsec],predata{a},time,'linear','extrap') ;  % interpolate the data
    stm(a,:) = interp1([.0001:.0001:SLsec],prestm(a,:),time,'linear','extrap') ;  % interpolate the stimulus assuming an original sampling interval of .0001sec
end

for a = 1:length(epochs2) ; % for each epoch
    data2(a,:) = interp1([SI2(a):SI2(a):SLsec],predata2{a},time,'linear','extrap') ;  % interpolate the data
    stm2(a,:) = interp1([.0001:.0001:SLsec],prestm2(a,:),time,'linear','extrap') ;  % interpolate the stimulus assuming an original sampling interval of .0001sec
end

% change points to be appropriate for new sample rate
PrePnts = round(Parameters.PrePnts*.0001/MSI) ;
StmPnts = round(Parameters.StmPnts*.0001/MSI) ;
PostPnts = round(Parameters.PostPnts*.0001/MSI) ;
STAPnts = round(Parameters.STAPnts*.0001/MSI) ;
SmoothPnts = round(Parameters.SmoothPnts*.0001/MSI) ;

% detect spikes in whole cell attached data
SpikePnts = SpikeDetection_WC(data,-20,(1/MSI)) ; % data,threshold,samplerate
SpikePnts2 = SpikeDetection_WC(data2,-20,(1/MSI)) ; % data,threshold,samplerate

% get spiketrains
[PSTH,SpikeTrain] = PSTHmaker(SpikePnts,Parameters.SmoothPnts,length(data),1/MSI) ;   % PSTHmaker(SpikePnts,SmoothPnts,RawDataLength,SampleRate)
[PSTH2,SpikeTrain2] = PSTHmaker(SpikePnts2,Parameters.SmoothPnts,length(data),1/MSI) ;   % PSTHmaker(SpikePnts,SmoothPnts,RawDataLength,SampleRate)

% prep matricies
ifrTrain = SpikeTrain*0 ;
ifrTrain2 = SpikeTrain2*0 ;

% get instantaneous firing rates plots
for a = 1:length(SpikePnts) ; % for each epoch
    ifr{a} = 1./(diff(SpikePnts{a})*MSI) ; % find the instantaneous firing rate (Hz) of each spike
    for b=1:length(ifr{a}) ; % for each spike
        ifrTrain(a,SpikePnts{a}(b):SpikePnts{a}(b+1)) = ifr{a}(b) ; % make the ifr the same until the next spike
    end
end

% get instantaneous firing rates plots
for a = 1:length(SpikePnts2) ; % for each epoch
    ifr2{a} = 1./(diff(SpikePnts2{a})*MSI) ; % find the instantaneous firing rate (Hz) of each spike
    for b=1:length(ifr2{a}) ; % for each spike
        ifrTrain2(a,SpikePnts2{a}(b):SpikePnts2{a}(b+1)) = ifr2{a}(b) ; % make the ifr the same until the next spike
    end
end


% 
for a = 1:min(length(SpikePnts),length(SpikePnts2)) ; % for each possible comparison between set 1 and 2
    ifrTrainDiff(a,:) = ifrTrain(a,:)-ifrTrain2(a,:) ; % find the difference between the two ifr train
    
    NewSpikes(a,:) = SpikeTrain2(a,:) ; % make a new matrix of spikes from the set 2
    FailSpikes(a,:) = SpikeTrain(a,:) ; % make a new matric of spikes from the set 1
    
    NewSpikes{a}(:,ifrTrainDiff(a,:)>=0) = 0 ; % get rid of the spikes that were in set 1 but not set 2
    FailSpikes{a}(:,ifrTrainDiff(a,:)<=0) = 0 ; % get rid of the spikes that were in set 2 but not set 1
    
    ifrTraingroup{a} = ifrTrain(group,:) ; % group ifrTrain
    ifrTraingroupMean{a} = mean(ifrTraingroup{a}) ; % mean ifrTrain
    SpikeTraingroup{a} = SpikeTrain(group,:) ; % group spiketrain
end

% group light stim, Spike response, and PSTH for each unique light stim
for a = 1:5 ; % there are five unique dynamic clamp sets (from 5 unique light stim)
    group2 = [a:5:length(SpikePnts2)] ; % these are repeated in sets so every 5th trial is from the same light stim
    stmgroup2{a} = stm2(group2,:) ; % group light stim
    ifrTraingroup2{a} = ifrTrain2(group2,:) ; % group ifrTrain
    ifrTraingroupMean2{a} = mean(ifrTraingroup2{a}) ; % mean ifrTrain
    SpikeTraingroup2{a} = SpikeTrain2(group2,:) ; % group spiketrain
    groupDiff = ifrTraingroupMean{a} - ifrTraingroupMean2{a} ; % subtract mean ifr trains
end

% compare PSTH from set 1 with PSTH from set 2

% find STA during times when there is a difference between the two PSTHs
for a=1:5 % for each unique light stim
    NewSpikes{a} = SpikeTraingroup2{a} ; % make a new matrix of spikes from the set 2
    FailSpikes{a} = SpikeTraingroup{a} ; % make a new matric of spikes from the set 1
    
    NewSpikes{a}(:,groupDiff{a}>=0) = 0 ; % keep only the spikes that were in set 2 but not set 1
    FailSpikes{a}(:,groupDiff{a}<=0) = 0 ; % keep only the spikes that were in set 1 but not set 2
    
    epochNSTAs{a} = ifft(fft(NewSpikes{a}(:,PrePnts+STAPnts:PrePnts+StmPnts),[],2).*conj(fft(stmgroup2{a}(:,PrePnts+STAPnts:PrePnts+StmPnts),[],2)),[],2) ; % correlation of each new spike train with the stim
    epochFSTAs{a} = ifft(fft(FailSpikes{a}(:,PrePnts+STAPnts:PrePnts+StmPnts),[],2).*conj(fft(stmgroup{a}(:,PrePnts+STAPnts:PrePnts+StmPnts),[],2)),[],2) ; % correlation of each failed spike train with the stim
    
        for b=1:size(epochNSTAs{a},1) ; % for each epochs sta for the new spikes
            epochNSTAs{a}(b,:) = epochNSTAs{a}(b,:)/sum(NewSpikes{a}(b,PrePnts+STAPnts:PrePnts+StmPnts)) ; %normalize each epoch STA by the number of spikes
        end
        
        for b=1:size(epochFSTAs{a},1) ; % for each epochs sta for the new spikes
            epochFSTAs{a}(b,:) = epochFSTAs{a}(b,:)/sum(FailSpikes{a}(b,PrePnts+STAPnts:PrePnts+StmPnts)) ; %normalize each epoch STA by the number of spikes
        end
        
    NewSpikeSTA(a,:) = mean(epochNSTAs{a}) ; % average STA across epochs
    FailSpikeSTA(a,:) = mean(epochFSTAs{a}) ; 
    
end

NewSpikeSTA = mean(NewSpikeSTA) ;
FailSpikeSTA = mean(FailSpikeSTA) ;


% prep export

identifier = ['NewSTA',id1,id2,num2str(A)] ;
ForIgor.(identifier) = NewSpikeSTA;

identifier = ['FailSTA',id1,id2,num2str(A)] ;
ForIgor.(identifier) = FailSpikeSTA;

figure
plot(data(2,:)) 
hold on
plot(SpikePnts{2},-20,'k*')
plot(data2(2,:),'r')
plot(SpikePnts2{2},-20,'k*')

figure
plot(time,PSTHgroupMean{1},'b')
hold on
plot(time,PSTHgroup2Mean{1},'r')
plot(time,PSTHgroupDiff{1},'g')

figure
plot(NewSpikeSTA)
hold on
plot(FailSpikeSTA,'r')
legend('-inh spikes','+inhspikes')
title('STA')

end


