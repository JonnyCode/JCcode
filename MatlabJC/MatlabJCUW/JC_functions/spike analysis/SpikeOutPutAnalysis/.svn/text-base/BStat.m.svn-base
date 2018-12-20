function ForIgor = BStat(Input,Parameters,id,A) ;

% this function will use my Burststatfinder to identify spikes inside and outside of defined bursts 

% get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each spike epoch
    [predata{a}, error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get cell attached data
    [SI(a), error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SI(a) = SI(a) * 1e-6; % Sampling interval in sec
end

% interpolate all signals so they are pseudo sampled at the same rate
MSI = min(SI) ; % find the signal with the highest sampling rate
SLsec = ((Parameters.PrePnts+Parameters.StmPnts+Parameters.PostPnts)/10000) ; % entire stimulus length in sec
time = [MSI:MSI:SLsec] ; % time vector in sec 
for a = 1:length(epochs) ; % for each epoch
    data(a,:) = interp1([SI(a):SI(a):SLsec],predata{a},time,'linear','extrap') ;  % interpolate the data
end
if ~strcmp(id,'CA') ; % if it is a whole cell recording
    data = data - 10 ; % account for 10mV liquid junction potential
end

% detect spikes in data
if strcmp(id,'CA') ;
    SpikePnts = SpikeDetection(data,15,(1/MSI)) ; % data,threshold,samplerate
else
    SpikePnts = SpikeDetection_WC(data,Parameters.WCSpikeThresh,(1/MSI)) ; % data,threshold,samplerate
end
    
% change points to be appropriate for new sample rate
PrePnts = round(Parameters.PrePnts*.0001/MSI) ;
StmPnts = round(Parameters.StmPnts*.0001/MSI) ;
PostPnts = round(Parameters.PostPnts*.0001/MSI) ;
STAPnts = round(Parameters.STAPnts*.0001/MSI) ;
DecimatePnts = round(Parameters.DecimatePnts*.0001/MSI) ;
QuietPnts = round(Parameters.QuietPnts*.0001/MSI) ; 


% find spikes during stim only
isi = [];
for a = 1:length(SpikePnts) ;                                                               % for each spike epoch ...
    Stm_spikes{a} = find(SpikePnts{a}>PrePnts+STAPnts & SpikePnts{a}<PrePnts+StmPnts) ;     % indicies of SpikePnts vector during the stimuli and far enough out to get the prespike wave form 
    NumStmSpikes(a) = length(Stm_spikes{a}) ;                                                   % number of spikes to be used per trial for STA
    isi=[isi,(MSI*1000)*diff(SpikePnts{a}(Stm_spikes{a}))] ;
    forburst{a} = SpikePnts{a}(Stm_spikes{a}) ;                                                 % this cell for burst detection
end % end epoch loop
NumStmspikesMean = mean(NumStmSpikes) ; % mean number of spikes per epoch as selected during the stimuli

% burst detection
[BurstPnts,LoneSpkPnts,NumSpikes,BurstISI,FracBurstSpikes] = BurstStatFinder(forburst,QuietPnts,Parameters.MinSpikeNumber) ; %(SpikePnts,QuietPnts,MinSpikeNumber)

for a = 1:length(SpikePnts) ; %per epoch
    NumBursts(a) = length(BurstPnts{a}) ;               % number burst 
    NumSpikeMean(a) = mean(NumSpikes{a}) ;              % mean number spikes per burst
    BurstISIMean(a) = (mean(BurstISI{a})*MSI*1000) ;    % mean isi in ms
end

NumBurstsMean = mean(NumBursts(a)) ;
NumSpikeMean2 = mean(NumSpikeMean) ;
BurstISIMean2 = mean(BurstISIMean) ;
FracBurstSpikeMean = mean(FracBurstSpikes) ; 

figure
rastor = gcf ;
for a = 1:length(SpikePnts) ;
    figure(rastor)
    subplot(5,1,a-5*(ceil(a/5)-1)) % spike rastors ploted in groups by same stim
    plot(SpikePnts{a}*MSI,epochs(a),'b.')
    hold on
    plot(forburst{a}*MSI,epochs(a),'g.')
    if NumBursts(a)~=0 ; % if there are any bursts
        plot(BurstPnts{a}*MSI,epochs(a),'ro')
        plot(LoneSpkPnts{a}*MSI,epochs(a),'co')
    end
end
title([id,num2str(A)])

% prep vectors for igor plotting
identifier = ['BurstStats',id,num2str(A)] ;
ForIgor.(identifier) = [NumStmspikesMean, FracBurstSpikeMean, NumBurstsMean, NumSpikeMean2,BurstISIMean2] ;
% mean per epoch: number of spikes, fraction in burst, number of bursts, mean number spikes per burst, mean isi per burst(ms),   

for a = 1:length(SpikePnts) ; % for each 
identifier = ['rastor',num2str(a),id,num2str(A)] ;
ForIgor.(identifier) = SpikePnts{a}*MSI ; % spike times in sec
end

identifier = ['time',id,num2str(A)] ;
ForIgor.(identifier) = time ;

figure
text(.1,.9,num2str([NumStmspikesMean, FracBurstSpikeMean, NumBurstsMean, NumSpikeMean2,BurstISIMean2]')) ;
title([id,num2str(A)])

figure % show 3rd example of spike detection
plot(time,data(3,:))
hold on
if strcmp(id,'CA') ; 
    plot(time(SpikePnts{3}(Stm_spikes{3})),15,'r*')
else
    plot(time(SpikePnts{3}(Stm_spikes{3})),Parameters.WCSpikeThresh,'r*')
end
title(['example',id,num2str(A)])

figure % show last example of spike detection
plot(time,data(end,:))
hold on
if strcmp(id,'CA') ; 
    plot(time(SpikePnts{length(Stm_spikes)}(Stm_spikes{length(Stm_spikes)})),15,'r*')
else
    plot(time(SpikePnts{length(Stm_spikes)}(Stm_spikes{length(Stm_spikes)})),Parameters.WCSpikeThresh,'r*')
end
title(['lastexample',id,num2str(A)])

end