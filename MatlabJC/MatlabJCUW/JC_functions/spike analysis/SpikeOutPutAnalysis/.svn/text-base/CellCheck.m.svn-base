function ForIgor = CellCheck(Input,Parameters,id,A) ;

% this function will help assess quality of CA, IC, and DC epochs and show
% spike detection

% get data
[fp, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
for a = 1:length(epochs) ; % for each spike epoch
    [predata{a}, error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data
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

% get power spectrum of raw data to assess oscillations etc
[powerspec_xvalues, mean_powerspec] = PowerSpectrumFinder(data,1/MSI) ;

% detect spikes 
if strcmp(id,'CA') ; % if cell attached spikes
    SpikePnts = SpikeDetection(data,10,(1/MSI)) ; % data,threshold,samplerate
else    
    SpikePnts = SpikeDetection_WC(data,Parameters.WCSpikeThresh,(1/MSI)) ; % data,threshold,samplerate
end

% change points to be appropriate for new sample rate
PrePnts = round(Parameters.PrePnts*.0001/MSI) ;
StmPnts = round(Parameters.StmPnts*.0001/MSI) ;
PostPnts = round(Parameters.PostPnts*.0001/MSI) ;
STAPnts = round(Parameters.STAPnts*.0001/MSI) ;
DecimatePnts = round(Parameters.DecimatePnts*.0001/MSI) ;

% find spikes rates during and before the stim
for a = 1:length(SpikePnts) ;                                                       % for each spike epoch ...
    Stm_spikes{a} = find(SpikePnts{a}>PrePnts+STAPnts & SpikePnts{a}<PrePnts+StmPnts) ; % indicies of SpikePnts vector during the stimuli and far enough out to get the prespike wave form (ie. to actually have been caused by the light) 
    RateStmSpikes(a) = length(Stm_spikes{a})/((StmPnts-STAPnts)*MSI) ;                  % rate of spikes (Hz) during the stim = number spikes/ number sample points * (length(sec) of each sample point)
    PreStm_spikes{a} = find(SpikePnts{a}<PrePnts & SpikePnts{a}>STAPnts) ;              % indicies of SpikePnts vector before the stim (but avoid early response because it is akin to lights turing on)
    RatePreStmSpikes(a) = length(PreStm_spikes{a})/((PrePnts-STAPnts)*MSI) ;                      % rate of spikes (Hz) before the stim
end % end epoch loop

meanStmSpikeRate = mean(RateStmSpikes) ;
meanPreStmSpikeRate = mean(RatePreStmSpikes) ; 

stdStmSpikeRate = std(RateStmSpikes) ;
stdPreStmSpikeRate = std(RatePreStmSpikes) ; 

% find subthreshold voltage stats before and during stim
if ~strcmp(id,'CA') ; % if this is a whole cell recording

lpdata = lowPassFilter(data, 1/MSI, 100) ; % low pass filter to get rid of spikes

meanPreStmVmean = mean(mean(lpdata(:,1:PrePnts),2)) ;                       % avearge across epochs of Voltage mean before stm
meanStmVmean = mean(mean(lpdata(:,PrePnts+STAPnts:PrePnts+StmPnts),2)) ;    % avearge across epochs of Voltage mean during stm
diffVmean = meanPreStmVmean - meanStmVmean ;                                % difference between mean prestm voltage and mean stm voltage

stdPreStmVmean = std(mean(lpdata(:,1:PrePnts),2)) ;                     % std across epochs of Voltage mean before stm
stdStmVmean = std(mean(lpdata(:,PrePnts+STAPnts:PrePnts+StmPnts),2)) ;  % std across epochs of Voltage mean during stm

meanPreStmVstd = mean(std(lpdata(:,1:PrePnts),[],2)) ;                     % avearge across epochs of Voltage std before stm
meanStmVstd = mean(std(lpdata(:,PrePnts+STAPnts:PrePnts+StmPnts),[],2)) ;  % avearge across epochs of Voltage std during stm

stdPreStmVstd = std(std(lpdata(:,1:PrePnts),[],2)) ;                       % std across epochs of Voltage std before stm
stdStmVstd = std(std(lpdata(:,PrePnts+STAPnts:PrePnts+StmPnts),[],2)) ;    % std across epochs of Voltage std during stm

figure % plot an example
plot(time,data(3,:)) ;
hold on
plot(time,lpdata(3,:),'r') ;
plot(time(1:PrePnts),ones(1,PrePnts)*meanPreStmVmean,'g-')
plot(time(PrePnts:PrePnts+StmPnts),ones(1,StmPnts+1)*meanStmVmean,'g-')
title([id,num2str(A)])

% prep vectors for export
identifier = ['meanPreStmVmean',id,num2str(A)] ;
ForIgor.(identifier) = meanPreStmVmean ;

identifier = ['meanStmVmean',id,num2str(A)] ;
ForIgor.(identifier) = meanStmVmean ;

identifier = ['diffVmean',id,num2str(A)] ;
ForIgor.(identifier) = diffVmean ;

identifier = ['stdPreStmVmean',id,num2str(A)] ;
ForIgor.(identifier) = stdPreStmVmean ;

identifier = ['stdStmVmean',id,num2str(A)] ;
ForIgor.(identifier) = stdStmVmean ;

identifier = ['meanPreStmVstd',id,num2str(A)] ;
ForIgor.(identifier) = meanPreStmVstd ;

identifier = ['meanStmVstd',id,num2str(A)] ;
ForIgor.(identifier) = meanStmVstd ;

identifier = ['stdPreStmVstd',id,num2str(A)] ;
ForIgor.(identifier) = stdPreStmVstd ;

identifier = ['stdStmVstd',id,num2str(A)] ;
ForIgor.(identifier) = stdStmVstd ;

end

identifier = ['meanPreStmSpikeRate',id,num2str(A)] ;
ForIgor.(identifier) = meanPreStmSpikeRate ;

identifier = ['meanStmSpikeRate',id,num2str(A)] ;
ForIgor.(identifier) = meanStmSpikeRate ;

identifier = ['stdPreStmSpikeRate',id,num2str(A)] ;
ForIgor.(identifier) = stdPreStmSpikeRate ;

identifier = ['stdStmSpikeRate',id,num2str(A)] ;
ForIgor.(identifier) = stdStmSpikeRate ;

identifier = ['xpowspec',id,num2str(A)] ;
ForIgor.(identifier) = powerspec_xvalues ;

identifier = ['powspec',id,num2str(A)] ;
ForIgor.(identifier) = mean_powerspec ;

% plot spike all spike detections
for a=1:length(SpikePnts) ; % for all examples
    figure 
    plot(time,data(a,:))
    title(num2str(epochs(a)))
    hold on
    if strcmp(id,'CA') ;
        plot(time(SpikePnts{a}),15,'k*') ;
        plot(time(SpikePnts{a}(Stm_spikes{a})),15,'r*') ;
    else
        plot(time(SpikePnts{a}),Parameters.WCSpikeThresh,'k*') ;
        plot(time(SpikePnts{a}(Stm_spikes{a})),Parameters.WCSpikeThresh,'r*') ;
    end
end

% plot power spectrum of raw data
figure
plot(powerspec_xvalues, mean_powerspec)
title([id,num2str(A)])
a = gca ;
set(a,'Xscale','log')
set(a,'Yscale','log')
