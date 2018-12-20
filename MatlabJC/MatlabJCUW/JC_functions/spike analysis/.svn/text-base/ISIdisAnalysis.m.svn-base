function ForIgor = ISIdisAnalysis(Input,Parameters,id,A) ;

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

% detect spikes 
if strcmp(id,'CA') ; % if cell attached spikes
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

% find spikes during stim only
isi = [];
for a = 1:length(SpikePnts) ;                                                               % for each spike epoch ...
    Stm_spikes{a} = find(SpikePnts{a}>PrePnts+STAPnts & SpikePnts{a}<PrePnts+StmPnts) ;    % indicies of SpikePnts vector during the stimuli and far enough out to get the prespike wave form 
    NumSpikes(a) = length(Stm_spikes{a}) ;                                                  % number of spikes during this period of time
    isi=[isi,(MSI*1000)*diff(SpikePnts{a}(Stm_spikes{a}))] ;                                % isi in ms as one long vector
end % end epoch loop
    
% isi distribution
[ISIdis,ISIbins] = hist(isi,[0:.5:1000]) ;
burstISI = find(isi<Parameters.isicut) ;
quietISI = find(isi>Parameters.isicut) ;
sumburst = length(burstISI) ; % number of isi below 20ms
sumquiet = length(quietISI) ; % number of isi above 20ms
meanburst = mean(isi(burstISI)) ; % mean of isi below 20ms
meanquiet = mean(isi(quietISI)) ; % mean of isi above 20ms

% get instantaneous firing rate plots
ifrTrain = data*0 ;

for a = 1:length(SpikePnts) ; % for each epoch
    ifr{a} = 1./(diff(SpikePnts{a})*MSI) ; % find the instantaneous firing rate (Hz) of each spike
    for b=1:length(ifr{a}) ; % for each spike
        ifrTrain(a,SpikePnts{a}(b):SpikePnts{a}(b+1)) = ifr{a}(b) ; % make the ifr the same until the next spike
    end
    %prep for export
    identifier = ['ifrTrain',num2str(a),id,num2str(A)] ;
    ForIgor.(identifier) = ifrTrain(a,:) ; 
end

% prep vectors for export
identifier = ['time',id,num2str(A)] ;
ForIgor.(identifier) = time ; 

identifier = ['ISIdis',id,num2str(A)] ;
ForIgor.(identifier) = ISIdis ;

identifier = ['ISIbins',id,num2str(A)] ;
ForIgor.(identifier) = ISIbins ;

identifier = ['meanburst',id,num2str(A)] ;
ForIgor.(identifier) = meanburst ;

identifier = ['meanquiet',id,num2str(A)] ;
ForIgor.(identifier) = meanquiet ;

identifier = ['sumburst',id,num2str(A)] ;
ForIgor.(identifier) = sumburst ;

identifier = ['sumquiet',id,num2str(A)] ;
ForIgor.(identifier) = sumquiet ;


figure % show example of spike detection
plot(time,data(3,:))
hold on
if strcmp(id,'CA') ; 
    plot(time(SpikePnts{3}(Stm_spikes{3})),15,'r*')
else
    plot(time(SpikePnts{3}(Stm_spikes{3})),Parameters.WCSpikeThresh,'r*')
end
title(['3example',id,num2str(A)])

figure % show last example of spike detection
plot(time,data(end,:))
hold on
if strcmp(id,'CA') ; 
    plot(time(SpikePnts{length(Stm_spikes)}(Stm_spikes{length(Stm_spikes)})),15,'r*')
else
    plot(time(SpikePnts{length(Stm_spikes)}(Stm_spikes{length(Stm_spikes)})),Parameters.WCSpikeThresh,'r*')
end
title(['lastexample',id,num2str(A)])

figure % show isi dis
plot(ISIbins,ISIdis) 
hold on
plot([meanburst,meanburst],[0,max(ISIdis)],'r')
plot([meanquiet,meanquiet],[0,max(ISIdis)],'r')
title(['isidis',id,num2str(A)])
axis([0 30 0 max(ISIdis)])

figure
plot(time,ifrTrain(3,:))
title(['3example',id,num2str(A)])

end % end function 