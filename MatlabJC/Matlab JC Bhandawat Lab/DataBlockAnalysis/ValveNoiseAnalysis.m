function ForIgor = ValveNoiseAnalysis(Input,id,A)

% this function will look at spike odor responses during a random flicker
% of the odor valve and extract stim to spike and stim to voltage filters
% JC 5/14/12

% parameters
sampleRate = 10000 ; % temp hard coded - should be saved in file 
driftCheckTime = 0.25 ; %(sec) time at begining and end which current injected is inspected for changes
absRefTime = 0.002 ; % (sec) 
minRise = 2 ; % (mV)
minFall = 1 ; % (mV)
PulseRecTime = 0.01 ; % (sec) time after current pulse that spontaneous firing rate can be assessed


% load data in matricies
rootdir = ['C:\Cafaro Data\', Input(A).cellname(1:6),'Data'];

odorRspTrials = str2num(Input(A).(id)) ;
numTrials = length(odorRspTrials) ; 

for a = 1:numTrials ;
    temp = load([rootdir,'\',Input(A).cellname,'\','voltage_',Input(A).cellname,'_',num2str(odorRspTrials(a))]) ;
    vData(a,:) = temp.voltage ;
    
    temp = load([rootdir,'\',Input(A).cellname,'\','current_',Input(A).cellname,'_',num2str(odorRspTrials(a))]) ;
    iData(a,:) = temp.current ;
    
    temp = load([rootdir,'\',Input(A).cellname,'\','Ao0_',Input(A).cellname,'_',num2str(odorRspTrials(a))]) ;
    ao0Data(a,:) = temp.Ao0 ; % odor
    
    temp = load([rootdir,'\',Input(A).cellname,'\','Ao1_',Input(A).cellname,'_',num2str(odorRspTrials(a))]) ;
    ao1Data(a,:) = temp.Ao1 ; % vext   
    
    temp = load([rootdir,'\',Input(A).cellname,'\','TrigTime_',Input(A).cellname,'_',num2str(odorRspTrials(a))]) ;
    tData(a) = temp.Trigtime ;
end

% time vector
time = [1:size(vData,2)]/sampleRate ;
tDataN = (tData - min(tData))*24*60^2 ; % convert to seconds since experiment began

% make odor valve pulse binary
ao0DataB = ao0Data ; 
ao0DataB(ao0Data>=5) = 1 ; 
ao0DataB(ao0Data<5) = 0 ;

% round Vext pulse to nearest 10 mV and get rid of single sample pulses
ao1DataR = round(ao1Data*100)/100 ; 
for a = 1:numTrials ;
    for b=1:length(ao1DataR(a,:))-2 ;
        if ao1DataR(a,b)~=ao1DataR(a,b+1) && ao1DataR(a,b)==ao1DataR(a,b+2);
            ao1DataR(a,b+1)= ao1DataR(a,b) ;
        end
    end
end

% make sure odor pulse was the same time
ao0DataBdiff = ao0DataB - repmat(ao0DataB(1,:),numTrials,1) ;
if sum(abs(ao0DataBdiff(:)))~=0 ;
    disp('odor pulse discrepancy') ;
end


% make sure the R input check was at the same time
ao1DataRdiff = ao1DataR - repmat(ao1DataR(1,:),numTrials,1) ;
if sum(abs(ao1DataRdiff(:)))~=0 ;
    disp('Vext pulse discrepancy') ;
end

% check that input current is not changing substantially during any of the trials
driftCheckPnts = driftCheckTime*sampleRate ;
for a = 1:numTrials ;
    Idrift(a) = mean(iData(a,1:driftCheckPnts)) - mean(iData(a,end-driftCheckPnts:end)) ;

    if Idrift(a)>1 ;
        disp(['significant I drift in trial',num2str(a)]) ;
    end
end

% index of current pulse
iipb = find(ao1DataR(1,:)~=0,1,'first') ; % current pulse beginging
iipe = find(ao1DataR(1,:)~=0,1,'last') ; % current pulse end

% detect spikes in voltage data
[spikePnt,SpikeData,NonSpikeData] = spikeFinder(vData,sampleRate,absRefTime,minRise,minFall) ;

% assess resting potential
for a = 1:numTrials ;
    Vrest(a) = mean(vData(a,1:iipb)) ; % g ohms
end

% spike trains
spikeTrains = zeros(size(vData)) ;
for a = 1:numTrials ;
    spikeTrains(a,spikePnt{a}) = 1 ;
end

% psth
psth = smooth(mean(spikeTrains,1),sampleRate*.01)' ; % 1 ms bins

% linear filter
ao0DataB_RS = reshape(ao0DataB(1,:)',sampleRate*2,length(ao0DataB(1,:))/(sampleRate*2))' ;
psth_RS = reshape(psth',sampleRate*2,length(ao0DataB(1,:))/(sampleRate*2))' ;

[LinearFilter] = LinFilterFinder(ao0DataB_RS,psth_RS, sampleRate, 5) ;
[LinearFilter] = LinFilterFinder(ao0DataB(1,:),mean(vData), sampleRate, 5) ;

% sta
sta = staFinder(ao0DataB,spikeTrains,10000) ;

% figures

figure % spike detection
for a = 1:numTrials ;
    plot(time,vData(a,:)) 
    hold on
    plot(time(spikePnt{a}),vData(a,spikePnt{a}),'r*')
    hold off
    pause
end

figure % voltage traces
subplot(5,1,1)
plot(time,mean(vData))

subplot(5,1,2:5) % spike raster
for a = 1:numTrials ;
    for b=1:length(spikePnt{a}) ;
        plot([1,1]*spikePnt{a}(b),[a-1,a])
        hold on
    end
end

figure
plot(time,ao0DataB(1,:))
hold on
plot(time,psth/max(psth),'r')

figure
plot(sta)

