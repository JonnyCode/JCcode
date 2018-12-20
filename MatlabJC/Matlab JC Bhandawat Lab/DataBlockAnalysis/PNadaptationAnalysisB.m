function ForIgor = PNadaptationAnalysisB(Input,id1,id2,A)

% this function will look at PN odor responses +/- background
% JC 4/24/12

% parameters
sampleRate = 10000 ; % (hz) temp hard coded - should be saved in file 
driftCheckTime = 0.25 ; %(sec) time at begining and end which current injected is inspected for changes
absRefTime = 0.002 ; % (sec) 
minRise = .4 ; % (mV)
minFall = .3 ; % (mV)
PulseRecTime = 0.01 ; % (sec) time after current pulse that spontaneous firing rate can be assessed

% load data in matricies
rootdir = ['C:\Cafaro Data\', Input(A).cellname(1:6),'Data'];

odorRspTrials = str2num(Input(A).(id1)) ;
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

% index of current pulse, odor pulse, and odor response times
iopb = find(ao0DataB(1,:)~=0,1,'first')-1 ; % odor pulse begining
iope = find(ao0DataB(1,:)~=0,1,'last') ; % odor pulse ending

iipb = find(ao1DataR(1,:)~=0,1,'first') ; % current pulse beginging
iipe = find(ao1DataR(1,1:iopb)~=0,1,'last') ; % current pulse end

vData_mean = mean(vData) ;
vData_mean_2thresh =mean(vData_mean(iipe:iopb))+2*(max(vData_mean(iipe:iopb))-mean(vData_mean(iipe:iopb))) ; % 2*max
iorb = find(vData_mean>vData_mean_2thresh,1,'first') ; % first point above thresh - odor response begining
iore = find(vData_mean>vData_mean_2thresh,1,'last') ; % last point above thresh - odor response ending

% detect spikes in voltage data
[spikePnt,SpikeData,NonSpikeData] = spikeFinder(vData,sampleRate,absRefTime,minRise,minFall) ;

% assess resting potential
for a = 1:numTrials ;
    Vrest(a) = mean(vData(a,(iipe+driftCheckPnts):iopb)) ; % g ohms
end

% odor response plateau potential
for a = 1:numTrials ;
    PlatPot(a) = mean(vData(a,iope-driftCheckPnts:iope)) ; %mV
    deltaPlatPot(a) = PlatPot(a)- Vrest(a) ; % mV
end

% spontaneous spike rate
PulseRecPnts = PulseRecTime*sampleRate ;
spontTime = (iopb-iipe)/sampleRate ;

for a = 1:numTrials ;
    SpontSpikeRate(a) = length(spikePnt{a}(spikePnt{a}>iipe+PulseRecPnts & spikePnt{a}<iopb))/spontTime ; % spikes/sec
end

% odor response spike rate
respTime = (iore-iorb)/sampleRate ;

for a = 1:numTrials ;
    RespSpikeRate(a) = length(spikePnt{a}(spikePnt{a}>iorb+PulseRecPnts & spikePnt{a}<iore))/respTime ; % spikes/sec
end

% delta spike rate
deltaSpikeRate = RespSpikeRate - SpontSpikeRate ;

% make background times second time stamps
for a=1:length(Input(A).(id2)) ;
    for b=1:length(Input(A).(id2){a}) ;
        bgTime{a}(b) = (datenum(Input(A).(id2){a}{b})- min(tData))*24*60^2 ;
    end
end

% assess input resisitance
if (iipe-driftCheckPnts)<iipb ; % if the points you want to check are before the begining of the pulse
    disp('I pulse is not that long') ;
else
    for a = 1:numTrials ;
        Ipulse = mean(iData(a,iipe-driftCheckPnts:iipe)) - mean(iData(a,1:iipb)) ;
        Vresp = mean(vData(a,iipe-driftCheckPnts:iipe)) - mean(vData(a,1:iipb)) ;
        Rin(a) = Vresp/Ipulse ; % g ohms
    end
end



% figures
% spike detection
figure 
for a = 1:numTrials ;
    plot(time,vData(a,:)) 
    hold on
    plot(time(spikePnt{a}),vData(a,spikePnt{a}),'r*')
    title(num2str(odorRspTrials(a)))
    hold off
    pause
end

% voltage traces
figure 
subplot(5,1,1)
plot(time,vData_mean)
hold on
plot(time(iorb:iore),vData_mean(iorb:iore),'r--')

% spike raster
subplot(5,1,2:5) 
for a = 1:numTrials ;
    for b=1:length(spikePnt{a}) ;
        plot([1,1]*spikePnt{a}(b),[a-1,a])
        hold on
    end
end

% spike rate (spont and evoked)
figure 
subplot(2,1,1)
plot(tDataN,SpontSpikeRate,'*')
hold on
plot(tDataN,RespSpikeRate,'r*')
xlabel('trig time (sec)')
ylabel('spike rate (hz)')
legend('spont','evoked')

for a=1:length(bgTime) ;
    for b=1:length(bgTime{a}) ;
        plot([bgTime{a}(b),bgTime{a}(b)],[0,max([SpontSpikeRate,RespSpikeRate])],'k')
    end
end

% spike rate (change from spont to evoked)
subplot(2,1,2) 
plot(tDataN,deltaSpikeRate,'+')
xlabel('trig time (sec)')
ylabel('delta spike rate (hz)')
hold on

for a=1:length(bgTime) ;
    for b=1:length(bgTime{a}) ;
        plot([bgTime{a}(b),bgTime{a}(b)],[min(deltaSpikeRate),max(deltaSpikeRate)],'k')
    end
end

% potential (spont and evoked)
figure 
subplot(2,1,1)
plot(tDataN,Vrest,'b*')
hold on
plot(tDataN,PlatPot,'r*')
ylabel('potential (mV)')
xlabel('trig time (sec)')
legend('rest','evoked')

for a=1:length(bgTime) ;
    for b=1:length(bgTime{a}) ;
        plot([bgTime{a}(b),bgTime{a}(b)],[min([Vrest,PlatPot]),max([Vrest,PlatPot])],'k')
    end
end

% pontential (change from rest to evoked)
subplot(2,1,2) 
plot(tDataN,deltaPlatPot,'+')
ylabel(' delta plateau potential (mV)')
xlabel('trig time (sec)')
hold on

for a=1:length(bgTime) ;
    for b=1:length(bgTime{a}) ;
        plot([bgTime{a}(b),bgTime{a}(b)],[min(deltaPlatPot),max(deltaPlatPot)],'k')
    end
end

% comparing spikerate and potential (spont and evoked)
figure 
subplot(2,1,1)
plot(Vrest,SpontSpikeRate,'*')
hold on
plot(PlatPot,RespSpikeRate,'r*')
xlabel('potential (mV)')
ylabel('spike rate (hz)')

% comparing spikerate and potential (change in spont to evoked)
subplot(2,1,2) 
plot(deltaPlatPot,deltaSpikeRate,'+')
xlabel('delta potential (mV)')
ylabel('delta spike rate (hz)')

% input resistance
figure 
plot(tDataN,Rin,'b*')
ylabel('Rinput G Ohms')
xlabel('trig time (sec)')
hold on
for a=1:length(bgTime) ;
    for b=1:length(bgTime{a}) ;
        plot([bgTime{a}(b),bgTime{a}(b)],[min(Rin),max(Rin)],'k')
    end
end


