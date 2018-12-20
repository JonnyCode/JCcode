function ForIgor = ValveNoiseAnalysisB(Input,A)

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

id = 'ValveNoise' ;

% load data in matricies
rootdir = ['Z:\Cafaro Data Backup\', Input(A).cellname(1:6),'Data'];

odorRspTrials = str2num(Input(A).(id)) ;
numTrials = length(odorRspTrials) ; 

for a = 1:numTrials ;
    temp = load([rootdir,'\',Input(A).cellname,'\','voltage_',Input(A).cellname,'_',num2str(odorRspTrials(a))]) ;
    vData{a} = temp.voltage ;
    
    temp = load([rootdir,'\',Input(A).cellname,'\','current_',Input(A).cellname,'_',num2str(odorRspTrials(a))]) ;
    iData{a} = temp.current ;
    
    temp = load([rootdir,'\',Input(A).cellname,'\','Ao0_',Input(A).cellname,'_',num2str(odorRspTrials(a))]) ;
    ao0Data{a} = temp.Ao0 ; % odor
    
    temp = load([rootdir,'\',Input(A).cellname,'\','Ao1_',Input(A).cellname,'_',num2str(odorRspTrials(a))]) ;
    ao1Data{a} = temp.Ao1 ; % vext   
    
    temp = load([rootdir,'\',Input(A).cellname,'\','TrigTime_',Input(A).cellname,'_',num2str(odorRspTrials(a))]) ;
    tData(a) = temp.Trigtime ;
end

% time vector
for a = 1:numTrials ;
    time{a} = [1:size(vData{a})]/sampleRate ;
end

tDataN = (tData - min(tData))*24*60^2 ; % convert to seconds since experiment began

% make odor valve pulse binary
ao0DataB = ao0Data ; 
for a = 1:numTrials ;
    ao0DataB{a}(ao0Data{a}>=5) = 1 ; 
    ao0DataB{a}(ao0Data{a}<5) = 0 ;
end

% round Vext pulse to nearest 10 mV and get rid of single sample pulses
for a = 1:numTrials ;    
    ao1DataR{a} = round(ao1Data{a}*100)/100 ; 
    for b=1:length(ao1DataR{a})-2 ;
        if ao1DataR{a}(b)~=ao1DataR{a}(b+1) && ao1DataR{a}(b)==ao1DataR{a}(b+2) ;
            ao1DataR{a}(b+1)= ao1DataR{a}(b) ;
        end
    end
end

% same noise seed or different?
ao0DataBdiff = ao0DataB{1}- ao0DataB{2} ;
if sum(abs(ao0DataBdiff(:)))~=0 ;
    disp('noise different seeds') ;
end

% was there a current pulse and was it at the same time always
for a=1:numTrials ;
    pulsei{a} = find(ao1DataR{a}~=ao1DataR{a}(1)) ;
    pulsei_length(a) = length(pulsei{a}) ;
end
    
if sum(pulsei_length)~=0 ;
    disp('Current is changed')
    if sum([pulsei_length - pulsei_length(1)])==0 ;
        if sum(sum(cell2mat(pulsei)-repmat(pulsei{1},numTrials,1))) ;
            disp('pulses at different times') ;
        end
    else
        disp('pulses at different times') ;
    end
end

% detect spikes in voltage data
for a=1:numTrials ;
    temp = spikeFinder(vData{a}',sampleRate,absRefTime,minRise,minFall) ;
    spikePnt{a} = temp{1} ;
end

% spike trains

for a = 1:numTrials ;
    spikeTrains{a} = zeros(size(vData{a}')) ;
    spikeTrains{a}(spikePnt{a}) = 1 ;
end

% sta
StaDataBin = 50000 ;
staLength = 7000 ;
for a = 1:numTrials ;
    sta{a} = staFinder(ao0DataB{a}'-mean(ao0DataB{a}'),spikeTrains{a},staLength) ;
end
    
% sta over time bins
for a = 1:numTrials ;
    NumSta(a) = length(spikeTrains{a})/StaDataBin ;
    for b=1:NumSta(a) ;
        staBlocks{a}(b,:) = staFinder(ao0DataB{a}(b*StaDataBin-StaDataBin+1:b*StaDataBin)'-mean(ao0DataB{a}(b*StaDataBin-StaDataBin+1:b*StaDataBin)'),...
            spikeTrains{a}(b*StaDataBin-StaDataBin+1:b*StaDataBin),staLength) ;
    end
end
    

% figures

figure % current pulses
for a = 1:numTrials ;
    plot(time{a},ao1DataR{a})
    hold on
end

figure % spike detection
for a = 1:numTrials ;
    plot(time{a},vData{a}) 
    hold on
    plot(time{a}(spikePnt{a}),vData{a}(spikePnt{a}),'r*')
    hold off
    pause
end

figure % spike rastor and odor stim
for a=1:numTrials ;
    for b=1:length(spikePnt{a}) ;
        plot([1,1]*time{a}(spikePnt{a}(b)),[a-.5,a],'r')
        hold on
    end
    plot(time{a},ao0DataB{a}/2+a-1,'k')  
end

figure % STAs
for a=1:numTrials ;
    plot(sta{a})
    hold on
end

figure
for a=1:numTrials ;
    plot(sta{a},'c')
    hold on
    
    plot([1:staLength],staBlocks{a})
    hold off
    pause
end

figure % STAs

for a=1:numTrials ;
    subplot(numTrials,max(NumSta)+1,a*(max(NumSta)+1)-(max(NumSta)+1)+1) ;
    plot([1:length(sta{a})]/10000,sta{a})
    axis tight
    ylim([-1,1])
    
    for b=1:NumSta(a) ;
        subplot(numTrials,max(NumSta)+1,a*(max(NumSta)+1)-(max(NumSta)+1)+1+b) ;
        plot([1:length(sta{a})]/10000,staBlocks{a}(b,:))
        axis tight
        ylim([-1,1])
    end
   
end


