function ForIgor = ORNanalysisB(Input,id,A) 

% hard coded parameters
numClusters = 2 ;
numPCs = 2 ;
NegDiffThreshStd = 1 ; 
sampleRate = 10000 ; %(hz)
SpikeThreshold = 10 ; % mv above base
spontTime = [0,13] ; % sec during which spontaneous rate should be assessed 
psthSmthTime = .1 ; %sec

% time to pnts
spontPnts = spontTime*sampleRate ;

% input
FilePath = ['C:\Cafaro Data\',Input(A).cellname(1:end-2),'DataSH\',Input(A).cellname,'\',Input(A).cellname,'.h5'] ; 
FrontEndString = ['voltage_',Input(A).cellname(1:end-2),'_',Input(A).cellname(end),'_']  ;


% get data and put in cell array
DataStruct = hdf52Struct(FilePath) ; % get data out of hdf5

Sweeps = str2num(Input(A).OdorRsp) ;
FieldNamesCell = FieldNameCellMaker(FrontEndString,Sweeps) ; % make cell with desired structure field names
for a = 1:length(Sweeps) ;
    Data{a} = DataStruct2Mat(DataStruct,FieldNamesCell(a)) ; % put vector data into cell array
end

% time
time = [1:length(Data{a})]/sampleRate ;

% detect and seperate spikes
for a = 1:length(Sweeps) ;
    [spikePntGroup] = spikeSorter(Data{a}, sampleRate, NegDiffThreshStd, numClusters, numPCs, false) ;
%     pause
%     close all
    for b = 1:numClusters ;
        spikePnts{a,b} = spikePntGroup{b} ;
    end
end

% spike trains and psth
for a = 1:numClusters ;
    for b = 1:length(Sweeps) ;
        SpikeTrain{a}(b,:) = zeros(1,length(Data{b})) ;
        SpikeTrain{a}(b,spikePnts{b,a}) = 1 ;
    end
    psth(a,:) = smooth(mean(SpikeTrain{a},1),psthSmthTime*sampleRate)*sampleRate ;
end
psth(:,[1:100,end-100:end]) = 0 ;

% slow pot change and odor response time pnts
for a = 1:length(Sweeps) ;
    temp(a,:) = lowPassFilter(Data{a},sampleRate,1) ;
end
DataLP = mean(temp,1) ;
DataLPthresh = mean(DataLP)-2*std(DataLP) ;
odorPnts(1) = find(DataLP<DataLPthresh,1,'first') ;
odorPnts(2) = find(DataLP<DataLPthresh,1,'last') ;

% spont spike number and rate
for a = 1:length(Sweeps) ;
    for b = 1:numClusters ; 
        spikeNum(a,b) = sum((spikePnts{a,b}>=spontPnts(1)).* (spikePnts{a,b}<=spontPnts(2))) ;
        SpontRate(a,b) = spikeNum(a,b)/(spontTime(2)-spontTime(1)) ;
    end
end
SpontRateMean = mean(SpontRate,1) ;
SpontRateStd = std(SpontRate,0,1) ;

% odor response number and rate
for a = 1:length(Sweeps) ;
    for b = 1:numClusters ; 
        spikeNum(a,b) = sum((spikePnts{a,b}>=odorPnts(1)).*(spikePnts{a,b}<=odorPnts(2))) ;
        OdorRate(a,b) = spikeNum(a,b)/((odorPnts(2)-odorPnts(1))/sampleRate) ;
    end
end
OdorRateMean = mean(OdorRate,1) ;
OdorRateStd = std(OdorRate,0,1) ;


% figures
colorVec = ['y','r']

% raw data and spike detection
figure
for a = 1:length(Sweeps) ;
    plot(time,Data{a})
    hold on
    for b=1:numClusters ;
        plot(time(spikePnts{a,b}),Data{a}(spikePnts{a,b}),[colorVec(b),'*'])
    end
    hold off
    pause
end

% spike data
figure 

% low pass data
subplot(numClusters*2+2,1,1)
plot(time,DataLP)
hold on
plot(time(odorPnts(1):odorPnts(2)),DataLP(odorPnts(1):odorPnts(2)),'r--')

% psth
subplot(numClusters*2+2,1,2)
plot(time,psth)

% spike raster
for a = 1:numClusters ;
    subplot(numClusters*2+2,1,((a-1)*numClusters)+[3:4]) 
    for b = 1:length(Sweeps) ;
        for c=1:length(spikePnts{b,a}) ;
            plot((1/sampleRate)*[1,1]*spikePnts{b,a}(c),[b-1,b])
            hold on
        end
    end
end

% for igor
identifier = ['SpontSpikeRate',num2str(A)] ;
ForIgor.(identifier) = SpontRateMean ;

identifier = ['SpontSpikeRateStd',num2str(A)] ;
ForIgor.(identifier) = SpontRateStd ;

identifier = ['OdorSpikeRate',num2str(A)] ;
ForIgor.(identifier) = OdorRateMean ;

identifier = ['OdorSpikeRateStd',num2str(A)] ;
ForIgor.(identifier) = OdorRateStd ;


