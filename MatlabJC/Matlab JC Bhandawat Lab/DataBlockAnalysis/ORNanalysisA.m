function ForIgor = ORNanalysisA(Input,id,A) 

numClusters = 1 ;
numPCs = 2 ;
NegDiffThreshStd = 2 ; 

% input
FilePath = ['C:\Cafaro Data\',Input(A).cellname(1:end-2),'Data\',Input(A).cellname,'\',Input(A).cellname,'.h5'] ; 
FrontEndString = ['voltage_',Input(A).cellname(1:end-2),'_',Input(A).cellname(end),'_']  ;

sampleRate = 10000 ; %(hz)
SpikeThreshold = 10 ; % mv above base

% get data and put in cell array
DataStruct = hdf52Struct(FilePath) ; % get data out of hdf5

SpontSweeps = str2num(Input(A).Spont) ;
FieldNamesCell = FieldNameCellMaker(FrontEndString,SpontSweeps) ; % make cell with desired structure field names
for a = 1:length(SpontSweeps) ;
    SpontData{a} = DataStruct2Mat(DataStruct,FieldNamesCell(a)) ; % put vector data into cell array
    SpontTime{a} = [1:length(SpontData{a})]/sampleRate ;
end

% detect and seperate spikes
for a = 1:length(SpontSweeps) ;
    [spikePntGroup] = spikeSorter(SpontData{a}, sampleRate, NegDiffThreshStd, numClusters, numPCs, false) ;
%     pause
%     close all
    for b = 1:numClusters ;
        spikePnts{a,b} = spikePntGroup{b} ;
    end
end

% spike number and rate
for a = 1:length(SpontSweeps) ;
    for b = 1:numClusters ; 
        spikeNum(a,b) = length(spikePnts{a,b}) ;
        SpontRate(a,b) = spikeNum(a,b)/SpontTime{a}(end) ;
    end
end

SpontRateMean = mean(SpontRate,1) ;
SpontRateStd = std(SpontRate,1) ;

% for igor
identifier = ['SpontSpikeRate',num2str(A)] ;
ForIgor.(identifier) = SpontRateMean ;

identifier = ['SpontSpikeRateStd',num2str(A)] ;
ForIgor.(identifier) = SpontRateStd ;

