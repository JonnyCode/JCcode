% quick ploting and basic analysis

% JC 4/13/12

% input
Cellname = '130606_1' ;
Trials = [94:113] ;
figNum = 1 ;
sampleRate = 10000 ; %(hz)
spikeDetectionParameters.absRefTime = .002 ; % (sec)
spikeDetectionParameters.minRise = .4 ; % (mv)
spikeDetectionParameters.minFall = .2 ; % (mv)
spikeDetectionParameters.filterOrder = 1
spikeDetectionParameters.lpfCutOff = 4000 ;
spikeDetectionParameters.minRiseTime = 0.001 ;
spikeDetectionParameters.minFallTime = 0.001 ;
spikeDetectionParameters.maxPlatTime = 0.001 ;
spikeDetectionParameters.NegDiffThreshStd = 1.5 ;
PsthBinTime = .5 ; % (sec)
dataType = 'vData' ;
spikeDetectionType = 'spikeSorter' ;

% output logicals (true or false)
plotOneAtAtime = true ; 
averageData = false ;
MakePowerSpectrum = false ;
MakeSpikeRaster = false ;
MakePsth = true ;
plotAo0Data = false ;

% get data and put in cell array
rootdir = ['Z:\Cafaro Data Backup\', Cellname(1:6),'Data'];

for a = 1:length(Trials) ;
    temp = load([rootdir,'\',Cellname,'\','voltage_',Cellname,'_',num2str(Trials(a))]) ;
    vData{a} = temp.voltage ;
    
    temp = load([rootdir,'\',Cellname,'\','current_',Cellname,'_',num2str(Trials(a))]) ;
    iData{a} = temp.current ;
    
    temp = load([rootdir,'\',Cellname,'\','Ao0_',Cellname,'_',num2str(Trials(a))]) ;
    ao0Data{a} = temp.Ao0 ; % odor
    
    temp = load([rootdir,'\',Cellname,'\','Ao1_',Cellname,'_',num2str(Trials(a))]) ;
    ao1Data{a} = temp.Ao1 ; % vext   
    
    temp = load([rootdir,'\',Cellname,'\','TrigTime_',Cellname,'_',num2str(Trials(a))]) ;
    tData(a) = temp.Trigtime ;
end

if strcmp(dataType,'vData') ;
    Data = vData ;
    for a = 1:length(Data) ;
        time{a} = [1:length(Data{a})]/sampleRate ;
    end
end

% plot data
figure 
for a = 1:length(Data) ;
    plot(time{a},Data{a})
    
    if plotOneAtAtime ;
        title(num2str(Trials(a)))
        hold off ;
        pause ;
    else
        title('all data')
        hold on ;
    end
end

% average data
if averageData ;
    DataMat = cell2mat(Data)' ;
    Data_mean = mean(DataMat) ;
    
    figure
    plot(time{1},Data_mean) ;
end
    
% power spectrum
if MakePowerSpectrum ;
    for a=1:length(Trials) ;
        [powerspec_xvalues(a,:), powerspec(a,:)] = PowerSpectrumFinder(Data{a},sampleRate) ;
    end
    
    figure
    plot(powerspec_xvalues(1,:),mean(powerspec,1)) ;
    h = gca ;
    set(h,'xscale','log')
    set(h,'yscale','log')
end

% spike rasters
if MakeSpikeRaster ; 
    DataMat = cell2mat(Data)' ;
    
    % spike detection
    if strcmp(spikeDetectionType, 'spikeFinder') ; 
        [SpikePnts,SpikeData,NonSpikeData] = spikeFinder(DataMat,sampleRate,absRefTime,minRise,minFall) ;
    elseif strcmp(spikeDetectionType, 'spikeSorter') ;
        for a=1:length(Trials) ;
            [spikePntGroup] = spikeSorter(DataMat(a,:), sampleRate, spikeDetectionParameters.NegDiffThreshStd, 1, 2, false) ;
            SpikePnts{a}= spikePntGroup{1} ;
        end
    end
            
        
    % spike detection check
    figure 
    for a=1:length(Trials) ;
        title(num2str(Trials(a)))
        plot(time{a},Data{a})
        hold on
        plot(time{a}(SpikePnts{a}),Data{a}(SpikePnts{a}),'r*')
        hold off
        pause
    end
        
    % spike raster
    for a=1:length(Trials) ;
        for b=1:length(SpikePnts{a}) ;
            plot([time{a}(SpikePnts{a}(b)),time{a}(SpikePnts{a}(b))],[a-1,a])
            hold on
        end
    end
    
end
      
% spike rasters
if MakePsth ; 
    DataMat = cell2mat(Data)' ;
    
    % spike detection
    if strcmp(spikeDetectionType, 'spikeFinder') ; 
        [SpikePnts,SpikeData,NonSpikeData] = spikeFinder(DataMat,sampleRate,absRefTime,minRise,minFall) ;
    elseif strcmp(spikeDetectionType, 'spikeSorter') ;
        for a=1:length(Trials) ;
            [spikePntGroup] = spikeSorter(DataMat(a,:), sampleRate, spikeDetectionParameters.NegDiffThreshStd, 1, 2, false) ;
            SpikePnts{a}= spikePntGroup{1} ;
        end
    end
            

    % spike detection check
    figure 
    for a=1:length(Trials) ;
        title(num2str(Trials(a)))
        plot(time{a},Data{a})
        hold on
        plot(time{a}(SpikePnts{a}),Data{a}(SpikePnts{a}),'r*')
        hold off
        pause
    end
    
    % spike trains and psth
    PsthBinPnts = PsthBinTime*sampleRate ;
    SpikeTrain = zeros(length(Trials),length(Data{1})) ;
    for a=1:length(Trials) ;
        SpikeTrain(a,SpikePnts{a}) = 1 ;
        SpikeTrain(a,:) = smooth(SpikeTrain(a,:),PsthBinPnts) ;
    end
    Psth = mean(SpikeTrain,1)*sampleRate ;
    
    plot(time{1},Psth)
    
end

% plot data and Ao0 together
if plotAo0Data == true ;
    figure 
    for a = 1:length(Data) ;
        plotyy(time{a},Data{a},time{a},ao0Data{a})
        hold on
        plot(time{a},Data{a})

        if plotOneAtAtime ;
            title(num2str(Trials(a)))
            hold off ;
            pause ;
        else
            title('all data')
            hold on ;
        end
    end
end










