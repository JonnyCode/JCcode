%quick ploting and basic analysis when data is from h5

% JC 1/3/12

% input
FilePath = 'C:\Cafaro Data\120314Data\120314#3\120314#3.h5' ; 
FrontEndString = 'voltage_120314_3_'  ;
Trials = [5] ;
figNum = 2 ;
sampleRate = 10000 ; %(hz)
SpikeThreshold = 10 ; % mv above base

% output logicals (true or false)
plotOneAtAtime = false ; 
averageData = false ;
MakePowerSpectrum = false ;
MakeSpikeRaster = true ;

% get data and put in cell array
DataStruct = hdf52Struct(FilePath) ; % get data out of hdf5
FieldNamesCell = FieldNameCellMaker(FrontEndString,Trials) ; % make cell with desired structure field names

for a = 1:length(Trials) ;
    Data{a} = DataStruct2Mat(DataStruct,FieldNamesCell(a)) ; % put vector data into cell array
    time{a} = [1:length(Data{a})]/sampleRate ;
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
    DataMat = cell2mat(Data') ;
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
    % spike detection
    cd Z:\MatlabJC\Matlab JC Bhandawat Lab\spike analysis
    load SpikeTemplate
    
    for a=1:length(Trials) ; 
        spikePnt = spikeFinder(Data{a},10000,.06,.001,.002) ;
        spikePnts{a} = spikePnt ;
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
        

