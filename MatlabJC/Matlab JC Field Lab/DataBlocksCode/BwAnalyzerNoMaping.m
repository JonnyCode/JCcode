function ForIgor = BwAnalyzerNoMaping(DataBlock, DB, Params)

% this function will analyze binary white stimuli response from array data

% JC 5/22/15
%DB=3; % TEMP

% parameters
Color_list = {'k','r','g','y','b'} ; % order of colors for each 
numDs = length(DataBlock(DB).BwPath) ; % number of data sets
RateBinT = 0.1 ; % (sec) time bin over which to calculate PSTH
StdBinT = 5 ; % (sec) time over which the std is calculated

% compare spatial and temporal receptive fields across data defined sets 
for ds = 1:numDs % for each data set

    % load data
    dataRun = load_data(DataBlock(DB).BwPath{ds}) ;
    dataRun = load_neurons(dataRun) ;

    marks_params.thresh = 3 ;
    dataRun = load_sta(dataRun) ;
    dataRun = load_params(dataRun) ;
    dataRun = get_sta_summaries(dataRun, 'all','marks_params', marks_params) ;
    dataRun = get_sta_fits_from_vision(dataRun) ;
    
    NumCells = length(dataRun.spikes) ;

    % temporal rf
    dataRun.stimulus.monitor_refresh = 60.35 ; % the monitor_refresh rate is not accurate
    trf_time(ds,:) = -[0:dataRun.stas.depth-1]*dataRun.stimulus.interval/dataRun.stimulus.monitor_refresh ; % time axis of temporal receptive field
    
    for c=1:NumCells ;
        Trf{ds}(c,:) = flipud(dataRun.stas.time_courses{c}) ;
        Trf_norm{ds}(c,:) = Trf{ds}(c,:)/norm(Trf{ds}(c,:)) ;
        
        [tcFit{ds}(c,:), final_params] = fit_time_course(fliplr(Trf{ds}(c,:))', 'verbose', false) ;
        tcFit_norm{ds}(c,:) = tcFit{ds}(c,:)/norm(tcFit{ds}(c,:)) ;
    end
        
    tcParams{ds} = tc_params_finder(fliplr(tcFit{ds}),trf_time(ds,:)) ;

    [firstPeakTimeHist{ds}, firstPeakTimeHistX{ds}] = hist(tcParams{ds}.firstPeakt,sort(trf_time(ds,:))) ;
    firstPeakTime_mean(ds) = nanmean(tcParams{ds}.firstPeakt) ;
    firstPeakTime_sem(ds) = nanstd(tcParams{ds}.firstPeakt)/sqrt(sum(~isnan(tcParams{ds}.firstPeakt))) ;

    [firstPeakHist{ds}, firstPeakHistX{ds}] = hist(tcParams{ds}.firstPeak,[-1:.01:1]) ;
    firstPeak_mean(ds) = nanmean(tcParams{ds}.firstPeak) ;
    firstPeak_sem(ds) = nanstd(tcParams{ds}.firstPeak)/sqrt(sum(~isnan(tcParams{ds}.firstPeak))) ;

    [DoTHist{ds}, DoTHistX{ds}] = hist(tcParams{ds}.DoT,[0:.25:10]) ;
    DoT_mean(ds) = nanmean(tcParams{ds}.DoT) ;
    DoT_sem(ds) = nanstd(tcParams{ds}.DoT)/sqrt(sum(~isnan(tcParams{ds}.DoT))) ; 
    
    % spike rate modulation
    PsthTime = [0:RateBinT:dataRun.duration] ;
    SpikeRate{ds} = nan(NumCells, length(PsthTime)) ;
    StdBinPnts = StdBinT/RateBinT ; % number of psth bins per std calc
    
    for c=1:NumCells ;
    	SpikeRate{ds}(c,:) = hist(dataRun.spikes{c}, PsthTime)/RateBinT ;
        SpikeRate_mean{ds}(c) = mean(SpikeRate{ds}(c,:)) ;
        
        SpikeRate_std{ds}(c) = 0 ; % preallocate
        for StdBins = [StdBinPnts:StdBinPnts:length(SpikeRate{ds}(c,:))] ;
            SpikeRate_std{ds}(c) = SpikeRate_std{ds}(c) + std(SpikeRate{ds}(c,StdBins-StdBinPnts+1:StdBins))/floor(length(SpikeRate{ds}(c,:))/StdBinPnts) ;
        end
    end   
    
    [SrStdHist{ds}, SrStdHistX{ds}] = hist(SpikeRate_std{ds},[0:1:100]) ;
    SrStd_mean(ds) = nanmean(SpikeRate_std{ds}) ;
    SrStd_sem(ds) = nanstd(SpikeRate_std{ds})/sqrt(sum(~isnan(SpikeRate_std{ds}))) ; 
end


%figures
figure
for ds = 1:numDs % for each data set
    subplot(4,1,1)
    plot(firstPeakTimeHistX{ds},firstPeakTimeHist{ds},Color_list{ds})
    xlabel('peak1 t')
    hold on

    subplot(4,1,2)
    plot(firstPeakHistX{ds},firstPeakHist{ds},Color_list{ds})
    xlabel('peak1 amp')
    hold on

    subplot(4,1,3)
    plot(DoTHistX{ds},DoTHist{ds},Color_list{ds})
    xlabel('DoT')
    hold on
    
    subplot(4,1,4)
    plot(SrStdHistX{ds},SrStdHist{ds},Color_list{ds})
    xlabel('Std spike rate')
    hold on
end

figure
for c=1:NumCells ; % for each cell
    hold off
    for ds = 1:numDs % for each data set
        plot(trf_time(ds,:),Trf_norm{ds}(c,:))
        hold on
    end
    title(num2str(c))
    pause
end
    

% for Igor
ForIgor = struct() ;
example_cell = 1 ; 

for ds = 1:numDs % for each data set  
    
    VecName = ['TrfTime','Cell',num2str(example_cell),'Ds',num2str(ds),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,trf_time(ds,:)) ; 

    VecName = ['TrfNorm','Cell',num2str(example_cell),'Ds',num2str(ds),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,Trf_norm{ds}(example_cell,:)) ; 
    
    
    VecName = ['StdHist','Ds',num2str(ds),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,SrStdHist{ds}) ; 
    
    VecName = ['StdHistX','Ds',num2str(ds),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,SrStdHistX{ds}) ; 
    
    
    VecName = ['fPeakTimeHist','Ds',num2str(ds),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,firstPeakTimeHist{ds}) ; 
    
    VecName = ['fPeakTimeHistX','Ds',num2str(ds),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,firstPeakTimeHistX{ds}) ; 
    
    
    VecName = ['fPeakAmpHist','Ds',num2str(ds),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,firstPeakHist{ds}) ; 
    
    VecName = ['fPeakAmpHistX','Ds',num2str(ds),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,firstPeakHistX{ds}) ; 
    
    
    VecName = ['DotHist','Ds',num2str(ds),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,DoTHist{ds}) ; 
    
    VecName = ['DotHistX','Ds',num2str(ds),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,DoTHistX{ds}) ; 
end







        
    
    
    


