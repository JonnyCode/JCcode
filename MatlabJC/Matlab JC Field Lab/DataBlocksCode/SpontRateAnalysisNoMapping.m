function ForIgor = SpontRateAnalysisNoMapping(DataBlock, DB, Params)

% modified from SpontRateAnalysis

% JC 7/6/15

% parameters
Color_list = {'c','k','r','g','y','b'} ; % order of colors for each 
RateBinT_long = 1 ; % sec
RateBinT_short = .001 ; % sec
MaxLag = 1 ; % sec
Ac_Trange = [.1,.2] ; % (sec) calculate mean of AC between these values 
numXrowsPlot = 5 ;
SpikeRateHistX = [0:10:50] ; % (Hz)
AcRatioHistX = [0:.01:1] ; % ratio
IsiDistX = [0:RateBinT_short:MaxLag] ;
pdfHistFlag = true ; 

% load data
for SpN = 1:length(DataBlock(DB).SpontPath) ;
    dataRun = load_data(DataBlock(DB).SpontPath{SpN}) ;
    dataRun = load_neurons(dataRun) ;
    dataRun = load_params(dataRun) ;
    dataRun = load_ei(dataRun, 'all') ; 
    
    NumCells = length(dataRun.spikes) ;
    
    % spike rate stats
    for c = 1:NumCells ; % for each cell of that type
         SpikeRate{SpN}(c,:) = hist(dataRun.spikes{c}, [0:RateBinT_long:dataRun.duration])/RateBinT_long ; % hz
         SpikeRate_mean{SpN}(c) = mean(SpikeRate{SpN}(c,:)) ;
         SpikeRate_std{SpN}(c) = mean(SpikeRate{SpN}(c,:)) ;   
    end
    SpikeRate_mean_mean(SpN) = mean(SpikeRate_mean{SpN}) ;
    SpikeRateHist{SpN} = hist(SpikeRate_mean{SpN},SpikeRateHistX) ;
    
    if pdfHistFlag ;
         SpikeRateHist{SpN} = SpikeRateHist{SpN}/sum(SpikeRateHist{SpN}) ;
    end
        
    
    % autocorrelation stats
    for c = 1:NumCells ; % for each cell of that type
        %IsiDist{SpN}{uc}(c,:) = hist(diff(dataRun.spikes{ci}),IsiDistX) ;         
         st = zeros(1,dataRun.duration * dataRun.sampling_rate) ; % make a spike train
         st(round(dataRun.spikes{c} * dataRun.sampling_rate)) = 1 ; % populate spike train
         Ac =  xcorr(st, MaxLag * dataRun.sampling_rate) ; % autocorrelation
         Ac = smooth(smooth(Ac,RateBinT_short* dataRun.sampling_rate),RateBinT_short * dataRun.sampling_rate) ; % double smoothed Ac
         Ac_smooth{SpN}(c,:) = Ac(MaxLag*dataRun.sampling_rate+1:end) ; % cut in half
         Ac_ratio{SpN}(c) = mean(Ac_smooth{SpN}(c,[Ac_Trange(1)*dataRun.sampling_rate:Ac_Trange(2)*dataRun.sampling_rate]))/Ac_smooth{SpN}(c,1) ; % peak ratio
    
%         [px,ps] = PowerSpectrumFinder(st,dataRun.sampling_rate) ; % power spectrum
%         pi(1) = find(px>=5 && px<=10) ;
        
    end
    Ac_ratio_mean(SpN) = mean(Ac_ratio{SpN}) ;
    AcRatioHist{SpN} = hist(Ac_ratio{SpN},AcRatioHistX) ;
end
AcX = [0:1/dataRun.sampling_rate:MaxLag] ; % autocorrelation X axis

% figures

figure
for c = 1:NumCells ; % for each cell of that type
    for SpN = 1:length(DataBlock(DB).SpontPath) ; 
        plot(AcX, Ac_smooth{SpN}(c,:), Color_list{SpN})
        hold on
        plot(Ac_Trange,[0,0],'r-')  
        text(.5,.5,num2str(Ac_ratio{SpN}(c)),'units','norm')
        xlim([.002,.2])
        
    end
    title(num2str(c))
    pause
    hold off
end

figure
for SpN = 1:length(DataBlock(DB).SpontPath) ;
    subplot(2,1,1)
    plot(SpikeRateHistX, SpikeRateHist{SpN}, Color_list{SpN})

    subplot(2,1,2)
    plot(AcRatioHistX, AcRatioHist{SpN}, Color_list{SpN})
end


% forIgor
Ac_ExampleCell = 17;
ForIgor = struct() ;

for SpN = 1:length(DataBlock(DB).SpontPath) ;        
    VecName = ['SpikeRateHistX','Spn',num2str(SpN),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,SpikeRateHistX) ; 

    VecName = ['SpikeRateHist','Spn',num2str(SpN),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,SpikeRateHist{SpN}) ; 

    VecName = ['AcRatioHistX','Spn',num2str(SpN),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,AcRatioHistX) ; %

    VecName = ['AcRatioHist','Spn',num2str(SpN),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,AcRatioHist{SpN}) ; %
    
    VecName = ['AcX','Cell',num2str(Ac_ExampleCell),'Spn',num2str(SpN),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,AcX) ; %

    VecName = ['Ac','Cell',num2str(Ac_ExampleCell),'Spn',num2str(SpN),'Db',num2str(DB)] ;
    ForIgor = setfield(ForIgor,VecName,Ac_smooth{SpN}(Ac_ExampleCell,:)) ; %
end





