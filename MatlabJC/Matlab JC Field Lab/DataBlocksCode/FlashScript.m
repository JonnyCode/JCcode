% script to analyze flash data and help me better understand data
% structure.

% JC 11/6/14

% data block

% cell from Greg - some trials have longer delays between flashblocks
Data(1).DataPath = '/Volumes/Berlioz/Analysis/2014-09-23-0/data000/data000' ;
Data(1).flashBlocks = {[1:21],[22:42],[43:63],[64:84],[85:105],[106:126],...
    [127:148],[149:169],[170:191],[192:212],[213:215],[216:236],[237:257],...
    [258:277],[278:297],[298:317],[318:337],[338:348]} ;
Data(1).flashDuration = [2,8,4,2,2,4,4,8,8,2,2,2,4,8,8,4,2,2] ;
Data(1).NDF = [5,5,5,4,5,5,4,4,5,4,4,3,4,4,3,3,2,1] ;

% load data
datarun = load_data(Data(1).DataPath);
datarun = load_neurons(datarun) ;

% psth for flash data
DB = 1; % data block
timeWindow = 2 ; % seconds
timeBin = .1 ; % seconds
time = [-timeWindow:timeBin:timeWindow] ;

for fb = 1:length(Data(DB).flashBlocks) ; % for each flash block

    for a = 1:length(datarun.spikes) ; % for each cell
        for b = Data(1).flashBlocks{fb} ; % for each flash
            Temp{b} = datarun.spikes{a}(datarun.spikes{a}>datarun.triggers(b)-timeWindow... 
                & datarun.spikes{a}<(datarun.triggers(b)+timeWindow))-datarun.triggers(b) ;
        end
        SpikesMat = cell2mat(Temp') ;
        histogram(a,:) = hist(SpikesMat,time)/length(Data(1).flashBlocks{fb}) ;
        clear Temp SpikesMat
    end
    psth_mean(fb,:) = mean(histogram/timeBin) ;
    psth_std(fb,:) = std(histogram/timeBin) ;
    clear histogram 

end

figure
for fb = 1:length(Data(DB).flashBlocks) ; % for each flash block
    subplot(3,6,fb)
    plot(time,psth_mean(fb,:),'k')
    hold on
    plot(time,psth_mean(fb,:)+psth_std(fb,:),'y:')
    plot(time,psth_mean(fb,:)-psth_std(fb,:),'y:')
    ylim([0,60])
    title(['NDF=',num2str(Data(1).NDF(fb)),',flash=',num2str(Data(1).flashDuration(fb)),'ms'])
end
