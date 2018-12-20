%% Quick plot
cellname = '030408Ac1' ;
epochs_str = {'[7:33]'} ; 
Amp = 0 ;
Color = {'b','r','g','k','y','c','r','g','k','y','b-.','r-.','g-.','k-.','y-.','k:'} ;
fig_num = 1;


average = 0 ; % average traces
spikePlot = 1 ; % spike raster and psth
threshold = 15 ;
spatial = 0 ;
dynamicClamp = 0 ; % get conductances injected
LEDstimulus = 1 ; % get stimulus

try
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',cellname]);
catch
    [fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',cellname]);
end

for a=1:length(epochs_str)
    epochs = str2num(epochs_str{a}) ;
    for b=1:length(epochs)
        [data{a}(b,:), error] = ITCReadEpoch(epochs(b), Amp, fp);
        [SI, error] = ITCGetSamplingInterval(epochs(b), fp);   
    
        if dynamicClamp==1;
            [excg{a}(b,:), inhg{a}(b,:), error] = ITCReadEpochStmGClamp(epochs(b), 0, fp);
        end
        
        if LEDstimulus==1;
             [stimulus{a}(b,:),error] = ITCReadEpochStm(epochs(b),0,fp) ;
        end
        
    end
    
    SI = SI*10^-6 ;
    time{a} = [1:length(data{a}(1,:))]*SI ;
    
    if spatial==1 ;
        try
            SpatialStimParams = hdf5load(['~/Data/mouse/',cellname,'_spatial.h5']) ;
        catch
            SpatialStimParams = hdf5load(['~/Data/primate/',cellname,'_spatial.h5']) ;
        end
        
        StrucString = ['params_epoch',num2str(epochs(1))] ;
        SpatialStimParams.(StrucString)       
    end       
end

figure(fig_num)
for a=1:length(epochs_str) 
    if average==1;
        plot(time{a},mean(data{a},1),Color{a})
    else
        plot(time{a},data{a},Color{a})
    end
    hold on
end

if dynamicClamp==1 ;
    figure
    for a=1:length(epochs_str) 
        if average==1;
            plot(time{a},mean(excg{a},1),Color{a})
            hold on
            plot(time{a},mean(inhg{a},1),[Color{a},'--'])
        else
            plot(time{a},excg{a},Color{a})
            hold on
            plot(time{a},inhg{a},[Color{a},'--'])
        end
        hold on
    end
end

if spikePlot==1 ;
    round = 0 ;
    for a=1:length(data) ;
        for b=1:size(data{a},1) ;
            round = round + 1 ;
            datahpf = highPassFilter(data{a}(b,:), 1/SI, 1) ;
             
            datahpfShift = circshift(datahpf,[0,1]) ; % shift
            SpikeTrain = nan(size(datahpf)) ; % prep spike trains
            if threshold>0 ;
                SpikeTrain(datahpf>=threshold & datahpfShift<threshold) = 0 ; % find indicies where change in current is above threshold
            else
                SpikeTrain(datahpf<=threshold & datahpfShift>threshold) = 0 ;
            end 
            SpikeTrain(1) = nan ; % cannot dectect first point as spike
            Train{a}(b,:) = SpikeTrain ;
            SpikeTimes{a}{b} = time{a}(Train{a}(b,:)==0) ;
            
            figure(fig_num+1)
            plot(datahpf)
            hold on
            plot(SpikeTrain+threshold,'r*')
            pause
            hold off
            
            clear SpikeTrain

            figure(fig_num)
            plot(time{a},Train{a}(b,:)+round,Color{a})
            hold on
        end     
    end
    set(gca,'ylim',[0 round+1])
    
    for a=1:length(SpikeTimes) ;
        figure
        for b=1:length(SpikeTimes{a}) ;
            for c=1:length(SpikeTimes{a}{b}) ;
                plot([SpikeTimes{a}{b}(c),SpikeTimes{a}{b}(c)],[b-1,b],Color{a})
                hold on
            end
        end
    end
    
    
end

%%
try
    SpatialStimParams = hdf5load(['~/Data/mouse/',cellname,'_spatial.h5']) ;
catch
    SpatialStimParams = hdf5load(['~/Data/primate/',cellname,'_spatial.h5']) ;
end


StrucString = ['params_epoch',] ;
SpatialStimParams.(StrucString)

            
        

            
            
            
            
    
    


