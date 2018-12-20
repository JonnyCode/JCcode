function ForIgor = SplitFieldAnalyzer(Input) ;

% this function will analyze data from presentation of the splitfield
% spatial stimulus (JC 10/28/09)

% get data
[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);

epochs = str2num(Input(A).SfCA) ;

for a = 1:length(epochs) ; % for each spike epoch

    [dataCA(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data

    [SICA(a), error] = ITCGetSamplingInterval(epochs(a), fp); % get sampling interval
    SICA(a) = SICA(a) * 1e-6; % Sampling interval in sec
end

if Input(A).ITC18flag == 1 ;
    SICA = SICA*1.25 ;
end

SpatialStimParams = hdf5load(['~/Data/primate/',Input(A).cellname,'_spatial.h5']) ; % load spatial stim params
frameRate = 60 ;
UnwrittenEpochs = [] ;
Unwritteni = [] ;
ZeroMeani = [] ;

for a = 1:length(epochs) ;
    StrucString = ['params_epoch',num2str(epochs(a))] ; 
    
    try
        Struct = SpatialStimParams.(StrucString) ;
        
        gammaTable = load('~/matlab/common/stim/gamma_table');
        gammaTable = gammaTable.gamma - min(gammaTable.gamma); %zero offset
        g_norm = gammaTable./max(gammaTable);
        
        intensityLeft = g_norm(Struct.intensityLeft+1) ;
        intensityRight = g_norm(Struct.intensityRight+1) ;

        ContrastLeft(a) = (intensityLeft-Struct.spatial_meanLevel)/Struct.spatial_meanLevel ;
        ContrastRight(a) = (intensityRight-Struct.spatial_meanLevel)/Struct.spatial_meanLevel ;
        DelayPnts(a) = floor((Struct.preFramesLeft-Struct.preFramesRight)/(frameRate*SICA(a))) ; % positive (delay left) negative (delay right)

        prePnts(a) = floor(Struct.spatial_prepts/(frameRate*SICA(a))) ;
        postPnts(a) = floor(Struct.spatial_postpts/(frameRate*SICA(a))) ;
        stimPnts(a) = floor(Struct.spatial_stimpts/(frameRate*SICA(a))) ;

        if Struct.preFramesLeft>=Struct.spatial_stimpts ; % if delay is longer than stim points
            ContrastLeft(a) = 0 ;
        end

        if Struct.preFramesRight>=Struct.spatial_stimpts ; % if delay is longer than stim points
            ContrastRight(a) = 0 ;
        end

        [Monitor(a,:), error] = ITCReadEpoch(epochs(a), 1, fp) ;    % get data
    catch ME
        UnwrittenEpochs = [UnwrittenEpochs,epochs(a)] ; % note epochs without spatial params recorded
        Unwritteni = [Unwritteni,a] ;
        
        ContrastLeft(a) = NaN ;
        ContrastRight(a) = NaN ;
        DelayPnts(a) = NaN ;
        prePnts(a) = NaN ;
        postPnts(a) = NaN ;
        stimPnts(a) = NaN ;
    end
    
    if Struct.spatial_meanLevel == 0 ;
        ZeroMeani = [ZeroMeani,a] ;
    end
    
end
writteni = setdiff([1:length(epochs)],[Unwritteni,ZeroMeani]) ;


% simple spike detection
SpikeTrain = zeros(size(dataCA)) ;
Threshold = -40 ;

hpData = highPassFilter(dataCA,1/SICA(1),2) ;
hpDataShift = circshift(hpData,[0,1]) ;
SpikeTrain(hpData<=Threshold & hpDataShift>=Threshold)=1 ;

% arrange spiketrain data by contrast and delay
ContrastLeftRange = unique(ContrastLeft(writteni)) ;
ContrastRightRange = unique(ContrastRight(writteni)) ;
DelayPntsRange = unique(DelayPnts(writteni)) ;

for a=1:length(DelayPntsRange) ;
    for b=1:length(ContrastLeftRange) ;
        for c=1:length(ContrastRightRange) ;
            SfArray_spikes{a}{b}{c} = SpikeTrain(DelayPnts==DelayPntsRange(a) & ContrastLeft==ContrastLeftRange(b) & ContrastRight==ContrastRightRange(c),:) ;
            sumSpikes{a}{b}{c} = sum(SfArray_spikes{a}{b}{c}) ;
            PSTHtemp{a}{b}{c} = sumSpikes{a}{b}{c}/(size(SfArray_spikes{a}{b}{c},1)*SICA(1)) ; %psth firing rate
            PSTHnorm{a}{b}{c} = PSTHtemp{a}{b}{c}/mean(PSTHtemp{a}{b}{c}(1:prePnts(1))) ; % fraction of firing rate above prestim mean
            PSTHnormSmooth1{a}{b}{c} = smooth(PSTHnorm{a}{b}{c},(1/SICA(1))*.001*1) ; % bin 1 ms
            PSTHnormSmooth10{a}{b}{c} = smooth(PSTHnorm{a}{b}{c},(1/SICA(1))*.001*10) ; % bin 10 ms
        end
    end
end


figure % strongest contrast delay responses
for a=1:length(DelayPntsRange) ;
    color=[1-(a-1)/(length(DelayPntsRange)-1) 0 (a-1)/(length(DelayPntsRange)-1)]
    if DelayPntsRange(a)>=0 ; % delay left, off left, on right
        subplot(4,2,1)
        plot(PSTHnormSmooth10{a}{1}{end},'color',color)
        hold on
        title('delay left, off left, on right')
    end
    if DelayPntsRange(a)<=0 ; % delay right, off left, on right
        subplot(4,2,2)
        plot(PSTHnormSmooth10{a}{1}{end},'color',color)
        hold on
        title('delay right, off left, on right')
    end
    if DelayPntsRange(a)>=0 ; % delay left, off left, off right 
        subplot(4,2,3)
        plot(PSTHnormSmooth10{a}{1}{1},'color',color)
        hold on
        title('delay left, off left, off right ')
    end
    if DelayPntsRange(a)<=0 ; % delay right, off left, off right 
        subplot(4,2,4)
        plot(PSTHnormSmooth10{a}{1}{1},'color',color)
        hold on
        title('delay right, off left, off right')
    end
    if DelayPntsRange(a)>=0 ; % delay left, on left, on right 
        subplot(4,2,5)
        plot(PSTHnormSmooth10{a}{end}{end},'color',color)
        hold on
        title('delay left, on left, on right')
    end
    if DelayPntsRange(a)<=0 ; % delay right, on left, on right 
        subplot(4,2,6)
        plot(PSTHnormSmooth10{a}{end}{end},'color',color)
        hold on
        title('delay right, on left, on right')
    end    
    if DelayPntsRange(a)>=0 ; % delay left, on left, off right 
        subplot(4,2,7)
        plot(PSTHnormSmooth10{a}{end}{1},'color',color)
        hold on
        title('delay left, on left, off right')
    end
    if DelayPntsRange(a)<=0 ; % delay right, on left, off right 
        subplot(4,2,8)
        plot(PSTHnormSmooth10{a}{end}{1},'color',color)
        hold on
        title('delay right, on left, off right')
    end
            
end

figure % delay = 0 ;
i = find(DelayPntsRange==min(abs(DelayPntsRange))) ;
for a=1:length(ContrastLeftRange)/2 + 1 ;
    color=[1-(a-1)/(length(ContrastLeftRange)-1) 0 (a-1)/(length(ContrastLeftRange)-1)]
    
    subplot(2,2,1) % off left on right
    plot(PSTHnormSmooth10{i}{a}{length(ContrastLeftRange)-a+1},'color',color)
    hold on
    title('off left on right')
    
    subplot(2,2,2) % off left off right
    plot(PSTHnormSmooth10{i}{a}{a},'color',color)
    hold on
    title('off left off right')
    
    subplot(2,2,3) % on left on right
    plot(PSTHnormSmooth10{i}{length(ContrastLeftRange)-a+1}{length(ContrastLeftRange)-a+1},'color',color)
    hold on
    title('on left on right') 
    
    subplot(2,2,4) % on left off right
    plot(PSTHnormSmooth10{i}{length(ContrastLeftRange)-a+1}{a},'color',color)
    hold on
    title('on left off right')
end
    
    
figure % comparing one side on and the other side on, same, or off
for a=1:floor(length(ContrastLeftRange)/2)
    subplot(floor(length(ContrastLeftRange)/2),2,a*2-1) % on left 
    plot(PSTHnormSmooth10{i}{ceil(length(ContrastLeftRange)/2)+a}{ceil(length(ContrastLeftRange)/2)+a},'b') ; % right on
    hold on
    plot(PSTHnormSmooth10{i}{ceil(length(ContrastLeftRange)/2)+a}{ceil(length(ContrastLeftRange)/2)},'g') ; % right no change
    plot(PSTHnormSmooth10{i}{ceil(length(ContrastLeftRange)/2)+a}{ceil(length(ContrastLeftRange)/2)-a},'r') ; % right off

    subplot(floor(length(ContrastLeftRange)/2),2,a*2) % on right
    plot(PSTHnormSmooth10{i}{ceil(length(ContrastLeftRange)/2)+a}{ceil(length(ContrastLeftRange)/2)+a},'b') ; % left on
    hold on
    plot(PSTHnormSmooth10{i}{ceil(length(ContrastLeftRange)/2)}{ceil(length(ContrastLeftRange)/2)+a},'g') ; % left no change
    plot(PSTHnormSmooth10{i}{ceil(length(ContrastLeftRange)/2)-a}{ceil(length(ContrastLeftRange)/2)+a},'r') ; % left off
end
    




