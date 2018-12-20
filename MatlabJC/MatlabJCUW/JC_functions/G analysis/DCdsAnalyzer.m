function ForIgor = DCdsAnalyzer(Input,Parameters,id,id2,A) ;

% this will analyze dynamic clamp data from on-off directional selective cells 
%injected with sim, shuffled, or sim repeated g obtained during moving bar.

% JC 11/12/09
% edditted 8/13/11 to include dynamic clamp through hardware flag


epochs = str2num(Input(A).(id)) ; %sim and shuffled g
epochsSimRepeat = str2num(Input(A).(id2)) ; % sim repeated

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/mouse/',Input(A).cellname]);

for a = 1:length(epochs) ;
    if Input(A).gHWflag ==1 ;
        [Data(a,:), error] = ITCReadEpoch(epochs(a), 2, fp) ; 
        [excg(a,:), inhg(a,:), error(a,:)] = ITCReadEpochStmGClamp(epochs(a), 3, fp);
    else 
        [Data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ; 
        [excg(a,:), inhg(a,:), error(a,:)] = ITCReadEpochStmGClamp(epochs(a), 0, fp);
    end
    
    I(a,:) = (Data(a,:)+80).*(inhg(a,:)) + Data(a,:).*excg(a,:) ; % added this line 8/16/10 to calculate current
end

for a = 1:length(epochsSimRepeat) ;
    if Input(A).gHWflag ==1 ;
        [dataSimRepeat(a,:), error] = ITCReadEpoch(epochsSimRepeat(a), 2, fp) ; 
        [excgSimRepeat(a,:), inhgSimRepeat(a,:), error(a,:)] = ITCReadEpochStmGClamp(epochsSimRepeat(a), 3, fp);
    else
        [dataSimRepeat(a,:), error] = ITCReadEpoch(epochsSimRepeat(a), 0, fp) ; 
        [excgSimRepeat(a,:), inhgSimRepeat(a,:), error(a,:)] = ITCReadEpochStmGClamp(epochsSimRepeat(a), 0, fp);
    end
end

[SI, error] = ITCGetSamplingInterval(epochs(1), fp); % get sampling interval    
SI = SI * 1e-6; % Sampling interval in sec
if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end
time = [SI:SI:SI*length(Data)] ;

SpikePnts = SpikeDetection_WC(Data,-30,10000) ;
SpikePntsSimRepeat = SpikeDetection_WC(dataSimRepeat,-30,10000) ;

spikeTrainSimRepeat = zeros(length(SpikePntsSimRepeat),length(Data)) ;
spikeTrainShuff = spikeTrainSimRepeat ;
spikeTrainSim = spikeTrainSimRepeat ;


round = 0 ;
for a = 1:2:length(epochs) ; % for each spike epoch
    round = round +1 ;

    dataSim(round,:) = Data(a,:) ;       
    dataShuff(round,:) = Data(a+1,:) ;       

    spikeTrainSim(round,SpikePnts{a}) = 1 ;
    spikeTrainShuff(round,SpikePnts{a+1}) = 1 ;

    spikeTimes_Sim{round} = time(SpikePnts{a}) ;
    spikeTimes_Shuff{round} = time(SpikePnts{a+1}) ;

    excgSim(round,:) = excg(a,:) ;
    inhgSim(round,:) = inhg(a,:) ;    
  
    excgShuff(round,:) = excg(a+1,:) ;
    inhgShuff(round,:) = inhg(a+1,:) ;
    
end

for a = 1:length(epochsSimRepeat)
    spikeTrainSimRepeat(a,SpikePntsSimRepeat{a}) = 1 ; 
    spikeTimes_SimRepeat{a} = time(SpikePntsSimRepeat{a}) ;
end

excgSim_Mean = mean(excgSim) ;
inhgSim_Mean = mean(inhgSim) ;

excgShuff_Mean = mean(excgShuff) ;
inhgShuff_Mean = mean(inhgShuff) ;

excgSimRepeat_Mean = mean(excgSimRepeat) ; %mean conductances
inhgSimRepeat_Mean = mean(inhgSimRepeat) ;


for a = 1:length(epochs)/2 ; % for each set of g
    
    excgSim_Res(a,:) = excgSim(a,:) - excgSim_Mean ;
    inhgSim_Res(a,:) = inhgSim(a,:) - inhgSim_Mean ;
    
    excgShuff_Res(a,:) = excgShuff(a,:) - excgShuff_Mean ;
    inhgShuff_Res(a,:) = inhgShuff(a,:) - inhgShuff_Mean ;

    ccSim(a,:) = xcov(excgSim_Res(a,:),inhgSim_Res(a,:),'coef') ;
    ccShuff(a,:) = xcov(excgShuff_Res(a,:),inhgShuff_Res(a,:),'coef') ;
end

ccSim_Mean = mean(ccSim) ;
ccShuff_Mean = mean(ccShuff) ;

corrCoefSim = ccSim_Mean(length(excgSim)) ;
corrCoefShuff = ccShuff_Mean(length(excgSim)) ;

ResVarexcgSim_Mean = mean(var(excgSim_Res,[],2)) ;
ResVarinhgSim_Mean = mean(var(inhgSim_Res,[],2)) ;

ResVarexcgShuff_Mean = mean(var(excgShuff_Res,[],2)) ;
ResVarinhgShuff_Mean = mean(var(inhgShuff_Res,[],2)) ;

time_cc = [SI(1)*([1:length(ccSim_Mean)] - (length(ccSim_Mean)+1)/2)] ;


% tuning curve with error bars
%groupTime = 3.0667 ; % (sec) time of each bar %OLD
groupTime = 3.0388 ; % 8/13/11 after g fixed
OnTime = .9 ; % (sec) time of on response

%prePnts = 2000 ; % pnts added to g before g stim %OLD
prePnts = 0 ; % 8/13/11 after g fixed
groupPnts = floor(groupTime/SI(1));
OnPnts = floor(OnTime/SI(1)) ;

for trial = 1:length(epochs)/2 ; % for each trial
    for a=1:8 ; %for each bar
        SpikeRateSim(trial,a) = sum(spikeTrainSim(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a))/groupTime ; % hz
        SpikeRateSim_ON(trial,a) = sum(spikeTrainSim(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a-(groupPnts-OnPnts)))/OnTime ; % hz
        SpikeRateSim_OFF(trial,a) = sum(spikeTrainSim(trial,prePnts+groupPnts*a-groupPnts+OnPnts+1:prePnts+groupPnts*a))/(groupTime-OnTime) ; % hz

        SpikeRateShuff(trial,a) = sum(spikeTrainShuff(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a))/groupTime ; % hz
        SpikeRateShuff_ON(trial,a) = sum(spikeTrainShuff(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a-(groupPnts-OnPnts)))/OnTime ; % hz
        SpikeRateShuff_OFF(trial,a) = sum(spikeTrainShuff(trial,prePnts+groupPnts*a-groupPnts+OnPnts+1:prePnts+groupPnts*a))/(groupTime-OnTime) ; % hz
        
        temp = xcov(excgSim_Res(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a),inhgSim_Res(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a),'coef') ;
        [m,peaki] = max(abs(temp)) ;
        ccPeakSim(trial,a) = temp(peaki) ;
        
        temp = xcov(excgShuff_Res(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a),inhgShuff_Res(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a),'coef') ;
        [m,peaki] = max(abs(temp)) ;
        ccPeakShuff(trial,a) = temp(peaki) ;
    end
end

meanccPeakSim = mean(ccPeakSim) ;
meanccPeakShuff = mean(ccPeakShuff) ;
    
for trial = 1:length(epochsSimRepeat) ; % for each trial
    for a=1:8 ; %for each bar
        SpikeRateSimRepeat(trial,a) = sum(spikeTrainSimRepeat(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a))/groupTime ; % hz
        SpikeRateSimRepeat_ON(trial,a) = sum(spikeTrainSimRepeat(trial,prePnts+groupPnts*a-groupPnts+1:prePnts+groupPnts*a-(groupPnts-OnPnts)))/OnTime ; % hz
        SpikeRateSimRepeat_OFF(trial,a) = sum(spikeTrainSimRepeat(trial,prePnts+groupPnts*a-groupPnts+OnPnts+1:prePnts+groupPnts*a))/(groupTime-OnTime) ; % hz
    end
end
 
% mean and std of firing rate for each bar angle
meanSpikeRateSim = mean(SpikeRateSim);
meanSpikeRateShuff = mean(SpikeRateShuff) ;
meanSpikeRateSimRepeat = mean(SpikeRateSimRepeat) ;

meanSpikeRateSim_norm = meanSpikeRateSim/max(mean(SpikeRateSim)) ;
meanSpikeRateShuff_norm = meanSpikeRateShuff/max(mean(SpikeRateShuff)) ;
meanSpikeRateSimRepeat_norm = meanSpikeRateSimRepeat/max(mean(SpikeRateSimRepeat)) ;

stdSpikeRateSim_norm = std(SpikeRateSim/max(mean(SpikeRateSim))) ;
stdSpikeRateShuff_norm = std(SpikeRateShuff/max(mean(SpikeRateShuff))) ;
stdSpikeRateSimRepeat_norm = std(SpikeRateSimRepeat/max(mean(SpikeRateSimRepeat))) ;

% vector strength and angle (% angle is not exactly right - see pair cads code)
BarAngles = [0:45:315] ;

for a=1:8 ; %for each bar
    
    vector_magX_Sim(a) = meanSpikeRateSim_norm (a)*cosd(BarAngles(a)) ; % magnitude of individual vector components
    vector_magY_Sim(a) = meanSpikeRateSim_norm (a)*sind(BarAngles(a)) ;
       
    vector_magX_Shuff(a) = meanSpikeRateShuff_norm(a)*cosd(BarAngles(a)) ;
    vector_magY_Shuff(a) = meanSpikeRateShuff_norm(a)*sind(BarAngles(a)) ;
           
end

vectorMagnitudeSim = sqrt(sum(vector_magX_Sim)^2+sum(vector_magY_Sim)^2) ; % magnitude of summed vector
vectorAngleSim = atand(sum(vector_magX_Sim)/sum(vector_magY_Sim)) ; % angle of summed vector
    
vectorMagnitudeShuff = sqrt(sum(vector_magX_Shuff)^2+sum(vector_magY_Shuff)^2) ;
vectorAngleShuff = atand(sum(vector_magX_Shuff)/sum(vector_magY_Shuff)) ;       

% calculate dsi (Talyor and Vaney 2002)
dsiSim = vectorMagnitudeSim/sum(meanSpikeRateSim_norm) ;
dsiShuff = vectorMagnitudeShuff/sum(meanSpikeRateShuff_norm) ;

% % figures
% figure
% %plot(BarAngles,meanSpikeRateSim)
% hold on
% errorbar(BarAngles,meanSpikeRateSim_norm,stdSpikeRateSim_norm,stdSpikeRateSim_norm,'k')
% 
% 
% figure
% set(gcf,'position',[101 25 1439 1064])
% subplot(3,2,1)
% errorbar(BarAngles,meanSpikeRateSim_norm,stdSpikeRateSim_norm,stdSpikeRateSim_norm,'k')
% hold on
% errorbar(BarAngles,meanSpikeRateShuff_norm,stdSpikeRateShuff_norm,stdSpikeRateShuff_norm,'g')
% errorbar(BarAngles,meanSpikeRateSimRepeat_norm,stdSpikeRateSimRepeat_norm,stdSpikeRateSimRepeat_norm,'y')
% % plot([0:45:315],SpikeRateSim,'k.')
% % plot([0:45:315],SpikeRateShuff,'g.')
% title('all spikes')
% xlabel('Bar angle')
% ylabel('spike rate norm (hz +/- std)')
% legend('sim','shuff','repeat')
% 
% subplot(3,2,2)
% plot(stdSpikeRateSim_norm,stdSpikeRateShuff_norm,'.')
% for a=1:8 ;
%     text(stdSpikeRateSim_norm(a),stdSpikeRateShuff_norm(a),num2str(BarAngles(a)))
% end
% line([0,max([stdSpikeRateSim_norm,stdSpikeRateShuff_norm])], [0,max([stdSpikeRateSim_norm,stdSpikeRateShuff_norm])])
% xlabel('std norm sim')
% ylabel('std norm shuff')
% 
% subplot(3,2,3)
% plot(time_cc,ccSim,'k')
% hold on
% plot(time_cc,ccShuff,'g')
% set(gca,'xlim',[-.2 .2])
% xlabel('lag time (sec)')
% ylabel('cc')
% 
% subplot(3,2,4)
% text([.1 ],[.9],['var exc sim = ',num2str(ResVarexcgSim_Mean)], 'Units','norm')
% text([.1 ],[.7],['var exc shuff = ',num2str(ResVarexcgShuff_Mean)], 'Units','norm')
% text([.1 ],[.5],['var inh sim = ',num2str(ResVarinhgSim_Mean)], 'Units','norm')
% text([.1 ],[.3],['var inh shuff = ',num2str(ResVarinhgShuff_Mean)], 'Units','norm')
% 
% subplot(3,2,5)
% errorbar([0:45:315],mean(SpikeRateSim_ON),-std(SpikeRateSim_ON),+std(SpikeRateSim_ON),'k')
% hold on
% errorbar([0:45:315],mean(SpikeRateShuff_ON),-std(SpikeRateShuff_ON),+std(SpikeRateShuff_ON),'g')
% errorbar([0:45:315],mean(SpikeRateSimRepeat_ON),-std(SpikeRateSimRepeat_ON),+std(SpikeRateSimRepeat_ON),'y')
% title('ON spikes')
% xlabel('Bar angle')
% ylabel('spike rate (hz +/- std)')
% 
% subplot(3,2,6)
% errorbar([0:45:315],mean(SpikeRateSim_OFF),-std(SpikeRateSim_OFF),+std(SpikeRateSim_OFF),'k')
% hold on
% errorbar([0:45:315],mean(SpikeRateShuff_OFF),-std(SpikeRateShuff_OFF),+std(SpikeRateShuff_OFF),'g')
% errorbar([0:45:315],mean(SpikeRateSimRepeat_OFF),-std(SpikeRateSimRepeat_OFF),+std(SpikeRateSimRepeat_OFF),'y')
% title('OFF spikes')
% xlabel('Bar angle')
% ylabel('spike rate (hz +/- std)')

% figure
% polar(0,max([vectorMagnitudeSim,vectorMagnitudeShuff,meanSpikeRateSim,meanSpikeRateShuff]))
% hold on
% 
% polar([vectorAngleSim,vectorAngleSim]*pi/180,[0,vectorMagnitudeSim],'k-')
% polar([vectorAngleShuff,vectorAngleShuff]*pi/180,[0,vectorMagnitudeShuff],'g-')
% 
% polar([0:45:315]*pi/180,meanSpikeRateSim,'k*')
% polar([0:45:315]*pi/180,meanSpikeRateShuff,'g*')

% figure(3)
% plot(meanccPeakSim,stdSpikeRateSim_norm/max(stdSpikeRateSim_norm),'k*')
% hold on
% plot(meanccPeakShuff,stdSpikeRateShuff_norm/max(stdSpikeRateShuff_norm),'g*')

% for Igor

% dsi for comparison
identifier = ['dsiSim',num2str(A)] ; 
ForIgor.(identifier) = dsiSim ;

identifier = ['dsiShuff',num2str(A)] ; 
ForIgor.(identifier) = dsiShuff ;


% std and cross corr coefs for comparsion plot
for a=1:8 ;
    identifier = ['stdFRSim',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = stdSpikeRateSim_norm(a) ; 

    identifier = ['stdFRShuff',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = stdSpikeRateShuff_norm(a) ;
    
%     identifier = ['ccPeakSim',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = meanccPeakSim(a) ; 
% 
%     identifier = ['ccPeakShuff',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = meanccPeakShuff(a) ;    
    
end
% 
% identifier = ['stdAll','cell',num2str(A)] ;
% ForIgor.(identifier) = [stdSpikeRateSim_norm,stdSpikeRateShuff_norm] ;
% 
% identifier = ['ccPeakAll','cell',num2str(A)] ;
% ForIgor.(identifier) = [meanccPeakSim,meanccPeakShuff] ;
% 
% 
% % vector magnitude and angle for polar plot
% identifier = ['polarvecMagSim',num2str(A)] ; 
% ForIgor.(identifier) = [0,vectorMagnitudeSim] ;
% 
% identifier = ['polarvecMagShuff',num2str(A)] ; 
% ForIgor.(identifier) = [0,vectorMagnitudeShuff] ;
% 
% identifier = ['polarvecAngleSim',num2str(A)] ; 
% ForIgor.(identifier) = [vectorAngleSim,vectorAngleSim] ;
% 
% identifier = ['polarvecAngleShuff',num2str(A)] ; 
% ForIgor.(identifier) = [vectorAngleShuff,vectorAngleShuff] ;
% 
% 
% 
% % mean spike rate and std for polar plots
% identifier = ['meanFRSim',num2str(A)] ;
% ForIgor.(identifier) = [meanSpikeRateSim_norm,meanSpikeRateSim_norm(1)] ;
% 
% identifier = ['meanFRShuff',num2str(A)] ;
% ForIgor.(identifier) = [meanSpikeRateShuff_norm,meanSpikeRateShuff_norm(1)] ;
% 
% identifier = ['polarBarAngles',num2str(A)] ;
% ForIgor.(identifier) = [BarAngles,BarAngles(1)] ;
% 
% for a=1:8 ;
%     identifier = ['polarStdSim',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = [meanSpikeRateSim_norm(a)-stdSpikeRateSim_norm(a),meanSpikeRateSim_norm(a)+stdSpikeRateSim_norm(a)] ;
%     
%     identifier = ['polarStdShuff',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = [meanSpikeRateShuff_norm(a)-stdSpikeRateShuff_norm(a),meanSpikeRateShuff_norm(a)+stdSpikeRateShuff_norm(a)] ;
%     
%     identifier = ['polarStdBarAngle',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = [BarAngles(a),BarAngles(a)] ;
% end

% mean spike rates raw and normalized linear
identifier = ['meanFRSimNormLin',num2str(A)] ;
ForIgor.(identifier) = meanSpikeRateSim_norm ;

% identifier = ['meanFRSimLin',num2str(A)] ;
% ForIgor.(identifier) = meanSpikeRateSim ;

identifier = ['stdFRSimNormLin',num2str(A)] ; % added 7/20/10
ForIgor.(identifier) = stdSpikeRateSim_norm ;

% mean spike rates raw and normalized linear % added 8/29/11
identifier = ['meanFRShuffNormLin',num2str(A)] ;
ForIgor.(identifier) = meanSpikeRateShuff_norm ;

identifier = ['stdFRShuffNormLin',num2str(A)] ; 
ForIgor.(identifier) = stdSpikeRateShuff_norm ;

% % individual voltage traces and spike rastors
stSim = nans(size(spikeTrainSim)) ;
stSim(spikeTrainSim==1)=1 ;

stShuff = nans(size(spikeTrainShuff)) ;
stShuff(spikeTrainShuff==1)=1 ;


for a=1:size(dataSim,1) ; % for each trial
    
%     identifier = ['dcvoltageSim',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = dataSim(a,:) ;
% 
%     identifier = ['dcvoltageShuff',num2str(a),'cell',num2str(A)] ;
%     ForIgor.(identifier) = dataShuff(a,:) ;
    
    identifier = ['SpikeTrainSim',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = stSim(a,:)*a ;

    identifier = ['SpikeTrainShuff',num2str(a),'cell',num2str(A)] ;
    ForIgor.(identifier) = stShuff(a,:)*a ;    

end
    
identifier = ['time',num2str(A)]
ForIgor.(identifier) = time ;

% % % current traces for reviewer
% rnd = 0 ;
% for a=1:2:size(I,1) ;
%     rnd = rnd+1 ;
%     identifier = ['DcDsISim',num2str(rnd),'cell',num2str(A)] ;
%     ForIgor.(identifier) = I(a,:) ;
%     
%     identifier = ['DcDsIShuff',num2str(rnd),'cell',num2str(A)] ;
%     ForIgor.(identifier) = I(a+1,:) ;
% end



