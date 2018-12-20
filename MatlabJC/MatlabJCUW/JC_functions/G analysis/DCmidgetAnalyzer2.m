function ForIgor = DCmidgetAnalyzer2(Input,Parameters,id,id2,A) ;

% this will analyze dynamic clamp data from midget cells injected with sim,
% shuffled, or sim repeated g .  This differs from DCmidgetAnalyzer because
% it does not include single hold epochs and sim repeated epochs are
% seperated from sim vs shuffled

% JC 9/23/09

% gaussianStd = 10.^[-3:.1:-1] ;

epochs = str2num(Input(A).(id)) ; %sim and shuffled g
epochsSimRepeat = str2num(Input(A).(id2)) ; % sim repeated

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);

for a = 1:length(epochs) ;
    [Data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ; 
    [excg(a,:), inhg(a,:), error(a,:)] = ITCReadEpochStmGClamp(epochs(a), 0, fp);
    %I(a,:) = (Data(a,:)+80).*(inhg(a,:)) + Data(a,:).*excg(a,:) ; % added this line 8/4/10 to calculate current
end

for a = 1:length(epochsSimRepeat) ;
    [dataSimRepeat(a,:), error] = ITCReadEpoch(epochsSimRepeat(a), 0, fp) ; 
    [excgSimRepeat(a,:), inhgSimRepeat(a,:), error(a,:)] = ITCReadEpochStmGClamp(epochsSimRepeat(a), 0, fp);
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
spikeTrainShuff = zeros(length(SpikePnts)/2,length(Data)) ;
spikeTrainSim = spikeTrainShuff ;


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

    ccSim(a,:) = xcorr(excgSim_Res(a,:),inhgSim_Res(a,:),'coef') ;
    ccShuff(a,:) = xcorr(excgShuff_Res(a,:),inhgShuff_Res(a,:),'coef') ;
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

% [meanPSTHsnr_SimRepeat,sumPSTHvar_SimRepeat,psth_SimRepeat] = PsthVar(spikeTrainSimRepeat,gaussianStd ,1/SI) ; 
% [meanPSTHsnr_Shuff,sumPSTHvar_Shuff,psth_Shuff] = PsthVar(spikeTrainShuff,gaussianStd ,1/SI) ; 
% [meanPSTHsnr_Sim,sumPSTHvar_Sim,psth_Sim] = PsthVar(spikeTrainSim,gaussianStd ,1/SI) ; 

% minIntTime_Sim = interp1(meanPSTHsnr_Sim,gaussianStd,1) ;
% minIntTime_Shuff = interp1(meanPSTHsnr_Shuff,gaussianStd,1) ;

[powerX_Sim,meanSpikeSpectrum_Sim,resSpikeSpectrum_Sim,meanSpikeSpectrum_smth_Sim,resSpikeSpectrum_smth_Sim,snrSpikeSpectrum_smth_Sim,VarSumResSpikeSpectrum_Sim] = snrSpikeSpectrum(spikeTrainSim,1/SI,.001) ;
[powerX_Shuff,meanSpikeSpectrum_Shuff,resSpikeSpectrum_Shuff,meanSpikeSpectrum_smth_Shuff,resSpikeSpectrum_smth_Shuff,snrSpikeSpectrum_smth_Shuff,VarSumResSpikeSpectrum_Shuff] = snrSpikeSpectrum(spikeTrainShuff,1/SI,.001) ;
[powerX_SimRepeat,meanSpikeSpectrum_SimRepeat,resSpikeSpectrum_SimRepeat,meanSpikeSpectrum_smth_SimRepeat,resSpikeSpectrum_smth_SimRepeat,snrSpikeSpectrum_smth_SimRepeat,VarSumResSpikeSpectrum_SimRepeat] = snrSpikeSpectrum(spikeTrainSimRepeat,1/SI,.001) ;

powerMin = find(powerX_Sim>=1,1,'first') ;
powerMax = find(powerX_Sim<=20,1,'last') ;

MeanSnrSim = mean(snrSpikeSpectrum_smth_Sim) ;
MeanSnrShuff = mean(snrSpikeSpectrum_smth_Shuff) ;
MeanSnrSimRepeat = mean(snrSpikeSpectrum_smth_SimRepeat) ;

PredSnrFactor = 1/sqrt((1-corrCoefShuff)/(1-corrCoefSim)) ;
SnrFactor = MeanSnrShuff/MeanSnrSim ;

ResVarSim = sum(resSpikeSpectrum_smth_Sim.PowerSpec(2:end)) ;
ResVarShuff = sum(resSpikeSpectrum_smth_Shuff.PowerSpec(2:end)) ;
ResVarSimRepeat = sum(resSpikeSpectrum_smth_SimRepeat.PowerSpec(2:end)) ;

PredResVarFactor = (1-corrCoefSim)/(1-corrCoefShuff) ;
ResVarFactor = ResVarSim/ResVarShuff ;

SumSnrSim = sum(meanSpikeSpectrum_Sim(powerMin:powerMax))/sum(resSpikeSpectrum_Sim(powerMin:powerMax)) ; 
SumSnrShuff = sum(meanSpikeSpectrum_Shuff(powerMin:powerMax))/sum(resSpikeSpectrum_Shuff(powerMin:powerMax)) ; 
SumSnrSimRepeat = sum(meanSpikeSpectrum_SimRepeat(powerMin:powerMax))/sum(resSpikeSpectrum_SimRepeat(powerMin:powerMax)) ; 


% % check spike precision with spike distance metric
% for a=1:length(spikeTimes_Sim) ;
%     stms_Sim{a} = spikeTimes_Sim{a}*1000 ;
%     stms_Shuff{a} = spikeTimes_Shuff{a}*1000 ;
% end
% 
% [DeltaT, tli_spike, tlj_spike, Percent_Pairs_Quantified, DeltaT_Histogram,XSim,DeltaT_CumProb_Sim] = DeltaT_Distribution_From_SCR(stms_Sim,.02) ;
% [DeltaT, tli_spike, tlj_spike, Percent_Pairs_Quantified, DeltaT_Histogram,XShuff,DeltaT_CumProb_Shuff] = DeltaT_Distribution_From_SCR(stms_Shuff,.02) ;


% figure 
% diffDshift = diff(circshift(Data,[0,1]),1,2) ;
% diffD = diff(Data,1,2) ;
% hist(Data(diffDshift>=0 & diffD<=0),[-80:40]); %finding local maxima
% xlabel('voltage')
% ylabel('number of observations')

% figure
% rasterPlot([spikeTimes_SimRepeat, spikeTimes_Shuff,spikeTimes_Sim])
% xlabel('time (sec)')
% ylabel('trial (not in order)')
% 
% figure
% plot(time_cc,ccSim_Mean,'k')
% hold on
% plot(time_cc,ccShuff_Mean,'g')
% plot(0,corrCoefSim,'ko')
% plot(0,corrCoefShuff,'go')
% axis([-.05 .05 -.2 .7])
% xlabel('lag (sec)')
% ylabel('cc')
% legend('sim','shuff')
% 
% figure
% set(gcf,'position',[358,24,1472,1060]) ;
% subplot(5,1,1)
% plot(meanSpikeSpectrum_smth_Sim.Freq,snrSpikeSpectrum_smth_Sim,'k')
% hold on
% plot(meanSpikeSpectrum_smth_Sim.Freq,snrSpikeSpectrum_smth_Shuff,'g')
% plot(meanSpikeSpectrum_smth_Sim.Freq,snrSpikeSpectrum_smth_Sim*PredSnrFactor,'c')
% %plot(meanSpikeSpectrum_smth_SimRepeat.Freq,snrSpikeSpectrum_smth_SimRepeat,'y')
% set(gca,'xscale','log')
% xlabel('frequency')
% ylabel('snr')
% legend('sim','shuff','prediction shuff')
% 
% subplot(5,1,2)
% plot(meanSpikeSpectrum_smth_Sim.Freq,snrSpikeSpectrum_smth_Sim./snrSpikeSpectrum_smth_Shuff,'k')
% set(gca,'xscale','log')
% set(gca,'ylim',[0 2])
% xlabel('frequency')
% ylabel('sim snr / shuff snr')
% 
% subplot(5,1,3)
% plot(meanSpikeSpectrum_smth_Sim.Freq,meanSpikeSpectrum_smth_Sim.PowerSpec,'k')
% hold on
% plot(meanSpikeSpectrum_smth_Sim.Freq,resSpikeSpectrum_smth_Sim.PowerSpec,'k--')
% plot(meanSpikeSpectrum_smth_Sim.Freq,meanSpikeSpectrum_smth_Shuff.PowerSpec,'g')
% plot(meanSpikeSpectrum_smth_Sim.Freq,resSpikeSpectrum_smth_Shuff.PowerSpec,'g--')
% plot(meanSpikeSpectrum_smth_Sim.Freq,resSpikeSpectrum_smth_Shuff.PowerSpec*PredResVarFactor,'c--')
% set(gca,'xscale','log','yscale','log')
% xlabel('frequency')
% ylabel('power')
% legend('mean sim','res sim','mean shuff', 'res shuff','pred res sim')
% 
% subplot(5,1,4)
% plot(powerX_Sim,meanSpikeSpectrum_Sim,'k')
% hold on
% plot(meanSpikeSpectrum_smth_Sim.Freq,meanSpikeSpectrum_smth_Sim.PowerSpec,'y')
% plot(powerX_Sim,resSpikeSpectrum_Sim,'k--')
% plot(resSpikeSpectrum_smth_Sim.Freq,resSpikeSpectrum_smth_Sim.PowerSpec,'y--')
% set(gca,'xscale','log','yscale','log')
% xlabel('frequency')
% ylabel('power')
% legend('full ps','smoothed ps')
% 
% subplot(5,1,5)
% plot(powerX_Shuff,meanSpikeSpectrum_Shuff,'g')
% hold on
% plot(meanSpikeSpectrum_smth_Shuff.Freq,meanSpikeSpectrum_smth_Shuff.PowerSpec,'y')
% plot(powerX_Shuff,resSpikeSpectrum_Shuff,'g--')
% plot(resSpikeSpectrum_smth_Shuff.Freq,resSpikeSpectrum_smth_Shuff.PowerSpec,'y--')
% set(gca,'xscale','log','yscale','log')
% xlabel('frequency')
% ylabel('power')
% legend('full ps','smoothed ps')

% figure
% plot(XSim,DeltaT_CumProb_Sim,'k')
% hold on
% plot(XShuff,DeltaT_CumProb_Shuff,'g')

% for igor

% identifier = ['PredSnrFactor',num2str(A)] ;
% ForIgor.(identifier) = PredSnrFactor ;
% 
% identifier = ['SnrFactor',num2str(A)] ;
% ForIgor.(identifier) = SnrFactor ;
% 
identifier = ['MeansnrSim',num2str(A)] ;
ForIgor.(identifier) = MeanSnrSim;

identifier = ['MeansnrShuff',num2str(A)] ;
ForIgor.(identifier) = MeanSnrShuff;
%
% identifier = ['MeansnrSimRepeat',num2str(A)] ;
% ForIgor.(identifier) = MeanSnrSimRepeat;

identifier = ['SumsnrSim',num2str(A)] ;
ForIgor.(identifier) = SumSnrSim;

identifier = ['SumsnrShuff',num2str(A)] ;
ForIgor.(identifier) = SumSnrShuff;

% identifier = ['SumsnrSimRepeat',num2str(A)] ;
% ForIgor.(identifier) = SumSnrSimRepeat;


% % whole mean spike spectrums, residual spike spectrums and snr spectrums
% identifier = ['powerSpecX',num2str(A)] ;
% ForIgor.(identifier) = meanSpikeSpectrum_smth_Sim.Freq ;
% 
% 
% identifier = ['meanSpikePowerSim',num2str(A)] ;
% ForIgor.(identifier) = meanSpikeSpectrum_smth_Sim.PowerSpec ;
% 
% identifier = ['meanSpikePowerShuff',num2str(A)] ;
% ForIgor.(identifier) = meanSpikeSpectrum_smth_Shuff.PowerSpec ;
% 
% 
% identifier = ['resSpikePowerSim',num2str(A)] ;
% ForIgor.(identifier) = resSpikeSpectrum_smth_Sim.PowerSpec ;
% 
% identifier = ['resSpikePowerShuff',num2str(A)] ;
% ForIgor.(identifier) = resSpikeSpectrum_smth_Shuff.PowerSpec ;
% 
% 
% identifier = ['snrSpikePowerSim',num2str(A)] ;
% ForIgor.(identifier) = snrSpikeSpectrum_smth_Sim ;
% 
% identifier = ['snrSpikePowerShuff',num2str(A)] ;
% ForIgor.(identifier) = snrSpikeSpectrum_smth_Shuff ;
% 
% % individual voltage traces and spike rastors
% stSim = nans(size(spikeTrainSim)) ;
% stSim(spikeTrainSim==1)=1 ;
% 
% stShuff = nans(size(spikeTrainShuff)) ;
% stShuff(spikeTrainShuff==1)=1 ;
% 
% stSimRepeat = nans(size(spikeTrainSimRepeat)) ;
% stSimRepeat(spikeTrainSimRepeat==1)=1 ;
% 
% r = 0 ;
% for a=1:2:size(Data,1) ; % for each trial
%     r=r+1 ;
%     
%     identifier = ['dcvoltageSim',num2str(r),'id',num2str(A)] ;
%     ForIgor.(identifier) = Data(a,:) ;
% 
%     identifier = ['dcvoltageShuff',num2str(r),'id',num2str(A)] ;
%     ForIgor.(identifier) = Data(a+1,:) ;
%        
%     identifier = ['spikeTrainSim',num2str(r),'id',num2str(A)] ;
%     ForIgor.(identifier) = stSim(r,:)*r ;
%     
%     identifier = ['spikeTrainShuff',num2str(r),'id',num2str(A)] ;
%     ForIgor.(identifier) = stShuff(r,:)*r ;
% end
% 
% for a=1:size(spikeTrainSimRepeat,1) ;
%     identifier = ['spikeTrainSimRepeat',num2str(a),'id',num2str(A)] ;
%     ForIgor.(identifier) = stSimRepeat(a,:)*a ;
% end
% 
% ForIgor.time = time ;
    


