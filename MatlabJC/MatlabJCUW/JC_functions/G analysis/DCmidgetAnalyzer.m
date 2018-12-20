function ForIgor = DCmidgetAnalyzer(Input,Parameters,id,A) ;

% this will analyze dynamic clamp data from midget cells injected with a 
gaussianStd = 10.^[-3:.1:-1] ;

epochs = str2num(Input(A).(id)) ;

[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);

for a = 1:length(epochs) ;
    [Data(a,:), error] = ITCReadEpoch(epochs(a), 0, fp) ; 
    [excg(a,:), inhg(a,:), error(a,:)] = ITCReadEpochStmGClamp(epochs(a), 0, fp);
end

[SI, error] = ITCGetSamplingInterval(epochs(1), fp); % get sampling interval    
SI = SI * 1e-6; % Sampling interval in sec
if Input(A).ITC18flag == 1 ;
    SI = SI*1.25 ;
end
time = [SI:SI:SI*length(Data)] ;

SpikePnts = SpikeDetection_WC(Data,-30,10000) ;

spikeTrainSimRepeat = zeros(length(epochs)/4,length(Data)) ;
spikeTrainShuff = spikeTrainSimRepeat ;
spikeTrainSim = spikeTrainSimRepeat ;
spikeTrainSh = spikeTrainSimRepeat ;

round = 0 ;
for a = 1:4:length(epochs) ; % for each spike epoch
    round = round +1 ;

    dataSimRepeat(round,:) = Data(a,:) ;    % get data
    dataShuff(round,:) = Data(a+1,:) ;       
    dataSim(round,:) = Data(a+2,:) ;       
    dataSh(round,:) = Data(a+3,:) ; 
    
    spikeTrainSimRepeat(round,SpikePnts{a}) = 1 ; % make spike trains
    spikeTrainShuff(round,SpikePnts{a+1}) = 1 ;
    spikeTrainSim(round,SpikePnts{a+2}) = 1 ;
    spikeTrainSh(round,SpikePnts{a+3}) = 1 ;   
    
    spikeTimes_SimRepeat{round} = time(SpikePnts{a}) ; %spike times
    spikeTimes_Shuff{round} = time(SpikePnts{a+1}) ;
    spikeTimes_Sim{round} = time(SpikePnts{a+2}) ;
    spikeTimes_Sh{round} = time(SpikePnts{a+3}) ;
    
    excgSimRepeat(round,:) = excg(a,:) ; % conductances
    inhgSimRepeat(round,:) = inhg(a,:) ;
    
    excgShuff(round,:) = excg(a+1,:) ;
    inhgShuff(round,:) = inhg(a+1,:) ;
    
    excgSim(round,:) = excg(a+2,:) ;
    inhgSim(round,:) = inhg(a+2,:) ;
    
    excgSh(round,:) = excg(a+3,:) ;
    inhgSh(round,:) = inhg(a+3,:) ;
    
end

excgSimRepeat_Mean = mean(excgSimRepeat) ; %mean conductances
inhgSimRepeat_Mean = mean(inhgSimRepeat) ;

excgShuff_Mean = mean(excgShuff) ;
inhgShuff_Mean = mean(inhgShuff) ;

excgSim_Mean = mean(excgSim) ;
inhgSim_Mean = mean(inhgSim) ;

excgSh_Mean = mean(excgSh) ;
inhgSh_Mean = mean(inhgSh) ;


for a = 1:length(epochs)/4 ; % for each set of g
        
    excgShuff_Res(a,:) = excgShuff(a,:) - excgShuff_Mean ;
    inhgShuff_Res(a,:) = inhgShuff(a,:) - inhgShuff_Mean ;

    excgSim_Res(a,:) = excgSim(a,:) - excgSim_Mean ;
    inhgSim_Res(a,:) = inhgSim(a,:) - inhgSim_Mean ;
    
    excgSh_Res(a,:) = excgSh(a,:) - excgSh_Mean ;
    inhgSh_Res(a,:) = inhgSh(a,:) - inhgSh_Mean ;
    
    ccShuff(a,:) = xcov(excgShuff_Res(a,:),inhgShuff_Res(a,:),'coef') ;
    ccSim(a,:) = xcov(excgSim_Res(a,:),inhgSim_Res(a,:),'coef') ;
    ccSh(a,:) = xcov(excgSh_Res(a,:),inhgSh_Res(a,:),'coef') ;
end

ccShuff_Mean = mean(ccShuff) ;
ccSim_Mean = mean(ccSim) ;
ccSh_Mean = mean(ccSh) ;

ResVarexcgShuff_Mean = mean(var(excgShuff_Res,[],2)) ;
ResVarinhgShuff_Mean = mean(var(inhgShuff_Res,[],2)) ;

ResVarexcgSim_Mean = mean(var(excgSim_Res,[],2)) ;
ResVarinhgSim_Mean = mean(var(inhgSim_Res,[],2)) ;

ResVarexcgSh_Mean = mean(var(excgSh_Res,[],2)) ;
ResVarinhgSh_Mean = mean(var(inhgSh_Res,[],2)) ;

time_cc = [SI(1)*([1:length(ccSh_Mean)] - (length(ccSh_Mean)+1)/2)] ;

[sumPSTHvar_SimRepeat,psth_SimRepeat] = PsthVar(spikeTrainSimRepeat,gaussianStd ,1/SI) ; 
[sumPSTHvar_Shuff,psth_Shuff] = PsthVar(spikeTrainShuff,gaussianStd ,1/SI) ; 
[sumPSTHvar_Sim,psth_Sim] = PsthVar(spikeTrainSim,gaussianStd ,1/SI) ; 
[sumPSTHvar_Sh,psth_Sh] = PsthVar(spikeTrainSh,gaussianStd ,1/SI) ; 

psthi = find(gaussianStd<=.002,1,'last') ; % 2 ms

psth_SimRepeat_mean = mean(psth_SimRepeat{psthi}) ; % mean psth
psth_Shuff_mean = mean(psth_Shuff{psthi}) ; % 
psth_Sim_mean = mean(psth_Sim{psthi}) ; % 
psth_Sh_mean = mean(psth_Sh{psthi}) ; % 

psth_SimRepeat_std = std(psth_SimRepeat{psthi}) ; % psth var
psth_Shuff_std = std(psth_Shuff{psthi}) ; % 
psth_Sim_std = std(psth_Sim{psthi}) ; % 
psth_Sh_std = std(psth_Sh{psthi}) ; % 

psth_SimRepeat_snr = nans(1,length(psth_SimRepeat_std)) ;
psth_Shuff_snr = nans(1,length(psth_SimRepeat_std)) ;
psth_Sim_snr = nans(1,length(psth_SimRepeat_std)) ;
psth_Sh_snr = nans(1,length(psth_SimRepeat_std)) ;

psth_SimRepeat_snr(psth_SimRepeat_mean>0) = psth_SimRepeat_mean(psth_SimRepeat_mean>0)./psth_SimRepeat_std(psth_SimRepeat_mean>0)  ; % time dependent signal to noise ratio 
psth_Shuff_snr(psth_Shuff_mean>0) = psth_Shuff_mean(psth_Shuff_mean>0)./psth_Shuff_std(psth_Shuff_mean>0)  ;
psth_Sim_snr(psth_Sim_mean>0) = psth_Sim_mean(psth_Sim_mean>0)./psth_Sim_std(psth_Sim_mean>0)  ;
psth_Sh_snr(psth_Sh_mean>0) = psth_Sh_mean(psth_Sh_mean>0)./psth_Sh_std(psth_Sh_mean>0)  ;

psth_SimRepeat_snr_Mean = nanmean(psth_SimRepeat_snr) ;  % mean time dependant signal to noise ratio
psth_Shuff_snr_Mean = nanmean(psth_Shuff_snr) ; 
psth_Sim_snr_Mean = nanmean(psth_Sim_snr) ;
psth_Sh_snr_Mean = nanmean(psth_Sh_snr) ;

figure 
diffDshift = diff(circshift(Data,[0,1]),1,2) ;
diffD = diff(Data,1,2) ;
hist(Data(diffDshift>=0 & diffD<=0),[-80:40]); %finding local maxima
xlabel('voltage')
ylabel('number of observations')

figure
rasterPlot([spikeTimes_SimRepeat, spikeTimes_Shuff,spikeTimes_Sim,spikeTimes_Sh])
xlabel('time (sec)')
ylabel('trial (not in order)')

figure
subplot(3,1,1)
plot(time,psth_SimRepeat_mean,'k')
hold on
plot(time,psth_Shuff_mean,'r')
plot(time,psth_Sim_mean,'b')
plot(time,psth_Sh_mean,'c')
xlabel('time')
ylabel('firing rate (hz)')
legend('repeated','shuffled','simultaneous','single hold')

subplot(3,1,2)
plot(time,psth_SimRepeat_std,'k')
hold on
plot(time,psth_Shuff_std,'r')
plot(time,psth_Sim_std,'b')
plot(time,psth_Sh_std,'c')
xlabel('time')
ylabel('std of firing rate (hz)')

subplot(3,1,3)
plot(time,psth_SimRepeat_snr,'k')
hold on
plot(time,psth_Shuff_snr,'r')
plot(time,psth_Sim_snr,'b')
plot(time,psth_Sh_snr,'c')

plot(time,ones(1,length(time))*psth_SimRepeat_snr_Mean,'k')
hold on
plot(time,ones(1,length(time))*psth_Shuff_snr_Mean,'r')
plot(time,ones(1,length(time))*psth_Sim_snr_Mean,'b')
plot(time,ones(1,length(time))*psth_Sh_snr_Mean,'c')
xlabel('time')
ylabel('snr of firing rate (hz)')

figure
subplot(2,1,1)
plot(gaussianStd,sumPSTHvar_SimRepeat,'k')
hold on
plot(gaussianStd,sumPSTHvar_Shuff,'r')
plot(gaussianStd,sumPSTHvar_Sim,'b')
plot(gaussianStd,sumPSTHvar_Sh,'c')

subplot(2,1,2)
text(0,.9,['Shuff exc mean res var=',num2str(ResVarexcgShuff_Mean)],'units','normalized')
text(0,.8,['Sim exc mean res var=',num2str(ResVarexcgSim_Mean)],'units','normalized')
text(0,.7,['Sh exc mean res var=',num2str(ResVarexcgSh_Mean)],'units','normalized')

text(0,.5,['Shuff inh mean res var=',num2str(ResVarinhgShuff_Mean)],'units','normalized')
text(0,.4,['Sim inh mean res var=',num2str(ResVarinhgSim_Mean)],'units','normalized')
text(0,.3,['Sh inh mean res var=',num2str(ResVarinhgSh_Mean)],'units','normalized')

figure
subplot(2,1,1)
plot(time_cc,ccSh,'c')
hold on
plot(time_cc,ccSh_Mean,'c','linewidth',2)
xlabel('lag (sec)')
ylabel('correlation coef')

subplot(2,1,2)
plot(time_cc,ccSim,'b')
hold on
plot(time_cc,ccSim_Mean,'b','linewidth',2)
plot(time_cc,ccShuff,'r')
plot(time_cc,ccShuff_Mean,'r','linewidth',2)
xlabel('lag (sec)')
ylabel('correlation coef')

% for igor

identifier = ['sumPSTHvar_SimRepeat',id,num2str(A)] ;
ForIgor.(identifier) = sumPSTHvar_SimRepeat ; 

identifier = ['sumPSTHvar_Shuff',id,num2str(A)] ;
ForIgor.(identifier) = sumPSTHvar_Shuff ; 

identifier = ['sumPSTHvar_Sim',id,num2str(A)] ;
ForIgor.(identifier) = sumPSTHvar_Sim ;

identifier = ['sumPSTHvar_Sh',id,num2str(A)] ;
ForIgor.(identifier) = sumPSTHvar_Sh ;

identifier = ['snr_SimRepeat',id,num2str(A)] ;
ForIgor.(identifier) = psth_SimRepeat_snr_Mean ; 

identifier = ['snr_Shuff',id,num2str(A)] ;
ForIgor.(identifier) = psth_Shuff_snr_Mean ; 

identifier = ['snr_Sim',id,num2str(A)] ;
ForIgor.(identifier) = psth_Sim_snr_Mean ;

identifier = ['snr_Sh',id,num2str(A)] ;
ForIgor.(identifier) = psth_Sh_snr_Mean ;
