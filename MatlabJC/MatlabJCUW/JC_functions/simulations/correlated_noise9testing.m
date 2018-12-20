% this simulation will assess effect of correlated noise of varaibility of
% the spike output.


%parameters
NoiseCom_Frac = [0:1] ; % range of common noise fraction of total
NoiseTot_Var = .5 ; % var of total noise in each conductance
NumTrials = 2 ; % number of repeats 

params.Einh = -80 ;
params.Eexc = 0 ;
params.Eleak = -60 ;
params.cap = .06 ;%nF
params.Vthresh = 10000 ;
params.AbsRef  = .002 ; %sec
params.RelRefTau = .008 ;%sec
params.RelRefAmp = 4 ; %mV
params.Gleak = 1 ;
params.SampleRate = 10000 ;

MaxDeltaTvector = 10.^[-3:.1:-1] ; % sec vector of times shifts for max delta t

si = .0001 ; % sample interval (sec)
time = [0:si:10] ; % sec

lt = length(time) ;

excLF = simFilter([0:si:.3],.050,.008,.12,.10,.05,.02) ; % linear filters(time,tpeak,peakRise,peakAmp,ttrough,troughDecay,troughAmp)
inhLF = simFilter([0:si:.3],.055,.009,.13,.11,.06,.02) ;

NFgainExc = 1 ;
NFgainInh = 1 ;

NFthreshExc = 0 ;
NFthreshInh = 0 ;

NFmaxExc = 300 ;
NFmaxInh = 300 ;


% set signal
signal = normrnd(2,2,1,lt) ;     % signal
signal = lowPassFilter(signal,1/si,60) ;    % low pass at 60hz

Trial = 0 ;
spikes = zeros(NumTrials*length(NoiseCom_Frac),lt) ;

for a = 1:length(NoiseCom_Frac) ; % for each new change in common noise amplitude
    
    NoiseCom_Std = sqrt(NoiseTot_Var * NoiseCom_Frac(a)) ; % standard deviation of common noise
    NoiseInd_Std = sqrt(NoiseTot_Var - NoiseCom_Std^2) ; % standard deviation of independant noise
    
    for b=1:NumTrials ; % for each trial
        Trial = Trial + 1 ;
        
        NoiseCommon = normrnd(0,NoiseCom_Std,1,lt) ; % create noise for this trial
        NoiseExc = normrnd(0,NoiseInd_Std,1,lt) ;
        NoiseInh = normrnd(0,NoiseInd_Std,1,lt) ;
        
        NoiseCommon = lowPassFilter(NoiseCommon,1/si,1000) ; % low pass at 1000hz
        NoiseExc = lowPassFilter(NoiseExc,1/si,1000) ; % low pass at 1000hz
        NoiseInh = lowPassFilter(NoiseInh,1/si,1000) ; % low pass at 1000hz
        
        GexcPreLN = signal + NoiseCommon + NoiseExc ; % add signal and noise
        GinhPreLN = signal + NoiseCommon + NoiseInh ;

        % pass signal and noise through linear and nonlinear filters

        GexcPreN = conv(GexcPreLN,excLF) ; % convolve linear filter
        GinhPreN = conv(GinhPreLN,inhLF) ; % convolve linear filter
        
        % pass through nonlinear filter(threshold, gain, saturate)
        GexcPre = GexcPreN*NFgainExc - NFgainExc*NFthreshExc ; % gain
        GinhPre = GinhPreN*NFgainInh - NFgainExc*NFthreshInh ;
        
        GexcPre(GexcPreN<NFthreshExc) = 0 ; %threshold
        GinhPre(GinhPreN<NFthreshInh) = 0 ; %threshold
        
        GexcPre(GexcPre>NFmaxExc) = NFmaxExc ; % saturate
        GinhPre(GinhPre>NFmaxInh) = NFmaxInh ; 
        
        Gexc(Trial,:) = GexcPre(1,(length(GexcPre)-lt)/2+1:((length(GexcPre)-lt)/2)+lt) ; % cut off convolution extra
        Ginh(Trial,:) = GinhPre(1,(length(GexcPre)-lt)/2+1:((length(GexcPre)-lt)/2)+lt) ;
        Gleak = 0 ; % leak conductance
                
        V_trace(Trial,:) = LIFmodelG_c(Gexc(Trial,:),Ginh(Trial,:),params) ; % run LIF model
        spikes(Trial,V_trace(Trial,:)==50) = 1 ; % spike train
        spiketimes{Trial} = time(V_trace(Trial,:)==50) ; % spike times
        
        V_filtered = V_trace(Trial,:) ; % filter out spikes to calculate synaptic current
        V_filtered(spikes(Trial,:) == 1) = params.Vthresh ;
        I(Trial,:) = Gexc(Trial,:).*(V_filtered-params.Eexc) + Ginh(Trial,:).*(V_filtered-params.Einh) + Gleak*(V_filtered-params.Eleak) ; % calculate synaptic current
    end
    GexcMean(a,:) = mean(Gexc(Trial-NumTrials+1:Trial,:)) ; % mean Gexc for this particular set of trials
    GinhMean(a,:) = mean(Ginh(Trial-NumTrials+1:Trial,:)) ; % mean Ginh for this particular set of trials
    V_traceMean(a,:) = mean(V_trace(Trial-NumTrials+1:Trial,:)) ; % mean voltage for this particular set of trials
    IMean(a,:) = mean(I(Trial-NumTrials+1:Trial,:)) ; %
end

% get cross correlations of estimated noise (residuals), estimate and predict total current noise, measure voltage noise (only applicable for nonspiking) 
lcc = length(GexcMean)*2-1 ; % length of cc
ccX = ([1:lcc] - (lcc+1)/2)*si ; % x values of cc
ccZeroi = (lcc+1)/2 ; % index of zero time lag
for a = 1:length(NoiseCom_Frac) ; % for each new change in common noise amplitude
    for b = a*NumTrials-NumTrials+1:a*NumTrials ; % for each trial         
     
        ccG(b,:) = xcov(GexcRes(b,:),GinhRes(b,:)) ; % cross corr of residuals 
        
        InoiseVar(b) = var(IRes(b,:)) ; % variance of current noise
        InoiseVarEst(b) = var(GexcRes(b,:)*(-45 - params.Eexc)+GinhRes(b,:)*(-45 - params.Einh)) ; % estimate variance of current noise
        InoiseVarPred(b) = mean(var(GexcRes(b,:)*(-45 - params.Eexc),[],2)) + mean(var(GinhRes(b,:)*(-45 - params.Einh),[],2)) + 2*ccIRes(b,ccZeroi)/lt ; % theoretical prediction based on cross correlation
        VnoiseVar(b) = var(VRes(b,:)) ; % variance of voltage
    end
    ccGResMean(a,:) = mean(ccGRes(a*NumTrials-NumTrials+1:a*NumTrials,:)) ; % cross correlation mean
    ccIResMean(a,:) = mean(ccIRes(a*NumTrials-NumTrials+1:a*NumTrials,:)) ; % cross correlation mean
    
    GexcVarMean(a) = mean(var(GexcRes(a*NumTrials-NumTrials+1:a*NumTrials,:),[],2)) ; % mean variance of gexc noise
    GinhVarMean(a) = mean(var(GinhRes(a*NumTrials-NumTrials+1:a*NumTrials,:),[],2)) ; % mean variance of ginh noise
    
    InoiseVarMean(a) = mean(InoiseVar(a*NumTrials-NumTrials+1:a*NumTrials)) ; % mean of current noise variance
    InoiseVarEstMean(a) = mean(InoiseVarEst(a*NumTrials-NumTrials+1:a*NumTrials)) ; % mean estimate of current noise variance
    
    InoiseVarPredMean(a) = mean(InoiseVarPred(a*NumTrials-NumTrials+1:a*NumTrials)) ; % mean estimate of current noise variance
    VnoiseVarMean(a) = mean(VnoiseVar(a*NumTrials-NumTrials+1:a*NumTrials)) ; % mean variance of voltage
end

ccGResCoefPeak = max(ccGResCoef,[],2) ; % peak of cc coef
ccGResCoefTrough = min(ccGResCoef,[],2) ; % trough of cc coef

% % check spike precision with spike match metric
% for a = 1:length(NoiseCom_Frac) ; % for each new change in common noise amplitude
%     [matrix,meanvector] = SpikeMatcher(spiketimes(a*NumTrials-NumTrials+1:a*NumTrials), MaxDeltaTvector) ;
%     SpkMIndPairs{a} = matrix ;
%     SpkMMean(a,:) = meanvector ; 
% end

% check spike precision with psth var metric
for a = 1:length(NoiseCom_Frac) ;
    [meanPSTHsnr(a,:),sumPSTHvar(a,:),psth{a}] = PsthVar(spikes(a*NumTrials-NumTrials+1:a*NumTrials,:),MaxDeltaTvector,1/si);
    %sumPSTHvar(a,:) = temp ;
end

% check signal to noise ratio of psth
psthi = find(MaxDeltaTvector<=.002,1,'last') ; % 2 ms
for a = 1:length(NoiseCom_Frac) ;
    mean_psth(a,:) = mean(psth{a}{psthi}) ;
    std_psth(a,:) = std(psth{a}{psthi}) ;
    snr(a,:)=nans(1,length(mean_psth(a,:))) ;
    snr(a,mean_psth(a,:)>0)=mean_psth(a,mean_psth(a,:)>0)./std_psth(a,mean_psth(a,:)>0) ;
    snrMean(a) = nanmean(snr(a,:)) ; 
end

% calculate infomation per spike and per second with psth (equation taken from Sincich et al 2009)
for a = 1:length(NoiseCom_Frac) ;
    InformationSpike(a) = 1/max(time)*sum(si*(mean_psth(a,mean_psth(a,:)>0)./mean(mean_psth(a,:))).*log2((mean_psth(a,mean_psth(a,:)>0)/mean(mean_psth(a,:))))) ;
    InformationSecond(a) = InformationSpike(a)*mean(mean_psth(a,:)) ; % psth is in hz thus it seemed reasonable to multiple by mean psth (spikes/sec) to get info/sec
end


% % check spike precision with spike distance metric 
% for a = 1:length(NoiseCom_Frac) ; % for each new change in common noise amplitude
%     for b = a*NumTrials-NumTrials+1:a*NumTrials ; % for each trial
%         for c = a*NumTrials-NumTrials+1:a*NumTrials ; % for each possible pair
%             sd{a}(b,c)=spkd_c(spiketimes{b},spiketimes{c},length(spiketimes{b}), length(spiketimes{c}), .02); % runs spike distance function at a given cost
%         end
%     end
%     sd_Mean(a) = mean(sd{a}(sd{a}~=0)) ; % get the mean spike distance for a particular amplitude of common noise
% end

% check spike precision with spike spectrum measure
for a = 1:length(NoiseCom_Frac) ;
    [powerX,meanSpikeSpectrum,resSpikeSpectrum,meanSpikeSpectrum_smth,resSpikeSpectrum_smth,snrSpikeSpectrum_smth,VarSumResSpikeSpectrum] = snrSpikeSpectrum(spikes(a*NumTrials-NumTrials+1:a*NumTrials,:),1/si,.001) ; 

    snrX(a,:) = meanSpikeSpectrum_smth.Freq ;
    snrSpectrum(a,:) = snrSpikeSpectrum_smth ;
    meanSnr(a) = mean(snrSpectrum(a,:)) ; %#ok<*AGROW>
end


% fit a line to Inoise as a function of cc(0) 
bITheory = GexcVarMean(1)*(-45 - params.Eexc)^2 + GinhVarMean(1)*(-45 - params.Einh)^2 ; % the y intecept (when no correlation) variances sum

figure
rasterPlot(spiketimes)
xlabel('time (sec)')
ylabel('trial')

color = ['b','g','r','k','c','y','b','g','r','k','c','y'] ;

figure
for a = 1:length(NoiseCom_Frac) ; % for each new change in common noise amplitude 
    subplot(2,1,1)
    plot(mean(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,ccZeroi)),InoiseVarMean(a),[color(a),'*'])
    hold on
    plot(mean(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,ccZeroi)),InoiseVarEstMean(a),[color(a),'+'])
    plot([0,1],[InoiseVarMean(1),0],'k--')
    plot([0,1],[InoiseVarEstMean(1),0],'b--')
    legend('from filtered I','from estimated I')
    xlabel('mean cc')
    ylabel('I variance')
    
    subplot(2,1,2)
    plot(mean(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,ccZeroi)),meanSnr(a),[color(a),'*'])
    hold on

    PredictedSnr(a) = meanSnr(1)*sqrt((1-mean(ccGResCoef(1*NumTrials-NumTrials+1:1*NumTrials,ccZeroi)))/(1-mean(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,ccZeroi)))) ;
    plot(mean(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,ccZeroi)),PredictedSnr(a),[color(a),'o'])
    
    legend('actual snr','predicted snr')
    xlabel('mean cc')
    ylabel('meanSNR')
end

figure
for a = 1:length(NoiseCom_Frac) ; % for each new change in common noise amplitude   
    subplot(4,2,1:2)
    plot(NoiseCom_Frac(a),ccGResCoefPeak(a*NumTrials-NumTrials+1:a*NumTrials),[color(a),'*'])
    hold on
    plot(NoiseCom_Frac(a),ccGResCoefTrough(a*NumTrials-NumTrials+1:a*NumTrials),[color(a),'o'])
    xlabel('fraction common noise')
    ylabel('peak or trough of cc coef')
    legend('peak','Trough')
        
    subplot(4,2,3:4)
    plot(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,ccZeroi),InoiseVar(a*NumTrials-NumTrials+1:a*NumTrials),[color(a),'*'])
    hold on
    plot(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,ccZeroi),InoiseVarPred(a*NumTrials-NumTrials+1:a*NumTrials),[color(a),'o'])
    xlabel('cc coef (0)')
    ylabel('I noise variance estimate')
    legend('actual estimate', 'predicted')

    if params.Vthresh>50 ;
        subplot(4,2,5:6)
        plot(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,ccZeroi),VnoiseVar(a*NumTrials-NumTrials+1:a*NumTrials),[color(a),'*'])
        hold on
        xlabel('cc coef (0)')
        ylabel('V noise variance')
    else
        subplot(4,2,5)
        plot(MaxDeltaTvector,sumPSTHvar(a,:),[color(a),'-'])
        hold on
        plot(MaxDeltaTvector(8),sumPSTHvar(a,8),[color(a),'o'])
        xlabel('gaussian var (sec)')
        ylabel('sum time dependent PSTH var')
        h=gca ;
        set(h,'xscale','log')
        
        subplot(4,2,6)
        plot(InoiseVarPredMean(a),sumPSTHvar(a,8),[color(a),'o'])
        hold on
%         plot(InoiseVarPredMean(a),sumPSTHvar(a,11),[color(a),'*'])
%         plot(InoiseVarPredMean(a),sumPSTHvar(a,21),[color(a),'o'])
        xlabel('Inoise var pred')
        ylabel('sum time dependent PSTH var')
    end
        
    subplot(4,1,4)
    plot(ccX,ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,:),[color(a)])
    hold on
    xlabel('lag (sec)')
    ylabel('cross corr coef')
end

figure
for a = 1:length(NoiseCom_Frac) ; % for each new change in common noise amplitude  
    subplot(5,2,1:2)
    plot(time,mean_psth(a,:),[color(a)])
    hold on
    plot(time,ones(1,length(time))*mean(mean_psth(a,:)),[color(a)])
    xlabel('time (sec)')
    ylabel('firing rate (Hz)')
    
    subplot(5,2,3:4)
    plot(time,std_psth(a,:),[color(a)])
    hold on
    plot(time,ones(1,length(time))*mean(std_psth(a,:)),[color(a)])
    xlabel('time (sec)')
    ylabel('std of firing rate (Hz)')
    
    subplot(5,2,5:6)
    plot(time,snr(a,:),[color(a)])
    hold on
    plot(time,ones(1,length(time))*mean(snr(a,:)),[color(a)])
    xlabel('time (sec)')
    ylabel('signal to noise ratio')
    
    subplot(5,2,7:8)
    plot(NoiseCom_Frac(a),snrMean(a),[color(a),'*'])
    hold on
    xlabel('fraction common noise')
    ylabel('mean snr')
    
    subplot(5,2,9)
    plot(NoiseCom_Frac(a),InformationSpike(a),[color(a),'*'])
    hold on
    xlabel('fraction common noise')
    ylabel('Information/Spike')
    
    subplot(5,2,10)
    plot(NoiseCom_Frac(a),InformationSecond(a),[color(a),'*'])
    hold on
    xlabel('fraction common noise')
    ylabel('Information/Second')
    
end

figure
for a = 1:length(NoiseCom_Frac) ; % for each new change in common noise amplitude   
    subplot(4,1,1)
    plot(NoiseCom_Frac(a),ccGResCoefPeak(a*NumTrials-NumTrials+1:a*NumTrials),[color(a),'*'])
    hold on
    plot(NoiseCom_Frac(a),ccGResCoefTrough(a*NumTrials-NumTrials+1:a*NumTrials),[color(a),'o'])
    xlabel('fraction common noise')
    ylabel('peak or trough of cc coef')
    legend('peak','Trough')
        
    subplot(4,1,2)
    plot(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,ccZeroi),InoiseVar(a*NumTrials-NumTrials+1:a*NumTrials),[color(a),'*'])
    hold on
    plot(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,ccZeroi),InoiseVarPred(a*NumTrials-NumTrials+1:a*NumTrials),[color(a),'o'])
    xlabel('cc coef (0)')
    ylabel('I noise variance estimate')
    legend('actual estimate', 'predicted')

    if params.Vthresh>50 ;
        subplot(4,1,3)
        plot(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,ccZeroi),VnoiseVar(a*NumTrials-NumTrials+1:a*NumTrials),[color(a),'*'])
        hold on
        xlabel('cc coef (0)')
        ylabel('V noise variance')
    else
        subplot(4,1,3)
        plot(MaxDeltaTvector,SpkMMean(a,:),[color(a),'-'])
        hold on
        xlabel('time shift (sec)')
        ylabel('percent of spikes that can be matched')
        h=gca ;
        set(h,'xscale','log')
    end
        
    subplot(4,1,4)
    plot(ccX,ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,:),[color(a)])
    hold on
    xlabel('lag (sec)')
    ylabel('cross corr coef')
end

subplot(4,1,2)
plot(ccGResCoef(:,ccZeroi),-bITheory*ccGResCoef(:,ccZeroi)+bITheory,'k--') 


figure
for a = 1:length(NoiseCom_Frac) ; % for each new change in common noise amplitude
    subplot(1,3,1)
    plot(NoiseCom_Frac(a),InoiseVarMean(a),[color(a),'*'])
    hold on
    xlabel('fraction of noise which is common')
    ylabel('mean variance I estimate')
    
    subplot(1,3,2)
    plot(NoiseCom_Frac(a),sd_Mean(a),[color(a),'*'])
    hold on
    xlabel('fraction of noise which is common')
    ylabel('mean spike distance')

    subplot(1,3,3)
    hold on
    plot(ccX,ccGResMean(a,:),color(a))
end



figure
plot(time,signal)
hold on
plot(time,NoiseCommon + NoiseExc,'r') ;

figure
plot(time,Gexc(1,:))
hold on
plot(time,Ginh(1,:),'r')

figure
plot(time,Gexc(end,:))
hold on
plot(time,Ginh(end,:),'r')

% for igor
stspikes = nans(size(spikes)) ;
stspikes(spikes==1)=1 ;

for a=1:size(Gexc,1) ;
    identifier = ['simulatedGexc',num2str(a)] ;
    ForIgor.(identifier) = Gexc(a,:) ;
    
    identifier = ['simulatedGinh',num2str(a)] ;
    ForIgor.(identifier) = Ginh(a,:) ;

    identifier = ['simulatedI',num2str(a)] ;
    ForIgor.(identifier) = I(a,:) ;
    
end

for a=1:length(NoiseCom_Frac) ;
    identifier = ['simulatedccX',num2str(a)] ;
    ForIgor.(identifier) = ccX ;
    
    identifier = ['simulatedcc',num2str(a)] ;
    ForIgor.(identifier) = mean(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,:)) ;
end

for a=1:size(Gexc,1) ;
    
    identifier = ['simulatedSpikeTrain',num2str(a)] ;
    ForIgor.(identifier) = stspikes(a,:) ;

end   


    
    

        

