% this simulation will assess effect of correlated noise of varaibility of
% the spike output.


%parameters
NoiseCom_Frac = 1 ; % range of common noise fraction of total
NoiseTot_Frac = [.05,.1,.2,.4] ; % fraction of total noise in each conductance
TotVar = 2 ; %variance of input (signal + noise)
NumTrials = 10 ; % number of repeats 
signal_mu = 0 ;

params.Einh = -80 ;
params.Eexc = 0 ;
params.Eleak = -60 ;
params.cap = .04 ;%nF
params.Vthresh = -45 ;
params.AbsRef  = .002 ; %sec
params.RelRefTau = .008 ;%sec
params.RelRefAmp = 4 ; %mV
params.Gleak = 1 ;
params.SampleRate = 10000 ;

stdGaussian = 10.^[-3:.1:-1] ; % std of gaussians to be convolved w/ spike trains

si = .0001 ; % sample interval (sec)
time = [0:si:10] ; % sec

lt = length(time) ;

excLF = simFilter([0:si:.3],.040,.009,.1,.10,.05,.02) ; % linear filters(time,tpeak,peakRise,peakAmp,ttrough,troughDecay,troughAmp)
inhLF = simFilter([0:si:.3],.05,.01,.1,.11,.06,.02) ;

NFgainExc = 1 ;
NFgainInh = 1 ;

NFthreshExc = 0 ;
NFthreshInh = 0 ;

NFmaxExc = 300 ;
NFmaxInh = 300 ;


% set signal
signal_orig = normrnd(signal_mu,TotVar,1,lt) ;     % signal
signal = lowPassFilter(signal_orig,1/si,60) ;    % low pass at 60hz

Trial = 0 ;
spikes = zeros(NumTrials*length(NoiseCom_Frac),lt) ;

for a = 1:length(NoiseTot_Frac) ; % for each new change in noise amplitude   

    NoiseCom_Std = sqrt(TotVar*NoiseTot_Frac(a) * NoiseCom_Frac) ; % standard deviation of common noise
    NoiseInd_Std = sqrt(TotVar*NoiseTot_Frac(a) - NoiseCom_Std^2) ; % standard deviation of independant noise
    
    for b=1:NumTrials ; % for each trial
        Trial = Trial + 1 ;
        
        NoiseCommon = normrnd(0,NoiseCom_Std,1,lt) ; % create noise for this trial
        NoiseExc = normrnd(0,NoiseInd_Std,1,lt) ;
        NoiseInh = normrnd(0,NoiseInd_Std,1,lt) ;
        
        NoiseCommon = lowPassFilter(NoiseCommon,1/si,1000) ; % low pass at 1000hz
        NoiseExc = lowPassFilter(NoiseExc,1/si,1000) ; % low pass at 1000hz
        NoiseInh = lowPassFilter(NoiseInh,1/si,1000) ; % low pass at 1000hz
        
        GnoiseExcPreLN = NoiseCommon + NoiseExc ; % add noise
        GnoiseInhPreLN = NoiseCommon + NoiseInh ;

        % pass signal and noise through linearfilters

        GnoiseExcPreN = conv(GnoiseExcPreLN,excLF) ; % convolve noise w/ linear filter
        GnoiseInhPreN = conv(GnoiseInhPreLN,inhLF) ; 
        
        GsignalExcPreN = conv(signal,excLF) ; % convolve signal w/ linear filter
        GsignalInhPreN = conv(signal,inhLF) ;
        
        % add noise and signal and set fraction of variance 
        NoiseAmpFactor_exc = sqrt((TotVar*NoiseTot_Frac(a))/(var(GnoiseExcPreN)+10^-100)) ;
        NoiseAmpFactor_inh = sqrt((TotVar*NoiseTot_Frac(a))/(var(GnoiseInhPreN)+10^-100)) ;        
        
        SignalAmpFactor_exc = sqrt((TotVar-(TotVar*NoiseTot_Frac(a)))/(var(GsignalExcPreN)+10^-100)) ;
        SignalAmpFactor_inh = sqrt((TotVar-(TotVar*NoiseTot_Frac(a)))/(var(GsignalInhPreN)+10^-100)) ;        
        
        GexcPreN = GsignalExcPreN*SignalAmpFactor_exc + GnoiseExcPreN*NoiseAmpFactor_exc ;
        GinhPreN = GsignalInhPreN*SignalAmpFactor_inh + GnoiseInhPreN*NoiseAmpFactor_inh ;
        
        % pass through nonlinear filter(threshold, gain, saturate)
        GexcPre = GexcPreN*NFgainExc - NFgainExc*NFthreshExc ; % gain
        GinhPre = GinhPreN*NFgainInh - NFgainExc*NFthreshInh ;
        
        GexcPre(GexcPreN<NFthreshExc) = 0 ; %threshold
        GinhPre(GinhPreN<NFthreshInh) = 0 ; %threshold
        
        GexcPre(GexcPre>NFmaxExc) = NFmaxExc ; % saturate
        GinhPre(GinhPre>NFmaxInh) = NFmaxInh ; 
        
        Gexc(Trial,:) = GexcPre(1,1:lt) ; % cut off convolution extra
        Ginh(Trial,:) = GinhPre(1,1:lt) ;
        Gleak = 0 ; % leak conductance
                
        V_trace(Trial,:) = LIFmodelG_c(Gexc(Trial,:),Ginh(Trial,:),params) ; % run LIF model
        spikes(Trial,V_trace(Trial,:)==50) = 1 ; % spike train
        spiketimes{Trial} = time(V_trace(Trial,:)==50) ; % spike times     
    end
    GexcMean(a,:) = mean(Gexc(Trial-NumTrials+1:Trial,:)) ; % mean Gexc for this particular set of trials
    GinhMean(a,:) = mean(Ginh(Trial-NumTrials+1:Trial,:)) ; % mean Ginh for this particular set of trials
    V_traceMean(a,:) = mean(V_trace(Trial-NumTrials+1:Trial,:)) ; % mean voltage for this particular set of trials
end

for a=1:length(NoiseTot_Frac)
    firingRate(a) = sum(sum(spikes(a*NumTrials-NumTrials+1:a*NumTrials,:)))/(max(time)*NumTrials);
end


% get cross correlations of estimated noise (residuals), estimate and predict total current noise, measure voltage noise (only applicable for nonspiking) 
lcc = length(GexcMean)*2-1 ; % length of cc
ccX = ([1:lcc] - (lcc+1)/2)*si ; % x values of cc
ccZeroi = (lcc+1)/2 ; % index of zero time lag
for a = 1:length(NoiseTot_Frac) ; % for each new change in common noise amplitude
    for b = a*NumTrials-NumTrials+1:a*NumTrials ; % for each trial         
        GexcRes(b,:) = Gexc(b,:) - GexcMean(a,:) ; % residual
        GinhRes(b,:) = Ginh(b,:) - GinhMean(a,:) ; % residual
        VRes(b,:) = V_trace(b,:) - V_traceMean(a,:) ; % residual
        
        ccGRes(b,:) = xcorr(GexcRes(b,:),GinhRes(b,:)) ; % cross corr of residuals
        ccGResCoef(b,:) = xcorr(GexcRes(b,:),GinhRes(b,:),'coef') ; % cross corr of residuals       
        ccIRes(b,:) = xcorr(GexcRes(b,:)*(-45 - params.Eexc),GinhRes(b,:)*(-45 - params.Einh)) ; % cross corr of current residuals
        
        InoiseVarEst(b) = var(GexcRes(b,:)*(-45 - params.Eexc)+GinhRes(b,:)*(-45 - params.Einh)) ; % estimate variance of current noise
        InoiseVarPred(b) = mean(var(GexcRes(b,:)*(-45 - params.Eexc),[],2)) + mean(var(GinhRes(b,:)*(-45 - params.Einh),[],2)) + 2*ccIRes(b,ccZeroi)/lt ; % theoretical prediction based on cross correlation
        VnoiseVar(b) = var(VRes(b,:)) ; % variance of voltage
    end
    ccGResCoefMean(a,:) = mean(ccGResCoef(a*NumTrials-NumTrials+1:a*NumTrials,:)) ;
    ccGResMean(a,:) = mean(ccGRes(a*NumTrials-NumTrials+1:a*NumTrials,:)) ; % cross correlation mean
    ccIResMean(a,:) = mean(ccIRes(a*NumTrials-NumTrials+1:a*NumTrials,:)) ; % cross correlation mean
    
    GexcResVarMean(a) = mean(var(GexcRes(a*NumTrials-NumTrials+1:a*NumTrials,:),[],2)) ; % mean variance of gexc noise
    GinhResVarMean(a) = mean(var(GinhRes(a*NumTrials-NumTrials+1:a*NumTrials,:),[],2)) ; % mean variance of ginh noise
    
    GexcVarMean(a) = mean(var(Gexc(a*NumTrials-NumTrials+1:a*NumTrials,:),[],2)) ; % mean variance of gexc
    GinhVarMean(a) = mean(var(Ginh(a*NumTrials-NumTrials+1:a*NumTrials,:),[],2)) ;
    
    InoiseVarEstMean(a) = mean(InoiseVarEst(a*NumTrials-NumTrials+1:a*NumTrials)) ; % mean estimate of current noise variance
    InoiseVarPredMean(a) = mean(InoiseVarPred(a*NumTrials-NumTrials+1:a*NumTrials)) ; % mean estimate of current noise variance
    VnoiseVarMean(a) = mean(VnoiseVar(a*NumTrials-NumTrials+1:a*NumTrials)) ; % mean variance of voltage
end

ccGResCoefPeak = max(ccGResCoef,[],2) ; % peak of cc coef
ccGResCoefTrough = min(ccGResCoef,[],2) ; % trough of cc coef

NoiseTot_FracActual = (GexcResVarMean+GinhResVarMean)./(GexcVarMean+GinhVarMean) ; % the fraction of the total variance that is residuals

% % check spike precision with spike match metric
% for a = 1:length(NoiseCom_Frac) ; % for each new change in common noise amplitude
%     [matrix,meanvector] = SpikeMatcher(spiketimes(a*NumTrials-NumTrials+1:a*NumTrials), stdGaussian) ;
%     SpkMIndPairs{a} = matrix ;
%     SpkMMean(a,:) = meanvector ; 
% end

% check spike precision with psth var metric
for a = 1:length(NoiseTot_Frac) ;
    [meanPSTHsnr(a,:),sumPSTHvar(a,:),psth{a}] = PsthVar(spikes(a*NumTrials-NumTrials+1:a*NumTrials,:),stdGaussian,1/si);
    %sumPSTHvar(a,:) = temp ;
end

% find gaussian width that provides a signal to noise ratio of 1
for a = 1:length(NoiseTot_Frac) ;
    minIntegrationTime(a) = interp1(meanPSTHsnr(a,:),stdGaussian,1) ;
end

% calculate infomation per spike and per second with psth (equation taken from Sincich et al 2009)
for a = 1:length(NoiseTot_Var) ;
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

% fit a line to Inoise as a function of cc(0) 

Color = [1-NoiseCom_Frac,0,NoiseCom_Frac] ;

figure
rasterPlot(spiketimes,[0:50],'r','r')
xlabel('time (sec)')
ylabel('trial')

figure
plot(stdGaussian,meanPSTHsnr,'*','color',Color)
hold on
xlabel('stdGuassian')
ylabel('meanPSTHsnr')
h=gca ;
set(h,'xscale','log')

figure
plot(NoiseTot_FracActual,minIntegrationTime,'*','color',Color)
hold on
xlabel('fraction of input variance that is noise')
ylabel('gaussian std to achiev snr of 1')

figure
plot(InoiseVarEstMean,minIntegrationTime,'*','color',Color)
hold on
xlabel('I noise input var')
ylabel('gaussian std to achiev snr of 1')

figure
plot(ccX,ccIResMean,'color',Color)
hold on
xlabel('lag (sec)')
ylabel('cross correlation')

figure
plot(ccX,ccGResCoefMean,'color',Color)
hold on
xlabel('lag (sec)')
ylabel('cross correlation coef')

figure
plot(InoiseVarEst,InoiseVarPred,'*','color',Color)
hold on
xlabel('I noise var est')
ylabel('I noise var pred')

figure
plot(NoiseTot_FracActual,firingRate,'*','color',Color) ;
hold on
xlabel('NoiseTot_FracActual')
ylabel('firingRate')

% ForIgor
% forIgor.time = time ;
% for a=1:NumTrials, iden=['spikes',num2str(a)], forIgor.(iden) = spikes(a,:)+a-1 ; end 
% forIgor.stdGaussian = stdGaussian ;
% for a=1:length(NoiseTot_Frac),iden = ['meanPSTHsnr',num2str(a)], forIgor.(iden) = meanPSTHsnr(a,:) ; end 
% forIgor.minIntTime = minIntegrationTime ; 
% forIgor.FracInputResVar = NoiseTot_FracActual ;
% forIgor.ccX = ccX ;
for a=1:length(NoiseTot_Frac), iden = ['ccGResCoefMean',num2str(a)], forIgor.(iden) = ccGResCoefMean(a,:) ; end
% forIgor.InoiseVarPredMean = InoiseVarPredMean ;
% forIgor.InoiseVarEstMean = InoiseVarEstMean ;
% forIgor.InoiseVarPred = InoiseVarPred ;
% forIgor.InoiseVarEst = InoiseVarEst ;

cd ~/Desktop
exportStructToHDF5(forIgor, 'noComm.h5', '/')