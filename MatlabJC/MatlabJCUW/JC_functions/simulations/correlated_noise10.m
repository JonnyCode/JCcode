% this simuluation is investigating nonshuffled (with noise correlations)
% vs nonshuffled as in experiment

%parameters
NoiseCom_Frac = 1 ; % common noise fraction of total
NoiseTot_Var = .5 ; % var of total noise in each conductance
NumTrials = 5 ; % number of repeats 

params.Einh = -80 ;
params.Eexc = 0 ;
params.Eleak = -60 ;
params.cap = .06 ;%nF
params.Vrest = -60;
params.Vthresh = 1000 ;
params.AbsRef  = .002 ; %sec
params.RelRefTau = .008 ;%sec
params.RelRefAmp = 4 ; %mV

MaxDeltaTvector = 10.^[-3:.1:0] ; % sec vector of times shifts for max delta t

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
signal = normrnd(.5,2,1,lt) ;     % signal
signal = lowPassFilter(signal,1/si,60) ;    % low pass at 60hz

% signal_exc = conv(GexcPreLN,excLF) ;
% se = signal_exc*NFgainExc - NFgainExc*NFthreshExc ;
% se(signal_exc<NFthreshExc) = 0 ;
% se(se>NFmaxExc) = NFmaxExc ;
% signal_exc = se(1:lt) ;
% 
% signal_inh = conv(GinhPreLN,inhLF) ;
% si = signal_inh*NFgainInh - NFgainInh*NFthreshInh ;
% si(signal_inh<NFthreshInh) = 0 ;
% si(si>NFmaxInh) = NFmaxInh ;
% signal_inh = si(1:lt) ;



NoiseCom_Std = sqrt(NoiseTot_Var * NoiseCom_Frac) ; % standard deviation of common noise
NoiseInd_Std = sqrt(NoiseTot_Var - NoiseCom_Std^2) ; % standard deviation of independant noise

for b=1:NumTrials ; % for each trial

    NoiseCommon = normrnd(0,NoiseCom_Std,1,lt) ; % create noise for this trial
    NoiseExc = normrnd(0,NoiseInd_Std,1,lt) ;
    NoiseInh = normrnd(0,NoiseInd_Std,1,lt) ;

    NoiseCommon = lowPassFilter(NoiseCommon,1/si,80) ; % low pass at 80hz
    NoiseExc = lowPassFilter(NoiseExc,1/si,80) ; % low pass at 80hz
    NoiseInh = lowPassFilter(NoiseInh,1/si,80) ; % low pass at 80hz

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

    Gexc(b,:) = GexcPre(1,1:lt) ; % cut off convolution extra
    Ginh(b,:) = GinhPre(1,1:lt) ;
    Gleak = 0 ; % leak conductance

    V_sim(b,:) = LIFmodelG(Gexc(b,:),Ginh(b,:),Gleak,1/si,params) ; % run LIF model for simultaneous g 
    spiketimes_sim{b} = time(V_sim(b,:)==50) ; % spike times for simultaneous 
    
    Isim(b,:) = Gexc(b,:).*(V_sim(b,:)-params.Eexc) + Ginh(b,:).*(V_sim(b,:)-params.Einh) ;
end
GexcMean = mean(Gexc,1) ; % mean Gexc 
GinhMean = mean(Ginh,1) ; % mean Gexc 

GinhShuf = circshift(Ginh,[4,0]) ;
V_shuf = LIFmodelG(Gexc,GinhShuf,Gleak,1/si,params) ; %circshift shifts inh G by 1 point to shuffle matches 
Ishuf = Gexc.*(V_shuf-params.Eexc) + GinhShuf.*(V_shuf-params.Einh) ;

for a = 1:size(V_shuf,1) ;
    spiketimes_shuf{a} = time(V_shuf(a,:)==50) ;
end

% get cross correlations of estimated noise (residuals), estimate and predict total current noise, measure voltage noise (only applicable for nonspiking) 
lcc = length(GexcMean)*2-1 ; % length of cc
ccX = ([1:lcc] - (lcc+1)/2)*si ; % x values of cc
ccZeroi = (lcc+1)/2 ; % index of zero time lag

index=[1:size(Gexc,1)] ;
for b = 1:size(Gexc,1) ; % for each trial
    
    Gres_exc(b,:) = Gexc(b,:)-GexcMean;
    Gres_inh(b,:) = Ginh(b,:)-GinhMean ;
    
    ccG_sim(b,:) = xcov(Gexc(b,:),Ginh(b,:)) ; % cross corr
    ccG_shuf(b,:) = xcov(Gexc(b,:),GinhShuf(b,:)) ;

    acGres_exc(b,:) = xcov(Gres_exc(b,:)) ;
    acGres_inh(b,:) = xcov(Gres_inh(b,:)) ;

    acG_exc(b,:) = xcov(Gexc(b,:)) ;
    acG_inh(b,:) = xcov(Ginh(b,:)) ;
    
end
ccMeanG = xcov(GexcMean,GinhMean) ;

varN = sqrt(mean(var(Gexc))*mean(var(Ginh))) ;
varN2 = sqrt(max(mean(acG_exc))*max(mean(acG_inh)))-max(ccMeanG) ;
varN3 = sqrt(max(mean(acGres_exc))*max(mean(acGres_inh))) ;

ccG_sc_sim = (mean(ccG_sim) - ccMeanG)./varN3 ; % shuffle corrected cross corr
ccG_sc_shuf = (mean(ccG_shuf) - ccMeanG)./(varG - max(ccMeanG)) ;

    InoiseVarEst_sim(b) = var(GexcRes(b,:)*(-45 - params.Eexc)+GinhRes(b,:)*(-45 - params.Einh)) ; % estimate variance of current noise
    InoiseVarPred_sim(b) = mean(var(GexcRes(b,:)*(-45 - params.Eexc),[],2)) + mean(var(GinhRes(b,:)*(-45 - params.Einh),[],2)) + 2*ccIRes_sim(b,ccZeroi)/lt ; % theoretical prediction based on cross correlation
end
ccGRes_simMean = mean(ccGRes_sim,1) ; % cross correlation mean
ccIRes_simMean = mean(ccIRes_sim,1) ; % cross correlation mean

ccGResCoef_simMean = mean(ccGResCoef_sim) ;
ccGResCoef_shufMean = mean(ccGResCoef_shuf) ;


% check spike precision with spike match metric
[spkM_sim,spkMmean_sim] = SpikeMatcher(spiketimes_sim, MaxDeltaTvector) ;
[spkM_shuf,spkMmean_shuf] = SpikeMatcher(spiketimes_shuf, MaxDeltaTvector) ;  

figure
for a = 1:size(V_shuf,1)
    if a<=size(V_sim,1) ;
        Vex = (V_sim(a,:)-min(V_sim(a,:)))/max(V_sim(a,:)-min(V_sim(a,:))) ;
        plot(time,a+Vex,'b')
        hold on
    end 
    Vex2 = (V_shuf(a,:)-min(V_shuf(a,:)))/max(V_shuf(a,:)-min(V_shuf(a,:))) ;
    plot(time,a+size(V_sim,1)+Vex2,'r')
end

figure
subplot(4,1,1)
plot(ccGRes_sim(:,ccZeroi),InoiseVarEst_sim,'b*')
hold on
plot(ccGRes_sim(:,ccZeroi),InoiseVarPred_sim,'bo')
plot(ccGRes_shuf(:,ccZeroi),InoiseVarEst_shuf,'r*')
plot(ccGRes_shuf(:,ccZeroi),InoiseVarPred_shuf,'ro')
xlabel('cc coef (0)')
ylabel('I noise variance estimate')
legend('actual estimate', 'predicted','shuffled')

subplot(4,1,2)
plot(MaxDeltaTvector,spkMmean_sim,'b')
hold on
plot(MaxDeltaTvector,spkMmean_shuf,'r')
xlabel('time shift (sec)')
ylabel('percent of spikes that can be matched')
h=gca ;
set(h,'xscale','log')

subplot(4,1,3)
plot(MaxDeltaTvector,spkM_sim,'b')
hold on
plot(MaxDeltaTvector,spkM_shuf,'r')
xlabel('time shift (sec)')
ylabel('percent of spikes that can be matched')
h=gca ;
set(h,'xscale','log')

subplot(4,1,4)
plot(ccX,ccGRes_simMean,'b')
hold on
plot(ccX,ccGRes_shufMean,'r')
xlabel('lag (sec)')
ylabel('cc (xcorr)')


subplot(4,1,2)
plot(ccGResCoef(:,ccZeroi),-bITheory*ccGResCoef(:,ccZeroi)+bITheory,'k--') 
