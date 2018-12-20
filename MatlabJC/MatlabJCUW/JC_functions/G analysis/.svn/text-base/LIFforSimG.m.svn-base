function ForIgor = LIFforSimG(Input,Parameters,id,A) ;

% this function will run a leaky integrate and fire model for simultaneous
% conductance recordings and assess spike output and subthresh voltage
% changes between simultaneous recordings and shuffled versions.

% get data
[fp, error] = ITCInitializeAnalysis(1000000, ['~/Data/primate/',Input(A).cellname]);

epochs = str2num(Input(A).(id)) ;
round = 0 ;
for a = 1:3:length(epochs) ; % for each spike epoch
    round = round +1 ;

    [dataExc(round,:), error] = ITCReadEpoch(epochs(a), 0, fp) ;    % get data

    [dataInh(round,:), error] = ITCReadEpoch(epochs(a+1), 0, fp) ;    % get data
    
    [dataAltV(round,:), error] = ITCReadEpoch(epochs(a+2), 0, fp) ;    % get data
    [SI(round), error] = ITCGetSamplingInterval(epochs(a+2), fp); % get sampling interval
    SI(round) = SI(round) * 1e-6; % Sampling interval in sec
end

SI = .0001 ; % if ITC 18

cyclepnts = 100 ; % number of sample points in a cycle, pnts between leaving hold1 and returning (Also gets rid of first cycle)    
FirstAltPnt = (cyclepnts/2) - 1 ; % first sample point you want to plot after begining of step from alternation 
%FirstAltPnt = 30 ;
LastAltPnt = FirstAltPnt ;  % last "                                                                    "

[prePnts, error] = ITCGetStmPrePts(epochs(1), 0, 0, fp) ; % points collected beyond which data is not worth analyzing
[postPnts, error] = ITCGetStmTailPts(epochs(1), 0, 0, fp) ;

samplerate = 1/SI(1) ; % Hz at which data was collected
time = [SI(1):SI(1):SI(1)*length(dataExc(1,:))] ; % time vector in seconds

% WINDOWING ALTERNATING V DATA
Alt_Exc = NaN(size(dataAltV)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold1
Alt_Inh = NaN(size(dataAltV)) ; % make a vector of NaNs that will serve as base for ploting alternating coductances at hold2_

for a = FirstAltPnt:LastAltPnt ;        % for each data point per cycle we want to plot
    Alt_Exc(:,[a:cyclepnts:end]) = dataAltV(:,[a:cyclepnts:end]) ; 
    Alt_Inh(:,[a+cyclepnts/2:cyclepnts:end]) = dataAltV(:,[a+cyclepnts/2:cyclepnts:end]) ; 
end 

% change currents into G and interplate the alt recordings
for a = 1:size(dataExc,1) ; % for each trial
    GExc(a,:) = dataExc(a,:)/-61 - min(dataExc(a,:)/-61) ; % get conductance from stable currrents
    GInh(a,:) = dataInh(a,:)/61 - min(dataInh(a,:)/61) ; 
    GAlt_Exc(a,:) = Alt_Exc(a,:)/-61 - min(Alt_Exc(a,:)/-61) ; % get conductance from alt current
    GAlt_Inh(a,:) = Alt_Inh(a,:)/61 - min(Alt_Inh(a,:)/61) ;

    b = find(isnan(GAlt_Exc(a,:)) == 0) ;     % find all the indices that are not nans
    GAlt_ExcInt(a,:) = interp1(b,GAlt_Exc(a,b),[1:length(GAlt_Exc)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b
    
    b = find(isnan(GAlt_Inh(a,:)) == 0) ;     % find all the indices that are not nans
    GAlt_InhInt(a,:) = interp1(b,GAlt_Inh(a,b),[1:length(GAlt_Inh)],'linear','extrap') ; % interpolate to find values that were not sampled
    clear b

end

% low pass filter to remove electrical crap
GExc_lpf = lowPassFilter(GExc, samplerate, 5000) ; %(signal,samplerate,cutoff frequ (hz))
GInh_lpf = lowPassFilter(GInh, samplerate, 5000) ; 
GAlt_ExcInt_lpf = lowPassFilter(GAlt_ExcInt, samplerate, 5000) ;
GAlt_InhInt_lpf = lowPassFilter(GAlt_InhInt, samplerate, 5000) ;

% high pass filter to remove slow drift 
GExc_hpf = highPassFilter(GExc_lpf, samplerate, 1) ; %(signal,samplerate,cutoff frequ (hz))
GInh_hpf = highPassFilter(GInh_lpf, samplerate, 1) ; 
GAlt_ExcInt_hpf = highPassFilter(GAlt_ExcInt_lpf, samplerate, 1) ;
GAlt_InhInt_hpf = highPassFilter(GAlt_InhInt_lpf, samplerate, 1) ;

% set up to deal with only smultaneous conductances
clearvars -except GAlt_ExcInt_hpf GAlt_InhInt_hpf id A Parameters prePnts postPnts time SI
for a= 1:size(GAlt_ExcInt_hpf,1) ;
    GExc(a,:) = GAlt_ExcInt_hpf(a,:) - mean(GAlt_ExcInt_hpf(a,prePnts:end-postPnts)) ;
    GInh(a,:) = GAlt_InhInt_hpf(a,:) - mean(GAlt_InhInt_hpf(a,prePnts:end-postPnts)) ; 
end
GExc = GExc - min(min(GAlt_ExcInt_hpf(:,prePnts:end-postPnts))) ;
GInh = GInh - min(min(GAlt_InhInt_hpf(:,prePnts:end-postPnts))) ;

    
% run LIF model for simultaneous and shuffled recodings
V_sim = LIFmodelG(GExc,GInh,Parameters.Gleak,10000,Parameters) ; %exc,inh,gleak,samplerate

V_shuf = LIFmodelG(GExc,circshift(GInh,[1,0]),Parameters.Gleak,10000,Parameters) ; %circshift shifts inh G by 1 to shuffle matches 


% get residuals
for a =1:size(GExc,1) ;
    GExcRes(a,:) = GExc(a,:) - mean(GExc) ;
    GInhRes(a,:) = GInh(a,:) - mean(GInh) ;
end

% get cross correlations of residuals
for a = 1:size(GExc,1) ;
    ccRes_sim(a,:) = xcorr(GExcRes(a,prePnts:end-postPnts),GInhRes(a,prePnts:end-postPnts),'coef') ;
    
    if a~=size(GExc,1) ;
        ccRes_shuf(a,:) = xcorr(GExcRes(a,prePnts:end-postPnts),GInhRes(a+1,prePnts:end-postPnts),'coef') ;
    else
         ccRes_shuf(a,:) = xcorr(GExcRes(a,prePnts:end-postPnts),GInhRes(1,prePnts:end-postPnts),'coef') ;
    end
       
end

ccRes_simMean = mean(ccRes_sim) ;
ccRes_shufMean = mean(ccRes_shuf) ;
time_cc = [SI(1)*([1:length(ccRes_simMean)] - (length(ccRes_simMean)+1)/2)] ;


if Parameters.Vthresh>10 ; % if this is no spiking model
% zero mean for each voltage response
V_sim = V_sim - repmat(mean(V_sim,2),1,length(V_sim)) ; 
V_shuf = V_shuf - repmat(mean(V_shuf,2),1,length(V_shuf)) ;

% mean voltage response and std 
V_sim_mean = mean(V_sim) ;
V_shuf_mean = mean(V_shuf) ;

V_sim_std = std(V_sim) ;
V_shuf_std = std(V_shuf) ;

% assess area above voltage distribution from model
for a=1:size(V_sim,1) ;
    hist_simPre(a,:) = hist(V_sim(a,30000:prePnts),[-40:40]) ; % histogram prestim  
    hist_sim(a,:) = hist(V_sim(a,prePnts:end-postPnts),[-40:40]) ; % histogram stim
end

for a=1:size(V_shuf,1) ;
    hist_shufPre(a,:) = hist(V_shuf(a,30000:prePnts),[-40:40]) ; % histogram
    hist_shuf(a,:) = hist(V_shuf(a,prePnts:end-postPnts),[-40:40]) ; % histogram
end

hist_simPre_mean = mean(hist_simPre) ;
hist_shufPre_mean = mean(hist_shufPre) ;

hist_sim_mean = mean(hist_sim) ;
hist_shuf_mean = mean(hist_shuf) ;

% assess time dependent variance of voltage from model
var_simPre = mean(var(V_sim(:,30000:prePnts))) ; 
var_sim = mean(var(V_sim(:,prePnts:end-postPnts))) ; 

var_shufPre = mean(var(V_shuf(:,30000:prePnts))) ; 
var_shuf = mean(var(V_shuf(:,prePnts:end-postPnts))) ;

% assess autocorr of V and its variance from trial to trial (hypothesis:
% auttocorrelation may be better preserved during simultatneous
% conductances)

for a=1:size(V_sim,1) ;
    autoCorr_sim(a,:) = xcov(V_sim(a,prePnts:end-postPnts)) ; %autocorr of signal
    autoCorr_simPre(a,:) = xcov(V_sim(a,30000:prePnts)) ; %autocorr of mean light
end

autoCorr_sim_mean = mean(autoCorr_sim) ;
autoCorr_sim_std = std(autoCorr_sim) ;

autoCorr_simPre_mean = mean(autoCorr_simPre) ;
autoCorr_simPre_std = std(autoCorr_simPre) ;

for a=1:size(V_shuf,1) ;
    autoCorr_shuf(a,:) = xcov(V_shuf(a,prePnts:end-postPnts)) ; %autocorr of signal
    autoCorr_shufPre(a,:) = xcov(V_shuf(a,30000:prePnts)) ; %autocorr of mean light
end

autoCorr_shuf_mean = mean(autoCorr_shuf) ;
autoCorr_shuf_std = std(autoCorr_shuf) ;

autoCorr_shufPre_mean = mean(autoCorr_shufPre) ;
autoCorr_shufPre_std = std(autoCorr_shufPre) ;

% figures
figure
plot(time, V_sim,'b')
hold on
plot(time,V_shuf,'r')
xlabel('time (sec)')
ylabel('voltage (mV)')
legend('simultaneous','shuffled')

figure
plot(time, V_sim_mean,'b')
hold on
plot(time,V_sim_mean+V_sim_std,'b--')
plot(time,V_sim_mean-V_sim_std,'b--')
plot(time,V_shuf_mean,'r')
plot(time,V_shuf_mean+V_shuf_std,'r--')
plot(time,V_shuf_mean-V_shuf_std,'r--')
xlabel('time (sec)')
ylabel('voltage (mV)')
title('mean and std Voltage')

figure
subplot(2,1,1)
plot([-40:40],hist_simPre_mean,'b')
hold on
plot([-40:40],hist_shufPre_mean,'r')
title('prestim')
xlabel('voltage (mV)')
ylabel('mean number of observations pre trial')

subplot(2,1,2)
plot([-40:40],hist_sim_mean,'b')
hold on
plot([-40:40],hist_shuf_mean,'r')
title('prestim')
xlabel('voltage (mV)')
ylabel('mean number of observations pre trial')

figure
subplot(2,1,1)
plot([1:length(autoCorr_simPre)],autoCorr_simPre,'b')
hold on
plot([1:length(autoCorr_shufPre)],autoCorr_shufPre,'r')
xlabel('sample points')
ylabel('voltage (mV)')
title('prestim auto correlation')

subplot(2,1,2)
plot([1:length(autoCorr_sim)],autoCorr_sim,'b')
hold on
plot([1:length(autoCorr_shuf)],autoCorr_shuf,'r')
xlabel('sample points')
ylabel('voltage (mV)')
title('stim auto correlation')

figure
subplot(2,1,1)
plot([1:length(autoCorr_simPre)],autoCorr_simPre_mean,'b')
hold on
plot([1:length(autoCorr_simPre)],autoCorr_simPre_mean+autoCorr_simPre_std,'b--')
plot([1:length(autoCorr_simPre)],autoCorr_simPre_mean-autoCorr_simPre_std,'b--')
plot([1:length(autoCorr_simPre)],autoCorr_shufPre_mean,'r')
plot([1:length(autoCorr_simPre)],autoCorr_shufPre_mean+autoCorr_shufPre_std,'r--')
plot([1:length(autoCorr_simPre)],autoCorr_shufPre_mean-autoCorr_shufPre_std,'r--')

subplot(2,1,2)
plot([1:length(autoCorr_sim)],autoCorr_sim_mean,'b')
hold on
plot([1:length(autoCorr_sim)],autoCorr_sim_mean+autoCorr_sim_std,'b--')
plot([1:length(autoCorr_sim)],autoCorr_sim_mean-autoCorr_sim_std,'b--')
plot([1:length(autoCorr_sim)],autoCorr_shuf_mean,'r')
plot([1:length(autoCorr_sim)],autoCorr_shuf_mean+autoCorr_shuf_std,'r--')
plot([1:length(autoCorr_sim)],autoCorr_shuf_mean-autoCorr_shuf_std,'r--')

else % if this is a spiking model
    for a=1:size(V_sim,1) ; % for each trial during the stimulus get spike times
        spiketimes_sim{a} = time(V_sim(a,prePnts:end-postPnts) ==50) ;
    end
    
    for a=1:size(V_shuf,1) ; % for each trial
        spiketimes_shuf{a} = time(V_shuf(a,prePnts:end-postPnts) ==50) ;
    end
    
    % check spike precision with spike match metric
    [allPairs_sim,meanPairs_sim] = SpikeMatcher(spiketimes_sim, 10.^[-3:.1:0]) ;
    [allPairs_shuf,meanPairs_shuf] = SpikeMatcher(spiketimes_shuf, 10.^[-3:.1:0]) ;
   
    % check spike precision with psth variance
    spiketrain_sim = zeros(size(V_sim,1),length(V_sim)) ;
    spiketrain_sim(V_sim==50) = 1 ;
    spiketrain_shuf = zeros(size(V_shuf,1),length(V_shuf)) ;
    spiketrain_shuf(V_shuf==50) = 1 ;
    
    [sumPSTHvar_sim,psth] = PsthVar(spiketrain_sim,10.^[-3:.1:0],1/SI(1)) ;
    [sumPSTHvar_shuf,psth] = PsthVar(spiketrain_shuf,10.^[-3:.1:0],1/SI(1)) ;
    
    figure
    subplot(3,1,1)
    plot(10.^[-3:.1:0],meanPairs_sim,'b')
    hold on
    plot(10.^[-3:.1:0],meanPairs_shuf,'r')
    h=gca
    set(h,'xscale','log')
    xlabel('max time shift (sec)')
    ylabel('percent spikes able to be matched')
    legend('simulataneous', 'shuffled')
    
    subplot(3,1,2)
    plot(10.^[-3:.1:0],allPairs_sim,'b')
    hold on
    plot(10.^[-3:.1:0],allPairs_shuf,'r')
    h=gca
    set(h,'xscale','log')
    xlabel('max time shift (sec)')
    ylabel('percent spikes able to be matched')
    legend('simulataneous', 'shuffled')
    
    subplot(3,1,3)
    plot(10.^[-3:.1:0],sumPSTHvar_sim,'b')
    hold on
    plot(10.^[-3:.1:0],sumPSTHvar_shuf,'r')
    h=gca
    set(h,'xscale','log')
    xlabel('gaussian variance (sec)')
    ylabel('sum psth var')
    legend('simulataneous', 'shuffled')
    
    
    figure
    subplot(2,1,1)
    for a=1:size(V_sim,1) ;
        Vex = (V_sim(a,:)-min(V_sim(a,:)))./max(V_sim(a,:)-min(V_sim(a,:))) ;
        plot(time, Vex+a,'b')
        hold on
    end
    xlabel('time (sec)')
    ylabel('trial')
  
    subplot(2,1,2)
    for a=1:size(V_shuf,1) ;
        Vex2 = (V_shuf(a,:)-min(V_shuf(a,:)))/max(V_shuf(a,:)-min(V_shuf(a,:))) ;
        hold on
        plot(time,Vex2+a,'r')
    end
    xlabel('time (sec)')
    ylabel('trial')
    
end % if spiking model

figure
subplot(2,1,1)
plot(time_cc,ccRes_simMean,'b')
hold on
plot(time_cc,ccRes_shufMean,'r')
xlabel('time (sec)')
ylabel('xcorr')

subplot(2,1,2)
plot(time_cc,ccRes_sim,'b')
hold on
plot(time_cc,ccRes_shuf,'r')
xlabel('time (sec)')
ylabel('xcorr')

figure
subplot(2,1,1)
plot(time,GExc)
title('exc')

subplot(2,1,2)
plot(time,GInh)
title('inh')

ForIgor.nothing = 0 ;

end % end function


