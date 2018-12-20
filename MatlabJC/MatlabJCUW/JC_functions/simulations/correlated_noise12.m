% this simulation will create conductances to a single cell and test the
% predictability of the variance of their spike number.

% 6/2/11

% parameters 
Noise_common_Frac_range = [0:.1:1] ; % fraction of common noise
Noise_std = 10 ;

Amp_exc = 40 ;
Amp_inh = 55 ;

SI = .0001 ; % sec sample interval
time_Max = 1 ; % sec
time_peak = .5 ; % sec
time_peakRise = .2 ; % sec

numTrials = 100 ;

params.Einh = -80 ;
params.Eexc = 0 ;
params.Eleak = -60 ;
params.cap = .04 ;%nF
params.Vthresh = -45 ;
params.AbsRef  = .002 ; %sec
params.RelRefTau = .008 ;%sec
params.RelRefAmp = 4 ; %mV
params.Gleak = 1 ;
params.Vrest = -60 ;

% create exc and inhibitory conductances (a discreet exc and inh g event difined by amplitude)
time = [0:SI:time_Max] ;
signal = exp(-((time-time_peak).^2)/(2*time_peakRise^2)) ; 

% for each Noise common Fraction
for a = 1:length(Noise_common_Frac_range) ;
    Noise_common_Frac = Noise_common_Frac_range(a) ;
    
    Noise_ind_std = sqrt(Noise_std^2 - Noise_std^2*Noise_common_Frac) ;
    Noise_common_std = sqrt(Noise_std^2*Noise_common_Frac) ;
    
    % for every trial
    G_exc{a} = nans(numTrials,length(time)) ; % preallocate
    G_inh{a} = nans(numTrials,length(time)) ;
    
    for b = 1:numTrials ; 
        Noise_exc = normrnd(0,Noise_ind_std) ;
        Noise_inh = normrnd(0,Noise_ind_std) ;
        Noise_common = normrnd(0,Noise_common_std) ;
    
        G_exc{a}(b,:) = signal*(Amp_exc + Noise_exc + Noise_common) ;
        G_inh{a}(b,:) = signal*(Amp_inh + Noise_inh + Noise_common) ;
        
        EoverI{a}(b,:) = G_exc{a}(b,:)./G_inh{a}(b,:) ;
        EminusI{a}(b,:) = G_exc{a}(b,:) - G_inh{a}(b,:) ;
    end
end

% generate spikes
spikeNumber = nans(numTrials,length(Noise_common_Frac_range)) ; % preallocate
for a = 1:length(Noise_common_Frac_range) ;
    voltage = LIFmodelG(G_exc{a},G_inh{a},1,1/SI,params) ;
    spikeNumber(:,a) = sum(voltage==50,2) ;
end

% variance of spikes and e/i
spikeNumber_var = var(spikeNumber,0,1) ;
for a = 1:length(Noise_common_Frac_range) ;
    EoverI_var_mean(a) = mean(var(EoverI{a},0,1)) ;
    EminusI_var_mean(a) = mean(var(EminusI{a},0,1)) ;
end

EoverI_var_mean_norm = EoverI_var_mean/max(EoverI_var_mean) ;
EminusI_var_mean_norm = EminusI_var_mean/max(EminusI_var_mean) ;    

% figure
figure
plot(mean(spikeNumber))

figure
plot(spikeNumber_var,EoverI_var_mean_norm,'r*')
hold on
plot(spikeNumber_var,EminusI_var_mean_norm,'b*')


% means, squared means, variances, residuals and noise covarrelations
for a = 1:length(Noise_common_Frac_range) ;
    spikeNumber_mean = mean(spikeNumber) ; % means
    G_exc_mean{a} = mean(G_exc{a}) ;
    G_inh_mean{a} = mean(G_exc{a}) ;

    G_exc_squared_mean{a} = mean(G_exc{a}.^2) ; % conductances squared
    G_inh_squared_mean{a} = mean(G_exc{a}.^2) ;
    
    spikeNumber_var = var(spikeNumber,0,2) ; % variances
    G_exc_var{a} = var(G_exc{a},0,2) ;
    G_inh_var{a} = var(G_exc{a},0,2) ;
    
    for b = 1:numTrials ;
        G_exc_res{a}(b,:) = G_exc{a}(b,:) - G_exc_mean{a} ; % residuals
        G_inh_res{a}(b,:) = G_inh{a}(b,:) - G_inh_mean{a} ;
    
        G_res_corr{a}(b,:) = xcorr(G_exc_res{a}(b,:),G_inh_res{a}(b,:),'unbiased') ; 
        G_exc_res_ac{a}(b,:) = xcorr(G_exc_res{a}(b,:),'unbiased') ;
        G_inh_res_ac{a}(b,:) = xcorr(G_inh_res{a}(b,:),'unbiased') ;
        
        G_exc_res_var(a,b) = var(G_exc_res{a}(b,:)) ;
        G_inh_res_var(a,b) = var(G_inh_res{a}(b,:)) ;
        
        
    end
    
    G_exc_res_var_mean(a) = mean(G_exc_res_var(a,:)) ;
    G_inh_res_var_mean(a) = mean(G_inh_res_var(a,:)) ;
    
    G_res_corr_mean_peak(a) = CCpeakFinder(mean(G_res_corr{a})) ;
    G_exc_res_ac_mean_peak(a) = CCpeakFinder(mean(G_exc_res_ac{a})) ;
    G_inh_res_ac_mean_peak(a) = CCpeakFinder(mean(G_inh_res_ac{a})) ;
end

% e minus i prediction
for a = 1:length(Noise_common_Frac_range) ;
    EminusI_pred(a) = G_exc_res_ac_mean_peak(a) + G_inh_res_ac_mean_peak(a) - 2*G_res_corr_mean_peak(a) ;
end

% figure
figure
plot(EminusI_var_mean,EminusI_pred,'*')

% e over i prediction
% for a = 1:length(Noise_common_Frac_range) ;
%     EoverI_pred(a) = var



%     
        
        
        
        
        
        
        
        
        
    
    