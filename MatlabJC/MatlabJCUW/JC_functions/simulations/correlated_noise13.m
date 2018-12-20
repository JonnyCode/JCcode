% this simulation will create conductances to a single cell and test the
% predictability of the variance of their spike number. modified from
% correlated_noise12.m

% JC 6/9/11

% parameters 
Noise_common_Frac_range = [0:.1:1] ; % fraction of common noise
Noise_std = 10 ;

Amp_exc = 40 ;
Amp_inh = 55 ;

SI = .0001 ; % sec sample interval
time_Max = 1 ; % sec
time_peak = .5 ; % sec
time_peakRise = .2 ; % sec

numTrials = 20 ;

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
    end
end

% generate spikes
spikeNumber = nans(numTrials,length(Noise_common_Frac_range)) ; % preallocate
for a = 1:length(Noise_common_Frac_range) ;
    voltage = LIFmodelG(G_exc{a},G_inh{a},1,1/SI,params) ;
    spikeNumber(:,a) = sum(voltage==50,2) ;
end

% means, varaiances and covariance
for a = 1:length(Noise_common_Frac_range) ;
    G_exc_mean(:,a) = mean(G_exc{a},2) ;
    G_inh_mean(:,a) = mean(G_inh{a},2) ;
    
    % variances and covariances
    spikeNumber_var(a) = var(spikeNumber(:,a)) ; 
    
    temp = cov(G_exc_mean(:,a),G_inh_mean(:,a)) ;
    G_exc_mean_var(a) = temp(1,1) ;
    G_inh_mean_var(a) = temp(2,2) ;
    G_mean_cov(a) = temp(1,2) ;
end

% e minus i
for a = 1:length(Noise_common_Frac_range) ;
    EminusI(:,a) = (G_exc_mean(:,a)-G_inh_mean(:,a)) ;
    EminusI_var(a) = var(EminusI(:,a)) ;
    
    EminusI_pred(a) = G_exc_mean_var(a) + G_inh_mean_var(a) - 2*G_mean_cov(a) ;
end

% figure
figure
plot(EminusI_var,EminusI_pred,'*')   

figure
plot(spikeNumber_var,EminusI_pred,'*')
        
        
        
        
        
        
        
    
    