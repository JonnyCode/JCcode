% this simulation will create conductances to a single cell and test the
% predictability of the variance of their spike number. modified from
% correlated_noise12.m

% JC 6/9/11

% parameters 
Noise_commonEI_Frac_range = [0:.25:1] ; % fraction of common noise between e and i
Noise_commonPair_Frac_range = [1] ; % fraction of common noise between ee and ii
Noise_std = 8 ;

Amp_exc = 40 ;
Amp_inh = 55 ;

SI = .0001 ; % sec sample interval
time_Max = 10 ; % sec
time_peak = .5 ; % sec
time_peakRise = .2 ; % sec

numTrials = 30 ;

params.Einh = -60 ;
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
signal = exp(-((time-time_peak).^2)/(2*time_peakRise^2)) + exp(-((time-time_peak -(time_peakRise*20) ).^2)/(2*time_peakRise^2)) + exp(-((time-time_peak-(time_peakRise*40)).^2)/(2*time_peakRise^2)) ; 

for c = 1:length(Noise_commonPair_Frac_range) ;
    % for each Noise common Fraction
    for a = 1:length(Noise_commonEI_Frac_range) ;
        Noise_commonEI_Frac = Noise_commonEI_Frac_range(a) ;
        Noise_commonPair_Frac = Noise_commonPair_Frac_range(c) ;

        Noise_var = Noise_std^2 ; % total noise in each e1, e2, i1, and i2
        Noise_Pair_var = Noise_var * Noise_commonPair_Frac * (1-Noise_commonEI_Frac) ; % variance of noise common to pair but e or i only (not ei)  
        Noise_EI_var = Noise_var * Noise_commonEI_Frac * (1-Noise_commonPair_Frac) ; % variance of noise common to ei only (not pair) 
        Noise_EIPair_var = Noise_var * Noise_commonPair_Frac * Noise_commonEI_Frac ; % varaince of noise common in ei to pair
        Noise_ind_var = Noise_var - Noise_Pair_var - Noise_EI_var - Noise_EIPair_var ; % variance of noise independant to a conductance in one cell (shared by no other conductance)

        % for every trial
        G_exc{a} = nans(numTrials,length(time)) ; % preallocate
        G_inh{a} = nans(numTrials,length(time)) ;

        for b = 1:numTrials ; 
            Noise_exc1 = normrnd(0,sqrt(Noise_ind_var)) ; % noise independant to a conductance in one cell (shared by no other conductance)
            Noise_inh1 = normrnd(0,sqrt(Noise_ind_var)) ;
            Noise_exc2 = normrnd(0,sqrt(Noise_ind_var)) ;
            Noise_inh2 = normrnd(0,sqrt(Noise_ind_var)) ;        

            Noise_ei1 = normrnd(0,sqrt(Noise_EI_var)) ; % noise common to ei only (not pair)
            Noise_ei2 = normrnd(0,sqrt(Noise_EI_var)) ;

            Noise_ePair = normrnd(0,sqrt(Noise_Pair_var)) ; % noise common to pair but e or i only (not ei) 
            Noise_iPair = normrnd(0,sqrt(Noise_Pair_var)) ;

            Noise_EIPair = normrnd(0,sqrt(Noise_EIPair_var)) ; % noise common in ei to pair

            G_exc1{a}(b,:) = signal*(Amp_exc + Noise_exc1 + Noise_ei1 + Noise_ePair + Noise_EIPair) ;
            G_inh1{a}(b,:) = signal*(Amp_inh + Noise_inh1 + Noise_ei1 + Noise_iPair + Noise_EIPair) ;
            G_exc2{a}(b,:) = signal*(Amp_exc + Noise_exc2 + Noise_ei2 + Noise_ePair + Noise_EIPair) ;
            G_inh2{a}(b,:) = signal*(Amp_inh + Noise_inh2 + Noise_ei2 + Noise_iPair + Noise_EIPair) ;               
        end
    end

    % generate spikes
    spikeNumber = nans(numTrials,length(Noise_commonEI_Frac_range)) ; % preallocate
    mean_voltage = 0 ;
    for a = 1:length(Noise_commonEI_Frac_range) ;
        voltage1 = LIFmodelG(G_exc1{a},G_inh1{a},1,1/SI,params) ;
        voltage2 = LIFmodelG(G_exc2{a},G_inh2{a},1,1/SI,params) ;
        
        mean_voltage = mean(voltage1(:))/length(Noise_commonEI_Frac_range) + mean_voltage ;  

        spikeNumber1(:,a) = sum(voltage1==50,2) ;
        spikeNumber2(:,a) = sum(voltage2==50,2) ;
    end

    % means, varaiances and covariance
    for a = 1:length(Noise_commonEI_Frac_range) ;
        G_exc1_mean(:,a) = mean(G_exc1{a},2) ;
        G_inh1_mean(:,a) = mean(G_inh1{a},2) ;

        G_exc2_mean(:,a) = mean(G_exc2{a},2) ;
        G_inh2_mean(:,a) = mean(G_inh2{a},2) ;

        % variances and covariances
        temp = cov(spikeNumber1(:,a),spikeNumber2(:,a)) ;
        spikeNumber1_var(a) = temp(1,1); 
        spikeNumber2_var(a) = temp(2,2) ;
        spikeNumber_cov(a) = temp(1,2) ;

        spikeNumber_corrCoef(a,c) = spikeNumber_cov(a) / sqrt(spikeNumber1_var(a)*spikeNumber2_var(a)) ;

        temp = cov(G_exc1_mean(:,a),G_exc2_mean(:,a)) ;
        G_exc1_mean_var(a) = temp(1,1) ;
        G_exc2_mean_var(a) = temp(2,2) ;
        G_exc_mean_cov(a) = temp(1,2) ;

        temp = cov(G_inh1_mean(:,a),G_inh2_mean(:,a)) ;
        G_inh1_mean_var(a) = temp(1,1) ;
        G_inh2_mean_var(a) = temp(2,2) ;
        G_inh_mean_cov(a) = temp(1,2) ;

        temp = cov(G_exc1_mean(:,a),G_inh2_mean(:,a)) ;
        G_exc1inh2_mean_cov(a) = temp(1,2) ;

        temp = cov(G_exc2_mean(:,a),G_inh1_mean(:,a)) ;
        G_exc2inh1_mean_cov(a) = temp(1,2) ; 

        temp = cov(G_exc1_mean(:,a),G_inh1_mean(:,a)) ;
        G_exc1inh1_mean_cov(a) = temp(1,2) ;

        temp = cov(G_exc2_mean(:,a),G_inh2_mean(:,a)) ;
        G_exc2inh2_mean_cov(a) = temp(1,2) ;
    end

    % e minus i
    for a = 1:length(Noise_commonEI_Frac_range) ;
        EminusI1(:,a) = (G_exc1_mean(:,a)-G_inh1_mean(:,a)) ;
        EminusI2(:,a) = (G_exc2_mean(:,a)-G_inh2_mean(:,a)) ;

        temp = cov(EminusI1(:,a),EminusI2(:,a)) ;
        EminusI1_var(a) = temp(1,1) ;
        EminusI2_var(a) = temp(2,2) ;
        EminusI_cov(a) = temp(1,2) ;

        EminusI_corrCoef(a,c) = EminusI_cov(a)/sqrt(EminusI1_var(a)*EminusI2_var(a)) ;

        EminusI1_var_pred(a) = (G_exc1_mean_var(a) + G_inh1_mean_var(a) - 2*G_exc1inh1_mean_cov(a)) ;
        EminusI2_var_pred(a) = (G_exc2_mean_var(a) + G_inh2_mean_var(a) - 2*G_exc2inh2_mean_cov(a)) ;

        EminusI_corrCoef_pred(a,c) = (G_exc_mean_cov(a) + G_inh_mean_cov(a) - G_exc1inh2_mean_cov(a) - G_exc2inh1_mean_cov(a))/...
            sqrt(EminusI1_var_pred(a)*EminusI2_var_pred(a)) ;
        
        % testing for best alpha
%         v_range = round(mean_voltage)-10:2:params.Vthresh+10 ;
%         for alphaRound = 1:length(v_range) ;
%             alpha = abs(params.Einh - v_range(alphaRound))/abs(params.Eexc - v_range(alphaRound)) ;
%             EminusI_corrCoef_pred_plusAlpha{alphaRound}(a,c) = (G_exc_mean_cov(a) + alpha^2*G_inh_mean_cov(a) - alpha*G_exc1inh2_mean_cov(a) - alpha*G_exc2inh1_mean_cov(a))/...
%                 sqrt((G_exc1_mean_var(a) + alpha^2*G_inh1_mean_var(a) - alpha*2*G_exc1inh1_mean_cov(a))*(G_exc2_mean_var(a) + alpha^2*G_inh2_mean_var(a) - alpha*2*G_exc2inh2_mean_cov(a))) ;
%         end
        alpha = abs(params.Einh - params.Vthresh)/abs(params.Eexc - params.Vthresh) ;
        EminusI_corrCoef_pred_plusAlpha(a,c) = (G_exc_mean_cov(a) + alpha^2*G_inh_mean_cov(a) - alpha*G_exc1inh2_mean_cov(a) - alpha*G_exc2inh1_mean_cov(a))/...
            sqrt((G_exc1_mean_var(a) + alpha^2*G_inh1_mean_var(a) - alpha*2*G_exc1inh1_mean_cov(a))*(G_exc2_mean_var(a) + alpha^2*G_inh2_mean_var(a) - alpha*2*G_exc2inh2_mean_cov(a))) ;
    end
end

% for a=1:alphaRound ;
%     temp = corrcoef(EminusI_corrCoef_pred_plusAlpha{a}(:),spikeNumber_corrCoef(:)) ;
%     cc(a) = temp(1,2) ;
%     
%     mse(a) = mean((EminusI_corrCoef_pred_plusAlpha{a}(:) - spikeNumber_corrCoef(:)).^2) ;
% end
% figure
% plot(v_range,mse)



% figure
figure
plot(spikeNumber1(:,1))
hold on
plot(spikeNumber2(:,2),'r')

figure
plot(EminusI1_var_pred,EminusI1_var,'*')
hold on
plot(EminusI2_var_pred,EminusI2_var,'ro')

figure
plot(EminusI_corrCoef_pred,EminusI_corrCoef,'*')   


figure
color = ['b','y','r','g','c'] ;
for a = 1:length(Noise_commonEI_Frac_range) ;
    plot(EminusI_corrCoef_pred(a,:),spikeNumber_corrCoef(a,:),[color(a),'*'])
    hold on
    plot(EminusI_corrCoef_pred_plusAlpha(a,:),spikeNumber_corrCoef(a,:),[color(a),'o'])
end
xlabel('predicted')
ylabel('actual spike corr')
        
figure
color = ['b','y','r','g','c'] ;
for a = 1:length(Noise_commonEI_Frac_range) ;
    plot(Noise_commonPair_Frac_range, spikeNumber_corrCoef(a,:),[color(a),'*'])
    hold on
end
xlabel('Noise commonPair Frac')
ylabel('actual spike corr') 

figure
color = ['b','y','r','g','c'] ;
for a = 1:length(Noise_commonEI_Frac_range) ;
    plot(Noise_commonEI_Frac_range, spikeNumber_corrCoef(:,a),[color(a),'*'])
    hold on
end        
xlabel('Noise commonEI Frac')
ylabel('actual spike corr')         
        
        
    
    