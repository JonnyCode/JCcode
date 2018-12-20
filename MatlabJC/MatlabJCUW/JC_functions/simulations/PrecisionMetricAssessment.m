% the script will compare different methods of assessing spike variability

%parameters
total_var = 70 ;
signal_mean = 0 ;
Gleak = 8 ;

GNoise_frac = [.16,.2] ; % signal var to noise var ratio
NumTrials = 3 ; % number of repeats 

si = .0001 ; % sample interval (sec)
time = [0:si:10] ; % sec
lt = length(time) ;

params.Einh = -80 ;
params.Eexc = 0 ;
params.Eleak = -60 ;
params.cap = .06 ;%nF
params.Vthresh = -45 ;
params.AbsRef  = .002 ; %sec
params.RelRefTau = .008 ;%sec
params.RelRefAmp = 4 ; %mV
params.Gleak = 1 ;
params.SampleRate = 10000 ;

MaxDeltaTvector = 10.^[-3:.1:0] ; % sec vector of times shifts for max delta t

signal = normrnd(signal_mean,sqrt(total_var),1,lt) ;     % signal
signal = lowPassFilter(signal,1/si,60) ;    % low pass at 60hz

for a = 1:length(GNoise_frac) ; % for each new change in current noise amplitude

    % set signal and noise std
    GNoise_std = sqrt(total_var * GNoise_frac(a)) ; 
    
    Gsignal = signal*sqrt(total_var/var(signal))*sqrt(1-GNoise_frac(a)) ; % first sqrt corrects for what is lost in low pass and second sqrt adjusts for noise amplitude

    spikeTrain{a} = zeros(NumTrials,lt) ;
    
    
    for Trial=1:NumTrials ; % for each trial
              
        GNoise = normrnd(0,GNoise_std,1,lt) ; % create noise for this trial
        GNoise = lowPassFilter(GNoise,1/si,80) ; 
        
        Gexc = Gsignal + GNoise ;
        Gexc(Gexc<0)=0 ;
        Ginh = zeros(1,lt) ;
        
        V = LIFmodelG_c(Gexc,Ginh,params) ; % run LIF model
        
        spikeTrain{a}(Trial,V==50)= 1 ;
        spiketimes{a}{Trial} = time(spikeTrain{a}(Trial,:)==1) ; % spike times   
        spiketimes_ms{a}{Trial}= time(spikeTrain{a}(Trial,:)==1)*1000 ; %spike times in ms
    end
    
    figure
    rasterPlot(spiketimes{a})
    
    % spike precision measures
    
    [temp,psth] = PsthVar(spikeTrain{a},MaxDeltaTvector,10000) ;
    sumPSTHvar(a,:) = temp ;
    
    [F,temp] = SpikeMatcher(spiketimes{a}, MaxDeltaTvector) ;
    Fmean(a,:) = temp ;
    
    
    [DeltaT, tli_spike, tlj_spike, Percent_Pairs_Quantified, DeltaT_Histogram,tempX,temp] = DeltaT_Distribution_From_SCR(spiketimes_ms{a},.02) ;
    DeltaT_CumProb(a,:) = temp ;
    X_Values(a,:) = tempX ;
    
end

figure
plot(MaxDeltaTvector,sumPSTHvar)
h=gca
set(h,'xscale','log')


figure
plot(MaxDeltaTvector,Fmean)
h=gca
set(h,'xscale','log')

figure
plot(X_Values(a,:),DeltaT_CumProb)
h=gca
set(h,'xscale','log')
