% correlated noise simulation 4/12/09
ts = .0001 ; %sec, time step
lParab = .01 ; %sec, length of parabola
maxt = 2 ; %sec, length of time
NoiseCom_Frac = [0:1] ; % range of common noise fraction of total
NoiseTot_Var = 2 ; % var of total noise in each conductance
NumTrials = 40 ; % number of repeats 
samplerate = 10000 ;
timeShift = .000 ; %sec inh is delayed
color = ['b*','r*','g*','k*','c*','b*','r*','g*','k*','c*','b*','r*','g*','k*'] ;

params.Einh = -80 ;
params.Eexc = 0 ;
params.Eleak = -60 ;
params.cap = .06 ;%nF
params.Vrest = -60;
params.Vthresh = -45 ;
params.AbsRef  = .002 ; %sec
params.RelRefTau = .008 ;%sec
params.RelRefAmp = 4 ; %mV


time = [ts:ts:maxt] ; %time in sec
lt = length(time) ;

S = zeros(1,lt) ; % signal

b=0 ;
for a= 1:lParab*2/ts:lt ; % for each parab that can fit
    b=b+.1 ; 
    S(a) = b ; % make a psuedo delta function that increases by .1
end

parab = (-([-lParab:ts:lParab].^2) +lParab^2)./lParab^2 ; % make a parabola

S = conv(S,parab) ;
S = S(1:lt) ;

Trial = 0 ;
spikes = zeros(NumTrials*length(NoiseCom_Frac),lt) ;

for a = 1:length(NoiseCom_Frac) ; % for each new change in common noise amplitude
    NoiseCom_Std = sqrt(NoiseTot_Var*NoiseCom_Frac(a)) ; % standard deviation of common noise
    NoiseInd_Std = sqrt(NoiseTot_Var - NoiseCom_Std^2) ; % standard deviation of independant noise
    
    for b=1:NumTrials ; % for each trial
        Trial = Trial + 1 ;
        
        NoiseCommon = normrnd(0,NoiseCom_Std,1,lt) ;
        NoiseExc = normrnd(0,NoiseInd_Std,1,lt) ;
        NoiseInh = normrnd(0,NoiseInd_Std,1,lt) ;
        
        GexcPre = S + NoiseCommon + NoiseExc ;
        GinhPre = S + NoiseCommon + NoiseInh ;
        
        pntShift = timeShift*samplerate ;
        GinhPre = [zeros(1,pntShift),GinhPre(1:end-pntShift)] ; %time shift inhibiton
        
        GexcPre = lowPassFilter(GexcPre, 10000, 60) ; % low pass filter both
        GinhPre = lowPassFilter(GinhPre, 10000, 60) ;
        
        GexcPre(GexcPre<0) = 0 ; %rectify
        GinhPre(GinhPre<0) = 0 ;
        
        Gexc(Trial,:) = GexcPre ;
        Ginh(Trial,:) = GinhPre ;
        Gleak = 0 ;
        
        V_trace(Trial,:) = LIFmodelG(Gexc(Trial,:),Ginh(Trial,:),Gleak,samplerate,params) ;
        spikes(Trial,V_trace(Trial,:)==50) = Trial ;
        spiketimes{Trial} = time(V_trace(Trial,:)==50) ; 
    end
    
end

figure
rasterPlot(spiketimes)
%plot(time,spikes,'b*') 

figure
plot(time,S,'k')
hold on
plot(time,Gexc(1,:),'b')
plot(time,Ginh(1,:),'r')

figure
plot(time,S,'k')
hold on
plot(time,Gexc(end,:),'b')
plot(time,Ginh(end,:),'r')

% check spike precision with spike distance metric 
for a = 1:length(NoiseCom_Frac) ; % for each new change in common noise amplitude
    for b = a*NumTrials-NumTrials+1:a*NumTrials ; % for each trial
        for c = a*NumTrials-NumTrials+1:a*NumTrials ; % for each possible pair
            sd{a}(b,c)=spkd_c(spiketimes{b},spiketimes{c},length(spiketimes{b}), length(spiketimes{c}), .02); % runs spike distance function at a given cost
        end
    end
    sd_Mean(a) = mean(sd{a}(sd{a}~=0)) ; % get the mean spike distance for a particular amplitude of common noise
end

figure
plot(NoiseCom_Frac,sd_Mean)
xlabel('fraction of noise which is common')
ylabel('mean spike distance')








