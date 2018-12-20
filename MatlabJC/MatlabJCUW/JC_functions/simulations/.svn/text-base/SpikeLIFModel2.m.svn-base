% This script will help model spike generation from measured conductances
% via a leaky integrate and fire model

% JC 5/16/08

% recordings to be used
Iexc_Epochs = {'[289]','[287]','[285]','[283]','[282]','[284]','[286]','[288]'} ; % All isolate exc epochs you want grouped
Iexc_Cells = {'050808Bc1','050808Bc1','050808Bc1','050808Bc1','050808Bc1','050808Bc1','050808Bc1','050808Bc1'} ; % The cells from which they originated

Iinh_Epochs = {'[401]','[399]','[397]','[395]','[394]','[396]','[398]','[400]'} ; %  
Iinh_Cells = {'050808Bc1','050808Bc1','050808Bc1','050808Bc1','050808Bc1','050808Bc1','050808Bc1','050808Bc1'} ; %

Leak_Epochs = {} ; % The current reponse to 5mV step to estimate Gleak
Leak_Cells = {} ;

Spike_Epochs = {} ; % Spike output you want to compare against the modeled spikes
Spike_Cells = {} ; %

prepnts = 5000 ;
stmpnts = 5000 ;
tailpnts = 5000 ;

% Model parameters
Eexc = 0 ; %mV reversal potential for exc
Einh = -80 ; 
Eleak = -90 ;
Eahp = -90 ;

Vrest = -50 ; %mV resting potential
Vthresh = -45 ; %mV spike threshold

GLeak = 20  ; %nS leak conductance (enter 0 if you want to extract the Gleak from the leak traces) 
Cap = .08 ; %nF (enter 0 if you want to extract the cells capacitance from the leak traces) 
% Ginh_offset = 2 ; %nS this will be added to Ginh to prevent high firing rates (enter NaN if you want this optimized see eq 5.12 "Theoretical Neuroscience" )

abs_ref = .0025 ; %Sec
AHPDecayTau = .02 ; %Sec decay time constant of G activated by spike 
AHPAmp = 5 ; %nS amplitude of AHP G

timestep = .00001 ; %sec time step used for differentiation

% LIF Model
for set = 1:length(Iexc_Epochs) ; % for each set of conductances you want to run
    
[fpExc, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',(Iexc_Cells{set})]);
[fpInh, error] = ITCInitializeAnalysis(500000, ['~/Data/Primate/',(Iinh_Cells{set})]);

[SamplingInterval, error] = ITCGetSamplingInterval(288, fpExc); % assume the same sampling interval for all traces
SamplingInterval = SamplingInterval * 1e-6;

Iexc_epochNums{set} = str2num(Iexc_Epochs{set}) ;
Iinh_epochNums{set} = str2num(Iinh_Epochs{set}) ;

for a = 1:length(Iexc_epochNums{set}) ; % get each individ exc recording...
    [Iexc{set}(a,:), error] = ITCReadEpoch(Iexc_epochNums{set}(a), 0, fpExc);
end

for b = 1:length(Iinh_epochNums{set}) ;
    [Iinh{set}(b,:), error] = ITCReadEpoch(Iinh_epochNums{set}(b), 0, fpInh);
end

Iexc_mean{set} = mean(Iexc{set},1) ; % take the mean of each individual recording
Iinh_mean{set} = mean(Iinh{set},1) ;

Gexc{set} = Iexc_mean{set}/(Einh-Eexc) ; % find conductances
Gexc{set} = Gexc{set} - mean(Gexc{set}(1:prepnts)) + 20  ;
% Gexc{set} = repmat(mean(Gexc{set}(1:prepnts)),1,length(Gexc{set})) ; % remove exc keep only pre step mean

% Gexc{set} = Gexc{set} - min(Gexc{set}) ;

Ginh{set} = Iinh_mean{set}/(Eexc-Einh) ;
%Ginh{set} = Ginh{set} - mean(Ginh{set}(1:prepnts)) ;
Ginh{set} = Ginh{set} - min(Ginh{set}) ;
Ginh{set} = repmat(mean(Ginh{set}(1:prepnts)),1,length(Ginh{set})) ; % remove inh keep only pre step mean

% adjust time parameters for samplerate
h = timestep/SamplingInterval ; %pnts    % adjust time step for samplerate
time = [1:h:length(Gexc{set})] ;%pnts    % create a time vector to interpolate for a smaller time step than what was sampled

% resample to make timestep
Gexc{set} = interp1([1:length(Gexc{set})],Gexc{set},time)  ;
Ginh{set} = interp1([1:length(Ginh{set})],Ginh{set},time)  ;

if GLeak == 0 | Cap == 0 ;
%--------------------------------------------------------------------------
% Extract cell leak conductance
%--------------------------------------------------------------------------

% read in leak traces (+5 mV steps at -60 mV typically)
StartLeak = 230;
EndLeak = 239;
LeakStep = 5e-3;            % in mV
clear LeakData;
for epoch = StartLeak:EndLeak
    [Data, error] = ITCReadEpoch(epoch, 0, fp);
    LeakData(epoch - StartLeak + 1, :) = Data;
end
Leak = mean(LeakData);
[LeakPrePts, error] = ITCGetStmPrePts(StartLeak, 0, 0, fp);
[LeakStmPts, error] = ITCGetStmPts(StartLeak, 0, 0, fp);
[LeakSamplingInterval, error] = ITCGetSamplingInterval(StartLeak, fp);
LeakSamplingInterval = LeakSamplingInterval * 1e-6;

Leak = Leak - mean(Leak(1:LeakPrePts));

if GLeak == 0 ;
% get leak conductance from second half of leak step
GLeak = mean(Leak(LeakPrePts + LeakStmPts/2:LeakPrePts + LeakStmPts)) * 1e-12 / LeakStep ;
GLeak = GLeak * 1e9 ;
end 

if Cap == 0 ;
% get capacitance by integrating charge in initial transient (c = q/V)
Leak = Leak - mean(Leak(LeakPrePts + LeakStmPts/2:LeakPrePts + LeakStmPts)) ;
Cap = sum(Leak(LeakPrePts:LeakPrePts + LeakStmPts)) * LeakSamplingInterval * 1e-12 / LeakStep ;
end

end % end leak and capacitance loop
end % end set loop

% MODEL
for set=1:length(Gexc) ;

Isyn{set} = time ;
Isyn{set}(1) = 0 ;

V{set} = time ;
V{set}(1) = Vrest ;

ref = 0 ;


for t = 2:length(time) ;

% calculate synaptic current
Isyn{set}(t) = Ginh{set}(t)*(V{set}(t-1)-Einh) + Gexc{set}(t)*(V{set}(t-1)-Eexc) + GLeak*(V{set}(t-1)-Eleak) ;%pA

% implement integrate and fire model c*dV/dt=Sum(g(V-E))
% dV/dt = -Ir/C , Ir = Sum(g(v-E))
if ~ref                    % if not within refrectory period... 
      V{set}(t) = V{set}(t-1) + timestep*(-Isyn{set}(t)/Cap) ;    %calculate the voltage

else
    V{set}(t) = Vrest ;       % reset voltage to rest
    ref = ref-1 ;
end

% Spike
if (V{set}(t) > Vthresh)      % if voltage is higher than threshold...
    V{set}(t) = 50 ;        % spike and make the next point at rest
    ref = round(abs_ref/timestep) ;
end

end % this ends t loop 

% Note times of spikes
Spiketimes{set}= find(V{set} == 50) ;

NumPreSpikes(set) = length(find(Spiketimes{set}<prepnts/h)) ; 
NumStimSpikes(set) = length(find(Spiketimes{set}>prepnts/h & Spiketimes{set}<(prepnts+stmpnts)/h)) ;
NumPostSpikes(set) = length(find(Spiketimes{set}>(prepnts+stmpnts)/h)) ;

NumSpikeDiff(set) = NumStimSpikes(set) - NumPreSpikes(set) ;

end % second set loop


for set = 1:length(Gexc) ;
figure(set)


subplot(3,1,1);
hold on
plot(V{set})
title('voltage')
xlabel('points')
ylabel('mV')

subplot(3,1,2)
plot(Gexc{set}, 'b')
hold on
plot(Ginh{set},'r')
title('synaptic conductances')
xlabel('points')
ylabel('nS')

subplot(3,1,3)
plot(Isyn{set})
title('I synaptic')
xlabel('points')
ylabel('pA')

end % figure set loop

figure
plot(NumSpikeDiff,'r*')
