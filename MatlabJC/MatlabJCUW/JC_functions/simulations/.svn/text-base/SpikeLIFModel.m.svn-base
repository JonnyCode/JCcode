% This script is a modified version of Fred's "SpikeGeneration" script and
% will help model spike generation from measured conductances

% JC 5/16/08

% recordings to be used
Iexc_Epochs = {'[283]','[285]','[287]','[289]'} ; % All isolate exc epochs you want grouped
Iexc_Cells = {'050808Bc1','050808Bc1','050808Bc1','050808Bc1'} ; % The cells from which they originated

Iinh_Epochs = {'[395]','[397]','[399]','[401]'} ; %  
Iinh_Cells = {'050808Bc1','050808Bc1','050808Bc1','050808Bc1'} ; %

Leak_Epochs = {} ; % The current reponse to 5mV step to estimate Gleak
Leak_Cells = {} ;

Spike_Epochs = {} ; % Spike output you want to compare against the modeled spikes
Spike_Cells = {} ; %

prepnts = 5000 ;
stmpnts = 5000 ;
tailpnts = 5000 ;

% Model parameters
Eexc = 0 ; %mV reversal potential for exc
Einh = -60 ; 
Eleak = -90 ;
Eahp = -90 ;

GLeak = 40 *1e-9 ; %S leak conductance (enter 0 if you want to extract the Gleak from the leak traces) 
Capacitance = .08 ; %nF (enter 0 if you want to extract the cells capacitance from the leak traces) 
% Ginh_offset = 2 ; %nS this will be added to Ginh to prevent high firing rates (enter NaN if you want this optimized see eq 5.12 "Theoretical Neuroscience" )

SpikeThreshold = -38 ; %mV
AbsRefTime = .0025 ; %Sec
AHPDecayTau = .02 ; %Sec decay time constant of G activated by spike 
AHPAmp = 5 ; %nS amplitude of AHP G

DecimatePts = 1 ;    % time step for numerical solution of differential equations = sample interval*decimatePnts

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
Gexc{set} = Gexc{set} - mean(Gexc{set}(1:prepnts)) +25 ;
% Gexc{set} = Gexc{set} - min(Gexc{set}) ;

Ginh{set} = Iinh_mean{set}/(Eexc-Einh) ;
%Ginh{set} = Ginh{set} - mean(Ginh{set}(1:prepnts)) ;
Ginh{set} = Ginh{set} - min(Ginh{set}) ;

Gexc{set} = decimate(Gexc{set}, DecimatePts)  ;
Ginh{set} = decimate(Ginh{set}, DecimatePts)  ;

if GLeak == 0 | Capacitance == 0 ;
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

if Capacitance == 0 ;
% get capacitance by integrating charge in initial transient (c = q/V)
Leak = Leak - mean(Leak(LeakPrePts + LeakStmPts/2:LeakPrePts + LeakStmPts)) ;
Cap = sum(Leak(LeakPrePts:LeakPrePts + LeakStmPts)) * LeakSamplingInterval * 1e-12 / LeakStep ;
end

end % end leak and capacitance loop
end % end set loop

for set=1:length(Gexc) ;
% model parameters
Parameters.VRevInh = Einh ;                             % reversal potential for inhibitory current
Parameters.VRevExc = Eexc ;                             % reversal potential for excitatory current
Parameters.VRevLeak = Eleak ;                           % reversal potential for leak
Parameters.GLeak = GLeak ;                              % leak conductance
Parameters.Cap = Capacitance ;                          % capacitance in nF
Parameters.TStep = SamplingInterval * DecimatePts ;     % time step for numerical solution of differential equations
Parameters.CurrentNoise = 0 ;                           % SD of Gaussian current noise
Parameters.NoiseTCon = 0.06 ;                           % time constant for noise correlation (in sec)
Parameters.Threshold = SpikeThreshold ;                 % spike threshold   
Parameters.AbsRefractTime = AbsRefTime ;                    % absolute refractory period in sec
Parameters.AHPDecay = AHPDecayTau ;                            % decay time constant (in sec) for AHP
Parameters.AHPAmp = AHPAmp ;                                 % amplitude (nS) of AHP conductance
Parameters.AHPVRev = Eahp ;                              % reversal potential of AHP conductance

[Spikes{set},Voltage{set}] = PredictRGCVoltageAndSpikes(Gexc{set}, Ginh{set}, Parameters);

end % second set loop

for a= 1:length(Gexc) ;
figure
tme = (1:length(Voltage{a})) * SamplingInterval * DecimatePts;
SpikeRaster{a} = ones(length(Voltage{a}), 1)*SpikeThreshold;
SpikeRaster{a}(round(Spikes{a}/(SamplingInterval * DecimatePts))) = 10;
subplot(2,1,1)
plot(tme, SpikeRaster{a}, 'r');
hold on
plot(tme, Voltage{a}, 'g');
axis tight;
subplot(2,1,2)
plot(Gexc{a})
hold on
plot(Ginh{a},'r')
end