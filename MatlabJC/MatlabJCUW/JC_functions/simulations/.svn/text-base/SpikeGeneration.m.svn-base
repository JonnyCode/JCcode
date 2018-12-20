% SpikeGeneration
%   Construct conductance based model to predict rgc spike response to
%   pattern of excitatory and inhibitory input.  Cell is modeled with
%   capacitance and leak conductance (taken from leak pulses).  Excitatory
%   and inhibitory conductances are taken from vc measurements.  These are
%   used to predict subthreshold voltage.  Potential spike times are
%   identified as points at which predicted voltage exceeds a threshold;
%   these are edited to account for (1) absolute refractory period; (2)
%   relative refractory period (modeled as AHP).  Most of the work is done
%   in PredictRGCVoltage.m and PredictRGCSpikeTimes.m
%
% Known limitations/issues
%   - kludge to compute spike distance in all cases (not doing all possible
%   comparisons)
%   - alter to use cell array of spike times rather than extracting
%
%   Created 7/26/05 FMR


% NOTE: this is setup for 072605c12 only!
[fp, error] = ITCInitializeAnalysis(500000, '/users/fred/data/072605c12');

%--------------------------------------------------------------------------
% Extract excitatory and inhibitory conductance traces
%--------------------------------------------------------------------------

% excitatory conductance data
StartEpoch = 87;        % first vc epoch at -60 mV
EndEpoch = 130;         % last at -60 mV
DrivingForce = 80;      % for conductance calculation

% get some information about epochs
[TargetLength, error] = ITCGetEpochSize(StartEpoch, fp);
[PrePts, error] = ITCGetStmPrePts(StartEpoch, 0, 0, fp);
[SamplingInterval, error] = ITCGetSamplingInterval(StartEpoch, fp);
SamplingInterval = SamplingInterval * 1e-6;

% read in vc data at -60 mV and scale to get excitatory conductance
clear ExcData;
GoodEpochs = 0;
for epoch = StartEpoch:EndEpoch
    [Length, error] = ITCGetEpochSize(epoch, fp);
    if (Length == TargetLength)
        GoodEpochs = GoodEpochs + 1;
        [Data, error] = ITCReadEpoch(epoch, 0, fp);
        ExcData(GoodEpochs, :) = -(Data - mean(Data(1:PrePts))) / DrivingForce;
    end
end

% inhibitory conductance data
StartEpoch = 153;       % first vc epoch at +20 mV
EndEpoch = 177;         % last at +20 mV
[TargetLength, error] = ITCGetEpochSize(StartEpoch, fp);

% read in vc data at +20 mV and scale to get inhibitory conductance
clear InhData;
GoodEpochs = 0;
for epoch = StartEpoch:EndEpoch
    [Length, error] = ITCGetEpochSize(epoch, fp);
    if (Length == TargetLength)
        GoodEpochs = GoodEpochs + 1;
        [Data, error] = ITCReadEpoch(epoch, 0, fp);
        InhData(GoodEpochs, :) = (Data - mean(Data(1:PrePts))) / DrivingForce;
    end
end

% calculate conductances
ExcG = mean(ExcData);
InhG = mean(InhData);
ExcG = ExcG - min(ExcG);
InhG = InhG - min(InhG);

% plot excitatory and inhibitory waveforms
figure(1)
subplot(1, 2, 1) 
plot(ExcG);
xlim([0 30000])
subplot(1, 2, 2) 
plot(InhG);
xlim([0 30000])

% decimate - DecimatePts variable will be used here and below to resample
% data at coarser sampling interval if desired
DecimatePts = 1;
ExcG = decimate(ExcG, DecimatePts);
InhG = decimate(InhG, DecimatePts);

%--------------------------------------------------------------------------
% Extract cell capacitance and leak conductance
%--------------------------------------------------------------------------

% read in leak traces (+5 mV steps at -60 mV typically)
StartLeak = 50;
EndLeak = 59;
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

% get leak conductance from second half of leak step
Leak = Leak - mean(Leak(1:LeakPrePts));
GLeak = mean(Leak(LeakPrePts + LeakStmPts/2:LeakPrePts + LeakStmPts)) * 1e-12 / LeakStep;

% get capacitance by integrating charge in initial transient (c = q/V)
Leak = Leak - mean(Leak(LeakPrePts + LeakStmPts/2:LeakPrePts + LeakStmPts));
Cap = sum(Leak(LeakPrePts:LeakPrePts + LeakStmPts)) * LeakSamplingInterval * 1e-12 / LeakStep;
fprintf(1, 'Leak Conductance = %d S (Resistance = %d ohms) Capacitance = %d F\n', GLeak, 1/GLeak, Cap);

%--------------------------------------------------------------------------
% Predict spike train for single set of parameters
%--------------------------------------------------------------------------

% model parameters
Parameters.VRevInh = -60;           % reversal potential for inhibitory current
Parameters.VRevExc = 0;             % reversal potential for excitatory current
Parameters.VRevLeak = -90;          % reversal potential for leak
Parameters.GLeak = GLeak * 1e9;     % leak conductance
Parameters.Cap = Cap * 1e9;         % capacitance in nF
Parameters.TStep = SamplingInterval * DecimatePts;    % time step for numerical solution of differential equations
Parameters.CurrentNoise = 150;      % SD of Gaussian current noise
Parameters.NoiseTCon = 0.06;        % time constant for noise correlation (in sec)
Parameters.Threshold = -55;         % spike threshold   
Parameters.AbsRefractTime = 0.0025;      % absolute refractory period in sec
Parameters.AHPDecay = 0.02;         % decay time constant (in sec) for AHP
Parameters.AHPAmp = 3;             % amplitude (nS) of AHP conductance
Parameters.AHPVRev = -90;           % reversal potential of AHP conductance

% plot window
Start = 5000 / DecimatePts;
Cutoff = 80000 / DecimatePts;

% predict and plot subthreshold voltage
[Spikes,Voltage] = PredictRGCVoltageAndSpikes(ExcG(1:Cutoff), InhG(1:Cutoff), Parameters);
figure(2); clf;
tme = (1:length(Voltage)) * SamplingInterval * DecimatePts;
plot(tme(Start:Cutoff), Voltage(Start:Cutoff), 'g');
axis tight;
hold on
SpikeRaster = zeros(length(Voltage), 1);
SpikeRaster(round(Spikes/(SamplingInterval * DecimatePts))) = 100;
plot(tme(Start:Cutoff), SpikeRaster(Start:Cutoff), 'r');
axis tight;

%--------------------------------------------------------------------------
% Read on cell spike responses and extra spike times
%--------------------------------------------------------------------------

Threshold = 25;
clear SpikeTime;
StartOnCell = 7;
EndOnCell = 25;

for epoch = 1:EndOnCell-StartOnCell+1
    [Data, error] = ITCReadEpoch(epoch+StartOnCell-1, 0, fp);
    Data = Data - mean(Data);
    PotentialSpikeTimes = find(Data(1:Cutoff*DecimatePts) > Threshold);
    SpikeTime{epoch}(1) = PotentialSpikeTimes(1) * SamplingInterval;
    NumSpikes = 1;
    for t = 1:length(PotentialSpikeTimes)
        TimeDiff = PotentialSpikeTimes(t) * SamplingInterval - SpikeTime{epoch}(NumSpikes);
        if (TimeDiff > Parameters.AbsRefractTime)
            NumSpikes = NumSpikes + 1;
            SpikeTime{epoch}(NumSpikes) = PotentialSpikeTimes(t) * SamplingInterval;
        end
    end
end

% plot single example
tme = (1:length(Voltage)) * SamplingInterval * DecimatePts;
SpikeRaster = zeros(length(Voltage), 1);
SpikeRaster(round(SpikeTime{1}/ (SamplingInterval * DecimatePts))) = 100;
plot(tme(Start:Cutoff), SpikeRaster(Start:Cutoff)-110);
hold off;

%--------------------------------------------------------------------------
% View different parts of record
%--------------------------------------------------------------------------
Start = 10000/DecimatePts;
Cutoff = 45000/DecimatePts;
figure(2); clf;
tme = (1:length(Voltage)) * SamplingInterval * DecimatePts;
SpikeRaster = zeros(length(Voltage), 1);
SpikeRaster(round(SpikeTime{1}/(SamplingInterval * DecimatePts))) = 100;
plot(SpikeRaster(Start:Cutoff)-110);
hold on;
SpikeRaster = zeros(length(Voltage), 1);
SpikeRaster(round(SpikeTime{EndOnCell-StartOnCell}/(SamplingInterval * DecimatePts))) = 100;
plot(SpikeRaster(Start:Cutoff)-220);
SpikeRaster = zeros(length(Voltage), 1);
SpikeRaster(round(Spikes / (SamplingInterval * DecimatePts))) = 100;
plot(SpikeRaster(Start:Cutoff), 'r');
plot(Voltage(Start:Cutoff), 'g');
hold off;
axis tight;

%--------------------------------------------------------------------------
% Search for best parameters using Victor distance as metric
%--------------------------------------------------------------------------

% model parameters
Parameters.VRevInh = -60;           % reversal potential for inhibitory current
Parameters.VRevExc = 0;             % reversal potential for excitatory current
Parameters.VRevLeak = -90;          % reversal potential for leak
Parameters.GLeak = GLeak * 1e9;     % leak conductance
Parameters.Cap = Cap * 1e9;         % capacitance in nF
Parameters.TStep = SamplingInterval * DecimatePts;    % time step for numerical solution of differential equations
Parameters.CurrentNoise = 150;      % SD of Gaussian current noise
Parameters.NoiseTCon = 0.06;        % time constant for noise correlation (in sec)
Parameters.Threshold = -55;         % spike threshold   
Parameters.AbsRefractTime = 0.0025;      % absolute refractory period in sec
Parameters.AHPDecay = 0.02;         % decay time constant (in sec) for AHP
Parameters.AHPAmp = 3;             % amplitude (nS) of AHP conductance
Parameters.AHPVRev = -90;           % reversal potential of AHP conductance

TestTimeShift = 0.004;              % time shift for spike distance measurement in sec
TestCost = 2/TestTimeShift;         % cost corresponding to time shift above
BestRatio = 1;                      % initial value for best ratio of spike distance 
                                    % to total spike count (i.e. sum of counts in two spike trains)

% loop through threshold and refractory parameters to find best match
for thresh = 1:5
    Parameters.Threshold = -60 + 2*thresh;
    fprintf(1, 'Threshold = %d\n', Parameters.Threshold);
    for rtime = 1:5
        Parameters.AbsRefractTime = rtime * 0.001;
        fprintf(1, '\n');
        for decay = 1:5
            Parameters.AHPDecay = 0.01 * decay;
            for amp = 1:5
                Parameters.AHPAmp =  amp;

                [Spikes,Voltage] = PredictRGCVoltageAndSpikes(ExcG(1:Cutoff), InhG(1:Cutoff), Parameters);
                d = 0;
                SpikeCount = 0;
                for epoch = 1:EndOnCell - StartOnCell
                    d = d + spkd_c(Spikes, SpikeTime{epoch}, length(Spikes), length(SpikeTime{epoch}), TestCost);
                    SpikeCount = SpikeCount + length(Spikes) + length(SpikeTime{epoch});
                end
                d = d / (EndOnCell - StartOnCell);
                SpikeCount = SpikeCount / (EndOnCell - StartOnCell);
                fprintf(1, '\tAbs = %d, Rel = %d %d , d = %d (%d) ratio = %d (%d)\n', Parameters.AbsRefractTime, Parameters.AHPDecay, Parameters.AHPAmp, d, SpikeCount, d/SpikeCount, BestRatio);
                if (d/SpikeCount < BestRatio)
                    BestRatio = d/SpikeCount;
                end
            end
        end
    end
end

% average distance between measured trains (does not do all
% comparisons, just all other trains with last)
SpikeCount = 0;
d = 0;
for epoch = 1:EndOnCell - StartOnCell - 1
    d = d + spkd_c(SpikeTime{EndOnCell - StartOnCell}, SpikeTime{epoch}, length(SpikeTime{EndOnCell - StartOnCell}), length(SpikeTime{epoch}), TestCost);
    SpikeCount = SpikeCount + length(SpikeTime{EndOnCell-StartOnCell}) + length(SpikeTime{epoch});
end
d = d / (EndOnCell - StartOnCell - 1);
SpikeCount = SpikeCount / (EndOnCell - StartOnCell - 1);
fprintf(1, 'distance for measured trains %d (%d) ratio = %d\n', d, SpikeCount, d/SpikeCount);

%--------------------------------------------------------------------------
% Victor distance for predicted trains
%--------------------------------------------------------------------------

Parameters.VRevInh = -60;           % reversal potential for inhibitory current
Parameters.VRevExc = 0;             % reversal potential for excitatory current
Parameters.VRevLeak = -90;          % reversal potential for leak
Parameters.GLeak = GLeak * 1e9;     % leak conductance
Parameters.Cap = Cap * 1e9;         % capacitance in nF
Parameters.TStep = SamplingInterval * DecimatePts;    % time step for numerical solution of differential equations
Parameters.CurrentNoise = 150;      % SD of Gaussian current noise
Parameters.NoiseTCon = 0.06;        % time constant for noise correlation (in sec)
Parameters.Threshold = -55;         % spike threshold   
Parameters.AbsRefractTime = 0.0025;      % absolute refractory period in sec
Parameters.AHPDecay = 0.02;         % decay time constant (in sec) for AHP
Parameters.AHPAmp = 3;             % amplitude (nS) of AHP conductance
Parameters.AHPVRev = -90;           % reversal potential of AHP conductance

NumTrains = 20;

clear Spikes;
for train = 1:NumTrains
    fprintf(1, 'starting train %d\n', train);
    [Spikes{train},Voltage] = PredictRGCVoltageAndSpikes(ExcG(1:Cutoff), InhG(1:Cutoff), Parameters);
    figure(1);
    plot(Voltage);
    pause(1);
end
d = 0;
SpikeCount =0;
for epoch = 2:NumTrains
    d = d + spkd_c(Spikes{1}, Spikes{epoch}, length(Spikes{1}), length(Spikes{epoch}), TestCost);
    SpikeCount = SpikeCount + length(Spikes{1}) + length(Spikes{epoch});
end
d = d / (NumTrains - 1);
SpikeCount = SpikeCount / (NumTrains-1);
fprintf(1, 'distance for predicted trains %d (%d) ratio = %d\n', d, SpikeCount, d/SpikeCount);


%--------------------------------------------------------------------------
% Play with parameters
%--------------------------------------------------------------------------

Parameters.VRevInh = -60;           % reversal potential for inhibitory current
Parameters.VRevExc = 0;             % reversal potential for excitatory current
Parameters.VRevLeak = -60;          % reversal potential for leak
Parameters.GLeak = GLeak * 1e9;     % leak conductance
Parameters.Cap = Cap * 1e9;         % capacitance in nF
Parameters.TStep = SamplingInterval;    % time step for numerical solution of differential equations
Parameters.CurrentNoise = 100;      % SD of Gaussian current noise
Parameters.NoiseTCon = 0.06;        % time constant for noise correlation (in sec)
Parameters.Threshold = -55;         % spike threshold   
Parameters.AbsRefractPts = 20;      % absolute refractory period
Parameters.AHPDecay = 0.01;         % decay time constant (in sec) for AHP
Parameters.AHPAmp = 20;             % amplitude (nS) of AHP conductance
Parameters.AHPVRev = -90;           % reversal potential of AHP conductance

ExcAmp = 1;
InhAmp = 4;

NumTrains = 5;

clear Spikes;
for train = 1:NumTrains
    fprintf(1, 'starting train %d\n', train);
    [Spikes{train},Voltage] = PredictRGCVoltageAndSpikes(ExcG(1:Cutoff), InhG(1:Cutoff), Parameters);
    figure(1);
    plot(Voltage);
    pause(1);
end
d = 0;
SpikeCount =0;
for epoch = 2:NumTrains
    d = d + spkd_c(Spikes{1}, Spikes{epoch}, length(Spikes{1}), length(Spikes{epoch}), TestCost);
    SpikeCount = SpikeCount + length(Spikes{1}) + length(Spikes{epoch});
end
d = d / (NumTrains - 1);
SpikeCount = SpikeCount / (NumTrains-1);
fprintf(1, 'distance for predicted trains %d (%d) ratio = %d\n', d, SpikeCount, d/SpikeCount);

%--------------------------------------------------------------------------
% Simple test cases - excitatory and inhibitory steps
%--------------------------------------------------------------------------

Parameters.VRevInh = -60;           % reversal potential for inhibitory current
Parameters.VRevExc = 0;             % reversal potential for excitatory current
Parameters.VRevLeak = -60;          % reversal potential for leak
Parameters.GLeak = GLeak * 1e9;     % leak conductance
Parameters.Cap = Cap * 1e9;         % capacitance in nF
Parameters.TStep = SamplingInterval;    % time step for numerical solution of differential equations
Parameters.CurrentNoise = 0;      % SD of Gaussian current noise
Parameters.NoiseTCon = 0.06;        % time constant for noise correlation (in sec)
Parameters.Threshold = -55;         % spike threshold   
Parameters.AbsRefractPts = 20;      % absolute refractory period
Parameters.AHPDecay = 0.05;         % decay time constant (in sec) for AHP
Parameters.AHPAmp = 1;             % amplitude (nS) of AHP conductance
Parameters.AHPVRev = -90;           % reversal potential of AHP conductance

clear ExcG InhG
Cutoff = 3000;
ExcG(1:1000) = 0;
ExcG(1000:2000) = 2;
ExcG(2000:Cutoff) = 0;

InhG(1:1000) = 0;
InhG(1000:2000) = 2;
InhG(2000:Cutoff) = 0;

clear Voltage Spikes
figure(1);clf
[Spikes,Voltage] = PredictRGCVoltageAndSpikes(ExcG(1:Cutoff), InhG(1:Cutoff), Parameters);
plot(Voltage);
hold on
SpikeRaster = zeros(length(Voltage), 1);
SpikeRaster(Spikes) = 100;
plot(SpikeRaster(Start:Cutoff), 'r');
axis tight;
hold off;
