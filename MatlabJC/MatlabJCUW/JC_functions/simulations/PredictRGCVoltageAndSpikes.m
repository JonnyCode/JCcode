function [Spikes,Voltage] = PredictRGCVoltageAndSpikes(ExcG, InhG, Params)

% integrate and fire model for spiking.  Subthreshold dyamics obey
% -CdV/dt = INoise + GExc(V - Vexc) - GInh(V - VInh) + GLeak(V - VLeak) +
%                   GAHP(V-VAHP) 
% generate spike everytime you cross specified threshold.
%
%   2/07 FMR

% generate noise current
if (Params.CurrentNoise > 0)
    CurrentNoise = normrnd(0, 1, length(ExcG), 1);
    Filter = ((1:length(ExcG)) * Params.TStep)';
    Filter = exp(-Filter / Params.NoiseTCon);
    FilterFFT = fft(Filter);
    CurrentNoiseFFT = fft(CurrentNoise);
    CurrentNoise = real(ifft(FilterFFT .* CurrentNoiseFFT));
    CurrentNoise = CurrentNoise - mean(CurrentNoise);
    CurrentNoise = CurrentNoise * Params.CurrentNoise / std(CurrentNoise);
else
    CurrentNoise(1:length(ExcG)) = 0;
end

% jam parameters into vector for call to c routine
Parameters(1) = Params.VRevInh;
Parameters(2) = Params.VRevExc;
Parameters(3) = Params.VRevLeak;
Parameters(4) = Params.GLeak;
Parameters(5) = Params.Cap;
Parameters(6) = Params.TStep;
Parameters(7) = Params.Threshold;
Parameters(8) = Params.AbsRefractTime;
Parameters(9) = Params.AHPDecay;
Parameters(10) = Params.AHPAmp;
Parameters(11) = Params.AHPVRev;

% call c routine to do dirty work - i.e. to solve above differential
% equation for subthreshold voltages, generate spikes when threshold
% exceeded and update AHP conductance following each spike
[Spikes,Voltage] = PredictRGCVoltageAndSpikes_c(ExcG, InhG, CurrentNoise, Parameters, length(ExcG));
