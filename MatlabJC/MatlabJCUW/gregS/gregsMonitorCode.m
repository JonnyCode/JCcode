rameMonitorData = LowPassFilter(frameMonitorData,20,1E-4);
%figure(3);
%plot(frameMonitorData);

stimIntervals = getThresCross(frameMonitorData,0.5,1);
periods = diff(stimIntervals);
%pause;

%first one seems off (I wonder if this is affecting my stim
meanPeriod = mean(periods(2:end));
pointsPerStimPt = meanPeriod*SampleInterval;

cycleFrames = SampleEpoch.stimuli.get(StimStreamName).get('parameters').get('CycleFrames');
F1 = 60/cycleFrames; %Hz
F2 = 2*F1; %Hz

function Ind = getThresCross(V,th,dir)
%dir 1 = up, -1 = down

Vorig = V(1:end-1);
Vshift = V(2:end);

if dir>0
    Ind = find(Vorig<th & Vshift>=th) + 1;
else
    Ind = find(Vorig>=th & Vshift<th) + 1;
end
