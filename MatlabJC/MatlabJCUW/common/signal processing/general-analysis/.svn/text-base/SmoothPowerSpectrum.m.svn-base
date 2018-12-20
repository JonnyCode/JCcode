function NewPowerSpec = SmoothPowerSpectrum(OriginalPowerSpec, OriginalFreq, SmoothFact, SkipPts)
%	NewPowerSpec = SmoothPowerSpectrum(OriginalPowerSpec, OriginalFreq, SmoothFact, SkipPts)
%
% This function takes a power spectrum and smooths it such that the points
% are roughly equally spaced on a log frequency scale.  Routine computes
% average and sd of points in window of width given by pnt^SmoothFact,
% where pnt is the pnt number in the original power spectrum.  Returns
% structure with three fields: average power spectrum values, sd, and
% frequency of new samples.  

% Created: FMR 10/03
% Revised 4/06 FMR
%   added skippts

NumOriginalPoints = length(OriginalPowerSpec) - SkipPts;
NumFinalPoints = (NumOriginalPoints/2).^(1/SmoothFact) + SkipPts;
StartWindow = SkipPts;

for cnt = 1:SkipPts
	NewPowerSpec.PowerSpec(cnt) = OriginalPowerSpec(cnt);
	NewPowerSpec.SDSpec(cnt) = nan;
	NewPowerSpec.Freq(cnt) = OriginalFreq(cnt);
end

for cnt = SkipPts+1:NumFinalPoints
	EndWindow = round(cnt^SmoothFact);
	if (EndWindow <= StartWindow)
		EndWindow = StartWindow + 1;
	end
	NewPowerSpec.PowerSpec(cnt) = mean(OriginalPowerSpec(StartWindow:EndWindow));
	NewPowerSpec.SDSpec(cnt) = sqrt(var(OriginalPowerSpec(StartWindow:EndWindow)));
	NewPowerSpec.Freq(cnt) = mean(OriginalFreq(StartWindow:EndWindow+1));
	StartWindow = EndWindow+1;
end
