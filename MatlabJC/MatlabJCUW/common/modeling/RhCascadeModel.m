% Calculate series of dim flash responses for multi-step shutoff
% model for rhodopsin activity.  First generate time course of 
% rhodopsin activity assuming rhodopsin shutoff is described as a
% series of independent first order reactions.  Each shutoff 
% reaction is assumed to produce an equal decrease in rhodopsin's 
% catalytic activity.  The rate constant for each shutoff reaction
% increases as well, such that each reaction controls an equal amount
% of rhodopsin's activity
%
% Created 7/01 FMR

function ReturnedCondition = RhFeedbackModel(NumSteps, InitialShutoffRate, Condition, Filter, NumResponses)

% put specified filter into frequency domain
filt = fft(Filter);
ReturnedCondition = Condition;

% generate series of responses
for resp = 1:NumResponses

	% initial settings for rhodopsin activity
	RhTimeCourse(1:Condition.EpochPts) = 0;
	CurrentStep = 1;
	CurrentCatalyticActivity = 1;
	ShutoffRate = InitialShutoffRate;
	
	% for each time point decide whether shutoff reaction has occurred.
	% if it has, update catalytic activity and shutoff rate.
	for cnt=1:Condition.EpochPts
		RhTimeCourse(cnt) = CurrentCatalyticActivity;
		% generate random number between 0 and 1 (uniform dist) and compare to shutoff rate
		if (rand(1) < ShutoffRate)
			ShutoffRate = InitialShutoffRate / CurrentStep;
			CurrentCatalyticActivity = 1 - CurrentStep / NumSteps;
			CurrentStep = CurrentStep + 1;
		end
		% are we done yet?
		if (CurrentStep == NumSteps)
			break;
		end
	end	

	% convolve rhodopsin time course with linear filter to generate modeled response
	RhTimeCourse = fft(RhTimeCourse);
	ReturnedCondition.EpochData(resp, :) = real(ifft(RhTimeCourse .* filt));
	
end
