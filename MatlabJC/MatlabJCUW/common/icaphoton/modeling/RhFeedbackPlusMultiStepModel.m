% Calculate series of dim flash responses for feedback model for 
% rhodopsin shutoff.  First generate time course of 
% rhodopsin activity assuming rhodopsin shutoff is controlled
% by a feedback signal that accumulates linearly over time
% and acts on rhodopsin shutoff cooperatively (with specified
% cooperativity).  Convolve with linear approximation to transduction
% cascade to generate modeled response.
%
% Created 7/01 FMR
%
% ReturnedCondition = RhFeedbackPlusMultiStepModel(Cooperativity, InitialShutoffRate, FeedbackGain, Condition, Filter, NumResponses)

function ReturnedCondition = RhFeedbackPlusMultiStepModel(NumSteps, Cooperativity, InitialShutoffRate, FeedbackGain, Condition, Filter, NumResponses)

filt = fft(Filter);
ReturnedCondition = Condition;

for resp = 1:NumResponses

	RhTimeCourse(1:Condition.EpochPts) = 0;
	CurrentStep = 1;
	CurrentCatalyticActivity = 1;
	ShutoffRate = InitialShutoffRate;

	for cnt=1:Condition.EpochPts
		RhTimeCourse(cnt) = 1;
		ShutoffRate = InitialShutoffRate * (FeedbackGain * cnt)^Cooperativity;
		if (rand(1) < ShutoffRate)
			break;
		end
	end	

	% for each time point decide whether shutoff reaction has occurred.
	% if it has, update catalytic activity and shutoff rate.
	for cnt=1:Condition.EpochPts
		RhTimeCourse(cnt) = CurrentCatalyticActivity;
		% generate random number between 0 and 1 (uniform dist) and compare to shutoff rate
		if (rand(1) < ShutoffRate)
			ShutoffRate = InitialShutoffRate * (FeedbackGain * cnt)^Cooperativity / CurrentStep;
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
