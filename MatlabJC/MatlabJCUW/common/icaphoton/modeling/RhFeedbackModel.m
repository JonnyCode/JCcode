% Calculate series of dim flash responses for feedback model for 
% rhodopsin shutoff.  First generate time course of 
% rhodopsin activity assuming rhodopsin shutoff is controlled
% by a feedback signal that accumulates linearly over time
% and acts on rhodopsin shutoff cooperatively (with specified
% cooperativity).  Convolve with linear approximation to transduction
% cascade to generate modeled response.
%
% Created 7/01 FMR

function ReturnedCondition = RhFeedbackModel(Cooperativity, InitialShutoffRate, FeedbackGain, Condition, Filter, NumResponses)

filt = fft(Filter);
ReturnedCondition = Condition;

for resp = 1:NumResponses

	RhTimeCourse(1:Condition.EpochPts) = 0;

	for cnt=1:Condition.EpochPts
		RhTimeCourse(cnt) = 1;
		ShutoffRate = InitialShutoffRate * (FeedbackGain * cnt)^Cooperativity;
		if (rand(1) < ShutoffRate)
			break;
		end
	end	

	RhTimeCourse = fft(RhTimeCourse);
	ReturnedCondition.EpochData(resp, :) = real(ifft(RhTimeCourse .* filt));
	
end
