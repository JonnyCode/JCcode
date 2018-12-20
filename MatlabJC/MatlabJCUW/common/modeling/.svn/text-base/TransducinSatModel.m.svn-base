% Calculate series of dim flash responses for model with
% single step rhodopsin shutoff and saturation at 
% transducin.  Thus amount of available transducin decreases
% exponentially over time.  
%
% Created 7/01 FMR

function ReturnedCondition = TransducinSatModel(TransducinSatRate, ShutoffRate, Condition, Filter, NumResponses)

filt = fft(Filter);
ReturnedCondition = Condition;

for resp = 1:NumResponses

	RhTimeCourse(1:Condition.EpochPts) = 0;

	for cnt=1:Condition.EpochPts
		RhTimeCourse(cnt) = exp(-TransducinSatRate * cnt);
		if (rand(1) < ShutoffRate)
			break;
		end
	end	

	RhTimeCourse = fft(RhTimeCourse);
	ReturnedCondition.EpochData(resp, :) = real(ifft(RhTimeCourse .* filt));
	
end
