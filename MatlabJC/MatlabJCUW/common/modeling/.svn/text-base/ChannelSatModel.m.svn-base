% Calculate series of dim flash responses for single step model for
% rhodopsin shutoff with saturation at channel.  First generate time course of 
% rhodopsin activity and convolve with linear filter approximation to transduction
% cascade to generate uncompressed response.  Pass this through exponential compression.
% General idea is that exp(-x) ~ 1-x for x small, hence no compression.  But as x 
% gets larger compression kicks in.

% Created 7/01 FMR

function ReturnedCondition = ChannelSatModel(ExpFact, ShutoffRate, Condition, Filter, NumResponses)

filt = fft(Filter);
ReturnedCondition = Condition;

for resp = 1:NumResponses

	RhTimeCourse(1:Condition.EpochPts) = 0;

	for cnt=1:Condition.EpochPts
		RhTimeCourse(cnt) = 1;
		if (rand(1) < ShutoffRate)
			break;
		end
	end	

	RhTimeCourse = fft(RhTimeCourse);
	UnCompressedResponse = -real(ifft(RhTimeCourse .* filt));
	UnCompressedResponse = UnCompressedResponse .* (sign(UnCompressedResponse) + 1)/2;
	ReturnedCondition.EpochData(resp, :) = -(1-exp(-UnCompressedResponse / ExpFact));
end
