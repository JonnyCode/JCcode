function DisplayEpochs(Condition)
%
% DisplayEpochs(Condition)
%
% DisplayEpochs.m
%	Superimposes all epochs in condition.
%
% Created 10/00 FMR
%  Revised GDF 12/01

NumEpochs = length(Condition.ExcludeEpochs)

hold on
for epoch = 1:NumEpochs
	if (Condition.ExcludeEpochs(epoch) == 0)
		plot(Condition.EpochData.Data(epoch, :));
	end
end
plot(Condition.AverageResponse, 'g');
hold off
