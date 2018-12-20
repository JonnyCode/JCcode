function NewEpochData = BaselineCorrectOvation(EpochData, StartBaseline, EndBaseline)
%
%	function NewEpochData = BaselineCorrectOvation(EpochData,StartBaseline, EndBaseline)
%
%	This function takes a small junk of baseline and subtracts off the average of that
%	baseline from the rest of the data.
%
%	Created: FMR 4/09
%
%%%SU 
%   load BaselineCorrectOvationTest/test
%   NewEpochData = BaselineCorrectOvation(EpochData,StartBaseline,EndBaseline);
%   for cnt = 1:size(NewEpochData, 1)
%       residual = mean(NewEpochData(cnt, StartBaseline:EndBaseline));
%   end
%   Indices = find(abs(residual) > 0);
%   Indices = find(abs(residual) > eps);
%%%TS isequal(length(Indices), 0)


[NumEpochs, EpochPts] = size(EpochData);

for epoch = 1:NumEpochs
    NewEpochData(epoch, :) = EpochData(epoch, :) - mean(EpochData(epoch, StartBaseline:EndBaseline));
end
