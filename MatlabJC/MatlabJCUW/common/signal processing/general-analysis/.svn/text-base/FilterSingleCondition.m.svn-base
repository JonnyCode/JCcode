function ReturnedCondition = FilterSingleCondition(CellInfo, Condition, TimeConstant, Frequency);
% FilterSingleCondition
%
% Low-pass filter single condition with filter specified.

EpochPts = length(Condition.AverageResponse);
Time = 1:EpochPts;
Time = Time * FindSearchPara(Condition, 'SampInterv') * Condition.DecimatePts * 1e-6;
Filter = exp(-Time / TimeConstant) .* cos(2 * 3.14159 * Time * Frequency);
NumConds = length(Condition);

for epoch = 1:length(Condition.ExcludeEpochs)
	resp = Condition.EpochData.Data(epoch, :);
	resp = conv(resp, Filter);
	Condition.EpochData.Data(epoch, 1:EpochPts) = resp(1:EpochPts);
end	
Condition = AverageAndVariance(CellInfo, Condition);

ReturnedCondition = Condition;
