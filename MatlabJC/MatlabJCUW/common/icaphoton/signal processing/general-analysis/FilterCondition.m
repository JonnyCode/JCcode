function ReturnedCondition = FilterCondition(CellInfo, Condition, TimeConstant, Frequency, Stages);
% FilterCondition
%
% Filter data in all members of specified EpochCondition or FamilyCondition with 
% filter specified.

EpochPts = length(Condition(1).AverageResponse);
Time = 1:EpochPts;
Time = Time * FindSearchPara(Condition(1), 'SampInterv') * Condition(1).DecimatePts * 1e-6;
Filter = (-Time/TimeConstant).^(Stages-1) .* exp(-Time / TimeConstant) .* cos(2 * 3.14159 * Time * Frequency);
NumConds = length(Condition);

for cond = 1:NumConds
	for epoch = 1:length(Condition(cond).EpochNumbers)
		resp = Condition(cond).EpochData.Data(epoch, :);
		resp = conv(resp, Filter);
		Condition(cond).EpochData.Data(epoch, 1:EpochPts) = resp(1:EpochPts);
	end	
	Condition(cond) = AverageAndVariance(CellInfo, Condition(cond));
end

ReturnedCondition = Condition;
