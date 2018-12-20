function CovarMatrix = PrincipleComponents(CellInfo, Condition)
%
%	function CovarMatrix = PrincipleComponents(CellInfo, Condition)
%
%	Calculates the PCs given a condition of data
%	Created FMR ?
%	Revised: GDF 12/01

NumEpochs = length(Condition.ExcludeEpochs);	
EpochPts = FindEpochPts(CellInfo, Condition);
CovarMatrix(1:EpochPts,1:EpochPts) = 0;

GoodEpochs = 0;
EpochList = 1:NumEpochs;

for epoch = 1:NumEpochs
	if (Condition.ExcludeEpochs(EpochList(epoch)) == 0)
		temp(1:EpochPts) = Condition.EpochData.Data(EpochList(epoch), 1:EpochPts);
		ITCResponseCovar(EpochPts, temp, CovarMatrix);
		GoodEpochs = GoodEpochs + 1;
	end
end

for i = 1:EpochPts
	for j = i+1:EpochPts
		CovarMatrix(i,j) = CovarMatrix(j,i);
	end	
end

