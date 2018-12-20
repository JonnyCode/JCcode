function DataMatrix = DataStruct2Mat(DataStruct,FieldNamesCell)

%this function will make a matrix from the data in a structure with only
%vectors("DataStruct") with fieldnames in a cell array ("FieldNamesCell")

%JC 1/2/12

DataLength = length(DataStruct.(FieldNamesCell{1})) ; % length of data in first listed field name   
DataMatrix = nan(length(FieldNamesCell),DataLength) ; % prep matrix

for a=1:length(FieldNamesCell) ; % for every field
    if size(DataStruct.(FieldNamesCell{a})) ~= [DataLength,1] ; % if the data is not the expected size
        error('data size error')
    else % if the data is the right size
        DataMatrix(a,:) = DataStruct.(FieldNamesCell{a}) ;
    end
end

        