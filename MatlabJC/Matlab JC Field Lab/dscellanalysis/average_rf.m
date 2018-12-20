function [datarun] = average_rf (datarun, cellids)

%Averages 4 pixels at a  time
%Stores it in datarun.stas.rf_average field
%Input: datarun structure and cell ids for which you want to average receptive field
%Output: datarun structure with new field

datarun.stas.rf_average = cell(length(datarun.cell_ids),1);

for i = 1: length(cellids)
    cellIndex = get_cell_indices(datarun, cellids(1,i));
    recField = get_rf(datarun, cellids(1,i));
    recFieldAve = zeros(size(recField,1)/2, size(recField,2)/2);
    b = 1;
    for j = 1:2:size(recField,1)
            a = 1;
        for k = 1:2:size(recField,2)
            recFieldAve(b,a) = (recField(j,k)+recField(j, k+1)+recField(j+1,k)+recField(j+1,k+1))/4;
            a = a+1;
        end
        b = b+1;
    end
    datarun.stas.rf_average{cellIndex,1} = recFieldAve;
end

end

    