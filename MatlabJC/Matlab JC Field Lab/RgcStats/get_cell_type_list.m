function cell_type_list = get_cell_type_list(datarun) 

% JC 4/11/2016
% output: cell_type_list gives cell_type for each cell

cell_type_list = cell(1,length(datarun.spikes)) ;
for ct= 1:length(datarun.cell_types) ; % for each cell type
    if ~isempty(datarun.cell_types{ct}.cell_ids) ;
        for cid = 1:length(datarun.cell_types{ct}.cell_ids) ; % for each cell
            ci = get_cell_indices(datarun,datarun.cell_types{ct}.cell_ids(cid)) ;
            cell_type_list{ci} = datarun.cell_types{ct}.name ;
        end
    end
end

    