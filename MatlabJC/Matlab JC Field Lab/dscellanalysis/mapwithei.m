function[masterind slaveind] = mapwithei(masterdata, slavedata, cellids, c)

[cell_list_map, failed_cells] = map_ei(masterdata, slavedata, 'master_cell_type', cellids, 'corr_threshold', c);
emptyCells = cellfun(@isempty,cell_list_map);
slaveind = get_cell_indices(slavedata, cell2mat(cell_list_map(1,(emptyCells==0))));
masterind =  find(emptyCells == 0);

end