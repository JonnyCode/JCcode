function ReturnedCellInfo = RmCIData(CellInfo)
% ReturnedCellInfo = RmCIData(CellInfo)
% Function to remove the field CellInfo.EpochData from the structure 
% CellInfo.  Structure must be loaded into the workspace.  Function
% removes the epoch data, and returns ReturnedCellInfo.  ReturnedCellInfo
% can also be called CellInfo.  To load data in the structure,
% use LoadCIData

% Created as rmSCData February, 2001 MKMK
% New Version RmCIData released August, 2001, MKMK

ReturnedCellInfo = rmfield(CellInfo,'EpochData');

