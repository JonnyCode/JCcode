function epochLists = getSelectedLeafEpochs(epochTree)
%Get lists of epochs from selected leaf nodes
%
%   epochLists is a cell array of auimodel.EpochList, one list for each
%   leaf node that's selected.
%
%   epochTree is an auimodel.EpochTree with all your data.
%
%   This is a utility used by epochTreeGUI and some analysis-filter-view
%   functions.

% benjamin.heasly@gmail.com
%   2 Feb. 2009

% To Do: change EpochTree.custom.isSelected to EpochTree.isSelected

epochLists = {};
if isobject(epochTree)
    leaves = epochTree.leafNodes;
    for ii = 1:length(leaves)
        if isfield(leaves{ii}.custom, 'isSelected') && leaves{ii}.custom.isSelected
            epochLists{end+1,1} = leaves{ii}.epochList;
        end
    end
end