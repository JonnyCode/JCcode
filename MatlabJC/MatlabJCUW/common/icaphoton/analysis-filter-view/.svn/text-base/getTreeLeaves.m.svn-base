function l = getTreeLeaves(epochTree, onlySelected)
%Get cell array of leaf nodes (optionally, only selected)
%
%   l is a cell array containing any leaf nodes of epochTree.  If
%   onlySelected is true, l contains only those leaf nodes where
%   EpochTree.custom.isSelected is true.
%
%   epochTree is an auimodel.EpochTree with all your data.
%
%   getTreeLeaves is a utility used by epochTreeGUI and some
%   analysis-filter-view functions.
%
%%%SU
%   tree = getFixtureTree;
%   for ii = 1:length(tree.leafNodes)
%       tree.leafNodes{ii}.custom.isSelected = false;
%   end
%
%   allLeaves = length(getTreeLeaves(tree));
%   noLeaves = length(getTreeLeaves(tree, true));
%   tree.leafNodes{1}.custom.isSelected = true;
%   oneLeaf = length(getTreeLeaves(tree, true));
%   clear tree;
%%%TS noLeaves == 0
%%%TS allLeaves > 0
%%%TS oneLeaf == 1

% benjamin.heasly@gmail.com
%   17 Feb. 2009

if nargin < 2 || isempty(onlySelected)
    onlySelected = false;
end

l = {};
if isobject(epochTree)
    if epochTree.isLeaf
        if ~onlySelected || (isfield(epochTree.custom, 'isSelected') && epochTree.custom.isSelected)
            % use a leaf node as-is
            l{1} = epochTree;
        end
    else
        % build a grand list of leaves
        for ii = 1:length(epochTree.leafNodes)
            if ~onlySelected || (isfield(epochTree.leafNodes{ii}.custom, 'isSelected') && epochTree.leafNodes{ii}.custom.isSelected)
                l{end+1,1} = epochTree.leafNodes{ii};
            end
        end
    end
end