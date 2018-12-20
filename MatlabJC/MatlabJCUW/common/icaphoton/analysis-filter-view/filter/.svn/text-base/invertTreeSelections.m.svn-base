function invertTreeSelections(epochTree, fig, doInit)
%Recursively invert isSelected noded of given EpochTree.
%
%   epochTree is your EpochTree or sub tree with awesome data.
%
%   fig is an optional handle to a Matlab figure.  Ignored.
%
%   doInit is an optional boolean to initialize this filter function.
%   Ignored.
%
%   invertTreeSelections() traverses the given epochTree's nodes and and
%   inverts every EpochTree.custom.isSelected.  When isSelected = [], or
%   DNE, invertTreeSelections() ignores it.
%
%%%SU
%   epochTree = getFixtureTree;
%   epochTree.custom.isSelected = true;
%   epochTree.children{1}.custom.isSelected = true;
%   epochTree.children{1}.children{1}.custom.isSelected = true;
%
%   invertTreeSelections(epochTree);
%   deselected(3) = epochTree.children{1}.children{1}.custom.isSelected;
%   deselected(2) = epochTree.children{1}.custom.isSelected;
%   deselected(1) = epochTree.custom.isSelected;
%
%   invertTreeSelections(epochTree);
%   reselected(3) = epochTree.children{1}.children{1}.custom.isSelected;
%   reselected(2) = epochTree.children{1}.custom.isSelected;
%   reselected(1) = epochTree.custom.isSelected;
%
%   clear epochTree
%%%TS all(deselected == false)
%%%TS all(reselected == true)

if nargin && isobject(epochTree)
    invertNode(epochTree);
end


function invertNode(epochTree)
if isfield(epochTree.custom, 'isSelected') && ~isempty(epochTree.custom.isSelected)
    epochTree.custom.isSelected = ~epochTree.custom.isSelected;
end

for ii = 1:length(epochTree.children)
    invertNode(epochTree.children{ii});
end