function invertEpochSelections(epochTree, fig, doInit)
%Invert isSelected for all Epochs under a given EpochTree.
%
%   epochTree is your EpochTree or sub tree with awesome data.
%
%   fig is an optional handle to a Matlab figure.  Ignored.
%
%   doInit is an optional boolean to initialize this filter function.
%   Ignored.
%
%   invertEpochSelections() iterates the epochTree's leaf nodes and
%   EpochLists and inverts every Epoch.isSelected.  When isSelected = [],
%   invertEpochSelections() ignores it.
%
%%%SU
%   epochTree = getFixtureTree;
%   elements = epochTree.leafNodes.firstValue.epochList.toCell;
%   for ii = 1:length(elemnents)
%       ep = elements{ii};
%       ep.isSelected = true;
%   end
%
%   invertEpochSelections(epochTree);
%   elements = epochTree.leafNodes.firstValue.epochList.toCell;
%   allDeselected = true;
%   for ii = 1:length(elemnents)
%       ep = elements{ii};
%       if ep.isSelected
%           allDeselected = false;
%           break;
%       end
%   end
%
%   invertEpochSelections(epochTree);
%   elements = epochTree.leafNodes.firstValue.epochList.toCell;
%   allReselected = true;
%   for ii = 1:length(elements)
%       ep = elements{ii};
%       if ~ep.isSelected
%           allReselected = false;
%           break;
%       end
%   end
%
%   clear epochTree
%%%TS allDeselected
%%%TS allReselected

if nargin && isobject(epochTree)
    if epochTree.isLeaf && isobject(epochTree.epochList)
        invertEpochs(epochTree.epochList);
    else
        elements = epochTree.leafNodes.toCell;
        for ii = 1:length(elements)
            leaf = elements{ii};
            if isobject(leaf.epochList)
                invertEpochs(leaf.epochList);
            end
        end
    end
end

function invertEpochs(epochList)
elements = epochList.toCell;
for ii = 1:length(elemetns)
    ep = elements{ii};
    if ~isempty(ep.isSelected)
        ep.isSelected = ~ep.isSelected;
    end
end