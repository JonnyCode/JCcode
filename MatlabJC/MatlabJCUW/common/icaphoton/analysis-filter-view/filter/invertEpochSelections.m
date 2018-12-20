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
%   eps = [epochTree.leafNodes{1}.epochList.elements{:}];
%   [eps.isSelected] = deal(true);
%
%   invertEpochSelections(epochTree);
%   deselected = [eps.isSelected];
%
%   invertEpochSelections(epochTree);
%   reselected = [eps.isSelected];
%
%   clear epochTree
%%%TS all(deselected == false)
%%%TS all(reselected == true)

if nargin && isobject(epochTree)
    if epochTree.isLeaf && isobject(epochTree.epochList)
        invertEpochs(epochTree.epochList);
    else
        for ii = 1:length(epochTree.leafNodes)
            if isobject(epochTree.leafNodes{ii}.epochList)
                invertEpochs(epochTree.leafNodes{ii}.epochList);
            end
        end
    end
end

function invertEpochs(epochList)
for ii = 1:epochList.length
    ep = epochList.elements{ii};
    if ~isempty(ep.isSelected)
        ep.isSelected = ~ep.isSelected;
    end
end