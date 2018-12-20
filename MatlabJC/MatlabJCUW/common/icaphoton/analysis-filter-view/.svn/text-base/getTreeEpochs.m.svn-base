function el = getTreeEpochs(epochTree, onlySelected)
%Get one EpochList with all Epochs in tree (optionally, only selected)
%
%   el is a single auimodel.EpochList containing all Epochs from epochTree.
%   If onlySelected is true, el contains only those Epochs where
%   Epoch.isSelected is true.
%
%   epochTree is an auimodel.EpochTree with all your data.
%
%   getTreeEpochs is a utility used by epochTreeGUI and some
%   analysis-filter-view functions.
%
%%%SU
%   tree = getFixtureTree;
%   el = getTreeEpochs(tree);
%   allEpochs = el.length;
%
%   for ii = 1:length(el.elements)
%       ep = el.elements{ii};
%       ep.isSelected = false;
%   end
%   sel = getTreeEpochs(tree, true);
%   noEpochs = sel.length;
%
%   ep = el.elements{1};
%   ep.isSelected = true;
%   sel = getTreeEpochs(tree, true);
%   oneEpoch = sel.length;
%
%   clear tree el sel ep
%%%TS allEpochs > 0
%%%TS noEpochs == 0
%%%TS oneEpoch == 1

% benjamin.heasly@gmail.com
%   2 Feb. 2009

if nargin < 2 || isempty(onlySelected)
    onlySelected = false;
end

import auimodel.EpochList;
el = EpochList;
if isobject(epochTree)

    if epochTree.isLeaf
        
        if onlySelected
            for jj = 1:length(epochTree.epochList.elements)
                epoch = epochTree.epochList.elements{jj};
                
                if ~isempty(epoch.isSelected) && epoch.isSelected
                    el.append(epoch);
                end
            end
        else

            % use the leaf node's epochList verbatim
            el = epochTree.epochList;
        end
    else

        % build a grand EpochList
        for ii = 1:length(epochTree.leafNodes)
            
            % iter = epochTree.leafNodes{ii}.epochList.iterator('cell');
            % while iter.hasNext
            %     epoch = iter.nextValue;
            
            for jj = 1:length(epochTree.leafNodes{ii}.epochList.elements)
                epoch = epochTree.leafNodes{ii}.epochList.elements{jj};
                
                if ~onlySelected || (~isempty(epoch.isSelected) && epoch.isSelected)
                    el.append(epoch);
                end
            end
        end
    end
end