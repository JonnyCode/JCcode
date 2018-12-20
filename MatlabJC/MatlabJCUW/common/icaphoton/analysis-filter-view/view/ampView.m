function ampView(epochTree, fig, doInit);
%Simple view function for superimposing amp data from selected nodes

% benjamin.heasly@gmail.com
%   2 Feb. 2009

if nargin && isobject(epochTree)
    if epochTree.isLeaf
        info = sprintf('(leaf)%s', getEpochTreeSplitString(epochTree));
    else
        info = sprintf('(%d leaves)%s', length(epochTree.leafNodes), getEpochTreeSplitString(epochTree));
    end
else
    
    info = 'no EpochTree given';
    return
end

if nargin < 2
    disp(sprintf('%s needs a figure', mfilename));
    return
end

if nargin < 3
    doInit = false;
end

figData = get(fig, 'UserData');

% init when told, or when this is not the 'current function' of the figure
if doInit || ~isfield(figData, 'currentFunction') || ~strcmp(figData.currentFunction, mfilename)
    
    % Don't clf!  This will clear the parent figure, which could be a whole GUI
    delete(get(fig, 'Children'));
    figData.currentFunction = mfilename;
    figData.ax = axes('Parent', fig);
    set(fig, 'UserData', figData);
end

cla(figData.ax);
title(figData.ax, 'getting data...');
drawnow

leaves = getTreeLeaves(epochTree, true);
streamName = 'no data';
for ii = 1:length(leaves)
    if ~isempty(leaves{ii}.epochList)
        col = dec2bin(mod(ii, 7),3)=='1';
        streamName = leaves{ii}.epochList.responseStreamNames{1};
        streamMean = mean(leaves{ii}.epochList.responsesByStreamName(streamName));
        line(1:length(streamMean), streamMean, 'Color', col, 'Parent', figData.ax);
        drawnow
    end
end
% just use the streamName from the last EpochList...might be wrong...
ylabel(figData.ax, sprintf('%s mean', streamName));
title(figData.ax, info);