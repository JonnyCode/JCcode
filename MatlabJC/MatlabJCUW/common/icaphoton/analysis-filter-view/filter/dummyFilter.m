function demoFilter(epochTree, fig, doInit);
%Placeholder for a filter function
%
%   Don't try to do any data work, just acknowledge the epochTree, the fig,
%   and the doInit flag.

if nargin && isobject(epochTree)
    info = sprintf('EpochTree with %d leaf nodes', length(epochTree.leafNodes));
else
    info = 'no EpochTree given';
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
    
    figData.treeInfo = uicontrol( ...
        'Parent', fig, ...
        'Style', 'text', ...
        'String', info, ...
        'Units', 'normalized', ...
        'Position', [0 .4 1 .1]);
    
    figData.callingInfo = uicontrol( ...
        'Parent', fig, ...
        'Style', 'text', ...
        'String', sprintf('%s initialized', mfilename), ...
        'Units', 'normalized', ...
        'Position', [0 .6 1 .1]);
    
    figData.timesCalled = 1;
    
    set(fig, 'UserData', figData);
else
    
    figData.timesCalled = figData.timesCalled + 1;
    set(figData.callingInfo, ...
        'String', sprintf('%s called %d times', mfilename, figData.timesCalled));
    set(fig, 'UserData', figData);
    
end