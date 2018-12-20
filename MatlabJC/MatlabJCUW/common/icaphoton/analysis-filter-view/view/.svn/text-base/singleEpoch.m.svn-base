function singleEpoch(epochTree, fig, doInit);
%Show one epoch at at time, for one leaf at a time

% benjamin.heasly@gmail.com
%   2 Feb. 2009

if ~nargin && ~isobject(epochTree)
    disp(sprintf('%s needs an EpochTree', mfilename));
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
    
    % create new panel and a slider
    delete(get(fig, 'Children'));
    figData.currentFunction = mfilename;
    figData.panel = uipanel('Parent', fig, ...
        'Units', 'normalized', ...
        'Position', [0 .1 1 .9]);
    figData.next = uicontrol('Parent', fig, ...
        'Units', 'normalized', ...
        'Position', [0 0 1 .1], ...
        'Style', 'slider', ...
        'Callback', {@plotNextEpoch, figData.panel});
    set(fig, 'UserData', figData);
end

% get all epochs from tree
el = getTreeEpochs(epochTree);
n = el.length;
if n
    set(figData.next, ...
        'Enable', 'on', ...
        'UserData', el, ...
        'Min',  1, ...
        'Max',  n+eps, ...
        'SliderStep', [1/n, 1/n], ...
        'Value', 1);
    plotNextEpoch(figData.next, [], figData.panel);
else
    delete(get(figData.panel, 'Children'));
    set(figData.next, 'Enable', 'off');
end

function plotNextEpoch(slider, event, panel)
set(panel, 'Title', 'getting data...');
drawnow;

ind = round(get(slider, 'Value'));
el = get(slider, 'UserData');

ep = el.elements{ind};
if ep.isSelected
    selStr = 'selected';
else
    selStr = 'not selected';
end
datStr = datestr(ep.startDate, 31);

resps = sort(fieldnames(ep.responses));
nResp = length(resps);
for ii = 1:nResp
    sp = subplot(nResp, 1, ii, 'Parent', panel);
    cla(sp);
    respData = ep.responses.(resps{ii}).data;
    if ischar(respData)
        % dont' break when lazyLoads disabled
        ylabel(sp, respData);
    else
        line(1:length(respData), respData, 'Parent', sp);
        ylabel(sp, resps{ii});
    end
end
set(panel, 'Title', sprintf('epoch %d of %d (%s)(%s)', ind, el.length, datStr, selStr));