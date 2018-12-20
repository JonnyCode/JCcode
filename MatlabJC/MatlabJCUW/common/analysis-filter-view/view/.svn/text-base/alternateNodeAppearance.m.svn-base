function alternateNodeAppearance(epochTree, fig, doInit);
%Specify alternate node appearance for tree browser

% benjamin.heasly@gmail.com
%   26 June 2009

if ~nargin && ~isobject(epochTree)
    return
end
if nargin < 2
    return
end
if nargin < 3
    doInit = false;
end

% init when told, or when this is not the 'current function' of the figure
figData = get(fig, 'UserData');
if doInit || ~isfield(figData, 'currentFunction') || ~strcmp(figData.currentFunction, mfilename)
    figData.currentFunction = mfilename;
    
    tip = sprintf('%s\n', ...
        'Specify appearance for the current tree node:', ...
        '    - alternate name as a string.', ...
        '    - alternate color as a Matlab ColorSpec (''[1 0 0]'', ''red'', ''r'').', ...
        '    - empty field to use default.', ...
        '', ...
        'Press return and click [refresh] in the tree browser.');
    
    figData.tips = uicontrol('Parent', fig', ...
        'Units', 'normalized', ...
        'Position', [.2 0, .6 .3], ...
        'Style', 'text', ...
        'String', tip, ...
        'HorizontalAlignment', 'left');
    
    figData.displayTable = uitable('Parent', fig', ...
        'Units', 'normalized', ...
        'Position', [.2 .4, .6 .4], ...
        'RowName', [], ...
        'Data', [], ...
        'ColumnName', {'appearance', 'default', 'alternate'}, ...
        'ColumnFormat', {'char', 'char', 'char'}, ...
        'ColumnEditable', [false, false, true]);
    
    set(fig, 'UserData', figData, 'ResizeFcn', {@canvasResizeFcn, figData});
    canvasResizeFcn(fig, [], figData);
end

% populate table with display data for current node
fields = fieldnames(epochTree.custom.display.alt);
alternates = struct2cell(epochTree.custom.display.alt);
defaults = struct2cell(rmfield(epochTree.custom.display, 'alt'));
for ii = 1:length(fields)
    if isnumeric(defaults{ii})
        defaults{ii} = num2str(defaults{ii});
    end
    if isnumeric(alternates{ii})
        alternates{ii} = num2str(alternates{ii});
    end
end
set(figData.displayTable, ...
    'Data', cat(2, fields, defaults, alternates), ...
    'CellEditCallback', {@enterAlternateAppearanceFcn, epochTree});

function canvasResizeFcn(fig, event, figData)
% set infoTable column widths proportionally
oldUnits = get(figData.displayTable, 'Units');
set(figData.displayTable, 'Units', 'pixels');
tablePos = get(figData.displayTable, 'Position');
set(figData.displayTable, 'Units', oldUnits);
set(figData.displayTable, 'ColumnWidth', ...
    {.3*tablePos(3), .3*tablePos(3), .3*tablePos(3)});

function enterAlternateAppearanceFcn(table, event, epochTree)
data = get(table, 'Data');
field = data{event.Indices(1), 1};
switch field
    case {'color', 'backgroundColor'}
        val = str2num(event.NewData);
        if isempty(val)
            val = event.NewData;
        end
        
    otherwise
        val = event.NewData;
end
epochTree.custom.display.alt.(field) = val;