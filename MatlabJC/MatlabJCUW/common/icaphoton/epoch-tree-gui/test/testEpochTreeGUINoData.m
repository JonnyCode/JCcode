function testEpochTreeGUINoData
%Run epochTreeGUI through some paces without linking to any EpochTree
%   Don't worry about actual graphical presentation

% should be constructable with no arguments
% should not open a figure
%%%SU
gui = epochTreeGUI();
%%%TS isobject(gui)
%%%TS isa(gui, 'epochTreeGUI')
%%%TS isempty(gui.fig)


% should be constructable with empty epochTree argument
% should open figure with interface widgets in default states
%%%SU
gui = epochTreeGUI([]);
madeFigure = ishandle(gui.fig);
madeBrowser = isfield(gui.treeBrowser, 'graphTree') && isobject(gui.treeBrowser.graphTree);
madeDBInterac = isfield(gui.databaseInteraction, 'existingTag') && ishandle(gui.databaseInteraction.existingTag);
madeCanvas = isfield(gui.plottingCanvas, 'panel') && ishandle(gui.plottingCanvas.panel);
madeAnalysis = isfield(gui.analysisTools, 'functionMenu') && ishandle(gui.analysisTools.functionMenu);
close all
%%%TS isobject(gui)
%%%TS isa(gui, 'epochTreeGUI')
%%%TS madeFigure
%%%TS madeBrowser
%%%TS madeDBInterac
%%%TS madeCanvas
%%%TS madeAnalysis


% should be able to invoke each callback/static method without error
%   even though behavior is probably trivial
%%%SU
gui = epochTreeGUI([]);
event = [];
bNode = gui.treeBrowser.graphTree.trunk;
callbackInvokations = { ...
    {@epochTreeGUI.updateViewCallback, bNode, event, gui}, ...
    {@epochTreeGUI.addTagCallback, gui.databaseInteraction.addTag, event, gui}, ...
    {@epochTreeGUI.existingTagCallback, gui.databaseInteraction.existingTag, event, gui}, ...
    {@epochTreeGUI.includeInAnalysisCallback, gui.databaseInteraction.includeInAnalysis, event, gui}, ...
    {@epochTreeGUI.newTagCallback, gui.databaseInteraction.newTag, event, gui}, ...
    {@epochTreeGUI.removeTagCallback, gui.databaseInteraction.removeTag, event, gui}, ...
    {@epochTreeGUI.mouseCheckedNodeCallback, bNode, event, gui}, ...
    {@epochTreeGUI.applyToolCallback, gui.analysisTools.applyTool, event, gui}, ...
    {@epochTreeGUI.functionMenuCallback, gui.analysisTools.functionMenu, event, gui}, ...
    {@epochTreeGUI.panTreeBrowserCallback, gui.treeBrowser.pan, event, gui}, ...
    {@epochTreeGUI.setEpochNodeIsSelectedCallback, bNode, gui}, ...
    {@epochTreeGUI.invertTreeSelectionsCallback, gui.treeBrowser.invertTreeSelections, event, gui}, ...
    {@epochTreeGUI.invertEpochSelectionsCallback, gui.treeBrowser.invertEpochSelections, event, gui}, ...
    };

callbacksPassed = true;
for cb = callbackInvokations
    try
        feval(cb{1}{:});
    catch e
        callbacksPassed = false;
        disp(sprintf('%s failed:\n\t"%s"', func2str(cb{1}{1}), e.message))
    end
end
close all
%%%TS callbacksPassed


% pan button in tree browser should toggle pan behavior of browser axes
%%%SU
gui = epochTreeGUI([]);
cb = get(gui.treeBrowser.pan, 'Callback');
set(gui.treeBrowser.pan, 'Value', true);
feval(cb{1}, gui.treeBrowser.pan, [], cb{2:end})
drawnow
panObj = pan(gui.fig);
panOn = panObj.Enable;

set(gui.treeBrowser.pan, 'Value', false);
feval(cb{1}, gui.treeBrowser.pan, [], cb{2:end})
drawnow
panObj = pan(gui.fig);
panOff = panObj.Enable;
clear panObj %% Matlab complains when this handle reloaded for TS eval
close all
%%%TS strcmp(panOn, 'on')
%%%TS strcmp(panOff, 'off')