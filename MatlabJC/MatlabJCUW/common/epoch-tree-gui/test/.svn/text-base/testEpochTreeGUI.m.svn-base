function testEpochTreeGUI
%Run epochTreeGUI through some paces without an EpochTree fixture
%   Don't worry about actual graphical presentation

% should be constructable with EpochTree argument
%   should read selections into checkboxes
%   should be able to refresh checkboxes after tree changes
%%%SU
tree = getFixtureTree;
tree.custom.isSelected = false;
child = tree.children.firstValue;
child.custom.isSelected = false;
grandChild = child.children.firstValue;
grandChild.custom.isSelected = true;

gui = epochTreeGUI(tree);
drawnow;
trunkUnchecked = gui.treeBrowser.graphTree.trunk.isChecked;
childUnchecked = gui.treeBrowser.graphTree.trunk.getChild(1).isChecked;
grandChildChecked = gui.treeBrowser.graphTree.trunk.getChild(1).getChild(1).isChecked;

tree.custom.isSelected = true;
child = tree.children.firstValue;
child.custom.isSelected = false;
grandChild = child.children.firstValue;
grandChild.custom.isSelected = false;

epochTreeGUI.refreshTreeBrowserCallback(gui.treeBrowser.refresh, [], gui);
drawnow;
trunkChecked = gui.treeBrowser.graphTree.trunk.isChecked;
childUncheckedStill = gui.treeBrowser.graphTree.trunk.getChild(1).isChecked;
grandChildUnchecked = gui.treeBrowser.graphTree.trunk.getChild(1).getChild(1).isChecked;

close('all');
drawnow;
clear('gui', 'tree');
%%%TS isequal(trunkUnchecked, false)
%%%TS isequal(childUnchecked, false)
%%%TS isequal(grandChildChecked, true)
%%%TS isequal(trunkChecked, true)
%%%TS isequal(childUncheckedStill, false)
%%%TS isequal(grandChildUnchecked, false)


% widget click should callout to toggle isSelected
%   and recur down the tree
%%%SU
tree = getFixtureTree;
tree.custom.isSelected = false;
child = tree.children.firstValue;
child.custom.isSelected = false;
grandChild = child.children.firstValue;
grandChild.custom.isSelected = true;

gui = epochTreeGUI(tree);
drawnow;

% % "click" a child node's checkbox
widget = gui.treeBrowser.graphTree.widgetList.getValue(2);
widget.checkBox.isChecked = true;
graphicalTree.checkWidgetCallback(widget.checkBox, [], widget);
drawnow;

trunkUnchecked = tree.custom.isSelected;
childChecked = tree.children.firstValue.custom.isSelected;
grandChildChecked = tree.children.firstValue.children.firstValue.custom.isSelected;

widget.checkBox.isChecked = false;
graphicalTree.checkWidgetCallback(widget.checkBox, [], widget);
drawnow;

trunkUncheckedStill = tree.custom.isSelected;
childUnchecked = tree.children.firstValue.custom.isSelected;
grandChildUnChecked = tree.children.firstValue.children.firstValue.custom.isSelected;

close('all');
drawnow;
clear('gui', 'tree');
%%%TS isequal(trunkUnchecked, false)
%%%TS isequal(childChecked, true)
%%%TS isequal(grandChildChecked, true)
%%%TS isequal(trunkUncheckedStill, false)
%%%TS isequal(childUnchecked, false)
%%%TS isequal(grandChildUnChecked, false)


% should be able to interact with database
%   should be able to add tag to epochs and find that epochs all have that tag
%   should be able to remove that tag and find that epochs no longer have that tag
%   should be able to set includeinanalysis
%   should be able to unset includeinanalysis
%%%SU
tree = getFixtureTree;
gui = epochTreeGUI(tree);

% % add a tag
newTag = 'laser_tag';
set(gui.databaseInteraction.newTag, 'String', newTag);
epochTreeGUI.addTagCallback(gui.databaseInteraction.addTag, [], gui);
epochTreeGUI.refreshTagCallback(gui.databaseInteraction.refreshTag, [], gui);
tags = get(gui.databaseInteraction.existingTag, 'String');
tagWasAdded = any(strcmp(tags, newTag));

% % remove tag
epochTreeGUI.removeTagCallback(gui.databaseInteraction.addTag, [], gui);
epochTreeGUI.refreshTagCallback(gui.databaseInteraction.refreshTag, [], gui);
tags = get(gui.databaseInteraction.existingTag, 'String');
tagWasRemoved = isempty(tags) || ~any(strcmp(tags, newTag));

% % include in analysis
epochTreeGUI.includeInAnalysisCallback(gui.databaseInteraction.includeInAnalysis, [], gui, true);
sel = getTreeEpochs(tree, true);
allIncluded = true;
elements = sel.toCell;
for ii = 1:length(elements)
    ep = elements{ii};
    if ~ep.includeInAnalysis
        allIncluded = false;
        break
    end
end

% % exclude from analysis
epochTreeGUI.includeInAnalysisCallback(gui.databaseInteraction.includeInAnalysis, [], gui, false);
sel = getTreeEpochs(tree, true);
allExcluded = true;
elements = sel.toCell;
for ii = 1:length(elements)
    ep = elements{ii};
    if ep.includeInAnalysis
        allExcluded = false;
        break
    end
end

close('all');
clear('gui', 'tree', 'el', 'sel');
% %%TS tagWasAdded
% %%TS tagWasRemoved
% %%TS allIncluded
% %%TS allExcluded


% should optionally not drill down to individual Epochs
%   assume fixture tree has some Epochs
%%%SU
tree = getFixtureTree;

% % do drill down, how many graphical nodes?
guiWithEpochs = epochTreeGUI(tree);
nodesWithEpochs = guiWithEpochs.treeBrowser.graphTree.nodeList.length;

% % do not drill down, how many graphical nodes?
guiWithoutEpochs = epochTreeGUI(tree, 'noEpochs');
nodesWithoutEpochs = guiWithoutEpochs.treeBrowser.graphTree.nodeList.length;

close('all');
clear('guiWithEpochs', 'guiWithoutEpochs', 'tree');
%%%TS nodesWithEpochs > nodesWithoutEpochs


% should check for, obey alternate color, bg color, name for each node
%%%SU
alt.name = 'altName';
alt.color = [1 1 1]*.53442;
alt.backgroundColor = [1 1 1]*.78832;

tree = getFixtureTree;
tree.custom.display.alt = alt;

gui = epochTreeGUI(tree);
preservedAlt = tree.custom.display.alt;

browserNode = gui.treeBrowser.graphTree.trunk;
browserNodeName = browserNode.name;
browserNodeColor = browserNode.textColor;
browserNodeBgColor = browserNode.textBackgroundColor;

close('all');
clear('gui', 'tree');
%%%TS isequal(alt.name, preservedAlt.name)
%%%TS isequal(alt.color, preservedAlt.color)
%%%TS isequal(alt.backgroundColor, preservedAlt.backgroundColor)
%%%TS isequal(alt.name, browserNodeName)
%%%TS isequal(alt.color, browserNodeColor)
%%%TS isequal(alt.backgroundColor, browserNodeBgColor)


% This test can run and pass, but causes Matlab to crash sometimes...?

% % % % pan button in tree browser should toggle pan behavior of browser axes
% % % %%%SU
% % % gui = epochTreeGUI([]);
% % % cb = get(gui.treeBrowser.pan, 'Callback');
% % % set(gui.treeBrowser.pan, 'Value', true);
% % % feval(cb{1}, gui.treeBrowser.pan, [], cb{2:end})
% % % drawnow
% % % panObj = pan(gui.fig);
% % % panOn = panObj.Enable;
% % % 
% % % set(gui.treeBrowser.pan, 'Value', false);
% % % feval(cb{1}, gui.treeBrowser.pan, [], cb{2:end})
% % % drawnow
% % % panObj = pan(gui.fig);
% % % panOff = panObj.Enable;
% % % delete('panObj');
% % % close('all');
% % % clear('gui', 'tree');
% % % %%%TS strcmp(panOn, 'on')
% % % %%%TS strcmp(panOff, 'off')