function testEpochTreeGUI
%Run epochTreeGUI through some paces without an EpochTree fixture
%   Don't worry about actual graphical presentation

% should be constructable with EpochTree argument
%   should read selections into checkboxes
%   should be able to refresh checkboxes after tree changes
%%%SU
tree = getFixtureTree;
tree.custom.isSelected = false;
tree.children{1}.custom.isSelected = false;
tree.children{1}.children{1}.custom.isSelected = true;

gui = epochTreeGUI(tree);
drawnow;
trunkUnchecked = gui.treeBrowser.graphTree.trunk.isChecked;
childUnchecked = gui.treeBrowser.graphTree.trunk.children(1).isChecked;
grandChildChecked = gui.treeBrowser.graphTree.trunk.children(1).children(1).isChecked;

tree.custom.isSelected = true;
tree.children{1}.custom.isSelected = false;
tree.children{1}.children{1}.custom.isSelected = false;

epochTreeGUI.refreshTreeBrowserCallback(gui.treeBrowser.refresh, [], gui);
drawnow;
trunkChecked = gui.treeBrowser.graphTree.trunk.isChecked;
childUncheckedStill = gui.treeBrowser.graphTree.trunk.children(1).isChecked;
grandChildUnchecked = gui.treeBrowser.graphTree.trunk.children(1).children(1).isChecked;

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
tree.children{1}.custom.isSelected = false;
tree.children{1}.children{1}.custom.isSelected = true;

gui = epochTreeGUI(tree);
drawnow;

% % "click" a child node's checkbox
widget = gui.treeBrowser.graphTree.drawWidgets(2);
widget.checkBox.isChecked = true;
graphicalTree.checkWidgetCallback(widget.checkBox, [], widget);
drawnow;

trunkUnchecked = tree.custom.isSelected;
childChecked = tree.children{1}.custom.isSelected;
grandChildChecked = tree.children{1}.children{1}.custom.isSelected;

widget.checkBox.isChecked = false;
graphicalTree.checkWidgetCallback(widget.checkBox, [], widget);
drawnow;

trunkUncheckedStill = tree.custom.isSelected;
childUnchecked = tree.children{1}.custom.isSelected;
grandUnchildChecked = tree.children{1}.children{1}.custom.isSelected;

close('all');
drawnow;
clear('gui', 'tree');
%%%TS isequal(trunkUnchecked, false)
%%%TS isequal(childChecked, true)
%%%TS isequal(grandChildChecked, true)
%%%TS isequal(trunkUncheckedStill, false)
%%%TS isequal(childUnchecked, false)
%%%TS isequal(grandUnchildChecked, false)


% should be able to interact with database
%   should be able to add tag to epochs and find that epochs all have that tag
%   should be able to remove that tag and find that epochs no longer have that tag
%   should be able to set includeinanalysis
%   should be able to unset includeinanalysis
%%%SU
tree = getFixtureTree;
gui = epochTreeGUI(tree);

% % for speed, select only a few Epochs
el = getTreeEpochs(tree);
for ii = 1:9
    ep = el.elements{ii};
    ep.isSelected = true;
end
for ii = 10:el.length
    ep = el.elements{ii};
    ep.isSelected = false;
end

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
for ii = 1:sel.length
    if ~sel.elements{ii}.includeInAnalysis
        allIncluded = false;
        break
    end
end

% % exclude from analysis
epochTreeGUI.includeInAnalysisCallback(gui.databaseInteraction.includeInAnalysis, [], gui, false);
sel = getTreeEpochs(tree, true);
allExcluded = true;
for ii = 1:sel.length
    if sel.elements{ii}.includeInAnalysis
        allExcluded = false;
        break
    end
end

close('all');
clear('gui', 'tree', 'el', 'sel');
%%%TS tagWasAdded
%%%TS tagWasRemoved
%%%TS allIncluded
%%%TS allExcluded