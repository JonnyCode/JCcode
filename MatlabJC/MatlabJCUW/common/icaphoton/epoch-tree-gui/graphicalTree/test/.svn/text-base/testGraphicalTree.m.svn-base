function testGraphicalTree
%Run graphicalTree through some paces.
%   Don't worry about actial graphical presentation
%   Implicitly test graphicalTreeNode and graphicalTreeNodeWidget

% should be constructable with no arguments
%%%SU
tree = graphicalTree();
%%%TS isobject(tree)
%%%TS isa(tree, 'graphicalTree')


% should be constructable with axes and name arguments
%   should precreate a bunch of drawable widgets and axes children
%   should create a trunk node with the same name
%%%SU
ax = axes;
name = 'larch';
tree = graphicalTree(ax, name);
widgets = get(ax, 'Children');
trunk = tree.trunk;
close all
%%%TS isobject(tree)
%%%TS isa(tree, 'graphicalTree')
%%%TS ~isempty(widgets)
%%%TS isobject(trunk)
%%%TS isa(trunk, 'graphicalTreeNode')
%%%TS strcmp(trunk.name, name)
%%%TS isobject(tree.drawWidgets)
%%%TS isa(tree.drawWidgets, 'graphicalTreeNodeWidget')


% tree's axes should adjust limits when figure resized or
fig = figure;
ax = axes('Parent', fig);
name = 'larch';
tree = graphicalTree(ax, name);
set(fig, 'Position', [1 1 1 1]*50);
drawnow
smallX = get(ax, 'XLim');
smallY = get(ax, 'YLim');
set(fig, 'Position', [1 1 1 1]*500);
drawnow
bigX = get(ax, 'XLim');
bigY = get(ax, 'YLim');
close all
%%%TS ~isequal(smallX, bigX)
%%%TS ~isequal(smallY, bigY)


% should be able to add loads of nodes to the tree
%   should create enough graphicalTreeNodeWidget to display expanded tree
%%%SU
ax = axes;
name = 'larch';
tree = graphicalTree(ax, name);
tree.trunk.isExpanded = true;
n = 3;
for ii = 1:n
    branch = tree.newNode(tree.trunk, 'branch');
    branch.isExpanded = true;
    for jj = 1:n
        twig = tree.newNode(branch, 'twig');
        twig.isExpanded = true;
        for kk = 1:n
            leaf = tree.newNode(twig, 'leaf');
            leaf.isExpanded = true;
        end
    end
end
tree.draw;
close all
%%%TS isequal(tree.trunk.numDescendants, n+n^2+n^3);
%%%TS length(tree.drawWidgets) > tree.trunk.numDescendants


% graphicalTreeNodeWidget.expandBox.box.ButtonDownFcn should invoke
% graphicalTree.expandWidgetCallback.  This should toggle node.isSelected,
% should set tree.CurrentNode, should redraw the tree.
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
tree.trunk.isExpanded = true;
tree.newNode(tree.trunk, 'leaf');
tree.draw;
expanded = tree.trunk.isExpanded;
twoDrawn = tree.drawCount;

box = tree.drawWidgets(1).expandBox.box;
cb = get(box, 'ButtonDownFcn');
feval(cb{1}, box, [], cb{2});
unexpanded = tree.trunk.isExpanded;
oneDrawn = tree.drawCount;

feval(cb{1}, box, [], cb{2});
reexpanded = tree.trunk.isExpanded;
twoDrawn2 = tree.drawCount;
close all
%%%TS isequal(expanded, true)
%%%TS isequal(unexpanded, false)
%%%TS isequal(reexpanded, true)
%%%TS isequal(tree.currentNode, tree.trunk)
%%%TS isequal(twoDrawn, 2)
%%%TS isequal(oneDrawn, 1)
%%%TS isequal(twoDrawn2, 2)


% graphicalTreeNodeWidget.expandBox.box.ButtonDownFcn should invoke
% graphicalTree.expandNodeButtonCallback, of the form {@foo}
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
cbInfo = 'callback happened';
foo = @(node, event) set(ax, 'UserData', cbInfo);
tree.expandNodeButtonCallback = {foo};
box = tree.drawWidgets(1).expandBox.box;
cb = get(box, 'ButtonDownFcn');
feval(cb{1}, box, [], cb{2});
trunkData = get(ax, 'UserData');
close all
%%%TS isequal(trunkData, cbInfo)


% graphicalTreeNodeWidget.expandBox.box.ButtonDownFcn should invoke
% graphicalTree.expandNodeButtonCallback, of the form {@foo, bar}
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
bar = 'callback happened';
foo = @(node, event, cbinfo) set(ax, 'UserData', cbinfo);
tree.expandNodeButtonCallback = {foo, bar};
box = tree.drawWidgets(1).expandBox.box;
cb = get(box, 'ButtonDownFcn');
feval(cb{1}, box, [], cb{2});
trunkData = get(ax, 'UserData');
close all
%%%TS isequal(trunkData, bar)


% graphicalTreeNodeWidget.checkBox.box.ButtonDownFcn should invoke
% graphicalTree.checkWidgetCallback.  This should toggle node.isChecked,
% should set tree.CurrentNode.
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
tree.trunk.isChecked = true;
tree.draw;
checked = tree.trunk.isChecked;

box = tree.drawWidgets(1).checkBox.box;
cb = get(box, 'ButtonDownFcn');
feval(cb{1}, box, [], cb{2});
unchecked = tree.trunk.isChecked;

feval(cb{1}, box, [], cb{2});
rechecked = tree.trunk.isChecked;
close all
%%%TS isequal(checked, true)
%%%TS isequal(unchecked, false)
%%%TS isequal(rechecked, true)
%%%TS isequal(tree.currentNode, tree.trunk)


% graphicalTreeNodeWidget.checkBox.box.ButtonDownFcn should invoke
% graphicalTree.checkNodeButtonCallback, of the form {@foo}
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
cbInfo = 'callback happened';
foo = @(node, event) set(ax, 'UserData', cbInfo);
tree.checkNodeButtonCallback = {foo};
box = tree.drawWidgets(1).checkBox.box;
cb = get(box, 'ButtonDownFcn');
feval(cb{1}, box, [], cb{2});
trunkData = get(ax, 'UserData');
close all
%%%TS isequal(trunkData, cbInfo)


% graphicalTreeNodeWidget.checkBox.box.ButtonDownFcn should invoke
% graphicalTree.checkNodeButtonCallback, of the form {@foo, bar}
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
bar = 'callback happened';
foo = @(node, event, cbinfo) set(ax, 'UserData', cbinfo);
tree.checkNodeButtonCallback = {foo, bar};
box = tree.drawWidgets(1).checkBox.box;
cb = get(box, 'ButtonDownFcn');
feval(cb{1}, box, [], cb{2});
trunkData = get(ax, 'UserData');
close all
%%%TS isequal(trunkData, bar)


% changing isChecked on a node should invoke
% graphicalTree.nodeBecameCheckedCallback, of the form {@foo}
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
cbInfo = 'callback happened';
foo = @(node) set(ax, 'UserData', cbInfo);
tree.nodeBecameCheckedCallback = {foo};
tree.trunk.setChecked(true);
tree.trunk.setChecked(false);
trunkData = get(ax, 'UserData');
close all
%%%TS isequal(trunkData, cbInfo)


% changing isChecked on a node should invoke
% graphicalTree.nodeBecameCheckedCallback, of the form {@foo, bar}
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
bar = 'callback happened';
foo = @(node, cbinfo) set(ax, 'UserData', cbinfo);
tree.nodeBecameCheckedCallback = {foo, bar};
tree.trunk.setChecked(true);
tree.trunk.setChecked(false);
trunkData = get(ax, 'UserData');
close all
%%%TS isequal(trunkData, bar)


% graphicalTreeNodeWidget.nameText.ButtonDownFcn should invoke
% graphicalTree.clickWidgetCallback.  This should set tree.CurrentNode.
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
text = tree.drawWidgets(1).nameText;
cb = get(text, 'ButtonDownFcn');
feval(cb{1}, text, [], cb{2});
close all
%%%TS isequal(tree.currentNode, tree.trunk)


% graphicalTreeNodeWidget.nameText.ButtonDownFcn should invoke
% graphicalTree.clickNodeButtonCallback, of the form {@foo}
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
cbInfo = 'callback happened';
foo = @(node, event) set(ax, 'UserData', cbInfo);
tree.clickNodeButtonCallback = {foo};
text = tree.drawWidgets(1).nameText;
cb = get(text, 'ButtonDownFcn');
feval(cb{1}, text, [], cb{2});
trunkData = get(ax, 'UserData');
close all
%%%TS isequal(trunkData, cbInfo)


% graphicalTreeNodeWidget.nameText.ButtonDownFcn should invoke
% graphicalTree.clickNodeButtonCallback, of the form {@foo, bar}
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
bar = 'callback happened';
foo = @(node, event, cbinfo) set(ax, 'UserData', cbinfo);
tree.clickNodeButtonCallback = {foo, bar};
text = tree.drawWidgets(1).nameText;
cb = get(text, 'ButtonDownFcn');
feval(cb{1}, text, [], cb{2});
trunkData = get(ax, 'UserData');
close all
%%%TS isequal(trunkData, bar)


% should be able to forget callbacks and make nodes forget userData
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
branch = tree.newNode(tree.trunk, 'branch');
leaf = tree.newNode(branch, 'leaf');
tree.trunk.userData = 'data';
branch.userData = 'data';
leaf.userData = 'data';
tree.expandNodeButtonCallback = {@disp};
tree.checkNodeButtonCallback = {@disp};
tree.clickNodeButtonCallback = {@disp};
tree.nodeBecameCheckedCallback = {@disp};
nodesHadData = ~isempty(tree.trunk.userData) ...
    && ~isempty(branch.userData) ...
    && ~isempty(leaf.userData);
treeHadCallbacks = ~isempty(tree.expandNodeButtonCallback) ...
    && ~isempty(tree.checkNodeButtonCallback) ...
    && ~isempty(tree.clickNodeButtonCallback) ...
    && ~isempty(tree.nodeBecameCheckedCallback);
tree.forgetExternalReferences;
nodesForgotData = isempty(tree.trunk.userData) ...
    && isempty(branch.userData) ...
    && isempty(leaf.userData);
treeForgotCallbacks = isempty(tree.expandNodeButtonCallback) ...
    && isempty(tree.checkNodeButtonCallback) ...
    && isempty(tree.clickNodeButtonCallback) ...
    && isempty(tree.nodeBecameCheckedCallback);
close all
%%%TS nodesHadData
%%%TS treeHadCallbacks
%%%TS nodesForgotData
%%%TS treeForgotCallbacks

% should be able to forget tree nodes and draw widgets,
% and make nodes forget each other
ax = axes;
tree = graphicalTree(ax, 'larch');
trunk = tree.trunk;
branch = tree.newNode(trunk, 'branch');
leaf = tree.newNode(branch, 'leaf');
childrenHadParent = isobject(branch.parent) && isobject(leaf.parent);
parentsHadChild = isobject(trunk.children) && isobject(branch.children);
treeHadTrunk = isobject(tree.trunk);
treeHadDrawWidgets = isobject(tree.drawWidgets);
tree.forgetInternalReferences;
childrenForgotParent = isempty(branch.parent) && isempty(leaf.parent);
parentsForgotChild = isempty(trunk.children) && isempty(branch.children);
treeForgotTrunk = isempty(tree.trunk);
treeForgotDrawWidgets = isempty(tree.drawWidgets);
close all
%%%TS childrenHadParent
%%%TS parentsHadChild
%%%TS treeHadTrunk
%%%TS treeHadDrawWidgets
%%%TS childrenForgotParent
%%%TS parentsForgotChild
%%%TS treeForgotTrunk
%%%TS treeForgotDrawWidgets

% graphicalTreeNodeWidgets should recolor their text to match the color of
% the model graphicalTreeNode
ax = axes;
tree = graphicalTree(ax, 'larch');

% red trunk
red = [1 0 0];
tree.trunk.textColor = red;
tree.draw;
trunkRedColor = get(tree.drawWidgets(1).nameText, 'Color');

% blue trunk
blue = [0 0 1];
tree.trunk.textColor = blue;
tree.draw;
trunkBlueColor = get(tree.drawWidgets(1).nameText, 'Color');

close all
%%%TS isequal(trunkRedColor, red)
%%%TS isequal(trunkBlueColor, blue)