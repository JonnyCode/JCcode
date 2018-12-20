function testGraphicalTree
%Run graphicalTree through some paces.
%   Don't worry about actial graphical presentation
%   Implicitly test graphicalTreeNode and graphicalTreeNodeWidget

% should be constructable with no arguments
%%%SU
tree = graphicalTree();
%%%TS isobject(tree)
%%%TS isa(tree, 'graphicalTree')


% tree should link up parent and child nodes
%%%SU
tree = graphicalTree();
child = tree.newNode(tree.trunk, 'child');
%%%TS isobject(child)
%%%TS isa(child, 'graphicalTreeNode')
%%%TS strcmp(child.name, 'child')
%%%TS isequal(tree.trunk.getChild(1), child)
%%%TS isequal(child.getParent, tree.trunk)


% parent nodes should account for descendants
% nodes should know how deep they are
%%%SU
tree = graphicalTree();
child1 = tree.newNode(tree.trunk, 'child1');
child2 = tree.newNode(tree.trunk, 'child2');
grandchild = tree.newNode(child1, 'grandchild');
%%%TS isequal(tree.trunk.numChildren, 2)
%%%TS isequal(tree.trunk.numDescendants, 3)
%%%TS isequal(child1.numDescendants, 1)
%%%TS isequal(tree.trunk.depth, 0)
%%%TS isequal(child1.depth, 1)
%%%TS isequal(child2.depth, 1)
%%%TS isequal(grandchild.depth, 2)

% nodes should recursively check descendants
%%%SU
tree = graphicalTree();
tree.trunk.isChecked = false;
child = tree.newNode(tree.trunk, 'child');
child.isChecked = false;
grandchild = tree.newNode(child, 'grandchild');
grandchild.isChecked = false;
tree.trunk.recursiveCheck = true;
tree.trunk.setChecked(true);
%%%TS tree.trunk.isChecked
%%%TS child.isChecked
%%%TS grandchild.isChecked

% nodes should recursively uncheck descendants
%%%SU
tree = graphicalTree();
tree.trunk.isChecked = true;
child = tree.newNode(tree.trunk, 'child');
child.isChecked = true;
grandchild = tree.newNode(child, 'grandchild');
grandchild.isChecked = true;
tree.trunk.recursiveCheck = true;
tree.trunk.setChecked(false);
%%%TS ~tree.trunk.isChecked
%%%TS ~child.isChecked
%%%TS ~grandchild.isChecked

% nodes should recursively account for checked descendants
%%%SU
tree = graphicalTree();
child = tree.newNode(tree.trunk, 'child');
grandchild = tree.newNode(child, 'grandchild');
tree.trunk.setChecked(false);
noChecks = tree.trunk.numCheckedDescendants;
grandchild.setChecked(true);
oneCheck = tree.trunk.numCheckedDescendants;
tree.trunk.setChecked(true);
twoChecks = tree.trunk.numCheckedDescendants;
child.setChecked(false);
noChecks2 = tree.trunk.numCheckedDescendants;
%%%TS isequal(noChecks, 0)
%%%TS isequal(oneCheck, 1)
%%%TS isequal(twoChecks, 2)
%%%TS isequal(noChecks2, 0)

% recursive accounting should be able to ignore nodes with
% isVisibleDescandant = false;
%%%SU
tree = graphicalTree();
child = tree.newNode(tree.trunk, 'child');
grandchild = tree.newNode(child, 'grandchild', false);
oneVisibleDescandant = tree.trunk.numDescendants;

tree.trunk.setChecked(false);
noChecks = tree.trunk.numCheckedDescendants;

grandchild.setChecked(true);
noChecks2 = tree.trunk.numCheckedDescendants;

tree.trunk.setChecked(true);
oneCheck = tree.trunk.numCheckedDescendants;

child.setChecked(false);
noChecks3 = tree.trunk.numCheckedDescendants;
%%%TS oneVisibleDescandant == 1
%%%TS noChecks == 0
%%%TS noChecks2 == 0
%%%TS oneCheck == 1
%%%TS noChecks3 == 0


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
%%%TS isobject(tree.widgetList.getValue(1))
%%%TS isa(tree.widgetList.getValue(1), 'graphicalTreeNodeWidget')


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
%%%TS isequal(tree.trunk.numDescendants, n+n^2+n^3)
%%%TS tree.widgetList.length > tree.trunk.numDescendants


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

box = tree.widgetList.getValue(1).expandBox.box;
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
%%%TS isequal(tree.getCurrentNode, tree.trunk)
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
box = tree.widgetList.getValue(1).expandBox.box;
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
box = tree.widgetList.getValue(1).expandBox.box;
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

box = tree.widgetList.getValue(1).checkBox.box;
cb = get(box, 'ButtonDownFcn');
feval(cb{1}, box, [], cb{2});
unchecked = tree.trunk.isChecked;

feval(cb{1}, box, [], cb{2});
rechecked = tree.trunk.isChecked;
close all
%%%TS isequal(checked, true)
%%%TS isequal(unchecked, false)
%%%TS isequal(rechecked, true)
%%%TS isequal(tree.getCurrentNode, tree.trunk)


% graphicalTreeNodeWidget.checkBox.box.ButtonDownFcn should invoke
% graphicalTree.checkNodeButtonCallback, of the form {@foo}
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
cbInfo = 'callback happened';
foo = @(node, event) set(ax, 'UserData', cbInfo);
tree.checkNodeButtonCallback = {foo};
box = tree.widgetList.getValue(1).checkBox.box;
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
box = tree.widgetList.getValue(1).checkBox.box;
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
text = tree.widgetList.getValue(1).nameText;
cb = get(text, 'ButtonDownFcn');
feval(cb{1}, text, [], cb{2});
close all
%%%TS isequal(tree.getCurrentNode, tree.trunk)


% graphicalTreeNodeWidget.nameText.ButtonDownFcn should invoke
% graphicalTree.clickNodeButtonCallback, of the form {@foo}
%%%SU
ax = axes;
tree = graphicalTree(ax, 'larch');
cbInfo = 'callback happened';
foo = @(node, event) set(ax, 'UserData', cbInfo);
tree.clickNodeButtonCallback = {foo};
text = tree.widgetList.getValue(1).nameText;
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
text = tree.widgetList.getValue(1).nameText;
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

% should be able to destroy all nodes and widgets
ax = axes;
tree = graphicalTree(ax, 'larch');
trunk = tree.trunk;
branch = tree.newNode(trunk, 'branch');
leaf = tree.newNode(branch, 'leaf');
childrenHadParent = isobject(branch.getParent) && isobject(leaf.getParent);
parentsHadChild = trunk.numChildren && branch.numChildren;
treeHadTrunk = isobject(tree.trunk);
treeHadDrawWidgets = isobject(tree.widgetList.getValue(1));
tree.forgetInternalReferences;
forgotWidgets = tree.widgetList.length == 0;
forgotNodes = tree.nodeList.length == 0;
forgotTrunk = isempty(tree.trunk);
close all
%%%TS childrenHadParent
%%%TS parentsHadChild
%%%TS treeHadTrunk
%%%TS treeHadDrawWidgets
%%%TS forgotWidgets
%%%TS forgotNodes
%%%TS forgotTrunk

% graphicalTreeNodeWidgets should recolor their text to match the color of
% the model graphicalTreeNode
ax = axes;
tree = graphicalTree(ax, 'larch');

% red trunk
red = [1 0 0];
tree.trunk.textColor = red;
tree.draw;
trunkRedColor = get(tree.widgetList.getValue(1).nameText, 'Color');

% blue trunk
blue = [0 0 1];
tree.trunk.textColor = blue;
tree.draw;
trunkBlueColor = get(tree.widgetList.getValue(1).nameText, 'Color');

close all
%%%TS isequal(trunkRedColor, red)
%%%TS isequal(trunkBlueColor, blue)