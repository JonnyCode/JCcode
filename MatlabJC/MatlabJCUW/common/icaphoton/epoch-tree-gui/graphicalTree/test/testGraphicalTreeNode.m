function testGraphicalTreeNode
%Run graphicalTreeNode through some paces, without any graphicalTree

% should be constructable with no arguments
%%%SU
node = graphicalTreeNode();
%%%TS isobject(node)
%%%TS isa(node, 'graphicalTreeNode')


% should be constructable with parent node and name arguments
%   nodes should refer to parent and child
%%%SU
parent = graphicalTreeNode([], [], 'parent');
child = graphicalTreeNode([], parent, 'child');
%%%TS isobject(parent) && isobject(child)
%%%TS isa(parent, 'graphicalTreeNode') && isa(child, 'graphicalTreeNode')
%%%TS strcmp(parent.name, 'parent') && strcmp(child.name, 'child')
%%%TS isequal(parent.children, child)
%%%TS isequal(child.parent, parent)


% parent nodes should account for descendants
% nodes should know how deep they are
%%%SU
parent = graphicalTreeNode([], [], 'parent');
child1 = graphicalTreeNode([], parent, 'child');
child2 = graphicalTreeNode([], parent, 'child');
grandchild = graphicalTreeNode([], child1, 'grandchild');
%%%TS isequal(length(parent.children), 2)
%%%TS isequal(parent.numDescendants, 3)
%%%TS isequal(child1.numDescendants, 1)
%%%TS isequal(parent.depth, 0)
%%%TS isequal(child1.depth, 1)
%%%TS isequal(child2.depth, 1)
%%%TS isequal(grandchild.depth, 2)

% nodes should recursively check descendants
%%%SU
parent = graphicalTreeNode([], [], 'parent');
parent.isChecked = false;
child = graphicalTreeNode([], parent, 'child');
child.isChecked = false;
grandchild = graphicalTreeNode([], child, 'grandchild');
grandchild.isChecked = false;
parent.recursiveCheck = true;
parent.setChecked(true);
%%%TS parent.isChecked
%%%TS child.isChecked
%%%TS grandchild.isChecked

% nodes should recursively uncheck descendants
%%%SU
parent = graphicalTreeNode([], [], 'parent');
parent.isChecked = true;
child = graphicalTreeNode([], parent, 'child');
child.isChecked = true;
grandchild = graphicalTreeNode([], child, 'grandchild');
grandchild.isChecked = true;
parent.recursiveCheck = true;
parent.setChecked(false);
%%%TS ~parent.isChecked
%%%TS ~child.isChecked
%%%TS ~grandchild.isChecked

% nodes should recursively account for checked descendants
%%%SU
parent = graphicalTreeNode([], [], 'parent');
child = graphicalTreeNode([], parent, 'child');
grandchild = graphicalTreeNode([], child, 'grandchild');
parent.setChecked(false);
noChecks = parent.numCheckedDescendants;
grandchild.setChecked(true);
oneCheck = parent.numCheckedDescendants;
parent.setChecked(true);
twoChecks = parent.numCheckedDescendants;
child.setChecked(false);
noChecks2 = parent.numCheckedDescendants;
%%%TS isequal(noChecks, 0);
%%%TS isequal(oneCheck, 1);
%%%TS isequal(twoChecks, 2);
%%%TS isequal(noChecks2, 0);

% recursive accounting should be able to ignore nodes with
% isVisibleDescandant = false;
%%%SU
parent = graphicalTreeNode([], [], 'parent');
child = graphicalTreeNode([], parent, 'child');
grandchild = graphicalTreeNode([], child, 'grandchild', false);
oneVisibleDescandant = parent.numDescendants;

parent.setChecked(false);
noChecks = parent.numCheckedDescendants;

grandchild.setChecked(true);
noChecks2 = parent.numCheckedDescendants;

parent.setChecked(true);
oneCheck = parent.numCheckedDescendants;

child.setChecked(false);
noChecks3 = parent.numCheckedDescendants;
%%%TS oneVisibleDescandant == 1
%%%TS noChecks == 0;
%%%TS noChecks2 == 0;
%%%TS oneCheck == 1;
%%%TS noChecks3 == 0;

% nodes should be able to recursively forget tree references
%%%SU
parent = graphicalTreeNode([], [], 'parent');
child = graphicalTreeNode([], parent, 'child');
grandchild = graphicalTreeNode([], child, 'grandchild');
childrenHadParent = isobject(child.parent) && isobject(grandchild.parent);
parentsHadChild = isobject(parent.children) && isobject(child.children);
parent.forgetTreeReferences;
allForgotChildren = isempty(parent.children) && isempty(child.children) && isempty(grandchild.children);
allForgotParents = isempty(parent.parent) && isempty(child.parent) && isempty(grandchild.parent);
%%%TS childrenHadParent
%%%TS parentsHadChild
%%%TS allForgotChildren
%%%TS allForgotParents