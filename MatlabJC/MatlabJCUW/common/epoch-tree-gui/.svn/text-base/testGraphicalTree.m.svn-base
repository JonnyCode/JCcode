% run my graphicalTree and graphicalTreeNode classes through some paces
%   I'm not sure how to test GUI features programatically, so try the
%   following:
%
%   should display a tree in an axes on the left side of a figure
%   tree should have a tunk with 10 branch nodes and 30 total leaf nodes
%   clicking a red square should show a node's children.
%   clicking a red circle should hide a node's children.
%   clicking a blue square should "select" a node and all its children
%   clicking a blue circle should de"select" a node and all its children
%   the last node clicked on should be highlighted in green
%   nodes should remain "select"ed and expanded even when buried under a
%       collapsed parent node
%   the pan(hand) tool in the figure tool bar should let you view different
%       parts of the tree (especially when its too big for the panel).

% benjamin.heasly@gmail.com
%   26 Jan 2009

close all
clear all
clear classes

f = figure;
pan = uipanel( ...
    'Parent',   f, ...
    'Units',    'normalized', ...
    'Position', [0 0 .5 1]);
ax = axes( ...
    'Parent',   pan, ...
    'Units',    'normalized', ...
    'Position', [0 0 1 .9]);

gTree = graphicalTree(ax, 'testTree');

for ii = 1:10
    branches(ii) = gTree.newNode(gTree.trunk, sprintf('branch %d', ii));
    branches(ii).isExpanded = false;
end
for ii = 1:30
    leaves(ii) = gTree.newNode(branches(ceil(ii/3)), sprintf('leaf %d', ii));
end

gTree.draw;