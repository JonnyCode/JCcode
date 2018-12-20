classdef graphicalTreeNode < handle
    
    properties
        tree = graphicalTree.empty;
        depth = 0;
        parent = graphicalTreeNode.empty;
        children = graphicalTreeNode.empty;
        numChildren = 0;
        numDescendants = 0;
        numCheckedDescendants = 0;
        
        name = 'no name';
        textColor = [0 0 0];
        userData = [];
        
        isExpanded = false;
        isChecked = false;
        recursiveCheck = true;
        isVisibleDescandant = true;
    end
    
    properties(Hidden=true)
        userDataIndex=nan;
    end
    
    methods
        function self = graphicalTreeNode(tree, parent, name, isVisibleDescandant)
            
            if nargin < 1
                return
            end
            if nargin > 3
                self.isVisibleDescandant = isVisibleDescandant;
            end
            
            self.tree = tree;
            self.name = name;
            self.parent = parent;
            if isobject(parent)
                self.depth = parent.depth+1;
                parent.addChild(self);
            end
        end
        
        function addChild(self, child)
            if isempty(self.children)
                self.children = child;
                self.numChildren = 1;
            else
                self.numChildren = self.numChildren + 1;
                self.children(self.numChildren) = child;
            end
            
            if child.isVisibleDescandant
                self.incrementDescendants;
            end
        end
        
        function allocateChildrenArray(self, size)
            if size > length(self.children)
                self.children(size) = graphicalTreeNode;
            end
        end
        
        function incrementDescendants(self)
            self.numDescendants = self.numDescendants + 1;
            if isobject(self.parent)
                self.parent.incrementDescendants;
            end
        end
        
        function incrementCheckedDescendants(self, diff)
            self.numCheckedDescendants = self.numCheckedDescendants + diff;
            if isobject(self.parent)
                self.parent.incrementCheckedDescendants(diff);
            end
        end
        
        function setChecked(self, isChecked)
            % recursive down tree
            checkDiff = isChecked - self.isChecked;
            ncdWas = self.numCheckedDescendants;
            ncd = self.becomeChecked(isChecked);
            
            diff = self.isVisibleDescandant*checkDiff + ncd - ncdWas;
            if isobject(self.parent) && diff ~= 0
                % recursive up tree
                self.parent.incrementCheckedDescendants(diff);
            end
        end
        
        function ncd = becomeChecked(self, isChecked)
            checkChanged = xor(self.isChecked, isChecked);
            self.isChecked = isChecked;
            
            if checkChanged && isobject(self.tree)
                % call out to user
                cb = self.tree.nodeBecameCheckedCallback;
                if length(cb) == 1
                    feval(cb{1}, self);
                elseif length(cb) > 1
                    feval(cb{1}, self, cb{2:end});
                end
            end
            
            % count up checked descendants
            ncd = 0;
            if isobject(self.children)
                if self.recursiveCheck
                    % count and set
                    ncd = sum(isChecked.*[self.children.isVisibleDescandant]);
                    for ii = 1:length(self.children)
                        ncd = ncd + self.children(ii).becomeChecked(isChecked);
                    end
                    self.numCheckedDescendants = ncd;
                else
                    % only count
                    ncd = countCheckedDescendants(self);
                end
            end
        end
        
        function ncd = countCheckedDescendants(self)
            if isobject(self.children)
                ncd = sum([self.children.isChecked] & [self.children.isVisibleDescandant]);
                for ii = 1:length(self.children)
                    ncd = ncd + self.children(ii).countCheckedDescendants;
                end
            else
                ncd = 0;
            end
            self.numCheckedDescendants = ncd;
        end
        
        function includeUnburied(self)
            
            % tree should draw me
            self.tree.includeInDraw(self);
            
            % recur on children
            if self.isExpanded && isobject(self.children)
                for ii = 1:length(self.children)
                    self.children(ii).includeUnburied;
                end
            end
        end
        
        function forgetUserData(self)
            for ii = 1:length(self.children)
                self.children(ii).forgetUserData;
            end
            self.userDataIndex = nan;
        end
        
        function userData = get.userData(self)
            if isobject(self.tree) && ~isnan(self.userDataIndex)
                userData = self.tree.getNodeUserData(self.userDataIndex);
            else
                userData = [];
            end
        end
        
        function set.userData(self, userData)
            if isobject(self.tree)
                if isnan(self.userDataIndex)
                    self.userDataIndex = self.tree.appendNodeUserData(userData);
                else
                    self.tree.setNodeUserData(self.userDataIndex, userData);
                end
            end
        end
        
        function forgetTreeReferences(self)
            for ii = 1:length(self.children)
                self.children(ii).forgetTreeReferences;
            end
            self.children = [];
            self.parent = [];
        end
    end
end