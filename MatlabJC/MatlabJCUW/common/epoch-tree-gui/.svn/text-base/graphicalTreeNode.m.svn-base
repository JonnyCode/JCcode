classdef graphicalTreeNode < handle
    
    properties
        tree;
        parent;
        children;
        numDescendants = 0;
        numCheckedDescendants = 0;
        
        name;
        userData;

        isBuried = false;
        isExpanded = false;
        isChecked = false;
        recursiveCheck = true;
        
        width = 20;
        height = 1;
        inset = 3;
        upset = .5;
        position;
        familyPosition;
        
        group;
        boxWidget;
        expandWidget;
        checkWidget;
        nameWidget;
    end
    
    methods
        function self = graphicalTreeNode(tree, parent, name)
            self.tree = tree;
            self.parent = parent;
            if isobject(parent)
                parent.appendChild(self);
            end
            self.name = name;

            self.group = hggroup('Parent', tree.ax);
            self.boxWidget = rectangle( ...
                'FaceColor',    'none', ...
                'LineStyle',    '-', ...
                'HitTest',      'on', ...
                'ButtonDownFcn', {@graphicalTree.clickWidgetCallback, self}, ...
                'Parent',       self.group);
            self.nameWidget = text( ...
                'BackgroundColor', 'none', ...
                'Color',    [0 0 0], ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'left', ...
                'String',   name, ...
                'Units',    'data', ...
                'VerticalAlignment', 'cap', ...
                'HitTest',  'off', ...
                'Parent',   self.group);
            self.expandWidget = graphicalCheckBox(self.group);
            self.expandWidget.faceColor = [1 0 0];
            self.expandWidget.callback = {@graphicalTree.expandWidgetCallback, self};
            self.checkWidget = graphicalCheckBox(self.group);
            self.checkWidget.faceColor = [0 0 1];
            self.checkWidget.callback = {@graphicalTree.checkWidgetCallback, self};
            
            % trigger mutators to get correct appearance
            self.position = [0, 0, self.width, self.height];
            self.isExpanded = self.isExpanded;
            self.isChecked = self.isChecked;
            self.showHighlight(false);
        end
        
        function delete(self)
            % try to speed up closing of parent figure
            self.userData = [];
        end
        
        function appendChild(self, child)
            self.children = cat(2, self.children, child);
            self.incrementDescendants;
        end
        
        function incrementDescendants(self)
            self.numDescendants = self.numDescendants + 1;
            if isobject(self.parent)
                self.parent.incrementDescendants;
            end
        end

        function setCheckedRecursive(self, isChecked)
            self.isChecked = isChecked;
            
            cb = self.tree.nodeIsCheckedCallback;
            if length(cb) == 1
                feval(cb{1}, self);
            elseif length(cb) > 1
                feval(cb{1}, self, cb{2:end});
            end

            if self.recursiveCheck
                for child = self.children
                    child.recursiveCheck = true;
                    child.setCheckedRecursive(isChecked);
                end
            end
        end
        
        function showChildren(self)
            % recursively manage real estate
            self.familyPosition = self.position;
            if isempty(self.children)
                % base case: do nothing?
            elseif self.isBuried
                for child = self.children
                    set(child.group, 'Visible', 'off');
                    child.isBuried = true;
                    child.showChildren;
                end
            elseif self.isExpanded
                for child = self.children
                    set(child.group, 'Visible', 'on');
                    child.isBuried = false;
                    child.position = self.getNextChildPosition(child);
                    child.showChildren;
                    self.increaseFamilyPosition(child);
                end
            else
                for child = self.children
                    set(child.group, 'Visible', 'off');
                    child.isBuried = true;
                    child.showChildren;
                end
            end
        end
        
        function nextChildPosition = getNextChildPosition(self, child)
            nextChildPosition = [ ...
                self.familyPosition(1)+self.inset, ...
                self.familyPosition(2)+self.familyPosition(4)+self.upset, ...
                child.width, ...
                child.height];
        end
        
        function increaseFamilyPosition(self, child)
            % familyPosition += child.familyPosition
            x = min(self.familyPosition(1), child.familyPosition(1));
            y = min(self.familyPosition(2), child.familyPosition(2));
            left = max(self.familyPosition(1)+self.familyPosition(3), child.familyPosition(1)+child.familyPosition(3));
            top = max(self.familyPosition(2)+self.familyPosition(4), child.familyPosition(2)+child.familyPosition(4));
            self.familyPosition = [x,y,left-x,top-y];
        end
        
        function showHighlight(self, isHighlighted)
            if isHighlighted
                set(self.boxWidget, ...
                    'EdgeColor', [0 1 0], 'LineWidth', 2);
            else
                set(self.boxWidget, ...
                    'EdgeColor', [1 1 1]*.5, 'LineWidth', 1);
            end
        end
        
        % property access methods:
        function self = set.position(self, position)
            set(self.boxWidget, 'Position', position);
            set(self.nameWidget, 'Position', position(1:2)+[3 0]);
            self.expandWidget.position = position.*[1 1 0 1]+[0 0 2 0];
            self.checkWidget.position = position.*[1 1 0 1]+[self.width-2 0 2 0];
            self.position = position;
        end
        function self = set.width(self, width)
            self.position(3) = width;
            self.width = width;
        end
        function self = set.height(self, height)
            self.position(4) = height;
            self.height = height;
        end
        function self = set.name(self, name)
            set(self.nameWidget, 'String', name);
            self.name = name;
        end
        function self = set.isExpanded(self, isExpanded)
            self.expandWidget.isChecked = isExpanded;
            self.isExpanded = isExpanded;
        end
        function self = set.isChecked(self, isChecked)
            self.checkWidget.isChecked = isChecked;
            self.isChecked = isChecked;
        end
    end
end