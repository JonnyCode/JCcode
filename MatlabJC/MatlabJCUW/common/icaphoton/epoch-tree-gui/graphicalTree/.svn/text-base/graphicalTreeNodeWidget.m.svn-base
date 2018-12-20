classdef graphicalTreeNodeWidget < handle
    %Reusable container for graphics used by graphicalTreeNode
    
    properties
        tree;
        ax;
        boundNode;
        drawIndex;
        
        inset = 3;
        downset = 1.5;
        position;
        
        group;
        expandBox;
        checkBox;
        nameText;
    end
    
    methods
        function self = graphicalTreeNodeWidget(tree)
            if nargin < 1
                return
            end

            % ax set method does significant config            
            self.tree = tree;
            self.ax = tree.ax;
        end
        
        function buildWidgets(self)
            self.group = hggroup('Parent', self.ax);
            self.nameText = text( ...
                'BackgroundColor',  'none', ...
                'Margin',           2, ...
                'Editing',          'off', ...
                'LineStyle',        '-', ...
                'LineWidth',        1, ...
                'Interpreter',      'none', ...
                'Units',            'data', ...
                'EraseMode',        'background', ...
                'Selected',         'off', ...
                'SelectionHighlight',   'off', ...
                'VerticalAlignment',    'middle', ...
                'HorizontalAlignment',  'left', ...
                'HitTest',          'on', ...
                'ButtonDownFcn',    {@graphicalTree.clickWidgetCallback, self}, ...
                'Parent',           self.group);
            self.expandBox = graphicalCheckBox(self.group);
            self.expandBox.color = [1 0 0];
            self.expandBox.callback = {@graphicalTree.expandWidgetCallback, self};
            self.checkBox = graphicalCheckBox(self.group);
            self.checkBox.color = [.5 .5 1];
            self.checkBox.callback = {@graphicalTree.checkWidgetCallback, self};
        end
        
        function bindNode(self, drawCount)
            self.drawIndex = drawCount;
            self.setPositions(drawCount);
            set(self.nameText, ...
                'String', self.boundNode.name, ...
                'Color', self.boundNode.textColor);
            self.expandBox.isChecked = self.boundNode.isExpanded;
            self.checkBox.isChecked = self.boundNode.isChecked;
            self.checkBox.isAlternateChecked = self.partialSelection;
            set(self.group, 'Visible', 'on');
        end
        
        function unBind(self)
            %self.boundNode = [];
            self.drawIndex = [];
            set(self.group, 'Visible', 'off');
        end
        
        function setPositions(self, drawCount)
            pos = [self.boundNode.depth*self.inset, (drawCount-.5)*self.downset];
            self.position = pos;
            self.expandBox.position = pos;
            set(self.nameText, 'Position', pos+[5 0]);
            self.checkBox.position = pos+[2 0];
        end
        
        function showHighlight(self, isHighlighted)
            if isHighlighted
                set(self.nameText, 'EdgeColor', [0 1 0]);
            else
                set(self.nameText, 'EdgeColor', 'none');
            end
        end
        
        function isPartial = partialSelection(self)
            %                   isChecked   ~isChecked
            % no descendants    /           []
            % some descendants  /           []
            % all escendants    X           []
            isPartial = self.boundNode.isChecked ...
                && self.boundNode.numCheckedDescendants/self.boundNode.numDescendants < 1;
        end
        
        function set.ax(self, ax)
            self.ax = ax;
            if ~isempty(self.group) && ishandle(self.group)
                set(self.group, 'Parent', ax);
            else
                self.buildWidgets;
            end
        end
    end
end