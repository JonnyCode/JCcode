classdef graphicalTree < handle
    
    properties
        ax;
        name;
        trunk;
        currentNode;
        
        % feval({callback{1}, node, event, callback{2:end}})
        expandNodeButtonCallback;
        checkNodeButtonCallback;
        clickNodeButtonCallback;
        
        % feval({callback{1}, node, callback{2:end}})
        nodeBecameCheckedCallback;
    end
    
    properties(Hidden=true)
        axParent;
        initialWidgets=25;
        drawCount=0;
        drawWidgets;
        currentWidget;
        
        nodeUserData;
    end
    
    methods
        function self = graphicalTree(ax, name)
            if nargin < 1
                return
            end
            
            self.name = name;
            self.trunk = self.newNode([], name);
            self.trunk.isExpanded = true;

            % ax set method does significant axes config.
            self.ax = ax;
            
            % assume a bunch of widgets are needed up front
            self.drawWidgets = graphicalTreeNodeWidget(self);
            for ii = self.initialWidgets:-1:2
                self.drawWidgets(ii) = graphicalTreeNodeWidget(self);
            end
            
            % begin with consistent appearance
            self.draw;
        end
        
        function node = newNode(self, parentNode, name)
            node = graphicalTreeNode(self, parentNode, name);
        end
        
        function index = appendNodeUserData(self, userData)
            if isempty(self.nodeUserData)
                index = 1;
                self.nodeUserData = containers.Map(index, userData, 'uniformvalues', true);
            else
                index = double(self.nodeUserData.Count)+1;
                self.nodeUserData(index) = userData;
            end
        end
        
        function nodeUserData = getNodeUserData(self, index)
            nodeUserData = self.nodeUserData(index);
        end
        
        function setNodeUserData(self, index, nodeUserData)
            self.nodeUserData(index) = nodeUserData;
        end
        
        function includeInDraw(self, node)
            % make/reuse graphical widget for this node
            self.drawCount = self.drawCount + 1;
            if self.drawCount > length(self.drawWidgets)
                % double up on widgets
                for ii = self.drawCount:(length(self.drawWidgets)*2)
                    self.drawWidgets(ii) = graphicalTreeNodeWidget(self);
                end
            end
            self.drawWidgets(self.drawCount).boundNode = node;
        end
        
        function draw(self);
            % figure out which nodes to draw
            self.drawCount = 0;
            self.trunk.includeUnburied;
            
            % set widget appearances
            for ii = 1:self.drawCount
                self.drawWidgets(ii).bindNode(ii);
            end
            
            % hide unused widgets
            for ii = self.drawCount+1:length(self.drawWidgets)
                self.drawWidgets(ii).unBind;
            end
        end
        
        function redrawChecks(self)
            % quickly(?) redraw only visible widgets
            for ii = 1:self.drawCount
                self.drawWidgets(ii).checkBox.isChecked = ...
                    self.drawWidgets(ii).boundNode.isChecked;
                self.drawWidgets(ii).checkBox.isAlternateChecked = ...
                    self.drawWidgets(ii).partialSelection;
            end
        end
        
        function setCurrentWidget(self, currentWidget)
            if isobject(self.currentWidget)
                self.currentWidget.showHighlight(false);
            end
            currentWidget.showHighlight(true);
            self.currentWidget = currentWidget;
        end
        
        function currentNode = get.currentNode(self)
            if isobject(self.currentWidget)
                currentNode = self.currentWidget.boundNode;
            else
                currentNode = [];
            end
        end
        
        function fig = resolveParentFigure(self)
            % get parent figure of self.axes
            fig = self.ax;
            while fig > 0 && ~strcmp(get(fig, 'Type'), 'figure')
                fig = get(fig, 'Parent');
            end
        end
        
        function forgetExternalReferences(self)
            self.expandNodeButtonCallback = [];
            self.checkNodeButtonCallback = [];
            self.clickNodeButtonCallback = [];
            self.nodeBecameCheckedCallback = [];
            self.trunk.forgetUserData;
            self.nodeUserData=[];
        end
        
        function forgetInternalReferences(self)
            self.trunk.forgetTreeReferences;
            self.trunk = [];
            for ii = 1:length(self.drawWidgets)
                self.drawWidgets(ii).boundNode = [];
            end
            self.drawWidgets = [];
        end
        
        function self = set.ax(self, ax)
            self.ax = ax;
            
            % make consistent
            self.axParent = get(self.ax, 'Parent');
            [xl, yl] = graphicalTree.getAxesLimsFromAxesSize(self.ax);
            set(self.ax, ...
                'Box',      'off', ...
                'Color',    [1 1 1], ...
                'Units',    'normalized', ...
                'XLim',     xl, ...
                'XScale',   'linear', ...
                'XTick',    [], ...
                'YLim',     yl, ...
                'YScale',   'linear', ...
                'YTick',    [], ...
                'YDir',     'reverse', ...
                'HitTest',  'on', ...
                'DrawMode', 'fast', ...
                'Visible',  'on');
            
            % parent figure and axParent might be the same object,
            %   or axParent might be e.g. a uipanel
            set(self.resolveParentFigure, ...
                'Renderer', 'painters', ...
                'Units',    'pixels', ...
                'WindowScrollWheelFcn',{@graphicalTree.scrollAxesWithMouseWheel, self}, ...
                'KeyPressFcn',  {@graphicalTree.navigateWithArrowKeys, self});
            set(self.axParent, ...
                'ResizeFcn', {@graphicalTree.adjustAxesOnContainerResize, self.ax});
            
            % inform graphical children
            if isobject(self.drawWidgets)
                for ii = 1:length(self.drawWidgets)
                    self.drawWidgets(ii).ax = ax;
                end
                self.draw;
            end
        end
    end
    
    methods(Static)
        
        % graphical containers callbacks
        function [xl, yl] = getAxesLimsFromAxesSize(ax)
            oldUnits = get(ax, 'Units');
            set(ax, 'Units', 'characters');
            pos = get(ax, 'Position');
            xl = [0, max(eps, pos(3))];
            yl = [0, max(eps, pos(4))];
            set(ax, 'Units', oldUnits);
        end
        function adjustAxesOnContainerResize(container, event, ax)
            [xl, yl] = graphicalTree.getAxesLimsFromAxesSize(ax);
            set(ax, 'XLim', xl, 'YLim', yl);
        end
        function scrollAxesWithMouseWheel(fig, event, gTree)
            % scroll when mouse is over the axes container
            mousePoint = get(fig, 'CurrentPoint');
            containerPos = getpixelposition(gTree.axParent);
            if mousePoint(1) >= containerPos(1) ...
                    && mousePoint(1) <= containerPos(1) + containerPos(3) ...
                    && mousePoint(2) >= containerPos(2) ...
                    && mousePoint(2) <= containerPos(2) + containerPos(4) ...
                    
                % keep axes on the tree
                treeBottom = gTree.drawWidgets(gTree.drawCount).position;
                yl = get(gTree.ax, 'YLim');
                inc = event.VerticalScrollCount;
                axTop = min(max(yl(1)+inc, 0), treeBottom(2));
                set(gTree.ax, 'YLim', [axTop, axTop+yl(2)-yl(1)]);
            end
        end
        function navigateWithArrowKeys(fig, event, gTree)
            switch event.Key
                case {'uparrow', 'downarrow'}
                    graphicalTree.highlightNodeWithArrowKeys(fig, event, gTree);
                case {'leftarrow', 'rightarrow'}
                    graphicalTree.expandNodeWithArrowKeys(fig, event, gTree);
                case {'space', 'x'}
                    graphicalTree.checkNodeWithSpaceAndX(fig, event, gTree);
            end
        end
        function checkNodeWithSpaceAndX(fig, event, gTree)
            if isobject(gTree.currentWidget)
                switch event.Key
                    case 'space'
                        check = false;
                    case 'x'
                        check = true;
                    otherwise
                        return
                end
                gTree.currentWidget.checkBox.isChecked = check;
                graphicalTree.checkWidgetCallback(gTree.currentWidget.checkBox, event, gTree.currentWidget);
            end
        end        
        function expandNodeWithArrowKeys(fig, event, gTree)
            if isobject(gTree.currentWidget)
                switch event.Key
                    case 'leftarrow'
                        expand = false;
                    case 'rightarrow'
                        expand = true;
                    otherwise
                        return
                end
                gTree.currentWidget.expandBox.isChecked = expand;
                graphicalTree.expandWidgetCallback(gTree.currentWidget.expandBox, event, gTree.currentWidget);
            end
        end
        function highlightNodeWithArrowKeys(fig, event, gTree)
            % highlight visible nodes with up/down arrow keys
            switch event.Key
                case 'uparrow'
                    increment = -1;
                case 'downarrow'
                    increment = +1;
                otherwise
                    return
            end
            
            if isobject(gTree.currentWidget)
                index = min(max(gTree.currentWidget.drawIndex+increment, 1), gTree.drawCount);
            else
                index = 1;
            end
            widget = gTree.drawWidgets(index);
            
            seekSelected = ~isempty(event.Modifier) && strcmp(event.Modifier, 'shift');
            while index>=1 && index<=gTree.drawCount && seekSelected
                widget = gTree.drawWidgets(index);
                seekSelected = ~widget.checkBox.isChecked;
                index = index+increment;
            end
            
            % behave just like mouseclick
            graphicalTree.clickWidgetCallback(widget.nameText, event, widget)
        end
        
        % node widget callbacks:
        function expandWidgetCallback(expandBox, event, nodeWidget)
            nodeWidget.boundNode.isExpanded = expandBox.isChecked;
            nodeWidget.tree.setCurrentWidget(nodeWidget);
            nodeWidget.tree.draw;
            drawnow;
            
            cb = nodeWidget.tree.expandNodeButtonCallback;
            if length(cb) == 1
                feval(cb{1}, nodeWidget.boundNode, event);
            elseif length(cb) > 1
                feval(cb{1}, nodeWidget.boundNode, event, cb{2:end});
            end
        end
        function checkWidgetCallback(checkBox, event, nodeWidget)
            nodeWidget.boundNode.setChecked(checkBox.isChecked);
            nodeWidget.tree.setCurrentWidget(nodeWidget);
            nodeWidget.tree.redrawChecks;
            drawnow;
            
            cb = nodeWidget.tree.checkNodeButtonCallback;
            if length(cb) == 1
                feval(cb{1}, nodeWidget.boundNode, event);
            elseif length(cb) > 1
                feval(cb{1}, nodeWidget.boundNode, event, cb{2:end});
            end
        end
        function clickWidgetCallback(nameText, event, nodeWidget)
            nodeWidget.tree.setCurrentWidget(nodeWidget);
            drawnow;
            
            cb = nodeWidget.tree.clickNodeButtonCallback;
            if length(cb) == 1
                feval(cb{1}, nodeWidget.boundNode, event);
            elseif length(cb) > 1
                feval(cb{1}, nodeWidget.boundNode, event, cb{2:end});
            end
        end
    end
end