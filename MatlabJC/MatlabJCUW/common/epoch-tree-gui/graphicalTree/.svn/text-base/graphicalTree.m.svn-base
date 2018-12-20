classdef graphicalTree < handle
    
    properties
        ax;
        name;
        trunk;

        nodeList;
        currentNodeKey;
        
        % feval({callback{1}, node, event, callback{2:end}})
        expandNodeButtonCallback;
        checkNodeButtonCallback;
        clickNodeButtonCallback;
        
        % feval({callback{1}, node, callback{2:end}})
        nodeBecameCheckedCallback;
        
        isBusy=false;
    end
    
    properties(Hidden = true)
        axParent;

        userDataList;

        drawCount = 0;
        initialWidgets = 100;
        widgetList;
        currentWidgetKey;
    end
    
    methods
        function self = graphicalTree(ax, name)
            
            self.nodeList = graphicalTreeList;
            self.userDataList = graphicalTreeList;
            self.widgetList = graphicalTreeList;
            
            if nargin > 1
                self.name = name;
            end
            self.trunk = self.newNode([], self.name);
            self.trunk.isExpanded = true;            

            if nargin > 0
                % set.ax method does significant axes config.
                self.ax = ax;
                
                % assume a bunch of widgets are needed up front
                for ii = 1:self.initialWidgets
                    widget = graphicalTreeNodeWidget(self);
                    self.widgetList.append(widget);
                    widget.selfKey = self.widgetList.length;
                end
                
                % begin with consistent appearance
                self.draw;
            end
        end
        
        function node = newNode(self, parent, name, isVisibleDescendant)
            if nargin < 2
                parent = [];
            end
            if nargin < 3
                name = '';
            end
            if nargin < 4
                isVisibleDescendant = true;
            end
            
            node = graphicalTreeNode(name, isVisibleDescendant);
            self.nodeList.append(node);
            
            node.selfKey = self.nodeList.length;
            node.tree = self;
            if isobject(parent) && isa(parent, 'graphicalTreeNode')
                node.parentKey = parent.selfKey;
                node.depth = parent.depth+1;
                parent.addChild(node);
            end
        end
        
        function index = appendNodeUserData(self, userData)
            self.userDataList.append(userData);
            index = self.userDataList.length;
        end
        
        function nodeUserData = getNodeUserData(self, index)
            nodeUserData = self.userDataList.getValue(index);
        end
        
        function setNodeUserData(self, index, nodeUserData)
            self.userDataList.setValue(index, nodeUserData);
        end
        
        function includeInDraw(self, node)
            % make sure there's a graphical widget available for caller
            self.drawCount = self.drawCount + 1;
            if self.drawCount > self.widgetList.length
                % double up on widgets
                for ii = self.drawCount:(self.widgetList.length*2)
                    widget = graphicalTreeNodeWidget(self);
                    self.widgetList.append(widget);
                    widget.selfKey = self.widgetList.length;
                end
            end
            widget = self.widgetList.getValue(self.drawCount);
            widget.boundNodeKey = node.selfKey;
        end
        
        function draw(self);
            % figure out which nodes to draw
            self.drawCount = 0;
            self.trunk.includeUnburied;
            
            % set widget appearances
            widgets = self.widgetList.getAllValues;
            for ii = 1:self.drawCount
                widgets{ii}.bindNode(ii);
            end
            
            % hide unused widgets
            for ii = self.drawCount+1:length(widgets)
                widgets{ii}.unBind;
            end
        end
        
        function redrawChecks(self)
            % quickly(?) redraw only visible widgets
            widgets = self.widgetList.getAllValues;
            for ii = 1:self.drawCount
                widget = widgets{ii};
                node = self.nodeList.getValue(widget.boundNodeKey);
                widget.checkBox.isChecked = node.isChecked;
                widget.checkBox.isAlternateChecked = widget.partialSelection;
            end
        end
        
        function setCurrentWidget(self, currentWidget)
            formerWidget = self.getCurrentWidget;
            if isobject(formerWidget)
                formerWidget.showHighlight(false);
            end
            currentWidget.showHighlight(true);
            self.currentWidgetKey = currentWidget.selfKey;
        end
        
        function currentWidget = getCurrentWidget(self)
            if ~isempty(self.currentWidgetKey)
                currentWidget = self.widgetList.getValue(self.currentWidgetKey);
            else
                currentWidget = [];
            end
        end
        
        function currentNode = getCurrentNode(self)
            widget = self.getCurrentWidget;
            if isobject(widget)
                currentNode = self.nodeList.getValue(widget.boundNodeKey);
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
            self.userDataList.removeAllValues;
        end
        
        function forgetInternalReferences(self)
            widgets = self.widgetList.getAllValues;
            for ii = 1:length(widgets)
                delete(widgets{ii});
            end
            self.widgetList.removeAllValues;
            
            nodes= self.nodeList.getAllValues;
            for ii = 1:length(nodes)
                delete(nodes{ii});
            end
            self.nodeList.removeAllValues;
            
            delete(self.trunk);
            self.trunk = [];
        end
        
        function set.ax(self, ax)
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
            if self.widgetList.length > 0
                widgets = self.widgetList.getAllValues;
                for ii = 1:length(widgets)
                    widgets{ii}.ax = ax;
                end
                self.draw;
            end
        end
        
        function set.isBusy(self, isBusy)
            self.isBusy = isBusy;
            widget = self.getCurrentWidget;
            if isobject(widget)
                widget.showBusy(isBusy);
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
                treeBottom = gTree.widgetList.getValue(gTree.drawCount).position;
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
            widget = gTree.getCurrentWidget;
            if isobject(widget)
                switch event.Key
                    case 'space'
                        check = false;
                    case 'x'
                        check = true;
                    otherwise
                        return
                end
                widget.checkBox.isChecked = check;
                graphicalTree.checkWidgetCallback(widget.checkBox, event, widget);
            end
        end
        function expandNodeWithArrowKeys(fig, event, gTree)
            widget = gTree.getCurrentWidget;
            if isobject(widget)
                switch event.Key
                    case 'leftarrow'
                        expand = false;
                    case 'rightarrow'
                        expand = true;
                    otherwise
                        return
                end
                widget.expandBox.isChecked = expand;
                graphicalTree.expandWidgetCallback(widget.expandBox, event, widget);
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

            currentWidget = gTree.getCurrentWidget;            
            if isobject(currentWidget)
                index = min(max(currentWidget.drawIndex+increment, 1), gTree.drawCount);
            else
                index = 1;
            end
            widget = gTree.widgetList.getValue(index);
            
            seekSelected = ~isempty(event.Modifier) && strcmp(event.Modifier, 'shift');
            while index>=1 && index<=gTree.drawCount && seekSelected
                widget = gTree.widgetList.getValue(index);
                seekSelected = ~widget.checkBox.isChecked;
                index = index+increment;
            end
            
            % behave just like mouseclick
            graphicalTree.clickWidgetCallback(widget.nameText, event, widget)
        end
        
        % node widget callbacks:
        function expandWidgetCallback(expandBox, event, nodeWidget)
            node = nodeWidget.tree.nodeList.getValue(nodeWidget.boundNodeKey);
            node.isExpanded = expandBox.isChecked;          

            nodeWidget.tree.setCurrentWidget(nodeWidget);
            nodeWidget.tree.isBusy = true;
            drawnow;
            
            nodeWidget.tree.draw;
            
            cb = nodeWidget.tree.expandNodeButtonCallback;
            if length(cb) == 1
                feval(cb{1}, node, event);
            elseif length(cb) > 1
                feval(cb{1}, node, event, cb{2:end});
            end
            nodeWidget.tree.isBusy = false;
            drawnow;
        end
        function checkWidgetCallback(checkBox, event, nodeWidget)
            node = nodeWidget.tree.nodeList.getValue(nodeWidget.boundNodeKey);
            node.setChecked(checkBox.isChecked);            

            nodeWidget.tree.setCurrentWidget(nodeWidget);
            nodeWidget.tree.isBusy = true;
            drawnow;
            
            nodeWidget.tree.redrawChecks;
            
            cb = nodeWidget.tree.checkNodeButtonCallback;
            if length(cb) == 1
                feval(cb{1}, node, event);
            elseif length(cb) > 1
                feval(cb{1}, node, event, cb{2:end});
            end
            nodeWidget.tree.isBusy = false;
            drawnow;
        end
        function clickWidgetCallback(nameText, event, nodeWidget)
            node = nodeWidget.tree.nodeList.getValue(nodeWidget.boundNodeKey);
            
            nodeWidget.tree.setCurrentWidget(nodeWidget);
            drawnow;
            
            cb = nodeWidget.tree.clickNodeButtonCallback;
            if length(cb) == 1
                feval(cb{1}, node, event);
            elseif length(cb) > 1
                feval(cb{1}, node, event, cb{2:end});
            end
        end
    end
end