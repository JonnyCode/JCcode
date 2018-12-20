classdef graphicalTree < handle
    
    properties
        ax;
        fig;
        trunk;
        currentNode;

        % feval({callback{1}, node, event, callback{2:end}})
        expandNodeButtonCallback;
        checkNodeButtonCallback;
        clickNodeButtonCallback;

        % feval({callback{1}, node, callback{2:end}})
        nodeIsCheckedCallback;
    end
    
    methods
        function self = graphicalTree(ax, name)
            self.ax = ax;
            self.fig = get(self.ax, 'Parent');
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
                'Visible',  'on');
            self.trunk = graphicalTreeNode(self, [], name);
            self.trunk.isExpanded = true;
            
            set(graphicalTree.resolveParentFigure(self.fig), ...
                'Renderer', 'painters', ...
                'WindowButtonDownFcn', {@graphicalTree.panAxesOnAltClick, self.ax});
            set(self.fig, ...
                'ResizeFcn', {@graphicalTree.adjustAxesOnFigureResize, self.ax});
            self.draw;
        end
        
        function node = newNode(self, parentNode, name)
            node = graphicalTreeNode(self, parentNode, name);
        end
        
        function self = set.currentNode(self, currentNode)
            if isobject(self.currentNode)
                self.currentNode.showHighlight(false);
            end
            currentNode.showHighlight(true);
            self.currentNode = currentNode;
        end
        
        function draw(self);
            self.trunk.showChildren;
        end
    end
    
    methods(Static)
        
        function fig = resolveParentFigure(fig)
            % fig might be figure, uibuttongroup, uipanel
            while ~strcmp(get(fig, 'Type'), 'figure')
                fig = get(fig, 'Parent');
            end
        end
        
        function [xl, yl] = getAxesLimsFromAxesSize(ax)
            oldUnits = get(ax, 'Units');
            set(ax, 'Units', 'characters');
            pos = get(ax, 'Position');
            xl = [0, max(eps, pos(3))];
            yl = [0, max(eps, pos(4))];
            set(ax, 'Units', oldUnits);
        end
        
        function adjustAxesOnFigureResize(fig, event, ax)
            [xl, yl] = graphicalTree.getAxesLimsFromAxesSize(ax);
            set(ax, 'XLim', xl, 'YLim', yl);
        end
        
        function panAxesOnAltClick(fig, event, ax)
            if strcmp(get(fig, 'SelectionType'), 'extend')
                pan(fig);
            end
        end
        
        % node widget callbacks:
        function expandWidgetCallback(widget, event, node)
            node.isExpanded = widget.isChecked;
            node.tree.currentNode = node;
            node.tree.draw;
            drawnow;
            
            cb = node.tree.expandNodeButtonCallback;
            if length(cb) == 1
                feval(cb{1}, node, event);
            elseif length(cb) > 1
                feval(cb{1}, node, event, cb{2:end});
            end
        end
        function checkWidgetCallback(widget, event, node)
            node.setCheckedRecursive(widget.isChecked);
            node.tree.currentNode = node;
            node.tree.draw;
            drawnow;

            cb = node.tree.checkNodeButtonCallback;
            if length(cb) == 1
                feval(cb{1}, node, event);
            elseif length(cb) > 1
                feval(cb{1}, node, event, cb{2:end});
            end
        end
        function clickWidgetCallback(widget, event, node)
            node.tree.currentNode = node;
            drawnow;
            
            cb = node.tree.clickNodeButtonCallback;
            if length(cb) == 1
                feval(cb{1}, node, event);
            elseif length(cb) > 1
                feval(cb{1}, node, event, cb{2:end});
            end
        end
    end
end