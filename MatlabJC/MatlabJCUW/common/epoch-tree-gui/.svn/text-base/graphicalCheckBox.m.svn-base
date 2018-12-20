classdef graphicalCheckBox < handle
    
    properties
        rect;
        parent;
        position;
        faceColor;
        
        % feval({callback{1}, node, event, callback{2:end}})
        callback;
        isChecked;
    end
    
    methods
        function self = graphicalCheckBox(parent)
            self.parent = parent;
            self.rect = rectangle( ...
                'LineStyle',    '-', ...
                'LineWidth',    1, ...
                'HitTest',      'on', ...
                'EdgeColor',    [0 0 0], ...
                'ButtonDownFcn',{@graphicalCheckBox.rectButtonDownFcn, self}, ...
                'Parent',       self.parent);
            self.position = [0 0 1 1];
            self.faceColor = [1 0 0];
            self.isChecked = false;
        end
        
        function delete(self)
            % try to speed up closing of parent figure
           self.callback = []; 
        end
        
        % property access courtesy methods:
        function self = set.parent(self, parent)
            set(self.rect, 'Parent', parent);
            self.parent = parent;
        end
        function self = set.position(self, position)
            set(self.rect, 'Position', position);
            self.position = position;
        end
        function self = set.faceColor(self, faceColor)
            set(self.rect, 'FaceColor', faceColor);
            self.faceColor = faceColor;
        end
        function self = set.isChecked(self, isChecked)
            set(self.rect, 'Curvature', [1 1]*isChecked);
            self.isChecked = isChecked;
        end
    end

    methods(Static)
        function rectButtonDownFcn(rect, event, gCheckBox)
            gCheckBox.isChecked = ~gCheckBox.isChecked;
            drawnow;

            cb = gCheckBox.callback;
            if length(cb) == 1
                feval(cb{1}, gCheckBox, event);
            elseif length(cb) > 1
                feval(cb{1}, gCheckBox, event, cb{2:end});
            end
        end
    end
end