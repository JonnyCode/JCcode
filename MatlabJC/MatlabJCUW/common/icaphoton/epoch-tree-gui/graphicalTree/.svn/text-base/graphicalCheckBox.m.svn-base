classdef graphicalCheckBox < handle
    
    properties
        parent;
        position;
        color;
        
        % feval({callback{1}, node, event, callback{2:end}})
        callback={};
        isChecked;
        isAlternateChecked;
    end
    
    properties(Hidden=true)
        box;
        checkedSymbol = 'X';
        uncheckedSymbol = '  ';
        altCheckedSymbol = '/';
    end
    
    methods
        function self = graphicalCheckBox(parent)
            if nargin < 1
                return
            end
            self.parent = parent;
            self.box = text( ...
                'Color',            [0 0 0], ...
                'EdgeColor',        'none', ...
                'Margin',           2, ...
                'Editing',          'off', ...
                'Interpreter',      'none', ...
                'EraseMode',        'background', ...
                'Units',            'data', ...
                'Selected',         'off', ...
                'SelectionHighlight',   'off', ...
                'HorizontalAlignment',  'left', ...
                'VerticalAlignment',    'middle', ...
                'HitTest',          'on', ...
                'ButtonDownFcn',    {@graphicalCheckBox.boxButtonDownFcn, self}, ...
                'Parent',           self.parent);
            
            self.position = [0 .5];
            self.color = [1 0 0];
            self.isChecked = false;
            self.isAlternateChecked = false;
        end
        
        % property access methods:
        function set.parent(self, parent)
            set(self.box, 'Parent', parent);
            self.parent = parent;
        end
        function set.position(self, position)
            set(self.box, 'Position', position);
            self.position = position;
        end
        function set.color(self, color)
            set(self.box, 'BackgroundColor', color);
            self.color = color;
        end
        function set.isChecked(self, isChecked)
            if isChecked
                set(self.box, 'String', self.checkedSymbol);
            else
                set(self.box, 'String', self.uncheckedSymbol);
            end
            self.isChecked = isChecked;
        end
        function set.isAlternateChecked(self, isAlternateChecked)
            if isAlternateChecked
                set(self.box, 'String', self.altCheckedSymbol);
            end
            self.isAlternateChecked = isAlternateChecked;
        end
    end
    
    methods(Static)
        function boxButtonDownFcn(box, event, gCheckBox)
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