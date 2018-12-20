% CODE DEBT: copied from MatlabAUIModel

% a simple, cell-backed list
classdef List < handle
  properties
    elements;
  end

  methods
    function self = List(startElements)
      if nargin > 0
        self.elements = startElements;
      else
        self.elements = {};
      end
    end

    function append(self, value)
      self.elements{ end + 1 } = value;
    end

    function result = contains(self, value)
      for i = 1 : length(self.elements)
        % CODE DEBT: factor out into a generalized equality method in Util
        if ischar(value)
          if strcmpi(value, self.elements{i})
            result = true;
            return;
          end
        elseif self.elements{i} == value
          result = true;
          return;
        end
      end

      result = false;
    end
    
    function remove(self, value)
      index = self.indexOf(value);

      if index > 0
        self.elements ...
          = {self.elements{1 : index - 1} self.elements{index + 1 : end}};
      end
    end

    function result = indexOf(self, value)
      for i = 1 : length(self.elements)
        if ischar(value)
          if strcmpi(value, self.elements{i})
            result = i;
            return;
          end
        elseif self.elements{i} == value
          result = i;
          return;
        end
      end

      result = -1;
    end

    function elemCount = count(self)
      elemCount = length(self.elements);
    end
  end
end

