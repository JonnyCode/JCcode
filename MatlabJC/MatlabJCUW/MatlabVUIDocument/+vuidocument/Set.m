% CODE DEBT: needs to be factored into a general Matlab utility package.
% copied from MatlabAUIModel.

% primitive, list-backed set.  add, contains, remove all O(n).
classdef Set < handle
  properties
    list;
  end

  methods
    function self = Set()
      import vuidocument.List;
      self.list = List();
    end

    function add(self, element)
      if ~ self.contains(element)
        self.list.append(element);
      end
    end

    function remove(self, element)
      self.list.remove(element);
    end

    function result = contains(self, element)
      result = self.list.contains(element);
    end

    function display(self)
      display(sprintf('\nset =\n'));
      display(self.list.elements);
    end
  end
end

