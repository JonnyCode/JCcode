% primitive, list-backed set.  add, contains, remove all O(n).
classdef Set < handle
  properties(SetAccess=private)
    elements;
  end

  properties(Access=private)
    map;
  end

  methods
    function self = Set()
      import containers.Map;
      self.map = Map();
    end

    function elems = get.elements(self)
      elems = self.map.keys;
    end

    function add(self, element)
      if ~ self.map.isKey(element)
        self.map(element) = [];
      end
    end

    function remove(self, element)
      self.map.remove(element);
    end

    function result = contains(self, element)
      result = self.map.isKey(element);
    end

    function display(self)
      display(sprintf('\nelements =\n'));
      display(self.map.keys);
    end
  end
end

