classdef LinkedListIterator < handle
  properties(SetAccess=protected)
    mode;
  end

  properties(Hidden)
    list;
    cursor;
    index;
  end

  methods
    % iterate a cached cell or by following links depending on what the user
    % wants.
    function set.mode(self, value)
      if strcmp(value, 'cell') || strcmp(value, 'links')
        self.mode = value;
      else
        error('Iterator mode must be either "cell" or "links".');
      end
    end

    function self = LinkedListIterator(list, mode)
      self.mode = mode;

      if strcmp(mode, 'cell')
        % precache cell array of elements.
        list.elements;
        self.list = list;
        self.index = 1;
      else
        self.cursor = list.headNode;
      end
    end

    function has = hasNext(self)
      % first handles case of init from empty list.
      switch self.mode
        case 'links'
          has = ~ isempty(self.cursor.nextNode);
        case 'cell'
          has = self.index < self.list.length + 1;
      end
    end

    function value = nextValue(self)
      if self.hasNext
        switch self.mode
          case 'links'
            value = self.nextNode.value;
          case 'cell'
            value = self.list.elements{self.index};
            self.index = self.index + 1;
        end
      else
        value = [];
      end
    end

    function node = nextNode(self)
      if self.hasNext
        switch self.mode
          case 'links'
            node = self.cursor.nextNode;
            self.cursor = self.cursor.nextNode;

          case 'cell'
            node = self.list.nodes{self.index};
            self.index = self.index + 1;
        end
      else
        node = [];
      end
    end
  end
end

