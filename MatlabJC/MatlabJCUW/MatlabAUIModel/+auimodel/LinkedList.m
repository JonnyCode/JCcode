classdef LinkedList < handle
  properties(SetAccess=protected)
    elements = {};
    length = 0;
  end

  properties(Hidden)
    headNode;
    nodes = {};
  end

  properties(Access=protected)
    tailNode;
    isDirty;
    elemsIsDirty = false;
    nodesIsDirty = false;
  end

  methods
    function set.isDirty(self, flag)
      self.elemsIsDirty = flag;
      self.nodesIsDirty = flag;
    end
  end

  methods(Hidden)
    function pushHelper(self, value)
      import auimodel.LinkedListNode;

      node = LinkedListNode(value);
      node.nextNode = self.headNode.nextNode;
      self.headNode.nextNode = node;

      self.length = self.length + 1;
      self.isDirty = true;
    end

    function appendHelper(self, value)
      import auimodel.LinkedListNode;

      node = LinkedListNode(value);
      self.tailNode.nextNode = node;
      self.tailNode = node;

      self.length = self.length + 1;
      self.isDirty = true;
    end
  end

  methods
    function self = LinkedList()
      import auimodel.LinkedListNode;
      % dummy first node
      self.headNode = LinkedListNode();
      self.tailNode = self.headNode;
    end

    function elems = get.elements(self)
      if self.elemsIsDirty
        self.elements = self.toCell();
        self.elemsIsDirty = false;
      end

      elems = self.elements;
    end

    function elems = get.nodes(self)
      if self.nodesIsDirty
        self.nodes = self.toCell(true);
        self.nodesIsDirty = false;
      end

      elems = self.nodes;
    end

    function append(self, value, iterateCell)
      if nargin > 2 && iscell(value) && iterateCell
        i = 1;
        for i = 1 : length(value)
          self.appendHelper(value{i});
        end
      else
        self.appendHelper(value);
      end
    end

    function value = firstValue(self)
      iter = self.iterator;
      value = iter.nextValue;
    end

    function push(self, value, iterateCell)
      if nargin > 2 && iscell(value) && iterateCell
        elems = length(value);
        for i = 1 : elems
          self.pushHelper(value{elems - i + 1});
        end
      else
        self.pushHelper(value);
      end
    end

    function iter = iterator(self, mode)
      import auimodel.LinkedListIterator;

      if nargin < 2
        mode = 'links';
      end

      iter = LinkedListIterator(self, mode);
    end

    function value = pop(self)
      if self.length > 0
        value = self.headNode.nextNode.value;

        newFirst = self.headNode.nextNode.nextNode;
        self.headNode.nextNode = newFirst;

        self.length = self.length - 1;
        self.isDirty = true;
      else
        value = [];
      end
    end

    function cellArray = toCell(self, forNodes)
      cellArray = cell(self.length, 1);

      i = 1;
      node = self.headNode.nextNode;

      while ~ isempty(node)
        if nargin > 1 && forNodes
          cellArray{i} = node;
        else
          cellArray{i} = node.value;
        end

        node = node.nextNode;
        i = i + 1;
      end
    end
  end
end

