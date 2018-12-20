classdef LinkedListNode < handle
  properties
    value = [];
    nextNode = [];
  end

  methods
    function self = LinkedListNode(value)
      if nargin > 0
        self.value = value;
      end
    end
  end
end

