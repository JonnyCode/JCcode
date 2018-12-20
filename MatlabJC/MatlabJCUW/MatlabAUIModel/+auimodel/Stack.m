% a simple, cell-backed stack
classdef Stack < auimodel.LinkedList
  methods
    function value = count(self)
      value = self.length;
    end
  end
end

