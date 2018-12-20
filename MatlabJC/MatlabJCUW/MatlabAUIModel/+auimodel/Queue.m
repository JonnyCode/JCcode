classdef Queue < auimodel.LinkedList
  methods
    function value = dequeue(self)
      value = self.pop;
    end

    function enqueue(self, value, iterateCell)
      self.append(value, iterateCell);
    end

    function value = count(self)
      value = self.length;
    end
  end
end

