classdef TreeDotVisitor < handle
  properties
    dotLines;
    dotNodes;
  end

  methods
    function self = TreeDotVisitor()
      import auimodel.LinkedList;
      import auimodel.Map;

      self.dotNodes = Map();
      self.dotLines = LinkedList();

      self.dotLines.append('digraph G {');
    end

    function finalize(self)
      self.dotLines.append('}');
    end

    function visit(self, node)
      if isobject(node.parent)
        self.writeDotConnection(node.parent, node);
      end
    end

    function declareDotNode(self, node)
      nextNodeNum = self.dotNodes.count();

      self.dotLines.append( sprintf('n%d [label = "%s"];', nextNodeNum, ...
        node.dotLabel) );

      self.dotNodes.set(node, nextNodeNum);
    end

    function writeDotConnection(self, src, dest)
      for node = [src, dest]
        if ~ self.dotNodes.hasKey(node)
          self.declareDotNode(node);
        end
      end

      self.dotLines.append(sprintf('n%d->n%d [label = "%s"];', ...
        self.dotNodes.get(src), self.dotNodes.get(dest), dest.dotEdgeLabel));
    end
  end
end

