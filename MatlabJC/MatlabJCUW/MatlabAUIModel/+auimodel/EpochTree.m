classdef EpochTree < handle
  properties
    splitKey;
    splitValue;
    splitValues;

    parent;
    children;

    isLeaf;

    epochList;
    leafNodes; % contains all direct or indirect descendant leaf nodes.

    custom = struct();
  end

  properties(Hidden)
    dotLabel;
    dotEdgeLabel;
    ovationExport;
  end

  methods(Static)
    function tree = build(source, splitKeyPaths)
      import auimodel.*;

      if ~ iscellstr(splitKeyPaths)
        MException('AUIModel:BadParam', ...
          'Second parameter must be a string cell of key paths.').throw();
      end

      if isa(source, 'auimodel.EpochList')
        epochList = source;
      elseif ExportLoader.isValidExport(source)
        epochList = EpochList(source);
      else
        error('Must supply an EpochList or OvationExport as first param.');
      end

      pathsStack = Stack();
      pathsStack.push(splitKeyPaths, true);

      tree = EpochTree(epochList, pathsStack, Stack(), epochList.ovationExport);
    end
  end

  methods(Access=private)
    function self = EpochTree(epochList, forwardKeyStack, backwardKeyStack, ...
        ovationExport)

      import auimodel.*;

      self.isLeaf = 0;

      % base case: there are no paths left to split on
      if forwardKeyStack.count() == 0
        self.isLeaf = 1;
        self.epochList = epochList;

        % parent will set our other properties
        return;
      end

      % advance to the next key and store it for rewinding later.
      self.splitKey = forwardKeyStack.pop();
      backwardKeyStack.push(self.splitKey);

      childMap = Map();

      pathGetter = epochList.firstValue.keyPathGetter(self.splitKey, false);

      % use the distinct values of the key path to construct buckets to
      % drop the corresponding epochs into.

      iter = epochList.iterator('cell');
      while iter.hasNext
        epoch = iter.nextValue;

        try
          splitValue = pathGetter(epoch);
        catch e
          splitValue = Null.instance;
        end

        if ~ childMap.hasKey(splitValue)
          subList = EpochList();
          subList.ovationExport = epochList.ovationExport;

          childMap.set(splitValue, subList);
        end

        childMap.get(splitValue).append(epoch);
      end

      splitValues = childMap.keys();
      child_count = length(splitValues);

      % TODO: consider possible benefits of converting to a LinkedList.
      children = cell(child_count, 1);

      leafNodes = LinkedList();
      for i = 1 : child_count
        splitValue = splitValues{i};

        list = childMap.get(splitValue);
        list.populateStreamNames();

        % recursive case: we traverse along our children depth first
        node = EpochTree(list, forwardKeyStack, backwardKeyStack, ...
          ovationExport);

        node.parent = self;
        node.ovationExport = ovationExport;
        node.splitValue = splitValue;

        % if we're a deepest branch, append our leaf children.  otherwise
        % aggregate from child branches.
        if node.isLeaf
          leafNodes.append(node);
        else
          leafNodes.append(node.leafNodes, true);
        end

        children{i} = node;

        % rewind the key stack for our parent.  leafs don't descend and so
        % don't require a rewind.
        if ~ node.isLeaf
          forwardKeyStack.push( backwardKeyStack.pop() );
        end
      end

      self.children = children;
      self.leafNodes = leafNodes.elements;
      self.ovationExport = epochList.ovationExport;
    end
  end

  methods
    function skpath = get.splitValues(self)
      import auimodel.Map;

      skpathMap = Map();
      self.splitKeyPathHelper(skpathMap);

      skpath = skpathMap;
    end

    function splitKeyPathHelper(self, skpathMap)
      if isobject(self.parent)
        % recursive case: there are still parents to transcend along
        self.parent.splitKeyPathHelper(skpathMap);
        skpathMap.set(self.parent.splitKey, self.splitValue);
      end
    end

    function accept(self, visitor)
      % base case: visitor only visits us
      visitor.visit(self);

      % recursive case: visitor should visit our children as well
      for i = 1 : length(self.children)
        self.children{i}.accept(visitor);
      end
    end

    % CODE DEBT: should be factored out into a general utility class 'tostr'
    % method.
    function str = get.dotEdgeLabel(self)
      import auimodel.Util;
      str = Util.tostr(self.splitValue);
    end

    function label = get.dotLabel(self)
      if self.isLeaf
        label = sprintf('[%d epochs]', self.epochList.length());
      else
        label = self.splitKey;
      end
    end

    function visualize(self)
      import auimodel.TreeDotBuilder;

      builder = TreeDotBuilder();
      builder.buildDot(self);
      builder.displayDiagram();
    end

    function saveTree(self, filename)
      for i = 1 : length(self.leafNodes)
        leaf = self.leafNodes{i};
        leaf.epochList.flush;
        leaf.epochList.disableLazyLoads;
      end

      % pseudo guid to avoid namespace collisions.
      v00bd9f99c3f4e7b42bfe6f1f5e240e77 = self;

      save(filename, 'v00bd9f99c3f4e7b42bfe6f1f5e240e77');
    end
  end

  methods(Static)
    function tree = loadTree(filename)
      import auimodel.ExportLoader;

      load(filename);
      tree = v00bd9f99c3f4e7b42bfe6f1f5e240e77;

      ExportLoader.configureStores(tree.ovationExport);

      for i = 1 : length(tree.leafNodes)
        leaf = tree.leafNodes{i};
        leaf.epochList.enableLazyLoads;
      end

      clear v00bd9f99c3f4e7b42bfe6f1f5e240e77;
    end
  end
end

