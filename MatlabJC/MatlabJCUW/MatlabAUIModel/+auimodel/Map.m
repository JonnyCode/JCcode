% HACK: map is not really a typeof linkedlist, but the "Hidden" attribute of
% properties doesn't appear to get inherited.  so we have to use
% access=protected instead.
classdef Map < auimodel.LinkedList

% can't use a cell property due to weird performance problems.  still
% O(n) but w/ a linked list instead.
% NOTE: still useful despite containers.Map because can key by objects.
% TODO: currently only tested with string, object, and float keys.

  methods
    % CODE DEBT: overriding property setters with something having different
    % semanics...  may cause weird behavior someplace.
    function set(self, key, value)
      node = self.getNode(key);

      if isempty(node)
        self.append({key, value});
      else
        node.value = {key, value};
      end
    end

    function node = getNode(self, key)
      import auimodel.Util;

      node = [];

      % currently optimized for setting keys rarely and reading them often.
      % makes sense to iterate a cell in this case.
      iter = self.iterator('cell');
      while iter.hasNext
        cursor = iter.nextNode;

        if Util.isEqual(cursor.value{1}, key)
          node = cursor;
          break;
        end
      end
    end

    function value = get(self, key)
      value = [];

      node = self.getNode(key);
      if ~ isempty(node)
        value = node.value{2};
      end
    end

    function result = hasKey(self, key)
      result = ~ isempty(self.getNode(key));
    end

    function values = values(self)
      values = cell(self.length, 1);

      i = 1;
      iter = self.iterator;
      while iter.hasNext
        pair = iter.nextValue;
        values{i} = pair{2};
        i = i + 1;
      end
    end

    function keys = keys(self)
      keys = cell(self.length, 1);

      i = 1;
      iter = self.iterator;
      while iter.hasNext
        pair = iter.nextValue;
        keys{i} = pair{1};
        i = i + 1;
      end
    end

    function keyCount = count(self)
      keyCount = self.length;
    end

    function display(self)
      import auimodel.Util;

      display(sprintf('\nmap =\n'));

      if self.length() > 0
        iter = self.iterator;
        while iter.hasNext
          pair = iter.nextValue;
          display(['    ' Util.tostr(pair{1}) ': ' Util.tostr(pair{2})]);
        end
        display(' ');
      else
        display(sprintf('[empty]\n'));
      end
    end
  end
end

