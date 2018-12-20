classdef Map < handle
% primitive O(n) Map implementation backed by a cell array.  optimizations would
% include a cell-backed hash table in matlab, then in a mex file using an array
% of mxArray pointers (could be quicker as hash code generators are probably
% already available for primitives?).

% TODO: currently only tested with string, object, and float keys.

% CODE DEBT: newer than one in auimodel.., need to be merged -> util lib.

  properties
    pairs;
  end

  methods
    function self = Map()
      self.pairs = {};
    end

    % CODE DEBT: overriding property setters with something having different
    % semanics...  may cause weird behavior someplace.
    function set(self, key, value)
      index = self.getIndex(key);

      if index == 0
        self.pairs{ length(self.pairs) + 1 } = {key, value};
      else
        self.pairs{index}{2} = value;
      end
    end

    function delete(self, key)
      index = self.getIndex(key);

      if index > 0
        self.pairs = {self.pairs{1 : index - 1} self.pairs{index + 1 : end}};
      end
    end

    function values = getValues(self)
      import vuidocument.List;

      list = List();
      for i = 1 : length(self.pairs)
        list.append(self.pairs{i}{2});
      end

      values = list.elements;
    end

    function index = getIndex(self, findKey)
      index = 0;

      for i = 1 : length(self.pairs)
        haveKey = self.pairs{i}{1};

        if strcmp(class(findKey), class(haveKey))
          if ischar(findKey)
            if strcmp(findKey, haveKey)
              index = i;
            end
          elseif findKey == haveKey
            index = i;
          end
        end
      end
    end

    % CODE DEBT: overriding property getters with something having different
    % semanics...  may cause weird behavior someplace.
    function value = get(self, key)
      value = [];
      index = self.getIndex(key);

      if ~ index == 0
        value = self.pairs{index}{2};
      end
    end

    function result = hasKey(self, key)
      result = ~ self.getIndex(key) == 0;
    end

    function keys = keys(self)
      import vuidocument.List;

      keys = List();

      for i = 1 : length(self.pairs)
        keys.append(self.pairs{i}{1});
      end

      keys = keys.elements;
    end

    function keyCount = count(self)
      keyCount = length(self.pairs);
    end

    function display(self)
      import vuidocument.Util;

      display(sprintf('\nmap =\n'));

      if self.count() > 0

        for i = 1 : length(self.pairs)
          pair = self.pairs{i};
          display(['    ' Util.tostr(pair{1}) ': ' Util.tostr(pair{2})]);
        end
        display(' ');
      else
        display(sprintf('[empty]\n'));
      end

    end

    % seems to get caught up w/ method invocation
%    function value = subsref(self, key)
%      value = self.get(key);
%    end
%
%    function subsasgn(self, key, value)
%      self.set(key, value);
%    end
  end
end

