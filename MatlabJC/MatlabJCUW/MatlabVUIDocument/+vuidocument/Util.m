classdef Util < handle
  methods(Static)
    function str = tostr(thing)
      if isnumeric(thing)
        str = num2str(thing);
      elseif ischar(thing)
        str = thing;
      elseif isobject(thing)
        % not sure how to check for a method, but this throws w/o, which is
        % desired.
        str = thing.toString();
      else
        disp(thing);
        MException('AUIModel:UnsupportedOperation', ...
          ['Do not know how to stringify preceding value.']).throw();
      end
    end
  end
end

