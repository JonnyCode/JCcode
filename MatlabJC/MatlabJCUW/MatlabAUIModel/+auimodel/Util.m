classdef Util < handle
  methods(Static)
    function result = isEqual(thing1, thing2)
      result = false;

      if strcmp(class(thing1), class(thing2))
        if ischar(thing1)
          if strcmp(thing1, thing2)
            result = true;
          end
        elseif thing1 == thing2
          result = true;
        end
      end
    end

    function s = join(d,varargin)
      %S=JOIN(D,L) joins a cell array of strings L by inserting string D in
      %            between each element of L.  Meant to work roughly like the
      %            PERL join function (but without any fancy regular expression
      %            support).  L may be any recursive combination of a list 
      %            of strings and a cell array of lists.
      %
      %For any of the following examples,
      %    >> join('_', {'this', 'is', 'a', 'string'} )
      %    >> join('_', 'this', 'is', 'a', 'string' )
      %    >> join('_', {'this', 'is'}, 'a', 'string' )
      %    >> join('_', {{'this', 'is'}, 'a'}, 'string' )
      %    >> join('_', 'this', {'is', 'a', 'string'} )
      %the result is:
      %    ans = 
      %        'this_is_a_string'
      %
      %Written by Gerald Dalley (dalleyg@mit.edu), 2004

      import auimodel.Util;

      if (length(varargin) == 0), 
        s = '';
      else
        if (iscell(varargin{1}))
            s = Util.join(d, varargin{1}{:});
        else
            s = varargin{1};
        end
        
        for ss = 2:length(varargin)
            s = [s d Util.join(d, varargin{ss})];
        end
      end
    end

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

