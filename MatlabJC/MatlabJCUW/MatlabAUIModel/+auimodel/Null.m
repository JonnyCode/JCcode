% null singleton to use as a splitValue for EpochTree when an epoch doesn't have
% a particular property being split on.
classdef Null < handle
  properties
    label = '[null]';
  end

  methods(Static)
    function instance = instance()
      persistent singleton;
      if isempty(singleton)
        singleton = auimodel.Null();
      end

      instance = singleton;
    end
  end

  methods
    function str = toString(self)
      str = '[null]';
    end
  end
end

