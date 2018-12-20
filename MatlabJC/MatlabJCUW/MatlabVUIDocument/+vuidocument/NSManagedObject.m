classdef NSManagedObject < handle
  properties
    objectID;
  end

  methods
    function self = NSManagedObject(entity)
      if nargin > 0
        self.objectID = entity.('objectID.URIRepresentation');
      end
    end
  end

  methods(Static)
    function attrs = recurseAttributes(self)
      attrs = {'objectID.URIRepresentation'};
    end
  end
end

