classdef TaggableCacheableEntity < auimodel.KeywordTaggableEntity
  properties
    isFromCache;
  end

  methods(Static)
    function attrs = recurseAttributes(self)
      import auimodel.*;
      attrs = KeywordTaggableEntity.recurseAttributes();
    end
  end

  methods
    function self = TaggableCacheableEntity(entity)
      self = self@auimodel.KeywordTaggableEntity(entity);

      persistent objectCache;

      if isempty(objectCache)
        objectCache = containers.Map();
      end

      if objectCache.isKey(self.objectID)
        self = objectCache(self.objectID);
        self.isFromCache = true;
      else
        objectCache(self.objectID) = self;
        self.isFromCache = false;
      end
    end
  end
end

