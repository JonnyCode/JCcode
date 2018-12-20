classdef KeywordTaggableEntity < auimodel.NSManagedObject
  properties
    keywords;
  end

  methods(Static)
    function attrs = recurseAttributes(self)
      import auimodel.*;
      attrs = NSManagedObject.recurseAttributes();
    end
  end

  methods
    function self = KeywordTaggableEntity(entity)
      self = self@auimodel.NSManagedObject(entity);
    end

    function loadKeywords(self)
      import auimodel.*;

      self.keywords = Set();
      proxy = CoreDataProxy.instance('auimodel');
      keywordEntities = proxy.getEntityProperty(self.objectID, 'keywords');

      self.keywords = Set();
      for i = 1 : length(keywordEntities)
        keywordTag = KeywordTag(keywordEntities{i});
        self.keywords.add(keywordTag.tag);
      end
    end

    function addKeywordTag(self, tag, saveContext)
      import auimodel.CoreDataProxy;

      if ~ self.keywords.contains(tag)
        proxy = CoreDataProxy.instance('auimodel');
        proxy.addKeywordToEntity(tag, self.objectID);
        self.keywords.add(tag);

        if nargin < 3 || saveContext
          proxy.saveContext();
        end
      end
    end

    function removeKeywordTag(self, tag, saveContext)
      import auimodel.CoreDataProxy;

      if self.keywords.contains(tag)
        proxy = CoreDataProxy.instance('auimodel');
        proxy.removeKeywordFromEntity(tag, self.objectID);
        self.keywords.remove(tag);

        if nargin > 2 && ~ saveContext
          proxy.saveContext();
        end
      end
    end
  end
end

