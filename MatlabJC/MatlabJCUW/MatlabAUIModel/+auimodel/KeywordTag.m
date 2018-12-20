classdef KeywordTag < auimodel.NSManagedObject
  properties
    tag;
  end

  methods(Static)
    function register()
      import auimodel.*;

      proxy = CoreDataProxy.instance('auimodel');

      proxy.registerClass('KeywordTag', ...
        [NSManagedObject.recurseAttributes(), {'tag'}]);
    end
  end

  methods
    function self = KeywordTag(entity)
      self = self@auimodel.NSManagedObject(entity);
      self.tag = entity.tag;
    end
  end
end

