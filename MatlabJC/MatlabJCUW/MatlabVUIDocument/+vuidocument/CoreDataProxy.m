classdef CoreDataProxy < handle
  properties
    modelName;
  end

  methods(Static)
    function instance = instance(myPackage)
      persistent instances;
      if isempty(instances)
        instances = containers.Map();
      end

      if ~ instances.isKey(myPackage)
        import(strcat(myPackage, '.CoreDataProxy'))

        instance = CoreDataProxy();
        instance.modelName = myPackage;

        instances(myPackage) = instance;
      else
        instance = instances(myPackage);
      end
    end
  end

  methods
    function loadModelBundle(self, bundlePath)
      CoreDataProxyMex(self.modelName, 'loadModelBundle', {bundlePath});
    end

    function setIsInitialized(self)
      CoreDataProxyMex(self.modelName, 'setIsInitialized', {});
    end

    function isInit = isInitialized(self)
      isInit = CoreDataProxyMex(self.modelName, 'getIsInitialized', {});
    end

    function loadPluginsAtPath(self, pluginsPath)
      CoreDataProxyMex(self.modelName, 'loadPluginsAtPath', {pluginsPath});
    end

    function result = performStaticSelector(self, class, selector, param)
      result = CoreDataProxyMex(self.modelName, 'performStaticSelector', ...
        {class, selector, param});
    end

    function performVoidStaticSelector(self, class, selector, param)
      CoreDataProxyMex(self.modelName, 'performVoidStaticSelector', ...
        {class, selector, param});
    end

    function performSelectorOnEntity(self, uri, selector)
      CoreDataProxyMex(self.modelName, 'performSelectorOnEntity', ...
        {uri, selector});
    end

    function loadBundle(self, bundlePath)
      CoreDataProxyMex(self.modelName, 'loadBundle', {bundlePath});
    end

    function outcome = addStores(self, storePaths, storeType)
      outcome = ...
        CoreDataProxyMex(self.modelName, 'addStores', {storePaths, storeType});
    end

    function registerClass(self, className, properties)
      CoreDataProxyMex(self.modelName, 'registerClass', ...
        [className, properties]);
    end

    function entity = getEntity(self, entityUri)
      entity = CoreDataProxyMex(self.modelName, 'getEntity', {entityUri});
    end

    function property = getEntityProperty(self, entityUri, propertyName)
      params = {entityUri, propertyName};
      property = CoreDataProxyMex(self.modelName, 'getEntityProperty', params);
    end

    function resetContext(self)
      CoreDataProxyMex(self.modelName, 'resetContext', {});
    end

    function saveContext(self)
      CoreDataProxyMex(self.modelName, 'saveContext', {});
    end

    function addKeywordToEntity(self, tag, entityUri)
      CoreDataProxyMex(self.modelName, 'addKeywordToEntity', {tag, entityUri});
    end

    function removeKeywordFromEntity(self, tag, entityUri)
      CoreDataProxyMex(self.modelName, 'removeKeywordFromEntity', ...
        {tag, entityUri});
    end

    function setEntityProperty(self, entityUri, propertyName, value)
      params = {entityUri, propertyName, value};
      CoreDataProxyMex(self.modelName, 'setEntityProperty', params);
    end
  end
end

