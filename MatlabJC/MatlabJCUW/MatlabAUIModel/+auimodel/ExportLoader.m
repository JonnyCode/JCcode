classdef ExportLoader < handle
  methods(Static)
    function isValid = isValidExport(ovationExport)
      isValid = isfield(ovationExport, 'persistentStoreURLs') ...
        && isfield(ovationExport, 'epochURIs');
    end

    function configureStores(ovationExport)
      import auimodel.*;

      proxy = CoreDataProxy.instance('auimodel');

      if ~ proxy.isInitialized()
        % locate myself
        root = strcat(fileparts(which('auimodel.ExportLoader')), '/');

        proxy.loadModelBundle(strcat(root, 'AUIModel.framework'));

        proxy.loadBundle(strcat(root, 'AUIProtocols.framework'));
        proxy.loadBundle(strcat(root, 'DAQFramework.framework'));

        proxy.loadPluginsAtPath(strcat(root, 'plugins'));

        proxy.performVoidStaticSelector('AUIIOController', ...
          'initializeControllers:', {}); % {} maps to nil

        proxy.setIsInitialized();

        Epoch.register();
      end

      proxy.addStores(ovationExport.persistentStoreURLs, 'sqlite');
    end
  end
end

