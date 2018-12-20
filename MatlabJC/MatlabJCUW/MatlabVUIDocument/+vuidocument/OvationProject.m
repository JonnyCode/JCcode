classdef OvationProject < handle
  methods(Static)
    % CODE DEBT: needs to be factored out into util pkg.
    function directory = basename(path)
      nodes = textscan(path, '%s', 'delimiter', '/');
      nodes = nodes{1}(2:length(nodes{1}));
      nodes = strcat(nodes, '/');

      directory = '/';
      for i = 1 : length(nodes) - 1
        directory = strcat(directory, nodes(i));
      end

      directory = char(directory);
    end

    function open(ovationExport)
      import vuidocument.*

      if ~ isfield(ovationExport, 'projectURL')
        error('Open must be passed an OvationExport to find the project.')
      end
      
      proxy = CoreDataProxy.instance('vuidocument');

      % need to be able to remove a store / unload a previous project for this
      % to work.
      if proxy.isInitialized()
        error(['Currently only one ovation project per Matlab session ' ...
          'supported.']);
      end

      root = OvationProject.basename(which('OvationProject'));

      proxy.loadModelBundle(strcat(root, 'VUIModel.framework'));

      % TODO: use pwd -> cd -> bd trick to resolve relative and ~ paths.
      % precond: abs proj path
      if ( proxy.addStores({['file://localhost/' ovationExport.projectURL]}, ...
          'xml') )
        proxy.setIsInitialized();
      end

      AnalysisRecord.register();
    end
  end
end

