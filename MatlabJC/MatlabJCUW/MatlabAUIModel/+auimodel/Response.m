classdef Response < auimodel.IOBase
  methods(Static)
    function register()
      import auimodel.*

      proxy = CoreDataProxy.instance('auimodel');
      proxy.registerClass('Response', IOBase.recurseAttributes());
    end
  end

  methods
    function self = Response(entity)
      self = self@auimodel.IOBase(entity);
    end

    function path = dataKeyPath(self)
      path = 'data';
    end
  end
end
 
