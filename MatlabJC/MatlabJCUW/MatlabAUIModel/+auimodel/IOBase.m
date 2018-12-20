classdef IOBase < auimodel.NSManagedObject
  properties
    channelID;
    externalDeviceGain;
    externalDeviceMode;
    streamName;
    data;
    lazyLoadsEnabled = true;
  end

  methods(Static)
    function attrs = recurseAttributes()
      import auimodel.*

      attrs = [NSManagedObject.recurseAttributes(), {'channelID', ...
        'externalDeviceGain', 'externalDeviceMode', ...
        'streamInfoUserDescription'}];
    end
  end

  methods
    function self = IOBase(entity)
      self = self@auimodel.NSManagedObject(entity);
      self.channelID = entity.channelID;
      self.externalDeviceGain = entity.externalDeviceGain;
      self.externalDeviceMode = entity.externalDeviceMode;

      self.data = 'Lazy loads disabled.';

      streamName = entity.streamInfoUserDescription;
      streamName = strrep(streamName, '.', '');
      streamName = strrep(streamName, ' ', '_');

      self.streamName = streamName;
    end

    function flush(self)
      self.data = 'Lazy loads disabled.';

      import auimodel.CoreDataProxy;
      proxy = CoreDataProxy.instance('auimodel');

      % NOTE: here we assume all subclasses of IOBase adhere to the CachedData
      % protocol.
      proxy.performSelectorOnEntity(self.objectID, 'flushCachedData');
    end

    function value = get.data(self)
      import auimodel.*

      if self.lazyLoadsEnabled && ischar(self.data)
        proxy = CoreDataProxy.instance('auimodel');
        self.data = proxy.getEntityProperty(self.objectID, self.dataKeyPath());
      end

      value = self.data;
    end
  end
end

