classdef Stimulus < auimodel.IOBase
  properties
    duration;
    parameters;
    stimulusID;
    version;
    description;
    sampleRate;
  end

  methods(Static)
    function register()
      import auimodel.*

      proxy = CoreDataProxy.instance('auimodel');

      proxy.registerClass('Stimulus', [IOBase.recurseAttributes(), ...
        {'duration', 'parameters', 'stimulusID', 'version', 'description', ...
         'sampleRate'}]);
    end
  end

  methods
    function self = Stimulus(entity)
      self = self@auimodel.IOBase(entity);
      self.objectID = entity.('objectID.URIRepresentation');
      self.duration = entity.duration;
      self.parameters = entity.parameters;
      self.stimulusID = entity.stimulusID;
      self.version = entity.version;
      self.description = entity.description;
      self.sampleRate = entity.sampleRate;
    end

    function path = dataKeyPath(self)
      path = 'reconstructedData';
    end
  end
end

