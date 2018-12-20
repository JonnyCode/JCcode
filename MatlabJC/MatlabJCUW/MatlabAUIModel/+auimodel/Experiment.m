classdef Experiment < auimodel.TaggableCacheableEntity
  properties
    purpose;
    otherNotes;
    startDate;
  end

  methods(Static)
    function register()
      import auimodel.*

      proxy = CoreDataProxy.instance('auimodel');
      proxy.registerClass('Experiment', ...
        [TaggableCacheableEntity.recurseAttributes(), ...
        {'purpose', 'otherNotes', 'startDate'}]);
    end
  end

  methods
    function self = Experiment(entity)
      self = self@auimodel.TaggableCacheableEntity(entity);
      self.purpose = entity.purpose;
      self.otherNotes = entity.otherNotes;
      self.startDate = datevec(entity.startDate);

      self.loadKeywords();
    end
  end
end

