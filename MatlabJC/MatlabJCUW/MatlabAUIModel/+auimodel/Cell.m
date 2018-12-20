classdef Cell < auimodel.TaggableCacheableEntity
  properties
    comment;
    label;
    startDate;
    experiment;
  end

  methods(Static)
    function register()
      import auimodel.*
      proxy = CoreDataProxy.instance('auimodel');
      proxy.registerClass('Cell', ...
        [TaggableCacheableEntity.recurseAttributes(), ...
        {'comment', 'label', 'startDate', 'experiment'}]);
    end
  end

  methods
    function self = Cell(entity)
      self = self@auimodel.TaggableCacheableEntity(entity);
      self.comment  = entity.comment;
      self.label  = entity.label;
      self.startDate  = datevec(entity.startDate);

      import auimodel.*

      proxy = CoreDataProxy.instance('auimodel');

      entity = proxy.getEntityProperty(self.objectID, 'experiment');
      self.experiment = Experiment(entity);

      self.loadKeywords();
    end

    function string = toString(self)
      string = datestr(self.startDate, 'mm/dd/yy HH:MM:SS AM');
    end
  end
end

