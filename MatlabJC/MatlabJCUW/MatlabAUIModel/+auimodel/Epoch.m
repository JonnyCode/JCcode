classdef Epoch < auimodel.KeywordTaggableEntity
  properties
    comment;
    duration;
    protocolSettings;
    protocolID;
    startDate;
    responses;
    stimuli;
    cell;
    includeInAnalysis;
    isSelected;
  end

  properties(Hidden)
    lazyLoadsEnabled;
  end

  methods(Static)
    function register()
      import auimodel.*;

      proxy = CoreDataProxy.instance('auimodel');

      proxy.registerClass('Epoch', ...
        [KeywordTaggableEntity.recurseAttributes(), ...
        {'comment', 'duration', 'protocolSettings', 'protocolID', ...
         'startDate', 'includeInAnalysis'}]);

      Stimulus.register();
      Response.register();
      KeywordTag.register(); % common dependency for epoch, experiment, and cell
      Experiment.register();
      Cell.register();
    end
  end

  methods
    function set.includeInAnalysis(self, flag)
      if flag
        value = 1.0;
      else
        value = 0.0;
      end

      % when loading epochs from disk we can't / don't want to round-trip and
      % produce a noop change.
      if self.lazyLoadsEnabled
        proxy = auimodel.CoreDataProxy.instance('auimodel');

        proxy.setEntityProperty(self.objectID, 'includeInAnalysis', value);
        proxy.saveContext();

        self.includeInAnalysis = value;
      else
        self.includeInAnalysis = value;
      end
    end

    function self = Epoch(entity)
      self = self@auimodel.KeywordTaggableEntity(entity);
      self.comment = entity.comment;
      self.duration = entity.duration;
      self.protocolSettings = entity.protocolSettings;
      self.protocolID = entity.protocolID;
      self.startDate = datevec(entity.startDate);
      self.includeInAnalysis = entity.includeInAnalysis;

      import auimodel.*

      % TODO: remember and document why we have to explicity get these 
      % properties instead of their faults being fired by MatKit's array
      % builder -> auto retrieval.

      proxy = CoreDataProxy.instance('auimodel');
      entities = proxy.getEntityProperty(self.objectID, 'responses');

      self.responses = struct();

      for i = 1 : length(entities)
        response = Response(entities{i});
        self.responses.(response.streamName) = response;
      end

      entities = proxy.getEntityProperty(self.objectID, 'stimuli');

      self.stimuli = struct();

      for i = 1 : length(entities)
        stimulus = Stimulus(entities{i});
        self.stimuli.(stimulus.streamName) = stimulus;
      end

      entity = proxy.getEntityProperty(self.objectID, 'cell');
      self.cell = Cell(entity);

      self.loadKeywords();
    end

    function setLazyLoads(self, prop, flag)
      streams = fields(self.(prop));

      for i = 1 : length(streams)
        self.(prop).(streams{i}).lazyLoadsEnabled = flag;
        i = i + 1;
      end
    end

    function disableLazyLoads(self)
      self.setLazyLoads('responses', false);
      self.setLazyLoads('stimuli', false);
      self.lazyLoadsEnabled = false;
    end

    function enableLazyLoads(self)
      self.setLazyLoads('responses', true);
      self.setLazyLoads('stimuli', true);
      self.lazyLoadsEnabled = true;
    end
  end

  % CODE DEBT: could be refactored into a generalized base class for "key coded"
  % objects.  currently only added to Epoch because that's where we need it for
  % EpochTree.
  methods
    function functor = keyPathGetter(self, keyPath, verify)
      import auimodel.Queue;
      import auimodel.Util;

      nodeCell = textscan(keyPath, '%s', 'delimiter', '.');
      nodeCell = nodeCell{1};

      if nargin > 2 && ~ verify 
        code = ['.(''' Util.join(''').(''', nodeCell) ''')'];
      else
        pathQueue = Queue();
        pathQueue.enqueue(nodeCell, true);

        code = self.codeForPath(self, pathQueue);
      end

      if ~ isempty(code)
        functor = eval(strcat('@(self) self', code));
      else
        functor = [];
      end
    end

    function isvalid = isValidKey(self, node, key)
      % TODO: the following doesn't work in the case where both are false.  a
      % matlab bug i think.
      %isvalid = isfield(self, key) || self.findprop(key).isvalid;

      isvalid = 0;
      if isstruct(node)
        isvalid = isfield(node, key);
      elseif isobject(node)
        isvalid = node.findprop(key).isvalid;
      else
        MException('AUIModel:BadKeyedType', ...
          ['Encountered unkeyable type "' class(node) '".']);
      end
    end

    function code = codeForPath(self, node, pathQueue)
      % base case: out of nodes, path complete.
      if pathQueue.count() < 1
        code = '';
        return;
      end

      % using strcat here effects some kind of type conversion...
      tryKey = strcat(pathQueue.dequeue());

      if self.isValidKey(node, tryKey)
        % recursive case: build out code below this node
        belowCode = self.codeForPath(node.(tryKey), pathQueue);
        
        if ischar(belowCode)
          code = strcat(self.codeForKey(tryKey), belowCode);

        % at some point an invalid key was encountered, and so the entire path
        % is bad.
        else
          code = [];
        end
      else
        code = [];
      end
    end

    function code = codeForKey(self, key)
      code = strcat('.(''', key, ''')');
    end
  end
end

