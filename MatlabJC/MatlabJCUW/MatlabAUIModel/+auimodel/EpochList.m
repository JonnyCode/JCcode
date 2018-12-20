classdef EpochList < auimodel.LinkedList
% CODE DEBT: should be a subclass of List -- ran into issues with constructor
% inheritance.  conforms to the same interface w/ code dupe for now.

  properties
    stimuliStreamNames;
    responseStreamNames;
    keywords;
  end

  properties(Hidden)
    ovationExport;
  end

  methods
    function self = EpochList(ovationExport)
      import auimodel.*

      self.keywords = Set();

      if nargin < 1
        return;
      end

      ExportLoader.configureStores(ovationExport);
      proxy = auimodel.CoreDataProxy.instance('auimodel');

      epoch_count = length(ovationExport.epochURIs);
      epochs = cell(epoch_count, 1);

      for i = 1 : epoch_count
        entity = proxy.getEntity(ovationExport.epochURIs{i});
        epoch = Epoch(entity);

        for tag = epoch.keywords.elements
          inner = tag{1};
          self.keywords.add(inner);
        end

        self.append(epoch);
      end

      self.ovationExport = ovationExport;

      self.isDirty = true;
      self.populateStreamNames();
    end

    function populateStreamNames(self)
      firstEpoch = self.firstValue();
      self.responseStreamNames = fields(firstEpoch.responses);
      self.stimuliStreamNames = fields(firstEpoch.stimuli);
    end

    % precondition: all responses have same duration / number of data points.
    % precondition: all epochs have a response from passed streamName.
    function matrix = dataMatrixByStreamName(self, streamType, streamName, ...
        selectedOnly)

      firstEpoch = self.firstValue();

      points = size(firstEpoch.(streamType).(streamName).data);
      points = points(2);

      matrix = zeros(self.length, points);

      i = 1;
      iter = self.iterator;
      while iter.hasNext
        epoch = iter.nextValue;

        % ??? Operands to the || and && operators must be convertible to logical
        % scalar values.
        % if selectedOnly && ~ epoch.isSelected

        if selectedOnly
          if isempty(epoch.isSelected)
            continue;
          end
        end

        thesePoints = size(epoch.(streamType).(streamName).data);
        thesePoints = thesePoints(2);

        if ~ (points == thesePoints)
          disp(sprintf(strcat('Inconsistent data length in %dth epoch ', ...
            ' (was %d, expected %d) -- terminating.'), i, thesePoints, points));
           return
        end

        matrix(i, 1 : points) = epoch.(streamType).(streamName).data;
        i = i + 1;
      end
    end

    function responses = responsesByStreamName(self, streamName, selectedOnly)
      if nargin < 3
        selectedOnly = false;
      end

      responses = self.dataMatrixByStreamName('responses', streamName, ...
        selectedOnly);
    end

    % clears the cached stimulus and response data for the epochs in this
    % list.
    function flush(self)
      iter = self.iterator;
      while iter.hasNext
        epoch = iter.nextValue;
        self.flushDataForType(epoch, 'responses');
        self.flushDataForType(epoch, 'stimuli');
      end
    end

    function flushDataForType(self, epoch, streamType)
      streamNames = fields(epoch.(streamType));
      for i = 1 : length(streamNames)
        epoch.(streamType).(streamNames{i}).flush();
      end
    end

    function stimuli = stimuliByStreamName(self, streamName, selectedOnly)
      if nargin < 3
        selectedOnly = false;
      end

      stimuli = self.dataMatrixByStreamName('stimuli', streamName, ...
        selectedOnly);
    end

    function addKeywordTag(self, tag)
      import auimodel.CoreDataProxy;
      proxy = CoreDataProxy.instance('auimodel');

      iter = self.iterator;
      while iter.hasNext
        epoch = iter.nextValue;
        epoch.addKeywordTag(tag, false);
      end

      proxy.saveContext();

      % round trip to verify persistence.
      iter = self.iterator;
      while iter.hasNext
        epoch = iter.nextValue;
        epoch.loadKeywords();
      end

      self.keywords.add(tag);
    end

    function removeKeywordTag(self, tag)
      import auimodel.CoreDataProxy;
      proxy = CoreDataProxy.instance('auimodel');

      iter = self.iterator;
      while iter.hasNext
        epoch = iter.nextValue;
        epoch.removeKeywordTag(tag, false);
      end

      proxy.saveContext();

      self.keywords.remove(tag);
    end

    function enableLazyLoads(self)
      iter = self.iterator;
      while iter.hasNext
        epoch = iter.nextValue;
        epoch.enableLazyLoads;
      end
    end

    function disableLazyLoads(self)
      iter = self.iterator;
      while iter.hasNext
        epoch = iter.nextValue;
        epoch.disableLazyLoads;
      end
    end
  end
end

