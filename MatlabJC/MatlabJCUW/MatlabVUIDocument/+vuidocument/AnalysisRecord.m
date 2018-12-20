classdef AnalysisRecord < vuidocument.NSManagedObject
  properties
    name;
    epochList;
    parameters;
    userDescription;
    svnRepositoryURI;
    codeRevision;
  end

  methods
    function set.svnRepositoryURI(self, svnRepositoryURI)
      % TODO: vfy svnRepositoryURI and active checkout here
      if ~ ischar(svnRepositoryURI)
        error('Repository URI must be a string / character array.');
      end

      self.svnRepositoryURI = svnRepositoryURI;
    end

    function uri = get.svnRepositoryURI(self)
      global AnalysisSvnRepositoryUri;
      if ~ isempty(self.svnRepositoryURI)
        uri = self.svnRepositoryURI;
      elseif exist('AnalysisSvnRepositoryUri') ...
          && ischar(AnalysisSvnRepositoryUri)
        uri = AnalysisSvnRepositoryUri;
      else
        uri = [];
      end
    end

    function set.codeRevision(self, codeRevision)
      if ~ isnumeric(codeRevision) || codeRevision < 0
        error('Revision numbers must be a positive number.');
      end

      self.codeRevision = codeRevision;
    end

    function set.name(self, name)
      if ~ ischar(name)
        error('Name must be a string / character array.');
      end

      self.name = name;
    end

    function set.userDescription(self, userDescription)
      if ~ (ischar(userDescription) || isempty(userDescription))
        error('Description must be a string or empty.');
      end

      self.userDescription = userDescription;
    end

    function set.parameters(self, parameters)
      if ~ isstruct(parameters)
        error('Parameters must be a struct (of analysis code parameters).')
      end

      self.parameters = parameters;
    end

    function set.epochList(self, epochList)
      if ~ isa(epochList, 'auimodel.EpochList')
        error('Epoch list must be an instance of auimodel.EpochList.');
      end

      self.epochList = epochList;
    end
  end

  properties(SetAccess=private)
    resources;
    keywords;
    timeStamp;
  end

  properties (GetAccess=private, SetAccess=private)
    resourceUris;
    autoTag;
  end

  methods(Static)
    function register()
      import vuidocument.*;

      proxy = CoreDataProxy.instance('vuidocument');

      proxy.registerClass('AnalysisRecord', ...
        [NSManagedObject.recurseAttributes(), ...
        {'name', 'userDescription', 'parameters', ... %'keywords', 'resources', 
         'codeRevision', 'svnRepositoryURI', 'timeStamp', 'autoTag'}]);
    end
  end

  methods
    function self = AnalysisRecord(name, parameters, epochList)
      import vuidocument.*

      if nargin < 3
        error('Name, analysis code parameters, and list of epochs required.');
      end

      proxy = CoreDataProxy.instance('vuidocument');

      if ~ proxy.isInitialized()
        error('Please open an Ovation project with OvationProject.open().');
      end

      self.name = name;
      self.parameters = parameters;
      self.epochList = epochList;

      self.keywords = Set();
      self.resources = Set();
      self.resourceUris = Map();
    end

    function addKeyWord(self, keyWord)
      self.keywords.add(keyWord);
    end

    function removeKeyWord(self, keyWord)
      self.keywords.remove(keyWord);
    end

    % precond: a resource with the filename of the one being added doesn't
    % already exist.
    function addResource(self, resource)
      % throws an error if the file doesn't exist
      ls(resource);

      fileUri = ['file://' pwd '/' resource];

      self.resources.add(resource);
      self.resourceUris.set(resource, fileUri);
    end

    function removeResource(self, resource)
      self.resources.remove(resource);
      self.resourceUris.delete(resource);
    end

    function entity = toEntity(self)
      entity = struct();
      entity.name = self.name;
      entity.userDescription = self.userDescription;
      entity.parameters = self.parameters;
      entity.keywords = self.keywords.list.elements;
      entity.resources = self.resourceUris.getValues();
      entity.svnRepositoryURI = self.svnRepositoryURI;
      entity.timeStamp = datevec(self.timeStamp);
      entity.codeRevision = self.codeRevision;
      entity.autoTag = self.autoTag;

      % TODO: make this more intelligent by having the user set it or
      % automatically determining a good root path.
      entity.resourceRootURL = 'file:///';
    end

    function populateFromEntity(self, entity)
      % CODE DEBT: breaking encapsulation here.  NSManagedObject should
      % probably have some kind of "populate" method other than the
      % constructor.
      self.objectID = entity.('objectID.URIRepresentation');

      % CODE DEBT: factor this pattern out into utility framework for populating
      % entities retrieved from the database (via some metaprogramming).  also
      % maybe marshaling into structs like the above.

      self.name = entity.name;
      self.userDescription = entity.userDescription;
      self.parameters = entity.parameters;
      self.codeRevision = entity.codeRevision;
      self.svnRepositoryURI = entity.svnRepositoryURI;
      self.timeStamp = entity.timeStamp;
      self.autoTag = entity.autoTag;

      % TODO: create proxied KeywordTag (should go in pending util lib) and
      % Resource classes, breaking out of them as appropriate to set these
      % properties.
%      self.keywords = entity.keywords;
%      self.resources = entity.resources;
    end

    function commit(self)
      import vuidocument.*

      if ~ isempty(self.objectID)
        error('This record has already been commited.');
      end

      if isempty(self.svnRepositoryURI)
        error(['SVN repository root not set.  Either set globally in ' ...
          'AnalysisSvnRepositoryUri or on svnRepositoryURI of this instance.']);
      end

      if isempty(self.codeRevision)
        error('Revision of analysis code not set in codeRevision.');
      end

      proxy = CoreDataProxy.instance('vuidocument');

      entity = proxy.performStaticSelector('AnalysisRecord', ...
        'insertAnalysisRecordForDictionary:intoManagedObjectContext:', ...
        self.toEntity());

      % round trip to allow easy verification of integrity.
      self.populateFromEntity(entity);

      self.epochList.addKeywordTag(self.autoTag);

      for tag = self.keywords.list.elements
        tag = tag{1};
        self.epochList.addKeywordTag(tag);
      end

      proxy.saveContext();
    end
  end
end

