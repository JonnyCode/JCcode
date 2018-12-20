classdef ModelExperiment < handle
    properties
        %input properties:
        %for all params stucts, if field is a cell array, it will be
        %considered a list of parameters to run in the model. All
        %combinations of the parameters will be tested. runsPerCondition is the
        %number of runs in each condition, nTrials (in the model itself)
        %will give an additional factor
        modelInitParams %struct
        modelFilterParams %cell array of structs, one for each layer
        modelNLParams %cell array of structs, one for each layer
        modelStimParams %struct
        runsPerCondition %random numbers in stim change each run
        
        modelInitFunction %function handle, returns Model (with layers and subunits set, transforms as well: ok if other stuff is not)
        setFilterFunctions %cell array of function handles (one for each layer), for setting filters
        setNLFunctions %cell array of function handles (one for each layer), for setting nonlinearities
        setStimFunction %function handle for setting stimulus
        
        %to be calculated during experiment:
        conditionMatrix %matrix (dimensioned by the number of variable params) of param structs. Once the model is run, a results field will be added to each struct        
        variableParams
    end
    
    properties (Hidden = true, Access = private)
    end
    
    methods        
        function buildconditionMatrix(self)
            %variableParams = struct;
            variableParams = [];
            z = 0; %counter for total number of variable params
            
            %init
            fnames = fieldnames(self.modelInitParams);
            for f=1:length(fnames)
                curFieldValue = self.modelInitParams.(fnames{f});
                if iscell(curFieldValue) && length(curFieldValue) > 1
                    z=z+1;
                    variableParams(z).name = fnames{f};
                    variableParams(z).type = 'init';
                    variableParams(z).values = curFieldValue;
                    variableParams(z).len = length(curFieldValue);
                end
            end
            
            %filter
            nLayers = length(self.modelFilterParams);
            for l=1:nLayers
                fnames = fieldnames(self.modelFilterParams{l});
                for f=1:length(fnames)
                    curFieldValue = self.modelFilterParams{l}.(fnames{f});
                    if iscell(curFieldValue) && length(curFieldValue) > 1
                        z=z+1;
                        variableParams(z).name = fnames{f};
                        variableParams(z).type = 'filt';
                        variableParams(z).values = curFieldValue;
                        variableParams(z).len = length(curFieldValue);
                        variableParams(z).layer = l;
                    end
                end
            end
            
            %NL
            for l=1:nLayers
                fnames = fieldnames(self.modelNLParams{l});
                for f=1:length(fnames)
                    curFieldValue = self.modelNLParams{l}.(fnames{f});
                    if iscell(curFieldValue) && length(curFieldValue) > 1
                        z=z+1;
                        variableParams(z).name = fnames{f};
                        variableParams(z).type = 'nlin';
                        variableParams(z).values = curFieldValue;
                        variableParams(z).len = length(curFieldValue);
                        variableParams(z).layer = l;
                    end
                end
            end
            
            %stim
            fnames = fieldnames(self.modelStimParams);
            for f=1:length(fnames)
                curFieldValue = self.modelStimParams.(fnames{f});
                if iscell(curFieldValue) && length(curFieldValue) > 1
                    z=z+1;
                    variableParams(z).name = fnames{f};
                    variableParams(z).type = 'stim';
                    variableParams(z).values = curFieldValue;
                    variableParams(z).len = length(curFieldValue);
                end
            end                    
            
           %got all variable params, now building cond. matrix
           L = 1; %length of each dimension
           for i=1:length(variableParams)               
            L(i) = variableParams(i).len;
           end          
           
           self.variableParams = variableParams;
                                 
           ND = length(variableParams); %number of dimensions (only 1,2,3 for now)
           
           if ND>3
               error('More than 3 variable params not supported');
           elseif ND == 0
               disp('No variable params');
               self.conditionMatrix = cell(1,1);
               self.conditionMatrix{1} = struct; %empty struct
           end
           
           %1D
           if ND == 1
              self.conditionMatrix = cell(1,L(1)); %initialize to correct size 
              for i=1:L(1)
                 self.conditionMatrix{i}.(variableParams(1).name) = variableParams(1).values{i};
                 self.conditionMatrix{i}.(['type_' variableParams(1).name]) = variableParams(1).type;
                 if isfield(variableParams(1), 'layer')
                     self.conditionMatrix{i}.(['layer_' variableParams(1).name]) = variableParams(1).layer;
                 end
              end
           end %1D
           
           %2D
           if ND == 2
               for i=1:L(1)
                   for j=1:L(2)
                       %param1
                       self.conditionMatrix{i,j}.(variableParams(1).name) = variableParams(1).values{i};
                       self.conditionMatrix{i,j}.(['type_' variableParams(1).name]) = variableParams(1).type;
                       if isfield(variableParams(1), 'layer')
                           self.conditionMatrix{i,j}.(['layer_' variableParams(1).name]) = variableParams(1).layer;
                       end
                       
                       %param2
                       self.conditionMatrix{i,j}.(variableParams(2).name) = variableParams(2).values{j};
                       self.conditionMatrix{i,j}.(['type_' variableParams(2).name]) = variableParams(2).type;
                       if isfield(variableParams(2), 'layer')
                           self.conditionMatrix{i,j}.(['layer_' variableParams(2).name]) = variableParams(2).layer;
                       end
                   end
               end
           end %2D
           
           %3D
           if ND == 3
               for i=1:L(1)
                   for j=1:L(2)
                       for k=1:L(3)
                           %param1
                           self.conditionMatrix{i,j,k}.(variableParams(1).name) = variableParams(1).values{i};
                           self.conditionMatrix{i,j,k}.(['type_' variableParams(1).name]) = variableParams(1).type;
                           if isfield(variableParams(1), 'layer')
                               self.conditionMatrix{i,j,k}.(['layer_' variableParams(1).name]) = variableParams(1).layer;
                           end
                           
                           %param2
                           self.conditionMatrix{i,j,k}.(variableParams(2).name) = variableParams(2).values{j};
                           self.conditionMatrix{i,j,k}.(['type_' variableParams(2).name]) = variableParams(2).type;
                           if isfield(variableParams(2), 'layer')
                               self.conditionMatrix{i,j,k}.(['layer_' variableParams(2).name]) = variableParams(2).layer;
                           end
                           
                           %param3
                           self.conditionMatrix{i,j,k}.(variableParams(3).name) = variableParams(3).values{k};
                           self.conditionMatrix{i,j,k}.(['type_' variableParams(3).name]) = variableParams(3).type;
                           if isfield(variableParams(3), 'layer')
                               self.conditionMatrix{i,j,k}.(['layer_' variableParams(3).name]) = variableParams(3).layer;
                           end
                       end
                   end
               end
           end %3D            
        end %buildCondition matrix
        
        function run(self) %run the experiment, storing results in conditionMatrix 
            if isempty(self.variableParams)
                allTypes = [];
            else
                allTypes = {self.variableParams.type};
            end
                
            modelState = cell(1,4);
            
            if isempty(strmatch('init',allTypes,'exact')) %no variable init params
                Model = self.modelInitFunction(self.modelInitParams);
                modelState{1} = 1;
                modelState{4} = 0;
                nLayers = length(Model.layers);
                modelState{2} = zeros(1,nLayers);
                modelState{3} = zeros(1,nLayers);
            end            
            
            %if initialized
            if modelState{1} == 1
                
                if isempty(strmatch('filt',allTypes,'exact')) %no variable filt params
                    for i=1:length(Model.layers)
                        curF = self.setFilterFunctions{i};
                        Model.layers(i) = curF(Model.layers(i),mergeStruct(self.modelInitParams, self.modelFilterParams{i}));
                        modelState{2}(i) = 1;
                    end
                end
                
                if isempty(strmatch('nlin',allTypes,'exact')) %no variable nlin params
                    for i=1:length(Model.layers)
                        curF = self.setNLFunctions{i};
                        Model.layers(i) = curF(Model.layers(i),mergeStruct(self.modelInitParams, self.modelNLParams{i}));
                        modelState{2}(i) = 1;
                    end
                end
                
                if isempty(strmatch('stim',allTypes,'exact')) %no variable stim params
                    Model = self.setStimFunction(Model,mergeStruct(self.modelInitParams, self.modelStimParams));
                    currentStimParams = mergeStruct(self.modelInitParams, self.modelStimParams); %for stim regeneration
                end
            end
            %if everything set, run model here, store result
            
            %get state of model at this point:
            %state = 4 element cell array, 1 = model part done
            %for elements with layers (filt and nlin), the cell has a 0 for
            %1 for each layer
            %[init, filt, nlin, stim]
                        
            N = numel(self.conditionMatrix);
            sZ = size(self.conditionMatrix);
            
            defaultModelState = modelState;
            
            disp(['Looping over ' num2str(N) ' conditions']);
            for i=1:N %loop over all variable conditions
                modelState = defaultModelState;
                currentInitParams = self.modelInitParams; %default
                
                curCondField = self.conditionMatrix{ind2sub(sZ,i)};
                fnames = fieldnames(curCondField);                                
                if isempty(modelState{1}) || modelState{1} == 0 %model not initialized
                    currentInitParams = mergeStruct(self.modelInitParams, curCondField); %overwrites the ones that need to change
                    Model = self.modelInitFunction(currentInitParams);                    
                    modelState{1} = 1;
                    modelState{4} = 0;
                    nLayers = length(Model.layers);
                    modelState{2} = zeros(1,nLayers);
                    modelState{3} = zeros(1,nLayers);
                end
                nLayers = length(Model.layers);
                                                
                typeFieldsInd = strmatch('type_',fnames);
                layerFieldsInd = strmatch('layer_',fnames);
                
                %initialized by now                
                %add filters if not added
                for l=1:nLayers
                    if modelState{2}(l) == 0
                        curF = self.setFilterFunctions{l};
                        for f=1:length(typeFieldsInd)
                            if ~isempty(strmatch('filt',curCondField.(fnames{typeFieldsInd(f)}))) && curCondField.(fnames{layerFieldsInd(f)}) == l
                                [null, paramName] = strtok(fnames{typeFieldsInd(f)},'_');
                                paramName = paramName(2:end); %remove leading '_'
                                tempFilterParams = self.modelFilterParams{l};
                                tempFilterParams.(paramName) = curCondField.(paramName);
                                Model.layers(l) = curF(Model.layers(l),mergeStruct(currentInitParams, tempFilterParams));
                                modelState{2}(l) = 1;
                                break; %although this should not be necessary is it should only find the right parameter once
                            end
                        end
                        if modelState{2}(l) == 0;
                            %set default (other layers)
                            Model.layers(l) = curF(Model.layers(l),mergeStruct(currentInitParams, self.modelFilterParams{l}));
                            modelState{2}(l) = 1;
                        end
                    end       
                end
                            
                %add nlins if not added
                for l=1:nLayers
                    if modelState{3}(l) == 0
                        curF = self.setNLFunctions{l};
                        for f=1:length(typeFieldsInd)
                            if ~isempty(strmatch('nlin',curCondField.(fnames{typeFieldsInd(f)}))) && curCondField.(fnames{layerFieldsInd(f)}) == l
                                %param is in current condition (to be set)
                                [null, paramName] = strtok(fnames{typeFieldsInd(f)},'_');
                                paramName = paramName(2:end); %remove leading '_'
                                tempNLParams = self.modelNLParams{l};
                                tempNLParams.(paramName) = curCondField.(paramName);
                                Model.layers(l) = curF(Model.layers(l),mergeStruct(currentInitParams, tempNLParams));
                                modelState{3}(l) = 1;                                
                                break; %although this should not be necessary is it should only find the right parameter once
                            end
                        end
                        if modelState{3}(l) == 0 %param not found
                            %set default (other layers)
                            Model.layers(l) = curF(Model.layers(l),mergeStruct(currentInitParams, self.modelNLParams{l}));
                            modelState{3}(l) = 1;
                        end
                    end
                end
                      
                %add stim if not added
                if modelState{4} == 0
                    curF = self.setStimFunction;
                    for f=1:length(typeFieldsInd)
                        if ~isempty(strmatch('stim',curCondField.(fnames{typeFieldsInd(f)})))
                            %param is in current condition (to be set)
                            [null, paramName] = strtok(fnames{typeFieldsInd(f)},'_');
                            paramName = paramName(2:end); %remove leading '_'
                            tempStimParams = self.modelStimParams;
                            tempStimParams.(paramName) = curCondField.(paramName);
                            Model = curF(Model,mergeStruct(currentInitParams, tempStimParams));
                            currentStimParams = mergeStruct(currentInitParams, tempStimParams); %for stim regeneration
                            modelState{4} = 1;
                            break; %although this should not be necessary is it should only find the right parameter once
                        end
                    end
                    if modelState{4} == 0;
                        %set default (other layers)
                        Model = curF(Model,mergeStruct(currentInitParams, self.modelStimParams));
                        currentStimParams = mergeStruct(currentInitParams, self.modelStimParams); %for stim regeneration
                        modelState{4} = 1;
                    end
                end               
                
                disp(['Running condition ' num2str(i) ' of ' num2str(N)]);
                for r=1:self.runsPerCondition                                     
                    if r>1 %regenerate stim            
                        disp(['generating stim ' num2str(r) ' of ' num2str(self.runsPerCondition)]);   
                        curF = self.setStimFunction;
                        Model = curF(Model,currentStimParams);
                    end    
                    disp('done');
                    Model.run;
                    self.conditionMatrix{ind2sub(sZ,i)}.Model = Model;
                    self.conditionMatrix{ind2sub(sZ,i)}.result{r} = Model.output;                
                end
                disp('done');
            end %for all N (conditions)                                    
        end
        
        
    end %methods    
end %class