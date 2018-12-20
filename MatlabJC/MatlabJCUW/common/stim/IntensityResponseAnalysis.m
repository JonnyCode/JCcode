classdef IntensityResponseAnalysis < AnalysisSuper
    properties
        epochNums
        parameterList %first is X axis, second is separate curves
        lowPassFilter
    end
    properties (Hidden = true)
        SpatialStimParams
        D
        Nepochs
        uniqueVals
        Ind
    end
    
    methods
        
        function self = IntensityResponseAnalysis(params)
            self.initAnalysisSuper(params);
            
            %process parameters
            temp = textscan(params.parameterList,'%s','Delimiter',',');
            self.parameterList = temp{1};
            self.epochNums = str2num(params.epochNums);
            
            %self.SpatialStimParams = params.tempVar;
            self.SpatialStimParams = hdf5load([self.DataFilePath,'_spatial.h5']) ; % load spatial stim params
            
            SampleEpoch = self.SpatialStimParams.(['params_epoch' num2str(self.epochNums(1))]);
            prepts = round(SampleEpoch.spatial_prepts / 60 * 10000);
            self.getProcessedData(prepts);
            if ~isempty(params.lowPassFilter) && params.lowPassFilter > 0
                self.D = LowPassFilter(self.D, params.lowPassFilter, 1E-4);
            end
            self.getParameterIndices;
            self.plotResponseCurves;
            %plot(mean(D));
            
        end
        
        %deal with missing epochs (in h5 file) here and in
        %getParameterIndices
        function getProcessedData(self,prepts)
            self.Nepochs = length(self.epochNums);
            [SampleData err] = ITCReadEpoch(self.epochNums(1),0,self.fp);
            L = length(SampleData);
            D = zeros(self.Nepochs,L);
            
            %get data for each epoch, skipping ones with no data in the h5
            %file
            z=1;
            goodEpochs = zeros(1,self.Nepochs);
            for i=1:self.Nepochs
                if isfield(self.SpatialStimParams,['params_epoch' num2str(self.epochNums(i))])
                    [D(z,:), err] = ITCReadEpoch(self.epochNums(i),0,self.fp);
                    z=z+1;
                    goodEpochs(i) = 1;
                else
                   disp(['skipping epoch ' num2str(self.epochNums(i))]);
                end
            end
            baseline = mean(D(:,1:prepts),2);
            self.D = D - repmat(baseline,1,size(D,2));
                        
            self.epochNums = self.epochNums(goodEpochs==1);
            self.Nepochs = length(self.epochNums);
            %self.Nepochs
        end
        
        function getParameterIndices(self)
            P = length(self.parameterList);
            allVals = zeros(self.Nepochs, P);
            
            %get data and params fr each epoch
            for i=1:self.Nepochs
                curParams = self.SpatialStimParams.(['params_epoch' num2str(self.epochNums(i))]);
                for p=1:P
                    if ischar(self.parameterList{p})
                        V = curParams.(self.parameterList{p});
                    elseif strcmp(class(self.parameterList{p}), 'function_handle')
                        V = self.parameterList{p}(curParams);
                    end
                    if isnumeric(V)
                        allVals(i,p) = V;
                    else
                        error('non-numeric param values not yet supported');
                    end
                    
                end
            end
            
            self.uniqueVals = unique(allVals,'rows');
            for i=1:length(self.uniqueVals)
                self.Ind{i} = find(ismember(allVals, self.uniqueVals(i,:), 'rows'));
            end
        end
        
        function plotResponseCurves(self)
            param2Vals = unique(self.uniqueVals(:,2));
            C = colormap;            
            figure;
            for i=1:length(param2Vals)
                curSectionInd = find(self.uniqueVals(:,2) == param2Vals(i));
                Xvalues = self.uniqueVals(curSectionInd,1);
                curMeans = zeros(1,length(curSectionInd));
                curErrs = zeros(1,length(curSectionInd));
                for j=1:length(curSectionInd)
                    curMeans(j) = mean(max(abs(mean((self.D(self.Ind{curSectionInd(j)},:))))));
                    curErrs(j) = std(abs(max((self.D(self.Ind{curSectionInd(j)},:)))))./sqrt(length(self.Ind{curSectionInd(j)}));                    
                end
                h = errorbar(Xvalues, curMeans,curErrs,'kx-');
                set(h,'color',C(ceil(size(C,1)./length(param2Vals)*i),:));
                hold on;                               
            end
            legend(num2str(param2Vals));            
            drawnow;
            pause;
        end
        
    end
end
