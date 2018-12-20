classdef PairedSpotsAnalysis < AnalysisSuper
    properties
        epochNums
        parameterList %splits other than position and pairSpots
        lowPassFilter
    end
    properties (Hidden = true)
        SpatialStimParams
        D
        Nepochs
        uniqueVals
        Ind
        position
        pairSpots
    end
    
    methods
        
        function self = PairedSpotsAnalysis(params)
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
            self.makeLinearityPlots;
            
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
            self.position = zeros(self.Nepochs,4);
            self.pairSpots = zeros(self.Nepochs,1);
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
                self.pairSpots(i) = curParams.pairSpots;
                if curParams.pairSpots
                    self.position(i,:) = [curParams.Xpos1 curParams.Ypos1 curParams.Xpos2 curParams.Ypos2];
                else
                    self.position(i,:) = [curParams.Xpos1 curParams.Ypos1 0 0];
                end
            end
            
            self.uniqueVals = unique(allVals,'rows');
            for i=1:length(self.uniqueVals)
                self.Ind{i} = find(ismember(allVals, self.uniqueVals(i,:), 'rows'));
            end
        end
        
        function makeLinearityPlots(self)
            C = colormap;
            L = size(self.D,2);
            for i=1:length(self.Ind) %for each segment
                figure(i);
                %segParamVals = self.allParamVals(self.Ind{i},:);
                segD = self.D(self.Ind{i},:);
                segPos = self.position(self.Ind{i},:);
                segPairSpots = self.pairSpots(self.Ind{i});
                
                pairInd = find(segPairSpots==1);
                pairD = segD(pairInd,:);
                pairPos = segPos(pairInd,:);
                uniquePos = unique(pairPos,'rows');
                
                singleInd = find(segPairSpots==0);
                singleD = segD(singleInd,:);
                singlePos = segPos(singleInd,:);
                
                LI = [];
                dist = [];
                for j=1:size(uniquePos,1)
                    posInd = [];
                    pos1Ind = [];
                    pos2Ind = [];
                    z=1;
                    for k=1:length(pairPos)
                        if isequal(pairPos(k,:), uniquePos(j,:))
                            posInd(z) = k;
                            z=z+1;
                        end
                    end
                    curD = pairD(posInd,:);
                    curMean = mean(curD,1);
                    %err = std(curD,[],1);
                    %get single components
                    pos1 = uniquePos(j,1:2);
                    pos2 = uniquePos(j,3:4);
                    z1=1;
                    z2=1;
                    for k=1:size(singlePos,1)
                        if isequal(singlePos(k,1:2), pos1)
                            pos1Ind(z1) = k;
                            z1=z1+1;
                        elseif isequal(singlePos(k,1:2), pos2)
                            pos2Ind(z2) = k;
                            z2=z2+1;
                        end
                    end
                    curD1 = singleD(pos1Ind,:);
                    curMean1 = mean(curD1,1);
                    curD2 = singleD(pos2Ind,:);
                    curMean2 = mean(curD2,1);
                    
                    LI(j) = (max(abs(curMean)) - max(abs(curMean1+curMean2))) / max(abs(curMean1+curMean2));
                    x1 = uniquePos(j,1); y1 = uniquePos(j,2);
                    x2 = uniquePos(j,3); y2 = uniquePos(j,4);
                    dist(j) = sqrt((x1-x2)^2+(y1-y2)^2)*1.8; %could add this pixelFactor (1.8) as a param
                    
                    subplot2(3,2,j);
                    plot(1:L,curMean,'r-');
                    hold on
                    plot(1:L,curMean1,'k-');
                    plot(1:L,curMean2,'k-');
                    plot(1:L,curMean1+curMean2,'b-');
                    title([num2str(self.uniqueVals(i)) ': ' num2str(uniquePos(j,:))]);
                    %hold off;
                end
                figure;                
                h = scatter(dist,LI);
                set(h,'MarkerEdgeColor','k');
                set(h,'MarkerFaceColor','k');
                xlabel('Distance (microns)');
                ylabel('Linearity Index');         
                maxLI = max(abs(LI));
                axis([min(dist)*0.9, max(dist)*1.1, -maxLI*1.1, maxLI*1.1]);
                title(num2str(self.uniqueVals(i)));
            end
            drawnow;
            pause;
        end
        
    end
end
