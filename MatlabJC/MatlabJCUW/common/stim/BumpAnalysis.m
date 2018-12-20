classdef BumpAnalysis < AnalysisSuper
    properties
        epochNums
        threshold % could be useful for choosing how to separate the singles
        s
    end
    
    methods

        % This analysis works on a cell's data
        % equivalently, on a tree, can use nodeEpochsMatrix on a cell tree
        % node to pull out data and then run the same analysis
        function self = BumpAnalysis(params)
            self.initAnalysisSuper(params) ; 
 
         % establish range of epochs
            self.epochNums = str2num(params.epochNums) ;
            
         % choose the right epochs
            self.bumpEpochInfo(self.epochNums);
            
         % gets r and epochs info
            self.EpochDataPuller(self.epochNums,0);
            
         % process r
            self.EpochData = lowPassIdealFilter(self.EpochData, 5);

            figure
            chartRecord(self.EpochData,self.s);
            drawnow
            pause
            
            figure
            bumps(self.EpochData,self.s,1);
            
        end
        
        function bumpEpochInfo(self,epochNums,channel)
            
            for e = 1:length(epochNums)
                [stims(e),error] = ITCGetStmAmp(epochNums(e), 0, 0, self.fp);
                [output(e),error] = ITCGetOutputChan(epochNums(e),0,self.fp);
            end
            
            % most common stms are the bumps
            stms = unique(stims);
            stmscnt = stms;
            for i = 1:length(stms)
                stmscnt(i) = sum(stims==stms(i));
            end
            % split the stms into two groups, assuming one much set much
            % outweighs the other
            clusters = clusterdata(stmscnt',2);
            clustcnt = unique(clusters);
            for c = 1:length(clustcnt)
                ecnt(c) = sum(stmscnt(clusters==clustcnt(c)));
            end
            bumpCluster = clustcnt(ecnt==max(ecnt));
            trashStims = stms(clusters~=bumpCluster);
            bumpEpochs = ones(size(epochNums));
            for i = 1:length(trashStims)
                bumpEpochs = logical(bumpEpochs.*(stims~=trashStims(i)));
            end
            self.epochNums = epochNums(bumpEpochs);
            self.s.e = self.epochNums;
            self.s.stims = stims(bumpEpochs);
            self.MakeTime(self.s.e(1));
            self.s.time = self.Time;
        end
        
    end
end