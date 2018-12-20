%%
% general setup stuff
clear

% define plot color sequence, axis fonts
PlotColors = 'bgrkymcbgrkymcbgrkymcbgrkymc';
set(0, 'DefaultAxesFontName','Helvetica')
set(0, 'DefaultAxesFontSize', 9)
colormap([0 0 0])
scrsz = get(0, 'ScreenSize');


%%
% database communication
import auimodel.*
import vuidocument.*

%%
% Ovation-index
cd ~/Analysis/Ovation-index/


%%  WT_DIMS
% test = 0;
% if test 
%     load wt_ames__dims_test
%     expName = 'wt_ames__dims_test';
%     list = EpochList(OvationExport);
% 
%     tree = EpochTree(list, {...
%         'protocolSettings.acquirino:cellBasename',...
%         'protocolSettings.stimuli:Green_LED:amp',...
%     });
% else
%     load wt_ames__dims
%     expName = 'wt_ames__dims';
%     list1 = EpochList(OvationExport);
%     
%     tree = EpochTree(list1, {...
%         'protocolSettings.acquirino:cellBasename',...
%         'protocolSettings.stimuli:Green_LED:amp',...
%     });
%     
%     load wt_ames__dims_thuys;
%     thuyslist = EpochList(OvationExport);
%     
%     thuystree = EpochTree(thuyslist, {...
%         'protocolSettings.acquirino:cellBasename',...
%         'protocolSettings.stimuli:Red_LED:amp',...
%     });
%     tree.rootNode.children(end+1:end+length(thuystree.rootNode.children)) = thuystree.rootNode.children(:);
% end
% 
% cd wt_ames

%tree.visualize;

%%  WT_L15_LOCKES_DIMS
% test = 0;
% 
% 
% load wt_l15_lockes_dims
% expName = 'wt_l15_lockes_dims';
% 
% list = EpochList(OvationExport);
% 
% tree = EpochTree(list, {...
%     'protocolSettings.acquirino:cellBasename',...
%     'protocolSettings.stimuli:Green_LED:amp',...
%     });
% 
% cd wt_l15_lockes

%%  ARRHET_AMES_DIMS
test = 0;

load ArrHet/arrhet_ames_dims
expName = 'arrhet_ames_dims';

list = EpochList(OvationExport);

tree = EpochTree(list, {...
    'protocolSettings.acquirino:cellBasename',...
    'protocolSettings.stimuli:Red_LED:amp',...
    });

cd wt_l15_lockes

%% Clean up
clearvars -except OvationExport list tree expName PlotColors

% standard time points; most epochs have 10ms stimuli
base = 100;
endp = 1010;

% Start with a cell
c = 1;

% Silence plots
plotindcell = 1; % plot individual cells

% Figure 1 is the profile of the cell
%  ___ ___ ___
% |___|___|___|
% |___|___|___|

if plotindcell
    figure(1)
end

for cell = tree.rootNode.children
    
    clear mu time phi peak ttpk intt
    alen = max(3,length(cell{1}.children));
    a = 1;
    for leaf = cell{1}.children

        clear l
        l = getLeafTraits(leaf);
        
        % Calibrate the LED for this condition
        % dphi
        LED = l.stin(1:end-4);
        switch LED
            case 'Red'
                fluxs = sprintf('dphi = rod%sLED%sFlux(l.dvec);',LED,num2str(2));
                eval(fluxs);
                cola = 1;
            case 'Green'
                fluxs = sprintf('dphi = rod%sLED%sFlux(l.dvec);',LED,num2str(2));
                eval(fluxs);
                if strfind(l.file,'072607')
                    fluxs = sprintf('dphi = rod%sLED%sFlux(l.dvec);',LED,num2str(1));
                    eval(fluxs);
                end
                cola = 1;
        end
        phi(a) = dphi*l.ampl*l.stmp*l.samp;
        %---

        % Process data, then shorten.
        data = -leaf{1}.epochList.responsesByStreamName(l.resn);
        data = baselineSubtract(data,l.prep);
        data = lowPassIdealFilter(data,20,l.samp);
        
        time = ((1:base+endp)-base)*l.samp;
        
        temp = data; clear data,
        data = temp(:,l.prep-base+1:l.prep+endp);

        mu(a,:) = mean(data);

        [mumx,mmxi] = max(mu(a,:));
        peak(a) = mumx; % only a single value
        ttpk(a) = (mmxi-base-l.stmp/2)*l.samp;
        intt(a) = max(cumsum(mu(a,:)))/peak(a)*l.samp;

        if plotindcell
            plotDimAmplitudeData(alen,a,time,data,mu,phi,peak,ttpk,intt,l.samp)
        end
        a = a+1;
    
    end

    if plotindcell
        %-- completing and saving cell figure(1)
        completeAndSaveDimFigure(alen,phi,peak,intt,ttpk,time,mu,cell{1},l.samp)
        figure(1), clf
        
    end

    tempmu = mu;
    for p = 1:length(phi)
        tempmu(p,:) = mu(p,:)/phi(p);
    end
    cellmu(c,:) = mean(tempmu); clear tempmu;
    
    experiment{c,:} = {l.file,time,mu,phi,peak,ttpk,intt};

    c = c+1;
end

plotDimsAcrossCells(cellmu,l.prep,l.stmp,l.samp,time,expName)

cd ~/Analysis/Ovation-index/wt
data_header = {'file','time','mu','phi','peak','ttpk','intt'};
eval(sprintf('save %s_data data_header experiment',expName));

%%
% [emax,maxi] = max(data,[],2);
        %peak(aind) = mean(emax)/phi(aind); % biased towards high values
        %ttpk(aind) = mean((maxi-prep-stmp/2)*samp);
        %intt(aind) = mean(max(cumsum(shortdata,2),[],2)./emax*samp);

