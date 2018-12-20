classdef epochTreeGUI < handle
    
    properties
        epochTree;
        
        fig;
        xDivLeft = .25;
        xDivRight = .1;
        yDiv = .2;
        
        treeBrowser = struct();
        databaseInteraction = struct();
        plottingCanvas = struct();
        analysisTools = struct();
        GUIControls = struct();
        
        isBusy=false;
    end
    
    properties(Hidden = true)
        title = 'Epoch Tree GUI';
        busyTitle = 'Epoch Tree GUI (busy...)';
    end
    
    methods
        function self = epochTreeGUI(epochTree)
            if nargin < 1
                return
            end
            
            self.epochTree = epochTree;
            
            % init (slowly)
            self.buildUIComponents;
            self.isBusy = true;
            self.initAnalysisTools;
            self.initTreeBrowser;
            self.isBusy = false;
        end
        
        function buildUIComponents(self)
            
            % clean figure
            if ~isempty(self.fig) && ishandle(self.fig)
                delete(self.treeBrowser.panel);
                delete(self.databaseInteraction.panel);
                delete(self.plottingCanvas.panel);
                delete(self.analysisTools.panel);
                delete(self.GUIControls.panel);
                clf(self.fig);
            else
                self.fig = figure;
            end
            set(self.fig, ...
                'Name',         self.title, ...
                'NumberTitle',  'off', ...
                'ToolBar',      'none', ...
                'DeleteFcn',    {@epochTreeGUI.figureDeleteCallback, self}, ...
                'HandleVisibility', 'on');
            
            % new panels and widgets
            self.buildTreeBrowserUI();
            self.buildDatabaseInteractionUI();
            self.buildPlottingCanvasUI();
            self.buildAnalysisToolsUI();
            self.buildGUIControlsUI();
        end
        
        function initTreeBrowser(self)
            
            if isfield(self.treeBrowser, 'graphTree') && isobject(self.treeBrowser.graphTree)
                delete(self.treeBrowser.graphTree);
            end
            
            % brand new graphicalTree object
            graphTree = graphicalTree(self.treeBrowser.treeAxes, 'epoch tree');
            graphTree.expandNodeButtonCallback = {@epochTreeGUI.updateViewCallback, self};
            graphTree.clickNodeButtonCallback = {@epochTreeGUI.updateViewCallback, self};
            graphTree.nodeBecameCheckedCallback = {@epochTreeGUI.setEpochNodeIsSelectedCallback, self};
            graphTree.checkNodeButtonCallback = {@epochTreeGUI.mouseCheckedNodeCallback, self};
            graphTree.draw;
            self.treeBrowser.graphTree = graphTree;
            
            % populate grahical tree with EpochTree and Epoch objects
            if isobject(self.epochTree)
                self.marryEpochNodesToWidgets(self.epochTree, graphTree.trunk);
            end
            self.refreshAllBrowserSelections;
        end
        
        function initAnalysisTools(self)
            self.analysisTools.selectedFunction = {};
            self.analysisTools.lastUsedFunction = {};
            set(self.analysisTools.view, 'Value', 1);
            self.populateFunctionMenuWithFunctions(get(self.analysisTools.view, 'String'));
        end
        
        function rebuildUIComponents(self)
            
            % save existing graphical tree object
            graphTree = self.treeBrowser.graphTree;
            
            % fresh figure
            self.buildUIComponents;
            self.isBusy = true;
            
            % rewire old graphicalTree to new axes
            self.treeBrowser.graphTree = graphTree;
            self.treeBrowser.graphTree.ax = self.treeBrowser.treeAxes;
            self.refreshAllBrowserSelections;
            
            % init analysis tools
            self.initAnalysisTools;
            
            self.isBusy = false;
        end
        
        function buildTreeBrowserUI(self)
            % main tree browser panel and axes:
            treeBrowser.panel = uipanel( ...
                'Title',    'tree browser', ...
                'Parent',   self.fig, ...
                'HandleVisibility', 'off', ...
                'Units',    'normalized', ...
                'Position', [0 self.yDiv self.xDivLeft 1-self.yDiv]);
            treeBrowser.treeAxes = axes( ...
                'Parent',	treeBrowser.panel, ...
                'HandleVisibility', 'on', ...
                'Units',    'normalized', ...
                'Position', [0 .05 1 .9]);
            
            % smaller widgets:
            treeBrowser.splitKeys = uicontrol( ...
                'Parent',   treeBrowser.panel, ...
                'Style',    'text', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   getEpochTreeSplitString(self.epochTree), ...
                'HorizontalAlignment', 'left', ...
                'Position', [0 .95 1 .05]);
            treeBrowser.invertTreeSelections = uicontrol( ...
                'Parent',   treeBrowser.panel, ...
                'Callback', {@epochTreeGUI.invertTreeSelectionsCallback, self}, ...
                'Style',    'pushbutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   '1 / tree', ...
                'HorizontalAlignment', 'left', ...
                'Position', [0 0 .2 .05]);
            treeBrowser.invertEpochSelections = uicontrol( ...
                'Parent',   treeBrowser.panel, ...
                'Callback', {@epochTreeGUI.invertEpochSelectionsCallback, self}, ...
                'Style',    'pushbutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   '1 / epoch', ...
                'HorizontalAlignment', 'left', ...
                'Position', [.2 0 .2 .05]);
            treeBrowser.pan = uicontrol( ...
                'Parent',   treeBrowser.panel, ...
                'Callback', {@epochTreeGUI.panTreeBrowserCallback, self}, ...
                'Style',    'togglebutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   'pan', ...
                'HorizontalAlignment', 'left', ...
                'Position', [.5 0 .2 .05]);
            treeBrowser.refresh = uicontrol( ...
                'Parent',   treeBrowser.panel, ...
                'Callback', {@epochTreeGUI.refreshTreeBrowserCallback, self}, ...
                'Style',    'pushbutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   'refresh', ...
                'HorizontalAlignment', 'left', ...
                'Position', [.8 0 .2 .05]);
            self.treeBrowser = treeBrowser;
        end
        
        function buildDatabaseInteractionUI(self)
            % main database interaction panel:
            databaseInteraction.panel = uipanel( ...
                'Title',    'db interaction', ...
                'Parent',   self.fig, ...
                'HandleVisibility', 'off', ...
                'Units',    'normalized', ...
                'Position', [0 0 self.xDivLeft self.yDiv]);
            
            % smaller widgets:
            databaseInteraction.existingTag = uicontrol( ...
                'Parent',   databaseInteraction.panel, ...
                'Callback', {@epochTreeGUI.existingTagCallback, self}, ...
                'Style',    'popupmenu', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   {'existing tags'}, ...
                'Position', [.1 .6 .7 .2]);
            databaseInteraction.refreshTag = uicontrol( ...
                'Parent',   databaseInteraction.panel, ...
                'Callback', {@epochTreeGUI.refreshTagCallback, self}, ...
                'Style',    'pushbutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   'refresh', ...
                'Position', [.8 .6 .2 .2]);
            databaseInteraction.newTag = uicontrol( ...
                'Parent',   databaseInteraction.panel, ...
                'Callback', {@epochTreeGUI.newTagCallback, self}, ...
                'Style',    'edit', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   'new tag', ...
                'HorizontalAlignment', 'left', ...
                'BackgroundColor', [1 1 1], ...
                'Position', [.1 .4 .7 .2]);
            databaseInteraction.addTag = uicontrol( ...
                'Parent',   databaseInteraction.panel, ...
                'Callback', {@epochTreeGUI.addTagCallback, self}, ...
                'Style',    'pushbutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   'add', ...
                'Position', [.8 .4 .2 .2]);
            databaseInteraction.removeTag = uicontrol( ...
                'Parent',   databaseInteraction.panel, ...
                'Callback', {@epochTreeGUI.removeTagCallback, self}, ...
                'Style',    'pushbutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   ' remove ', ...
                'Position', [.8 .2 .2 .2]);
            databaseInteraction.includeInAnalysis = uicontrol( ...
                'Parent',   databaseInteraction.panel, ...
                'Callback', {@epochTreeGUI.includeInAnalysisCallback, self, true}, ...
                'Style',    'pushbutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   'include', ...
                'Position', [.1 0 .3 .2]);
            databaseInteraction.excludeFromAnalysis = uicontrol( ...
                'Parent',   databaseInteraction.panel, ...
                'Callback', {@epochTreeGUI.includeInAnalysisCallback, self, false}, ...
                'Style',    'pushbutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   'exclude', ...
                'Position', [.4 0 .3 .2]);
            self.databaseInteraction = databaseInteraction;
        end
        
        function buildPlottingCanvasUI(self)
            % big panel for plotting or custom tools
            self.plottingCanvas.panel = uipanel( ...
                'Title',    'plotting canvas', ...
                'Parent',   self.fig, ...
                'HandleVisibility', 'on', ...
                'Units',    'normalized', ...
                'Position', [self.xDivLeft self.yDiv 1-self.xDivLeft 1-self.yDiv]);
        end
        
        function buildAnalysisToolsUI(self)
            % main analysis tools panel:
            analysisTools.panel = uipanel( ...
                'Title',    'analysis tools', ...
                'HandleVisibility', 'off', ...
                'Parent',   self.fig, ...
                'Units',    'normalized', ...
                'Position', [self.xDivLeft 0 1-self.xDivLeft-self.xDivRight self.yDiv]);
            
            % mutually exclusive buttons
            analysisTools.afvGroup = uibuttongroup( ...
                'Parent',   analysisTools.panel, ...
                'SelectionChangeFcn', {@epochTreeGUI.afvGroupCallback, self}, ...
                'Units',    'normalized', ...
                'Position', [.1 .6 .8 .3]);
            analysisTools.analysis = uicontrol( ...
                'Parent',   analysisTools.afvGroup, ...
                'Style',    'togglebutton', ...
                'Units',    'normalized', ...
                'UserData', 1, ...
                'FontSize', 10, ...
                'String',   'analysis', ...
                'Position', [0 0 .25 1]);
            analysisTools.filter = uicontrol( ...
                'Parent',   analysisTools.afvGroup, ...
                'Style',    'togglebutton', ...
                'Units',    'normalized', ...
                'UserData', 1, ...
                'FontSize', 10, ...
                'String',   'filter', ...
                'Position', [.25 0 .25 1]);
            analysisTools.view = uicontrol( ...
                'Parent',   analysisTools.afvGroup, ...
                'Style',    'togglebutton', ...
                'Units',    'normalized', ...
                'UserData', 1, ...
                'FontSize', 10, ...
                'String',   'view', ...
                'Position', [.5 0 .25 1]);
            analysisTools.noUpdate = uicontrol( ...
                'Parent',   analysisTools.afvGroup, ...
                'Style',    'togglebutton', ...
                'Units',    'normalized', ...
                'UserData', 1, ...
                'FontSize', 10, ...
                'String',   'no update', ...
                'Position', [.75 0 .25 1]);
            
            % other widgets:
            analysisTools.functionMenu = uicontrol( ...
                'Parent',   analysisTools.panel, ...
                'Callback', {@epochTreeGUI.functionMenuCallback, self}, ...
                'Style',    'popupmenu', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   {'pick a tool'}, ...
                'Position', [.2 .2 .4 .2]);
            analysisTools.applyTool = uicontrol( ...
                'Parent',   analysisTools.panel, ...
                'Callback', {@epochTreeGUI.applyToolCallback, self}, ...
                'Style',    'pushbutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   'apply', ...
                'Position', [.6 .2 .2 .2]);
            self.analysisTools = analysisTools;
        end
        
        function buildGUIControlsUI(self)
            % GUI-wide control buttons
            GUIControls.panel = uipanel( ...
                'Title',    'GUI', ...
                'HandleVisibility', 'off', ...
                'Parent',   self.fig, ...
                'Units',    'normalized', ...
                'Position', [1-self.xDivRight 0 self.xDivRight self.yDiv]);
            
            GUIControls.rebuild = uicontrol( ...
                'Parent',   GUIControls.panel, ...
                'Callback', {@epochTreeGUI.rebuildCallback, self}, ...
                'Style',    'pushbutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   'reload', ...
                'Position', [.1 .2 .8 .2]);
            GUIControls.saveTree = uicontrol( ...
                'Parent',   GUIControls.panel, ...
                'Callback', {@epochTreeGUI.saveTreeCallback, self}, ...
                'Style',    'pushbutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   'save tree', ...
                'Position', [.1 .4 .8 .2]);
            GUIControls.loadTree = uicontrol( ...
                'Parent',   GUIControls.panel, ...
                'Callback', {@epochTreeGUI.loadTreeCallback, self}, ...
                'Style',    'pushbutton', ...
                'Units',    'normalized', ...
                'FontSize', 10, ...
                'String',   'load tree', ...
                'Position', [.1 .6 .8 .2]);
            
            self.GUIControls = GUIControls;
        end
        
        function set.isBusy(self, isBusy)
            self.isBusy = isBusy;
            if isBusy
                set(self.fig, 'Name', self.busyTitle);
            else
                set(self.fig, 'Name', self.title);
            end
            drawnow;
        end
        
        function marryEpochNodesToWidgets(self, epochNode, browserNode)
            browserNode.userData = epochNode;
            if ~isfield(epochNode.custom, 'isSelected')
                epochNode.custom.isSelected = false;
            end
            epochNode.custom.isCapsule = false;
            
            if epochNode.isLeaf
                % base case: new node for each leaf epoch
                import auimodel.*;
                iter = epochNode.epochList.iterator('cell');
                while iter.hasNext
                    ep = iter.nextValue;
                    if isempty(ep.isSelected)
                        ep.isSelected = true;
                    end
                    
                    % encapsulate epoch in an EpochTree node
                    %   which is not really part of the main EpochTree
                    capsuleList = EpochList;
                    capsuleList.append(ep);
                    capsuleList.populateStreamNames;
                    
                    capsuleNode = EpochTree.build(capsuleList, {});
                    capsuleNode.custom.isCapsule = true;
                    capsuleNode.parent = epochNode;
                    
                    capsuleName = sprintf('%3d: %s', iter.index-1, datestr(ep.startDate, 31));
                    capsuleWidget = browserNode.tree.newNode(browserNode, capsuleName, false);
                    capsuleWidget.userData = capsuleNode;
                end
                
                % do not propagate checking to epochs
                browserNode.recursiveCheck = false;
                
            else
                % recur: new browserNode for each child node
                for ii = 1:length(epochNode.children)
                    splitValue = epochNode.children{ii}.splitValue;
                    
                    if isobject(splitValue)
                        childName = splitValue.toString();
                    else
                        childName = num2str(epochNode.children{ii}.splitValue);
                    end
                    
                    childWidget = browserNode.tree.newNode(browserNode, childName);
                    self.marryEpochNodesToWidgets(epochNode.children{ii}, childWidget);
                end
            end
        end
        
        function refreshAllBrowserSelections(self, startNode)
            self.isBusy = true;
            
            if nargin < 2
                % default to whole tree
                startNode = self.treeBrowser.graphTree.trunk;
            end
            
            % read the EpochTree selections into the browser tree
            self.readEpochTreeSelections(startNode);
            
            % let entire(?) tree do accounting and update appearance
            %   might be able to start at startNode, would need to
            %   propagate any change up the tree
            self.treeBrowser.graphTree.trunk.countCheckedDescendants;
            self.treeBrowser.graphTree.draw;
            
            self.isBusy = false;
        end
        
        function readEpochTreeSelections(self, browserNode)
            node = browserNode.userData;
            if isobject(node)
                
                % deal with encapsulated single epochs
                if node.custom.isCapsule
                    browserNode.isChecked = node.epochList.elements{1}.isSelected;
                else
                    browserNode.isChecked = node.custom.isSelected;
                end
            end
            
            % recur: set child selections
            for ii = 1:length(browserNode.children)
                self.readEpochTreeSelections(browserNode.children(ii));
            end
        end
        
        function invokeAnalysisTool(self)
            if isempty(self.analysisTools.selectedFunction)
                return
            end
            
            self.isBusy = true;
            
            % initialize the current analysis tool?
            doInit = ~isequal(self.analysisTools.selectedFunction, self.analysisTools.lastUsedFunction);
            if doInit
                delete(get(self.plottingCanvas.panel, 'Children'));
            end
            
            % invoke analysis/filter/view function
            [epochNode, browserNode] = self.getCurrentNode;
            feval(self.analysisTools.selectedFunction, ...
                epochNode, ...
                self.plottingCanvas.panel, ...
                doInit);
            self.analysisTools.lastUsedFunction = self.analysisTools.selectedFunction;
            
            % refresh selections on return from filter function
            if get(self.analysisTools.filter, 'Value')
                self.refreshAllBrowserSelections(browserNode);
            end
            
            self.isBusy = false;
        end
        
        function [epochNode, browserNode] = getCurrentNode(self)
            % try to find the "EpochTree" highlighted in the tree browser
            if isobject(self.treeBrowser.graphTree.currentNode)
                browserNode = self.treeBrowser.graphTree.currentNode;
            else
                browserNode = self.treeBrowser.graphTree.trunk;
            end
            epochNode = browserNode.userData;
        end
        
        function reinvokeView(self)
            % if a view function is current, reinvoke it
            if get(self.analysisTools.view, 'Value')
                self.invokeAnalysisTool;
            end
        end
        
        function populateTagMenuWithTags(self)
            
            self.isBusy = true;
            
            menuWidget = self.databaseInteraction.existingTag;
            tags = self.getSelectedKeywordsIntersection;
            if isempty(tags)
                set(menuWidget, ...
                    'String',   {'no common tags'}, ...
                    'Enable',   'inactive', ...
                    'Value',    1);
            else
                set(menuWidget, ...
                    'String',   tags, ...
                    'Enable',	'on', ...
                    'Value',    1);
            end
            
            self.isBusy = false;
        end
        
        function commonKeys = getSelectedKeywordsIntersection(self)
            commonKeys = {};
            
            sel = getTreeEpochs(self.getCurrentNode, true);
            if isempty(sel)
                return
            end
            
            iter = sel.iterator('cell');
            if iter.hasNext
                % start with some list of keys
                ep = iter.nextValue;
                commonKeys = ep.keywords.elements;
                while iter.hasNext && ~isempty(commonKeys)
                    % pare down the list
                    ep = iter.nextValue;
                    commonKeys = intersect(commonKeys, ep.keywords.elements);
                end
            end
        end
        
        function populateFunctionMenuWithFunctions(self, subDir)
            global ANALYSIS_FILTER_VIEW_FOLDER
            funcFolder = fullfile(ANALYSIS_FILTER_VIEW_FOLDER, subDir);
            funcFiles = getMFiles(funcFolder, false);
            funcNames = cell(1, length(funcFiles)+1);
            funcNames{1} = '---pick a function---';
            for ii = 1:length(funcFiles)
                [pat, funcNames{ii+1}] = fileparts(funcFiles{ii});
            end
            set(self.analysisTools.functionMenu, 'String', funcNames);
        end
        
        function delete(self)
            % try to close the GUI figure
            if ~isempty(self.fig) && ishandle(self.fig) && strcmp(get(self.fig, 'BeingDeleted'), 'off')
                close(self.fig);
            end
        end
    end
    
    methods(Static)
        % tree browser callbacks:
        function setEpochNodeIsSelectedCallback(browserNode, gui)

            node = browserNode.userData;
            if isobject(node)
                if node.custom.isCapsule
                    % set selection to encapsulated Epoch
                    epoch = node.epochList.elements{1};
                    epoch.isSelected = browserNode.isChecked;
                end
                % set selction to EpochTree node (redundant for capsule)
                node.custom.isSelected = browserNode.isChecked;
            end
        end
        function mouseCheckedNodeCallback(browserNode, event, gui)
            gui.treeBrowser.graphTree.redrawChecks;
            gui.reinvokeView;
        end
        function updateViewCallback(browserNode, event, gui)
            gui.reinvokeView;
        end
        function invertTreeSelectionsCallback(widget, event, gui)
            % invoke filter function directly, on current node.
            % don't behave like menu-selected analysis tools function
            if exist('invertTreeSelections', 'file')
                self.isBusy = true;
                invertTreeSelections(gui.getCurrentNode);
                gui.refreshAllBrowserSelections;
                self.isBusy = false;
            end
        end
        function invertEpochSelectionsCallback(widget, event, gui)
            % invoke filter function directly, on current node.
            % don't behave like menu-selected analysis tools function
            if exist('invertEpochSelections', 'file')
                self.isBusy = true;
                invertEpochSelections(gui.getCurrentNode);
                gui.refreshAllBrowserSelections;
                self.isBusy = false;
            end
        end

        function panTreeBrowserCallback(widget, event, gui)
            p = pan(gui.fig);
            if get(widget, 'Value')
                % STUPID, axes will not pan when HandleVisibility=off
                %   even for the pan button in builtin figure toolbar
                
                set(gui.treeBrowser.treeAxes, 'HandleVisibility', 'on');
                set(p, 'Enable', 'on');
                setAllowAxesPan(p, gui.treeBrowser.treeAxes, true);
            else
                set(gui.treeBrowser.treeAxes, 'HandleVisibility', 'off');
                set(p, 'Enable', 'off');
                setAllowAxesPan(p, gui.treeBrowser.treeAxes, false);
            end
        end
        function refreshTreeBrowserCallback(widget, event, gui)
            if isfield(gui.treeBrowser, 'graphTree')
                gui.refreshAllBrowserSelections;
            end
        end
        
        % database interaction callbacks:
        function addTagCallback(widget, event, gui)
            gui.isBusy = true;
            tag = get(gui.databaseInteraction.newTag, 'String');
            sel = getTreeEpochs(gui.getCurrentNode, true);
            sel.addKeywordTag(tag);
            gui.populateTagMenuWithTags;
            gui.isBusy = false;
        end
        function removeTagCallback(widget, event, gui)
            gui.isBusy = true;
            tag = get(gui.databaseInteraction.newTag, 'String');
            sel = getTreeEpochs(gui.getCurrentNode, true);
            sel.removeKeywordTag(tag);
            gui.populateTagMenuWithTags;
            gui.isBusy = false;
        end
        function newTagCallback(widget, event, gui)
            % test validity?
        end
        function existingTagCallback(widget, event, gui)
            tags = get(widget, 'String');
            if ~isempty(tags)
                selectedTag = tags{get(widget, 'Value')};
                set(gui.databaseInteraction.newTag, 'String', selectedTag);
            end
        end
        function refreshTagCallback(widget, event, gui)
            gui.populateTagMenuWithTags;
        end
        function includeInAnalysisCallback(widget, event, gui, isIncluded)
            gui.isBusy = true;
            sel = getTreeEpochs(gui.getCurrentNode, true);
            iter = sel.iterator('cell');
            while iter.hasNext
                iter.nextValue.includeInAnalysis = isIncluded;
            end
            gui.isBusy = false;
        end
        
        % analysis tools callbacks:
        function afvGroupCallback(widget, event, gui)
            
            % remember last function used
            set(event.OldValue, 'UserData', get(gui.analysisTools.functionMenu, 'Value'));
            
            subdir = get(event.NewValue, 'String');
            if strcmp(subdir, 'no update')
                % don't do much
                set(gui.analysisTools.functionMenu, 'Enable', 'off');
            else
                % populate function menu with selected function type,
                %   recall last function of this type used
                gui.populateFunctionMenuWithFunctions(subdir);
                set(gui.analysisTools.functionMenu, ...
                    'Enable', 'on', ...
                    'Value', get(event.NewValue, 'UserData'));
                epochTreeGUI.functionMenuCallback(gui.analysisTools.functionMenu, [], gui);
            end
        end
        function functionMenuCallback(widget, event, gui)
            functionNames = get(widget, 'String');
            functionIndex = get(widget, 'Value');
            if ~isempty(functionNames) && functionIndex > 1
                gui.analysisTools.selectedFunction = str2func(functionNames{functionIndex});
            end
        end
        function applyToolCallback(widget, event, gui)
            gui.invokeAnalysisTool;
        end
        
        % GUI controls callbacks:
        function rebuildCallback(widget, event, gui)
            gui.rebuildUIComponents;
        end
        function saveTreeCallback(widget, event, gui)
            if isobject(gui.epochTree)
                [n, p] = uiputfile({'*.mat';'*.*'}, 'Save EpochTree as');
                if ischar(n)
                    gui.isBusy = true;
                    gui.epochTree.saveTree(fullfile(p,n));
                    gui.isBusy = false;
                end
            end
        end
        function loadTreeCallback(widget, event, gui)
            
            [n, p] = uigetfile({'*.mat';'*.*'}, 'Load EpochTree from file');
            if ischar(n)
                
                gui.isBusy = true;
                
                % forget old tree
                if isobject(gui.epochTree) && isfield(gui.treeBrowser, 'graphTree')
                    % destroy handle references to reduce closing time
                    gui.treeBrowser.graphTree.forgetExternalReferences;
                    gui.treeBrowser.graphTree.forgetInternalReferences;
                    delete(gui.treeBrowser.graphTree);
                end
                
                % get new tree
                import auimodel.EpochTree;
                gui.epochTree = EpochTree.loadTree(fullfile(p,n));
                
                % init (slowly)
                gui.buildUIComponents;
                gui.isBusy = true;
                gui.initAnalysisTools;
                gui.initTreeBrowser;
                gui.isBusy = false;
            end
        end
        
        % when the GUI figure closes, try to delete this object
        function figureDeleteCallback(fig, event, gui)
            gui.isBusy = true;
            if isobject(gui)
                if isfield(gui.treeBrowser, 'graphTree')
                    % destroy handle references to reduce closing time
                    gui.treeBrowser.graphTree.forgetExternalReferences;
                    gui.treeBrowser.graphTree.forgetInternalReferences;
                end
                delete(gui);
            end
        end
    end
end