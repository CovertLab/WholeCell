% Creates GUI window which allows the user to quickly view plots of
% simulation. User can either select plots from list box, or scroll through
% plots using the previous and next page buttons.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/9/2011
classdef StateBrowser < handle
    properties (SetAccess = protected)
        simulationIdx
        
        selectedPlotNumber
        plotConfigurations
        plotData
        
        simulationMetadata
        states
        stateStrings
        stateSubsets
        timeMin
        timeMax
        timeStep
        
        timeFrom = 1;
        timeTo = 1e9;
        timeBy = 1;
        
        verbosity = 0;
    end
    
    %handles
    properties (Access = protected)
        figureHandle
        optionPanelHandle
        plotConfigurationSubpanelHandle
        plotConfigurationTabButtonPanelHandle
        plotConfigurationTabPanelHandle
        plotConfigurationTabHandles
        plotConfigurationPanelHandles
        dataDirectorySubpanelHandle
        dataDirectoryListboxHandle
        dataDirectorySelectButtonHandle
        simulationMetadataGridHandle
        timeSubpanelHandle
        timeFromLabelHandle
        timeToLabelHandle
        timeByLabelHandle
        timeFromHandle
        timeToHandle
        timeByHandle
        actionSubpanelHandle
        plotButtonHandle
        saveButtonHandle
        plotPanelHandle
        axesHandles
    end
    
    methods
        function this = StateBrowser(dataDirectory, plotConfigurations, timeFrom, timeTo, timeBy, verbosity)
            %capture simulation metadata
            this.simulationMetadata = edu.stanford.covert.cell.sim.util.SimulationDiskUtil.getSimulations();
            
            %option window
            this.open();
            
            %select simulation
            if nargin == 0
                dataDirectory = this.simulationMetadata(1).directory;
            end
            this.selectDirectory(dataDirectory);
            
            %load plot configuration
            this.selectedPlotNumber = 1;
            if nargin < 2
                this.plotConfigurations = this.createPlotConfiguration();
                this.createPlotConfigurationTabPanel();
            else
                this.plotConfigurations = plotConfigurations;
                for i = 1:numel(plotConfigurations)
                    [tfs, this.plotConfigurations(i).state] = ismember(this.plotConfigurations(i).state, this.stateStrings);
                    if ~all(tfs)
                        throw(MException('StateBrowser:error', 'Undefined state %d', i));
                    end
                    
                    if ~isempty(this.plotConfigurations(i).stateSubset)
                        [tfs, this.plotConfigurations(i).stateSubset] = ismember(this.plotConfigurations(i).stateSubset, this.stateSubsets{this.plotConfigurations(i).state});
                        if ~all(tfs)
                            throw(MException('StateBrowser:error', 'Undefined state subset %d', i));
                        end
                    end
                end
                this.createPlotConfigurationTabPanel();
            end
            
            %set time range
            if nargin >= 3
                this.timeFrom = timeFrom;
            end
            if nargin >= 4
                this.timeTo = timeTo;
            end
            if nargin >= 5
                this.timeBy = timeBy;
            end
            
            %plot
            if nargin >= 2
                this.plot();
            end
            
            %set verbosity
            if nargin >= 6
                this.verbosity = verbosity;
            end
        end
        
        function open(this)
            if ~isempty(this.figureHandle)
                return;
            end
            
            %% figure
            this.figureHandle = figure('Name', 'Whole Cell Simulation :: State Browser', ...
                'Position', [100 100 1200 800], ...
                'Units', 'pixels', ...
                'Color', [0.9412 0.9412 0.9412], ...
                'NumberTitle', 'off', ...
                'MenuBar', 'none', ...
                'Toolbar', 'none', ...
                'PaperPositionMode', 'auto', ...
                'ResizeFcn', @(hObject, eventData) this.layout());
            
            %% option panel
            this.optionPanelHandle = uipanel(...
                'Parent', this.figureHandle, ...
                'Title', [], ...
                'BorderType', 'none');
            
            %plot configure
            this.plotConfigurationSubpanelHandle = uipanel(...
                'Parent', this.optionPanelHandle, ...
                'Title', [], ...
                'BorderType', 'none');
            this.plotConfigurationTabButtonPanelHandle = uipanel(...
                'Parent', this.plotConfigurationSubpanelHandle, ...
                'Title', [], ...
                'BorderType', 'none');
            this.plotConfigurationTabPanelHandle = uipanel('Parent', this.plotConfigurationSubpanelHandle, 'Title', []);
            this.createPlotConfigurationTabPanel();
            
            %data directory
            this.dataDirectorySubpanelHandle = uipanel(...
                'Parent', this.optionPanelHandle, ...
                'Title', 'Simulation', ...
                'BorderType', 'etchedin');
            this.dataDirectoryListboxHandle = uicontrol(this.dataDirectorySubpanelHandle, ...
                'Style', 'popupmenu',...
                'String', {this.simulationMetadata.name}, ...
                'Value', 1, ...
                'HorizontalAlignment', 'Left', ...
                'Units', 'normalized');
            this.dataDirectorySelectButtonHandle = uicontrol(this.dataDirectorySubpanelHandle, ...
                'Style', 'pushbutton',...
                'String', 'Select', ...
                'HorizontalAlignment', 'Center', ...
                'Units', 'normalized', ...
                'Callback', @(hObject, eventData) this.selectDirectory());
            if ~isdeployed
                properties = [
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Simulation', 'DisplayName', 'Name')
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Simulation', 'DisplayName', 'Description')
                    PropertyGridField('integer', uint32(0), 'Type', PropertyType('uint32', 'scalar'), 'Description', '', 'Category', 'Simulation', 'DisplayName', 'Length (s)')
                    PropertyGridField('integer', uint32(0), 'Type', PropertyType('uint32', 'scalar'), 'Description', '', 'Category', 'Metadata',   'DisplayName', 'Revision')
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Metadata',   'DisplayName', 'Differences From Revision')
                    PropertyGridField('integer', uint32(0), 'Type', PropertyType('uint32', 'scalar'), 'Description', '', 'Category', 'Metadata',   'DisplayName', 'Knowledge Base WID')
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Metadata',   'DisplayName', 'Group ID')
                    PropertyGridField('integer', uint32(0), 'Type', PropertyType('uint32', 'scalar'), 'Description', '', 'Category', 'Metadata',   'DisplayName', 'Group Index')
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Metadata',   'DisplayName', 'Start Time')
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Metadata',   'DisplayName', 'End Time')
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Metadata',   'DisplayName', 'Output Directory')
                    PropertyGridField('integer', uint32(0), 'Type', PropertyType('uint32', 'scalar'), 'Description', '', 'Category', 'Metadata',   'DisplayName', 'Segment Size (s)')
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Reseacher',  'DisplayName', 'Name')
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Reseacher',  'DisplayName', 'Affiliation')
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Reseacher',  'DisplayName', 'Username')
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Reseacher',  'DisplayName', 'Email')
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Reseacher',  'DisplayName', 'Hostname')
                    PropertyGridField('string',  '',        'Type', PropertyType('char', 'row'),      'Description', '', 'Category', 'Reseacher',  'DisplayName', 'IP')
                    ];
                properties = properties.GetHierarchy();
                this.simulationMetadataGridHandle = PropertyGrid(this.dataDirectorySubpanelHandle, ...
                    'Properties', properties, ...
                    'Position', [0 0 0.5 1]);
            end
            
            %time
            this.timeSubpanelHandle = uipanel(...
                'Parent', this.optionPanelHandle, ...
                'Title', 'Time', ...
                'BorderType', 'etchedin');
            
            this.timeFromLabelHandle = uicontrol(...
                'Parent', this.timeSubpanelHandle, ...
                'Style', 'text', ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'right', ...
                'String', 'From'); %#ok<*NASGU>
            this.timeFromHandle = uicontrol(...
                'Parent', this.timeSubpanelHandle, ...
                'Style', 'edit', ...
                'Units', 'normalized', ...
                'String', this.timeFrom, ...
                'Callback', @(hObject, eventData) this.timeFromCallback());
            
            this.timeToLabelHandle = uicontrol(...
                'Parent', this.timeSubpanelHandle, ...
                'Style', 'text', ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'right', ...
                'String', 'To'); %#ok<*NASGU>
            this.timeToHandle = uicontrol(...
                'Parent', this.timeSubpanelHandle, ...
                'Style', 'edit', ...
                'Units', 'normalized', ...
                'String', this.timeTo, ...
                'Callback', @(hObject, eventData) this.timeToCallback());
            
            this.timeByLabelHandle = uicontrol(...
                'Parent', this.timeSubpanelHandle, ...
                'Style', 'text', ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'right', ...
                'String', 'By'); %#ok<*NASGU>
            this.timeByHandle = uicontrol(...
                'Parent', this.timeSubpanelHandle, ...
                'Style', 'edit', ...
                'Units', 'normalized', ...
                'String', this.timeBy, ...
                'Callback', @(hObject, eventData) this.timeByCallback());
            
            %action
            this.actionSubpanelHandle = uipanel(...
                'Parent', this.optionPanelHandle, ...
                'Title', [], ...
                'BorderType', 'none');
            this.plotButtonHandle = uicontrol(this.actionSubpanelHandle, ...
                'Style', 'pushbutton',...
                'String', 'Plot', ...
                'Units', 'normalized', ...
                'Callback', @(hObject, eventData) this.plot);
            this.saveButtonHandle = uicontrol(this.actionSubpanelHandle, ...
                'Style', 'pushbutton',...
                'String', 'Save', ...
                'Units', 'normalized', ...
                'Callback', @(hObject, eventData) this.selectFileToSave);
            
            %% plot panel
            this.plotPanelHandle = uipanel('Parent', this.figureHandle, 'Title', 'Plots');
            
            %% layout window
            this.layout();
        end
        
        function layout(this)
            figureSize = get(this.figureHandle, 'Position');
            figureWidth = figureSize(3);
            figureHeight = figureSize(4);
            
            optionPanelWidth = 300;
            plotPanelWidth = figureWidth - optionPanelWidth - 20;
            if ~isdeployed
                dataDirectorySubpanelHeight = 210;
                simulationMetadataGridHeight = 178;
            else
                dataDirectorySubpanelHeight = 40;
                simulationMetadataGridHeight = 0;
            end
            actionSubpanelHeight = 30;
            timeSubpanelHeight = 75;
            plotConfigurationSubpanelHeight = figureHeight - dataDirectorySubpanelHeight - actionSubpanelHeight - timeSubpanelHeight - 20-12;
            tabButtonHeight = 25;
            
            %main panels
            set(this.optionPanelHandle, 'Position', [0 0 optionPanelWidth/figureWidth 1]);
            set(this.plotPanelHandle, 'Position', [(optionPanelWidth+20)/figureWidth 4/figureHeight (plotPanelWidth-4)/figureWidth 1-8/figureHeight]);
            
            %option subpanels
            set(this.plotConfigurationSubpanelHandle, 'Position', ...
                [4/optionPanelWidth (actionSubpanelHeight+timeSubpanelHeight)/figureHeight+16/figureHeight 1-8/optionPanelWidth plotConfigurationSubpanelHeight/figureHeight]);
            set(this.dataDirectorySubpanelHandle, 'Position', ...
                [8/optionPanelWidth (actionSubpanelHeight+plotConfigurationSubpanelHeight+timeSubpanelHeight)/figureHeight+27/figureHeight 1-16/optionPanelWidth dataDirectorySubpanelHeight/figureHeight]);
            set(this.timeSubpanelHandle, 'Position', ...
                [8/optionPanelWidth actionSubpanelHeight/figureHeight+12/figureHeight 1-16/optionPanelWidth timeSubpanelHeight/figureHeight]);
            set(this.actionSubpanelHandle, 'Position', ...
                [4/optionPanelWidth 4/figureHeight 1-8/optionPanelWidth actionSubpanelHeight/figureHeight]);
            
            %option controls
            set(this.plotConfigurationTabButtonPanelHandle, 'Position', ...
                [4/optionPanelWidth 1-tabButtonHeight/plotConfigurationSubpanelHeight 1-8/optionPanelWidth tabButtonHeight/plotConfigurationSubpanelHeight]);
            set(this.plotConfigurationTabPanelHandle, 'Position', ...
                [4/optionPanelWidth 4/plotConfigurationSubpanelHeight 1-8/optionPanelWidth 1-(tabButtonHeight+2)/plotConfigurationSubpanelHeight]);
            
            set(this.dataDirectoryListboxHandle, 'Position', ...
                [0.01 (2+6+simulationMetadataGridHeight)/dataDirectorySubpanelHeight 0.7 1-(8+simulationMetadataGridHeight)/dataDirectorySubpanelHeight])
            set(this.dataDirectorySelectButtonHandle, 'Position', ...
                [0.73 (2+6+simulationMetadataGridHeight)/dataDirectorySubpanelHeight 0.26 1-(8+simulationMetadataGridHeight)/dataDirectorySubpanelHeight])
            if ~isdeployed
                set(this.simulationMetadataGridHandle, 'Position', ...
                    [0.01 2/dataDirectorySubpanelHeight 0.98 simulationMetadataGridHeight/dataDirectorySubpanelHeight])
            end
            
            set(this.timeFromLabelHandle, 'Position', ...
                [0 0.71 0.1 0.25]);
            set(this.timeToLabelHandle, 'Position', ...
                [0 0.39 0.1 0.25]);
            set(this.timeByLabelHandle, 'Position', ...
                [0 0.07 0.1 0.25]);
            set(this.timeFromHandle, 'Position', ...
                [0.11 0.71 0.88 0.25]);
            set(this.timeToHandle, 'Position', ...
                [0.11 0.39 0.88 0.25]);
            set(this.timeByHandle, 'Position', ...
                [0.11 0.07 0.88 0.25]);
            
            set(this.plotButtonHandle, 'Position', ...
                [0.08 6/actionSubpanelHeight 0.4 1-10/actionSubpanelHeight])
            set(this.saveButtonHandle, 'Position', ...
                [0.52 6/actionSubpanelHeight 0.4 1-10/actionSubpanelHeight])
            
            this.createPlotConfigurationTabPanel();
        end
        
        function close(this)
            delete(this.figureHandle);
            
            this.figureHandle = [];
            this.optionPanelHandle = [];
            this.plotConfigurationSubpanelHandle = [];
            this.plotConfigurationTabButtonPanelHandle = [];
            this.plotConfigurationTabPanelHandle = [];
            this.plotConfigurationTabHandles = [];
            this.plotConfigurationPanelHandles = [];
            this.dataDirectorySubpanelHandle = [];
            this.dataDirectorySelectButtonHandle = [];
            this.dataDirectoryListboxHandle = [];
            this.simulationMetadataGridHandle = [];
            this.timeSubpanelHandle = [];
            this.timeFromHandle = [];
            this.timeToHandle = [];
            this.timeByHandle = [];
            this.actionSubpanelHandle = [];
            this.plotButtonHandle = [];
            this.saveButtonHandle = [];
            this.plotPanelHandle = [];
            this.axesHandles = [];
        end
        
        %% layout window fullscreen
        function maximizeWindow(this)
            if ispc
                maxfig(this.figureHandle, 1);
            else
                warning('WholeCell:warning', 'Full screen not supported on this platform');
            end
        end
        
        function selectDirectory(this, dataDirectory)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            %set variables
            if nargin > 1
                [~, this.simulationIdx] = ismember(dataDirectory, {this.simulationMetadata.directory});
                set(this.dataDirectoryListboxHandle, 'Value', this.simulationIdx);
            else
                this.simulationIdx = get(this.dataDirectoryListboxHandle, 'Value');
            end
            
            %load simulation
            dataDirectory = this.simulationMetadata(this.simulationIdx).directory;
            this.states = DiskLogger.getAvailableStates(dataDirectory);
            this.stateStrings = cellfun(@(x, y) [x ' - ' y], this.states(:, 1), this.states(:, 2), 'UniformOutput', false);
            this.stateSubsets = cell(size(this.stateStrings));
            tmp = load('src_test/+edu/+stanford/+covert/+cell/+sim/fixtures/Simulation.mat');
            sim = tmp.fixture;
            for i = 1:numel(sim.states)
                s = sim.states{i};
                if isa(s, 'edu.stanford.covert.cell.sim.MoleculeCountState')
                    j = find(strcmp(this.states(:,1), s.wholeCellModelID(7:end)) & strcmp(this.states(:, 2), 'counts'));
                    this.stateSubsets{j} = s.wholeCellModelIDs; %#ok<FNDSB>
                end
            end
            clear tmp smp s;
            metadata = DiskLogger.loadMetadata(dataDirectory);
            options = DiskLogger.loadOptions(dataDirectory);
            this.timeMin = 1;
            this.timeMax = metadata.lengthSec;
            this.timeStep = options.stepSizeSec;
            
            %setup metadata panel
            if ~isdeployed
                properties = this.simulationMetadataGridHandle.Properties;
                md = this.simulationMetadata(this.simulationIdx);
                tmp = findobj(properties, 'Category', 'Simulation', '-and', 'DisplayName', 'Name'); tmp.Value = md.shortDescription;
                tmp = findobj(properties, 'Category', 'Simulation', '-and', 'DisplayName', 'Description'); tmp.Value = md.longDescription;
                tmp = findobj(properties, 'Category', 'Simulation', '-and', 'DisplayName', 'Length (s)'); tmp.Value = uint32(md.lengthSec);
                tmp = findobj(properties, 'Category', 'Metadata', '-and', 'DisplayName', 'Revision'); tmp.Value = uint32(md.revision);
                tmp = findobj(properties, 'Category', 'Metadata', '-and', 'DisplayName', 'Differences From Revision'); tmp.Value = md.differencesFromRevision;
                tmp = findobj(properties, 'Category', 'Metadata', '-and', 'DisplayName', 'Knowledge Base WID'); tmp.Value = uint32(md.knowledgeBaseWID);
                tmp = findobj(properties, 'Category', 'Metadata', '-and', 'DisplayName', 'Group ID'); tmp.Value = md.simGroup;
                tmp = findobj(properties, 'Category', 'Metadata', '-and', 'DisplayName', 'Group Index'); tmp.Value = uint32(md.simIdx);
                tmp = findobj(properties, 'Category', 'Metadata', '-and', 'DisplayName', 'Start Time'); tmp.Value = md.startTime;
                tmp = findobj(properties, 'Category', 'Metadata', '-and', 'DisplayName', 'End Time'); tmp.Value = md.endTime;
                tmp = findobj(properties, 'Category', 'Metadata', '-and', 'DisplayName', 'Output Directory'); tmp.Value = md.outputDirectory;
                tmp = findobj(properties, 'Category', 'Metadata', '-and', 'DisplayName', 'Segment Size (s)'); tmp.Value = uint32(md.segmentSizeStep);
                tmp = findobj(properties, 'Category', 'Reseacher', '-and', 'DisplayName', 'Name'); tmp.Value = [md.firstName ' ' md.lastName];
                tmp = findobj(properties, 'Category', 'Reseacher', '-and', 'DisplayName', 'Affiliation'); tmp.Value = md.affiliation;
                tmp = findobj(properties, 'Category', 'Reseacher', '-and', 'DisplayName', 'Username'); tmp.Value = md.userName;
                tmp = findobj(properties, 'Category', 'Reseacher', '-and', 'DisplayName', 'Email'); tmp.Value = md.email;
                tmp = findobj(properties, 'Category', 'Reseacher', '-and', 'DisplayName', 'Hostname'); tmp.Value = md.hostName;
                tmp = findobj(properties, 'Category', 'Reseacher', '-and', 'DisplayName', 'IP'); tmp.Value = md.ipAddress;
                for i = 1:numel(properties)
                    if ischar(properties(i).Value)
                        properties(i).Description = properties(i).Value;
                    else
                        properties(i).Description = num2str(properties(i).Value);
                    end
                end
                this.simulationMetadataGridHandle.Properties = properties;
                this.simulationMetadataGridHandle.DisableToolBar();
                set(this.simulationMetadataGridHandle.Control, 'BorderType', 'none')
            end
            
            this.timeFromCallback();
            this.timeToCallback();
            this.timeByCallback();
        end
        
        function config = createPlotConfiguration(this)
            this.selectedPlotNumber = 1;
            stateSubset = this.stateSubsets{1};
            if ~isempty(stateSubset)
                stateSubset = 1;
            end
            
            config = struct(...
                'state', 1, ...
                'stateSubset', stateSubset);
        end
        
        function createPlotConfigurationTabPanel(this)
            delete(this.plotConfigurationTabHandles);
            delete(this.plotConfigurationPanelHandles);
            
            this.plotConfigurationTabHandles = [];
            this.plotConfigurationPanelHandles = [];
            
            for i = 1:numel(this.plotConfigurations)
                [tabHandle, panelHandle] = this.createPlotConfigurationTab(i);
                this.plotConfigurationTabHandles = [this.plotConfigurationTabHandles; tabHandle];
                this.plotConfigurationPanelHandles = [this.plotConfigurationPanelHandles; panelHandle];
            end
            
            [tabHandle, panelHandle] = this.createPlotConfigurationNewTab();
            this.plotConfigurationTabHandles = [this.plotConfigurationTabHandles; tabHandle];
            this.plotConfigurationPanelHandles = [this.plotConfigurationPanelHandles; panelHandle];
            
            set(this.plotConfigurationTabHandles(this.selectedPlotNumber), 'BackgroundColor', [0.7 0.7 0.7]);
            set(this.plotConfigurationPanelHandles(this.selectedPlotNumber), 'Visible', 'On');
        end
        
        function [tabHandle, panelHandle] = createPlotConfigurationTab(this, plotNumber)
            tabW = min(0.15, 1/(1+numel(this.plotConfigurations)));
            tabHandle = uicontrol(...
                'Parent', this.plotConfigurationTabButtonPanelHandle, ...
                'style', 'pushbutton', ...
                'String', sprintf('Plot-%d', plotNumber), ...
                'BackgroundColor', [0.9412 0.9412 0.9412], ...
                'SelectionHighlight', 'off', ...
                'HorizontalAlignment', 'left', ...
                'Units', 'normalized',...
                'Position', [tabW * (plotNumber-1) 0 tabW 1], ...
                'Callback', @(hObject, eventData) this.selectPlotConfiguration(plotNumber));
            panelHandle = uipanel(...
                'Parent', this.plotConfigurationTabPanelHandle,...
                'Visible', 'off',...
                'Units', 'normalized',...
                'Position', [0 0 1 1], ...
                'BorderType', 'none');
            
            figureSize = get(this.figureHandle, 'Position');
            figureWidth = figureSize(3);
            figureHeight = figureSize(4);
            panelHeight = figureHeight - 30 - 30 - 4 - 25;
            panelWidth = 300 - 8;
            stateHeight = 20;
            subsetHeight = panelHeight - stateHeight - 44;
            
            uicontrol(...
                'Parent', panelHandle, ...
                'Style', 'text', ...
                'Units', 'normalized', ...
                'Position', [4/panelWidth 1-(stateHeight+12)/panelHeight 41/panelWidth stateHeight/panelHeight], ...
                'HorizontalAlignment', 'right', ...
                'String', 'State'); %#ok<*NASGU>
            stateHandle = uicontrol(...
                'Parent', panelHandle, ...
                'Style', 'popupmenu', ...
                'Units', 'normalized', ...
                'Position', [49/panelWidth 1-(stateHeight+8)/panelHeight 1-53/panelWidth 20/panelHeight], ...
                'String', this.stateStrings, ...
                'Value', this.plotConfigurations(plotNumber).state, ...
                'Callback', @stateCallback); %#ok<*NASGU>
            
            if ~isempty(this.stateSubsets{this.plotConfigurations(plotNumber).state})
                enable = 'on';
                nEl = numel(this.stateSubsets{this.plotConfigurations(plotNumber).state});
            else
                val2 = [];
                enable = 'off';
                nEl = 0;
            end
            uicontrol(...
                'Parent', panelHandle, ...
                'Style', 'text', ...
                'Units', 'normalized', ...
                'Position', [4/panelWidth 1-(stateHeight+8+20+10)/panelHeight 41/panelWidth 20/panelHeight], ...
                'HorizontalAlignment', 'right', ...
                'String', 'Subset'); %#ok<*NASGU>
            stateSubsetHandle = uicontrol(...
                'Parent', panelHandle, ...
                'Style', 'listbox', ...
                'Units', 'normalized', ...
                'Position', [49/panelWidth 1-(stateHeight+subsetHeight+36)/panelHeight 1-53/panelWidth subsetHeight/panelHeight], ...
                'String', this.stateSubsets{this.plotConfigurations(plotNumber).state}, ...
                'Value', this.plotConfigurations(plotNumber).stateSubset, ...
                'Min', 1, ...
                'Max', nEl, ...
                'Enable', enable, ...
                'Callback', @stateSubsetCallback);
            
            function stateCallback(~, ~)
                this.plotConfigurations(plotNumber).state = get(stateHandle, 'Value');
                this.plotConfigurations(plotNumber).stateSubset = 1;
                if ~isempty(this.stateSubsets{get(stateHandle, 'Value')})
                    enable = 'on';
                else
                    enable = 'off';
                end
                set(stateSubsetHandle, ...
                    'String', this.stateSubsets{get(stateHandle, 'Value')}, ...
                    'Value', this.plotConfigurations(plotNumber).stateSubset, ...
                    'Enable', enable, ...
                    'Min', 1, ...
                    'Max', numel(this.stateSubsets{get(stateHandle, 'Value')}));
            end
            
            function stateSubsetCallback(~, ~)
                this.plotConfigurations(plotNumber).stateSubset = get(stateSubsetHandle, 'Value');
            end
        end
        
        function [tabHandle, panelHandle] = createPlotConfigurationNewTab(this)
            tabW = min(0.15, 1/(1+numel(this.plotConfigurations)));
            tabHandle = uicontrol(...
                'Parent', this.plotConfigurationTabButtonPanelHandle, ...
                'style', 'pushbutton', ...
                'String', '+', ...
                'BackgroundColor', [0.9412 0.9412 0.9412], ...
                'SelectionHighlight', 'off', ...
                'HorizontalAlignment', 'center', ...
                'Units', 'normalized',...
                'Position', [tabW * numel(this.plotConfigurations) 0 min(tabW, 0.1) 1],...
                'Callback', @(hObject, eventData) this.addPlotConfigurationTab());
            panelHandle = uipanel(...
                'Parent', this.plotConfigurationTabPanelHandle,...
                'Visible', 'off',...
                'Units', 'normalized',...
                'Position', [0 0 1 1], ...
                'BorderType', 'none');
        end
        
        function selectPlotConfiguration(this, plotNumber)
            set(this.plotConfigurationTabHandles(this.selectedPlotNumber), ...
                'BackgroundColor', [0.9412 0.9412 0.9412]);
            set(this.plotConfigurationPanelHandles(this.selectedPlotNumber), ...
                'Visible', 'Off');
            
            set(this.plotConfigurationTabHandles(plotNumber), ...
                'BackgroundColor', [0.7 0.7 0.7]);
            set(this.plotConfigurationPanelHandles(plotNumber), ...
                'Visible', 'On');
            
            this.selectedPlotNumber = plotNumber;
        end
        
        function addPlotConfigurationTab(this)
            this.plotConfigurations = [
                this.plotConfigurations
                this.createPlotConfiguration];
            this.selectedPlotNumber = numel(this.plotConfigurations);
            this.createPlotConfigurationTabPanel();
        end
        
        function timeFromCallback(this)
            this.timeFrom = ...
                ceil((min(this.timeMax, max(this.timeMin, str2double(get(this.timeFromHandle, 'String')))) -this.timeMin) / this.timeBy) * this.timeBy + this.timeMin;
            this.timeTo = ...
                floor((min(this.timeMax, max(this.timeMin, max(this.timeFrom, str2double(get(this.timeToHandle, 'String'))))) -this.timeMin) / this.timeBy) * this.timeBy + this.timeMin;
            set(this.timeFromHandle, 'String', this.timeFrom);
            set(this.timeToHandle, 'String', this.timeTo);
        end
        
        function timeToCallback(this)
            this.timeTo = ...
                floor((min(this.timeMax, max(this.timeMin, str2double(get(this.timeToHandle, 'String')))) -this.timeMin) / this.timeBy) * this.timeBy + this.timeMin;
            this.timeFrom = ...
                ceil((min(this.timeMax, max(this.timeMin, min(this.timeTo, str2double(get(this.timeFromHandle, 'String'))))) -this.timeMin) / this.timeBy) * this.timeBy + this.timeMin;
            set(this.timeToHandle, 'String', this.timeTo);
            set(this.timeFromHandle, 'String', this.timeFrom);
        end
        
        function timeByCallback(this)
            this.timeBy = ...
                round(min(this.timeMax - this.timeMin,  max(this.timeStep, str2double(get(this.timeByHandle, 'String')))) / this.timeStep) * this.timeStep;
            set(this.timeByHandle, 'String', this.timeBy);
            this.timeFromCallback();
            this.timeToCallback();
        end
        
        function plot(this)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            nRows = ceil(sqrt(numel(this.plotConfigurations)));
            nCols = max(nRows, ceil(numel(this.plotConfigurations) / nRows));
            nRows = ceil(numel(this.plotConfigurations) / nCols);
            
            delete(this.axesHandles);
            this.axesHandles = [];
            
            this.plotData = repmat(struct('time', [], 'data', []), numel(this.plotConfigurations), 1);
            for i = 1:numel(this.plotConfigurations)
                axesHandle = subplot(nRows, nCols, i, 'Parent', this.plotPanelHandle);
                cla(axesHandle);
                this.axesHandles = [
                    this.axesHandles
                    axesHandle];
                
                state = [this.states(this.plotConfigurations(i).state, :) ':' ':'];
                labels = {};
                if ~isempty(this.stateSubsets{this.plotConfigurations(i).state})
                    state{3} = this.plotConfigurations(i).stateSubset;
                    labels = this.stateSubsets{this.plotConfigurations(i).state}(this.plotConfigurations(i).stateSubset);
                end
                if this.verbosity
                    fprintf('Loading data ... ');
                end
                stateVals = DiskLogger.load(...
                    this.simulationMetadata(this.simulationIdx).directory, ...
                    [{'Time' 'values' ':' ':'}; state], ...
                    this.timeFrom, ...
                    this.timeTo, ...
                    this.timeBy, ...
                    'extract');
                if this.verbosity
                    fprintf('done\n');
                end
                time = permute(stateVals.Time.values, [1 3 2]);
                data = permute(sum(stateVals.(state{1}).(state{2}), 2), [1 3 2]);
                
                this.plotData(i).time = time;
                this.plotData(i).data = data;
                
                if isa(data, 'edu.stanford.covert.util.SparseMat')
                    data = full(sum(data ~= 0, 1));
                end
                
                if time(end) > 2 * 3600
                    time = time / 3600;
                    timeLabel = 'Time (h)';
                elseif time(end) > 2 * 60
                    time = time / 60;
                    timeLabel = 'Time (m)';
                else
                    timeLabel = 'Time (s)';
                end
                plot(axesHandle, time, data);
                xlabel(axesHandle, timeLabel, 'FontSize', 12);
                ylabel(axesHandle, state(1:2), 'FontSize', 12);
                xlim(axesHandle, [min(time) max(time)]);
                if ~isempty(labels)
                    legend(axesHandle, labels, 'Interpreter', 'None');
                end
            end
        end
        
        function selectFileToSave(this)
            if isdeployed
                fileFormats = {
                    '*.mat', 'MATLAB workspace (*.mat)'
                    };
                defualtFile = 'plots.mat';
            else
                fileFormats = {
                    '*.eps', 'EPS Level 1 (*.ai)'
                    '*.jpg', 'JPEG image (*.jpg)'
                    '*.fig', 'MATLAB figure (*.fig)'
                    '*.mat', 'MATLAB workspace (*.mat)'
                    '*.pdf', 'Portable Document Format (*.pdf)'
                    '*.png', 'Portable Network Graphics (*.png)'
                    '*.tif', 'TIFF image, compressed (*.tif)'
                    };
                defualtFile = 'plots.pdf';
            end
            
            %prompt for file
            fileName = uiputfile(fileFormats, 'Save plots', defualtFile);
            
            %canceled
            if fileName == 0
                return;
            end
            
            %save
            this.save(fileName);
        end
        
        function save(this, fileName)
            if isdeployed && (length(fileName) < 4 || ~isequal(fileName(end-3:end), '.mat'))
                warning('WholeCell:warning', 'Figure printing not supported on this platform');
                return;
            end
            
            if this.verbosity
                fprintf('Saving ...');
            end
            if length(fileName) >= 4 && isequal(fileName(end-3:end), '.mat')
                configurations = this.plotConfigurations;
                data = this.plotData;
                save(fileName, 'data', 'configurations');
            else
                saveas(this.figureHandle, fileName);
            end
            if this.verbosity
                fprintf(' done\n');
            end
        end
    end
end