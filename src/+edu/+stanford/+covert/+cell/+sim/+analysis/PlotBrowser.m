%PlotBrowser
% Creates GUI window which allows the user to quickly view plots of
% simulation. User can either select plots from list box, or scroll through
% plots using the previous and next page buttons.
%
%Uses meta class to find names of all public static  plot methods in
%analysis classes. Their signature should be:
%   function plot*(sim, axes, time, compartments)
%   end
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/6/2011
classdef PlotBrowser < handle
    properties (SetAccess = protected)
        simulation
        
        plotNames
        plotFunctions
        
        figure
        buttonsPanel
        compartmentListBoxTitle
        compartmentListBox
        plotListBoxTitle
        plotListBox
        plotButton
        prevButton
        nextButton
        saveButton
        printButton
        plotPanel
    end
    
    properties (Constant = true)
        plotClasses = {
            'Chromosome'
            'MetabolicReaction'
            'Transcription'
            'Translation'
            }
    end
    
    methods
        function this = PlotBrowser(sim, selectedPlotNames, selectedCompartmentNames)
            this = this@handle;
            this.simulation = sim;
            
            %find plots
            this.discoverPlots();
            
            %open window
            this.open();
            
            %layout window fullscreen
            screenSize = get(0, 'ScreenSize');
            screenSize(2) = 30;
            screenSize(4) = screenSize(4) - screenSize(2);
            screenSize([3 4]) = max(1, screenSize([3 4]));
            set(this.figure, 'OuterPosition', screenSize);
            
            %plot
            if nargin < 2
                selectedPlotNames = this.plotNames{1};
            end
            if nargin < 3
                selectedCompartmentNames = 'Sum';
            end
            this.selectPlots(selectedPlotNames, selectedCompartmentNames);
        end
    end
    
    methods
        function selectPlots(this, selectedPlotNames, selectedCompartmentNames)
            %select plots
            [~, selectedPlots] = ismember(selectedPlotNames, this.plotNames);
            set(this.plotListBox, 'Value', selectedPlots);
            
            %select compartments
            if nargin < 3
                selectedCompartmentNames = 'Sum';
            end
            [~, selectedCompartments] = ismember(selectedCompartmentNames, get(this.compartmentListBox, 'String'));
            set(this.compartmentListBox, 'Value', selectedCompartments);
            
            %plot
            this.graphSelectedPlots();
        end
        
        function open(this)
            %figure
            this.figure = figure('Name', 'Whole Cell Simulation :: Plot Browser', ...
                'NumberTitle', 'off', 'ResizeFcn', @this.layoutWindow); %#ok<CPROP,PROP>
            
            %list box panel
            this.buttonsPanel = uipanel('Title', 'Configure Plots');
            
            %select compartment
            this.compartmentListBoxTitle = uicontrol(this.buttonsPanel, 'Style', 'text',...
                'String', 'Select Compartments', 'HorizontalAlignment', 'Center');
            this.compartmentListBox = uicontrol(this.buttonsPanel, 'Style', 'listbox',...
                'String', [{'Sum'; 'Mean'; 'Median'; 'Std Dev'}; this.simulation.compartment.names],...
                'Min', 1, 'Max', 4 + this.simulation.compartment.count, 'Value', 1);
            
            %plot list box
            this.plotListBoxTitle = uicontrol(this.buttonsPanel,'Style', 'text',...
                'String','Select Plots', 'HorizontalAlignment', 'Center');
            this.plotListBox = uicontrol(this.buttonsPanel, 'Style', 'listbox',...
                'String', this.plotNames, 'Min', 1, 'Max', length(this.plotNames), 'Value', 1);
            
            %plot button
            this.plotButton = uicontrol(this.buttonsPanel, 'Style', 'pushbutton',...
                'String', 'Plot', 'Callback', @this.graphSelectedPlots);
            
            %previous page button
            this.prevButton = uicontrol(this.buttonsPanel, 'Style', 'pushbutton',...
                'String', 'Prev', 'Callback', @this.graphPrevPlots);
            
            %next page button
            this.nextButton = uicontrol(this.buttonsPanel, 'Style', 'pushbutton',...
                'String', 'Next', 'Callback', @this.graphNextPlots);
            
            %save button
            this.saveButton = uicontrol(this.buttonsPanel, 'Style', 'pushbutton',...
                'String', 'Save', 'Callback', @this.save);
            
            %print button
            this.printButton = uicontrol(this.buttonsPanel, 'Style', 'pushbutton',...
                'String', 'Print', 'Callback', @this.print);
            
            %plot panel
            this.plotPanel = uipanel('Title','Plots');
        end
        
        function close(this)
            close(this.figure);
            
            this.figure = [];
            this.buttonsPanel = [];
            this.compartmentListBoxTitle = [];
            this.compartmentListBox = [];
            this.plotListBoxTitle = [];
            this.plotListBox = [];
            this.plotButton = [];
            this.prevButton = [];
            this.nextButton = [];
            this.saveButton = [];
            this.printButton = [];
            this.plotPanel = [];
        end
    end
    
    %helper methods
    methods (Access = protected)
        %Return list of names of plotting functions
        function discoverPlots(this)
            this.plotNames = cell(0, 1);
            this.plotFunctions = cell(0, 1);
            
            metaPackage = meta.package.fromName('edu.stanford.covert.cell.sim.analysis');
            for i = 1:numel(metaPackage.Classes)
                metaClass = metaPackage.Classes{i};
                
                %ignore test clases, handle clases
                if any(ismember(cellfun(@(class) class.Name, metaClass.SuperClasses, 'UniformOutput', false), {'TestCase', 'handle'}))
                    continue;
                end
                
                %ignore unsupported classes
                className = metaClass.Name(find(metaClass.Name == '.', 1, 'last') + 1:end);
                if ~ismember(className, this.plotClasses)
                    continue;
                end
                
                %find plotting methods
                for j = 1:numel(metaClass.Methods)
                    metaMethod = metaClass.Methods{j};
                    
                    %filter out methods to match required signature
                    if ...
                            length(metaMethod.Name) < 5 || ...
                            ~strcmp(metaMethod.Name(1:4), 'plot') || ...
                            ~metaMethod.Static || ...
                            ~strcmp(metaMethod.Access, 'public') || ...
                            length(metaMethod.InputNames) < 4
                        continue;
                    end
                    
                    this.plotNames = [this.plotNames; className ':' strrep(metaMethod.Name(5:end), '_', ' ')];
                    this.plotFunctions = [this.plotFunctions; metaClass.Name '.' metaMethod.Name];
                end
                
                %sort plots by class, name
                [~, order] = sort(this.plotNames);
                this.plotNames = this.plotNames(order);
                this.plotFunctions = this.plotFunctions(order);
            end
        end
        
        function graphSelectedPlots(this, ~, ~)
            sim = this.simulation;
            selectedPlots = sort(get(this.plotListBox, 'Value'));
            time = permute(sim.state('Time').values, [3 1 2]);
            compartments = sort(get(this.compartmentListBox, 'Value') - 4);
            
            rows = ceil(sqrt(length(selectedPlots)));
            cols = ceil(length(selectedPlots) / rows);
            for i = 1:length(selectedPlots)
                axes = subplot(rows, cols, i, 'Parent', this.plotPanel);
                cla(axes);
                legend(axes, 'off');
                set(axes, 'Box', 'on');
                
                feval(this.plotFunctions{selectedPlots(i)}, sim, axes, time, compartments);
                
                set(axes, 'fontSize', 6);
                title(axes, this.plotNames{selectedPlots(i)}, 'fontSize', 10);
            end
        end
    end
    
    %Callbacks
    methods (Access = protected)
        function layoutWindow(this, ~, ~)
            %figure size
            figureSize = get(this.figure, 'Position');
            figureWidth = figureSize(3);
            figureHeight = figureSize(4);
            
            %list box panel
            panelWidth = 300;
            set(this.buttonsPanel, 'Position', [0 0 panelWidth / figureWidth 1]);
            
            %select compartment
            compartmentListBoxWidth = panelWidth - 20;
            compartmentListBoxHeight = 100;
            compartmentListBoxSize = [10 120 compartmentListBoxWidth compartmentListBoxHeight];
            set(this.compartmentListBoxTitle, 'Position', [10 120+compartmentListBoxHeight compartmentListBoxWidth 20]);
            set(this.compartmentListBox, 'Position', compartmentListBoxSize);
            
            %plot list box
            plotListBoxWidth = compartmentListBoxWidth;
            plotListBoxHeight = figureHeight-140-compartmentListBoxHeight-35-10-20;
            plotListBoxSize = [10 120+compartmentListBoxHeight+10+35 plotListBoxWidth plotListBoxHeight];
            set(this.plotListBoxTitle, 'Position', [10 120+compartmentListBoxHeight+10+35+plotListBoxHeight compartmentListBoxWidth 20]);
            set(this.plotListBox, 'Position', plotListBoxSize);
            
            %plot button
            buttonWidth = panelWidth-20;
            buttonHeight = 20;
            buttonSize = [10 70 buttonWidth buttonHeight];
            set(this.plotButton, 'Position', buttonSize);
            
            %previous page button
            buttonWidth = (panelWidth-20-10)/2;
            buttonSize = [10 40 buttonWidth buttonHeight];
            set(this.prevButton, 'Position', buttonSize);
            
            %next page button
            buttonSize = [10+buttonWidth+10 40 buttonWidth buttonHeight];
            set(this.nextButton, 'Position', buttonSize);
            
            %save button
            buttonSize = [10 10 buttonWidth buttonHeight];
            set(this.saveButton, 'Position', buttonSize);
            
            %print button
            buttonSize = [10 + buttonWidth + 10 10 buttonWidth buttonHeight];
            set(this.printButton, 'Position', buttonSize);
            
            %plot panel
            set(this.plotPanel, 'Position', [panelWidth / figureWidth 0 1 - panelWidth/figureWidth 1]);
        end
        
        function graphPrevPlots(this, ~, ~)
            selectedPlots = get(this.plotListBox, 'Value');
            
            unselectedPlots = setdiff(1:length(this.plotNames), selectedPlots);
            
            %select next plots
            if ...
                    max(selectedPlots) - min(selectedPlots) + 1 == length(selectedPlots) || ...
                    max(unselectedPlots) - min(unselectedPlots) + 1 == length(unselectedPlots)
                selectedPlots = mod(selectedPlots - length(selectedPlots) - 1, length(this.plotNames)) + 1;
            else
                selectedPlots = (1:length(selectedPlots)) + min(selectedPlots) - 1;
            end
            set(this.plotListBox ,'Value', sort(selectedPlots));
            
            %graph
            this.graphSelectedPlots();
        end
        
        function graphNextPlots(this, ~, ~)
            selectedPlots = get(this.plotListBox, 'Value');
            
            unselectedPlots = setdiff(1:length(this.plotNames), selectedPlots);
            
            %select next plots
            if ...
                    max(selectedPlots) - min(selectedPlots) + 1 == length(selectedPlots) || ...
                    max(unselectedPlots) - min(unselectedPlots) + 1 == length(unselectedPlots)
                selectedPlots = mod(selectedPlots+length(selectedPlots)-1,length(this.plotNames))+1;
            else
                selectedPlots = (1:length(selectedPlots)) + min(selectedPlots) - 1;
            end
            set(this.plotListBox, 'Value', sort(selectedPlots));
            
            %graph
            this.graphSelectedPlots();
        end
        
        function save(this, ~, ~)
            filemenufcn(this.figure, 'FileSaveAs')
        end
        
        function print(this, ~, ~)
            printdlg(this.figure);
        end
    end
end