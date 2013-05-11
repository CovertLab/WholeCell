%CellState
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef CellState < handle
    %constants
    properties
        wholeCellModelID  %ID of the instance of this class
        name              %name of the instance of this class
    end
    
    %constants
    properties
        parameterNames    %names of properties that are parameters
        parameterIndexs   %names of properties that are indices into parameters
    end
    
    %helper objects
    properties
        randStream        %random stream
    end
    
    %options
    properties
        verbosity = 0;    %0 = no output, 5 = maximum output
        seed = [];        %set to any number for reproducible random streams
    end
    
    %constructor
    methods
        function this = CellState(wholeCellModelID, name)
            this.wholeCellModelID = wholeCellModelID;
            this.name = name;
            
            this.constructRandStream();
        end
    end
    
    %random numbers
    methods
        function constructRandStream(this)
            this.randStream = edu.stanford.covert.util.RandStream('mcg16807');
        end
        
        function seedRandStream(this)
            if isempty(this.seed)
                this.seed = round(mod(now, 1) * 1e7);
            end
            this.randStream.reset(this.seed);
        end
    end
        
    methods
        %Build object graph by storing references to other state objects and
        %processes
        function storeObjectReferences(this, simulation) %#ok<MANU,INUSD>
        end
    end
    
    methods
        function initializeConstants(this, knowledgeBase, simulation)
            %copy options from simulation
            this.verbosity   = simulation.verbosity;
            this.seed        = simulation.seed;
            
            %set values of parameter properties
            %store annotation of parameter properties
            state = findobj(knowledgeBase.states,'wholeCellModelID',this.wholeCellModelID);
            this.parameterNames = {state.parameters.name}';
            this.parameterIndexs = {state.parameters.index}';
            for i = 1:length(state.parameters)
                p = state.parameters(i);
                if ~isempty(p.index)
                    this.(p.name)(this.(p.index),1) = p.defaultValue;
                else
                    this.(p.name) = p.defaultValue;
                end
            end
        end
    end
    
    %get/set methods of input/output properties. Creates structs with
    %fields equal to that of properties of the process.
    methods
        %options
        function value = getOptions(this)
            value = this.getAsStruct(this.optionNames);
        end
        
        function setOptions(this, value)
            this.setFromStruct(value, this.optionNames);
        end
        
        %parameters
        function value = getParameters(this)
            value = this.getAsStruct(this.parameterNames);
        end
        
        function setParameters(this, value)
            this.setFromStruct(value, this.parameterNames);
        end
        
        %fixed constants
        function value = getFixedConstants(this)
            value = this.getAsStruct(this.fixedConstantNames);
        end
        
        function setFixedConstants(this, value)
            this.setFromStruct(value, this.fixedConstantNames);
        end
        
        %fit constants
        function value = getFittedConstants(this)
            value = this.getAsStruct(this.fittedConstantNames);
        end
        
        function setFittedConstants(this, value)
            this.setFromStruct(value, this.fittedConstantNames);
        end
        
        %time courses
        function value = getTimeCourses(this)
            value = this.getAsStruct(this.stateNames);
        end
        
        function setTimeCourses(this, value)
            this.setFromStruct(value, this.stateNames);
        end
    end
    
    methods (Access = private)
        %gets process property values in a struct
        function value = getAsStruct(this, fields)
            value = struct;
            validNames = fieldnames(this);
            for i = 1:length(fields)
                if ismember(fields{i}, validNames)
                    value.(fields{i}) = this.(fields{i});
                end
            end
        end
        
        %sets process property values according to struct values
        function setFromStruct(this, value, fields)
            validNames = fieldnames(this);
            for i = 1:length(fields)
                if isfield(value, fields{i}) && ismember(fields{i}, validNames)
                    try %#ok<TRYNC>
                        this.(fields{i}) = value.(fields{i});
                    end
                end
            end
        end
    end
    
    %printing
    methods
        function disp(this)
            metadata = metaclass(this);
            fprintf('%s %s CellState object with ', this.name, metadata.Name);
            if ~isempty(this.stateNames)
                fprintf('state properties:\n - %s', strjoin(sprintf('\n - '),this.stateNames{:}));
            else
                fprintf('no state properties');
            end
            fprintf('\n\n');
        end
    end
end
