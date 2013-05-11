%Database logger. Implements simulation logger interface.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/10/2011
classdef DatabaseLogger < edu.stanford.covert.cell.sim.util.Logger
    %options
    properties (SetAccess = protected)
        metadata          %metadata
        downsampleStepSec %max state history to cache in memory during a run (s)
        database          %database instance
        simulationWID     %simulation database ID
    end
    
    %saved data
    properties (SetAccess = protected)
        log
        randStreamStates
    end
    
    %indices into simulation state
    properties (Access = protected)
        stateIndex_time
    end
    
    methods
        function this = DatabaseLogger(database, downsampleStepSec, metadata)
            if ~exist('metadata', 'var')
                metadata = struct;
            end
            
            this.database = database;
            this.downsampleStepSec = downsampleStepSec;
            this.metadata = metadata;
        end
    end
    
    methods
        function this = addMetadata(this, metadata)
            names = fieldnames(metadata);
            for i = 1:length(names)
                this.metadata.(names{i}) = metadata.(names{i});
            end
        end
        
        function this = initialize(this, sim)
            %% validate options
            %segment length
            validateattributes(sim.lengthSec / this.downsampleStepSec, {'numeric'}, {'real', 'nonnegative', 'integer'});
            
            %downsample step
            validateattributes(sim.lengthSec / this.downsampleStepSec, {'numeric'}, {'real', 'nonnegative', 'integer'});
            
            %metadata
            this.validateMetadata();
            
            %% indices
            this.stateIndex_time = sim.stateIndex('Time');
            
            %% metdata
            %start time
            this.metadata.startTime = datestr(clock, 'yyyy-mm-dd HH:MM:SS');
            this.metadata.endTime = [];
            this.metadata.lengthSec = [];
            
            %output directory
            this.metadata.outputDirectory = '';
            
            %segment step
            this.metadata.segmentSizeStep = this.downsampleStepSec;
            
            %down sample step
            this.metadata.downsampleStepSec = this.downsampleStepSec;
            
            %% allocate memory
            %cell state
            metaStates = this.getMetaStates(sim);
            this.log = this.allocateMemory(metaStates, sim.lengthSec / this.downsampleStepSec);
            
            %rand stream state
            this.randStreamStates = struct('simulation', [], 'states', struct(), 'processes', struct());
            tmp = sim.getRandStreamStates();
            this.randStreamStates.simulation = zeros(numel(tmp.simulation), sim.lengthSec / this.downsampleStepSec);
            for i = 1:numel(sim.states)
                o = sim.states{i};
                this.randStreamStates.states.(o.wholeCellModelID(7:end)) = zeros(numel(tmp.states.(o.wholeCellModelID(7:end))), sim.lengthSec / this.downsampleStepSec);
            end
            for i = 1:numel(sim.processes)
                o = sim.processes{i};
                this.randStreamStates.processes.(o.wholeCellModelID(9:end)) = zeros(numel(tmp.processes.(o.wholeCellModelID(9:end))), sim.lengthSec / this.downsampleStepSec);
            end
        end
        
        function this = append(this, sim)
            t = sim.states{this.stateIndex_time};
            time = t.values;
            i = time / this.downsampleStepSec;
            if mod(i, 1) ~= 0
                return;
            end
            
            this.copyFromState(sim, i);
            this.copyRandStreamStates(sim, i);
        end
        
        function this = finalize(this, sim)
            %reconnect to server
            this.database.reopen();
            
            %meta data -- end time
            this.metadata.endTime = datestr(clock, 'yyyy-mm-dd HH:MM:SS');
            this.metadata.lengthSec = sim.states{this.stateIndex_time}.values;
            
            %for convenience
            md = this.metadata;
            
            %get contact, prompt for new contact if necessary
            contactWID = edu.stanford.covert.cell.kb.KnowledgeBaseUtil.getContact(...
                md, md.knowledgeBaseWID, this.database);
            
            %create simulation object in database, store meta data
            this.database.setNullValue(0);
            this.database.prepareStatement('CALL set_simulation("{S}","{S}","{Si}","{Si}","{S}","{S}","{S}","{S}","{S}","{Si}","{Si}","{Si}","{S}","{S}","{Si}")', ...
                md.shortDescription, md.longDescription, contactWID,...
                md.revision, md.differencesFromRevision, ...
                md.userName, md.hostName, md.ipAddress, md.outputDirectory, ...
                md.lengthSec, md.segmentSizeStep, md.downsampleStepSec,...
                md.startTime, md.endTime, ...
                md.knowledgeBaseWID);
            this.simulationWID = this.database.query().WID;
            
            %append last data bits and trim
            t = sim.states{this.stateIndex_time};
            time = t.values;
            i = ceil(time / this.downsampleStepSec);
            this.copyFromState(sim, i);
            
            states = sim.states;
            for j = 1:length(states)
                state = states{j};
                stateID = state.wholeCellModelID(7:end);
                names = state.stateNames;
                for k = 1:length(names)
                    name = names{k};
                    this.log.(stateID).(name) = this.log.(stateID).(name)(:, :, 1:i);
                end
            end
            
            %append last rand stream state and trim
            this.copyRandStreamStates(sim, i);
            this.randStreamStates.simulation = this.randStreamStates.simulation(:, 1:i);
            for k = 1:numel(sim.states)
                o = sim.states{k};
                this.randStreamStates.states.(o.wholeCellModelID(7:end)) = this.randStreamStates.states.(o.wholeCellModelID(7:end))(:, 1:i);
            end
            for k = 1:numel(sim.processes)
                o = sim.processes{k};
                this.randStreamStates.processes.(o.wholeCellModelID(9:end)) = this.randStreamStates.processes.(o.wholeCellModelID(9:end))(:, 1:i);
            end
            
            %store to options, parameters, fitting constants, state to database
            this.saveOptions(sim);
            this.saveParameters(sim);
            this.saveFittedConstants(sim);
            this.saveRandStreamStates(sim);
            this.saveTimecourses(sim);
            
            % set load error to false
            this.database.prepareStatement('CALL set_simulation_success("{Si}")', this.simulationWID);
            this.database.query();
            
            %clear log
            this.log = [];
        end
    end
    
    methods
        function validateMetadata(this)
            if ...
                    ~isfield(this.metadata, 'shortDescription') || ...
                    ~isfield(this.metadata, 'longDescription') || ...
                    ~isfield(this.metadata, 'email') || ...
                    ~isfield(this.metadata, 'firstName') || ...
                    ~isfield(this.metadata, 'lastName') || ...
                    ~isfield(this.metadata, 'affiliation') || ...
                    ~isfield(this.metadata, 'knowledgeBaseWID') || ...
                    ~isfield(this.metadata, 'revision') || ...
                    ~isfield(this.metadata, 'differencesFromRevision') || ...
                    ~isfield(this.metadata, 'userName') || ...
                    ~isfield(this.metadata, 'hostName') || ...
                    ~isfield(this.metadata, 'ipAddress')
                throw(MException('DatabaseLogger:invalidMetadata', 'Metadata missing'));
            end
        end
        
        %Copies the current state to a particular index in a segment.
        function copyFromState(this, sim, k)
            import edu.stanford.covert.cell.sim.util.DatabaseLogger;
            
            states = sim.states;
            for i = 1:length(states)
                state = states{i};
                stateID = state.wholeCellModelID(7:end);
                names = state.stateNames;
                for j = 1:length(names)
                    name = names{j};
                    data = state.(name);
                    if      size(data, 1) > size(this.log.(stateID).(name), 1) || ...
                            size(data, 2) > size(this.log.(stateID).(name), 2)
                        tmp = this.log.(stateID).(name);
                        this.log.(stateID).(name) = DatabaseLogger.allocateData(class(tmp), [size(data) size(tmp, 3)]);
                    end
                    
                    if isa(data, 'edu.stanford.covert.util.SparseMat')
                        if k == 1
                            this.log.(stateID).(name) = data;
                        else
                            this.log.(stateID).(name) = cat(3, this.log.(stateID).(name), data);
                        end
                    else
                        this.log.(stateID).(name)(1:size(data, 1), 1:size(data, 2), k) = data;
                    end
                end
            end
        end
        
        function copyRandStreamStates(this, sim, k)
            tmp = sim.getRandStreamStates();
            this.randStreamStates.simulation(:, k) = tmp.simulation;
            for i = 1:numel(sim.states)
                o = sim.states{i};
                this.randStreamStates.states.(o.wholeCellModelID(7:end))(:, k) = tmp.states.(o.wholeCellModelID(7:end));
            end
            for i = 1:numel(sim.processes)
                o = sim.processes{i};
                this.randStreamStates.processes.(o.wholeCellModelID(9:end))(:, k) = tmp.processes.(o.wholeCellModelID(9:end));
            end
        end
    end
    
    methods (Static)
        function metaStates = getMetaStates(sim, includeDependentStates)
            metaStates.names = cell(0, 1);
            metaStates.properties = cell(0, 4);
            for i = 1:length(sim.states)
                state = sim.states{i};
                stateID = strrep(state.wholeCellModelID, 'State_', '');
                metaStates.names = [
                    metaStates.names;
                    stateID];
                names = state.stateNames;
                if nargin >= 2 && includeDependentStates
                    names = [names; state.dependentStateNames]; %#ok<AGROW>
                end
                for j = 1:length(names)
                    metaStates.properties = [
                        metaStates.properties;
                        {stateID names{j} class(state.(names{j})) size(state.(names{j}))}];
                end
            end
        end
        
        %Allocates space for one segment's worth of time course data.
        function log = allocateMemory(metaStates, segmentSizeStep)
            import edu.stanford.covert.cell.sim.util.DatabaseLogger;
            
            log = struct;
            
            for i = 1:size(metaStates.names, 1)
                log.(metaStates.names{i, 1}) = struct;
            end
            
            for i = 1:size(metaStates.properties, 1)
                log.(metaStates.properties{i, 1}).(metaStates.properties{i, 2}) = DatabaseLogger.allocateData(...
                    metaStates.properties{i, 3}, [metaStates.properties{i, 4}, segmentSizeStep]);
            end
        end
        
        function data = allocateData(dataType, sz)
            import edu.stanford.covert.cell.sim.util.DatabaseLogger;
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.util.SparseMat;
            
            switch dataType
                case 'edu.stanford.covert.util.SparseMat', data = SparseMat([], [], sz);
                case 'edu.stanford.covert.util.CircularSparseMat',  data = CircularSparseMat([], [], sz, 1);
                case 'char',     data = char(zeros(sz));
                case 'logical',  data = false(sz);
                case 'cell',     data = cell(sz);
                otherwise,
                    try
                        data = zeros(sz, dataType);
                    catch exception
                        throw(MException('DatabaseLogger:invalidDataType', sprintf('%s is not a numeric data type', dataType)));
                    end
            end
        end
    end
    
    %save to database
    methods (Access = protected)
        %Saves options to database
        function saveOptions(this, sim)
            %create option entries
            this.saveOptionsHelper(rmfield(sim.getOptions(), {'states','processes'}), []);
            
            %states
            for i = 1:length(sim.states)
                this.saveOptionsHelper(sim.states{i}.getOptions(), sim.states{i});
            end
            
            %processes
            for i = 1:length(sim.processes)
                this.saveOptionsHelper(sim.processes{i}.getOptions(), sim.processes{i});
            end
        end
        
        function saveOptionsHelper(this, options, module)
            optionNames = fieldnames(options);
            for i = 1:length(optionNames)
                value = edu.stanford.covert.io.jsonFormat(options.(optionNames{i}));

                if ~isempty(module)
                    this.database.prepareStatement('CALL set_simulation_option("{S}","{S}",NULL,"{S}","{Si}","{Si}")',...
                        module.wholeCellModelID, optionNames{i}, value, this.simulationWID, this.metadata.knowledgeBaseWID);
                else
                    this.database.prepareStatement('CALL set_simulation_option(NULL,"{S}",NULL,"{S}","{Si}","{Si}")',...
                        optionNames{i}, value, this.simulationWID, this.metadata.knowledgeBaseWID);
                end
                this.database.query();
            end
        end
        
        % Saves parameters to database
        function saveParameters(this, sim)
            %states
            for i = 1:length(sim.states)
                this.saveParametersHelper(sim.states{i});
            end
            
            %processes
            for i = 1:length(sim.processes)
                this.saveParametersHelper(sim.processes{i});
            end
        end
        
        function saveParametersHelper(this, module)
            knowledgeBaseWID = this.metadata.knowledgeBaseWID;

            objectFieldNames = fieldnames(module);
            for i = 1:length(module.parameterNames)
                if ~any(strcmp(module.parameterNames{i}, objectFieldNames))
                    continue;
                end

                if isempty(module.parameterIndexs{i})
                    value = module.(module.parameterNames{i});
                else
                    value = module.(module.parameterNames{i})(module.(module.parameterIndexs{i}));
                end

                if iscell(value)
                    value = strjoin(';',value{:});
                elseif ~isempty(value)
                    strValue = [];
                    for j = 1:length(value)
                        strValue = [strValue num2str(value(j)) ';']; %#ok<AGROW>
                    end
                    value = strValue(1:end-1);
                end

                this.database.prepareStatement('CALL set_simulation_parameter("{S}","{S}","{S}","{S}","{Si}","{Si}")', ...
                    module.wholeCellModelID, module.parameterNames{i}, module.parameterIndexs{i}, value, this.simulationWID, knowledgeBaseWID);
                this.database.query();
            end
        end
        
        %Saves fit contants to database
        function saveFittedConstants(this, sim)
            %create fit constants entries -- states
            for i = 1:length(sim.states)
                this.saveFittedConstantsHelper(sim.states{i});
            end
            
            %create fit constants entries -- processes
            for i = 1:length(sim.processes)
                this.saveFittedConstantsHelper(sim.processes{i});
            end
        end
        
        function saveFittedConstantsHelper(this, module)
            import edu.stanford.covert.io.jsonFormat;
            
            knowledgeBaseWID = this.metadata.knowledgeBaseWID;
            
            fittedConstantNames = module.fittedConstantNames;
            for i = 1:length(fittedConstantNames)
                fittedConstantName = fittedConstantNames{i};
                this.database.prepareStatement(...
                    'CALL set_simulation_fittedconstant("{S}","{S}",NULL,"{S}",NULL,"{Si}","{Si}")',...
                    module.wholeCellModelID, fittedConstantName, jsonFormat(module.(fittedConstantName)),...
                    this.simulationWID, knowledgeBaseWID);
                this.database.query();
            end
        end
        
        %Saves rand stream states to database
        function saveRandStreamStates(this, sim)
            import edu.stanford.covert.io.jsonFormat;
            
            knowledgeBaseWID = this.metadata.knowledgeBaseWID;
            
            %create rand stream states entries -- simulation
            this.database.prepareStatement(...
                'CALL set_simulation_randstreamstate(NULL,NULL,NULL,"{S}",NULL,"{Si}","{Si}")',...
                jsonFormat(this.randStreamStates.simulation),...
                this.simulationWID, knowledgeBaseWID);
            this.database.query();
            
            %create rand stream states entries -- states
            for i = 1:length(sim.states)
                o = sim.states{i};
                this.database.prepareStatement(...
                    'CALL set_simulation_randstreamstate("{S}",NULL,NULL,"{S}",NULL,"{Si}","{Si}")',...
                    o.wholeCellModelID, jsonFormat(this.randStreamStates.states.(o.wholeCellModelID(7:end))),...
                    this.simulationWID, knowledgeBaseWID);
                this.database.query();
            end
            
            %create rand stream states entries -- processes
            for i = 1:length(sim.processes)
                o = sim.processes{i};
                this.database.prepareStatement(...
                    'CALL set_simulation_randstreamstate("{S}",NULL,NULL,"{S}",NULL,"{Si}","{Si}")',...
                    o.wholeCellModelID, jsonFormat(this.randStreamStates.processes.(o.wholeCellModelID(9:end))),...
                    this.simulationWID, knowledgeBaseWID);
                this.database.query();
            end
        end
        
        %Saves time courses to database
        function saveTimecourses(this, sim)
            import edu.stanford.covert.io.jsonFormat;
            
            knowledgeBaseWID = this.metadata.knowledgeBaseWID;
            
            for i = 1:length(sim.states)
                state = sim.states{i};
                stateID = state.wholeCellModelID(7:end);
                names = state.stateNames;
                for j = 1:length(names)
                    this.database.prepareStatement(...
                        'CALL set_simulation_timecourse("{S}","{S}",NULL,"{S}",NULL,"{Si}","{Si}")',...
                        state.wholeCellModelID, names{j}, jsonFormat(this.log.(stateID).(names{j})), ...
                        this.simulationWID, knowledgeBaseWID);
                    this.database.query();
                end
            end
        end
    end
    
    %load from database
    methods (Static)
        %load simulation from database (already downsampled)
        function [states, metadata, options, parameters, fittedConstants, randStreamStates] = load(sim, database, simulationWID)
            import edu.stanford.covert.cell.sim.util.DatabaseLogger;
            
            metadata = DatabaseLogger.loadMetadata(database, simulationWID);
            options = DatabaseLogger.loadOptions(sim, database, simulationWID);
            parameters = DatabaseLogger.loadParameters(sim, database, simulationWID);
            fittedConstants = DatabaseLogger.loadFittedConstants(sim, database, simulationWID);
            randStreamStates = DatabaseLogger.loadRandStreamStates(sim, database, simulationWID);
            states = DatabaseLogger.loadTimecourses(sim, database, simulationWID);
        end
        
        %Loads simulation meta data from database
        function value = loadMetadata(database, simulationWID)
            database.setNullValue(0);
            database.prepareStatement('CALL get_simulation_summary("{Si}")', simulationWID);
            result = database.query();
            
            value = struct;
            value.knowledgeBaseWID        = uint32(result.KnowledgeBaseWID(1));
            value.shortDescription        = result.ShortDescription{1};
            value.longDescription         = result.LongDescription{1};
            value.revision                = result.Revision(1);
            value.differencesFromRevision = result.DifferencesFromRevision{1};
            value.lengthSec               = result.Length(1);
            value.downsampleStepSec       = result.SampleStep(1);
            value.segmentSizeStep         = result.SegmentStep(1);
            value.firstName               = result.FirstName{1};
            value.lastName                = result.LastName{1};
            value.affiliation             = result.Affiliation{1};
            value.email                   = result.Email{1};
            value.userName                = result.UserName{1};
            value.hostName                = result.HostName{1};
            value.ipAddress               = result.IPAddress{1};
            value.outputDirectory         = result.OutputDirectory{1};
            value.startTime               = result.StartDate{1}(1:19);
            value.endTime                 = result.EndDate{1}(1:19);
        end
        
        %Loads simulation options from database
        function value = loadOptions(sim, database, simulationWID)
            database.setNullValue(0);
            database.prepareStatement('CALL get_simulation_options("{Si}")', simulationWID);
            result = database.query();
            
            value = struct('states', struct, 'processes', struct);
            for i = 1:numel(sim.states)
                value.states.(sim.states{i}.wholeCellModelID(7:end)) = struct;
            end
            for i = 1:numel(sim.processes)
                value.processes.(sim.processes{i}.wholeCellModelID(9:end)) = struct;
            end
            
            for i = 1:length(result.WID)
                if isempty(result.Module{i})
                    value.(result.Name{i}) = edu.stanford.covert.io.jsonParse(result.Value{i});
                elseif result.Module{i}(1) == 'S'
                    value.states.(result.Module{i}(7:end)).(result.Name{i}) = edu.stanford.covert.io.jsonParse(result.Value{i});
                else
                    value.processes.(result.Module{i}(9:end)).(result.Name{i}) = edu.stanford.covert.io.jsonParse(result.Value{i});
                end
            end
        end
        
        %Loads simulation parameters from database
        function value = loadParameters(sim, database, simulationWID)
            database.setNullValue(0);
            database.prepareStatement('CALL get_simulation_parameters("{Si}")', simulationWID);
            result = database.query();
            
            value = struct('states', struct, 'processes', struct);
            for i = 1:numel(sim.states)
                value.states.(sim.states{i}.wholeCellModelID(7:end)) = struct;
            end
            for i = 1:numel(sim.processes)
                value.processes.(sim.processes{i}.wholeCellModelID(9:end)) = struct;
            end
            
            for i = 1:length(result.WID)
                if isempty(result.Value{i})
                    tmp = {};
                else
                    tmp = strsplit(';', result.Value{i})';
                    tmpNum = zeros(size(tmp));
                    for j = 1:length(tmp)
                        tmpNum(j) = str2double(tmp{j});
                    end
                    if ~any(isnan(tmpNum))
                        tmp = tmpNum;
                    end
                end
                
                if result.Module{i}(1) == 'S'
                    if isempty(result.Index{i})
                        value.states.(result.Module{i}(7:end)).(result.Name{i}) = tmp;
                    else
                        o = sim.state(result.Module{i}(7:end));
                        value.states.(result.Module{i}(7:end)).(result.Name{i})(o.(result.Index{i}), 1) = tmp;
                    end
                else
                    if isempty(result.Index{i})
                        value.processes.(result.Module{i}(9:end)).(result.Name{i}) = tmp;
                    else
                        o = sim.process(result.Module{i}(9:end));
                        value.processes.(result.Module{i}(9:end)).(result.Name{i})(o.(result.Index{i}), 1) = tmp;
                    end
                end
            end
        end
        
        %Loads fitted constants from database
        function value = loadFittedConstants(sim, database, simulationWID)
            database.setNullValue(0);
            database.prepareStatement('CALL get_simulation_fittedconstants("{Si}")', simulationWID);
            result = database.query();
            
            value = struct('states', struct, 'processes', struct);
            for i = 1:numel(sim.states)
                value.states.(sim.states{i}.wholeCellModelID(7:end)) = struct;
            end
            for i = 1:numel(sim.processes)
                value.processes.(sim.processes{i}.wholeCellModelID(9:end)) = struct;
            end
            
            for i = 1:length(result.WID)
                if result.Module{i}(1) == 'S'
                    value.states.(result.Module{i}(7:end)).(result.Name{i}) = edu.stanford.covert.io.jsonParse(result.Value{i});
                else
                    value.processes.(result.Module{i}(9:end)).(result.Name{i}) = edu.stanford.covert.io.jsonParse(result.Value{i});
                end
            end
        end
        
        %Loads rand stream states from database
        function value = loadRandStreamStates(~, database, simulationWID)
            database.setNullValue(0);
            database.prepareStatement('CALL get_simulation_randstreamstates("{Si}")', simulationWID);
            result = database.query();
            
            value = struct('simulation', [], 'states', struct, 'processes', struct);
            for i = 1:length(result.WID)
                if isempty(result.Module{i})
                    value.simulation = edu.stanford.covert.io.jsonParse(result.Value{i});
                elseif result.Module{i}(1) == 'S'
                    value.states.(result.Module{i}(7:end)) = edu.stanford.covert.io.jsonParse(result.Value{i});
                else
                    value.processes.(result.Module{i}(9:end)) = edu.stanford.covert.io.jsonParse(result.Value{i});
                end
            end
        end
        
        %Loads time courses from database
        function states = loadTimecourses(sim, database, simulationWID)
            import edu.stanford.covert.cell.sim.util.DatabaseLogger;
            
            database.setNullValue(0);
            database.prepareStatement('CALL get_simulation_timecourses("{Si}")', simulationWID);
            result = database.query();
            
            metaStates = DatabaseLogger.getMetaStates(sim);
            states = DatabaseLogger.allocateMemory(metaStates, 0);
            
            for i = 1:length(result.WID)
                states.(result.Module{i}(7:end)).(result.Name{i}) = edu.stanford.covert.io.jsonParse(result.Value{i});
            end
        end
    end
end