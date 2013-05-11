% Simulation disk logger.
%   Stores and retrieves simulated dynamics and simulation meta data
%   to/from disk. Generally, each batch of simulations is stored as folder
%   containing 1 subfolder for each simulation. Each simulation subfolder
%   contains the simulated dynamics of a single cell. The simulated
%   dynamics of each single cell are stored in two ways:
%   1) Simulated dynamics are stored indexed primarily by time with the
%      file name pattern state-(\d+).mat
%   2) Simulated dynamics are stored indexed primarily by state with the
%      file name pattern state-(state name)-(property name).mat
%   During the execution of each simulation results are stored in the first
%   form. However this form is inefficient for most analyses. Consequently,
%   after the completion of each simulation, simulated dynamics are
%   reindexed by state.
%
%   The SimuationEnsemble class provides additional methods for retrieving
%   the simulated dynamics of multiple cells and entire cellular
%   populations. The SimulationDiskUtil class provides additional methods
%   for retrieving simulation meta data.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/10/2011
classdef DiskLogger < edu.stanford.covert.cell.sim.util.Logger
    %options
    properties (SetAccess = protected)
        metadata          %metadata
        outputDirectory   %output directory
        segmentSizeStep   %max state history to cache in memory during a run (s)
        verbosity         %verbosity
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
        function this = DiskLogger(outputDirectory, segmentSizeStep, metadata, verbosity)
            if ~exist(outputDirectory, 'dir')
                mkdir(outputDirectory);
            end
            if ~exist('metadata', 'var')
                metadata = struct;
            end
            if ~exist('verbosity', 'var')
                verbosity = 0;
            end
            
            this.outputDirectory = outputDirectory;
            this.segmentSizeStep = segmentSizeStep;
            this.metadata = metadata;
            this.verbosity = verbosity;
        end
    end
    
    methods
        function setOptions(this, varargin)
            if isstruct(varargin{1})
                options = varargin{1}; %#ok<*PROP>
            else
                options = struct(varargin{:});
            end
            
            metaClass = metaclass(this);
            fields = intersect(fieldnames(options), cellfun(@(x) x.Name, metaClass.Properties, 'UniformOutput',false));
            for i = 1:numel(fields)
                this.(fields{i}) = options.(fields{i});
            end
        end
        
        function this = addMetadata(this, varargin)
            if isstruct(varargin{1})
                metadata = varargin{1}; %#ok<*PROP>
            else
                metadata = struct(varargin{:});
            end
            
            names = fieldnames(metadata);
            for i = 1:length(names)
                this.metadata.(names{i}) = metadata.(names{i});
            end
        end
        
        function this = initialize(this, sim)
            %% validate options
            %segment length
            validateattributes(sim.lengthSec / sim.stepSizeSec / this.segmentSizeStep, {'numeric'}, {'real', 'nonnegative', 'integer'});
            
            %downsample step
            validateattributes(sim.lengthSec / sim.stepSizeSec / this.segmentSizeStep, {'numeric'}, {'real', 'nonnegative', 'integer'});
            
            %metadata
            this.validateMetadata();
            
            %% indices
            this.stateIndex_time = sim.stateIndex('Time');
            
            %% metdata
            %states
            this.metadata.stateNames = cell(0, 2);
            this.metadata.dependentStateNames = cell(0, 2);
            for j = 1:length(sim.states)
                state = sim.states{j};
                stateID = {state.wholeCellModelID(7:end)};
                this.metadata.stateNames = [
                    this.metadata.stateNames
                    repmat(stateID, numel(state.stateNames), 1) state.stateNames
                    ];
                this.metadata.dependentStateNames = [
                    this.metadata.dependentStateNames
                    repmat(stateID, numel(state.dependentStateNames), 1) state.dependentStateNames
                    ];
            end
            
            %start time
            this.metadata.startTime = datestr(clock, 'yyyy-mm-dd HH:MM:SS');
            this.metadata.endTime = [];
            this.metadata.lengthSec = [];
            
            %output directory
            this.metadata.outputDirectory = strrep(this.outputDirectory, '\', '/');
            if ~isempty(this.outputDirectory) && ~exist(this.outputDirectory, 'dir')
                mkdir(this.outputDirectory);
            end
            
            %segment step
            this.metadata.segmentSizeStep = this.segmentSizeStep;
            
            %down sample step
            this.metadata.downsampleStepSec = sim.stepSizeSec;
            
            %% store initial state as "segment zero"
            metaStates = this.getMetaStates(sim, true);
            this.log = this.allocateMemory(metaStates, 1);
            this.copyFromState(sim, 1);
            this.saveSegmentToDisk(0);
            
            %% allocate memory
            %cell state
            this.log = this.allocateMemory(metaStates, this.segmentSizeStep);
            
            %rand stream states
            this.randStreamStates = struct('simulation', [], 'states', struct(), 'processes', struct());
            tmp = sim.getRandStreamStates();
            this.randStreamStates.simulation = zeros(numel(tmp.simulation), sim.lengthSec + 1);
            for i = 1:numel(sim.states)
                o = sim.states{i};
                this.randStreamStates.states.(o.wholeCellModelID(7:end)) = zeros(numel(tmp.states.(o.wholeCellModelID(7:end))), sim.lengthSec + 1);
            end
            for i = 1:numel(sim.processes)
                o = sim.processes{i};
                this.randStreamStates.processes.(o.wholeCellModelID(9:end)) = zeros(numel(tmp.processes.(o.wholeCellModelID(9:end))), sim.lengthSec + 1);
            end
            this.copyRandStreamStates(sim);
        end
        
        function this = append(this, sim)
            %cell state
            nSteps = sim.states{this.stateIndex_time}.values / sim.stepSizeSec;
            i = mod(nSteps - 1, this.segmentSizeStep) + 1;
            this.copyFromState(sim, i);
            
            if i == this.segmentSizeStep
                this.saveSegmentToDisk(nSteps / this.segmentSizeStep);
            end
            
            %rand stream states
            this.copyRandStreamStates(sim);
        end
        
        function this = finalize(this, sim)
            %references
            states = sim.states;
            
            %meta data -- end time
            this.metadata.endTime = datestr(clock, 'yyyy-mm-dd HH:MM:SS');
            this.metadata.lengthSec = sim.states{this.stateIndex_time}.values;
            
            %append last data bits, trim last segment, save, and clear
            nSteps = sim.states{this.stateIndex_time}.values / sim.stepSizeSec;
            i = mod(nSteps - 1, this.segmentSizeStep) + 1;
            this.copyFromState(sim, i);
            
            for j = 1:length(states)
                state = states{j};
                stateID = state.wholeCellModelID(7:end);
                names = [state.stateNames; state.dependentStateNames];
                for k = 1:length(names)
                    name = names{k};
                    this.log.(stateID).(name) = this.log.(stateID).(name)(:, :, 1:i);
                end
            end
            
            this.saveSegmentToDisk(ceil(nSteps / this.segmentSizeStep));
            
            %append rand stream states, trim, and save
            this.copyRandStreamStates(sim);
            this.randStreamStates.simulation = this.randStreamStates.simulation(:, 1:nSteps+1);
            for i = 1:numel(sim.states)
                o = sim.states{i};
                this.randStreamStates.states.(o.wholeCellModelID(7:end)) = this.randStreamStates.states.(o.wholeCellModelID(7:end))(:, 1:nSteps+1);
            end
            for i = 1:numel(sim.processes)
                o = sim.processes{i};
                this.randStreamStates.processes.(o.wholeCellModelID(9:end)) = this.randStreamStates.processes.(o.wholeCellModelID(9:end))(:, 1:nSteps+1);
            end
            
            %store to disk
            this.saveMetadata(sim);
            this.saveOptions(sim);
            this.saveParameters(sim);
            this.saveFittedConstants(sim);
            this.saveRandStreamStates(sim);
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
                throw(MException('DiskLogger:invalidMetadata', 'Metadata missing'));
            end
        end
        
        %Copies the current state to a particular index in a segment.
        function copyFromState(this, sim, k)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            states = sim.states;
            for i = 1:length(states)
                state = states{i};
                stateID = state.wholeCellModelID(7:end);
                names = [state.stateNames; state.dependentStateNames];
                for j = 1:length(names)
                    name = names{j};
                    data = state.(name);
                    if      size(data, 1) > size(this.log.(stateID).(name), 1) || ...
                            size(data, 2) > size(this.log.(stateID).(name), 2)
                        tmp = this.log.(stateID).(name);
                        this.log.(stateID).(name) = DiskLogger.allocateData(class(tmp), [size(data) size(tmp, 3)]);
                        this.log.(stateID).(name)(1:size(tmp,1), 1:size(tmp, 2), :) = tmp;
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
        
        function copyRandStreamStates(this, sim)
            nSteps = sim.states{this.stateIndex_time}.values / sim.stepSizeSec;
            tmp = sim.getRandStreamStates();
            
            this.randStreamStates.simulation(:, nSteps + 1) = tmp.simulation;
            for i = 1:numel(sim.states)
                o = sim.states{i};
                this.randStreamStates.states.(o.wholeCellModelID(7:end))(:, nSteps + 1) = tmp.states.(o.wholeCellModelID(7:end));
            end
            for i = 1:numel(sim.processes)
                o = sim.processes{i};
                this.randStreamStates.processes.(o.wholeCellModelID(9:end))(:, nSteps + 1) = tmp.processes.(o.wholeCellModelID(9:end));
            end
        end
        
        function saveMetadata(this, ~)
            this.saveStructToDisk(this.outputDirectory, 'metadata.mat', this.metadata);
        end
        
        function saveOptions(this, sim)
            this.saveStructToDisk(this.outputDirectory, 'options.mat', sim.getOptions);
        end
        
        function saveParameters(this, sim)
            this.saveStructToDisk(this.outputDirectory, 'parameters.mat', sim.getParameters);
        end
        
        function saveFittedConstants(this, sim)
            this.saveStructToDisk(this.outputDirectory, 'fittedConstants.mat', sim.getFittedConstants);
        end
        
        function saveRandStreamStates(this, ~)
            this.saveStructToDisk(this.outputDirectory, 'randStreamStates.mat', this.randStreamStates);
        end
        
        function saveSegmentToDisk(this, segmentIdx)
            this.saveStructToDisk(this.outputDirectory, ['state-' num2str(segmentIdx) '.mat'], this.log);
        end
        
        %number of time steps per segment
        function value = getSegmentLength(this, sim)
            value = this.segmentSizeSec / sim.stepSizeSec;
        end
        
        %number of time segments
        function value = getNumSegments(this, sim)
            value = sim.lengthSec / this.getSegmentLength;
        end
        
        function clearLog(this)
            this.log = [];
        end
    end
    
    %load from disk
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
                if nargin >= 2 && exist('includeDependentStates', 'var') && includeDependentStates
                    names = [names; state.dependentStateNames]; %#ok<AGROW>
                end
                for j = 1:length(names)
                    metaStates.properties = [
                        metaStates.properties;
                        {stateID names{j} class(state.(names{j})) size(state.(names{j}))}];
                end
            end
         end
        
        function [stateNames, stateDataTypes, isStateDataBuiltinDenseMatrix] = getAvailableStates(outputDirectory)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            log0 = DiskLogger.loadStructFromDisk(sprintf('%s/state-%d.mat', outputDirectory, 0));
            stateNames = cell(0, 2);
            stateDataTypes = cell(0, 1);
            isStateDataBuiltinDenseMatrix = false(0, 1);
            fields = fieldnames(log0);
            for i = 1:numel(fields)
                fieldsi = fieldnames(log0.(fields{i}));
                stateNames = [
                    stateNames;
                    repmat(fields(i), size(fieldsi, 1), 1) fieldsi(:)]; %#ok<AGROW>
                for j = 1:numel(fieldsi)
                    stateDataTypes = [
                        stateDataTypes;
                        class(log0.(fields{i}).(fieldsi{j}))
                        ]; %#ok<AGROW>
                    isStateDataBuiltinDenseMatrix = [
                        isStateDataBuiltinDenseMatrix;
                        isnumeric(log0.(fields{i}).(fieldsi{j})) && ~issparse(log0.(fields{i}).(fieldsi{j}))
                        ]; %#ok<AGROW>
                end
            end
        end
        
        %load simulation from disk
        %- stateNames is a n x 2 cell array; 1st col->state ID, 2nd col->property name
        %- downsampleType can be 'extract' or 'mean'
        function [states, metadata, options, parameters, fittedConstants, randStreamStates] = load(...
                outputDirectory, stateNames, ...
                initTime, finTime, downsampleStepSec, ...
                downsampleType)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
                        
            metadata = DiskLogger.loadMetadata(outputDirectory); metadata.downsampleStepSec = downsampleStepSec;
            options = DiskLogger.loadOptions(outputDirectory);
            if nargout >= 4
                parameters = DiskLogger.loadParameters(outputDirectory);
            end
            if nargout >= 5
                fittedConstants = DiskLogger.loadFittedConstants(outputDirectory);
            end
            if nargout >= 6
                randStreamStates = DiskLogger.loadRandStreamStates(outputDirectory);
            end
            states = DiskLogger.loadTimecourses(outputDirectory, stateNames, initTime, finTime, downsampleStepSec, downsampleType, options, metadata);
        end
        
        function value = loadMetadata(outputDirectory, varargin)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            value = DiskLogger.loadStructFromDisk([outputDirectory '/metadata.mat'], varargin{:});
        end
        
        function value = loadOptions(outputDirectory, varargin)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            value = DiskLogger.loadStructFromDisk([outputDirectory '/options.mat'], varargin{:});
        end
        
        function value = loadParameters(outputDirectory, varargin)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            value = DiskLogger.loadStructFromDisk([outputDirectory '/parameters.mat'], varargin{:});
        end
        
        function value = loadFittedConstants(outputDirectory, varargin)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            value = DiskLogger.loadStructFromDisk([outputDirectory '/fittedConstants.mat'], varargin{:});
        end
        
        function value = loadRandStreamStates(outputDirectory, varargin)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            value = DiskLogger.loadStructFromDisk([outputDirectory '/randStreamStates.mat'], varargin{:});
        end
        
        function states = loadTimecourses(outputDirectory, stateNames, initTime, finTime, downsampleStepSec, ...
                downsampleType, options, metadata, method)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            %process options
            if ~exist('downsampleStepSec', 'var') || isempty(downsampleStepSec)
                downsampleStepSec = options.stepSizeSec;
            elseif mod(downsampleStepSec, options.stepSizeSec) ~= 0
                throw(MException('DiskLogger:invalid', 'downsampleStepSec must be a multiple of stepSizeSec'));
            end
            if ~exist('downsampleType', 'var') || isempty(downsampleType)
                downsampleType = 'extract';
            elseif ~ismember(downsampleType, {'extract', 'mean'})
                throw(MException('DiskLogger:invalid', 'downsampleType must be one of ''extract'' or ''mean'''));
            end
            if strcmp(downsampleType, 'mean') && downsampleStepSec == options.stepSizeSec
                downsampleType = 'extract';
            end
            if isempty(initTime)
                initTime = 1;
            else
                validateattributes(initTime, {'numeric'}, {'>=', 1, '<=', metadata.lengthSec});
            end
            if isempty(finTime)
                finTime = metadata.lengthSec;
            else
                validateattributes(finTime, {'numeric'}, {'>=', initTime, '<=', metadata.lengthSec});
            end
            if ~exist('method', 'var')
                method = 'byState';
            end
            
            %allocate memory
            numSegments = ceil(metadata.lengthSec / options.stepSizeSec / metadata.segmentSizeStep);
            if ischar(stateNames) && stateNames(1) == '-'
                logN = DiskLogger.loadStructFromDisk(sprintf('%s/state-%d.mat', outputDirectory, numSegments));
                metaStates.names = fieldnames(logN);
                switch stateNames
                    case '-all'
                        if isfield(metadata, 'stateNames')
                            stateNames = [metadata.stateNames; metadata.dependentStateNames];
                        else
                            %for backwards compatability with data from before dependent
                            %states were stored (Revision < 1733)
                            stateNames = DiskLogger.getAvailableStates(outputDirectory);
                        end
                    case '-independent'
                        stateNames = metadata.stateNames;
                    case '-dependent'
                        stateNames = metadata.dependentStateNames;
                    otherwise
                        throw(MException('DiskLogger:error', 'Invalid stateNames option ''%s'', stateNames'));
                end
            else
                logN = DiskLogger.loadStructFromDisk(sprintf('%s/state-%d.mat', outputDirectory, numSegments), stateNames(:, 1));                
                metaStates.names = intersect(fieldnames(logN), unique(stateNames(:, 1)));
            end
            metaStates.properties = [stateNames cell(size(stateNames, 1), 2)];
            for i = 1:size(stateNames, 1)
                tmp = logN.(stateNames{i, 1}).(stateNames{i, 2});
                metaStates.properties{i, 3} = class(tmp);
                metaStates.properties{i, 4} = [size(tmp, 1) size(tmp, 2)];
                if size(stateNames, 2) >= 3
                    if ischar(stateNames{i, 3}) && ismember(stateNames{i, 3}, {'-sum'; '-nnz'})
                        metaStates.properties{i, 3} = 'double';
                        metaStates.properties{i, 4}(1) = 1;
                    elseif isnumeric(stateNames{i, 3})
                        metaStates.properties{i, 4}(1) = size(stateNames{i, 3}, 1);
                    end
                end
                if size(stateNames, 2) >= 4
                    if ischar(stateNames{i, 4}) && ismember(stateNames{i, 4}, {'-sum'; '-nnz'})
                        metaStates.properties{i, 3} = 'double';
                        metaStates.properties{i, 4}(2) = 1;
                    elseif isnumeric(stateNames{i, 4})
                        metaStates.properties{i, 4}(2) = size(stateNames{i, 4}, 1);
                    end
                end
            end
            states = DiskLogger.allocateMemory(metaStates, ceil((finTime-initTime) / downsampleStepSec));
            
            %downsample timecourse
            if strcmp(method, 'byState') && exist([outputDirectory filesep 'state-Time-values.mat'], 'file') && finTime - initTime > metadata.segmentSizeStep
                switch downsampleType
                    case 'extract', states = DiskLogger.downSample_extract_byState(states, options, metadata, outputDirectory, stateNames, initTime, finTime, downsampleStepSec);
                    case 'mean',    states = DiskLogger.downSample_mean_byState(   states, options, metadata, outputDirectory, stateNames, initTime, finTime, downsampleStepSec);
                end
            else
                switch downsampleType
                    case 'extract', states = DiskLogger.downSample_extract_byTime(states, options, metadata, outputDirectory, stateNames, initTime, finTime, downsampleStepSec);
                    case 'mean',    states = DiskLogger.downSample_mean_byTime(   states, options, metadata, outputDirectory, stateNames, initTime, finTime, downsampleStepSec);
                end
            end
        end
        
        function states = downSample_extract_byState(states, options, metadata, outputDirectory, stateNames, initTime, finTime, downsampleStepSec)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            times = ceil(initTime/options.stepSizeSec : downsampleStepSec/options.stepSizeSec : ceil(finTime/downsampleStepSec)*downsampleStepSec/options.stepSizeSec) * options.stepSizeSec;
            times(end) = finTime;
            
            data = DiskLogger.loadStructFromDisk(sprintf('%s/state-%s-%s.mat', outputDirectory, 'Time', 'values'));
            data = data.data;
            allTimes = isequal(times, permute(data, [2 3 1]));
            if allTimes
                idxs = ':';
            else
                idxs = ismembc2(times, permute(data, [2 3 1]));
            end
            
            sNames = unique(stateNames(:, 1));
            for i = 1:numel(sNames);
                stateName = sNames{i};
                propNames = stateNames(strcmp(stateNames(:, 1), stateName), 2:end);
                for j = 1:size(propNames, 1)
                    fileName = sprintf('%s/state-%s-%s.mat', outputDirectory, stateName, propNames{j, 1});
                    if ~exist(fileName, 'file')
                        states = DiskLogger.downSample_extract_byTime(...
                            states, options, metadata, outputDirectory, [{stateName} propNames(j, :)], initTime, finTime, downsampleStepSec);
                        continue;
                    end
                    
                    tmp = DiskLogger.loadStructFromDisk(fileName);
                    tmp = tmp.data;
                    
                    calcMargins = false(1, 2);
                    calcSums = false(1, 2);
                    calcNnzs = false(1, 2);
                    s.type = '()';
                    s.subs = {':' ':' ':'};
                    if size(stateNames, 2) > 2
                        s.subs(1:2) = propNames(j, 2:3);
                        calcMargins = cellfun(@(x) isnumeric(x) && size(x, 2) > 1, propNames(j, 2:3));
                        calcSums = strcmp(propNames(j, 2:3), '-sum');
                        calcNnzs = strcmp(propNames(j, 2:3), '-nnz');
                        s.subs(calcSums | calcNnzs) = {':'};
                    end
                    s.subs{3} = idxs;
                    
                    tmp = subsref(tmp, s);
                    
                    if calcMargins(1)
                        siz = size(tmp);
                        siz = [size(s.subs{1}) siz(2:end)];
                        tmp = permute(sum(reshape(tmp, siz), 2), [1 3 4 2]);
                    elseif calcSums(1)
                        tmp = sum(tmp, 1);
                    elseif calcNnzs(1)
                        if iscell(tmp)
                            tmp = ~cellfun(@isempty, tmp);
                        end
                        tmp = sum(tmp ~= 0, 1);
                    end
                    if calcMargins(2)
                        siz = size(tmp);
                        siz = [siz(1) size(s.subs{2}) siz(3:end)];
                        tmp = permute(sum(reshape(tmp, siz), 3), [1 2 4 3]);
                    elseif calcSums(2)
                        tmp = sum(tmp, 2);
                    elseif calcNnzs(2)
                        if iscell(tmp)
                            tmp = ~cellfun(@isempty, tmp);
                        end
                        tmp = sum(tmp ~= 0, 2);
                    end
                    
                    if  isa(tmp, 'edu.stanford.covert.util.SparseMat') && ...
                            ~isa(states.(stateName).(propNames{j, 1}), 'edu.stanford.covert.util.SparseMat')
                        tmp = full(tmp);
                    end
                    
                    states.(stateName).(propNames{j, 1}) = tmp;
                    
                    clear tmp;
                end
            end
        end
        
        function states = downSample_mean_byState(states, options, metadata, outputDirectory, stateNames, initTime, finTime, downsampleStepSec)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            times = ceil(initTime/options.stepSizeSec : downsampleStepSec/options.stepSizeSec : ceil(finTime/downsampleStepSec)*downsampleStepSec/options.stepSizeSec) * options.stepSizeSec;
            
            data = DiskLogger.loadStructFromDisk(sprintf('%s/state-%s-%s.mat', outputDirectory, 'Time', 'values'));
            data = data.data;
            idxs = ismembc2(times, permute(data, [2 3 1]));
            
            sNames = unique(stateNames(:, 1));
            for i = 1:numel(sNames);
                stateName = sNames{i};
                propNames = stateNames(strcmp(stateNames(:, 1), stateName), 2:end);
                
                for j = 1:size(propNames, 1)
                    fileName = sprintf('%s/state-%s-%s.mat', outputDirectory, stateName, propNames{j, 1});
                    if ~exist(fileName, 'file')
                        states = DiskLogger.downSample_mean_byTime(...
                            states, options, metadata, outputDirectory, [{stateName} propNames(j, :)], initTime, finTime, downsampleStepSec);
                        continue;
                    end
                    data = DiskLogger.loadStructFromDisk(fileName);
                    data = data.data;
                    
                    calcMargins = false(1, 2);
                    calcSums = false(1, 2);
                    calcNnzs = false(1, 2);
                    s.type = '()';
                    s.subs = {':' ':' ':'};
                    if size(stateNames, 2) > 2
                        s.subs(1:2) = propNames(j, 2:3);
                        calcMargins = cellfun(@(x) isnumeric(x) && size(x, 2) > 1, propNames(j, 2:3));
                        calcSums = strcmp(propNames(j, 2:3), '-sum');
                        calcNnzs = strcmp(propNames(j, 2:3), '-nnz');
                        s.subs(calcSums | calcNnzs) = {':'};
                    end
                    if iscell(states.(stateName).(propNames{j, 1}))
                        s.subs{3} = idxs;
                        states.(stateName).(propNames{j, 1}) = subsref(data, s);
                    else
                        for k = 1:numel(idxs)
                            s.subs{3} = min(finTime, idxs(k)+(0:downsampleStepSec-1));
                            tmp = mean(subsref(data, s), 3);
                            
                            if calcMargins(1)
                                siz = size(tmp);
                                siz = [size(s.subs{1}) siz(2:end)];
                                tmp = permute(sum(reshape(tmp, siz), 2), [1 3 4 2]);
                            elseif calcSums(1)
                                tmp = sum(tmp, 1);
                            elseif calcNnzs(1)
                                if iscell(tmp)
                                    tmp = ~cellfun(@isempty, tmp);
                                end
                                tmp = sum(tmp ~= 0, 1);
                            end
                            if calcMargins(2)
                                siz = size(tmp);
                                siz = [siz(1) size(s.subs{2}) siz(3:end)];
                                tmp = permute(sum(reshape(tmp, siz), 3), [1 2 4 3]);
                            elseif calcSums(2)
                                tmp = sum(tmp, 2);
                            elseif calcNnzs(2)
                                if iscell(tmp)
                                    tmp = ~cellfun(@isempty, tmp);
                                end
                                tmp = sum(tmp ~= 0, 2);
                            end
                            
                            states.(stateName).(propNames{j, 1})(:, :, k) = tmp;
                            
                            clear tmp;
                        end
                    end
                    
                    clear data;
                end
            end
        end
        
        function states = downSample_extract_byTime(states, options, metadata, outputDirectory, stateNames, initTime, finTime, downsampleStepSec)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            numSegments = ceil(metadata.lengthSec / options.stepSizeSec / metadata.segmentSizeStep);
            
            times = ceil(initTime/options.stepSizeSec : downsampleStepSec/options.stepSizeSec : ceil(finTime/downsampleStepSec)*downsampleStepSec/options.stepSizeSec) * options.stepSizeSec;
            times(end) = finTime;
            
            tmp2 = cell(numSegments, size(stateNames, 1));
            for i = 1:numSegments
                if ~any(ismembc(times, (i-1) * metadata.segmentSizeStep + (1:metadata.segmentSizeStep)))
                    continue;
                end
                try
                    logI = DiskLogger.loadStructFromDisk(sprintf('%s/state-%d.mat', outputDirectory, i), [stateNames(:, 1); 'Time']);
                catch exception
                    exception.addCause(MException('DiskLogger:error', ...
                        'Unable to load segment %d of states {%s}', ...
                        i, strjoin(', ', stateNames{:, 1}))).rethrow();
                end
                allTimes = all(ismembc(logI.Time.values, times));
                segIdxs = ismembc2(times, logI.Time.values);
                gblIdxs = find(segIdxs > 0);
                segIdxs = segIdxs(gblIdxs);
                for j = 1:size(stateNames, 1)
                    try
                        calcMargins = false(1, 2);
                        calcSums = false(1, 2);
                        calcNnzs = false(1, 2);
                        s.type = '()';
                        s.subs = {':' ':' ':'};
                        if size(stateNames, 2) > 2
                            s.subs(1:2) = stateNames(j, 3:4);
                            calcMargins = cellfun(@(x) isnumeric(x) && size(x, 2) > 1, stateNames(j, 3:4));
                            calcSums = strcmp(stateNames(j, 3:4), '-sum');
                            calcNnzs = strcmp(stateNames(j, 3:4), '-nnz');
                            s.subs(calcSums | calcNnzs) = {':'};
                        end
                        
                        if ~allTimes
                            s.subs{3} = segIdxs;
                        end
                        tmp = subsref(logI.(stateNames{j, 1}).(stateNames{j, 2}), s);
                        
                        if calcMargins(1)
                            siz = size(tmp);
                            siz = [size(s.subs{1}) siz(2:end)];
                            tmp = permute(sum(reshape(tmp, siz), 2), [1 3 4 2]);
                        elseif calcSums(1)
                            tmp = sum(tmp, 1);
                        elseif calcNnzs(1)
                            if iscell(tmp)
                                tmp = ~cellfun(@isempty, tmp);
                            end
                            tmp = sum(tmp ~= 0, 1);
                        end
                        if calcMargins(2)
                            siz = size(tmp);
                            siz = [siz(1) size(s.subs{2}) siz(3:end)];
                            tmp = permute(sum(reshape(tmp, siz), 3), [1 2 4 3]);
                        elseif calcSums(2)
                            tmp = sum(tmp, 2);
                        elseif calcNnzs(2)
                            if iscell(tmp)
                                tmp = ~cellfun(@isempty, tmp);
                            end
                            tmp = sum(tmp ~= 0, 2);
                        end
                        
                        if isa(states.(stateNames{j, 1}).(stateNames{j, 2}), 'edu.stanford.covert.util.SparseMat')
                            if ~isa(tmp, 'edu.stanford.covert.util.SparseMat')
                                tmp = edu.stanford.covert.util.SparseMat(tmp);
                            end
                            tmp2{i, j} = tmp;
                        else
                            states.(stateNames{j, 1}).(stateNames{j, 2})(1:size(tmp, 1), 1:size(tmp, 2), gblIdxs) = tmp;
                        end
                    catch exception
                        exception.addCause(MException('DiskLogger:error', ...
                            'Unable to load segment %d of state %s.%s', ...
                            i, stateNames{j, 1}, stateNames{j, 2})).rethrow();
                    end
                    
                    clear tmp;
                end
                
                clear logI;
            end
            
            for j = 1:size(stateNames, 1)
                if isa(states.(stateNames{j, 1}).(stateNames{j, 2}), 'edu.stanford.covert.util.SparseMat')
                    tfs = ~cellfun(@isempty, tmp2(:, j));
                    tmp = tmp2(tfs, j);
                    
                    siz = [0 0];
                    for i = 1:numel(tmp)
                        siz(1) = max(siz(1), size(tmp{i}, 1));
                        siz(2) = max(siz(2), size(tmp{i}, 2));
                    end
                    for i = 1:numel(tmp)
                        if size(tmp{i}, 1) ~= siz(1)
                            tmp{i} = [tmp{i}; edu.stanford.covert.util.SparseMat([], [], [siz(1)-size(tmp{i}, 1) size(tmp{i}, 2) size(tmp{i}, 3)])];
                        end
                        if size(tmp{i}, 2) ~= siz(2)
                            tmp{i} = [tmp{i}  edu.stanford.covert.util.SparseMat([], [], [size(tmp{i}, 1) siz(2)-size(tmp{i}, 2) size(tmp{i}, 3)])];
                        end
                    end
                    
                    states.(stateNames{j, 1}).(stateNames{j, 2}) = cat(3, tmp{:});
                    
                    clear tmp;
                end
            end
        end
        
        function states = downSample_mean_byTime(states, options, metadata, outputDirectory, stateNames, initTime, finTime, downsampleStepSec)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            numSegments = ceil(metadata.lengthSec / options.stepSizeSec / metadata.segmentSizeStep);
            
            times = ceil(initTime/options.stepSizeSec : downsampleStepSec/options.stepSizeSec : ceil(finTime/downsampleStepSec)*downsampleStepSec/options.stepSizeSec) * options.stepSizeSec;
            
            log2 = DiskLogger.loadStructFromDisk(sprintf('%s/state-%d.mat', outputDirectory, 1), [stateNames(:,1); 'Time']);
            for i = 1:numSegments
                log1 = log2;
                if i < numSegments
                    try
                        log2 = DiskLogger.loadStructFromDisk(sprintf('%s/state-%d.mat', outputDirectory, i+1), [stateNames(:,1); 'Time']);
                    catch exception
                        exception.addCause(MException('DiskLogger:error', ...
                            'Unable to load segment %d of states {%s}', ...
                            i, strjoin(', ', stateNames{:, 1}))).rethrow();
                    end
                end
                segIdxs = ismembc2(times, log1.Time.values);
                gblIdxs = find(segIdxs > 0);
                segIdxs = segIdxs(gblIdxs);
                
                for k = 1:numel(segIdxs)
                    for j = 1:size(stateNames, 1)
                        try
                            tmpSegIdxs = segIdxs(k)+(1:downsampleStepSec/options.stepSizeSec)-1;
                            segmentLength = size(log1.(stateNames{j, 1}).(stateNames{j, 2}), 3);
                            calcMargins = false(1, 2);
                            calcSums = false(1, 2);
                            calcNnzs = false(1, 2);
                            s.type = '()';
                            s.subs = {':' ':' ':'};
                            if size(stateNames, 2) > 2
                                s.subs(1:2) = stateNames(j, 3:4);
                                calcMargins = cellfun(@(x) isnumeric(x) && size(x, 2) > 1, stateNames(j, 3:4));
                                calcSums = strcmp(stateNames(j, 3:4), '-sum');
                                calcNnzs = strcmp(stateNames(j, 3:4), '-nnz');
                                s.subs(calcSums | calcNnzs) = {':'};
                            end
                            if tmpSegIdxs(end) <= segmentLength
                                s.subs{3} = tmpSegIdxs;
                                tmp = subsref(log1.(stateNames{j, 1}).(stateNames{j, 2}), s);
                            else
                                s1 = s;
                                s2 = s;
                                s1.subs{3} = tmpSegIdxs(1):size(log1.(stateNames{j, 1}).(stateNames{j, 2}), 3);
                                s2.subs{3} = 1:tmpSegIdxs(end)-segmentLength;
                                tmp = cat(3, ...
                                    subsref(log1.(stateNames{j, 1}).(stateNames{j, 2}), s1), ...
                                    subsref(log2.(stateNames{j, 1}).(stateNames{j, 2}), s2));
                            end
                            
                            if calcMargins(1)
                                siz = size(tmp);
                                siz = [size(s.subs{1}) siz(2:end)];
                                tmp = permute(sum(reshape(tmp, siz), 2), [1 3 4 2]);
                            elseif calcSums(1)
                                tmp = sum(tmp, 1);
                            elseif calcNnzs(1)
                                if iscell(tmp)
                                    tmp = ~cellfun(@isempty, tmp);
                                end
                                tmp = sum(tmp ~= 0, 1);
                            end
                            if calcMargins(2)
                                siz = size(tmp);
                                siz = [siz(1) size(s.subs{2}) siz(3:end)];
                                tmp = permute(sum(reshape(tmp, siz), 3), [1 2 4 3]);
                            elseif calcSums(2)
                                tmp = sum(tmp, 2);
                            elseif calcNnzs(2)
                                if iscell(tmp)
                                    tmp = ~cellfun(@isempty, tmp);
                                end
                                tmp = sum(tmp ~= 0, 2);
                            end
                            
                            switch class(tmp)
                                case 'cell', tmp = tmp(:, :, 1);
                                otherwise, tmp = mean(tmp, 3);
                            end
                            
                            if isa(states.(stateNames{j, 1}).(stateNames{j, 2}), 'edu.stanford.covert.util.SparseMat')
                                if gblIdxs(k) == 1
                                    if ~isa(tmp, 'edu.stanford.covert.util.SparseMat')
                                        tmp = edu.stanford.covert.util.SparseMat(tmp);
                                    end
                                    states.(stateNames{j, 1}).(stateNames{j, 2}) = tmp;
                                else
                                    states.(stateNames{j, 1}).(stateNames{j, 2}) = cat(3, ...
                                        states.(stateNames{j, 1}).(stateNames{j, 2}), ...
                                        tmp);
                                end
                            else
                                states.(stateNames{j, 1}).(stateNames{j, 2})(1:size(tmp,1), 1:size(tmp,2), gblIdxs(k)) = tmp;
                            end
                        catch exception
                            exception.addCause(MException('DiskLogger:error', ...
                                'Unable to load segment %d %s.%s', ...
                                i, stateNames{j, 1}, stateNames{j, 2})).rethrow();
                        end
                        
                        clear tmp;
                    end
                end
                
                clear log1;
            end
            
            clear log2;
        end
        
        function saveStructToDisk(outputDirectory, path, data, structTF, appendTf)
            options = {'-v7'};
            if nargin >= 5 && appendTf
                options = [options; '-append'];
            end
            if (nargin < 4 && isstruct(data)) || (nargin >= 4 && structTF)
                options = [options; '-struct'];
            end
            
            lastwarn('');
            save([outputDirectory filesep path],  options{:}, 'data');
            
            [~, id] = lastwarn();
            if isequal(id, 'MATLAB:save:sizeTooBigForMATFile')
                tmp = whos('data');
                if tmp.bytes < 5e9
                    options{1} = '-v7.3';
                    save([outputDirectory filesep path],  options{:}, 'data');
                else
                    warning('WholeCell:warning', 'Data too large to store as mat file ''%s''', [outputDirectory filesep path]);
                end
            end
        end
        
        function structure = loadStructFromDisk(path, fields, varargin)
            if nargin == 1
                structure = load(path);
            elseif nargin == 2 && iscell(fields)
                fields = unique(fields);
                structure = load(path, fields{:});
            elseif nargin == 2
                structure = load(path, fields);
            elseif nargin > 2
                structure = load(path, fields, varargin{:});
            end
        end
        
        %Allocates space for one segment's worth of time course data.
        function log = allocateMemory(metaStates, segmentSizeStep)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            log = struct;
            
            for i = 1:size(metaStates.names, 1)
                log.(metaStates.names{i, 1}) = struct;
            end
            
            for i = 1:size(metaStates.properties, 1)
                log.(metaStates.properties{i, 1}).(metaStates.properties{i, 2}) = DiskLogger.allocateData(...
                    metaStates.properties{i, 3}, [metaStates.properties{i, 4}, segmentSizeStep]);
            end
        end
        
        function data = allocateData(dataType, sz, sparsity)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.util.SparseMat;
            
            if prod(sz) > 10^8 && (nargin < 3 || sparsity < 1/3) && ...
                    ismember(dataType, {'logical', 'double', 'uint8', 'uint16', 'uint32', 'uint64', 'int8', 'int16' ,'int32', 'int64'})                    
                dataType = 'edu.stanford.covert.util.SparseMat';
            end
            
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
                        exception.addCause(MException('DiskLogger:invalidDataType', '%s is not a numeric data type', dataType)).rethrow();
                    end
            end
        end
    end
    
    methods (Static = true)
        function reindexTimeCourses(outputDirectory)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            %% get meta data, options
            metadata = DiskLogger.loadMetadata(outputDirectory);
            options = DiskLogger.loadOptions(outputDirectory);
            stateNames = DiskLogger.getAvailableStates(outputDirectory);
            
            %% re-index simulation, organized by state
            for i = 1:size(stateNames, 1)
                try
                    data = DiskLogger.loadTimecourses(outputDirectory, ...
                        stateNames(i, :), ...
                        1, metadata.lengthSec, metadata.downsampleStepSec, 'extract', ...
                        options, metadata, 'byTime');
                    DiskLogger.saveStructToDisk(outputDirectory, ['state-' stateNames{i, 1} '-' stateNames{i, 2} '.mat'], ...
                        data.(stateNames{i, 1}).(stateNames{i, 2}), false);                    
                catch exception
                    if ismember(exception.identifier, {'MATLAB:save:errorClosingFile'; 'MATLAB:nomem'})
                        warning('WholeCell:warning:reindexing', 'Unable to reindex %s.%s\n%s', stateNames{i, 1}, stateNames{i, 2}, exception.getReport());
                    else
                        exception.addCause(MException('DiskLogger:error', 'Unable to reindex %s.%s', stateNames{i, 1}, stateNames{i, 2})).rethrow();
                    end
                end
                
                %clean up
                clear data;
            end
        end
    end
end
