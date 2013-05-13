% Whole cell simulation class.
% - Runs simulations
% - Downsamples simulations
% - Stores and loads simulation to/from local disk or a database
% - Simulations are stored both on local disk and in the database with the
%   meta data, parameters, and time courses list below. Several dependent
%   time courses are also stored in the database.
% - Opens simulations in dashboard
%
% Options
% - lengthSec
% - stepSizeSec
% - verbosity
% - seed
% - geneticKnockouts
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/7/2010
classdef Simulation < handle
    properties (Constant)
        optionNames = {
            'lengthSec'
            'stepSizeSec'
            'verbosity'
            'seed'
            'macromoleculeStateInitialization'
            'geneticKnockouts'
            'stimulus'
            'media'
            };
    end
    
    properties (SetAccess = private)
        compartment             %compartments
        gene                    %indices of genes
        
        stateMetadata           %all state IDs, classes, names, etc.
        states                  %cell array of states in alphabetical order
        state_time              %reference to time state
        state_metabolite        %reference to metabolite state
        state_stimulus          %reference to stimulus state
        
        processMetadata         %all process IDs, classes, names, etc.
        processes               %cell array of processes in alphabetical order
        processInitOrderIndexs  %indices of processes in initialization order
        processEvalOrderIndexs  %indices of processes in state evolution order
        processesInInitOrder    %processes in initialization order
        processesInEvalOrder    %processes in state evolution order
        processIndex_translation        %index of translation process
        processIndex_tRNAAminoacylation %index of tRNA aminoacylation process
    end
    
    properties (Access = private)
        randStream              %random number stream
    end
    
    %options. Use applyOptions(name1, value1, ...) to change one or more
    %options before a run.
    properties (SetAccess = private)
        lengthSec   = 50000; %length of simulation (s); ~1.5x 9 h doubling time
        stepSizeSec = 1;     %time scale of simulation (s)
        verbosity   = 1;     %0 = no output, 5 = maximum output
        seed        = [];    %set to any number for reproducible random streams
        macromoleculeStateInitialization = 'multinomial';  %Toggles how inital state is calculated; see initializeState
        
        geneticKnockouts = {}; %Whole Cell Model IDs of genes to be knocked out
        stimulus         = []; %stimulus set values
        media            = []; %media set values
    end
        
    methods
        function this = Simulation(kbStates, kbProcesses)
            %construct components
            this.constructRandStream();
            
            %construct constants
            this.compartment = edu.stanford.covert.cell.sim.constant.Compartment();
            this.gene = edu.stanford.covert.cell.sim.constant.Gene();
            
            %construct states
            this.states = {};
            if exist('kbStates', 'var')
                sMetadata.wholeCellModelIDs = {kbStates.wholeCellModelID}';
                sMetadata.names             = {kbStates.name}';
                sMetadata.classes           = {kbStates.class}';
                
                this.stateMetadata = sMetadata;
                this.stateMetadata.indexs = ...
                    this.buildNameIndex(sMetadata.wholeCellModelIDs);
                this.constructStates(sMetadata.wholeCellModelIDs);
            end
            
            %construct processes
            this.processes = {};
            if exist('kbProcesses', 'var')
                mMetadata.wholeCellModelIDs = {kbProcesses.wholeCellModelID}';
                mMetadata.names             = {kbProcesses.name}';
                mMetadata.classes           = {kbProcesses.class}';
                mMetadata.initOrder         = [kbProcesses.initializationOrder]';
                mMetadata.evalOrder         = [kbProcesses.evaluationOrder]';
                
                this.processMetadata = mMetadata;
                this.processMetadata.indexs = ...
                    this.buildNameIndex(mMetadata.wholeCellModelIDs);
                
                this.constructProcesses(mMetadata.wholeCellModelIDs);
            end
            
            %link processes-states
            for i = 1:numel(this.states)
                this.states{i}.storeObjectReferences(this);
            end
            
            %link processes-states
            for i = 1:numel(this.processes)
                this.processes{i}.storeObjectReferences(this);
            end
        end
    end
    
    methods
        function state = state(this, wholeCellModelID)
            i = this.stateIndex(wholeCellModelID);
            if i
                state = this.states{i};
            else
                state = [];
            end
        end
        
        function process = process(this, wholeCellModelID)
            i = this.processIndex(wholeCellModelID);
            if i
                process = this.processes{i};
            else
                process = [];
            end
        end
        
        function i = stateIndex(this, wholeCellModelID)
            [~, i] = ismember(wholeCellModelID, this.stateWholeCellModelIDs);
            if i == 0
                [~, i] = ismember(['State_' wholeCellModelID], this.stateWholeCellModelIDs);
            end
        end
        
        function i = processIndex(this, wholeCellModelID)
            [~,i] = ismember(wholeCellModelID, this.processWholeCellModelIDs);
            if i == 0
                [~, i] = ismember(['Process_' wholeCellModelID], this.processWholeCellModelIDs);
            end
        end
        
        function result = stateWholeCellModelIDs(this)
            result = cellfun(@(s) s.wholeCellModelID, this.states, 'UniformOutput', false);
        end
        
        function result = processWholeCellModelIDs(this)
            result = cellfun(@(m) m.wholeCellModelID, this.processes, 'UniformOutput', false);
        end
        
        function this = applyAllParameters(this, varargin)
            if nargin >= 2 && isstruct(varargin{1})
                values = varargin{1};
            else
                values = cell2struct(varargin(2:2:end), varargin(1:2:end-1), 2);
            end
            
            this.applyOptions(values);
            this.applyParameters(values);
            this.applyFittedConstants(values);
            this.applyFixedConstants(values);
        end
        
        function this = applyOptions(this, varargin)
            if nargin >= 2 && isstruct(varargin{1})
                options = varargin{1};
            else
                options = cell2struct(varargin(2:2:end), varargin(1:2:end-1), 2);
            end
            this.setFromStruct(options, 'options', this.optionNames);
            if mod(this.lengthSec, this.stepSizeSec)
                throw(MException('Simulation:error',...
                    sprintf('Simulation length (%d) must be a multiple of step size (%d)',...
                    this.lengthSec, this.stepSizeSec)));
            end
        end
        
        function this = applyParameters(this, parameters)
            this.setFromStruct(parameters, 'parameters', {});
        end
        
        %apply perturbations
        %- genetic: adjust constants so that effectively the genes'
        %  products won't be expressed. This is done by setting the gene
        %  products' decay rates to be very high
        function this = applyPerturbations(this)
            %modify constants
            this.applyPerturbationsToConstants();
            
            %modify initial conditions
            this.applyPerturbationsToState()
        end
        
        function applyPerturbationsToConstants(this)
            %% genetic knockouts
            %references
            g = this.gene;
            r = this.state('Rna');
            m = this.state('ProteinMonomer');
            c = this.state('ProteinComplex');
            
            %indices
            [~, geneIdxs] = ismember(this.geneticKnockouts, this.gene.wholeCellModelIDs);
            if ~all(geneIdxs)
                throw(MException('Simulation:error', 'undefined genes'));
            end
            matureRNAIdxs = find(any(r.matureRNAGeneComposition(geneIdxs(~ismember(geneIdxs, g.mRNAIndexs)), :), 1));
            monIdxs = find(ismember(this.gene.mRNAIndexs, geneIdxs));
            cpxIdxs = find(any(any(c.proteinComplexComposition(geneIdxs, :, :), 3), 1));
            
            %set constants to that effectively the knocked out genes aren't
            %expressed. More specifically we allow genes to be expressed,
            %but then decay their gene products quickly.
            decayRate = 1e6;
            
            r.decayRates(r.processedIndexs(matureRNAIdxs), :)     = decayRate;
            r.decayRates(r.matureIndexs(matureRNAIdxs), :)        = decayRate;
            r.decayRates(r.boundIndexs(matureRNAIdxs), :)         = decayRate;
            r.decayRates(r.misfoldedIndexs(matureRNAIdxs), :)     = decayRate;
            r.decayRates(r.damagedIndexs(matureRNAIdxs), :)       = decayRate;
            r.decayRates(r.aminoacylatedIndexs(matureRNAIdxs), :) = decayRate;
            
            m.decayRates(m.nascentIndexs(monIdxs), :)             = decayRate;
            m.decayRates(m.signalSequenceIndexs(monIdxs), :)      = decayRate;
            m.decayRates(m.processedIIndexs(monIdxs), :)          = decayRate;
            m.decayRates(m.processedIIIndexs(monIdxs), :)         = decayRate;
            m.decayRates(m.foldedIndexs(monIdxs), :)              = decayRate;
            m.decayRates(m.matureIndexs(monIdxs), :)              = decayRate;
            m.decayRates(m.inactivatedIndexs(monIdxs), :)         = decayRate;
            m.decayRates(m.boundIndexs(monIdxs), :)               = decayRate;
            m.decayRates(m.misfoldedIndexs(monIdxs), :)           = decayRate;
            m.decayRates(m.damagedIndexs(monIdxs), :)             = decayRate;
            
            c.decayRates(c.nascentIndexs(cpxIdxs), :)             = decayRate;
            c.decayRates(c.matureIndexs(cpxIdxs), :)              = decayRate;
            c.decayRates(c.inactivatedIndexs(cpxIdxs), :)         = decayRate;
            c.decayRates(c.boundIndexs(cpxIdxs), :)               = decayRate;
            c.decayRates(c.misfoldedIndexs(cpxIdxs), :)           = decayRate;
            c.decayRates(c.damagedIndexs(cpxIdxs), :)             = decayRate;
        end
        
        function applyPerturbationsToState(this)
            import edu.stanford.covert.cell.sim.constant.Condition;
            
            %% genetic knockouts
            %references
            g = this.gene;
            r = this.state('Rna');
            m = this.state('ProteinMonomer');
            c = this.state('ProteinComplex');
            
            %indices
            [~, geneIdxs] = ismember(this.geneticKnockouts, this.gene.wholeCellModelIDs);
            if ~all(geneIdxs)
                throw(MException('Simulation:error', 'undefined genes'));
            end
            matureRNAIdxs = find(any(r.matureRNAGeneComposition(geneIdxs(~ismember(geneIdxs, g.mRNAIndexs)), :), 1));
            monIdxs = find(ismember(this.gene.mRNAIndexs, geneIdxs));
            cpxIdxs = find(any(any(c.proteinComplexComposition(geneIdxs, :, :), 3), 1));
            
            %modify initial conditions:
            %- set RNAs, protein counts to zero
            %- synchronize molecule counts with other states
            %  - Chromosome: bound proteins
            %  - RNA polmyerase, transcript: No active transcription
            %  - Ribosome, polypeptide: No active translation
            %  - FtsZ ring: No formed ring
            initRNACounts = r.counts;
            r.counts(r.processedIndexs(matureRNAIdxs), :)     = 0;
            r.counts(r.matureIndexs(matureRNAIdxs), :)        = 0;
            r.counts(r.boundIndexs(matureRNAIdxs), :)         = 0;
            r.counts(r.misfoldedIndexs(matureRNAIdxs), :)     = 0;
            r.counts(r.damagedIndexs(matureRNAIdxs), :)       = 0;
            r.counts(r.aminoacylatedIndexs(matureRNAIdxs), :) = 0;
            if any(any(r.updateExternalState(-(initRNACounts - r.counts), true)))
                throw(MException('Simulation:error', 'Unable to properly knockout gene'));
            end
            
            initMonCounts = m.counts;
            m.counts(m.nascentIndexs(monIdxs), :)        = 0;
            m.counts(m.signalSequenceIndexs(monIdxs), :) = 0;
            m.counts(m.processedIIndexs(monIdxs), :)     = 0;
            m.counts(m.processedIIIndexs(monIdxs), :)    = 0;
            m.counts(m.foldedIndexs(monIdxs), :)         = 0;
            m.counts(m.matureIndexs(monIdxs), :)         = 0;
            m.counts(m.inactivatedIndexs(monIdxs), :)    = 0;
            m.counts(m.boundIndexs(monIdxs), :)          = 0;
            m.counts(m.misfoldedIndexs(monIdxs), :)      = 0;
            m.counts(m.damagedIndexs(monIdxs), :)        = 0;
            if any(any(m.updateExternalState(-(initMonCounts - m.counts), true)))
                throw(MException('Simulation:error', 'Unable to properly knockout gene'));
            end
            
            initCpxCounts = c.counts;
            c.counts(c.nascentIndexs(cpxIdxs), :)     = 0;
            c.counts(c.matureIndexs(cpxIdxs), :)      = 0;
            c.counts(c.inactivatedIndexs(cpxIdxs), :) = 0;
            c.counts(c.boundIndexs(cpxIdxs), :)       = 0;
            c.counts(c.misfoldedIndexs(cpxIdxs), :)   = 0;
            c.counts(c.damagedIndexs(cpxIdxs), :)     = 0;
            if any(any(c.updateExternalState(-(initCpxCounts - c.counts), true)))
                throw(MException('Simulation:error', 'Unable to properly knockout gene'));
            end
            
            %% stimulus
            stim = this.state('Stimulus');
            newVals = this.stimulus;
            vals = stim.setValues;
            tfs = true(size(vals, 1), 1);
            for i = 1:size(newVals, 1)
                tfs = tfs & ~any(...
                    ismember(vals(:, Condition.objectCompartmentIndexs), newVals(:, Condition.objectCompartmentIndexs)) & ...
                    ((vals(:, Condition.initialTimeIndexs) >= newVals(:, Condition.initialTimeIndexs) & ...
                    vals(:, Condition.initialTimeIndexs) <= newVals(:, Condition.finalTimeIndexs)) | ...
                    (vals(:, Condition.finalTimeIndexs) >= newVals(:, Condition.initialTimeIndexs) & ...
                    vals(:, Condition.finalTimeIndexs) <= newVals(:, Condition.finalTimeIndexs)) | ...
                    (vals(:, Condition.finalTimeIndexs) <= newVals(:, Condition.initialTimeIndexs) & ...
                    vals(:, Condition.finalTimeIndexs) >= newVals(:, Condition.finalTimeIndexs))));
            end
            stim.setValues = [vals(tfs, :); newVals];
            
            %% media
            m = this.state('Metabolite');
            newVals = this.media;
            vals = m.setCounts;
            tfs = true(size(vals, 1), 1);
            for i = 1:size(newVals, 1)
                tfs = tfs & ~any(...
                    ismember(vals(:, Condition.objectCompartmentIndexs), newVals(:, Condition.objectCompartmentIndexs)) & ...
                    ((vals(:, Condition.initialTimeIndexs) >= newVals(:, Condition.initialTimeIndexs) & ...
                    vals(:, Condition.initialTimeIndexs) <= newVals(:, Condition.finalTimeIndexs)) | ...
                    (vals(:, Condition.finalTimeIndexs) >= newVals(:, Condition.initialTimeIndexs) & ...
                    vals(:, Condition.finalTimeIndexs) <= newVals(:, Condition.finalTimeIndexs)) | ...
                    (vals(:, Condition.finalTimeIndexs) <= newVals(:, Condition.initialTimeIndexs) & ...
                    vals(:, Condition.finalTimeIndexs) >= newVals(:, Condition.finalTimeIndexs))));
            end
            m.setCounts = [vals(tfs, :); newVals];
        end
        
        function this = applyFittedConstants(this, value)
            this.setFromStruct(value, 'fittedConstants', {});
        end
        
        function this = applyFixedConstants(this, value)
            this.setFromStruct(value, 'fixedConstants', {});
        end
        
        function this = applyRandStreamStates(this, value)
            if size(value.simulation(:, 1), 2) > 1
                warning('WholeCell:warning', 'Applying first rand stream state');
            end
            
            this.randStream.state = value.simulation(:, 1);
            for i = 1:numel(this.states)
                o = this.states{i};
                o.randStream.state = value.states.(o.wholeCellModelID(7:end))(:, 1);
            end
            for i = 1:numel(this.processes)
                o = this.processes{i};
                o.randStream.state = value.processes.(o.wholeCellModelID(9:end))(:, 1);
            end
        end
        
        function this = applyTimeCourses(this, value)
            this.setFromStruct(value, 'timeCourses', {});
        end
        
        function constructRandStream(this)
            %seed private stream
            this.randStream = edu.stanford.covert.util.RandStream('mcg16807');
            
            %seed state private streams
            for i = 1:length(this.states)
                o = this.states{i};
                o.constructRandStream();
            end
            
            %seed process private streams
            for i = 1:length(this.processes)
                o = this.processes{i};
                o.constructRandStream();
            end
        end
        
        function this = seedRandStream(this)
            %pick a seed
            if isempty(this.seed)
                this.seed = round(mod(now, 1) * 1e7);
            end
            
            %seed default stream
            if verLessThan('matlab', '7.12')
                globalStream = RandStream.getDefaultStream(); %#ok<GETRS>
            else
                globalStream = RandStream.getGlobalStream();
            end
            globalStream.reset(1);
            
            %seed private stream
            this.randStream.reset(this.seed)
            
            %seed state private streams
            for i = 1:length(this.states)
                o = this.states{i};
                o.seed = this.seed;
                o.seedRandStream();
            end
            
            %seed process private streams
            for i = 1:length(this.processes)
                o = this.processes{i};
                o.seed = this.seed;
                o.seedRandStream();
            end
        end
        
        function value = getAllParameters(this)
            import edu.stanford.covert.util.StructUtil;
            
            value = struct;
            value = StructUtil.catstruct(value, this.getOptions());
            value = StructUtil.catstruct(value, this.getParameters());
            value = StructUtil.catstruct(value, this.getFittedConstants());
            value = StructUtil.catstruct(value, this.getFixedConstants());
        end
        
        function value = getOptions(this)
            value = this.getAsStruct('options', this.optionNames);
        end
        
        function value = getParameters(this)
            value = this.getAsStruct('parameters', {});
        end
        
        function value = getFittedConstants(this)
            value = this.getAsStruct('fittedConstants', {});
        end
        
        function value = getFixedConstants(this)
            value = this.getAsStruct('fixedConstants', {});
        end
        
        function value = getTimeCourses(this)
            value = this.getAsStruct('timeCourses', {}, true, false);
        end
        
        function value = getRandStreamStates(this)
            value = struct('simulation', [], 'states', struct(), 'processes', struct());
            value.simulation = this.randStream.state;
            for i = 1:numel(this.states)
                o = this.states{i};
                value.states.(o.wholeCellModelID(7:end)) = o.randStream.state;
            end
            for i = 1:numel(this.processes)
                o = this.processes{i};
                value.processes.(o.wholeCellModelID(9:end)) = o.randStream.state;
            end
        end
        
        %number of simulation steps
        function value = getNumSteps(this)
            value = this.lengthSec / this.stepSizeSec;
        end
        
        % Allocate arrays to store values of physical properties over course of
        % simulation. Each of these arrays has one column.
        function this = allocateMemoryForState(this, numTimePoints)
            %states
            for i = 1:numel(this.states)
                this.states{i}.allocateMemory(numTimePoints);
            end
            
            %processes
            for i = 1:numel(this.processes)
                this.processes{i}.allocateMemoryForState(numTimePoints);
            end
        end
    end
    
    %loading, allocation, initialization
    methods
        this = initializeConstants(this, knowledgeBase)
        this = initializeState(this)
        this = run(this, varargin)
        [this, requirements, allocations, usages] = evolveState(this, outputDirectory)
        [daughters, state_daughts] = divideState(this)
    end
    
    %protected helper methods
    methods (Access = protected)
        function this = constructStates(this, stateWholeCellModelIDs)
            if ~exist('stateWholeCellModelIDs', 'var')
                stateWholeCellModelIDs = this.stateMetadata.wholeCellModelIDs;
            end
            
            this.states = cell(numel(stateWholeCellModelIDs), 1);
            [~, map1] = ismember(stateWholeCellModelIDs, this.stateMetadata.wholeCellModelIDs);
            [~, map2] = ismember(cellfun(@(x) ['State_' x], stateWholeCellModelIDs, 'UniformOutput', false), this.stateMetadata.wholeCellModelIDs);
            map = map1 + map2;
            
            %construct states
            for i = 1:numel(stateWholeCellModelIDs)
                constructor = str2func(...
                    ['edu.stanford.covert.cell.sim.state.' this.stateMetadata.classes{map(i)}]);
                this.states{i} = constructor(...
                    this.stateMetadata.wholeCellModelIDs{map(i)}, this.stateMetadata.names{map(i)});
            end
            
            %construct references to states
            this.state_time       = this.state('Time');
            this.state_metabolite = this.state('Metabolite');
            this.state_stimulus   = this.state('Stimulus');
        end
        
        function this = constructProcesses(this, processWholeCellModelIDs)
            if ~exist('processWholeCellModelIDs', 'var')
                processWholeCellModelIDs = this.processMetadata.wholeCellModelIDs;
            end
            
            [~, map1] = ismember(processWholeCellModelIDs, this.processMetadata.wholeCellModelIDs);
            [~, map2] = ismember(cellfun(@(x) ['Process_' x], processWholeCellModelIDs, 'UniformOutput', false), this.processMetadata.wholeCellModelIDs);
            map = map1 + map2;
            
            if ~isempty(map)
                this.processes{length(map),1} = [];
                for i = 1:length(map)
                    constructor = str2func(...
                        ['edu.stanford.covert.cell.sim.process.' this.processMetadata.classes{map(i)}]);
                    this.processes{i} = constructor(...
                        this.processMetadata.wholeCellModelIDs{map(i)}, this.processMetadata.names{map(i)});
                end
            else
                this.processes(:) = [];
            end
            
            [this.processesInInitOrder, this.processInitOrderIndexs] = this.createProcessOrdering(this.processMetadata.initOrder);
            [this.processesInEvalOrder, this.processEvalOrderIndexs] = this.createProcessOrdering(this.processMetadata.evalOrder);
            
            this.processIndex_translation = this.processIndex('Translation');
            this.processIndex_tRNAAminoacylation = this.processIndex('tRNAAminoacylation');
        end
        
        function [processes, indexs] = createProcessOrdering(this, completeOrdering)
            ids = this.processWholeCellModelIDs;
            [~,order] = sort(completeOrdering);
            indexs = zeros(length(this.processes),1,'int8');
            j = 1;
            for i = 1:length(order)
                if completeOrdering(order(i))
                    [~,k] = ismember(this.processMetadata.wholeCellModelIDs(order(i)), ids);
                    if k
                        indexs(j) = k;
                        j = j + 1;
                    end
                end
            end
            indexs = indexs(1:j-1);
            
            processes = this.processes(indexs);
        end
        
        function result = buildNameIndex(~, names)
            result = struct;
            for i = 1:length(names)
                result.(names{i}) = i;
            end
        end
        
        function value = getAsStruct(this, field, fields, statesFlag, processesFlag)
            if nargin < 4
                statesFlag = true;
            end
            if nargin < 5
                processesFlag = true;
            end
            
            value = struct();
            if statesFlag
                value.states = struct;
            end
            if processesFlag
                value.processes = struct;
            end
            
            propertyNames = fieldnames(this);
            for i = 1:length(fields)
                if ismember(fields{i}, propertyNames)
                    value.(fields{i}) = this.(fields{i});
                end
            end
            
            if statesFlag
                for i = 1:length(this.states)
                    id = strrep(this.states{i}.wholeCellModelID, 'State_', '');
                    value.states.(id) = struct;
                    mValue = this.states{i}.(['get' upper(field(1)) field(2:end)]);
                    if ~isempty(fieldnames(mValue))
                        value.states.(id) = mValue;
                    end
                end
            end
            
            if processesFlag
                for i = 1:length(this.processes)
                    id = strrep(this.processes{i}.wholeCellModelID, 'Process_', '');
                    value.processes.(id) = struct;
                    mValue = this.processes{i}.(['get' upper(field(1)) field(2:end)]);
                    if ~isempty(fieldnames(mValue))
                        value.processes.(id) = mValue;
                    end
                end
            end
        end
        
        function this = setFromStruct(this, value, field, fields)
            propertyNames = fieldnames(this);
            for i = 1:length(fields)
                if isfield(value, fields{i}) && ismember(fields{i}, propertyNames)
                    try %#ok<TRYNC>
                        this.(fields{i}) = value.(fields{i});
                    end
                end
            end
            
            if isfield(value, 'states')
                for i = 1:length(this.states)
                    state = this.states{i};
                    id = strrep(state.wholeCellModelID, 'State_', '');
                    if isfield(value.states, id)
                        state.(['set' upper(field(1)) field(2:end)])(value.states.(id));
                    end
                end
            end
            
            if isfield(value, 'processes')
                for i = 1:length(this.processes)
                    process = this.processes{i};
                    id =  strrep(process.wholeCellModelID, 'Process_', '');
                    if isfield(value.processes, id)
                        process.(['set' upper(field(1)) field(2:end)])(value.processes.(id));
                    end
                end
            end
        end
    end
    
    %helper methods for testing
    methods
        function value = getForTest(this, propName)
            value = this.(propName);
        end
        
        function this = setForTest(this, propName, value)
            this.(propName) = value;
        end
    end
end
