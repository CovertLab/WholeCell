% Defines a simulation process interface:
% - ID, name, index
% - Flags for initialized
% - Options
% - Parameters
% - Locally stored state -- stimuli, substrate, enzymes, bound enzymes
%   involved in process
%   - IDs
%   - names
%   - indices
%   - mass
%   - current value (during simulation)
% - Effects of process on simulation biomass, metabolism production
% - Minimum enzyme expression consistent with growth rate
% - Contribution of process to initial state, time evolution of cell
% - Resource requirements of process
% - Independent random stream -- all generation of random numbers must
%   either be through method of randStream property, or the randomStream
%   property must be set to the default random stream prior to execution
%   and the default random stream reset after execution using
%   randStream.setDefault() and randStream.resetDefault()
%
% Provides support for:
% - Mapping network components -- stimuli, substrates, enzymes, to/from
%   simulation. Mapping is computed by initializeConstants and is used by
%   copyStateFromStimulation and copyToState
% - Initializing process constants from knowledge base
% - Retrieving/sending constrants to/from simulation
% - Retrieving/sending state to/from simulation
% - Computing gene composition of process enzymes
%
% Intialization code is executed as follows in the simulation
% initializeState method:
%   1. process.copyFromState()
%   2. process.initializeState()
%   3. process.copyToState()
%
% Time evolution code is executed as follows in the simulation evolveState
% method:
%   1. process.copyFromState()
%   2. process.evolveState()
%   3. process.copyToState()
%
% Differences between handling of stimuli, substrates, enzymes, and
% boundEnzymes:
% - All are copied from the simulation by copyFromState prior to
%   initialize and evolve state
% - Substrates are copied back to simulation by copyToState after
%   initialize and evolve state
% - Stimuli are not copied back to simulation. That is the process cannot
%   edit the stimuli of the simulation
%
% Dimensions of locally store state
%  - stimuli      length(stimuliWholeCellModelIDs)   x compartments x time
%  - substrates   length(substrateWholeCellModelIDs) x compartments x time
%  - enzymes      length(enzymeWholeCellModelIDs)    x compartments x time
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/14/2010
classdef Process < handle
    %metadata
    properties
        wholeCellModelID        %Whole Cell Model ID of process
        name                    %name of process
        initialized = false;    %whether constants of process have been initialized
    end

    properties (SetAccess = protected)
        randStream              %random number stream
    end

    %Property annotations -- subclasses should set these properties to
    %cell arrays containing names of process properties
    properties (Abstract, Constant)
        optionNames__           %names of process properties that are considered options, and can be overridden through the simulation's option property setter
        fixedConstantNames__    %names of process properties that are considered fixed constants
        fittedConstantNames__   %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        localStateNames__       %names of process properties that are part of the simulation's state and owned either by the simulation, other processes, or other properties of this process; these properties will be stored to disk/database when the simulation is stored
    end

    %Public interface to the property annotations. Base process classes add
    %properties that they define to the above private property annotations
    properties (Dependent, SetAccess = protected)
        optionNames             %names of process properties that are considered options, and can be overridden through the simulation's option property setter
        fixedConstantNames      %names of process properties that are considered fixed constants
        fittedConstantNames     %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        localStateNames         %names of process properties that are part of the simulation's state and owned either by the simulation, other processes, or other properties of this process; these properties will be stored to disk/database when the simulation is stored
    end   

    %options
    properties
        stepSizeSec = 1;        %time scale of simulation (seconds)
        verbosity   = 1;        %0 = no output, 5 = maximum output
        seed        = [];       %set to any number for reproducible random streams
    end

    %enumerations
    properties (Constant)
    end

    %IDs
    %Subclasses should set these properties to cell arrays containing the
    %Whole Cell Model IDs of stimuli, metabolites, RNAs, and proteins which
    %should be considered stimuli, substrates, or enzymes within the
    %process. initializeConstants computes a mapping between the
    %simulation's stimuli, metabolites, RNAs, and proteins properties and
    %the stimuli, substrates, and enzymes properties of the process
    properties (Abstract)
        stimuliWholeCellModelIDs     %cell array of whole cell model IDs of stimuli, metabolites, RNAs, and proteins which should be considered stimuli within the process
        substrateWholeCellModelIDs   %cell array of whole cell model IDs of stimuli, metabolites, RNAs, and proteins which should be considered substrates within the process
        enzymeWholeCellModelIDs      %cell array of whole cell model IDs of stimuli, metabolites, RNAs, and proteins which should be considered enzymes within the process
    end

    %Names of stimuli, substrates, and enzymes. Set by the
    %initializeConstants method
    properties
        stimuliNames     %cell array of stimuli names
        substrateNames   %cell array of substrate names
        enzymeNames      %cell array of enzyme names
    end

    properties
        parameterNames   %Names of properties that are considered parameters in the knowledge base and can be overriden via simulation.applyParameters
        parameterIndexs  %Names of properties that represent indices over parameter properties
    end

    %Indices which hold the mapping between the simulation's stimuli,
    %metabolites, RNAs, and proteins properties and the stimuli,
    %substrates, and enzymes properties of the process. These are set by the
    %initializeConstants method from the cell arrays above.
    %
    %GlobalIndexs properties store the index of the component within the
    %simulation's property. GlobalIndexs properties have size
    %[numComponents of particular type within simulation (eg stimuli,
    %metabolites, RNA, monomer, complex) X compartments].
    %
    %LocalIndexs properties store the index of the component within the
    %processes property. LocalIndexs properties have size
    %[numComponents of particular type within simulation (eg stimuli,
    %metabolites, RNA, monomer, complex) X 1].
    %
    %CompartmentIndexs properties store the index of the component
    %compartment within the simulation's property. CompartmentIndexs properties have size
    %[numComponents of particular type within simulation (eg stimuli,
    %metabolites, RNA, monomer, complex) X compartments]
    properties
        stimulusStimulusLocalIndexs               = zeros(0,1);
        stimulusStimulusGlobalIndexs              = zeros(0,1);
        stimulusStimulusCompartmentIndexs         = zeros(0,1);
        stimulusStimulusGlobalCompartmentIndexs   = zeros(0,1);
        stimulusMetaboliteLocalIndexs             = zeros(0,1);
        stimulusMetaboliteGlobalIndexs            = zeros(0,1);
        stimulusMetaboliteCompartmentIndexs       = zeros(0,1);
        stimulusMetaboliteGlobalCompartmentIndexs = zeros(0,1);
        stimulusRNALocalIndexs                    = zeros(0,1);
        stimulusRNAGlobalIndexs                   = zeros(0,1);
        stimulusRNACompartmentIndexs              = zeros(0,1);
        stimulusRNAGlobalCompartmentIndexs        = zeros(0,1);
        stimulusMonomerLocalIndexs                = zeros(0,1);
        stimulusMonomerGlobalIndexs               = zeros(0,1);
        stimulusMonomerCompartmentIndexs          = zeros(0,1);
        stimulusMonomerGlobalCompartmentIndexs    = zeros(0,1);
        stimulusComplexLocalIndexs                = zeros(0,1);
        stimulusComplexGlobalIndexs               = zeros(0,1);
        stimulusComplexCompartmentIndexs          = zeros(0,1);
        stimulusComplexGlobalCompartmentIndexs    = zeros(0,1);
        stimulusGlobalIndexs                      = zeros(0,1);
        
        substrateStimulusLocalIndexs              = zeros(0,1);
        substrateStimulusGlobalIndexs             = zeros(0,1);
        substrateStimulusCompartmentIndexs        = zeros(0,1);
        substrateStimulusGlobalCompartmentIndexs  = zeros(0,1);
        substrateMetaboliteLocalIndexs            = zeros(0,1);
        substrateMetaboliteGlobalIndexs           = zeros(0,1);
        substrateMetaboliteCompartmentIndexs      = zeros(0,1);
        substrateMetaboliteGlobalCompartmentIndexs= zeros(0,1);
        substrateRNALocalIndexs                   = zeros(0,1);
        substrateRNAGlobalIndexs                  = zeros(0,1);
        substrateRNACompartmentIndexs             = zeros(0,1);
        substrateRNAGlobalCompartmentIndexs       = zeros(0,1);
        substrateMonomerLocalIndexs               = zeros(0,1);
        substrateMonomerGlobalIndexs              = zeros(0,1);
        substrateMonomerCompartmentIndexs         = zeros(0,1);
        substrateMonomerGlobalCompartmentIndexs   = zeros(0,1);
        substrateComplexLocalIndexs               = zeros(0,1);
        substrateComplexGlobalIndexs              = zeros(0,1);
        substrateComplexCompartmentIndexs         = zeros(0,1);
        substrateComplexGlobalCompartmentIndexs   = zeros(0,1);
        substrateGlobalIndexs                     = zeros(0,1);
        
        enzymeStimulusLocalIndexs                 = zeros(0,1);
        enzymeStimulusGlobalIndexs                = zeros(0,1);
        enzymeStimulusCompartmentIndexs           = zeros(0,1);
        enzymeStimulusGlobalCompartmentIndexs     = zeros(0,1);
        enzymeMetaboliteLocalIndexs               = zeros(0,1);
        enzymeMetaboliteGlobalIndexs              = zeros(0,1);
        enzymeMetaboliteCompartmentIndexs         = zeros(0,1);
        enzymeMetaboliteGlobalCompartmentIndexs   = zeros(0,1);
        enzymeRNALocalIndexs                      = zeros(0,1);
        enzymeRNAGlobalIndexs                     = zeros(0,1);
        enzymeRNACompartmentIndexs                = zeros(0,1);
        enzymeRNAGlobalCompartmentIndexs          = zeros(0,1);
        enzymeMonomerLocalIndexs                  = zeros(0,1);
        enzymeMonomerGlobalIndexs                 = zeros(0,1);
        enzymeMonomerCompartmentIndexs            = zeros(0,1);
        enzymeMonomerGlobalCompartmentIndexs      = zeros(0,1);
        enzymeComplexLocalIndexs                  = zeros(0,1);
        enzymeComplexGlobalIndexs                 = zeros(0,1);
        enzymeComplexCompartmentIndexs            = zeros(0,1);
        enzymeComplexGlobalCompartmentIndexs      = zeros(0,1);
        enzymeBoundRNAGlobalCompartmentIndexs     = zeros(0,1);
        enzymeBoundMonomerGlobalCompartmentIndexs = zeros(0,1);
        enzymeBoundComplexGlobalCompartmentIndexs = zeros(0,1);
        enzymeGlobalIndexs                        = zeros(0,1);
    end

    %Indices of the global compartments of stimuli, substrates, and
    %enzymes. Computed by getter methods from the compartment indices in
    %the above section. These properties have size [numComponents of particular type within process (eg stimuli,
    %substrates, enzymes) X compartments]
    properties (Dependent)
        stimuliCompartments
        substrateCompartments
        enzymeCompartments
    end

    %constants
    properties
        substrateMolecularWeights  %molecular weights of substrates, set by initialize constants
        enzymeMolecularWeights     %molecular weights of enzymes, set by initialize constants
    end

    %global state (stored locally for convenience)
    properties
        compartment
        gene
        states
    end
    
    %global state side effects (changes that simulation should make to parts of
    %its state that the process doesn't know about/listen to/care about)
    properties
        simulationStateSideEffects %list of changes simulation should make to its state, SimulationStateSideEffect array
    end

    %local state copied from simulation
    properties
        stimuli                   %values of stimuli, updated before each evolveState call by copyFromState, size: [length(stimuliWholeCellModelIDs) X compartments X time]
        substrates                %values of substrates available to all processes, updated before each evolveState call by copyFromState, size: [length(substrateWholeCellModelIDs) X compartments X time]
        enzymes                   %values of enzymes, updated before each evolveState call by copyFromState, size: [length(enzymeWholeCellModelIDs) X compartments X time]
        boundEnzymes              %values of boundEnzymes, updated before each evolveState call by copyFromState, size: [length(enzymeWholeCellModelIDs) X compartments X time]
    end

    %dependent local state (implemented as dependent property for convenience)
    properties (Dependent)
        dryWeight            %dry weight represented by state variables accessible to process including stimuli, substrates, and enzymes, subclasses may need to override the getDryWeight method which calculates this        
    end
    
    %links to other parts of cell state
    properties
        geometry
        stimulus
        metabolite
        rna
        monomer
        complex
    end

    %constructor
    methods
        %sets process meta data
        function this = Process(wholeCellModelID, name)
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

    %communication between process/simulation
    methods
        function storeObjectReferences(this, simulation)
            this.compartment = simulation.compartment;
            this.gene = simulation.gene;
            
            this.geometry = simulation.state('Geometry');
            this.stimulus = simulation.state('Stimulus');
            this.metabolite = simulation.state('Metabolite');
            this.rna = simulation.state('Rna');
            this.monomer = simulation.state('ProteinMonomer');
            this.complex = simulation.state('ProteinComplex');
            this.states = {this.geometry; this.metabolite; this.rna; this.monomer; this.complex};
        end
        
        %Initializes process's constants from knowledge base and simulation
        %including
        %- computes mapping between process's locally stored state (stimuli,
        %  substrates, enzymes) and the simulation's state using the
        %  stimiliWholeCellModelIDs, substrateWholeCellModelIDs, and
        %  enzymeWholeCellModelIDs properties that subclasses must define
        %- fetches molecular weights of substrates and enzymes from
        %  simulation
        %- set values of parameter properties of process
        %- copies option values from simulation:
        %  - stepSizeSec
        %  - verbosity
        %  - seed
        %
        %Sets initialized property to true
        %
        %By default sets the mapping between the process's locally stored
        %state and the simulation state to only copy over 1 compartment for
        %each component. Subclasses can set the retainStimuliCompartments,
        %retainSubstrateCompartments, retainEnzymeCompartments options to
        %true to setup the mapping to copy over all compartments of a
        %component to the simulation's locally stored state.
        function initializeConstants(this, knowledgeBase, simulation, options)
            %default options
            if ~exist('options','var') || ~isstruct(options)
                options = struct;
            end
            if ~isfield(options, 'retainStimuliCompartments')
                options.retainStimuliCompartments = false;
            end
            if ~isfield(options, 'retainSubstrateCompartments')
                options.retainSubstrateCompartments = false;
            end
            if ~isfield(options, 'retainEnzymeCompartments')
                options.retainEnzymeCompartments = false;
            end

            %process
            process = findobj(knowledgeBase.processes, 'wholeCellModelID', this.wholeCellModelID);

            %set values of parameter properties
            %store annotation of parameter properties
            this.parameterNames = {process.parameters.name}';
            this.parameterIndexs = {process.parameters.index}';

            parameters = struct;
            validNames = fieldnames(this);
            for i = 1:length(process.parameters)
                p = process.parameters(i);
                if ~isempty(p.index)
                    if ismember(p.index, validNames)
                        parameters.(p.name)(this.(p.index),1) = p.defaultValue;
                    end
                else
                    parameters.(p.name) = p.defaultValue;
                end
            end
            this.setParameters(parameters);

            %copy options from simulation
            this.stepSizeSec = simulation.stepSizeSec;
            this.verbosity   = simulation.verbosity;
            this.seed        = simulation.seed;

            %- compute mapping of stimuli to simulation's state
            %- fetch names of stimuli from simulation
            %- fetch molecular weights of stimuli from simulation
            [this.stimuliNames...
             this.stimulusStimulusLocalIndexs...
             this.stimulusStimulusGlobalIndexs...
             this.stimulusStimulusCompartmentIndexs...
             this.stimulusStimulusGlobalCompartmentIndexs...
             this.stimulusMetaboliteLocalIndexs...
             this.stimulusMetaboliteGlobalIndexs...
             this.stimulusMetaboliteCompartmentIndexs...
             this.stimulusMetaboliteGlobalCompartmentIndexs...
             this.stimulusRNALocalIndexs...
             this.stimulusRNAGlobalIndexs...
             this.stimulusRNACompartmentIndexs...
             this.stimulusRNAGlobalCompartmentIndexs...
             this.stimulusMonomerLocalIndexs...
             this.stimulusMonomerGlobalIndexs...
             this.stimulusMonomerCompartmentIndexs...
             this.stimulusMonomerGlobalCompartmentIndexs...
             this.stimulusComplexLocalIndexs...
             this.stimulusComplexGlobalIndexs...
             this.stimulusComplexCompartmentIndexs...
             this.stimulusComplexGlobalCompartmentIndexs...
             ~, ~, ~, ...
             this.stimulusGlobalIndexs] = this.initializeConstantsHelper(...
                this.stimuliWholeCellModelIDs, options.retainStimuliCompartments);
           
            %- compute mapping of substrates to simulation's state
            %- fetch names of substrates from simulation
            %- fetch molecular weights of substrates from simulation
            [this.substrateNames...
             this.substrateStimulusLocalIndexs...
             this.substrateStimulusGlobalIndexs...
             this.substrateStimulusCompartmentIndexs...
             this.substrateStimulusGlobalCompartmentIndexs...
             this.substrateMetaboliteLocalIndexs...
             this.substrateMetaboliteGlobalIndexs...
             this.substrateMetaboliteCompartmentIndexs...
             this.substrateMetaboliteGlobalCompartmentIndexs...
             this.substrateRNALocalIndexs...
             this.substrateRNAGlobalIndexs...
             this.substrateRNACompartmentIndexs...
             this.substrateRNAGlobalCompartmentIndexs...
             this.substrateMonomerLocalIndexs...
             this.substrateMonomerGlobalIndexs...
             this.substrateMonomerCompartmentIndexs...
             this.substrateMonomerGlobalCompartmentIndexs...
             this.substrateComplexLocalIndexs...
             this.substrateComplexGlobalIndexs...
             this.substrateComplexCompartmentIndexs...
             this.substrateComplexGlobalCompartmentIndexs...
             ~, ~, ~, ...
             this.substrateGlobalIndexs...
             this.substrateMolecularWeights] = this.initializeConstantsHelper(...
                this.substrateWholeCellModelIDs, options.retainSubstrateCompartments);
            
            %- compute mapping of enzymes to simulation's state
            %- fetch names of enzymes from simulation
            %- fetch molecular weights of enzymes from simulation
            [this.enzymeNames...
             this.enzymeStimulusLocalIndexs...
             this.enzymeStimulusGlobalIndexs...
             this.enzymeStimulusCompartmentIndexs...
             this.enzymeStimulusGlobalCompartmentIndexs...
             this.enzymeMetaboliteLocalIndexs...
             this.enzymeMetaboliteGlobalIndexs...
             this.enzymeMetaboliteCompartmentIndexs...
             this.enzymeMetaboliteGlobalCompartmentIndexs...
             this.enzymeRNALocalIndexs...
             this.enzymeRNAGlobalIndexs...
             this.enzymeRNACompartmentIndexs...
             this.enzymeRNAGlobalCompartmentIndexs...
             this.enzymeMonomerLocalIndexs...
             this.enzymeMonomerGlobalIndexs...
             this.enzymeMonomerCompartmentIndexs...
             this.enzymeMonomerGlobalCompartmentIndexs...
             this.enzymeComplexLocalIndexs...
             this.enzymeComplexGlobalIndexs...
             this.enzymeComplexCompartmentIndexs...
             this.enzymeComplexGlobalCompartmentIndexs...
             this.enzymeBoundRNAGlobalCompartmentIndexs...
             this.enzymeBoundMonomerGlobalCompartmentIndexs...
             this.enzymeBoundComplexGlobalCompartmentIndexs...
             this.enzymeGlobalIndexs...
             this.enzymeMolecularWeights] = this.initializeConstantsHelper(...
                this.enzymeWholeCellModelIDs, options.retainEnzymeCompartments);

            this.initialized = true;
        end

        %Retrieve state from simulation/other processes
        %- cell volume
        %- stimuli
        %- substrates
        %- enzymes
        %- bound enzymes
        %
        %Uses mappings computed by initializeConstants to retrieve state
        %from simulation
        %
        %Called before initializeState, evolveState by simulation
        %initializeState and evolveState methods
        function copyFromState(this)
            this.stimuli = this.copyStimuliFromState(...
                this.stimulus.values, this.metabolite.counts, this.rna.counts, this.monomer.counts, this.complex.counts);
            this.substrates = this.copySubstratesFromState(...
                this.stimulus.values, this.metabolite.counts, this.rna.counts, this.monomer.counts, this.complex.counts);
            [this.enzymes, this.boundEnzymes] = this.copyEnzymesFromState(...
                this.stimulus.values, this.metabolite.counts, this.rna.counts, this.monomer.counts, this.complex.counts);
        end
        
        function stimuli = copyStimuliFromState(this, stimulusValues, metaboliteCounts, rnaCounts, monomerCounts, complexCounts)
            numTimePoints = size(stimulusValues, 3);
            stimuli = zeros(numel(this.stimuliWholeCellModelIDs), size(this.stimulusStimulusGlobalCompartmentIndexs, 2), numTimePoints);
            
            if isempty(stimuli)
                return;
            end
                        
            if numTimePoints == 1
                if ~isempty(this.stimulusStimulusGlobalCompartmentIndexs)
                    stimuli(this.stimulusStimulusLocalIndexs, :) = ...
                        stimulusValues(this.stimulusStimulusGlobalCompartmentIndexs);
                end
                if ~isempty(this.stimulusMetaboliteGlobalCompartmentIndexs)
                    stimuli(this.stimulusMetaboliteLocalIndexs, :) = ...
                        metaboliteCounts(this.stimulusMetaboliteGlobalCompartmentIndexs);
                end
                if ~isempty(this.stimulusRNAGlobalCompartmentIndexs)
                    stimuli(this.stimulusRNALocalIndexs, :) = ...
                        rnaCounts(this.stimulusRNAGlobalCompartmentIndexs);
                end
                if ~isempty(this.stimulusMonomerGlobalCompartmentIndexs)
                    stimuli(this.stimulusMonomerLocalIndexs, :) = ...
                        monomerCounts(this.stimulusMonomerGlobalCompartmentIndexs);
                end
                if ~isempty(this.stimulusComplexGlobalCompartmentIndexs)
                    stimuli(this.stimulusComplexLocalIndexs, :) = ...
                        complexCounts(this.stimulusComplexGlobalCompartmentIndexs);
                end
            else
                if ~isempty(this.stimulusStimulusLocalIndexs)
                    stimuli(this.stimulusStimulusLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(stimulusValues, ...
                        this.stimulusStimulusGlobalIndexs, ...
                        this.stimulusStimulusCompartmentIndexs);
                end
                if ~isempty(this.stimulusMetaboliteLocalIndexs)
                    stimuli(this.stimulusMetaboliteLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(metaboliteCounts, ...
                        this.stimulusMetaboliteGlobalIndexs,...
                        this.stimulusMetaboliteCompartmentIndexs);
                end
                if ~isempty(this.stimulusRNALocalIndexs)
                    stimuli(this.stimulusRNALocalIndexs, :, :) = ...
                        this.copyFromStateHelper(rnaCounts, ...
                        this.rna.matureIndexs(this.stimulusRNAGlobalIndexs),...
                        this.stimulusRNACompartmentIndexs);
                end
                if ~isempty(this.stimulusMonomerLocalIndexs)
                    stimuli(this.stimulusMonomerLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(monomerCounts, ...
                        this.monomer.matureIndexs(this.stimulusMonomerGlobalIndexs),...
                        this.stimulusMonomerCompartmentIndexs);
                end
                if ~isempty(this.stimulusComplexLocalIndexs)
                    stimuli(this.stimulusComplexLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(complexCounts, ...
                        this.complex.matureIndexs(this.stimulusComplexGlobalIndexs),...
                        this.stimulusComplexCompartmentIndexs);
                end
            end
        end
        
        function substrates = copySubstratesFromState(this, stimulusValues, metaboliteCounts, rnaCounts, monomerCounts, complexCounts)            
            numTimePoints = size(metaboliteCounts, 3);
            substrates = zeros(numel(this.substrateWholeCellModelIDs), size(this.substrateMetaboliteGlobalCompartmentIndexs, 2), numTimePoints);
            
            if isempty(substrates)
                return;
            end
                        
            if numTimePoints == 1
                if ~isempty(this.substrateStimulusGlobalCompartmentIndexs)
                    substrates(this.substrateStimulusLocalIndexs, :) = ...
                        stimulusValues(this.substrateStimulusGlobalCompartmentIndexs);
                end
                if ~isempty(this.substrateMetaboliteGlobalCompartmentIndexs)
                    substrates(this.substrateMetaboliteLocalIndexs, :) = ...
                        metaboliteCounts(this.substrateMetaboliteGlobalCompartmentIndexs);
                end
                if ~isempty(this.substrateRNAGlobalCompartmentIndexs)
                    substrates(this.substrateRNALocalIndexs, :) = ...
                        rnaCounts(this.substrateRNAGlobalCompartmentIndexs);
                end
                if ~isempty(this.substrateMonomerGlobalCompartmentIndexs)
                    substrates(this.substrateMonomerLocalIndexs, :) = ...
                        monomerCounts(this.substrateMonomerGlobalCompartmentIndexs);
                end
                if ~isempty(this.substrateComplexGlobalCompartmentIndexs)
                    substrates(this.substrateComplexLocalIndexs, :) = ...
                        complexCounts(this.substrateComplexGlobalCompartmentIndexs);
                end
            else
                if ~isempty(this.substrateStimulusLocalIndexs)
                    substrates(this.substrateStimulusLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(stimulusValues, ...
                        this.substrateStimulusGlobalIndexs, ...
                        this.substrateStimulusCompartmentIndexs);
                end
                if ~isempty(this.substrateMetaboliteLocalIndexs)
                    substrates(this.substrateMetaboliteLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(metaboliteCounts, ...
                        this.substrateMetaboliteGlobalIndexs,...
                        this.substrateMetaboliteCompartmentIndexs);
                end
                if ~isempty(this.substrateRNALocalIndexs)
                    substrates(this.substrateRNALocalIndexs, :, :) = ...
                        this.copyFromStateHelper(rnaCounts, ...
                        this.rna.matureIndexs(this.substrateRNAGlobalIndexs),...
                        this.substrateRNACompartmentIndexs);
                end
                if ~isempty(this.substrateMonomerLocalIndexs)
                    substrates(this.substrateMonomerLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(monomerCounts, ...
                        this.monomer.matureIndexs(this.substrateMonomerGlobalIndexs),...
                        this.substrateMonomerCompartmentIndexs);
                end
                if ~isempty(this.substrateComplexLocalIndexs)
                    substrates(this.substrateComplexLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(complexCounts, ...
                        this.complex.matureIndexs(this.substrateComplexGlobalIndexs),...
                        this.substrateComplexCompartmentIndexs);
                end
            end
        end
        
        function [enzymes, boundEnzymes] = copyEnzymesFromState(this, stimulusValues, metaboliteCounts, rnaCounts, monomerCounts, complexCounts)
            numTimePoints = size(complexCounts, 3);
            enzymes = zeros(numel(this.enzymeWholeCellModelIDs), size(this.enzymeBoundMonomerGlobalCompartmentIndexs, 2), numTimePoints);
            boundEnzymes = zeros(size(enzymes));
            
            if isempty(enzymes)
                return;
            end
            
            if numTimePoints == 1
                if ~isempty(this.enzymeStimulusGlobalCompartmentIndexs)
                    enzymes(this.enzymeStimulusLocalIndexs, :) = ...
                        stimulusValues(this.enzymeStimulusGlobalCompartmentIndexs);
                end
                if ~isempty(this.enzymeMetaboliteGlobalCompartmentIndexs)
                    enzymes(this.enzymeMetaboliteLocalIndexs, :) = ...
                        metaboliteCounts(this.enzymeMetaboliteGlobalCompartmentIndexs);
                end
                if ~isempty(this.enzymeRNAGlobalCompartmentIndexs)
                    enzymes(this.enzymeRNALocalIndexs, :) = ...
                        rnaCounts(this.enzymeRNAGlobalCompartmentIndexs);
                    boundEnzymes(this.enzymeRNALocalIndexs, :) = ...
                        rnaCounts(this.enzymeBoundRNAGlobalCompartmentIndexs);
                end
                if ~isempty(this.enzymeMonomerGlobalCompartmentIndexs)
                    enzymes(this.enzymeMonomerLocalIndexs, :) = ...
                        monomerCounts(this.enzymeMonomerGlobalCompartmentIndexs);
                    boundEnzymes(this.enzymeMonomerLocalIndexs, :) = ...
                        monomerCounts(this.enzymeBoundMonomerGlobalCompartmentIndexs);
                end
                if ~isempty(this.enzymeComplexGlobalCompartmentIndexs)
                    enzymes(this.enzymeComplexLocalIndexs, :) = ...
                        complexCounts(this.enzymeComplexGlobalCompartmentIndexs);
                    boundEnzymes(this.enzymeComplexLocalIndexs, :) = ...
                        complexCounts(this.enzymeBoundComplexGlobalCompartmentIndexs);
                end
            else
                if ~isempty(this.enzymeStimulusLocalIndexs)
                    enzymes(this.enzymeStimulusLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(stimulusValues, ...
                        this.enzymeStimulusGlobalIndexs, ...
                        this.enzymeStimulusCompartmentIndexs);
                end
                if ~isempty(this.enzymeMetaboliteLocalIndexs)
                    enzymes(this.enzymeMetaboliteLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(metaboliteCounts, ...
                        this.enzymeMetaboliteGlobalIndexs,...
                        this.enzymeMetaboliteCompartmentIndexs);
                end
                if ~isempty(this.enzymeRNALocalIndexs)
                    enzymes(this.enzymeRNALocalIndexs, :, :) = ...
                        this.copyFromStateHelper(rnaCounts, ...
                        this.rna.matureIndexs(this.enzymeRNAGlobalIndexs),...
                        this.enzymeRNACompartmentIndexs);
                end
                if ~isempty(this.enzymeMonomerLocalIndexs)
                    enzymes(this.enzymeMonomerLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(monomerCounts, ...
                        this.monomer.matureIndexs(this.enzymeMonomerGlobalIndexs),...
                        this.enzymeMonomerCompartmentIndexs);
                end
                if ~isempty(this.enzymeComplexLocalIndexs)
                    enzymes(this.enzymeComplexLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(complexCounts, ...
                        this.complex.matureIndexs(this.enzymeComplexGlobalIndexs),...
                        this.enzymeComplexCompartmentIndexs);
                end
                
                if nargin < 2
                    return;
                end
                if ~isempty(this.enzymeRNALocalIndexs)
                    boundEnzymes(this.enzymeRNALocalIndexs, :, :) = ...
                        this.copyFromStateHelper(rnaCounts, ...
                        this.rna.boundIndexs(this.enzymeRNAGlobalIndexs),...
                        this.enzymeRNACompartmentIndexs);
                end
                if ~isempty(this.enzymeMonomerLocalIndexs)
                    boundEnzymes(this.enzymeMonomerLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(monomerCounts, ...
                        this.monomer.boundIndexs(this.enzymeMonomerGlobalIndexs),...
                        this.enzymeMonomerCompartmentIndexs);
                end
                if ~isempty(this.enzymeComplexLocalIndexs)
                    boundEnzymes(this.enzymeComplexLocalIndexs, :, :) = ...
                        this.copyFromStateHelper(complexCounts, ...
                        this.complex.boundIndexs(this.enzymeComplexGlobalIndexs),...
                        this.enzymeComplexCompartmentIndexs);
                end
            end
        end

        %Send state to simulation/other processes
        %- substrates
        %- enzymes
        %- bound enzymes
        %
        %Note stimuli are not stored back to the simulation.
        %
        %Uses mappings computed by initializeConstants to send state
        %from simulation.
        %
        %Called after initializeState, evolveState by simulation
        %initializeState and evolveState methods
        %
        %Note: this function will behave incorrectly if a process touches the
        %same molecular species as both a substrate and enzyme, and modifies the
        %substrate and enzymes counts of this species. As of 1/15/2011 there
        %were no processes which did this, so the function behaved correctly.
        function copyToState(this)
            numTimePoints = size(this.substrates, 3);
            
            if ~isempty(this.enzymes)
                if numTimePoints == 1
                    if ~isempty(this.enzymeStimulusGlobalCompartmentIndexs)
                        this.stimulus.values(this.enzymeStimulusGlobalCompartmentIndexs) = ...
                            this.enzymes(this.enzymeStimulusLocalIndexs, :);
                    end
                    if ~isempty(this.enzymeMetaboliteGlobalCompartmentIndexs)
                        this.metabolite.counts(this.enzymeMetaboliteGlobalCompartmentIndexs) = ...
                            this.enzymes(this.enzymeMetaboliteLocalIndexs, :);
                    end
                    if ~isempty(this.enzymeRNAGlobalCompartmentIndexs)
                        this.rna.counts(this.enzymeRNAGlobalCompartmentIndexs) = ...
                            this.enzymes(this.enzymeRNALocalIndexs, :);
                        this.rna.counts(this.enzymeBoundRNAGlobalCompartmentIndexs) = ...
                            this.boundEnzymes(this.enzymeRNALocalIndexs, :);
                    end
                    if ~isempty(this.enzymeMonomerGlobalCompartmentIndexs)
                        this.monomer.counts(this.enzymeMonomerGlobalCompartmentIndexs) = ...
                            this.enzymes(this.enzymeMonomerLocalIndexs, :);
                        this.monomer.counts(this.enzymeBoundMonomerGlobalCompartmentIndexs) = ...
                            this.boundEnzymes(this.enzymeMonomerLocalIndexs, :);
                    end
                    if ~isempty(this.enzymeComplexGlobalCompartmentIndexs)
                        this.complex.counts(this.enzymeComplexGlobalCompartmentIndexs) = ...
                            this.enzymes(this.enzymeComplexLocalIndexs, :);
                        this.complex.counts(this.enzymeBoundComplexGlobalCompartmentIndexs) = ...
                            this.boundEnzymes(this.enzymeComplexLocalIndexs, :);
                    end
                else
                    this.copyToStateHelper(this.stimulus, ...
                        'values', ...
                        'enzymes', ...
                        this.enzymeStimulusGlobalIndexs, ...
                        this.enzymeStimulusCompartmentIndexs, ...
                        this.enzymeStimulusLocalIndexs);
                    this.copyToStateHelper(this.metabolite, ...
                        'counts', ...
                        'enzymes',...
                        this.enzymeMetaboliteGlobalIndexs,...
                        this.enzymeMetaboliteCompartmentIndexs, ...
                        this.enzymeMetaboliteLocalIndexs);
                    this.copyToStateHelper(this.rna, ...
                        'counts', ...
                        'enzymes',...
                        this.rna.matureIndexs(this.enzymeRNAGlobalIndexs),...
                        this.enzymeRNACompartmentIndexs, ...
                        this.enzymeRNALocalIndexs);
                    this.copyToStateHelper(this.monomer, ...
                        'counts', ...
                        'enzymes',...
                        this.monomer.matureIndexs(this.enzymeMonomerGlobalIndexs),...
                        this.enzymeMonomerCompartmentIndexs, ...
                        this.enzymeMonomerLocalIndexs);
                    this.copyToStateHelper(this.complex, ...
                        'counts', ...
                        'enzymes',...
                        this.complex.matureIndexs(this.enzymeComplexGlobalIndexs),...
                        this.enzymeComplexCompartmentIndexs, ...
                        this.enzymeComplexLocalIndexs);
                    
                    this.copyToStateHelper(this.rna, ...
                        'counts', ...
                        'boundEnzymes',...
                        this.rna.boundIndexs(this.enzymeRNAGlobalIndexs),...
                        this.enzymeRNACompartmentIndexs, ...
                        this.enzymeRNALocalIndexs);
                    this.copyToStateHelper(this.monomer, ...
                        'counts', ...
                        'boundEnzymes',...
                        this.monomer.boundIndexs(this.enzymeMonomerGlobalIndexs),...
                        this.enzymeMonomerCompartmentIndexs, ...
                        this.enzymeMonomerLocalIndexs);
                    this.copyToStateHelper(this.complex, ...
                        'counts', ...
                        'boundEnzymes',...
                        this.complex.boundIndexs(this.enzymeComplexGlobalIndexs),...
                        this.enzymeComplexCompartmentIndexs, ...
                        this.enzymeComplexLocalIndexs);
                end
            end
            
            if ~isempty(this.substrates)
                if numTimePoints == 1
                    if ~isempty(this.substrateStimulusGlobalCompartmentIndexs)
                        this.stimulus.values(this.substrateStimulusGlobalCompartmentIndexs) = ...
                            this.substrates(this.substrateStimulusLocalIndexs, :);
                    end
                    if ~isempty(this.substrateMetaboliteGlobalCompartmentIndexs)
                        this.metabolite.counts(this.substrateMetaboliteGlobalCompartmentIndexs) = ...
                            this.substrates(this.substrateMetaboliteLocalIndexs, :);
                    end
                    if ~isempty(this.substrateRNAGlobalCompartmentIndexs)
                        this.rna.counts(this.substrateRNAGlobalCompartmentIndexs) = ...
                            this.substrates(this.substrateRNALocalIndexs, :);
                    end
                    if ~isempty(this.substrateMonomerGlobalCompartmentIndexs)
                        this.monomer.counts(this.substrateMonomerGlobalCompartmentIndexs) = ...
                            this.substrates(this.substrateMonomerLocalIndexs, :);
                    end
                    if ~isempty(this.substrateComplexGlobalCompartmentIndexs)
                        this.complex.counts(this.substrateComplexGlobalCompartmentIndexs) = ...
                            this.substrates(this.substrateComplexLocalIndexs, :);
                    end
                else
                    this.copyToStateHelper(this.stimulus, ...
                        'values', ...
                        'substrates', ...
                        this.substrateStimulusGlobalIndexs, ...
                        this.substrateStimulusCompartmentIndexs, ...
                        this.substrateStimulusLocalIndexs);
                    this.copyToStateHelper(this.metabolite, ...
                        'counts', ...
                        'substrates',...
                        this.substrateMetaboliteGlobalIndexs,...
                        this.substrateMetaboliteCompartmentIndexs, ...
                        this.substrateMetaboliteLocalIndexs);
                    this.copyToStateHelper(this.rna, ...
                        'counts', ...
                        'substrates',...
                        this.rna.matureIndexs(this.substrateRNAGlobalIndexs),...
                        this.substrateRNACompartmentIndexs, ...
                        this.substrateRNALocalIndexs);
                    this.copyToStateHelper(this.monomer, ...
                        'counts', ...
                        'substrates',...
                        this.monomer.matureIndexs(this.substrateMonomerGlobalIndexs),...
                        this.substrateMonomerCompartmentIndexs, ...
                        this.substrateMonomerLocalIndexs);
                    this.copyToStateHelper(this.complex, ...
                        'counts', ...
                        'substrates',...
                        this.complex.matureIndexs(this.substrateComplexGlobalIndexs),...
                        this.substrateComplexCompartmentIndexs, ...
                        this.substrateComplexLocalIndexs);
                end
            end
        end
    end

    %Allocate memory for state
    %- stimuli
    %- substrates
    %- enzymes
    %- bound enzymes
    %
    %Called by simulation allocateMemoryForState method which is called at the
    %beginning of simulation.initializeConstants
    methods
        function allocateMemoryForState(this, numTimePoints)
            stimuliCompartments = max([
                size(this.stimulusStimulusGlobalIndexs, 2)
                size(this.stimulusMetaboliteGlobalIndexs, 2)
                size(this.stimulusRNAGlobalIndexs, 2)
                size(this.stimulusMonomerGlobalIndexs, 2)
                size(this.stimulusComplexGlobalIndexs, 2)]);
            substrateCompartments = max([
                size(this.substrateStimulusGlobalIndexs, 2)
                size(this.substrateMetaboliteGlobalIndexs, 2)
                size(this.substrateRNAGlobalIndexs, 2)
                size(this.substrateMonomerGlobalIndexs, 2)
                size(this.substrateComplexGlobalIndexs, 2)]);
            enzymeCompartments = max([
                size(this.enzymeStimulusGlobalIndexs, 2)
                size(this.enzymeMetaboliteGlobalIndexs, 2)
                size(this.enzymeRNAGlobalIndexs, 2)
                size(this.enzymeMonomerGlobalIndexs, 2)
                size(this.enzymeComplexGlobalIndexs, 2)]);
            
            this.stimuli      = zeros(length(this.stimuliWholeCellModelIDs),   stimuliCompartments, numTimePoints);
            this.substrates   = zeros(length(this.substrateWholeCellModelIDs), substrateCompartments, numTimePoints);
            this.enzymes      = zeros(length(this.enzymeWholeCellModelIDs), enzymeCompartments, numTimePoints);
            this.boundEnzymes = zeros(length(this.enzymeWholeCellModelIDs), enzymeCompartments, numTimePoints);
        end
    end

    %model. Subclasses must implement these methods
    methods (Abstract = true)
        %Calculate contribution to FBA objective using
        %- weight fractions: dnaWeight, rnaWeight, proteinWeight, remainingWeight
        %- macromolecular steady-state: RNAs, monomers, complexs
        %- macromolecular decays: rnaDecays, monomerDecays, complexDecays
        %- macromolecular productions (steady-state + decays): rnaProductions, monomerProducts, complexProductions
        %
        %Calculate minimum expression (at beginning of cell cycle) consistent with
        %cell cycle length given current estimates of
        % - Cell weight
        % - Cell cycle length, area under cell growth curve
        % - Composition: dNMP, NMP, AA, other
        % - Expression: RNA, gene, protein monomers
        % - Decay rates: RNA, protein monomers
        %
        %Called by simulation fitConstants method prior to
        %after simulation.initializeConstants and after
        %simulation.initializeState
        %
        %This function can access the cell state and local process state, but
        %shouldn't modify this state or any constants
        [biomassProduction, byproducts, minimumEnzymeExpression, maximumEnzymeExpression] = ...
            calcResourceRequirements_LifeCycle(this, constants, states)

        %initialization
        %
        %called once, before the simulation begins
        initializeState(this)

        %resource requirements
        %
        %called once per time step, prior to evolveState
        result = calcResourceRequirements_Current(this)

        %simulation
        %
        %called once per time step during the simulation
        evolveState(this)
    end
    
    %get/set methods of annotation properties
    methods
        %options
        function value = get.optionNames(this)
            value = [this.optionNames__; {'stepSizeSec'; 'verbosity'; 'seed'}];
        end

        %fixed constants
        function value = get.fixedConstantNames(this)
            value = this.computeFixedConstantsNames();
        end
        
        %provides way for subclasses to override getter for
        %fixedConstantNames property. Annotates as fixed constants:
        %- compartments
        %- molecular weights
        function value = computeFixedConstantsNames(this)
            value = [{
                'stimuliCompartments';
                'substrateCompartments';
                'enzymeCompartments';
                'substrateMolecularWeights';
                'enzymeMolecularWeights'};
                this.fixedConstantNames__];
        end
        
        %fit constants
        function value = get.fittedConstantNames(this)
            value = this.computeFittedConstantsNames();
        end

        %provides way for subclasses to override getter for
        %fittedConstantNames property. Annotates as fitted constants:
        %- biomassProduction
        %- byproducts
        %- minimumEnzymeExpression
        function value = computeFittedConstantsNames(this)
            value = this.fittedConstantNames__;
        end

        %dependent time courses
        function value = get.localStateNames(this)
            value = this.computeLocalStateNames();
        end

        %provides way for subclasses to override getter for
        %localStateNames property. Annotates as
        %dependent time courses.
        %- stimuli
        %- substrates
        %- enzymes
        %- boundEnzymes
        function value = computeLocalStateNames(this)
            value = [{
                'stimuli'; 'substrates'; 'enzymes'; 'boundEnzymes'}; 
                this.localStateNames__];
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

    %Compute indices of the global compartments of stimuli, substrates, and
    %enzymes. Computed by getter methods from the compartment indices in
    %the above section. These properties have size [numComponents of particular type within process (eg stimuli,
    %substrates, enzymes) X compartments]
    methods
        %Compute indices of the global compartments of stimuli.
        function value = get.stimuliCompartments(this)
            value = this.getComponentCompartments('stimulus', 'stimuli');
        end

        %Compute indices of the global compartments of substrates.
        function value = get.substrateCompartments(this)
            value = this.getComponentCompartments('substrate', 'substrate');
        end

        %Compute indices of the global compartments of enzymes.
        function value = get.enzymeCompartments(this)
            value = this.getComponentCompartments('enzyme', 'enzyme');
        end

        %Compute indices of the global compartments of components.
        function value = getComponentCompartments(this, type, types)
            numCompartments = max([
                size(this.([type 'StimulusGlobalIndexs']),2)
                size(this.([type 'MetaboliteGlobalIndexs']),2)
                size(this.([type 'RNAGlobalIndexs']),2)
                size(this.([type 'MonomerGlobalIndexs']),2)
                size(this.([type 'ComplexGlobalIndexs']),2)]);

            value = zeros(size(this.([types 'WholeCellModelIDs']),1),  numCompartments);

            value(this.([type 'StimulusLocalIndexs']),:) = ...
                this.([type 'StimulusCompartmentIndexs']);

            value(this.([type 'MetaboliteLocalIndexs']),:) = ...
                this.([type 'MetaboliteCompartmentIndexs']);

            value(this.([type 'RNALocalIndexs']),:) = ...
                this.([type 'RNACompartmentIndexs']);

            value(this.([type 'MonomerLocalIndexs']),:) = ...
                this.([type 'MonomerCompartmentIndexs']);

            value(this.([type 'ComplexLocalIndexs']),:) = ...
                this.([type 'ComplexCompartmentIndexs']);
        end
    end

    %get methods of dependent local state
    methods
        %dry weight represented by state variables accessible to process
        %including stimuli, substrates, and enzymes, subclasses may need to
        %override the getDryWeight method which calculates this
        function value = get.dryWeight(this)
            value = this.getDryWeight;
            for i = 1:numel(this.states)
                value = value + this.states{i}.dryWeight;
            end
        end

        %Provides way for subclasses to override getter method of
        %dryWeight property
        function value = getDryWeight(this)
            if size(this.substrates, 3) == 1
                value = ...
                    this.substrateMolecularWeights' * sum(this.substrates,   2) + ...
                    this.enzymeMolecularWeights'    * sum(this.enzymes,      2) + ...
                    this.enzymeMolecularWeights'    * sum(this.boundEnzymes, 2);
            else
                value = permute(...
                    this.substrateMolecularWeights' * permute(sum(this.substrates,   2),[1 3 2]) + ...
                    this.enzymeMolecularWeights'    * permute(sum(this.enzymes,      2),[1 3 2]) + ...
                    this.enzymeMolecularWeights'    * permute(sum(this.boundEnzymes, 2),[1 3 2]), ...
                    [1 3 2]);
            end
            value = value / edu.stanford.covert.util.ConstantUtil.nAvogadro;
        end
    end

    %helper methods of initializeConstants
    methods (Access = protected)
        %- computes mapping between process's locally stored state (stimuli,
        %  substrates, enzymes) and the simulation's state using the
        %  stimuliWholeCellModelIDs, substrateWholeCellModelIDs, and
        %  enzymeWholeCellModelIDs properties that subclasses must define
        %- fetches molecular weights of substrates and enzymes from
        %  simulation
        %
        %By default sets the mapping between the process's locally stored
        %state and the simulation state to only copy over 1 compartment for
        %each component. Subclasses can set the retainStimuliCompartments,
        %retainSubstrateCompartments, retainEnzymeCompartments options to
        %true to setup the mapping to copy over all compartments of a
        %component to the simulation's locally stored state.
        function [...
            names, ...
            stimulusLocalIndexs, stimulusGlobalIndexs, stimulusCompartmentIndexs, stimulusGlobalCompartmentIndexs, ...
            metaboliteLocalIndexs, metaboliteGlobalIndexs, metaboliteCompartmentIndexs, metaboliteGlobalCompartmentIndexs, ...
            rnaLocalIndexs, rnaGlobalIndexs, rnaCompartmentIndexs, rnaGlobalCompartmentIndexs,...
            monomerLocalIndexs, monomerGlobalIndexs, monomerCompartmentIndexs, monomerGlobalCompartmentIndexs,...
            complexLocalIndexs, complexGlobalIndexs, complexCompartmentIndexs, complexGlobalCompartmentIndexs,...            
            bound_rnaGlobalCompartmentIndexs, bound_monomerGlobalCompartmentIndexs, bound_complexGlobalCompartmentIndexs, ...
            globalIndexs, molecularWeights] = ...
                initializeConstantsHelper(...
                    this, wholeCellModelIDs, retainCompartments)

            if ~exist('retainCompartments','var')
                retainCompartments = false;
            end

            names             = cell(length(wholeCellModelIDs),1);
            globalIndexs      = zeros(length(wholeCellModelIDs),1);
            compartmentIndexs = zeros(length(wholeCellModelIDs),1);
            stimulis          = false(length(wholeCellModelIDs),1);
            metabolites       = false(length(wholeCellModelIDs),1);
            rnas              = false(length(wholeCellModelIDs),1);
            monomers          = false(length(wholeCellModelIDs),1);
            complexs          = false(length(wholeCellModelIDs),1);

            molecularWeights  = zeros(length(wholeCellModelIDs),1);

            stimuliIDs        = this.stimulus.wholeCellModelIDs;
            metaboliteIDs     = this.metabolite.wholeCellModelIDs;
            rnaIDs            = this.rna.wholeCellModelIDs(this.rna.matureIndexs);
            monomerIDs        = this.monomer.wholeCellModelIDs(this.monomer.matureIndexs);
            complexIDs        = this.complex.wholeCellModelIDs(this.complex.matureIndexs);

            stimNames         = this.stimulus.names;
            metaboliteNames   = this.metabolite.names;
            rnaNames          = this.rna.names(this.rna.matureIndexs);
            monomerNames      = this.monomer.names(this.monomer.matureIndexs);
            complexNames      = this.complex.names(this.complex.matureIndexs);

            stimuliCompartments         = repmat(this.compartment.extracellularIndexs, size(stimuliIDs));
            metaboliteCompartments      = repmat(this.compartment.cytosolIndexs, size(metaboliteIDs));
            metaboliteCompartments(this.metabolite.hydrophobicIndexs) = this.compartment.membraneIndexs;
            rnaCompartments             = repmat(this.compartment.cytosolIndexs, size(rnaIDs));
            monomerCompartments         = this.monomer.compartments(this.monomer.matureIndexs);
            complexCompartments         = this.complex.compartments(this.complex.matureIndexs);

            for i = 1:length(wholeCellModelIDs)
                index = find(strcmp(stimuliIDs,wholeCellModelIDs{i}),1,'first');
                if ~isempty(index)
                    names{i} = stimNames{index};
                    globalIndexs(i) = index;
                    compartmentIndexs(i) = stimuliCompartments(index);
                    molecularWeights(i) = 0;
                    stimulis(i) = true;
                    continue;
                end

                index = find(strcmp(metaboliteIDs,wholeCellModelIDs{i}),1,'first');
                if ~isempty(index)
                    names{i} = metaboliteNames{index};
                    globalIndexs(i) = index;
                    compartmentIndexs(i) = metaboliteCompartments(index);
                    molecularWeights(i) = this.metabolite.molecularWeights(index);
                    metabolites(i) = true;
                    continue;
                end

                index = find(strcmp(rnaIDs,wholeCellModelIDs{i}),1,'first');
                if ~isempty(index)
                    names{i} = rnaNames{index};
                    globalIndexs(i) = index;
                    compartmentIndexs(i) = rnaCompartments(index);
                    molecularWeights(i) = this.rna.molecularWeights(this.rna.matureIndexs(index));
                    rnas(i) = true;
                    continue;
                end

                index = find(strcmp(monomerIDs,wholeCellModelIDs{i}),1,'first');
                if ~isempty(index)
                    names{i} = monomerNames{index};
                    globalIndexs(i) = index;
                    compartmentIndexs(i) = monomerCompartments(index);
                    molecularWeights(i) = this.monomer.molecularWeights(this.monomer.matureIndexs(index));
                    monomers(i) = true;
                    continue;
                end

                index = find(strcmp(complexIDs,wholeCellModelIDs{i}),1,'first');
                if ~isempty(index)
                    names{i} = complexNames{index};
                    globalIndexs(i) = index;
                    compartmentIndexs(i) = complexCompartments(index);
                    molecularWeights(i) = this.complex.molecularWeights(this.complex.matureIndexs(index));
                    complexs(i) = true;
                    continue;
                end
            end
            
            invalidIdxs = find(~(...
                stimulis    | ...
                metabolites | ...
                rnas        | ...
                monomers    | ...
                complexs));
            if ~isempty(invalidIdxs)
                throw(MException('Process:invalidIDs','%s not found',strjoin(', ', wholeCellModelIDs{invalidIdxs})));
            end

            stimulusLocalIndexs   = reshape(find(stimulis),[],1);
            metaboliteLocalIndexs = reshape(find(metabolites),[],1);
            rnaLocalIndexs        = reshape(find(rnas),[],1);
            monomerLocalIndexs    = reshape(find(monomers),[],1);
            complexLocalIndexs    = reshape(find(complexs),[],1);

            if retainCompartments
                nCmp = this.compartment.count;

                stimulusGlobalIndexs        = repmat(globalIndexs(stimulusLocalIndexs),   1, nCmp);
                metaboliteGlobalIndexs      = repmat(globalIndexs(metaboliteLocalIndexs), 1, nCmp);
                rnaGlobalIndexs             = repmat(globalIndexs(rnaLocalIndexs),        1, nCmp);
                monomerGlobalIndexs         = repmat(globalIndexs(monomerLocalIndexs),    1, nCmp);
                complexGlobalIndexs         = repmat(globalIndexs(complexLocalIndexs),    1, nCmp);

                stimulusCompartmentIndexs   = repmat(1:nCmp, length(stimulusLocalIndexs),   1);
                metaboliteCompartmentIndexs = repmat(1:nCmp, length(metaboliteLocalIndexs), 1);
                rnaCompartmentIndexs        = repmat(1:nCmp, length(rnaLocalIndexs),        1);
                monomerCompartmentIndexs    = repmat(1:nCmp, length(monomerLocalIndexs),    1);
                complexCompartmentIndexs    = repmat(1:nCmp, length(complexLocalIndexs),    1);
            else
                stimulusGlobalIndexs        = globalIndexs(stimulusLocalIndexs);
                metaboliteGlobalIndexs      = globalIndexs(metaboliteLocalIndexs);
                rnaGlobalIndexs             = globalIndexs(rnaLocalIndexs);
                monomerGlobalIndexs         = globalIndexs(monomerLocalIndexs);
                complexGlobalIndexs         = globalIndexs(complexLocalIndexs);

                stimulusCompartmentIndexs   = compartmentIndexs(stimulusLocalIndexs);
                metaboliteCompartmentIndexs = compartmentIndexs(metaboliteLocalIndexs);
                rnaCompartmentIndexs        = compartmentIndexs(rnaLocalIndexs);
                monomerCompartmentIndexs    = compartmentIndexs(monomerLocalIndexs);
                complexCompartmentIndexs    = compartmentIndexs(complexLocalIndexs);
            end
            
            stimulusGlobalCompartmentIndexs = sub2ind(...
                [numel(this.stimulus.wholeCellModelIDs) this.compartment.count], ...
                stimulusGlobalIndexs, ...
                stimulusCompartmentIndexs);
            metaboliteGlobalCompartmentIndexs = sub2ind(...
                [numel(this.metabolite.wholeCellModelIDs) this.compartment.count], ...
                metaboliteGlobalIndexs, ...
                metaboliteCompartmentIndexs);
            rnaGlobalCompartmentIndexs = sub2ind(...
                [numel(this.rna.wholeCellModelIDs) this.compartment.count], ...
                this.rna.matureIndexs(rnaGlobalIndexs), ...
                rnaCompartmentIndexs);
            monomerGlobalCompartmentIndexs = sub2ind(...
                [numel(this.monomer.wholeCellModelIDs) this.compartment.count], ...
                this.monomer.matureIndexs(monomerGlobalIndexs), ...
                monomerCompartmentIndexs);
            complexGlobalCompartmentIndexs = sub2ind(...
                [numel(this.complex.wholeCellModelIDs) this.compartment.count], ...
                this.complex.matureIndexs(complexGlobalIndexs), ...
                complexCompartmentIndexs);
            
            bound_rnaGlobalCompartmentIndexs = sub2ind(...
                [numel(this.rna.wholeCellModelIDs) this.compartment.count], ...
                this.rna.boundIndexs(rnaGlobalIndexs), ...
                rnaCompartmentIndexs);
            bound_monomerGlobalCompartmentIndexs = sub2ind(...
                [numel(this.monomer.wholeCellModelIDs) this.compartment.count], ...
                this.monomer.boundIndexs(monomerGlobalIndexs), ...
                monomerCompartmentIndexs);
            bound_complexGlobalCompartmentIndexs = sub2ind(...
                [numel(this.complex.wholeCellModelIDs) this.compartment.count], ...
                this.complex.boundIndexs(complexGlobalIndexs), ...
                complexCompartmentIndexs);
        end
        
        function initializeConstants_overrideSubstrates(this, wholeCellModelIDs)
            %set whole cell model IDs
            this.substrateWholeCellModelIDs = wholeCellModelIDs;

            %recompute mapping of substrates onto simulation
            [this.substrateNames...
             this.substrateStimulusLocalIndexs...
             this.substrateStimulusGlobalIndexs...
             this.substrateStimulusCompartmentIndexs...
             this.substrateStimulusGlobalCompartmentIndexs...
             this.substrateMetaboliteLocalIndexs...
             this.substrateMetaboliteGlobalIndexs...
             this.substrateMetaboliteCompartmentIndexs...
             this.substrateMetaboliteGlobalCompartmentIndexs...
             this.substrateRNALocalIndexs...
             this.substrateRNAGlobalIndexs...
             this.substrateRNACompartmentIndexs...
             this.substrateRNAGlobalCompartmentIndexs...
             this.substrateMonomerLocalIndexs...
             this.substrateMonomerGlobalIndexs...
             this.substrateMonomerCompartmentIndexs...
             this.substrateMonomerGlobalCompartmentIndexs...
             this.substrateComplexLocalIndexs...
             this.substrateComplexGlobalIndexs...
             this.substrateComplexCompartmentIndexs...
             this.substrateComplexGlobalCompartmentIndexs...
             ~, ~, ~, ...
             this.substrateGlobalIndexs...
             this.substrateMolecularWeights] = ...
                this.initializeConstantsHelper(wholeCellModelIDs, size(this.substrates,2)>1);
        end

        function initializeConstants_overrideEnzymes(this, wholeCellModelIDs)
            %set whole cell model IDs
            this.enzymeWholeCellModelIDs  = wholeCellModelIDs;

            %recompute mapping of enzymes onto simulation
            [this.enzymeNames...
             this.enzymeStimulusLocalIndexs...
             this.enzymeStimulusGlobalIndexs...
             this.enzymeStimulusCompartmentIndexs...
             this.enzymeStimulusGlobalCompartmentIndexs...
             this.enzymeMetaboliteLocalIndexs...
             this.enzymeMetaboliteGlobalIndexs...
             this.enzymeMetaboliteCompartmentIndexs...
             this.enzymeMetaboliteGlobalCompartmentIndexs...
             this.enzymeRNALocalIndexs...
             this.enzymeRNAGlobalIndexs...
             this.enzymeRNACompartmentIndexs...
             this.enzymeRNAGlobalCompartmentIndexs...
             this.enzymeMonomerLocalIndexs...
             this.enzymeMonomerGlobalIndexs...
             this.enzymeMonomerCompartmentIndexs...
             this.enzymeMonomerGlobalCompartmentIndexs...
             this.enzymeComplexLocalIndexs...
             this.enzymeComplexGlobalIndexs...
             this.enzymeComplexCompartmentIndexs...
             this.enzymeComplexGlobalCompartmentIndexs...
             this.enzymeBoundRNAGlobalCompartmentIndexs...
             this.enzymeBoundMonomerGlobalCompartmentIndexs...
             this.enzymeBoundComplexGlobalCompartmentIndexs...
             this.enzymeGlobalIndexs...
             this.enzymeMolecularWeights] = ...
                this.initializeConstantsHelper(wholeCellModelIDs, size(this.enzymes,2)>1);
        end
    end
    
    methods
        %Computes local indices of stimuli. That is, this method returns
        %the indices of wholeCellModelIDs within this.stimulusWholeCellModelIDs.
        function value = stimulusIndexs(this, wholeCellModelIDs)
            value = this.componentIndexs(wholeCellModelIDs, 'stimulus');
        end

        %Computes local indices of substrates. That is, this method returns
        %the indices of wholeCellModelIDs within this.substrateWholeCellModelIDs.
        function value = substrateIndexs(this, wholeCellModelIDs)
            value = this.componentIndexs(wholeCellModelIDs, 'substrate');
        end

        %Computes local indices of enzymes. That is, this method returns
        %the indices of wholeCellModelIDs within this.enzymeWholeCellModelIDs.
        function value = enzymeIndexs(this, wholeCellModelIDs)
            value = this.componentIndexs(wholeCellModelIDs, 'enzyme');
        end

        %Computes local indices of components. That is, this method returns
        %the indices of wholeCellModelIDs within this.*WholeCellModelIDs.
        function value = componentIndexs(this, wholeCellModelIDs, type)
            [tf,value] = ismember(wholeCellModelIDs, this.([type 'WholeCellModelIDs']));
            value = value(tf);
            if ~all(tf)
                missingComponents = setdiff(wholeCellModelIDs,this.([type 'WholeCellModelIDs']));
                throw(MException('Process:error', '%s not found', strjoin(', ', missingComponents{:})));
            end
        end
    end

    %helper functions of get/set state
    methods
        %Use simulation-process mapping constructed by initializeConstants
        %to copy component values from simulation. Called by
        %copyFromState.
        function lclState = copyFromStateHelper(~, gblState, gblIdxs, cmpIdxs)
            numTimePoints = size(gblState, 3);
            lclState = zeros([size(gblIdxs) numTimePoints]);
            if isempty(gblIdxs)
                return;
            end
            
            times = permute(1:numTimePoints, [1 3 2]);
            if ~isa(gblState, 'edu.stanford.covert.util.SparseMat');
                lclState = gblState(sub2ind(...
                    size(gblState), ...
                    gblIdxs(:, :, ones(numTimePoints, 1)),...
                    cmpIdxs(:, :, ones(numTimePoints, 1)),...
                    times(ones(size(gblIdxs, 1), 1), ones(size(gblIdxs, 2), 1), :)));
            else
                lclState = reshape(gblState([...
                    reshape(gblIdxs(:, :, ones(numTimePoints, 1)), [], 1) ...
                    reshape(cmpIdxs(:, :, ones(numTimePoints, 1)), [], 1) ...
                    reshape(times(ones(size(gblIdxs, 1), 1), ones(size(gblIdxs, 2), 1), :), [], 1)
                    ]), [size(gblIdxs) numTimePoints]);
            end
        end
        
        %Use simulation-process mapping constructed by initializeConstants
        %to copy component values to simulation. Called by
        %copyToState.
        function copyToStateHelper(this, state, statePropertyName, lclPropertyName, stateIdxs, cmpIdxs, lclIdxs)
            if isempty(lclIdxs)
                return;
            end
            
            numTimePoints = size(state.(statePropertyName), 3);
            state.(statePropertyName)(sub2ind(...
                size(state.(statePropertyName)), ...
                repmat(stateIdxs, [1 1 numTimePoints]),...
                repmat(cmpIdxs, [1 1 numTimePoints]),...
                permute(repmat((1:numTimePoints)', [1 size(stateIdxs)]), [2 3 1]))) = ...
                this.(lclPropertyName)(lclIdxs, :, :);
        end
        
        %Converts units of process state for example from numbers of
        %components to their concentration. type and types inputs are
        %strings representing the process property that needs to be scaled.
        %Scale is a string representing what units the property should be
        %scaled to.
        function scaledValues = scaleState(this, type, types, scale)
            %unscaled values
            unscaledValues = this.(types);
            
            %scale values
            switch scale
                case 'mM'
                    scaledValues = unscaledValues / ...
                        (edu.stanford.covert.util.ConstantUtil.nAvogadro/1000) / ...
                        this.geometry.volume;
                case 'M'
                    scaledValues = unscaledValues / ...
                        edu.stanford.covert.util.ConstantUtil.nAvogadro / ...
                        this.geometry.volume;
                otherwise
                    scaledValues = unscaledValues;
            end
            
            %don't scale stimuli
            scaledValues(this.([type 'StimulusLocalIndexs']),:) = ...
                unscaledValues(this.([type 'StimulusLocalIndexs']),:);
        end
    end

    %Methods which return gene composition of components as adjacency
    %matrix (genes X components). Used by applyPerturbations method of
    %simulation and testGeneEssentiality method of process test classes.
    methods
        %Returns gene composition of stimuli as adjacency matrix (genes X stimuli)
        function value = stimuliGeneComposition(this)
            value = this.componentGeneComposition('stimulus', 'stimuli');
        end

        %Returns gene composition of substrate as adjacency matrix (genes X substrate)
        function value = substrateGeneComposition(this)
            value = this.componentGeneComposition('substrate', 'substrate');
        end

        %Returns gene composition of enzyme as adjacency matrix (genes X enzyme)
        function value = enzymeGeneComposition(this)
            value = this.componentGeneComposition('enzyme', 'enzyme');
        end

        %methods which return gene composition of components as adjacency
        %matrix (genes X components). Called by
        %- stimuliGeneComposition
        %- substrateGeneComposition
        %- enzymeGeneComposition
        function value = componentGeneComposition(this, type, types)
            value = zeros(length(this.gene.wholeCellModelIDs), length(this.([types 'WholeCellModelIDs'])));

            value(:, this.([type 'RNALocalIndexs'])) = ...
                this.rna.matureRNAGeneComposition(:, this.([type 'RNAGlobalIndexs'])(:,1));

            value(sub2ind(...
                size(value),...
                this.gene.mRNAIndexs(this.([type 'MonomerGlobalIndexs'])(:,1)),...
                this.([type 'MonomerLocalIndexs']))) = ...
                1;

            value(:, this.([type 'ComplexLocalIndexs'])) = ...
                sum(this.complex.proteinComplexComposition(:, this.([type 'ComplexGlobalIndexs'])(:,1),:),3);
        end
    end

    %methods for selecting and converting indices
    methods
        function i = substrateToMetabolite(this, i)
            [tf,idxs] = ismember(i, this.substrateMetaboliteLocalIndexs);
            i = this.substrateMetaboliteGlobalIndexs(idxs(tf));
        end
        
        function i = metaboliteToSubstrate(this, i)
            [tf,idxs] = ismember(i, this.substrateMetaboliteGlobalIndexs);
            i = this.substrateMetaboliteLocalIndexs(idxs(tf));
        end
        
        function i = enzymeToMonomer(this, i)
            [tf,idxs] = ismember(i, this.enzymeMonomerLocalIndexs);
            i = this.enzymeMonomerGlobalIndexs(idxs(tf));
        end
        
        function i = monomerToEnzyme(this, i)
            [tf,idxs] = ismember(i, this.enzymeMonomerGlobalIndexs);
            i = this.enzymeMonomerLocalIndexs(idxs(tf));
        end
        
        function i = enzymeToComplex(this, i)
            [tf,idxs] = ismember(i, this.enzymeComplexLocalIndexs);
            i = this.enzymeComplexGlobalIndexs(idxs(tf));
        end
        
        function i = complexToEnzyme(this, i)
            [tf,idxs] = ismember(i, this.enzymeComplexGlobalIndexs);
            i = this.enzymeComplexLocalIndexs(idxs(tf));
        end
    end
end
