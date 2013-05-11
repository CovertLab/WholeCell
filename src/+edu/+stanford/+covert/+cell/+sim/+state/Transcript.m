%Transcripts
%
% @wholeCellModelID State_Transcript
% @name             Transcripts
% @description
%
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 12/16/2010
classdef Transcript < edu.stanford.covert.cell.sim.CellState
    %property annotations
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity'
            'seed'
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'genomeLength'
            'transcriptionUnitFivePrimeCoordinates'
            'transcriptionUnitDirections'
            'transcriptionUnitSequences'
            'transcriptionUnitLengths'
            };
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such  
        stateNames              = {   %names of properties which are part of the simulation's state
            'boundTranscriptionUnits'
            'boundTranscriptProgress'
            'boundTranscriptChromosome'
            'abortedTranscripts'
            };
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            %totalBaseCounts
            %'abortedSequences'
            %'rnaMaxRNAPolymeraseState'
            'rnaBoundRNAPolymerases'
            };
    end
       
    %constants
    properties (Constant)
        nullTranscriptValue = 0; %value for transcripts that do not exist; needed to sync with rnaPolymerases
        newTranscriptValue  = 1; %value for transcripts that do not exist; needed to sync with rnaPolymerases
    end

    %constants
    properties
        genomeLength                                  %length of genome
        
        transcriptionUnitWholeCellModelIDs            %whole cell model ids of transcription units
        transcriptionUnitFivePrimeCoordinates         %transcription unit 5' coordinates
        transcriptionUnitDirections                   %transcription unit directions
        transcriptionUnitSequences                    %transcription unit sequences
        transcriptionUnitLengths                      %transcription unit sequence lengths
        
        nmpMolecularWeights                           %molecular weights of NMPs
    end
    
    %state
    properties
        %these first three properties describing transcripts currently
        %being made need to be kept in sync with each other. 
        boundTranscriptionUnits     %transcription unit of all transcripts currently being made
        boundTranscriptProgress     %progress on transcription unit of all transcripts currently being made
        boundTranscriptChromosome   %chromosome of all transcripts currently being made
        
        abortedTranscripts          %identity of aborted transcripts (transcripts X [transcription Unit, length])
        totalBaseCounts             %total base counts of aborted transcripts
    end   
    
    %dependent state
    properties (Dependent = true)
        abortedSequences            %sequences of aborted transcripts
        rnaMaxRNAPolymeraseState    %Maximum nascent transcript length of each transcription unit
        rnaBoundRNAPolymerases      %Number of RNA polymerases actively transcribing each gene
    end

    properties (Dependent = true)
        dryWeight                   %dry weight of this class' state properties
    end
    
    %references to other parts of cell state
    properties   
        chromosome        
    end

    %constructor
    methods
        function this = Transcript(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);            
        end
    end
    
    methods
        %build object graph
        function storeObjectReferences(this, simulation)
            this.chromosome = simulation.state('Chromosome');
        end
        
        function initializeConstants(this, knowledgeBase, simulation)
		    this.initializeConstants@edu.stanford.covert.cell.sim.CellState(knowledgeBase, simulation);
			
            this.genomeLength                           = knowledgeBase.genome.sequenceLength;
            
            this.transcriptionUnitWholeCellModelIDs     = {knowledgeBase.transcriptionUnits.wholeCellModelID}';
            this.transcriptionUnitFivePrimeCoordinates  = [knowledgeBase.transcriptionUnits.fivePrimeCoordinate]';
            this.transcriptionUnitDirections            = [knowledgeBase.transcriptionUnits.direction]';
            this.transcriptionUnitSequences             = {knowledgeBase.transcriptionUnits.sequence}';
            this.transcriptionUnitLengths               = [knowledgeBase.transcriptionUnits.sequenceLength]';
            
            this.nmpMolecularWeights = simulation.state('Metabolite').molecularWeights(simulation.state('Metabolite').nmpIndexs);
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemory(this, numTimePoints)   
            this.boundTranscriptionUnits   = zeros(0, 1, numTimePoints);
            this.boundTranscriptProgress   = zeros(0, 1, numTimePoints);
            this.boundTranscriptChromosome = zeros(0, 1, numTimePoints);
            
            this.abortedTranscripts = zeros(0, 2, numTimePoints);
        end
    end
    
    methods
        function initialize(this)
            this.abortedTranscripts = zeros(0, 2);
        end
    end
    
    %getters
    methods
        function value = get.dryWeight(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            nSeq = 0;
            seq = [];
            
            %actively transcribing transcripts
            for i = 1:size(this.boundTranscriptionUnits, 1)
                if this.boundTranscriptProgress(i) > this.nullTranscriptValue
                    nSeq = nSeq + 1;
                    seq = [seq this.transcriptionUnitSequences{this.boundTranscriptionUnits(i)}(1:(this.boundTranscriptProgress(i) - 1))]; %#ok<AGROW>
                end
            end
            
            %aborted
            for i = 1:size(this.abortedTranscripts, 1)
                if this.abortedTranscripts(i, 1) ~= 0
                    nSeq = nSeq + 1;
                    seq = [seq this.transcriptionUnitSequences{this.abortedTranscripts(i, 1)}(1:this.abortedTranscripts(i, 2))]; %#ok<AGROW>
                end
            end
            
            %mass
            value = (...
                + [sum(seq == 'A') sum(seq == 'C') sum(seq == 'G') sum(seq == 'U')] * this.nmpMolecularWeights ...
                - (numel(seq) - nSeq) * (ConstantUtil.elements.H + ConstantUtil.elements.O) ...
                ) / ConstantUtil.nAvogadro;
        end
        
        function value = get.abortedSequences(this)
            value = cell(size(this.abortedTranscripts, 1), 1);
            for i = 1:size(this.abortedTranscripts, 1)
                if this.abortedTranscripts(i, 1) == 0
                    continue;
                end
                value{i} = this.transcriptionUnitSequences{this.abortedTranscripts(i, 1)}(1:this.abortedTranscripts(i, 2));
            end
        end
        
        function value = get.totalBaseCounts(this)
            import edu.stanford.covert.cell.kb.ssRNA;
            
            value = zeros(1, 4);
            
            for i = 1:size(this.boundTranscriptionUnits, 1)
                if this.boundTranscriptProgress(i) > this.nullTranscriptValue
                    seq = this.transcriptionUnitSequences{this.boundTranscriptionUnits(i)}(1:(this.boundTranscriptProgress(i) - 1));
                    value = value + ssRNA.computeBaseCount(seq, 4, 1:4);
                end
            end
            
            abortedSeqs = this.abortedSequences;
            for i = 1:numel(abortedSeqs)
                value = value + ssRNA.computeBaseCount(abortedSeqs{i}, 4, 1:4);
            end
        end
        
        %Maximum nascent transcript length of each transcription unit
        function value = get.rnaMaxRNAPolymeraseState(this)
            value = zeros(length(this.transcriptionUnitLengths), 1, size(this.boundTranscriptionUnits, 3));
            for i = 1:size(this.boundTranscriptionUnits, 1)
                for k = 1:size(this.boundTranscriptionUnits, 3)
                    boundTranscriptionUnit = this.boundTranscriptionUnits(i, 1, k);
                    if boundTranscriptionUnit == 0
                        continue;
                    end
                    value(boundTranscriptionUnit, 1, k) = max(...
                        value(boundTranscriptionUnit, 1, k),...
                        this.transcriptionUnitLengths(i, 1, k));
                end
            end
        end

        %Number of RNA polymerases actively transcribing each gene
        function value = get.rnaBoundRNAPolymerases(this)
            nTUs = numel(this.transcriptionUnitLengths);
            nTimePoints = size(this.boundTranscriptionUnits, 3);
            value = zeros(nTUs, 1, nTimePoints);
            for k = 1:nTimePoints
                value(:, 1, k) = histc(this.boundTranscriptionUnits(:, 1, k), 1:nTUs);
            end
        end
    end
end
