%RNA Polymerase
%
% @wholeCellModelID State_RNAPolymerase
% @name             RNA Polymerases
% @description
%
%   states represents the current state / pseudostate (actively
%   transcribing, specifically bound, non-specifically bound, free,
%   non-existent) of each RNA polymerase, where each state is indicated by the
%   enumeration:
%   - rnaPolymeraseActivelyTranscribingValue
%   - rnaPolymeraseSpecificallyBoundValue
%   - rnaPolymeraseNonSpecificallyBoundValue
%   - rnaPolymeraseFreeValue
%   - rnaPolymeraseNotExistValue (state exists as a way to account for memory
%     allocated for future RNA polymerases)
%
% Information about positions of the polymerases on the DNA and the
% progress of RNA polymerases transcribing specific transcrips is all
% contained within the chromosomeState class and newTranscriptState class.
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 12/13/2010

classdef RNAPolymerase < edu.stanford.covert.cell.sim.CellState
    %property annotations
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity'
            'seed'
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'stateExpectations'
            };
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        stateNames              = {   %names of properties which are part of the simulation's state
            'states'
            'positionStrands'
            'transcriptionFactorBindingProbFoldChange'
            'supercoilingBindingProbFoldChange'
            };
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            'stateOccupancies'
            'nActive'
            'nSpecificallyBound'
            'nNonSpecificallyBound'
            'nFree'
            };
    end
    
    %constants
    properties (Constant)
        activelyTranscribingIndex   = 1; %index within rnaPolymeraseStateOccupancies
        specificallyBoundIndex      = 2; %index within rnaPolymeraseStateOccupancies
        nonSpecificallyBoundIndex   = 3; %index within rnaPolymeraseStateOccupancies
        freeIndex                   = 4; %index within rnaPolymeraseStateOccupancies
        
        activelyTranscribingValue   = 1; %value within states
        specificallyBoundValue      = -3; %value within states
        nonSpecificallyBoundValue   = -1; %value within states
        freeValue                   = -2; %value within states
        notExistValue               = 0; %value within states
        
        stateValues                 = [1; -3; -1; -2]; %values of states
    end
    
    %fixed biological constants
    properties
        stateExpectations       %Expected fractional occupancies of RNA polymerase states
    end
    
    %state
    properties
        states                                    %RNA polymerase state (see Transcription process)
        positionStrands                           %the position on the chromosome(s) for each RNApolymerase
        transcriptionFactorBindingProbFoldChange  %fold change effect of currently active transcription factors on RNA polymerase binding probabilities [nTUs X 2]
        supercoilingBindingProbFoldChange         %fold change in transcription probabilities due to supercoiling sigma, [# transcUnits x 2 chromosomes]
    end
    
    properties (Dependent = true)
        stateOccupancies        %number of RNA polymerases in various states
        nActive                 %number of actively transcribing RNA polymerases
        nSpecificallyBound      %number of specifically bound RNA polymerases
        nNonSpecificallyBound   %number of non-specifically bound RNA polymerases
        nFree                   %number of free RNA polymerases
    end
    
    %dependent state
    properties (Constant)
        dryWeight = 0;          %dry weight of this class' state properties
    end
    
    %references to other parts of cell state
    properties
        chromosome
        transcripts
    end
    
    %constructor
    methods
        function this = RNAPolymerase(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    methods
        %build object graph
        function storeObjectReferences(this, simulation)
            this.chromosome = simulation.state('Chromosome');
            this.transcripts = simulation.state('Transcript');
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemory(this, numTimePoints)
            nTUs = numel(this.chromosome.transcriptionUnitWholeCellModelIDs);
            
            this.states                                   = zeros(0, 1, numTimePoints);
            this.positionStrands                          = zeros(0, 2, numTimePoints);
            this.transcriptionFactorBindingProbFoldChange = zeros(nTUs, 2, numTimePoints);
            this.supercoilingBindingProbFoldChange        = zeros(nTUs, 2, numTimePoints);
        end
    end
    
    %initialize state
    methods
        function initialize(~)
        end
    end
    
    %external interface
    methods
        function degradeFreePolymerase(this, nPolymerases)
            idxs = find(this.states == this.freeValue, nPolymerases, 'first');
            this.states(idxs) = this.notExistValue;
        end
        
        function releasePolymerase(this, posStrnds, proteinIsDegraded)
            import edu.stanford.covert.util.CircularSparseMat;
            
            c = this.chromosome;
            t = this.transcripts;
            
            tfs = CircularSparseMat.ismember_subs(this.positionStrands, posStrnds, [c.sequenceLen c.nCompartments]);
            idxs = find(tfs);
            if isempty(idxs)
                return;
            end
            
            tfs(tfs) = this.states(tfs) > this.activelyTranscribingValue;
            idxs2 = find(tfs);
            for i = 1:numel(idxs2)
                if t.boundTranscriptProgress(idxs2(i)) <= 1
                    continue;
                end
                t.abortedTranscripts = [
                    t.abortedTranscripts
                    t.boundTranscriptionUnits(idxs2(i)) t.boundTranscriptProgress(idxs2(i))-1
                    ];
            end
            if proteinIsDegraded
                this.states(idxs) = this.notExistValue;
            else
                this.states(idxs) = this.freeValue;
            end
            this.positionStrands(idxs, :) = 0;
            t.boundTranscriptionUnits(idxs) = t.nullTranscriptValue;
            t.boundTranscriptProgress(idxs) = t.nullTranscriptValue;
            t.boundTranscriptChromosome(idxs) = t.nullTranscriptValue;
        end
    end
    
    %getters
    methods
        %number of actively transcribing RNA polymerases
        function value = get.stateOccupancies(this)
            numTimePoints = size(this.states,3);
            
            value = zeros(4,1,numTimePoints);
            value(this.activelyTranscribingIndex,:,:) = this.nActive;
            value(this.specificallyBoundIndex,:,:)    = this.nSpecificallyBound;
            value(this.nonSpecificallyBoundIndex,:,:) = this.nNonSpecificallyBound;
            value(this.freeIndex,:,:)                 = this.nFree;
            
            value = value ./ repmat(sum(sum(value,1),2), [4,1,1]);
        end
        
        %number of actively transcribing RNA polymerases
        function value = get.nActive(this)
            value = sum(this.states >= this.activelyTranscribingValue);
        end
        
        %number of specifically bound RNA polymerases
        function value = get.nSpecificallyBound(this)
            value = sum(this.states == this.specificallyBoundValue);
        end
        
        %number of non-specifically bound RNA polymerases
        function value = get.nNonSpecificallyBound(this)
            value = sum(this.states == this.nonSpecificallyBoundValue);
        end
        
        %number of free RNA polymerases
        function value = get.nFree(this)
            value = sum(this.states == this.freeValue);
        end
    end
end
