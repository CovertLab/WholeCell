%Ribosome
%
% @wholeCellModelID State_Ribosome
% @name             Ribosome
% @description
%
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/4/2011
classdef Ribosome < edu.stanford.covert.cell.sim.CellState
    %property annotations
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {}; %names of process properties that are considered fixed constants
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        stateNames              = {   %names of properties which are part of the simulation's state
            'states';
            'boundMRNAs';
            'mRNAPositions';
            'tmRNAPositions';
            };
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            'stateOccupancies'
            'nActive'
            'nNotExist'
            'nStalled'
            %'nMRNAsBound'
            };
    end
    
    %enumerations
    properties (Constant)
        activeIndex   =  1; %index in stateOccupancies
        notExistIndex =  2; %index in stateOccupancies
        stalledIndex  =  3; %index in stateOccupancies
        
        activeValue   =  1; %value within states
        notExistValue =  0; %value within states
        stalledValue  = -1; %value within states
    end
        
    %state
    properties
        states           %ribosome state
        boundMRNAs       %index of gene each ribosome is bound to
        mRNAPositions    %length of nascent monomer
        tmRNAPositions   %length of nascent proteolysis tag
    end
    
    %dependent local state (implemented as dependent property for convenience)
    properties (Dependent = true)
        stateOccupancies %numbers of ribosomes in various states
        
        nActive          %number of actively translating
        nNotExist        %number of non-existent ribosomes
        nStalled         %number of stalled ribosomes
        
        nMRNAsBound      %Number of ribosomes actively translating each mRNA
    end
    
    %dependent state
    properties (Constant)
        dryWeight = 0;                 %dry weight of this class' state properties
    end
        
    %references to other parts of cell state
    properties
        gene        %gene
        polypeptide %polypeptides
        chromosome  %chromosomes
        rna         %RNA
        complex     %macromolecular complexes
    end
    
    %constructor
    methods
        function this = Ribosome(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    methods
        %set references to state objects
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.CellState(simulation);
            
            this.gene = simulation.gene;
            this.polypeptide = simulation.state('Polypeptide');
            this.chromosome = simulation.state('Chromosome');
            this.rna = simulation.state('Rna');
            this.complex = simulation.state('ProteinComplex');
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemory(this, numTimePoints)
            this.states         = zeros(0, 1, numTimePoints);
            this.boundMRNAs     = zeros(0, 1, numTimePoints);
            this.mRNAPositions  = zeros(0, 1, numTimePoints);
            this.tmRNAPositions = zeros(0, 1, numTimePoints);
        end
    end
    
    %initialize state
    methods
        function initialize(~)
        end
    end
    
    methods
        function releaseRibosome(this, nRibs, nTMRNA, mRNAs)
            poly = this.polypeptide;
            
            ribIdxs = [];
            if nRibs > 0
                ribIdxs = this.randStream.randomlySelectNRows(find(this.states ~= this.notExistValue), nRibs); %#ok<FNDSB>
                if numel(ribIdxs) < nRibs
                    throw(MException('Ribosome:error', 'Unable release ribosomes'));
                end
            end
            
            tmRNAIdxs = [];
            if nTMRNA > 0
                tmRNAIdxs = this.randStream.randomlySelectNRows(find(this.states == this.stalledValue), nTMRNA); %#ok<FNDSB>
                if numel(tmRNAIdxs) < nTMRNA
                    throw(MException('Ribosome:error', 'Unable release tmRNA'));
                end
            end
            
            mRNAIdxs = [];
            if any(mRNAs)
                mRNAGenes = this.rna.matureRNAGeneComposition(this.gene.mRNAIndexs, this.rna.matureMRNAIndexs) * mRNAs;
                totMRNAGenes = this.rna.matureRNAGeneComposition(this.gene.mRNAIndexs, this.rna.matureMRNAIndexs) * ...
                    this.rna.counts(this.rna.matureIndexs(this.rna.matureMRNAIndexs), this.rna.compartment.cytosolIndexs);
                tmp = find(mRNAGenes);
                mRNAIdxs = [];
                for i = 1:numel(tmp)
                    tmpIdxs = find(this.boundMRNAs == tmp(i));
                    if isempty(tmpIdxs)
                        continue;
                    end
                    mRNAIdxs = [mRNAIdxs;
                        tmpIdxs(this.randStream.rand(size(tmpIdxs)) <= mRNAGenes(tmp(i)) / totMRNAGenes(tmp(i)));
                        ]; %#ok<AGROW>
                end
            end
            
            if ~isempty(ribIdxs) + ~isempty(mRNAIdxs) + ~isempty(tmRNAIdxs) > 1
                idxs = unique([ribIdxs; mRNAIdxs; tmRNAIdxs]);
            else
                idxs = [ribIdxs; mRNAIdxs; tmRNAIdxs];
            end
            
            %abort peptides
            for i = 1:numel(idxs)
                if this.mRNAPositions(idxs(i)) == 0 && this.tmRNAPositions(idxs(i)) == 0
                    continue;
                end
                poly.abortedPolypeptides = [
                    poly.abortedPolypeptides
                    this.boundMRNAs(idxs(i)) this.mRNAPositions(idxs(i)) 0
                    ];
                
                if this.states(idxs(i)) == this.stalledValue
                    poly.abortedPolypeptides(end, 3) = this.tmRNAPositions(idxs(i));
                end
            end
            
            %release ribosomes
            if numel(idxs) - nRibs > 0
                this.complex.counts(this.complex.boundIndexs(this.complex.ribosome70SIndexs), this.complex.compartment.cytosolIndexs) = ...
                    + this.complex.counts(this.complex.boundIndexs(this.complex.ribosome70SIndexs), this.complex.compartment.cytosolIndexs) ...
                    - (numel(idxs) - nRibs);
                this.complex.counts(this.complex.matureIndexs([this.complex.ribosome30SIndexs; this.complex.ribosome50SIndexs]), this.complex.compartment.cytosolIndexs) = ...
                    + this.complex.counts(this.complex.matureIndexs([this.complex.ribosome30SIndexs; this.complex.ribosome50SIndexs]), this.complex.compartment.cytosolIndexs) ...
                    + (numel(idxs) - nRibs);
            end
            
            %release tmRNA
            if sum(this.states(idxs) == this.stalledValue) - nTMRNA > 0
                this.rna.counts(this.rna.boundIndexs(this.rna.matureTMRNAIndexs), this.rna.compartment.cytosolIndexs) = ...
                    + this.rna.counts(this.rna.boundIndexs(this.rna.matureTMRNAIndexs), this.rna.compartment.cytosolIndexs) ...
                    - (sum(this.states(idxs) == this.stalledValue) - nTMRNA);
                this.rna.counts(this.rna.matureIndexs(this.rna.matureTMRNAIndexs), this.rna.compartment.cytosolIndexs) = ...
                    + this.rna.counts(this.rna.matureIndexs(this.rna.matureTMRNAIndexs), this.rna.compartment.cytosolIndexs) + ...
                    + (sum(this.states(idxs) == this.stalledValue) - nTMRNA);
            end
            
            %update polypeptide, ribosome states
            this.states(idxs) = this.notExistValue;
            this.boundMRNAs(idxs) = this.notExistValue;
            this.mRNAPositions(idxs) = this.notExistValue;
            this.tmRNAPositions(idxs) = this.notExistValue;
            poly.boundMRNAs(idxs) = this.boundMRNAs(idxs);
            poly.nascentMonomerLengths(idxs) = this.mRNAPositions(idxs);
            poly.proteolysisTagLengths(idxs) = this.tmRNAPositions(idxs);
        end
    end
    
    %getters
    methods
        %occupancies of ribosome states
        %fraction of all 30S Ribosomes that are active or are in the
        %30SIF3 state.
        function value = get.stateOccupancies(this)
            value = zeros(3, 1, size(this.states, 3));
            value(this.activeIndex,   :, :) = this.nActive;
            value(this.notExistIndex, :, :) = this.nNotExist;
            value(this.stalledIndex,  :, :) = this.nStalled;
            value = value ./ repmat(sum(value, 1), [3 1 1]);
        end

        %number of actively translating ribosomes
        function value = get.nActive(this)
            value = sum(this.states == this.activeValue);
        end

        %number of non-existent ribosomes
        function value = get.nNotExist(this)
            value = sum(this.states == this.notExistValue);
        end

        %number of stalled ribosomes
        function value = get.nStalled(this)
            value = sum(this.states == this.stalledValue);
        end        

        %Number of ribosomes actively translating each mRNA
        function value = get.nMRNAsBound(this)
            value = zeros(length(this.polypeptide.monomerLengths), 1, size(this.states, 3));
            for i = 1:size(this.states, 1)
                for k = 1:size(this.states, 3)
                    boundMRNA = this.boundMRNAs(i, 1, k);
                    if boundMRNA < 1
                        continue; 
                    end
                    value(boundMRNA, 1, k) = ...
                        value(boundMRNA, 1, k) + 1;
                end
            end
        end
    end
end
