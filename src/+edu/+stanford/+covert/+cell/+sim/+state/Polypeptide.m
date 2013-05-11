%Polypeptide
%
% @wholeCellModelID State_Polypeptide
% @name             Polypeptide
% @description
%
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 12/16/2010
classdef Polypeptide < edu.stanford.covert.cell.sim.CellState
    %property annotations
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'monomerLengths';
            'monomerTRNASequences';
            'monomerAASequences';
            'proteolysisTagLength';
            'proteolysisTagMolecularWeight';
            'proteolysisTagTRNASequence';
            'proteolysisTagAASequence';
            };
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        stateNames              = {   %names of properties which are part of the simulation's state
            'boundMRNAs';
            'nascentMonomerLengths';
            'proteolysisTagLengths';
            'abortedPolypeptides'
            };
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            %'totalBaseCounts'
            %'mRNAMaxNascentMonomerLength'
            %'abortedSequences'
            %'abortedSequenceLengths'
            };
    end
    
    %fixed biological constants
    properties
        aminoAcidMolecularWeights          %amino acid molecular weights
        
        monomerWholeCellModelIDs           %protein monomer molecular weights
        monomerLengths                     %protein monomer lengths
        monomerTRNASequences               %tRNA sequences of mRNAs
        monomerAASequences                 %amino acid sequences of mRNAs

        proteolysisTagLength               %amino acid length of proteolysis tag
        proteolysisTagMolecularWeight      %amino acid weight of proteolysis tag
        proteolysisTagTRNASequence         %tRNA sequence of proteolysis tag
        proteolysisTagAASequence           %amino acid sequence of proteolysis tag
    end
    
    %state
    properties
        boundMRNAs                        %index of gene each ribosome is bound to
        nascentMonomerLengths             %length of nascent monomer
        proteolysisTagLengths             %length of nascent proteolysis tag
        abortedPolypeptides               %identify of aborted polypeptides aborted polypeptide X [monomer, mRNA length, tmRNA length]
    end
    
    %dependent local state (implemented as dependent property for convenience)
    properties (Dependent = true)
        totalBaseCounts                   %counts of amino acids in nascent polypeptides
        mRNAMaxNascentMonomerLength       %Maximum translation progress of each mRNA
        abortedSequences                  %sequences of proteolysis tagged polypeptides
        abortedSequenceLengths            %lengths of proteolysis tagged polypeptides
        dryWeight                         %dry weight of this class' state properties
    end   
    
    %constructor
    methods
        function this = Polypeptide(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    methods
        %build object graph
        function storeObjectReferences(~, ~)
        end
        
        function initializeConstants(this, knowledgeBase, simulation)
            import edu.stanford.covert.cell.kb.ProteinMonomer;

            this.initializeConstants@edu.stanford.covert.cell.sim.CellState(knowledgeBase, simulation);
            
            this.aminoAcidMolecularWeights = [knowledgeBase.metabolites(knowledgeBase.aminoAcidIndexs).molecularWeight]';
            
            %protein monomers
            this.monomerWholeCellModelIDs = {knowledgeBase.proteinMonomers.wholeCellModelID}';
            this.monomerLengths           = [knowledgeBase.proteinMonomers.sequenceLength]';
            this.monomerTRNASequences     = knowledgeBase.proteinMonomerTRNASequences;
            this.monomerAASequences       = {knowledgeBase.proteinMonomers.sequence}';

            %tmRNA proteolysis tag
            tmRNAGene = findobj(knowledgeBase.genes,'wholeCellModelID', 'MG_0004');
            tmRNAProteolysisTag = findobj(tmRNAGene.features,'type','tmRNA proteolysis tag');
            tmRNAProteolysisTagDNASequence = {seqrcomplement(tmRNAProteolysisTag.sequence)};
            tmRNAProteolysisTagTRNASequence = knowledgeBase.computeTRNASequences(tmRNAProteolysisTagDNASequence, false);
            this.proteolysisTagTRNASequence = [0;tmRNAProteolysisTagTRNASequence{1}];
            this.proteolysisTagLength = length(this.proteolysisTagTRNASequence);
            this.proteolysisTagAASequence = [...
                ProteinMonomer.bases(knowledgeBase.aminoAcidIndexs == tmRNAGene.aminoAcid.idx)...
                nt2aa(seqrcomplement(tmRNAProteolysisTag.sequence), 'GeneticCode', knowledgeBase.translationTable)];
            baseCount = ProteinMonomer.computeBaseCount(this.proteolysisTagAASequence, 21, 1:21, false);
            this.proteolysisTagMolecularWeight  = ProteinMonomer.computeMolecularWeight(baseCount, this.proteolysisTagLength, this.aminoAcidMolecularWeights);
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemory(this, numTimePoints)
            this.boundMRNAs            = zeros(0, 1, numTimePoints);
            this.nascentMonomerLengths = zeros(0, 1, numTimePoints);
            this.proteolysisTagLengths = zeros(0, 1, numTimePoints);
            this.abortedPolypeptides   = zeros(0, 3, numTimePoints);
        end
    end
    
    methods
        function initialize(this)
            this.abortedPolypeptides = zeros(0, 3);
        end
    end
    
    %getters
    methods
        function value = get.totalBaseCounts(this)
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            
            value = zeros(1, 21);
            
            for i = 1:size(this.boundMRNAs, 1)
                sequence = [];
                if this.nascentMonomerLengths(i) > 0
                    sequence = [sequence this.monomerAASequences{this.boundMRNAs(i)}(1:this.nascentMonomerLengths(i))]; %#ok<*AGROW>
                end
                if this.proteolysisTagLengths(i) > 0
                    sequence = [sequence this.proteolysisTagAASequence(1:this.proteolysisTagLengths(i))]; %#ok<*AGROW>
                end
                if ~isempty(sequence)
                    value = value + ProteinMonomer.computeBaseCount(sequence, 21, 1:21, true);
                end
            end
            
            abortedSeqs = this.abortedSequences;
            for i = 1:numel(abortedSeqs)
                value = value + ProteinMonomer.computeBaseCount(abortedSeqs{i}, 21, 1:21, true);
            end
        end
        
        %Maximum translation progress of each mRNA
        function value = get.mRNAMaxNascentMonomerLength(this)
            value = zeros(length(this.monomerLengths), 1, size(this.boundMRNAs, 3));
            for i = 1:size(this.boundMRNAs, 1)
                for k = 1:size(this.boundMRNAs, 3)
                    boundMRNA = this.boundMRNAs(i,1,k);
                    if boundMRNA < 1
                        continue; 
                    end
                    value(boundMRNA, 1, k) = max(...
                        value(boundMRNA, 1, k), ...
                        this.nascentMonomerLengths(i, 1, k));
                end
            end
        end
        
        function value = get.abortedSequences(this)
            value = cell(size(this.abortedPolypeptides, 1), 1);
            for i = 1:size(this.abortedPolypeptides, 1)
                if this.abortedPolypeptides(i, 1) == 0
                    continue;
                end
                value{i} = [
                    this.monomerAASequences{this.abortedPolypeptides(i, 1)}(1:this.abortedPolypeptides(i, 2)) ...
                    this.proteolysisTagAASequence(1:this.abortedPolypeptides(i, 3))
                    ];
            end
        end
        
        function value = get.abortedSequenceLengths(this)
            value = sum(this.abortedPolypeptides(:, 2:3), 2);
        end
        
        function value = get.dryWeight(this)
            % import classes
            import edu.stanford.covert.util.ConstantUtil;
            
            seq = [];
            nSeq = 0;
            
            %weight of nascent polypeptides
            for i = 1:size(this.boundMRNAs, 1)
                tmpSeq = [];
                if this.nascentMonomerLengths(i) > 0
                    tmpSeq = [tmpSeq ...
                        this.monomerAASequences{this.boundMRNAs(i)}(1:this.nascentMonomerLengths(i))]; %#ok<*AGROW>
                end
                if this.proteolysisTagLengths(i) > 0
                    tmpSeq = [tmpSeq ...
                        this.proteolysisTagAASequence(1:this.proteolysisTagLengths(i))]; %#ok<*AGROW>
                end
                if ~isempty(tmpSeq)
                    seq = [seq tmpSeq];
                    nSeq = nSeq + 1;
                end
            end
            
            % weight of proteolysis tagged monomers
            for i = 1:size(this.abortedPolypeptides, 1)
                if this.abortedPolypeptides(i, 1) ~= 0
                    nSeq = nSeq + 1;
                    seq = [seq ...
                        this.monomerAASequences{this.abortedPolypeptides(i, 1)}(1:this.abortedPolypeptides(i, 2)) ...
                        this.proteolysisTagAASequence(1:this.abortedPolypeptides(i, 3))
                        ];
                end
            end
            
            value = (...
                + [
                sum(seq == 'A')
                sum(seq == 'R')
                sum(seq == 'N')
                sum(seq == 'D')
                sum(seq == 'C')
                sum(seq == 'Q')
                sum(seq == 'E')
                sum(seq == 'G')
                sum(seq == 'H')
                sum(seq == 'I')
                sum(seq == 'L')
                sum(seq == 'K')
                sum(seq == 'M')
                sum(seq == 'F')
                sum(seq == 'P')
                sum(seq == 'S')
                sum(seq == 'T')
                sum(seq == 'W')
                sum(seq == 'Y')
                sum(seq == 'V')
                ]' * this.aminoAcidMolecularWeights(1:20) ...
                - (numel(seq) - nSeq) * (ConstantUtil.elements.H + 2 * ConstantUtil.elements.O) ...
                ) / ConstantUtil.nAvogadro;
        end
    end
end