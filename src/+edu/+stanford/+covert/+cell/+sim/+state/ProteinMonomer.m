%ProteinMonomer
%- nascent
%- processed
%- folded
%- mature
%- bound
%
%Translation    NTPs->nascent
%Processing     nascent->processedI->processedII->folded->mature
%Modification   mature->mature
%Misfolding     mature->misfolded
%               inactivated->misfolded
%Refolding      misfolded->mature
%Activation     inactivated->mature
%Inactivation   mature->inactivated
%Damage         Damaged complex->damaged monomer
%Enzymatic Use  mature->bound
%               bound->mature (upon completion)
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 1/5/2011
classdef ProteinMonomer < edu.stanford.covert.cell.sim.MoleculeCountState
    %Constants
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'molecularWeights';
            'baseCounts';
            'lengths';            
            'compartments'
            'minimumAverageExpression'
            'macromoleculeStateInitializationVariation'
            };
        fittedConstantNames     = {   %names of process properties that are considered fitted constants, and should be stored with the simulation as such
            'halfLives';
            };
        stateNames              = {   %names of properties which are part of the simulation's state
            'counts'
            };
        dependentStateNames     = {}; %names of properties which can be calculated from the simulation's state
    end
	    
    %indices
    properties
        nascentIndexs           %index within monomers
        processedIIndexs        %index within monomers
        processedIIIndexs       %index within monomers
        signalSequenceIndexs    %index within monomers
        foldedIndexs            %index within monomers
        matureIndexs            %index within monomers
        inactivatedIndexs       %index within monomers
        boundIndexs             %index within monomers
        misfoldedIndexs         %index within monomers
        damagedIndexs           %index within monomers
        
        translationFactorIndexs %index within matureIndexs
    end
    
    %constants
    properties
        macromoleculeStateInitializationVariation          %Toggles amount of variation in initial state; see intializeState
        minimumAverageExpression                           %minimum average monomer expression
    end
    
    %references to objects
    properties
        chromosome
        ribosome
    end
    
    %constructor
    methods
        function this = ProteinMonomer(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.MoleculeCountState(wholeCellModelID, name);
        end
    end
    
    methods
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.MoleculeCountState(simulation);
            this.chromosome = simulation.state('Chromosome');
            this.ribosome = simulation.state('Ribosome');
        end
    end
    
    methods
        function initializeConstants(this, knowledgeBase, simulation)
            this.initializeConstants@edu.stanford.covert.cell.sim.MoleculeCountState(knowledgeBase, simulation);
            
            numMonomers = knowledgeBase.numProteinMonomers;
            
            this.nascentIndexs        = (1:numMonomers)';
            this.processedIIndexs     = (1:numMonomers)' + this.nascentIndexs(end);
            this.processedIIIndexs    = (1:numMonomers)' + this.processedIIndexs(end);
            this.signalSequenceIndexs = (1:numMonomers)' + this.processedIIIndexs(end);
            this.foldedIndexs         = (1:numMonomers)' + this.signalSequenceIndexs(end);
            this.matureIndexs         = (1:numMonomers)' + this.foldedIndexs(end);
            this.inactivatedIndexs    = (1:numMonomers)' + this.matureIndexs(end);
            this.boundIndexs          = (1:numMonomers)' + this.inactivatedIndexs(end);
            this.misfoldedIndexs      = (1:numMonomers)' + this.boundIndexs(end);
            this.damagedIndexs        = (1:numMonomers)' + this.misfoldedIndexs(end);
            
            this.wholeCellModelIDs = repmat({knowledgeBase.proteinMonomers.wholeCellModelID}', 10, 1);
            this.names = repmat({knowledgeBase.proteinMonomers.name}', 10, 1);
            this.molecularWeights = [...
                knowledgeBase.proteinMonomers.molecularWeight ...
                knowledgeBase.proteinMonomers.processedISequenceMolecularWeight  ...
                knowledgeBase.proteinMonomers.processedIISequenceMolecularWeight  ...
                knowledgeBase.proteinMonomers.signalSequenceMolecularWeight ...
                knowledgeBase.proteinMonomers.foldedSequenceMolecularWeight ...
                repmat([knowledgeBase.proteinMonomers.matureSequenceMolecularWeight], 1, 5)]';
            this.baseCounts = [
                reshape([knowledgeBase.proteinMonomers.baseCount],                      [], numMonomers)';
                reshape([knowledgeBase.proteinMonomers.processedISequenceBaseCount],    [], numMonomers)';
                reshape([knowledgeBase.proteinMonomers.processedIISequenceBaseCount],   [], numMonomers)';
                reshape([knowledgeBase.proteinMonomers.signalSequenceBaseCount],        [], numMonomers)';
                reshape([knowledgeBase.proteinMonomers.foldedSequenceBaseCount],        [], numMonomers)';
                repmat(reshape([knowledgeBase.proteinMonomers.matureSequenceBaseCount], [], numMonomers)', 5, 1)];
            this.lengths = [...
                knowledgeBase.proteinMonomers.sequenceLength  ...
                knowledgeBase.proteinMonomers.processedISequenceLength  ...
                knowledgeBase.proteinMonomers.processedIISequenceLength  ...
                knowledgeBase.proteinMonomers.signalSequenceLength  ...
                knowledgeBase.proteinMonomers.foldedSequenceLength  ...
                repmat([knowledgeBase.proteinMonomers.matureSequenceLength], 1, 5)]';
            this.compartments = double(repmat(knowledgeBase.proteinMonomerCompartments, 10, 1));
            this.compartments(this.signalSequenceIndexs) = this.compartment.cytosolIndexs;
            this.halfLives = [...
                knowledgeBase.proteinMonomers.halfLife ...
                knowledgeBase.proteinMonomers.processedISequenceHalfLife ...
                knowledgeBase.proteinMonomers.processedIISequenceHalfLife ...
                knowledgeBase.proteinMonomers.signalSequenceHalfLife ...
                knowledgeBase.proteinMonomers.foldedSequenceHalfLife ...
                repmat([knowledgeBase.proteinMonomers.matureSequenceHalfLife], 1, 3) ...
                zeros(1, numMonomers)  ...
                zeros(1, numMonomers)]';
            
            this.halfLives(~ismember(this.compartments, [this.compartment.cytosolIndexs; this.compartment.terminalOrganelleCytosolIndexs])) = Inf;
            this.halfLives(this.signalSequenceIndexs) = 0;
            
            this.translationFactorIndexs = this.getIndexs('MG_026_MONOMER');
        end
    end
    
    methods
        function notUpdatingProteins = updateExternalState(this, deltaProteins, proteinIsDegraded)
            c = this.chromosome;
            
            notUpdatingProteins = zeros(size(deltaProteins));
            deltaBoundProteins = -deltaProteins(this.boundIndexs, this.compartment.cytosolIndexs);
            
            %translation factors
            deltaBoundProteins(this.translationFactorIndexs) = 0;
            
            %update chromosome state
            idxs = find(deltaBoundProteins);
            [posStrnds, proteins] = find(c.monomerBoundSites);
            
            chrReleasePosStrnds = zeros(0, 2);
            for i = 1:numel(idxs)
                if sum(proteins == idxs(i)) < deltaBoundProteins(idxs(i))
                    throw(MException('ProteinMonomer:error', 'Error updating external state'))
                end
                chrReleasePosStrnds = [
                    chrReleasePosStrnds;
                    this.randStream.randomlySelectNRows(posStrnds(proteins == idxs(i), :), deltaBoundProteins(idxs(i)))
                    ]; %#ok<AGROW>
            end
            
            releasedMonomers = c.setRegionProteinUnbound(...
                chrReleasePosStrnds, 1, idxs, [], ...
                false, false, false, proteinIsDegraded);
            if ~isequal(deltaBoundProteins(idxs), releasedMonomers)
                throw(MException('ProteinMonomer:error', 'Chromosomally bound proteins impropely released'));
            end
        end
    end
    
    %helper methods
    methods
        function value = getIndexs(this, wholeCellModelIDs)
            [~, value] = ismember(wholeCellModelIDs, this.wholeCellModelIDs(this.matureIndexs));
        end
    end
end
