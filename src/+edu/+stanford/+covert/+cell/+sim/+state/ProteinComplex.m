%ProteinComplex
%- nascent
%- mature
%- bound
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 1/5/2011
classdef ProteinComplex < edu.stanford.covert.cell.sim.MoleculeCountState
    
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity'
            'seed'
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'molecularWeights'
            'baseCounts'
            'lengths'
            'compartments'
            'proteinComplexComposition'
            'minimumAverageExpression'
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
        nascentIndexs           %index within complexs
        matureIndexs            %index within complexs
        inactivatedIndexs       %index within complexs
        boundIndexs             %index within complexs
        misfoldedIndexs         %index within complexs
        damagedIndexs           %index within complexs
        
        ribosome30SIndexs       %index within matureIndexs
        ribosome50SIndexs       %index within matureIndexs
        ribosome70SIndexs       %index within matureIndexs
        translationFactorIndexs %index within matureIndexs
        ftsZGTPIndexs           %index within matureIndexs
        ftsZGDPIndexs           %index within matureIndexs
        rnaPolymeraseIndexs     %index within matureIndexs
        replisomeIndexs         %index within matureIndexs
        dnaPolymeraseIndexs     %index within matureIndexs
        dnaAPolymerIndexs       %index within matureIndexs
    end
    
    properties
        proteinComplexComposition %protein complex composition (monomers X complexes X compartments)
        formationProcesses        %indices of proceses where each complex is formed
        minimumAverageExpression  %minimum average expression
    end
    
    %references to objects
    properties
        chromosome
        rnaPolymerase
        ribosome
        ftsZRing
    end
    
    %constructor
    methods
        function this = ProteinComplex(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.MoleculeCountState(wholeCellModelID, name);
        end
    end
    
    methods
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.MoleculeCountState(simulation);            
            this.chromosome = simulation.state('Chromosome');
            this.rnaPolymerase = simulation.state('RNAPolymerase');
            this.ribosome = simulation.state('Ribosome');
            this.ftsZRing = simulation.state('FtsZRing');
        end
    end
    
    methods
        function initializeConstants(this, knowledgeBase, simulation)
            this.initializeConstants@edu.stanford.covert.cell.sim.MoleculeCountState(knowledgeBase, simulation);
            
            numComplexs = knowledgeBase.numProteinComplexs;
            
            this.nascentIndexs     = (1:numComplexs)';
            this.matureIndexs      = (1:numComplexs)' + this.nascentIndexs(end);
            this.inactivatedIndexs = (1:numComplexs)' + this.matureIndexs(end);
            this.boundIndexs       = (1:numComplexs)' + this.inactivatedIndexs(end);
            this.misfoldedIndexs   = (1:numComplexs)' + this.boundIndexs(end);
            this.damagedIndexs     = (1:numComplexs)' + this.misfoldedIndexs(end);
            
            this.wholeCellModelIDs = repmat({knowledgeBase.proteinComplexs.wholeCellModelID}', 6, 1);
            this.names = repmat({knowledgeBase.proteinComplexs.name}', 6, 1);
            this.baseCounts = [
                reshape([knowledgeBase.proteinComplexs.baseCount], [], numComplexs)';
                repmat(reshape([knowledgeBase.proteinComplexs.matureBaseCount], [], numComplexs)', 5, 1)];
            this.molecularWeights = [...
                knowledgeBase.proteinComplexs.molecularWeight  ...
                repmat([knowledgeBase.proteinComplexs.matureMolecularWeight], 1, 5)]';
            this.compartments = double(repmat(knowledgeBase.proteinComplexCompartments, 6, 1));
            this.halfLives = [...
                knowledgeBase.proteinComplexs.halfLife  ...
                repmat([knowledgeBase.proteinComplexs.matureHalfLife], 1, 3) ...
                zeros(1, knowledgeBase.numProteinComplexs) ...
                zeros(1, knowledgeBase.numProteinComplexs)]';
            
            this.proteinComplexComposition = knowledgeBase.proteinComplexAllRNAComposition;
            this.proteinComplexComposition([knowledgeBase.mRNAGenes.idx], :, :) = knowledgeBase.proteinComplexAllMonomerComposition;
            
            tmp = [knowledgeBase.proteinComplexs.complexFormationProcess];
            this.formationProcesses = repmat([tmp.idx]', 6, 1);
            
            this.halfLives(~ismember(this.compartments, [this.compartment.cytosolIndexs; this.compartment.terminalOrganelleCytosolIndexs])) = Inf;
            
            this.ribosome30SIndexs = this.getIndexs('RIBOSOME_30S');
            this.ribosome50SIndexs = this.getIndexs('RIBOSOME_50S');
            this.ribosome70SIndexs = this.getIndexs('RIBOSOME_70S');
            this.translationFactorIndexs = this.getIndexs({
                'MG_089_DIMER'
                'MG_433_DIMER'
                'MG_451_DIMER'
                });
            this.ftsZGTPIndexs = this.getIndexs('MG_224_9MER_GTP');
            this.ftsZGDPIndexs = this.getIndexs('MG_224_9MER_GDP');
            this.rnaPolymeraseIndexs = this.getIndexs({
                'RNA_POLYMERASE'
                'RNA_POLYMERASE_HOLOENZYME'
                });
            this.replisomeIndexs = this.getIndexs({
                'MG_094_HEXAMER'
                'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE'
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX'
                'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE'
                });
            this.dnaPolymeraseIndexs = this.getIndexs({
                'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE'
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX'
                'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE'
                });
            this.dnaAPolymerIndexs = this.getIndexs({
                'MG_469_2MER_1ATP_ADP' %DnaA 2mer-(1)ATP-(1)ADP
                'MG_469_2MER_ATP'      %DnaA-ATP 2mer
                'MG_469_3MER_2ATP_ADP' %DnaA 3mer-(2)ATP-(1)ADP
                'MG_469_3MER_ATP'      %DnaA-ATP 3mer
                'MG_469_4MER_3ATP_ADP' %DnaA 4mer-(3)ATP-(1)ADP
                'MG_469_4MER_ATP'      %DnaA-ATP 4mer
                'MG_469_5MER_4ATP_ADP' %DnaA 5mer-(4)ATP-(1)ADP
                'MG_469_5MER_ATP'      %DnaA-ATP 5mer
                'MG_469_6MER_5ATP_ADP' %DnaA 6mer-(5)ATP-(1)ADP
                'MG_469_6MER_ATP'      %DnaA-ATP 6mer
                'MG_469_7MER_6ATP_ADP' %DnaA 7mer-(6)ATP-(1)ADP
                'MG_469_7MER_ATP'      %DnaA-ATP 7mer
                });
        end
    end
            
    methods
        function notUpdatingProteins = updateExternalState(this, deltaProteins, proteinIsDegraded)
            c = this.chromosome;
            
            notUpdatingProteins = zeros(size(deltaProteins));
            deltaFreeProteins = deltaProteins(this.matureIndexs, this.compartment.cytosolIndexs);
            deltaBoundProteins = deltaProteins(this.boundIndexs, this.compartment.cytosolIndexs);
            
            %update ribosome state
            this.ribosome.releaseRibosome(-deltaBoundProteins(this.ribosome70SIndexs), 0, 0);
            deltaBoundProteins(this.ribosome70SIndexs) = 0;
            
            %update bound translation elongation factors
            deltaBoundProteins(this.translationFactorIndexs) = 0;
            
            %update bound FtsZ state
            notUpdatingProteins(this.boundIndexs([this.ftsZGTPIndexs; this.ftsZGDPIndexs]), this.compartment.cytosolIndexs) = ...
                this.ftsZRing.releaseFtsZ(-deltaBoundProteins([this.ftsZGTPIndexs; this.ftsZGDPIndexs]));
            deltaBoundProteins(this.ftsZGTPIndexs) = 0;
            deltaBoundProteins(this.ftsZGDPIndexs) = 0;
            
            %prevent changes to bound DnaA complexes to that replication
            %initiation not disrupted
            if any(deltaBoundProteins(this.dnaAPolymerIndexs) < 0)
                warning('WholeCell:warning', 'DnaA complex not decayed');
                notUpdatingProteins(this.boundIndexs(this.dnaAPolymerIndexs), this.compartment.cytosolIndexs) = ...
                    -deltaBoundProteins(this.dnaAPolymerIndexs);
                deltaBoundProteins(this.dnaAPolymerIndexs) = 0;
            end
            
            %prevent changes to bound DNA polymerase so that replication
            %not disrupted; issue warning
            if any(deltaBoundProteins(this.replisomeIndexs) < 0)
                warning('WholeCell:warning', 'DNA polymerase not decayed');
                notUpdatingProteins(this.boundIndexs(this.replisomeIndexs), this.compartment.cytosolIndexs) = ...
                    -deltaBoundProteins(this.replisomeIndexs);
                deltaBoundProteins(this.replisomeIndexs) = 0;
            end
            
            %update chromosomally bound proteins
            idxs = find(deltaBoundProteins < 0);
            [posStrnds, proteins] = find(c.complexBoundSites);
            
            chrReleasePosStrnds = zeros(0, 2);
            for i = 1:numel(idxs)
                if sum(proteins == idxs(i)) < -deltaBoundProteins(idxs(i))
                    throw(MException('ProteinComplex:error', 'Error updating external state'))
                end
                chrReleasePosStrnds = [
                    chrReleasePosStrnds;
                    this.randStream.randomlySelectNRows(posStrnds(proteins == idxs(i), :), -deltaBoundProteins(idxs(i)))
                    ]; %#ok<AGROW>
            end
            
            [~, releasedComplexs] = c.setRegionProteinUnbound(...
                chrReleasePosStrnds, 1, [], idxs, ...
                false, false, false, proteinIsDegraded);
            if ~isequal(-deltaBoundProteins(idxs), releasedComplexs)
                throw(MException('ProteinComplex:error', 'Chromosomally bound proteins impropely released'));
            end
            
            %update free RNA polymerase
            if deltaFreeProteins(this.rnaPolymeraseIndexs(1)) && proteinIsDegraded
                this.rnaPolymerase.degradeFreePolymerase(-deltaFreeProteins(this.rnaPolymeraseIndexs(1)));
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
