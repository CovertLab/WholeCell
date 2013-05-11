%Metabolite
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 1/5/2011
classdef Metabolite < edu.stanford.covert.cell.sim.MoleculeCountState
    %Constants
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity'
            'seed'
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'molecularWeights'
            'baseCounts'
            'lengths'
            'halfLives'
            'compartments'
            'meanNTPConcentration'
            'meanNDPConcentration'
            'meanNMPConcentration'
            };
        fittedConstantNames     = {   %names of process properties that are considered fitted constants, and should be stored with the simulation as such
            'biomassComposition'
            'biomassProduction'
            'byproducts'
            };
        stateNames              = {   %names of properties which are part of the simulation's state
            'counts'
            'processRequirements'   
            'processAllocations'
            'processUsages'
            };
        dependentStateNames     = {}; %names of properties which can be calculated from the simulation's state
    end
    
    %indices
    properties
        hydrophobicIndexs    %indices of hydrophobic counts
        aminoAcidIndexs      %amino acids indices within counts
        dnmpIndexs           %dAMP, dCMP, dGMP, dTMP indices within counts
        dndpIndexs           %dADP, dCDP, dGDP, dTDP indices within counts
        dntpIndexs           %dATP, dCTP, dGTP, dTTP indices within counts
        nmpIndexs            %AMP, CMP, GMP, UMP indices within counts
        ndpIndexs            %ADP, CDP, GDP, UDP indices within counts
        ntpIndexs            %ATP, CTP, GTP, UTP indices within counts
        ampIndexs            %AMP index within counts
        adpIndexs            %ADP index within counts
        atpIndexs            %ATP index within counts
        phosphateIndexs      %Pi index within counts
        diphosphateIndexs    %PPi index within counts
        waterIndexs          %H2O index within counts
        hydrogenIndexs       %H index within counts
        dr5pIndexs           %DR5P index within counts
        unmodifiedBaseIndexs %AD, CSN, GN, THY indices within counts
        m6ADIndexs           %m6AD index within counts        
        lipidIndexs          %lipid indices within counts
        polyamineIndexs      %polyamine indices within counts
        carbohydrateIndexs   %carbohydrate indices within counts
        vitaminIndexs        %vitamin indices within counts
        ionIndexs            %ion indices within counts
        antibioticIndexs     %antibiotic indices within counts
    end
    
    %constants
    properties
        setCounts                                     %experimentally set metabolic conditions (concentration, time when concentration is applied)
        
        experimentalBiomassComposition                %experimentally observed biomass composition (molecules / cell)
        experimentalBiomassCompositionWeightFractions %experimentally observed biomass composition (weight fraction)
        experimentalBiomassCompositionMolFractions    %experimentally observed biomass composition (mol fraction)
        biomassComposition                            %metabolite composition of biomass (mollecules / cell) (metabolites X compartments)
        biomassProduction                             %biomass production
        byproducts                                    %expected energy and other byproducts of networks (ADP, AMP, PPi, Pi, H)
        processBiomassProduction                      %expected metabolic demand of each process
        processByproduct                              %expected metabolic byproducts produced by each process
        
        meanNTPConcentration                          %mean ATP, GTP concentration (mM)
        meanNDPConcentration                          %mean ADP, GDP concentration (mM)
        meanNMPConcentration                          %mean AMP, CMP, GMP, UMP concentration (mM)
        
        processWholeCellModelIDs                      %IDs of processes
    end
    
    %constants
    properties (Dependent)
        metabolismProduction       %biomass production
        
        experimentalNMPComposition %experimentally determined mol fractions of NMPs
        experimentalAAComposition  %experimentally determined mol fractions of AAs
        
        nmpComposition             %mol fractions of NMPs
        aaComposition              %mol fractions of AAs
        
        wetWeight                  %water weight (g)
    end
    
    %state
    properties
        processRequirements        %instantaneous metabolic requirements of each process
        processAllocations         %instantaneous metabolites allocated to each process
        processUsages              %instantaneous metabolic usages of each process
    end
    
    %constructor
    methods
        function this = Metabolite(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.MoleculeCountState(wholeCellModelID, name);
        end
    end
    
    methods
        function initializeConstants(this, knowledgeBase, simulation)
            import edu.stanford.covert.cell.sim.constant.Condition;
            import edu.stanford.covert.util.ConstantUtil;
            
            this.initializeConstants@edu.stanford.covert.cell.sim.MoleculeCountState(knowledgeBase, simulation);
            
            %ids, names, weights
            this.wholeCellModelIDs = {knowledgeBase.metabolites.wholeCellModelID}';
            this.names             = {knowledgeBase.metabolites.name}';
            this.molecularWeights  = [knowledgeBase.metabolites.molecularWeight]';
            
            %indices
            this.hydrophobicIndexs    = find([knowledgeBase.metabolites.hydrophobic]);
            this.aminoAcidIndexs      = this.getIndexs({'ALA'; 'ARG'; 'ASN'; 'ASP'; 'CYS'; 'GLN'; 'GLU'; 'GLY'; 'HIS'; 'ILE'; 'LEU'; 'LYS'; 'MET'; 'PHE'; 'PRO'; 'SER'; 'THR'; 'TRP'; 'TYR'; 'VAL'; 'FMET'});                                                           
            this.dnmpIndexs           = this.getIndexs({'DAMP'; 'DCMP'; 'DGMP'; 'DTMP'});
            this.dndpIndexs           = this.getIndexs({'DADP'; 'DCDP'; 'DGDP'; 'DTDP'});
            this.dntpIndexs           = this.getIndexs({'DATP'; 'DCTP'; 'DGTP'; 'DTTP'});
            this.nmpIndexs            = this.getIndexs({'AMP';  'CMP';  'GMP';  'UMP'});
            this.ndpIndexs            = this.getIndexs({'ADP';  'CDP';  'GDP';  'UDP'});
            this.ntpIndexs            = this.getIndexs({'ATP';  'CTP';  'GTP';  'UTP'});
            this.ampIndexs            = this.getIndexs({'AMP'});
            this.adpIndexs            = this.getIndexs({'ADP'});
            this.atpIndexs            = this.getIndexs({'ATP'});
            this.phosphateIndexs      = this.getIndexs({'PI'});
            this.diphosphateIndexs    = this.getIndexs({'PPI'});
            this.waterIndexs          = this.getIndexs({'H2O'});
            this.hydrogenIndexs       = this.getIndexs({'H'});
            this.dr5pIndexs           = this.getIndexs({'DR5P'});
            this.unmodifiedBaseIndexs = this.getIndexs({'AD'; 'CSN'; 'GN'; 'THY'});
            this.m6ADIndexs           = this.getIndexs({'m6AD'});
            this.lipidIndexs          = find(strcmp({knowledgeBase.metabolites.category}', 'lipid'));
            this.polyamineIndexs      = find(strcmp({knowledgeBase.metabolites.category}', 'polyamine'));
            this.carbohydrateIndexs   = find(strcmp({knowledgeBase.metabolites.category}', 'carbohydrate'));
            this.vitaminIndexs        = find(strcmp({knowledgeBase.metabolites.category}', 'vitamin'));
            this.ionIndexs            = find(strcmp({knowledgeBase.metabolites.category}', 'ion'));
            this.antibioticIndexs     = find(strcmp({knowledgeBase.metabolites.category}', 'antibiotic'));
            
            %media
            this.setCounts = zeros(0, 6);
            for i = 1:knowledgeBase.numMetabolites
                if isempty(knowledgeBase.metabolites(i).mediaCompartments); continue; end;
                
                value = zeros(length(knowledgeBase.metabolites(i).mediaCompartments), 6);
                value(:, Condition.objectIndexs)      = knowledgeBase.metabolites(i).idx;
                value(:, Condition.compartmentIndexs) = [knowledgeBase.metabolites(i).mediaCompartments.idx];
                value(:, Condition.valueIndexs)       = knowledgeBase.metabolites(i).mediaConcentrations;
                value(:, Condition.initialTimeIndexs) = knowledgeBase.metabolites(i).mediaInitialTimes;
                value(:, Condition.finalTimeIndexs)   = knowledgeBase.metabolites(i).mediaFinalTimes;
                value(:, Condition.objectCompartmentIndexs) = sub2ind([knowledgeBase.numMetabolites knowledgeBase.numCompartments],...
                    value(:, Condition.objectIndexs),...
                    value(:, Condition.compartmentIndexs));
                
                this.setCounts = [this.setCounts; value];
            end
            
            this.setCounts(:, Condition.valueIndexs) = round(...
                this.setCounts(:, Condition.valueIndexs) * ...
                edu.stanford.covert.util.ConstantUtil.nAvogadro / 1000 * simulation.state('Geometry').chamberVolume);
            
            %biomass
            this.experimentalBiomassComposition = -knowledgeBase.biomassComposition;
            this.experimentalBiomassCompositionMolFractions = sum(this.experimentalBiomassComposition, 2) / sum(this.experimentalBiomassComposition(:));
            this.experimentalBiomassCompositionWeightFractions = ...
                sum(this.experimentalBiomassCompositionMolFractions, 2) .* this.molecularWeights / ...
                (sum(this.experimentalBiomassCompositionMolFractions, 2)' * this.molecularWeights);
            
            %processes
            this.processWholeCellModelIDs = {knowledgeBase.processes.wholeCellModelID}';
        end
        
        function allocateMemory(this, numTimePoints)
            import edu.stanford.covert.util.SparseMat;
            
            this.allocateMemory@edu.stanford.covert.cell.sim.MoleculeCountState(numTimePoints);
            
            nProcesses = numel(this.processWholeCellModelIDs);
            this.processRequirements = SparseMat([], [], [numel(this.wholeCellModelIDs) * this.compartment.count  nProcesses  numTimePoints]);
            this.processAllocations = SparseMat([], [], [numel(this.wholeCellModelIDs) * this.compartment.count  nProcesses  numTimePoints]);
            this.processUsages = SparseMat([], [], [numel(this.wholeCellModelIDs) * this.compartment.count  nProcesses  numTimePoints]);
        end
    end
    
    %getters
    methods
        function value = get.metabolismProduction(this)
            value = this.biomassProduction - this.byproducts;
        end
        
        %experimental composition
        function value = get.experimentalNMPComposition(this)
            value = this.experimentalBiomassCompositionMolFractions(this.nmpIndexs);
            value = value / sum(value);
        end
        
        function value = get.experimentalAAComposition(this)
            value = this.experimentalBiomassCompositionMolFractions(this.aminoAcidIndexs);
            value = value / sum(value);
        end
        
        %composition
        function value = get.nmpComposition(this)
            value = this.biomassComposition(this.nmpIndexs, :);
            value = value / sum(value(:));
        end
        
        function value = get.aaComposition(this)
            value = sum(this.biomassComposition(this.aminoAcidIndexs, :), 2);
            value = value / sum(value(:));
        end
        
        %weight
        function value = calcDryWeight(this)
            if size(this.counts, 3) == 1
                value = this.molecularWeights' * this.counts;
            else
                value = multiprod(this.molecularWeights', this.counts, [1 2], [1 2]);
            end
            value = ...
                + value / edu.stanford.covert.util.ConstantUtil.nAvogadro ...
                - this.wetWeight;
        end
        
        function value = get.wetWeight(this)
            value = this.counts(this.waterIndexs, :, :) * this.molecularWeights(this.waterIndexs) / edu.stanford.covert.util.ConstantUtil.nAvogadro;
        end
    end
end
