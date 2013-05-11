%Chromosome
% Integration point for processes which interact with specific
% positions/strands of the cell's chromosome(s).
% - Represents the portion of chromosome(s) accessible to enzymes. That is
%   positions/strands which are NOT
%   - damaged in any way (no gap sites, abasic sites, damaged
%     sugar-phosphates, damaged bases, cross links, strand breaks,
%     or Holliday junctions)
%   - stably bound by enzymes
%   - single stranded
%
% Terminology:
% ==================
%         Site  single base/bond of chromosomes, indicated by strand index and
%               number of bases/bonds along 5'->3' strand from ORI [position X
%               strand]
%       Region  contiguous set of bases/bonds of chromosomes, indicated by start
%               and end positions (bases/bonds along 5'->3' strand from ORI and
%               strand (positive/negative)
%   Accessible  polymerized, not bound by protein, and not damaged
% Inaccessible  not polymerized, bound by protein, or damaged
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010

%TODO
%- more precise ignoreDamageFilter
%- include gapsites in isRegionPolymerized
classdef Chromosome < edu.stanford.covert.cell.sim.CellState
    %constants
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'doubleStrandBreakSeparation';
            'strandBreakClassification_doubleStrandBreakSeparation';
            'strandBreakClassification_segmentLength';
            'sequence';
            'sequenceLen';
            'oriCPosition';
            'terCPosition';
            'transcriptionUnitStartCoordinates';
            'transcriptionUnitLengths';
            'transcriptionUnitStrands';
            'monomerDNAFootprints';
            'complexDNAFootprints';
            'monomerDNAFootprintBindingStrandedness';
            'complexDNAFootprintBindingStrandedness';
            'monomerDNAFootprintRegionStrandedness';
            'complexDNAFootprintRegionStrandedness';
            'reactionBoundMonomer';
            'reactionBoundComplex';
            'reactionMonomerCatalysisMatrix';
            'reactionComplexCatalysisMatrix';
            'reactionThresholds';
            'relaxedBasesPerTurn';
            'equilibriumSuperhelicalDensity';
            'supercoiledSuperhelicalDensityTolerance';
            };
        fittedConstantNames = {};  %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        stateNames = {             %names of properties which are part of the simulation's state
            'polymerizedRegions';
            'linkingNumbers';
            'monomerBoundSites';
            'complexBoundSites';
            'gapSites';
            'abasicSites';
            'damagedSugarPhosphates';
            'damagedBases';
            'intrastrandCrossLinks';
            'strandBreaks';
            'hollidayJunctions';
            'segregated';
            };
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            %'unpolymerizedRegions'
            'singleStrandedRegions'
            'doubleStrandedRegions'
            'geneCopyNumbers'
            %'polymerizedGenes'
            'transcriptionUnitCopyNumbers'
            %'polymerizedTranscriptionUnits'
            %'geneCopyNumbers_Accessible'
            %'transcriptionUnitCopyNumbers_Accessible'
            %'accessibleGenes'
            %'accessibleTranscriptionUnits'
            'ploidy'
            %'linkingNumbers_minFreeEnergy'
            %'supercoils'
            'superhelicalDensity'
            %'supercoiled'
            %'damagedSites'
            %'damagedSites_shifted_incm6AD'
            %'damagedSites_nonRedundant'
            %'damagedSites_excm6AD'
            %'gapSites3'
            %'gapSites5'
            %'abasicSites3'
            %'abasicSites5'
            %'damagedSugarPhosphates3'
            %'damagedSugarPhosphates5'
            %'damagedBases3'
            %'damagedBases5'
            %'strandBreaks3'
            %'strandBreaks5'
            %'intrastrandCrossLinks3'
            %'intrastrandCrossLinks5'
            %'hollidayJunctions3'
            %'hollidayJunctions5'
            %'singleStrandBreaks'
            %'doubleStrandBreaks'
            %'strandBreakClassification'
            %'restrictableMunIRMSites'
            %'hemiunmethylatedMunIRMSites'
            };
        
        nCompartments = 4;              %number of strands
        nChromosomes = 2;               %number of chromosomes
        
        strandIndexs_positive    = [1; 3]; %5'->3' strands
        strandIndexs_negative    = [2; 4]; %3'->5' strands
        
        strandIndexs_ch1         = [1; 2]; %chromosome 1 strands; sort to "mother" cell during division
        strandIndexs_ch2         = [3; 4]; %chromosome 2 strands; sort to "daughter" cell during division
        
        strandIndexs_template    = [1; 4]; %strands which are serving as templates for replication
        strandIndexs_nonTemplate = [2; 3]; %strands which are not serving as templates for replication
        
        strandIndexs_old         = [1; 4]; %strands which (at end the end of the cell cycle) were present at the beginning of the cell cycle
        strandIndexs_new         = [2; 3]; %strands which (at end the end of the cell cycle) were synthesized during the cell cycle
        
        dnaStrandedness_ssDNA = 1; %values of ssDNA footprints within *DNAFootprintStrandedness
        dnaStrandedness_dsDNA = 2; %values of dsDNA footprints within *DNAFootprintStrandedness
        dnaStrandedness_xsDNA = 3; %values of ssDNA/dsDNA footprints within *DNAFootprintStrandedness
        
        strandBreakClassification_index_NB    = 1; %index within strandBreakClassification of NB
        strandBreakClassification_index_SSB   = 2; %index within strandBreakClassification of SSB
        strandBreakClassification_index_SSB_  = 3; %index within strandBreakClassification of SSB+
        strandBreakClassification_index_2SSB  = 4; %index within strandBreakClassification of 2SSB
        strandBreakClassification_index_DSB   = 5; %index within strandBreakClassification of DSB
        strandBreakClassification_index_DSB_  = 6; %index within strandBreakClassification of DSB+
        strandBreakClassification_index_DSB__ = 7; %index within strandBreakClassification of DSB++
    end
    
    %computed ids, names, indices
    properties
        transcriptionUnitWholeCellModelIDs     %whole cell model IDs of transcription units
        transcriptionUnitNames                 %names of transcription units
        
        monomerIndexs_ligase                   %index within ProteinMonomer.matureIndexs of DNA ligase
        complexIndexs_dnaPolymerase            %index within ProteinComplex.matureIndexs of DNA polymerase
        complexIndexs_DisA                     %index within ProteinComplex.matureIndexs of DisA
        complexIndexs_rnaPolymerase            %index within ProteinComplex.matureIndexs of RNA polymerase
        monomerIndexs_replisome                %index within ProteinMonomer.matureIndexs of replication machinery
        complexIndexs_replisome                %index within ProteinComplex.matureIndexs of replication machinery
        
        reactionWholeCellModelIDs              %IDs of bound protein release reactions
        reactionNames                          %names of bound protein release reactions
    end
    
    %constants
    properties
        doubleStrandBreakSeparation                           = 1;   %max separtion in bases between SSB's to be consider a DBSB
        strandBreakClassification_doubleStrandBreakSeparation = 10;  %maximum separation of single strand breaks which is considered double strand break [PUB_0486]
        strandBreakClassification_segmentLength               = 216; %length of segments for which SSB/SSB+/2SSB/DSB/DSB+/DSB++ classification applies [PUB_0486]        
        
        sequence                                %chromosome sequence
        sequenceLen                             %chromosome sequence length (number of bases)
        sequenceGCContent                       %chromosome G/C content
        oriCPosition                            %oriC position
        terCPosition                            %terC position
        
        transcriptionUnitStartCoordinates       %genomic coordinates of transcription units
        transcriptionUnitLengths                %genomic coordinates of transcription units
        transcriptionUnitStrands                %genomic direction of transcription units
        
        monomerDNAFootprints                    %number of bases of DNA each monomer occupies when bound to DNA
        complexDNAFootprints                    %number of bases of DNA each complex occupies when bound to DNA
        monomerDNAFootprintBindingStrandedness  %enumeration of binding strandedness of DNA footprints (ssDNA, dsDNA, ssDNA/dsDNA) using this.dnaStrandedness_*sDNA
        complexDNAFootprintBindingStrandedness  %enumeration of binding strandedness of DNA footprints (ssDNA, dsDNA, ssDNA/dsDNA) using this.dnaStrandedness_*sDNA
        monomerDNAFootprintRegionStrandedness   %enumeration of region strandedness of DNA footprints (ssDNA, dsDNA, ssDNA/dsDNA) using this.dnaStrandedness_*sDNA
        complexDNAFootprintRegionStrandedness   %enumeration of region strandedness of DNA footprints (ssDNA, dsDNA, ssDNA/dsDNA) using this.dnaStrandedness_*sDNA
        
        reactionBoundMonomer                    %protein monomer released by each reaction [reactions X 1]
        reactionBoundComplex                    %protein complex released by each reaction [reactions X 1]
        reactionMonomerCatalysisMatrix          %monomers required to catalyze the release of a bound protein [reactions X monomers]
        reactionComplexCatalysisMatrix          %complexs required to catalyze the release of a bound protein [reactions X complexs]
        reactionThresholds                      %number of proteins required to catalyze each release reaction [reactions X 1]
        
        relaxedBasesPerTurn                     %Number of dna bases per turn for the relaxed LK calculation (10.5)
        equilibriumSuperhelicalDensity          %equilibrium superhelical density; also known as specific linking difference / \sigma_{sp} [-0.06; PUB_0749]
        supercoiledSuperhelicalDensityTolerance %tolerance in superhelical density to be considered supercoiled (0.1)
    end
    
    %state
    properties
        polymerizedRegions       %integers [positions x strands] indicating the start positions of polymerized regions of strands and their lengths
        linkingNumbers           %integers [positions x strands] indicating the current linking number of each double-stranded region
        
        monomerBoundSites        %indices [positions x strands] indicating start positions of protein monomers bound to DNA bases
        complexBoundSites        %indices [positions x strands] indicating start positions of macromolecular complexes bound to DNA bases
        
        gapSites                 %boolean [positions x strands] indicating positions of gap sites
        abasicSites              %boolean [positions x strands] indicating positions of abasic sites
        damagedSugarPhosphates   %indices [positions x strands] indicating metabolite identity of damaged sugar-phosphates
        damagedBases             %indices [positions x strands] indicating metabolite identity of damaged bases
        intrastrandCrossLinks    %boolean [positions x strands] indicating metabolite identity of intrastrand cross links in DNA
        strandBreaks             %boolean [positions x strands] indicating positions of strand breaks in strands of DNA
        hollidayJunctions        %boolean [positions x strands] indicating positions of holliday junctions
        
        segregated               %boolean indicating whether or not the chromsomes are segregated
    end
    
    %properties to keep track of whether or not the dependent properties need to
    %be recomputed
    properties (SetAccess = protected)
        validated
        
        validated_polymerizedRegions
        validated_linkingNumbers
        validated_proteinBoundSites
        validated_damaged
        validated_gapSites
        validated_abasicSites
        validated_damagedSugarPhosphates
        validated_damagedBases
        validated_intrastrandCrossLinks
        validated_strandBreaks
        validated_hollidayJunctions
        validated_segregated
        
        validated_unpolymerizedRegions
        validated_singleStrandedRegions
        validated_doubleStrandedRegions
        validated_geneCopyNumbers
        validated_ploidy
        validated_polymerizedGenes
        validated_transcriptionUnitCopyNumbers
        validated_polymerizedTranscriptionUnits
        validated_geneCopyNumbers_Accessible
        validated_transcriptionUnitCopyNumbers_Accessible
        validated_accessibleGenes
        validated_accessibleTranscriptionUnits
        validated_linkingNumbers_minFreeEnergy
        validated_supercoils
        validated_superhelicalDensity
        validated_supercoiled
        validated_damagedSites
        validated_damagedSites_shifted_incm6AD
        validated_damagedSites_nonRedundant
        validated_damagedSites_excm6AD
        validated_gapSites3
        validated_gapSites5
        validated_abasicSites3
        validated_abasicSites5
        validated_damagedSugarPhosphates3
        validated_damagedSugarPhosphates5
        validated_damagedBases3
        validated_damagedBases5
        validated_strandBreaks3
        validated_strandBreaks5
        validated_intrastrandCrossLinks3
        validated_intrastrandCrossLinks5
        validated_hollidayJunctions3
        validated_hollidayJunctions5
        validated_singleStrandBreaks
        validated_doubleStrandBreaks
        validated_strandBreakClassification
        validated_munIRMSiteMethylationStatus
        validated_munIRMSiteRestrictionStatus
        validated_dryWeight
    end
    
    %dependent local state (implemented as dependent property for convenience)
    properties (SetAccess = protected)
        unpolymerizedRegions         %integers indicating the start positions of unpolymerized regions (ie. not yet replicated) of strands and their lengths
        singleStrandedRegions
        doubleStrandedRegions
        geneCopyNumbers
        polymerizedGenes
        transcriptionUnitCopyNumbers
        polymerizedTranscriptionUnits
        geneCopyNumbers_Accessible
        transcriptionUnitCopyNumbers_Accessible
        accessibleGenes
        accessibleTranscriptionUnits
        ploidy
        
        linkingNumbers_minFreeEnergy %free energy mininum linking number
        supercoils                   %difference between linkingNumbers and free energy mininum linking number
        superhelicalDensity          %supercoils / free energy minimum linking number
        supercoiled                  %boolean (1 x 2) indicating whether or not each chromosome is properly supercoiled (eg. within some tolerance of the free energy minimum)
        
        damagedSites                 %integers (genome length x 4) indicating identities of bases which are damaged or which are adjacent to damaged bonds
        damagedSites_shifted_incm6AD %integers (genome length x 4) indicating identities of damaged bases/sites
        damagedSites_nonRedundant    %integers (genome length x 4) indicating identities of damaged bases/sites
        damagedSites_excm6AD         %integers (genome length x 4) indicating identities of damaged bases/sites
        gapSites3                    %boolean (genome length x 4) indicating positions of gap sites 3' to bases
        gapSites5                    %boolean (genome length x 4) indicating positions of gap sites 5' to bases
        abasicSites3                 %boolean (genome length x 4) indicating positions of abasic sites 3' to bases
        abasicSites5                 %boolean (genome length x 4) indicating positions of abasic sites 5' to bases
        damagedSugarPhosphates3      %integers (genome length x 4) indicating indices of damaged sugar-phosphate 3' to bases
        damagedSugarPhosphates5      %integers (genome length x 4) indicating indices of damaged sugar-phosphate 5' to bases
        damagedBases3                %integers (genome length x 4) indicating indices of damaged bases 3' to bases
        damagedBases5                %integers (genome length x 4) indicating indices of damaged bases 5' to bases
        strandBreaks3                %boolean (genome length x 4) indicating positions of strand breaks 3' to bases
        strandBreaks5                %boolean (genome length x 4) indicating positions of strand breaks 5' to bases
        intrastrandCrossLinks3       %boolean (genome length x 4) indicating positions of intrastrand cross links 3' to bases
        intrastrandCrossLinks5       %boolean (genome length x 4) indicating positions of intrastrand cross links 5' to bases
        hollidayJunctions3           %boolean (genome length x 4) indicating positions of holliday junctions 3' to bases
        hollidayJunctions5           %boolean (genome length x 4) indicating positions of holliday junctions 5' to bases
        
        singleStrandBreaks           %boolean (genome length x 4) indicating positions of single strand breaks (strand breaks excluding double strand breaks and strand breaks adjacent to gapSites)
        doubleStrandBreaks           %boolean (genome length x 4) indicating positions of double strand breaks        
        strandBreakClassification
        
        restrictableMunIRMSites      %boolean (genome length x 4) indicating positions of restrictable MunI R/M sites because they aren't methylated
        hemiunmethylatedMunIRMSites  %boolean (genome length x 4) indicating positions of hemi-unmethylated MunI R/M sites
        
        dryWeight                    %dry weight of this class' state properties
    end
    
    %references to objects
    properties
        compartment   %compartments
        gene          %genes
        
        metabolite    %metabolites
        transcript    %transcript
        rnaPolymerase %RNA polymerase
        
        dnaRepair     %DNA repair process
    end
    
    %constructor
    methods
        function this = Chromosome(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    methods
        function storeObjectReferences(this, simulation)
            this.compartment = simulation.compartment;
            this.gene = simulation.gene;
            this.metabolite = simulation.state('Metabolite');
            this.transcript = simulation.state('Transcript');
            this.rnaPolymerase = simulation.state('RNAPolymerase');
            this.dnaRepair = simulation.process('DNARepair');
        end
    end
    
    methods
        function initializeConstants(this, knowledgeBase, simulation)
            this.initializeConstants@edu.stanford.covert.cell.sim.CellState(knowledgeBase, simulation);
            
            %% import classes
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            
            %% sequence
            this.sequence = ChromosomeSequence(knowledgeBase.genome.sequence);
            this.sequenceLen = size(this.sequence, 1);
            this.sequenceGCContent = getGCContent(this.sequence);
            oriC = findobj(knowledgeBase.genomeFeatures, 'wholeCellModelID', 'oriC');
            terC = findobj(knowledgeBase.genomeFeatures, 'wholeCellModelID', 'terC');
            this.oriCPosition = oriC.startCoordinate;
            this.terCPosition = terC.startCoordinate;
            
            %Ensure that relaxedBasesPerTurn is defined such that the relaxed
            %linking number will be a whole number for this organism.
            this.relaxedBasesPerTurn = length(this.sequence) * (1 + this.equilibriumSuperhelicalDensity) / ...
                ceil(length(this.sequence) / this.relaxedBasesPerTurn * (1 + this.equilibriumSuperhelicalDensity));
            
            %% transcription units
            this.transcriptionUnitWholeCellModelIDs = {knowledgeBase.transcriptionUnits.wholeCellModelID}';
            this.transcriptionUnitNames             = {knowledgeBase.transcriptionUnits.name}';
            this.transcriptionUnitStartCoordinates  = [knowledgeBase.transcriptionUnits.startCoordinate]';
            this.transcriptionUnitLengths           = [knowledgeBase.transcriptionUnits.sequenceLength]';
            this.transcriptionUnitStrands           = 2-[knowledgeBase.transcriptionUnits.direction]';
            
            %% proteins
            this.monomerIndexs_ligase        = simulation.state('ProteinMonomer').getIndexs({'MG_254_MONOMER'});
            this.complexIndexs_dnaPolymerase = simulation.state('ProteinComplex').getIndexs({'DNA_POLYMERASE_CORE'});
            this.complexIndexs_DisA          = simulation.state('ProteinComplex').getIndexs({'MG_105_OCTAMER'});
            this.complexIndexs_rnaPolymerase = sort(simulation.state('ProteinComplex').getIndexs({'RNA_POLYMERASE'; 'RNA_POLYMERASE_HOLOENZYME'}));
            this.monomerIndexs_replisome     = zeros(0, 1);
            this.complexIndexs_replisome     = sort(simulation.state('ProteinComplex').getIndexs({
                'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE'
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX'
                'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE'
                'MG_001_DIMER'
                'MG_094_HEXAMER'
                }));
            
            this.monomerDNAFootprints = ceil([knowledgeBase.proteinMonomers.dnaFootprint]');
            this.complexDNAFootprints = ceil([knowledgeBase.proteinComplexs.dnaFootprint]');
            
            strandedValues = cell(3, 1);
            strandedValues{this.dnaStrandedness_ssDNA} = 'ssDNA';
            strandedValues{this.dnaStrandedness_dsDNA} = 'dsDNA';
            strandedValues{this.dnaStrandedness_xsDNA} = 'Either';
            
            monomerStrandedness = {knowledgeBase.proteinMonomers.dnaFootprintBindingStrandedness}';
            complexStrandedness = {knowledgeBase.proteinComplexs.dnaFootprintBindingStrandedness}';
            monomerStrandedness(cellfun(@isempty, monomerStrandedness)) = {'dsDNA'};
            complexStrandedness(cellfun(@isempty, complexStrandedness)) = {'dsDNA'};
            [~, this.monomerDNAFootprintBindingStrandedness] = ismember(monomerStrandedness, strandedValues);
            [~, this.complexDNAFootprintBindingStrandedness] = ismember(complexStrandedness, strandedValues);
            
            monomerStrandedness = {knowledgeBase.proteinMonomers.dnaFootprintRegionStrandedness}';
            complexStrandedness = {knowledgeBase.proteinComplexs.dnaFootprintRegionStrandedness}';
            monomerStrandedness(cellfun(@isempty, monomerStrandedness)) = {'dsDNA'};
            complexStrandedness(cellfun(@isempty, complexStrandedness)) = {'dsDNA'};
            [~, this.monomerDNAFootprintRegionStrandedness] = ismember(monomerStrandedness, strandedValues);
            [~, this.complexDNAFootprintRegionStrandedness] = ismember(complexStrandedness, strandedValues);
            
            %% proteins which have ability to release other proteins bound to chromosomes
            state = findobj(knowledgeBase.states, 'wholeCellModelID', this.wholeCellModelID);
            reactions = state.reactions;
            reactionGlobalIndexs = [reactions.idx]';
            this.reactionWholeCellModelIDs = {reactions.wholeCellModelID}';
            this.reactionNames = {reactions.name}';
            
            this.reactionBoundMonomer = zeros(size(reactions));
            this.reactionBoundComplex = zeros(size(reactions));
            [i, j] = find(max(0, knowledgeBase.reactionProteinMonomerStoichiometryMatrix(:, reactionGlobalIndexs, this.compartment.cytosolIndexs)));
            this.reactionBoundMonomer(j) = i;
            [i, j] = find(max(0, knowledgeBase.reactionProteinComplexStoichiometryMatrix(:, reactionGlobalIndexs, this.compartment.cytosolIndexs)));
            this.reactionBoundComplex(j) = i;
            
            this.reactionMonomerCatalysisMatrix = max(0, -knowledgeBase.reactionProteinMonomerStoichiometryMatrix(:, reactionGlobalIndexs, this.compartment.cytosolIndexs)');
            this.reactionComplexCatalysisMatrix = max(0, -knowledgeBase.reactionProteinComplexStoichiometryMatrix(:, reactionGlobalIndexs, this.compartment.cytosolIndexs)');
            
            this.reactionThresholds = ...
                sum(this.reactionMonomerCatalysisMatrix, 2) + ...
                sum(this.reactionComplexCatalysisMatrix, 2);
        end
    end
    
    methods
        %allocate memory
        function allocateMemory(this, numTimePoints)
            import edu.stanford.covert.util.CircularSparseMat;
            
            this.polymerizedRegions       = CircularSparseMat([], [], [this.sequenceLen, this.nCompartments, numTimePoints], 1);
            this.linkingNumbers           = CircularSparseMat([], [], [this.sequenceLen, this.nCompartments, numTimePoints], 1);

            this.monomerBoundSites        = CircularSparseMat([], [], [this.sequenceLen, this.nCompartments, numTimePoints], 1);
            this.complexBoundSites        = CircularSparseMat([], [], [this.sequenceLen, this.nCompartments, numTimePoints], 1);
            
            this.gapSites                 = CircularSparseMat([], false(0,1), [this.sequenceLen, this.nCompartments, numTimePoints], 1);
            this.abasicSites              = CircularSparseMat([], false(0,1), [this.sequenceLen, this.nCompartments, numTimePoints], 1);
            this.damagedSugarPhosphates   = CircularSparseMat([], [],         [this.sequenceLen, this.nCompartments, numTimePoints], 1);
            this.damagedBases             = CircularSparseMat([], [],         [this.sequenceLen, this.nCompartments, numTimePoints], 1);
            this.intrastrandCrossLinks    = CircularSparseMat([], [],         [this.sequenceLen, this.nCompartments, numTimePoints], 1);
            this.strandBreaks             = CircularSparseMat([], false(0,1), [this.sequenceLen, this.nCompartments, numTimePoints], 1);
            this.hollidayJunctions        = CircularSparseMat([], false(0,1), [this.sequenceLen, this.nCompartments, numTimePoints], 1);
            
            this.segregated               = false(1, 1, numTimePoints);
            
            this.invalidate();
        end
    end
    
    methods
        function initialize(this)
            this.allocateMemory(1);
            
            this.polymerizedRegions(1, this.strandIndexs_ch1) = this.sequenceLen;
            this.linkingNumbers(1, this.strandIndexs_ch1) = this.sequenceLen / this.relaxedBasesPerTurn * (1 + this.equilibriumSuperhelicalDensity);
            
            this.invalidate();
        end
    end

    %general methods which query, but do not modify state
    methods
        %- positionsStrands must be a 2 column vector
        %- lengths must have same number of rows as positionsStrands
        function [tfs, idxs, positionsStrands, lengths] = sampleAccessibleRegions(this, ...
                nSites, weights, positionsStrands, lengths, ...
                bindingMonomers, bindingComplexs, ...
                isPositionsStrandFootprintCentroid, ignoreDamageFilter, returnOverlappingRegions, returnExtentAccessible, ...
                checkRegionSupercoiled)

            if nargin < 12
                checkRegionSupercoiled = false;
            end
            
            if isempty(weights)
                weights = ones(size(positionsStrands, 1), 1);
            end
            
            [footprint, footprint3Prime, footprint5Prime, footprintBindingStrandedness] = this.getDNAFootprint(bindingMonomers, bindingComplexs);
            
            idxs = [];
            while any(weights) && numel(idxs) < nSites
                %sample sites
                nMoreSites = min(max(2 * (nSites - numel(idxs)), 10), nnz(weights));
                selectedSites = this.randStream.randsample(numel(weights), nMoreSites, false, weights);
                weights(selectedSites) = 0;
                
                %determine if selected sites are accessible
                [tmpTfs, ~, ~, extents] = this.isRegionAccessible(positionsStrands(selectedSites, :), lengths(selectedSites), ...
                    bindingMonomers, bindingComplexs, isPositionsStrandFootprintCentroid, ignoreDamageFilter, ...
                    returnExtentAccessible, true, checkRegionSupercoiled);
                newIdxs = selectedSites(tmpTfs);
                lengths(selectedSites(~tmpTfs)) = 0;
                lengths(selectedSites(tmpTfs)) = extents;
                
                %if proteins binding, make sure proteins won't sterically overlap
                if ~returnOverlappingRegions
                    [idxs, newIdxs] = this.excludeOverlappingRegions(...
                        idxs, newIdxs, positionsStrands, lengths, ...
                        footprint, footprint3Prime, footprint5Prime, isPositionsStrandFootprintCentroid, ...
                        footprintBindingStrandedness == this.dnaStrandedness_dsDNA);
                end
                
                if numel(newIdxs) > nSites - numel(idxs)
                    newIdxs = newIdxs(1:nSites - numel(idxs));
                end
                idxs = [idxs; newIdxs];
            end
            
            %format output
            idxs = sort(idxs);
            tfs = false(size(positionsStrands, 1), 1);
            tfs(idxs) = true;
            positionsStrands = positionsStrands(idxs, :);
            lengths = lengths(idxs, :);
        end
        
        %randomly select among accessible sites (which if sequence seq is
        %specified, have this sequence) with probability probOrNSites (or if
        %probOrNSites > 1, randomly select probOrNSites sites)
        function positionsStrands = sampleAccessibleSites(this, prob, nSites, seq)
            %for convenience
            dnaLength = this.sequenceLen;
            posStrnds = find(this.polymerizedRegions);
            nStrands = max(posStrnds(:, 2));
            seqLen = numel(seq);
            
            %estimate total number of accessible sites
            [~, boundMonomers] = find(this.monomerBoundSites);
            [~, boundComplexs] = find(this.complexBoundSites);
            
            nAccessibleSites = ...
                collapse(this.polymerizedRegions) ...
                - sum(this.monomerDNAFootprints(boundMonomers, 1)) ...
                - sum(this.complexDNAFootprints(boundComplexs, 1)) ...
                - nnz(this.damagedSites);
            
            %calculate number of sites to select
            if ~isempty(prob)
                nGC = sum(seq == 'G' | seq == 'C');
                nSites = min(nSites, this.randStream.stochasticRound(...
                    nAccessibleSites * prob * ...
                    (this.sequenceGCContent/2)^nGC * ((1-this.sequenceGCContent)/2)^(seqLen-nGC)));
            end
            if nSites == 0
                positionsStrands = zeros(0, 2);
                return;
            end
            
            %randomly select nSites undamaged sites with sequence equal to seq
            positionsStrands = zeros(0, 2);
            iter = 0;
            while iter < 15
                %iterations
                iter = iter + 1;
                
                %more sites
                nMoreSites = nSites - size(positionsStrands, 1);
                nMoreSites = max(2 * nMoreSites, nMoreSites + 10);
                
                %pick random positions and strand
                positions = ceil(dnaLength * this.randStream.rand(nMoreSites, 1));
                strands   = ceil(nStrands * this.randStream.rand(nMoreSites, 1));
                
                %throw away positions that aren't inaccessible
                if isempty(seq)
                    [~, idxs] = this.isRegionAccessible([positions strands], 1, [], [], true, [], false, false);
                else
                    dir = 2 * (mod(strands, 2) == 1) - 1;
                    pos = 0:seqLen - 1;
                    
                    subsequences = this.sequence.subsequence(...
                        positions(:, ones(1, seqLen)) + ...
                        dir(:, ones(1, seqLen)) .* pos(ones(nMoreSites, 1), :), ...
                        strands);
                    
                    if size(seq, 1) == 1
                        if isscalar(seq)
                            idx = find(subsequences == seq);
                        elseif size(seq, 2) == 2
                            idx = find(subsequences(:, 1) == seq(1) & subsequences(:, 2) == seq(2));
                        else
                            idx = find(all(subsequences == seq(ones(size(subsequences, 1), 1), :), 2));
                        end
                    else
                        idx = find(ismember(subsequences, seq, 'rows'));
                    end
                    positions = positions(idx, :);
                    strands = strands(idx, :);
                    
                    [~, idxs] = this.isRegionAccessible([positions strands], seqLen, [], [], true, [], false, false);
                end
                
                %append to list of valid positions and strands
                positionsStrands = edu.stanford.covert.util.SparseMat.unique_subs(...
                    [positionsStrands; positions(idxs) strands(idxs)], [dnaLength this.nCompartments]);
                
                %if more sites that requested, randomly select
                if size(positionsStrands, 1) >= nSites
                    positionsStrands = this.randStream.randomlySelectNRows(positionsStrands, nSites);
                    break;
                end
            end
            
            %sort
            if nSites > 1
                positionsStrands = edu.stanford.covert.util.SparseMat.sort_subs(positionsStrands, [dnaLength this.nCompartments]);
            end
        end
        
        %Filters a list of sites (positionsStrands -- strands and positions along
        %strands), returning only those sites which are accessible to the
        %queried protein monomers and complexes (that is site which have been
        %polymerized, which aren't damaged, and which are either not bound by
        %protein, or are bound by proteins which at least one of the query
        %protein can release from the chromosome). Query proteins are indicated
        %by their indices within simulation.matureMatureIndexs
        %(monomerIndexs) and simulation.matureIndexs (complexIndexs).
        %
        %If monomerIndexs and complexIndexs are null, or aren't defined, this
        %function only returns sites that have been polymerized, aren't damaged,
        %and aren't bound by any proteins.
        %
        %Also returns the indices (idxs) of the query sites which are
        %accessible, and a boolean vector (tfs) indicated whether or not each
        %query site is accessible
        %
        %Returns true/false if regions defined by positionsStrands and lengths
        %are accessible/inaccessible.
        %
        %true  ==> region is accessible
        %false ==> region is not accessible
        function [tfs, idxs, positionsStrands, lengths] = isRegionAccessible(this, ...
                positionsStrands, lengths, ...
                bindingMonomers, bindingComplexs, ...
                isPositionsStrandFootprintCentroid, ignoreDamageFilter, ...
                returnExtent, returnOverlappingRegions, checkRegionSupercoiled)
            if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                return;
            end
            
            %% default values of inputs
            if ~all(lengths)
                throw(MException('Chromosome:invalidInput', 'Lengths must be non-zero integers'));
            elseif numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            if isempty(ignoreDamageFilter)
                ignoreDamageFilter = zeros(size(positionsStrands, 1), 0);
            elseif size(ignoreDamageFilter, 1) <= 1
                ignoreDamageFilter = ignoreDamageFilter(ones(size(positionsStrands, 1), 1), :);
            end
            if nargin < 9
                returnOverlappingRegions = true;
            end
            if nargin < 10
                checkRegionSupercoiled = false;
            end
            
            %% indices of complexs and monomers which at least one of the the
            %query proteins can release, and therefore shouldn't be cause for a
            %site to be filtered out if they are bound to a site
            [releasableMonomerIndexs, releasableComplexIndexs] = this.getReleasableProteins(bindingMonomers, bindingComplexs);
            
            %% DNA footprint of binding proteins,
            [footprint, footprint3Prime, footprint5Prime, footprintBindingStrandedness, footprintRegionStrandedness] = this.getDNAFootprint(bindingMonomers, bindingComplexs);
            
            %adding query positions if input positions are the centroid of
            %the binding proteins' DNA footpring
            queryPositionsStrands = positionsStrands;
            if isPositionsStrandFootprintCentroid
                strnd = queryPositionsStrands(:, 2);
                
                queryPositionsStrands(mod(strnd, 2) == 1 & lengths > 0, 1) = queryPositionsStrands(mod(strnd, 2) == 1 & lengths > 0, 1) - footprint5Prime;
                queryPositionsStrands(mod(strnd, 2) == 0 & lengths > 0, 1) = queryPositionsStrands(mod(strnd, 2) == 0 & lengths > 0, 1) - footprint3Prime;
                queryPositionsStrands(mod(strnd, 2) == 1 & lengths < 0, 1) = queryPositionsStrands(mod(strnd, 2) == 1 & lengths < 0, 1) + footprint3Prime;
                queryPositionsStrands(mod(strnd, 2) == 0 & lengths < 0, 1) = queryPositionsStrands(mod(strnd, 2) == 0 & lengths < 0, 1) + footprint5Prime;
                
                ignoreDamageFilter(mod(strnd, 2) == 1 & lengths > 0, :) = ignoreDamageFilter(mod(strnd, 2) == 1 & lengths > 0, :) + footprint5Prime;
                ignoreDamageFilter(mod(strnd, 2) == 0 & lengths > 0, :) = ignoreDamageFilter(mod(strnd, 2) == 0 & lengths > 0, :) + footprint3Prime;
                ignoreDamageFilter(mod(strnd, 2) == 1 & lengths < 0, :) = ignoreDamageFilter(mod(strnd, 2) == 1 & lengths < 0, :) - footprint3Prime;
                ignoreDamageFilter(mod(strnd, 2) == 0 & lengths < 0, :) = ignoreDamageFilter(mod(strnd, 2) == 0 & lengths < 0, :) - footprint5Prime;
            elseif lengths < 0
                queryPositionsStrands(:, 1) = queryPositionsStrands(:, 1) + footprint - 1;
            end
            queryLengths = lengths + sign(lengths) * (footprint - 1);
            
            %% filter query sites
            switch footprintRegionStrandedness
                case this.dnaStrandedness_dsDNA
                    [~, ~, ~, polymerized] = this.isRegionDoubleStranded(queryPositionsStrands, queryLengths, true, checkRegionSupercoiled);
                case this.dnaStrandedness_ssDNA
                    [~, ~, ~, polymerized] = this.isRegionSingleStranded(queryPositionsStrands, queryLengths, true);                
                case this.dnaStrandedness_xsDNA
                    %TODO: calculate extent
                    polymerized = queryLengths;
                otherwise
                    throw(MException('Chromosome:invalidInput','Invalid footprint strandedness'));
            end
            
            [~, ~, ~, proteinFree] = this.isRegionProteinFree(queryPositionsStrands, queryLengths, true, ...
                releasableMonomerIndexs, releasableComplexIndexs, ...
                footprintBindingStrandedness == this.dnaStrandedness_dsDNA, ...
                footprintRegionStrandedness == this.dnaStrandedness_xsDNA);
            [~, ~, ~, undamaged] = this.isRegionUndamaged(queryPositionsStrands, queryLengths, ...
                footprintBindingStrandedness == this.dnaStrandedness_dsDNA, ignoreDamageFilter, true, ...
                footprintRegionStrandedness == this.dnaStrandedness_xsDNA);
            
            extents = sign(lengths) .* max(0, (min(abs([proteinFree polymerized undamaged]), [], 2) - (footprint - 1)));
            
            %% Make sure regions don't sterically overlap
            if ~returnOverlappingRegions
                if returnExtent
                    tfs = find(extents ~= 0);
                else
                    tfs = find(lengths == extents);
                end
                [~, idxs] = this.excludeOverlappingRegions([], tfs, positionsStrands, extents, ...
                    footprint, footprint3Prime, footprint5Prime, isPositionsStrandFootprintCentroid, ...
                    footprintBindingStrandedness == this.dnaStrandedness_dsDNA);
                extents(setdiff(1:end, idxs)) = 0;
            end
            
            %% format output
            %- extract sites which pass filter
            %- reshape indices of sites which pass filter
            %- construct boolean indicating which sites pass filter
            if returnExtent
                tfs = true(size(positionsStrands, 1), 1);
            else
                tfs = lengths == extents;
            end
            idxs = find(tfs);
            idxs = idxs(:);
            positionsStrands = positionsStrands(idxs, :);
            lengths = extents(idxs, :);
        end
    end
    
    %private methods for query state by region
    methods
        function [tfs, idxs, positionsStrands, lengths] = isRegionSingleStranded(this, positionsStrands, lengths, returnExtent)
            if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                return;
            end
            if numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            
            [polPosStrnds, polLens] = find(this.polymerizedRegions);
            
            %only 1 chromosome
            if size(positionsStrands, 2) == 2 && all(polLens == this.sequenceLen) && mod(numel(polLens), 2) == 0 && all(polPosStrnds(:, 2) == (1:numel(polLens))')
                if returnExtent
                    tfs = true(size(positionsStrands, 1), 1);
                    idxs = (1:numel(tfs))';
                    lengths(:) = 0;
                else
                    tfs = false(size(lengths));
                    idxs = zeros(0, 1);
                    positionsStrands = zeros(0, 2);
                    lengths = zeros(0, 1);
                end
                return;
            end
            
            %all other cases
            [tfs, idxs, positionsStrands, lengths] = this.isRegionPolymerized(positionsStrands, lengths, returnExtent, this.dnaStrandedness_ssDNA);
        end
        
        function [tfs, idxs, positionsStrands, lengths] = isRegionDoubleStranded(this, positionsStrands, lengths, returnExtent, checkRegionSupercoiled)
            if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                return;
            end
            if nargin < 5
                checkRegionSupercoiled = false;
            end
            
            if numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            
            [polPosStrnds, polLens] = find(this.polymerizedRegions);
            if checkRegionSupercoiled
                supPosStrnds = find(this.supercoiled);
            end
            
            %only fully formed chromosomes
            if size(positionsStrands, 2) == 2 && all(polLens == this.sequenceLen) && mod(numel(polLens), 2) == 0 && all(polPosStrnds(:, 2) == (1:numel(polLens))')
                if numel(polLens) == 2
                    if returnExtent
                        tfs = true(size(positionsStrands, 1), 1);
                        idxs = (1:numel(tfs))';
                        lengths(positionsStrands(:, 2) > 2) = 0;
                        if checkRegionSupercoiled
                            lengths(~ismembc(positionsStrands(:, 2), supPosStrnds(:, 2))) = 0;
                        end
                    else
                        tfs = positionsStrands(:, 2) <= 2;
                        if checkRegionSupercoiled
                            tfs = tfs & ismembc(positionsStrands(:, 2), supPosStrnds(:, 2));
                        end
                        idxs = find(tfs);
                        positionsStrands = positionsStrands(tfs, :);
                        lengths = lengths(tfs, 1);
                    end
                else
                    if returnExtent
                        tfs = true(size(positionsStrands, 1), 1);
                        idxs = (1:numel(tfs))';
                        if checkRegionSupercoiled
                            lengths(~ismembc(positionsStrands(:, 2), supPosStrnds(:, 2))) = 0;
                        end
                    else
                        tfs = true(size(positionsStrands, 1), 1);
                        if checkRegionSupercoiled
                            tfs = tfs & ismembc(positionsStrands(:, 2), supPosStrnds(:, 2));
                        end
                        idxs = find(tfs);
                        positionsStrands = positionsStrands(tfs, :);
                        lengths = lengths(tfs, 1);
                    end
                end
                return;
            end
            
            %all other cases
            [tfs, idxs, positionsStrands, lengths] = this.isRegionPolymerized(positionsStrands, lengths, returnExtent, this.dnaStrandedness_dsDNA, checkRegionSupercoiled);
        end
        
        function [tfs, idxs, positionsStrands, lengths] = isRegionPolymerized(this, ...
                positionsStrands, lengths, returnExtent, checkRegionStrandedness, checkRegionSupercoiled, polymerizedRegions)
            if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                return;
            end
            if ~all(lengths)
                throw(MException('Chromosome:invalidInput','Lengths must be non-zero integers'));
            elseif numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            if nargin < 7
                polymerizedRegions = this.polymerizedRegions;
            end
            
            if nargin >= 5 && ~isempty(checkRegionStrandedness) && checkRegionStrandedness == this.dnaStrandedness_dsDNA
                [regions, regionLengths] = find(this.doubleStrandedRegions);
                if nargin >= 6 && ~isempty(checkRegionSupercoiled) && checkRegionSupercoiled
                    tfs = this.supercoiled(regions);
                    regions = regions(tfs, :);
                    regionLengths = regionLengths(tfs, :);
                end
            elseif nargin >= 5 && ~isempty(checkRegionStrandedness) && checkRegionStrandedness == this.dnaStrandedness_ssDNA
                [regions, regionLengths] = find(this.singleStrandedRegions);
            else
                [regions, regionLengths] = find(polymerizedRegions);
            end
            
            %all strands either fully, or not at all synthesized
            if all(regionLengths == this.sequenceLen) && size(regions, 2) == 2
                if size(regions, 1) == this.nCompartments
                    %all chromosomes fully synethesized
                    tfs = true(size(positionsStrands,1), 1);
                    idxs = (1:size(positionsStrands,1))';
                elseif isequal(regions(:,2), [1; 2])
                    %first chromosome fully synethesized
                    if returnExtent
                        tfs = true(size(positionsStrands,1), 1);
                        idxs = (1:size(positionsStrands,1))';
                        lengths(~ismembc(positionsStrands(:,2), regions(:, 2))) = 0;
                    else
                        tfs = ceil(positionsStrands(:,2)/2) == 1;
                        idxs = reshape(find(tfs), [], 1);
                        positionsStrands = positionsStrands(idxs, :);
                        lengths = lengths(idxs, :);
                    end
                else
                    %all strands either fully, or not at all synthesized
                    if returnExtent
                        tfs = true(size(positionsStrands,1), 1);
                        idxs = (1:size(positionsStrands,1))';
                        lengths(~ismembc(positionsStrands(:,2), regions(:, 2))) = 0;
                    else
                        tfs = ismembc(positionsStrands(:,2), regions(:, 2));
                        idxs = reshape(find(tfs), [], 1);
                        positionsStrands = positionsStrands(idxs, :);
                        lengths = lengths(idxs, :);
                    end
                end
                return;
            end
            
            regionStartCoors = regions(:, 1);
            regionEndCoors = regionStartCoors  + regionLengths - 1;
            regionStrandTimes = [regions(:, 2:end) ones(size(regions, 1), size(positionsStrands,2) - size(regions, 2))];
            
            startCoors = mod(positionsStrands(:, 1) + min(0, lengths+1) - 1, this.sequenceLen) +1;
            endCoors = startCoors + abs(lengths) - 1;
            strandTimes = [positionsStrands(:, 2:end) ones(size(positionsStrands, 1), size(regions, 2) - size(positionsStrands,2))];
            idxs = (1:size(positionsStrands, 1))';
            
            tmpIdxs = find(endCoors > this.sequenceLen);
            strandTimes = [strandTimes; strandTimes(tmpIdxs, :)];
            startCoors = [startCoors; max(1, startCoors(tmpIdxs) - this.sequenceLen)];
            endCoors = [endCoors; endCoors(tmpIdxs) - this.sequenceLen];
            endCoors(tmpIdxs) = this.sequenceLen;
            idxs = [idxs; idxs(tmpIdxs)];
            
            tmpIdxs = find(startCoors < 0);
            strandTimes = [strandTimes; strandTimes(tmpIdxs, :)];
            startCoors = [startCoors; startCoors(tmpIdxs) + this.sequenceLen];
            endCoors = [endCoors; min(this.sequenceLen, endCoors(tmpIdxs) + this.sequenceLen)];
            startCoors(tmpIdxs) = 1;
            idxs = [idxs; idxs(tmpIdxs)];
            
            extents = zeros(size(positionsStrands, 1), 1);
            if ndims(polymerizedRegions) == 2
                regionStrandTimeSiz = size(polymerizedRegions, 2);
            else
                regionStrandTimeSiz = [size(polymerizedRegions, 2) size(polymerizedRegions, 3)];
            end
            
            tfs = ...
                ((max(endCoors) >= regionStartCoors & max(endCoors) <= regionEndCoors) | ...
                (min(startCoors) >= regionStartCoors & min(startCoors) <= regionEndCoors) | ...
                (min(startCoors) < regionStartCoors & max(endCoors) > regionEndCoors)) & ...
                edu.stanford.covert.util.SparseMat.ismember_subs(regionStrandTimes, strandTimes, regionStrandTimeSiz);
            regionStartCoors = regionStartCoors(tfs);
            regionEndCoors = regionEndCoors(tfs);
            regionStrandTimes = regionStrandTimes(tfs, :);
            
            if size(positionsStrands, 1) == numel(idxs)
                [~, idxs1, idxs2] = edu.stanford.covert.util.SparseMat.unique_subs([startCoors endCoors], [this.sequenceLen; this.sequenceLen]);
                for j = 1:numel(idxs1)
                    idxs3 = find(idxs2 == j);
                    startCoor = startCoors(idxs3(1));
                    endCoor = endCoors(idxs3(1));
                    len = lengths(idxs(idxs3(1)));
                    if len > 0
                        idxs3 = idxs3(idxs3 <= size(positionsStrands, 1) | extents(idxs(idxs3)) >= this.sequenceLen - startCoor + 1);
                        tfs = startCoor >= regionStartCoors & startCoor <= regionEndCoors;
                        tmpregionEndCoors = regionEndCoors(tfs);
                        [tfs, idxs4] = edu.stanford.covert.util.SparseMat.ismember_subs(strandTimes(idxs3, :), regionStrandTimes(tfs, :), regionStrandTimeSiz);
                        extents(idxs(idxs3(tfs))) = extents(idxs(idxs3(tfs))) + min(tmpregionEndCoors(idxs4(tfs)), endCoor) - startCoor + 1;
                    else
                        idxs3 = idxs3(idxs3 <= size(positionsStrands, 1) | -extents(idxs(idxs3)) >= endCoor);
                        tfs = endCoor >= regionStartCoors & endCoor <= regionEndCoors;
                        tmpregionStartCoors = regionStartCoors(tfs);
                        [tfs, idxs4] = edu.stanford.covert.util.SparseMat.ismember_subs(strandTimes(idxs3, :), regionStrandTimes(tfs, :), regionStrandTimeSiz);
                        extents(idxs(idxs3(tfs))) = extents(idxs(idxs3(tfs))) - (endCoor - max(tmpregionStartCoors(idxs4(tfs)), startCoor) + 1);
                    end
                end
            else
                for i = 1:numel(startCoors)
                    if lengths(idxs(i)) > 0
                        if i > size(positionsStrands, 1) && extents(idxs(i)) < this.sequenceLen - startCoors(idxs(i)) + 1
                            continue;
                        end
                        
                        tmpIdxs = find(...
                            startCoors(i) >= regionStartCoors & ...
                            startCoors(i) <= regionEndCoors & ...
                            ismember(regionStrandTimes, strandTimes(i, :), 'rows'), 1, 'first');
                        
                        if ~isempty(tmpIdxs)
                            extents(idxs(i)) = extents(idxs(i)) + min(regionEndCoors(tmpIdxs), endCoors(i)) - startCoors(i) + 1;
                        end
                    else
                        if i > size(positionsStrands, 1) && -extents(idxs(i)) < endCoors(idxs(i))
                            continue;
                        end
                        
                        tmpIdxs = find(...
                            endCoors(i) >= regionStartCoors & ...
                            endCoors(i) <= regionEndCoors & ...
                            ismember(regionStrandTimes, strandTimes(i, :), 'rows'), 1, 'first');
                        
                        if ~isempty(tmpIdxs)
                            extents(idxs(i)) = extents(idxs(i)) - (endCoors(i) - max(regionStartCoors(tmpIdxs), startCoors(i)) + 1);
                        end
                    end
                end
            end
            
            if returnExtent
                tfs = true(size(positionsStrands, 1), 1);
            else
                tfs = extents == lengths;
            end
            
            idxs = reshape(find(tfs), [], 1);
            positionsStrands = positionsStrands(idxs, :);
            lengths = extents(idxs, :);
        end
        
        function [tfs, idxs, positionsStrands, lengths] = isRegionNotPolymerized(this, positionsStrands, lengths, returnExtent)
            if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                return;
            end
            
            if ~all(lengths)
                throw(MException('Chromosome:invalidInput','Lengths must be non-zero integers'));
            elseif numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            
            unpolRgns = this.unpolymerizedRegions;
            if nnz(unpolRgns) == 0
                if returnExtent
                    tfs = true(size(positionsStrands,1), 1);
                    idxs = (1:size(positionsStrands,1))';
                    lengths =  zeros(size(positionsStrands,1), 1);
                else
                    tfs = false(size(positionsStrands,1), 1);
                    idxs = zeros(0, 1);
                    positionsStrands = zeros(0, ndims(positionsStrands));
                    lengths = zeros(0, 1);
                end
                return;
            end
            
            [regions, regionLengths] = find(unpolRgns);
            regionStartCoors = regions(:, 1);
            regionEndCoors = regionStartCoors  + regionLengths - 1;
            regionStrandTimes = [regions(: ,2:end) ones(size(regions, 1), size(positionsStrands,2) - size(regions, 2))];
            
            startCoors = mod(positionsStrands(:, 1) + min(0, lengths+1) - 1, this.sequenceLen) +1;
            endCoors = startCoors + abs(lengths) - 1;
            strandTimes = [positionsStrands(: ,2:end) ones(size(positionsStrands, 1), size(regions, 2) - size(positionsStrands,2))];
            idxs = (1:size(positionsStrands, 1))';
            
            tmpIdxs = find(endCoors > this.sequenceLen);
            strandTimes = [strandTimes; strandTimes(tmpIdxs, :)];
            startCoors = [startCoors; max(1, startCoors(tmpIdxs) - this.sequenceLen)];
            endCoors = [endCoors; endCoors(tmpIdxs) - this.sequenceLen];
            endCoors(tmpIdxs) = this.sequenceLen;
            idxs = [idxs; idxs(tmpIdxs)];
            
            tmpIdxs = find(startCoors < 0);
            strandTimes = [strandTimes; strandTimes(tmpIdxs, :)];
            startCoors = [startCoors; startCoors(tmpIdxs) + this.sequenceLen];
            endCoors = [endCoors; min(this.sequenceLen, endCoors(tmpIdxs) + this.sequenceLen)];
            startCoors(tmpIdxs) = 1;
            idxs = [idxs; idxs(tmpIdxs)];
            
            extents = zeros(size(positionsStrands, 1), 1);
            for i = 1:numel(startCoors)
                if lengths(idxs(i)) > 0
                    if i > size(positionsStrands, 1) && extents(idxs(i)) < this.sequenceLen - startCoors(idxs(i)) + 1
                        continue;
                    end
                    
                    tmpIdxs = find(...
                        startCoors(i) >= regionStartCoors & ...
                        startCoors(i) <= regionEndCoors & ...
                        ismember(regionStrandTimes, strandTimes(i, :), 'rows'), 1, 'first');
                    
                    if ~isempty(tmpIdxs)
                        extents(idxs(i)) = extents(idxs(i)) + min(regionEndCoors(tmpIdxs), endCoors(i)) - startCoors(i) + 1;
                    end
                else
                    if i > size(positionsStrands, 1) && -extents(idxs(i)) < endCoors(idxs(i))
                        continue;
                    end
                    
                    tmpIdxs = find(...
                        endCoors(i) >= regionStartCoors & ...
                        endCoors(i) <= regionEndCoors & ...
                        ismember(regionStrandTimes, strandTimes(i, :), 'rows'), 1, 'first');
                    
                    if ~isempty(tmpIdxs)
                        extents(idxs(i)) = extents(idxs(i)) - (endCoors(i) - max(regionStartCoors(tmpIdxs), startCoors(i)) + 1);
                    end
                end
            end
            
            if returnExtent
                tfs = true(size(positionsStrands, 1), 1);
            else
                tfs = extents == lengths;
            end
            
            idxs = reshape(find(tfs), [], 1);
            positionsStrands = positionsStrands(idxs, :);
            lengths = extents(idxs, :);
        end

        %lengths must be non-negative integer
        %releasableMonomerIndexs and releasableComplexIndexs must be sorted
        function [tfs, idxs, positionsStrands, lengths, monomers, complexs] = isRegionProteinFree(this, ...
                positionsStrands, lengths, returnExtent, releasableMonomerIndexs, releasableComplexIndexs, bothStrands, allStrands)
            if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                monomers = zeros(0, 0);
                complexs = zeros(0, 0);
                return;
            end
            
            if numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            if nargin < 8
                allStrands = false;
            end
            
            startCoors = mod(positionsStrands(:, 1) + min(0, lengths+1) - 1, this.sequenceLen) + 1;
            endCoors = startCoors + abs(lengths) - 1;
            strnds = positionsStrands(:, 2);
            if bothStrands
                strnds = ceil(strnds/2);
            end
            if allStrands
                strnds(:) = 1;
            end
            idxs = (1:size(positionsStrands, 1))';
            
            tmpIdxs = find(endCoors > this.sequenceLen);
            startCoors = [startCoors; max(1, startCoors(tmpIdxs) - this.sequenceLen)];
            endCoors = [endCoors; endCoors(tmpIdxs) - this.sequenceLen];
            endCoors(tmpIdxs) = this.sequenceLen;
            strnds = [strnds; strnds(tmpIdxs, :)];
            idxs = [idxs; idxs(tmpIdxs)];
            
            tmpIdxs = find(startCoors < 0);
            startCoors = [startCoors; startCoors(tmpIdxs) + this.sequenceLen];
            endCoors = [endCoors; min(this.sequenceLen, endCoors(tmpIdxs) + this.sequenceLen)];
            startCoors(tmpIdxs) = 1;
            strnds = [strnds; strnds(tmpIdxs, :)];
            idxs = [idxs; idxs(tmpIdxs)];
            
            %initialize output
            extents  = abs(lengths);
            monomers = sparse(size(positionsStrands, 1), max(abs([lengths; 0])));
            complexs = sparse(size(positionsStrands, 1), max(abs([lengths; 0])));
            
            %bound monomers
            [subs, vals] = find(this.monomerBoundSites);
            if ~bothStrands
                dsTfs = this.monomerDNAFootprintBindingStrandedness(vals) == this.dnaStrandedness_dsDNA;
                subs = [
                    subs(~dsTfs, :)
                    subs( dsTfs, 1) 2*ceil(subs(dsTfs, 2)/2)-1
                    subs( dsTfs, 1) 2*ceil(subs(dsTfs, 2)/2)
                    ];
                vals = [
                    vals(~dsTfs)
                    vals( dsTfs)
                    vals( dsTfs)
                    ];
                [subs, order] = edu.stanford.covert.util.SparseMat.sort_subs(subs, [this.sequenceLen this.nCompartments]);
                vals = vals(order, :);
            end
            if allStrands
                tfs = ~ismembc(vals, this.monomerIndexs_replisome);
                subs = subs(tfs, :);
                vals = vals(tfs);
            end
            monomerStarts = subs(:, 1);
            monomerEnds   = monomerStarts + this.monomerDNAFootprints(vals, :) - 1;
            monomerStrnds = subs(:, 2);
            if bothStrands
                monomerStrnds = ceil(monomerStrnds/2);
            end
            if allStrands
                monomerStrnds(:) = 1;
            end
            
            tmpIdxs = find(monomerEnds > this.sequenceLen);
            monomerStarts = [monomerStarts; ones(numel(tmpIdxs), 1)];
            monomerEnds   = [monomerEnds;   monomerEnds(tmpIdxs, :) - this.sequenceLen];
            monomerStrnds = [monomerStrnds; monomerStrnds(tmpIdxs, :)];
            vals = [vals; vals(tmpIdxs, :)];
            
            for i = 1:numel(startCoors)
                monMask = ...
                    (monomerStarts >= startCoors(i) & monomerStarts <= endCoors(i)) | ...
                    (monomerEnds >= startCoors(i) & monomerEnds <= endCoors(i)) | ...
                    (monomerStarts <= startCoors(i) & monomerEnds >= endCoors(i));
                tmpIdxs = find(monMask);
                if isempty(tmpIdxs); continue; end;
                
                tmpIdxs = tmpIdxs(monomerStrnds(monMask, 1) ==  strnds(i, 1));
                if isempty(tmpIdxs); continue; end;
                
                if lengths(idxs(i)) > 0
                    if i <= size(positionsStrands, 1)
                        monomerPos = monomerStarts(tmpIdxs) - startCoors(i) + 1;
                    else
                        monomerPos = monomerStarts(tmpIdxs) - startCoors(i) + 1 + (this.sequenceLen-positionsStrands(idxs(i), 1));
                    end
                else
                    if i <= size(positionsStrands, 1)
                        monomerPos = endCoors(i) - monomerEnds(tmpIdxs) + 1;
                    else
                        monomerPos = endCoors(i) - monomerEnds(tmpIdxs) + 1 + positionsStrands(idxs(i), 1);
                    end
                end
                monomerPos = max(monomerPos, 1);
                
                extents(idxs(i)) = min([
                    extents(idxs(i));
                    monomerPos(~ismembc(vals(tmpIdxs), releasableMonomerIndexs))-1]);
                
                monomers(idxs(i), monomerPos) = vals(tmpIdxs); %#ok<SPRIX>
            end
            
            %bound complexes
            [subs, vals] = find(this.complexBoundSites);
            if ~bothStrands
                dsTfs = this.complexDNAFootprintBindingStrandedness(vals) == this.dnaStrandedness_dsDNA;
                subs = [
                    subs(~dsTfs, :)
                    subs( dsTfs, 1) 2*ceil(subs(dsTfs, 2)/2)-1
                    subs( dsTfs, 1) 2*ceil(subs(dsTfs, 2)/2)
                    ];
                vals = [
                    vals(~dsTfs)
                    vals( dsTfs)
                    vals( dsTfs)
                    ];
                [subs, order] = edu.stanford.covert.util.SparseMat.sort_subs(subs, [this.sequenceLen this.nCompartments]);
                vals = vals(order, :);
            end
            if allStrands
                tfs = ~ismembc(vals, this.complexIndexs_replisome);
                subs = subs(tfs, :);
                vals = vals(tfs);
            end
            complexStarts = subs(:, 1);
            complexEnds   = complexStarts + this.complexDNAFootprints(vals, :) - 1;
            complexStrnds = subs(:, 2);
            if bothStrands
                complexStrnds = ceil(complexStrnds/2);
            end
            if allStrands
                complexStrnds(:) = 1;
            end
            
            tmpIdxs = find(complexEnds > this.sequenceLen);
            complexStarts = [complexStarts; ones(numel(tmpIdxs), 1)];
            complexEnds   = [complexEnds;   complexEnds(tmpIdxs, :) - this.sequenceLen];
            complexStrnds = [complexStrnds; complexStrnds(tmpIdxs, :)];
            vals = [vals; vals(tmpIdxs, :)];
            
            for i = 1:numel(startCoors)
                comMask = ...
                    (complexStarts >= startCoors(i) & complexStarts <= endCoors(i)) | ...
                    (complexEnds >= startCoors(i) & complexEnds <= endCoors(i)) | ...
                    (complexStarts <= startCoors(i) & complexEnds >= endCoors(i));
                tmpIdxs = find(comMask);
                if isempty(tmpIdxs); continue; end;
                
                tmpIdxs = tmpIdxs(complexStrnds(comMask, 1) ==  strnds(i, 1));
                if isempty(tmpIdxs); continue; end;
                
                if lengths(idxs(i)) > 0
                    if i <= size(positionsStrands, 1)
                        complexPos = complexStarts(tmpIdxs) - startCoors(i) + 1;
                    else
                        complexPos = complexStarts(tmpIdxs) - startCoors(i) + 1 + (this.sequenceLen-positionsStrands(idxs(i), 1));
                    end
                else
                    if i <= size(positionsStrands, 1)
                        complexPos = endCoors(i) - complexEnds(tmpIdxs) + 1;
                    else
                        complexPos = endCoors(i) - complexEnds(tmpIdxs) + 1 + positionsStrands(idxs(i), 1);
                    end
                end
                complexPos = max(complexPos, 1);
                
                extents(idxs(i)) = min([
                    extents(idxs(i));
                    complexPos(~ismembc(vals(tmpIdxs), releasableComplexIndexs))-1]);
                
                complexs(idxs(i), complexPos) = vals(tmpIdxs); %#ok<SPRIX>
            end
            
            %format output
            extents = sign(lengths) .* extents;
            if returnExtent
                tfs = true(size(positionsStrands, 1), 1);
            else
                tfs = extents == lengths;
            end
            idxs = reshape(find(tfs), [], 1);
            positionsStrands = positionsStrands(idxs, :);
            lengths = extents(idxs, :);
        end
        
        %lengths must be non-negative integer
        function [tfs, idxs, positionsStrands, lengths] = isRegionUndamaged(this, ...
                positionsStrands, lengths, isEitherStrandDamaged, ignoreDamageFilter, returnExtent, isAnyStrandDamaged)
            if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                return;
            end
            
            if numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            if nargin < 7
                isAnyStrandDamaged = false;
            end
            
            [damagedPosStrnds, damages] = find(this.damagedSites);
            if isempty(damages)
                tfs = true(size(lengths));
                idxs = (1:numel(lengths))';
                return;
            end
            
            if isempty(ignoreDamageFilter)
                ignoreDamageFilter = zeros(size(positionsStrands, 1), 0);
            elseif size(ignoreDamageFilter, 1) <= 1
                ignoreDamageFilter = reshape(ignoreDamageFilter(ones(size(positionsStrands, 1), 1), :), size(positionsStrands, 1), []);
            end
                        
            damagedPositions   = damagedPosStrnds(:, 1);
            damagedStrandTimes = [damagedPosStrnds(:, 2:end) ones(size(damagedPosStrnds, 1), size(positionsStrands, 2) - size(damagedPosStrnds, 2))];
            strandTimes        = [positionsStrands(:, 2:end) ones(size(positionsStrands, 1), size(damagedPosStrnds, 2) - size(positionsStrands, 2))];
            
            if isAnyStrandDamaged
                strandTimes(:, 1) = 1;
                damagedStrandTimes(:, 1) = 1;
            elseif isEitherStrandDamaged
                strandTimes(:, 1) = ceil(strandTimes(:, 1) / 2);
                damagedStrandTimes(:, 1) = ceil(damagedStrandTimes(:, 1) / 2);
            end
                        
            startPositions = mod(positionsStrands(:, 1) + min(0, lengths+1) - 1, this.sequenceLen) + 1;
            endPositions   = startPositions + abs(lengths) - 1;
            idxs = (1:size(positionsStrands,1))';
            
            tmpIdxs = find(endPositions > this.sequenceLen);
            strandTimes = [strandTimes; strandTimes(tmpIdxs)];
            startPositions = [startPositions; max(1, startPositions(tmpIdxs) - this.sequenceLen)];
            endPositions = [endPositions; endPositions(tmpIdxs) - this.sequenceLen];
            endPositions(tmpIdxs) = this.sequenceLen;
            idxs = [idxs; idxs(tmpIdxs)];
            
            tmpIdxs = find(startPositions < 0);
            strandTimes = [strandTimes; strandTimes(tmpIdxs)];
            startPositions = [startPositions; startPositions(tmpIdxs) + this.sequenceLen];
            endPositions = [endPositions; min(this.sequenceLen, endPositions(tmpIdxs) + this.sequenceLen)];
            startPositions(tmpIdxs) = 1;
            idxs = [idxs; idxs(tmpIdxs)];
            
            extents = zeros(size(positionsStrands,1), 1);
            for i = 1:numel(startPositions)
                tmpIdxs = ...
                    damagedPositions >= startPositions(i, 1) & ...
                    damagedPositions <= endPositions(i, 1) & ...
                    ismember(damagedStrandTimes, strandTimes(i, :), 'rows');
                
                if lengths(idxs(i)) > 0
                    if i > size(positionsStrands, 1) && extents(idxs(i)) < this.sequenceLen - startPositions(idxs(i)) + 1
                        continue;
                    end
                    
                    if i <= size(positionsStrands, 1)
                        tmpIgnoreDamageFilter = ignoreDamageFilter(idxs(i), :);
                    else
                        tmpIgnoreDamageFilter = ignoreDamageFilter(idxs(i), :) - (this.sequenceLen - startPositions(idxs(i)) + 1);
                    end
                    
                    extent = min(setdiff(damagedPositions(tmpIdxs) - startPositions(i) + 1, tmpIgnoreDamageFilter)) - 1;
                    if isempty(extent)
                        extent = endPositions(i) - startPositions(i) + 1;
                    end
                    extents(idxs(i)) = extents(idxs(i)) + extent;
                else
                    if i <= size(positionsStrands, 1)
                        tmpIgnoreDamageFilter = ignoreDamageFilter(idxs(i), :);
                    else
                        tmpIgnoreDamageFilter = ignoreDamageFilter(idxs(i), :) + positionsStrands(idxs(i), 1);
                    end
                    
                    extent = min(setdiff(endPositions(i) - damagedPositions(tmpIdxs) + 1, tmpIgnoreDamageFilter - lengths(idxs(i)) - 1)) - 1;
                    if isempty(extent)
                        extents(idxs(i)) = extents(idxs(i)) -(endPositions(i) - startPositions(i) + 1);
                    else
                        extents(idxs(i)) = extents(idxs(i)) - extent;
                    end
                end
            end
            
            if returnExtent
                tfs = true(size(positionsStrands, 1), 1);
            else
                tfs = extents == lengths;
            end
            idxs = reshape(find(tfs), [], 1);
            positionsStrands = positionsStrands(idxs, :);
            lengths = extents(idxs, :);
        end
    end
    
    %additional private methods used to query state
    methods
        %Converts from the base-pair centric view of the chromosomes that this
        %class uses to a strand-centric view where each column represents a
        %single strand.
        function varargout = getStrandView(this, outputs)
            import edu.stanford.covert.util.CircularSparseMat;
            
            if nargin < 2
                outputs = {
                    'polymerizedRegions';
                    'monomerBoundSites';
                    'complexBoundSites';
                    'gapSites';
                    'abasicSites';
                    'damagedSugarPhosphates';
                    'damagedBases';
                    'intrastrandCrossLinks';
                    'strandBreaks';
                    'hollidayJunctions';
                    };
            end
            if ~iscell(outputs)
                outputs = {outputs};
            end
            if numel(outputs) < nargout
                throw(MException('Chromosome:tooManyOutputs', 'Too many outputs requested'));
            end
            
            varargout = cell(nargout, 1);
            for i = 1:nargout
                tmp = this.(outputs{i});
                switch outputs{i}
                    case 'polymerizedRegions'
                        [positionsStrands, lengths] = find(this.polymerizedRegions);
                        positionsStrands = [positionsStrands ones(size(positionsStrands,1), 3-size(positionsStrands,2))];
                        idx1 = find(positionsStrands(:,2) == this.strandIndexs_negative(1));
                        idx2 = find(positionsStrands(:,2) == this.strandIndexs_negative(2));
                        
                        starts1 = positionsStrands(idx1, 1);
                        starts2 = positionsStrands(idx2, 1);
                        ends1 = starts1 + lengths(idx1) - 1;
                        ends2 = starts2 + lengths(idx2) - 1;
                        times1 = positionsStrands(idx1, 3:end);
                        times2 = positionsStrands(idx2, 3:end);
                        
                        tmpStarts2 = zeros(0,1);
                        tmpEnds2 = zeros(0,1);
                        tmpTimes2 = zeros(0,size(times2,2));
                        for j = 1:numel(idx2)
                            idx = find(...
                                ((starts1 >= starts2(j) & starts1 <= ends2(j)) | ...
                                (ends1    >= starts2(j) & ends1   <= ends2(j)) | ...
                                (starts1  <= starts2(j) & ends1   >= ends2(j))) & ...
                                all(times1 == times2(j, :), 2));                                                                                   
                            
                            tmpStarts2 = [tmpStarts2;
                                max(starts1(idx), starts2(j))];
                            tmpEnds2 = [tmpEnds2;
                                min(ends1(idx), ends2(j))];
                            tmpTimes2 = [tmpTimes2;
                                times1(idx, :)];
                            
                            starts1 = [starts1(setdiff(1:end, idx));
                                min([starts1(idx); starts2(j)])];
                            ends1 = [ends1(setdiff(1:end, idx));
                                max([ends1(idx); ends2(j)])];
                            times1 = [times1(setdiff(1:end, idx), :);
                                times2(j)];
                        end
                        
                        tmpPosStrnds = [positionsStrands(setdiff(1:end, [idx1; idx2]), :);
                            starts1 this.strandIndexs_negative(ones(numel(starts1), 1)) times1;
                            tmpStarts2 this.strandIndexs_negative(2 * ones(numel(tmpStarts2), 1)) tmpTimes2];
                        tmpLengths = [lengths(setdiff(1:end, [idx1; idx2]), :);
                            ends1 - starts1 + 1;
                            tmpEnds2 - tmpStarts2 + 1;];                        
                        tmp = CircularSparseMat(tmpPosStrnds, tmpLengths, size(tmp), 1);
                        tmp = this.mergeAdjacentRegions(tmp);
                    case 'linkingNumbers'
                        throw(MException('Chromosome:tooLazy', 'You''re going to have to write this yourself.'));
                    otherwise
                        [positionsStrands, lengths] = find(this.polymerizedRegions(:, this.strandIndexs_negative(2), :));            
                        [tmpPosStrnds, tmpVals] = find(tmp);
                        
                        positionsStrands = [positionsStrands ones(size(positionsStrands,1), 3-size(positionsStrands,2))];
                        tmpPosStrnds = [tmpPosStrnds ones(size(tmpPosStrnds,1), 3-size(tmpPosStrnds,2))];
                        for j = size(positionsStrands, 1)
                            idx1 = ...
                                tmpPosStrnds(:, 1) >= positionsStrands(j, 1) & ...
                                tmpPosStrnds(:, 1) <= positionsStrands(j, 1) + lengths(j)-1 & ...
                                tmpPosStrnds(:, 2) == this.strandIndexs_negative(1) & ...
                                all(tmpPosStrnds(:, 3:end) == positionsStrands(j, 3:end), 2);
                            idx2 = ...
                                tmpPosStrnds(:, 1) >= positionsStrands(j, 1) & ...
                                tmpPosStrnds(:, 1) <= positionsStrands(j, 1) + lengths(j)-1 & ...
                                tmpPosStrnds(:, 2) == this.strandIndexs_negative(2) & ...
                                all(tmpPosStrnds(:, 3:end) == positionsStrands(j, 3:end), 2);
                            tmpPosStrnds(idx1, 2) = this.strandIndexs_negative(2);
                            tmpPosStrnds(idx2, 2) = this.strandIndexs_negative(1);
                        end
                        tmp = CircularSparseMat(tmpPosStrnds, tmpVals, size(tmp), 1);
                end
                varargout{i} = tmp;
            end
        end
        
        %Calculates the footprint of a group of proteins as the maximum of their
        %individual footprints. Also calculates the number of bases that the
        %footprint spans 3'- and 5'- to the centroid base of the footprint.
        function [footprint, footprint3Prime, footprint5Prime, bindingStrandedness, regionStrandedness] = getDNAFootprint(this, monomers, complexs)
            footprint = max([1;
                this.monomerDNAFootprints(monomers);
                this.complexDNAFootprints(complexs)]);
            
            [footprint3Prime, footprint5Prime] = this.calculateFootprintOverhangs(footprint);
            
            bindingStrandedness = [this.monomerDNAFootprintBindingStrandedness(monomers); this.complexDNAFootprintBindingStrandedness(complexs)];
            if isempty(bindingStrandedness)
                bindingStrandedness = this.dnaStrandedness_dsDNA;
            elseif ~isscalar(bindingStrandedness)
                if any(bindingStrandedness ~= bindingStrandedness(1))
                    throw(MException('Chromosome:invalidInput','Footprints must be of a single strandedness'));
                else
                    bindingStrandedness = bindingStrandedness(1);
                end
            end
            
            regionStrandedness = [this.monomerDNAFootprintRegionStrandedness(monomers); this.complexDNAFootprintRegionStrandedness(complexs)];
            if isempty(regionStrandedness)
                regionStrandedness = this.dnaStrandedness_dsDNA;
            elseif ~isscalar(regionStrandedness)
                if any(regionStrandedness ~= regionStrandedness(1))
                    throw(MException('Chromosome:invalidInput','Footprints must be of a single strandedness'));
                else
                    regionStrandedness = regionStrandedness(1);
                end
            end
        end
        
        %return sorted releasableMonomerIndexs, releasableComplexIndexs
        function [releasableMonomerIndexs, releasableComplexIndexs] = getReleasableProteins(this, bindingMonomers, bindingComplexs)
            releasableMonomerIndexs = this.reactionBoundMonomer .* (...
                sum(this.reactionMonomerCatalysisMatrix(:, bindingMonomers), 2) + ...
                sum(this.reactionComplexCatalysisMatrix(:, bindingComplexs), 2) >= ...
                this.reactionThresholds);
            
            if any(releasableMonomerIndexs)
                releasableMonomerIndexs = sort(releasableMonomerIndexs);
                if ~isempty(releasableMonomerIndexs)
                    releasableMonomerIndexs = releasableMonomerIndexs([diff(releasableMonomerIndexs) ~= 0; true], 1);
                end
            else
                releasableMonomerIndexs = 0;
            end
            
            releasableComplexIndexs = this.reactionBoundComplex .* (...
                sum(this.reactionMonomerCatalysisMatrix(:, bindingMonomers), 2) + ...
                sum(this.reactionComplexCatalysisMatrix(:, bindingComplexs), 2) >= ...
                this.reactionThresholds);
            if any(releasableComplexIndexs)
                releasableComplexIndexs = sort(releasableComplexIndexs);
                if ~isempty(releasableComplexIndexs)
                    releasableComplexIndexs = releasableComplexIndexs([diff(releasableComplexIndexs) ~= 0; true], 1);
                end
            else
                releasableComplexIndexs = 0;
            end
        end
        
        %Estimates regions that are accessible to a given protein species. That is
        %regions:
        %- that are polymerized, and have the strandedness that the protein requires
        %- that are not already bound by the specified protein species
        function [rgnPosStrnds, rgnLens] = getAccessibleRegions(this, monomerIdx, complexIdx, checkRegionSupercoiled)
            if nargin < 4
                checkRegionSupercoiled = false;
            end
            
            %find lengths and type of footprint
            [dnaFtpt, ~, ~, dnaFtptBindingStrandedness, dnaFtptRegionStrandedness] = this.getDNAFootprint(monomerIdx, complexIdx);
            
            %initialize excluded regions
            excPosStrnds = [];
            excLens = [];
            
            %exclude regions which are bound by proteins
            [releasableMonomerIndexs, releasableComplexIndexs] = this.getReleasableProteins(monomerIdx, complexIdx);
            
            [monPosStrnds, mons] = find(this.monomerBoundSites);
            [cpxPosStrnds, cpxs] = find(this.complexBoundSites);
            
            iMonSS = find(~ismembc(mons, releasableMonomerIndexs) & this.monomerDNAFootprintBindingStrandedness(mons) == this.dnaStrandedness_ssDNA);
            iMonDS = find(~ismembc(mons, releasableMonomerIndexs) & this.monomerDNAFootprintBindingStrandedness(mons) == this.dnaStrandedness_dsDNA);
            iCpxSS = find(~ismembc(cpxs, releasableComplexIndexs) & this.complexDNAFootprintBindingStrandedness(cpxs) == this.dnaStrandedness_ssDNA);
            iCpxDS = find(~ismembc(cpxs, releasableComplexIndexs) & this.complexDNAFootprintBindingStrandedness(cpxs) == this.dnaStrandedness_dsDNA);
            
            excPosStrnds = [
                excPosStrnds
                monPosStrnds(iMonSS, :)
                monPosStrnds(iMonDS, 1) 2*ceil(monPosStrnds(iMonDS, 2)/2)-1
                monPosStrnds(iMonDS, 1) 2*ceil(monPosStrnds(iMonDS, 2)/2)
                cpxPosStrnds(iCpxSS, :)
                cpxPosStrnds(iCpxDS, 1) 2*ceil(cpxPosStrnds(iCpxDS, 2)/2)-1
                cpxPosStrnds(iCpxDS, 1) 2*ceil(cpxPosStrnds(iCpxDS, 2)/2)
                ];
            excLens = [
                excLens
                this.monomerDNAFootprints(mons(iMonSS))
                this.monomerDNAFootprints(mons(iMonDS))
                this.monomerDNAFootprints(mons(iMonDS))
                this.complexDNAFootprints(cpxs(iCpxSS))
                this.complexDNAFootprints(cpxs(iCpxDS))
                this.complexDNAFootprints(cpxs(iCpxDS))
                ];
            
            %exclude regions which are damaged
            dmgPosStrnds = find(this.damagedSites);
            
            excPosStrnds = [
                excPosStrnds;
                dmgPosStrnds];
            excLens = [
                excLens;
                ones(size(dmgPosStrnds, 1), 1)];
            
            %exclude other strand if proteins by both strands
            if dnaFtptBindingStrandedness == this.dnaStrandedness_dsDNA
                excPosStrnds(:, 2) = ceil(excPosStrnds(:, 2) / 2);
            end
            
            %find polymerized regions
            switch dnaFtptRegionStrandedness
                case this.dnaStrandedness_ssDNA
                    if dnaFtptBindingStrandedness == this.dnaStrandedness_ssDNA
                        [polRgnPosStrnds, polRgnLens] = find(this.singleStrandedRegions);
                    else
                        throw(MException('Chromosome:error', 'unsupported strandedness: %s', dnaFtptRegionStrandedness));
                    end
                case this.dnaStrandedness_dsDNA
                    if dnaFtptBindingStrandedness == this.dnaStrandedness_dsDNA
                        [polRgnPosStrnds, polRgnLens] = find(this.doubleStrandedRegions(:, 1:2:end));
                        if checkRegionSupercoiled
                            tfs = this.supercoiled([polRgnPosStrnds(:, 1) 2 * polRgnPosStrnds(:, 2)]);
                            polRgnPosStrnds = polRgnPosStrnds(tfs, :);
                            polRgnLens = polRgnLens(tfs, :);
                        end
                    else
                        [polRgnPosStrnds, polRgnLens] = find(this.doubleStrandedRegions);
                        if checkRegionSupercoiled
                            tfs = this.supercoiled(polRgnPosStrnds);
                            polRgnPosStrnds = polRgnPosStrnds(tfs, :);
                            polRgnLens = polRgnLens(tfs, :);
                        end
                    end
                otherwise
                    throw(MException('Chromosome:error', 'unsupported strandedness: %s', dnaFtptRegionStrandedness));
            end
            
            %find difference of polymerized and excluded regions
            [rgnPosStrnds, rgnLens] = this.excludeRegions(polRgnPosStrnds, polRgnLens, excPosStrnds, excLens);
            
            %return strands
            if dnaFtptBindingStrandedness == this.dnaStrandedness_dsDNA
                rgnPosStrnds(:, 2) = 2 * rgnPosStrnds(:, 2) - 1;
            end
            
            %return only regions with length at least dnaFtpt
            idx = find(rgnLens >= dnaFtpt);
            rgnPosStrnds = rgnPosStrnds(idx, :);
            rgnLens = rgnLens(idx, :);
        end
        
        function value = getDamagedSites(this, includeBases, includeBonds, includeBase5Prime, includeBond5Prime, includeBase3Prime, includeBond3Prime, includeM6AD)
            import edu.stanford.covert.util.CircularSparseMat;
            
            if ~includeBases
                value = CircularSparseMat([], [], size(this.gapSites), 1);
            else
                %damaged based
                if includeM6AD
                    value = this.damagedBases;
                else
                    [posStrnds, dmgs] = find(this.damagedBases);
                    posStrnds(dmgs == this.metabolite.m6ADIndexs, :) = [];
                    dmgs(dmgs == this.metabolite.m6ADIndexs, :) = [];
                    value = CircularSparseMat(posStrnds, dmgs, [this.sequenceLen this.nCompartments], 1);
                end
                
                if nnz(this.gapSites)
                    value = value + this.gapSites;
                end
                if nnz(this.abasicSites)
                    value = value + this.abasicSites;
                end
                if nnz(this.damagedSugarPhosphates)
                    value = value + this.damagedSugarPhosphates;
                end
                if nnz(this.intrastrandCrossLinks)
                    value = value + this.intrastrandCrossLinks;
                    if includeBase5Prime
                        value = value + this.intrastrandCrossLinks5;
                    end
                    if includeBase3Prime
                        value = value + this.intrastrandCrossLinks3;
                    end
                end
            end
            
            if includeBonds
                if nnz(this.strandBreaks)
                    value = value + this.strandBreaks;
                    if includeBond5Prime
                        value = value + this.strandBreaks5;
                    end
                    if includeBond3Prime
                        value = value + this.strandBreaks3;
                    end
                end
                if nnz(this.hollidayJunctions)
                    value = value + this.hollidayJunctions;
                    if includeBond5Prime
                        value = value + this.hollidayJunctions5;
                    end
                    if includeBond3Prime
                        value = value + this.hollidayJunctions3;
                    end
                end
            end
        end
    end
    
    %additional methods which query, but do not modify state
    %unlike those above, these methods are very focused; they return information
    %about the Chromosome useful to only specific processes
    methods
        function [unmethylatedSites, hemimethylatedSites, methylatedSites, cleavedSites, inaccessibleRegions] = rmStatus(this, ...
                sites, methylatedPositions, restrictionPositions, bindingMonomers, bindingComplexs)
            import edu.stanford.covert.util.SparseMat;
            
            warningState = warning('query', 'SparseMat:inefficient');
            warning('off', 'SparseMat:inefficient');
            
            if nargin < 5
                bindingMonomers = [];
            end
            if nargin < 6
                bindingComplexs = [];
            end
            
            nonmethylatedPositions = [
                1:methylatedPositions(1)-1 methylatedPositions(1)+1:size(sites,2);
                1:methylatedPositions(2)-1 methylatedPositions(2)+1:size(sites,2)]';
            nonRestrictionPositions = [
                1:restrictionPositions(1)-1 restrictionPositions(1)+1:size(sites,2);
                1:restrictionPositions(2)-1 restrictionPositions(2)+1:size(sites,2)]';
            
            sitesChromosomesPositions = [
                size(sites, 1)
                size(sites, 2)
                size(this.damagedBases, 2)/2
                size(this.damagedBases, 3)
                ]';
            
            damagedRegions = ...
                (                this.damagedBases(                  sites(:,methylatedPositions(1)),                   1:2:end,:)~=this.metabolite.m6ADIndexs >0 & ...
                this.damagedBases(                  sites(:,methylatedPositions(1)),                   1:2:end,:)~=0                                       >0)   | ...
                (                this.damagedBases(                  sites(:,methylatedPositions(2)),                   2:2:end,:)~=this.metabolite.m6ADIndexs >0 & ...
                this.damagedBases(                  sites(:,methylatedPositions(2)),                   2:2:end,:)~=0                                       >0)   | ...
                xor(             this.strandBreaks(                  sites(:,restrictionPositions(1)),                  1:2:end,:), ...
                this.strandBreaks(                  sites(:,restrictionPositions(2)),                  2:2:end,:))                                               | ...
                collapse(reshape(this.damagedBases(          reshape(sites(:,nonmethylatedPositions(:,1)),        [],1), 1:2:end,:), sitesChromosomesPositions-[0 1 0 0]),2) >0   | ...
                collapse(reshape(this.damagedBases(          reshape(sites(:,nonmethylatedPositions(:,2)),        [],1), 2:2:end,:), sitesChromosomesPositions-[0 1 0 0]),2) >0   | ...
                this.damagedBases(                  sites(:,methylatedPositions(1)),                   2:2:end,:) ~=0                                            | ...
                this.damagedBases(                  sites(:,methylatedPositions(2)),                   1:2:end,:) ~=0                                            | ...
                collapse(reshape(this.strandBreaks(          reshape(sites(:,nonRestrictionPositions(1:end-1,1)),[],1), 1:2:end,:), sitesChromosomesPositions-[0 2 0 0]),2) >0    | ...
                collapse(reshape(this.strandBreaks(          reshape(sites(:,nonRestrictionPositions(2:end  ,2)),[],1), 2:2:end,:), sitesChromosomesPositions-[0 2 0 0]),2) >0    | ...
                this.strandBreaks(                  sites(:,restrictionPositions(1)),                  2:2:end,:) ~=0                                            | ...
                this.strandBreaks(                  sites(:,restrictionPositions(2)),                  1:2:end,:) ~=0                                            | ...
                collapse(reshape(this.intrastrandCrossLinks( reshape(sites,                                      [],1), 1:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.intrastrandCrossLinks( reshape(sites,                                      [],1), 2:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.damagedSugarPhosphates(reshape(sites,                                      [],1), 1:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.damagedSugarPhosphates(reshape(sites,                                      [],1), 2:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.abasicSites(           reshape(sites,                                      [],1), 1:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.abasicSites(           reshape(sites,                                      [],1), 2:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.gapSites(              reshape(sites,                                      [],1), 1:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.gapSites(              reshape(sites,                                      [],1), 2:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.hollidayJunctions(     reshape(sites(:,1:end-1),                           [],1), 1:2:end,:), sitesChromosomesPositions-[0 1 0 0]),2) >0    | ...
                collapse(reshape(this.hollidayJunctions(     reshape(sites(:,2:end),                             [],1), 2:2:end,:), sitesChromosomesPositions-[0 1 0 0]),2) >0;
            
            methylatedSites = ...
                this.damagedBases(sites(:, methylatedPositions(1)), 1:2:end,:) == this.metabolite.m6ADIndexs & ...
                this.damagedBases(sites(:, methylatedPositions(2)), 2:2:end,:) == this.metabolite.m6ADIndexs;
            
            hemimethylatedSites = xor(...
                this.damagedBases(sites(:, methylatedPositions(1)), 1:2:end,:) == this.metabolite.m6ADIndexs, ...
                this.damagedBases(sites(:, methylatedPositions(2)), 2:2:end,:) == this.metabolite.m6ADIndexs);
            
            unmethylatedSites = ...
                this.damagedBases(sites(:, methylatedPositions(1)), 1:2:end,:) ~= this.metabolite.m6ADIndexs & ...
                this.damagedBases(sites(:, methylatedPositions(2)), 2:2:end,:) ~= this.metabolite.m6ADIndexs;
            
            cleavedSites = ...
                this.strandBreaks(sites(:, restrictionPositions(1)), 1:2:end,:) & ...
                this.strandBreaks(sites(:, restrictionPositions(2)), 2:2:end,:);
            
            siteLength = size(sites, 2);
            ignoreDamageFilter =  1:siteLength;
            sitesStrands = [
                sites(:, 1)   ones(size(sites, 1), 1)
                sites(:, 1) 3*ones(size(sites, 1), 1)];
            inaccessibleRegions = damagedRegions | ...
                (cleavedSites & ~unmethylatedSites) | ...
                ~reshape(this.isRegionAccessible(sitesStrands, siteLength, bindingMonomers, bindingComplexs, true, ignoreDamageFilter, false, true), [], this.nCompartments / 2);
                        
            methylatedSites     = ~inaccessibleRegions & methylatedSites;
            hemimethylatedSites = ~inaccessibleRegions & hemimethylatedSites;
            unmethylatedSites   = ~inaccessibleRegions & unmethylatedSites & ~cleavedSites;
            cleavedSites        = ~inaccessibleRegions & cleavedSites;            
            
            unmethylatedSites   = SparseMat(unmethylatedSites);
            hemimethylatedSites = SparseMat(hemimethylatedSites);
            methylatedSites     = SparseMat(methylatedSites);
            cleavedSites        = SparseMat(cleavedSites);
            inaccessibleRegions = SparseMat(inaccessibleRegions);
            
            if strcmp(warningState.state, 'on'); warning('on', 'SparseMat:inefficient'); end;
        end
    end

    %methods which modify the state of this class, and possibly request
    %modifications to other parts of the simulation's state
    methods
        function sideEffects = setRegionUnwound(this, positions, lengths)
            sideEffects = edu.stanford.covert.cell.sim.SimulationStateSideEffect.empty(0, 1);

            if ~isequal(size(positions, 1), size(lengths, 1))
                throw(MException('Chromosome:invalidInput', 'positions and lengths must have same number of rows'));
            end
            if size(positions, 2) ~= 1
                throw(MException('Chromosome:invalidInput', 'positions must have 1 columns'));
            end
            if size(lengths, 2) ~= 1
                throw(MException('Chromosome:invalidInput', 'lengths must have 1 column'));
            end
            if ~any(lengths)
                return;
            end
            L = this.sequenceLen;
            if any(positions > L)
                throw(MException('Chromosome:invalidInput', 'positions cannot wrap ORI'));
            end
            
            positions = positions(logical(lengths));
            lengths = lengths(logical(lengths));
            n = size(positions, 1);
            if ~all(this.isRegionDoubleStranded([positions ones(n, 1)], lengths, false))
                throw(MException('Chromosome:invalidInput', 'regions must be double-stranded'));
            end
            if ~all(positions == 1 | positions == L | ...
                    this.isRegionSingleStranded([max(1, positions - 1) ones(n, 1)], 1, false) | ...
                    this.isRegionSingleStranded([min(L, positions + 1) ones(n, 1)], 1, false))
                throw(MException('Chromosome:invalidInput',...
                    'unwinding must begin at either end of dsDNA or continue where it left off'));
            end
            if ~all(this.isRegionNotPolymerized([
                    positions 3*ones(n, 1);
                    positions 4*ones(n, 1)],...
                    [lengths; lengths], false))
                throw(MException('Chromosome:invalidInput','chromosome 2 region cannot be polymerized'));
            end
            
            oldStrd = this.strandIndexs_ch1(2);
            newStrd = this.strandIndexs_ch2(2);
            
            directions = sign(lengths);
            positions = positions + min(0, lengths + 1);
            lengths = abs(lengths);
            
            for i = 1:n
                len = lengths(i);
                pos = positions(i,1);
                dir = directions(i);
                
                %if necessary, move region of initial negative strand of
                %chromosome 1 to chromosome 2
                this.monomerBoundSites      = this.shiftStrandToNewChromosome(this.monomerBoundSites,      pos, len, oldStrd, newStrd);
                this.complexBoundSites      = this.shiftStrandToNewChromosome(this.complexBoundSites,      pos, len, oldStrd, newStrd);
                this.gapSites               = this.shiftStrandToNewChromosome(this.gapSites,               pos, len, oldStrd, newStrd);
                this.abasicSites            = this.shiftStrandToNewChromosome(this.abasicSites,            pos, len, oldStrd, newStrd);
                this.damagedSugarPhosphates = this.shiftStrandToNewChromosome(this.damagedSugarPhosphates, pos, len, oldStrd, newStrd);
                this.damagedBases           = this.shiftStrandToNewChromosome(this.damagedBases,           pos, len, oldStrd, newStrd);
                this.intrastrandCrossLinks  = this.shiftStrandToNewChromosome(this.intrastrandCrossLinks,  pos, len, oldStrd, newStrd);
                this.strandBreaks           = this.shiftStrandToNewChromosome(this.strandBreaks,           pos, len, oldStrd, newStrd);
                this.hollidayJunctions      = this.shiftStrandToNewChromosome(this.hollidayJunctions,      pos, len, oldStrd, newStrd);                
                
                %update region of chromosomes that have been polymerized                
                [regionStartPositions, regionLengths] = find(this.polymerizedRegions(:, oldStrd));
                regionStartPositions = regionStartPositions(:, 1);
                idx = find(regionStartPositions <= pos & regionStartPositions + regionLengths > pos);                
                this.polymerizedRegions(regionStartPositions(idx), oldStrd) = pos - regionStartPositions(idx);
                if dir == 1
                    if pos ~= regionStartPositions(idx)
                        throw(MException('Chromosome:error', 'programmer error: unwinding bad region'));
                    end
                    this.linkingNumbers(pos + len, 1:2) = this.linkingNumbers(pos, oldStrd);
                    this.linkingNumbers(pos, 1:2) = 0;
                end
                if (dir == 1 && len == regionLengths(idx) || ...
                    dir == -1 && pos == regionStartPositions(idx)) && ...
                    abs(this.linkingNumbers([pos oldStrd])) > 1e-6
                    throw(MException('Chromosome:invalidInput',...
                        'cannot completely unwind region with nonzero linking number'));
                end
                if regionLengths(idx) - (pos - regionStartPositions(idx)) - len ~= 0
                    this.polymerizedRegions(pos+len, oldStrd) = ...
                        regionLengths(idx) - (pos - regionStartPositions(idx)) - len;
                end
                this.polymerizedRegions(pos, newStrd) = len;
                
                this.mergeOwnAdjacentRegions();
            end
        end
                
        function spmat = shiftStrandToNewChromosome(~, spmat, pos, len, oldStrd, newStrd)
            import edu.stanford.covert.util.CircularSparseMat;
            
            [subs, vals] = find(spmat);
            tfs = ...
                subs(:, 1) >= pos & ...
                subs(:, 1) <= pos + len - 1 & ...
                subs(:, 2) == oldStrd;
            subs(tfs, 2) = newStrd;
            spmat = CircularSparseMat(subs, vals, size(spmat), 1);
        end
        
        function sideEffects = setRegionPolymerized(this, positionsStrands, lengths)
            import edu.stanford.covert.cell.sim.SimulationStateSideEffect;

            L = this.sequenceLen;            
            if ~all(positionsStrands(:,2) == 1 | positionsStrands(:,2) == 2)
                throw(MException('Chromosome:invalidInput', 'positionsStrands must be a valid template strand (eg. 1 or 2)'));
            end
            if any(positionsStrands(:,1) + lengths -1 > L)
                throw(MException('Chromosome:invalidInput', 'positionsStrands cannot wrap ORI'));
            end
            if ~isequal(size(positionsStrands,1), size(lengths,1))
                throw(MException('Chromosome:invalidInput', 'positionsStrands and lengths must have same number of rows'));
            end
            if size(positionsStrands, 2) ~= 2
                throw(MException('Chromosome:invalidInput', 'positionsStrands must have 2 columns'));
            end
            if size(lengths, 2) ~= 1
                throw(MException('Chromosome:invalidInput', 'lengths must have 1 column'));
            end
            
            positionsStrands(:, 1) = positionsStrands(: ,1) + min(0, lengths + 1);
            lengths = abs(lengths);
            
            for i = 1:size(positionsStrands, 1)
                len = lengths(i);
                pos = positionsStrands(i,1);
                tmpStrd = this.strandIndexs_template(positionsStrands(i,2), :);
                nonTmpStrd = this.strandIndexs_nonTemplate(positionsStrands(i,2), :);
                
                %if no polymerization requested, do nothing
                if len == 0
                    continue;
                end

                %check that template exists
                if ~this.isRegionPolymerized([pos tmpStrd], len, false)
                    throw(MException('Chromosome:error','cannot polymerize a region without a template'));
                end

                %check that strand hasn't been polymerized
                if ~this.isRegionNotPolymerized([pos nonTmpStrd], len, false)
                    throw(MException('Chromosome:error','cannot polymerize a region that''s already been polymerized'));
                end

                %set region polymerized
                this.polymerizedRegions([pos nonTmpStrd]) = len;
                this.linkingNumbers([pos tmpStrd; pos nonTmpStrd]) = len / this.relaxedBasesPerTurn;
                this.mergeOwnAdjacentRegions();
            end

            %side effects
            sideEffects = SimulationStateSideEffect.empty(0, 1);
        end

        %lengths must be non-negative integers
        function [releasedMonomers, releasedComplexs, sideEffects, tfs, idxs, positionsStrands, lengths] = ...
                setSiteProteinBound(this, positionsStrands, maxBindings, weights, binding_monomerIndexs, binding_complexIndexs, ...
                mainEffects_monomerIndexs, mainEffects_complexIndexs, ...
                isBindingStable, isPositionsStrandFootprintCentroid, lengths, isBindingProcessive, ignoreDamageFilter, checkRegionSupercoiled)
                
            if nargin < 14
                checkRegionSupercoiled = false;
            end
            
            if isBindingStable && numel(binding_monomerIndexs) + numel(binding_complexIndexs) > 1
                throw(MException('Chromosome:invalidInput', 'Can only bind 1 protein at a time'));
            end
            if isBindingStable && ~isBindingProcessive && ~all(lengths == 1)
                throw(MException('Chromosome:invalidInput','If any(lengths~=1) then if binding is stable it must be processive'));
            end
            if numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            
            %randomly select among accessible sites
            [tfs, idxs, positionsStrands, lengths] = this.sampleAccessibleRegions(maxBindings, weights, positionsStrands, lengths, ...
                binding_monomerIndexs, binding_complexIndexs, isPositionsStrandFootprintCentroid, ignoreDamageFilter, ~isBindingStable, isBindingProcessive, checkRegionSupercoiled);
            
            %if positions are centroid, shift positions to start coordinate view
            [footprint, footprint3Prime, footprint5Prime, footprintBindingStrandedness, footprintRegionStrandedness] = this.getDNAFootprint(binding_monomerIndexs, binding_complexIndexs);
            
            releasePositionsStrands = positionsStrands;
            releasePositionsStrands(:, 1) = releasePositionsStrands(:, 1) + min(0, lengths + 1);
            if isPositionsStrandFootprintCentroid
                releasePositionsStrands(mod(releasePositionsStrands(:,2),2)==1, 1) = releasePositionsStrands(mod(releasePositionsStrands(:,2),2)==1, 1) - footprint5Prime;
                releasePositionsStrands(mod(releasePositionsStrands(:,2),2)==0, 1) = releasePositionsStrands(mod(releasePositionsStrands(:,2),2)==0, 1) - footprint3Prime;
            end
            
            %remove proteins currently bound to selected accessible sites
            if footprintRegionStrandedness == this.dnaStrandedness_xsDNA
                [monomerPosStrands, monomers] = find(this.monomerBoundSites);
                [complexPosStrands, complexs] = find(this.complexBoundSites);
                [releasableMonomerIndexs, releasableComplexIndexs] = this.getReleasableProteins(binding_monomerIndexs, binding_complexIndexs);
                
                tmpIdxs = find(~ismembc(monomers, releasableMonomerIndexs));
                monomerPosStrands = monomerPosStrands(tmpIdxs, :);
                monomers = monomers(tmpIdxs, :);
                
                tmpIdxs = find(~ismembc(complexs, releasableComplexIndexs));
                complexPosStrands = complexPosStrands(tmpIdxs, :);
                complexs = complexs(tmpIdxs, :);
                
                [releasePositionsStrands, releaseLens] = this.excludeRegions([
                    releasePositionsStrands(:, 1)     ones(size(releasePositionsStrands, 1), 1)
                    releasePositionsStrands(:, 1) 2 * ones(size(releasePositionsStrands, 1), 1)
                    releasePositionsStrands(:, 1) 3 * ones(size(releasePositionsStrands, 1), 1)
                    releasePositionsStrands(:, 1) 4 * ones(size(releasePositionsStrands, 1), 1)], ...
                    repmat(footprint + abs(lengths) - 1, 4, 1), ...
                    [monomerPosStrands; complexPosStrands], [this.monomerDNAFootprints(monomers); this.complexDNAFootprints(complexs)]);
            else
                releaseLens = footprint + abs(lengths) - 1;
            end
            [releasedMonomers, releasedComplexs, sideEffects] = ...
                this.setRegionProteinUnbound(releasePositionsStrands, releaseLens, ...
                mainEffects_monomerIndexs, mainEffects_complexIndexs, ...                
                footprintRegionStrandedness == this.dnaStrandedness_dsDNA, ...
                footprintBindingStrandedness == this.dnaStrandedness_dsDNA, ...
                false, false);
            
            %if binding is stable, bind protein to selected accessible sites
            if isBindingStable
                bindPositionsStrands = positionsStrands;
                if isBindingProcessive
                    bindPositionsStrands(:,1) = bindPositionsStrands(:,1) + lengths - sign(lengths);
                end
                if isPositionsStrandFootprintCentroid
                    bindPositionsStrands(mod(bindPositionsStrands(:,2),2)==1, 1) = bindPositionsStrands(mod(bindPositionsStrands(:,2),2)==1, 1) - footprint5Prime;
                    bindPositionsStrands(mod(bindPositionsStrands(:,2),2)==0, 1) = bindPositionsStrands(mod(bindPositionsStrands(:,2),2)==0, 1) - footprint3Prime;
                end

                if ~isempty(binding_monomerIndexs)
                    this.monomerBoundSites(bindPositionsStrands) = binding_monomerIndexs;
                    releasedMonomers(mainEffects_monomerIndexs == binding_monomerIndexs) = ...
                        releasedMonomers(mainEffects_monomerIndexs == binding_monomerIndexs) - ...
                        numel(idxs);
                else
                    this.complexBoundSites(bindPositionsStrands) = binding_complexIndexs;
                    releasedComplexs(mainEffects_complexIndexs == binding_complexIndexs) = ...
                        releasedComplexs(mainEffects_complexIndexs == binding_complexIndexs) - ...
                        numel(idxs);
                end
            end
        end
        
        %change the identity of bound proteins
        function [releasedMonomers, releasedComplexs, sideEffects] = ...
                modifyBoundProtein(this, positionsStrands, newMonomers, newComplexs, ...
                mainEffects_monomerIndexs, mainEffects_complexIndexs)
             
            import edu.stanford.covert.cell.sim.SimulationStateSideEffect;
            import edu.stanford.covert.cell.sim.SimulationStateSideEffectItem;
            import edu.stanford.covert.util.countUnique;
            
            if any(all([newMonomers newComplexs], 2))
                throw(MException('Only 1 protein can be bound at each site'));
            end
            
            %get identities of currently bound proteins
            oldMonomers = this.monomerBoundSites(positionsStrands);
            oldComplexs = this.complexBoundSites(positionsStrands);
            
            %check that changing protein identify doesn't increase footprint
            %size (otherwise setSiteProteinBound should be used)
            oldFootprints = zeros(size(positionsStrands, 1));
            oldFootprints(oldMonomers ~= 0) = this.monomerDNAFootprints(oldMonomers(oldMonomers~=0));
            oldFootprints(oldComplexs ~= 0) = this.complexDNAFootprints(oldComplexs(oldComplexs~=0));
            
            newFootprints = zeros(size(positionsStrands, 1));
            newFootprints(newMonomers ~= 0) = this.monomerDNAFootprints(newMonomers(newMonomers~=0));
            newFootprints(newComplexs ~= 0) = this.complexDNAFootprints(newComplexs(newComplexs~=0));
            
            if any(newFootprints > oldFootprints)
                throw(MException('New protein footprints cannot be larger that old ones'));
            end
            
            %modify identities of bound proteins
            this.monomerBoundSites(positionsStrands(newMonomers~=0, :)) = newMonomers(newMonomers~=0, :);
            this.complexBoundSites(positionsStrands(newComplexs~=0, :)) = newComplexs(newComplexs~=0, :);
            
            %summarize bound/released monomers/complexs
            releasedMonomers = zeros(size(mainEffects_monomerIndexs));
            releasedComplexs = zeros(size(mainEffects_complexIndexs));
            
            oldMonomers = oldMonomers(oldMonomers~=0);
            [tfs, idxs] = ismember(oldMonomers, mainEffects_monomerIndexs);
            if ~all(tfs)
                throw(MException('Chromosome:invalidInput','All binding/releasing proteins must be mainEffects'));
            end
            [idxs, counts] = countUnique(idxs);
            releasedMonomers(idxs) = releasedMonomers(idxs) + counts;
            
            newMonomers = newMonomers(newMonomers~=0);
            [tfs, idxs] = ismember(newMonomers, mainEffects_monomerIndexs);
            if ~all(tfs)
                throw(MException('Chromosome:invalidInput','All binding/releasing proteins must be mainEffects'));
            end
            [idxs, counts] = countUnique(idxs);
            releasedMonomers(idxs) = releasedMonomers(idxs) - counts;
            
            oldComplexs = oldComplexs(oldComplexs~=0);
            [tfs, idxs] = ismember(oldComplexs, mainEffects_complexIndexs);
            if ~all(tfs)
                throw(MException('Chromosome:invalidInput','All binding/releasing proteins must be mainEffects'));
            end
            [idxs, counts] = countUnique(idxs);
            releasedComplexs(idxs) = releasedComplexs(idxs) + counts;
            
            newComplexs = newComplexs(newComplexs~=0);
            [tfs, idxs] = ismember(newComplexs, mainEffects_complexIndexs);
            if ~all(tfs)
                throw(MException('Chromosome:invalidInput','All binding/releasing proteins must be mainEffects'));
            end
            [idxs, counts] = countUnique(idxs);
            releasedComplexs(idxs) = releasedComplexs(idxs) - counts;
            
            sideEffects = SimulationStateSideEffect.empty(0, 1);
        end
        
        function [releasedMonomers, releasedComplexs, sideEffects] = setRegionProteinUnbound(this, ...
                positionsStrands, lengths, mainEffects_monomerIndexs, mainEffects_complexIndexs, ...
                regionBothStrands, bindingBothStrands, suspendExternalStateUpdating, proteinIsDegraded)
            import edu.stanford.covert.cell.sim.SimulationStateSideEffect;
            import edu.stanford.covert.cell.sim.SimulationStateSideEffectItem;
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.util.countUnique;
            
            startCoors = mod(positionsStrands(:, 1) + min(0, lengths+1) - 1, this.sequenceLen) + 1;
            endCoors = startCoors + abs(lengths) - 1;
            strndTimes = [positionsStrands(:, 2:end) ones(size(positionsStrands, 1), size(this.monomerBoundSites, 3) - size(positionsStrands,2))];
            if bindingBothStrands
                strndTimes(:, 1) = ceil(strndTimes(:, 1) / 2); 
            end
            
            tmpIdxs = find(endCoors > this.sequenceLen);
            startCoors = [startCoors; max(1, startCoors(tmpIdxs) - this.sequenceLen)];
            endCoors = [endCoors; endCoors(tmpIdxs) - this.sequenceLen];
            endCoors(tmpIdxs) = this.sequenceLen;
            strndTimes = [strndTimes; strndTimes(tmpIdxs, :)];
            
            tmpIdxs = find(startCoors < 0);
            startCoors = [startCoors; startCoors(tmpIdxs) + this.sequenceLen];
            endCoors = [endCoors; min(this.sequenceLen, endCoors(tmpIdxs) + this.sequenceLen)];
            startCoors(tmpIdxs) = 1;
            strndTimes = [strndTimes; strndTimes(tmpIdxs, :)];
            
            %bound monomers
            [subs, vals] = find(this.monomerBoundSites);
            unbindingMonomers = false(size(vals));
            monomerStarts = subs(:, 1);
            monomerEnds   = monomerStarts + this.monomerDNAFootprints(vals, :) - 1;
            monomerStrndTimes = [subs(:, 2:end) ones(size(subs,1), size(strndTimes,2)-size(subs,2)+1)];
            idxs = (1:numel(vals))';
            if bindingBothStrands
                monomerStrndTimes(:, 1) = ceil(monomerStrndTimes(:, 1) / 2);
            elseif regionBothStrands
                tfs = this.monomerDNAFootprintBindingStrandedness(vals) == this.dnaStrandedness_dsDNA;
                idxs = [idxs(~tfs); idxs(tfs); idxs(tfs)];
                monomerStarts = [monomerStarts(~tfs); monomerStarts(tfs); monomerStarts(tfs)];
                monomerEnds = [monomerEnds(~tfs); monomerEnds(tfs); monomerEnds(tfs)];
                monomerStrndTimes = [
                    monomerStrndTimes(~tfs, :)
                    2*ceil(monomerStrndTimes(tfs, 1)/2)-1 monomerStrndTimes(tfs, 2:end)
                    2*ceil(monomerStrndTimes(tfs, 1)/2) monomerStrndTimes(tfs, 2:end)
                    ];
            end
            
            tmpIdxs = find(monomerEnds > this.sequenceLen);
            monomerStarts = [monomerStarts; ones(numel(tmpIdxs), 1)];
            monomerEnds   = [monomerEnds;   monomerEnds(tmpIdxs, :) - this.sequenceLen];
            monomerStrndTimes = [monomerStrndTimes; monomerStrndTimes(tmpIdxs, :)];
            idxs = [idxs; idxs(tmpIdxs, :)];
            
            for i = 1:numel(startCoors)
                monMask = ...
                    (monomerStarts >= startCoors(i) & monomerStarts <= endCoors(i)) | ...
                    (monomerEnds >= startCoors(i) & monomerEnds <= endCoors(i)) | ...
                    (monomerStarts <= startCoors(i) & monomerEnds >= endCoors(i));
                tmpIdxs = find(monMask);
                if isempty(tmpIdxs); continue; end;
                
                if isvector(strndTimes)
                    tmpIdxs = tmpIdxs(monomerStrndTimes(monMask, 1) ==  strndTimes(i, 1));
                elseif isscalar(tmpIdxs)
                    tmpIdxs = tmpIdxs(all(monomerStrndTimes(monMask, :) == strndTimes(i, :)));
                else
                    tmpIdxs = tmpIdxs(ismember(monomerStrndTimes(monMask, :), strndTimes(i, :), 'rows'));
                end
                if isempty(tmpIdxs); continue; end;
                
                unbindingMonomers(idxs(tmpIdxs)) = true;
            end
            
            releasedMonomers = zeros(numel(mainEffects_monomerIndexs), 1);
            if isempty(unbindingMonomers) || ~any(unbindingMonomers)
                gblIdxs = [];
                counts = [];
            elseif isscalar(unbindingMonomers)
                gblIdxs = vals(unbindingMonomers);
                lclIdxs = find(gblIdxs == mainEffects_monomerIndexs);
                if isempty(lclIdxs)
                    counts = 1;
                else
                    releasedMonomers(lclIdxs, 1) = 1;
                    counts = zeros(0, 1);
                    gblIdxs = zeros(0, 1);
                end
            else
                [gblIdxs, counts] = countUnique(vals(unbindingMonomers, 1));
                [tfs, lclIdxs] = ismember(gblIdxs, mainEffects_monomerIndexs);
                
                releasedMonomers(lclIdxs(tfs)) = counts(tfs);
                
                tfs = tfs | gblIdxs == 0;
                gblIdxs = gblIdxs(~tfs);
                counts = counts(~tfs);
            end

            if isempty(gblIdxs)
                sideEffects_releasedMonomers = SimulationStateSideEffect.empty(0, 1);
            else
                sideEffects_releasedMonomers = SimulationStateSideEffect.empty(numel(gblIdxs), 0);
            end
            for i = 1:numel(gblIdxs)
                sideEffects_releasedMonomers(i, 1) = SimulationStateSideEffect([...
                    SimulationStateSideEffectItem('ProteinMonomer', 'counts', 'matureIndexs', gblIdxs(i), this.compartment.cytosolIndexs,  counts(i)); ...
                    SimulationStateSideEffectItem('ProteinMonomer', 'counts', 'boundIndexs',  gblIdxs(i), this.compartment.cytosolIndexs, -counts(i))]);
            end
            
            unbindingMonomerSubs = subs(unbindingMonomers, :);
            unbindingMonomerVals = vals(unbindingMonomers, :);
            this.monomerBoundSites = CircularSparseMat(subs(~unbindingMonomers, :), vals(~unbindingMonomers, 1), [this.sequenceLen this.nCompartments], 1);
            
            %bound complexes
            [subs, vals] = find(this.complexBoundSites);
            unbindingComplexs = false(size(vals));
            complexStarts = subs(:, 1);
            complexEnds   = complexStarts + this.complexDNAFootprints(vals, :) - 1;
            complexStrndTimes = [subs(:, 2:end) ones(size(subs,1), size(strndTimes,2)-size(subs,2)+1)];
            idxs = (1:numel(vals))';
            if bindingBothStrands
                complexStrndTimes(:, 1) = ceil(complexStrndTimes(:, 1) / 2);
            elseif regionBothStrands
                tfs = this.complexDNAFootprintBindingStrandedness(vals) == this.dnaStrandedness_dsDNA;
                idxs = [idxs(~tfs); idxs(tfs); idxs(tfs)];
                complexStarts = [complexStarts(~tfs); complexStarts(tfs); complexStarts(tfs)];
                complexEnds = [complexEnds(~tfs); complexEnds(tfs); complexEnds(tfs)];
                complexStrndTimes = [
                    complexStrndTimes(~tfs, :)
                    2*ceil(complexStrndTimes(tfs, 1)/2)-1 complexStrndTimes(tfs, 2:end)
                    2*ceil(complexStrndTimes(tfs, 1)/2) complexStrndTimes(tfs, 2:end)
                    ];
            end
            
            tmpIdxs = find(complexEnds > this.sequenceLen);
            complexStarts = [complexStarts; ones(numel(tmpIdxs), 1)];
            complexEnds   = [complexEnds;   complexEnds(tmpIdxs, :) - this.sequenceLen];
            complexStrndTimes    = [complexStrndTimes; complexStrndTimes(tmpIdxs, :)];
            idxs = [idxs; idxs(tmpIdxs, :)];
            
            for i = 1:numel(startCoors)
                comMask = ...
                    (complexStarts >= startCoors(i) & complexStarts <= endCoors(i)) | ...
                    (complexEnds >= startCoors(i) & complexEnds <= endCoors(i)) | ...
                    (complexStarts <= startCoors(i) & complexEnds >= endCoors(i));
                tmpIdxs = find(comMask);
                if isempty(tmpIdxs); continue; end;
                
                if isvector(strndTimes)
                    tmpIdxs = tmpIdxs(complexStrndTimes(comMask, 1) ==  strndTimes(i, 1));
                elseif isscalar(tmpIdxs)
                    tmpIdxs = tmpIdxs(all(complexStrndTimes(comMask, :) == strndTimes(i, :)));
                else
                    tmpIdxs = tmpIdxs(ismember(complexStrndTimes(comMask, :), strndTimes(i, :), 'rows'));
                end
                if isempty(tmpIdxs); continue; end;
                
                unbindingComplexs(idxs(tmpIdxs)) = true;
            end
            
            releasedComplexs = zeros(numel(mainEffects_complexIndexs), 1);
            if isempty(unbindingComplexs) || ~any(unbindingComplexs)
                gblIdxs = [];
                counts = [];
            elseif isscalar(unbindingComplexs)
                gblIdxs = vals(unbindingComplexs);
                lclIdxs = find(gblIdxs == mainEffects_complexIndexs);
                if isempty(lclIdxs)
                    counts = 1;
                else
                    releasedComplexs(lclIdxs, 1) = 1;
                    counts = zeros(0, 1);
                    gblIdxs = zeros(0, 1);
                end
            else
                [gblIdxs, counts] = countUnique(vals(unbindingComplexs, 1));
                [tfs, lclIdxs] = ismember(gblIdxs, mainEffects_complexIndexs);
                
                releasedComplexs(lclIdxs(tfs)) = counts(tfs);
                tfs = tfs | gblIdxs == 0;
                gblIdxs = gblIdxs(~tfs);
                counts = counts(~tfs);
            end
            
            if isempty(gblIdxs)
                sideEffects_releasedComplexs = SimulationStateSideEffect.empty(0, 1);
            else
                sideEffects_releasedComplexs = SimulationStateSideEffect.empty(numel(gblIdxs), 0);
            end
            for i = 1:numel(gblIdxs)
                sideEffects_releasedComplexs(i, 1) = SimulationStateSideEffect([...
                    SimulationStateSideEffectItem('ProteinComplex', 'counts', 'matureIndexs', gblIdxs(i), this.compartment.cytosolIndexs,  counts(i)); ...
                    SimulationStateSideEffectItem('ProteinComplex', 'counts', 'boundIndexs',  gblIdxs(i), this.compartment.cytosolIndexs, -counts(i))]);
            end
            
            unbindingComplexSubs = subs(unbindingComplexs, :);
            unbindingComplexVals = vals(unbindingComplexs, :);
            this.complexBoundSites = CircularSparseMat(subs(~unbindingComplexs, :), vals(~unbindingComplexs, 1), [this.sequenceLen this.nCompartments], 1);
            
            %effects on other states
            if ~suspendExternalStateUpdating
                this.updateExternalState(unbindingMonomerSubs, unbindingMonomerVals, unbindingComplexSubs, unbindingComplexVals, proteinIsDegraded);
            end
            
            %side effects
            sideEffects = [
                sideEffects_releasedMonomers;
                sideEffects_releasedComplexs];
        end
        
        function [positionsStrands, sideEffects] = stochasticallySetProteinUnbound(this, monomerIndex, complexIndex, ...
                rate, protectedPositionsStrands, protectedLengths, suspendExternalStateUpdating, proteinIsDegraded)
            import edu.stanford.covert.cell.sim.SimulationStateSideEffect;
            import edu.stanford.covert.util.CircularSparseMat;
            
            sideEffects = SimulationStateSideEffect.empty(0, 1);
            
            if ~isempty(monomerIndex)
                if ~isempty(complexIndex)
                    throw(MException('Chromosome:invalidInput', 'Can only unbind one protein at a time'));
                end
                [positionsStrands, proteins] = find(this.monomerBoundSites);
                unbinding = proteins == monomerIndex;
            elseif ~isempty(complexIndex)
                [positionsStrands, proteins] = find(this.complexBoundSites);
                unbinding = proteins == complexIndex;
            else
                positionsStrands = zeros(0, 2);
                return;
            end
            if ~any(unbinding)
                positionsStrands = zeros(0, 2);
                return;
            end
            
            %randomly select bound protein to release at specified rate
            if isfinite(rate)
                unbinding(unbinding) = this.randStream.rand(sum(unbinding), 1) < rate;
                if ~any(unbinding)
                    positionsStrands = zeros(0, 2);
                    return;
                end
            end
            
            %exclude protected sites
            if ~isempty(protectedPositionsStrands)
                footprint = this.getDNAFootprint(monomerIndex, complexIndex);
                for i = 1:size(protectedPositionsStrands, 1)
                    unbndPosStrnds = positionsStrands(unbinding, :);
                    
                    tfs = ...
                        (unbndPosStrnds(:, 1)                 >= protectedPositionsStrands(i, 1) & ...
                         unbndPosStrnds(:, 1)                 <= protectedPositionsStrands(i, 1) + protectedLengths(i) - 1) | ...
                        (unbndPosStrnds(: ,1) + footprint - 1 >= protectedPositionsStrands(i, 1) & ...
                         unbndPosStrnds(:, 1) + footprint - 1 <= protectedPositionsStrands(i, 1) + protectedLengths(i) - 1) | ...
                        (unbndPosStrnds(:, 1)                 <= protectedPositionsStrands(i, 1) & ...
                         unbndPosStrnds(:, 1) + footprint - 1 >= protectedPositionsStrands(i, 1) + protectedLengths(i) - 1);
                    
                    tfs(tfs) = unbndPosStrnds(tfs, 2) == protectedPositionsStrands(i, 2);
                    
                    unbinding(unbinding) = ~tfs;
                    
                    if all(tfs)
                        positionsStrands = zeros(0, 2);
                        return;
                    end
                end
            end
            
            %update bound proteins
            if ~isempty(monomerIndex)
                unbindingMonomerSubs = positionsStrands(unbinding, :);
                unbindingMonomerVals = proteins(unbinding, 1);
                unbindingComplexSubs = zeros(0, 2);
                unbindingComplexVals = zeros(0, 1);
                this.monomerBoundSites = CircularSparseMat(positionsStrands(~unbinding, :), proteins(~unbinding, 1), [this.sequenceLen this.nCompartments], 1);
            else
                unbindingMonomerSubs = zeros(0, 2);
                unbindingMonomerVals = zeros(0, 1);
                unbindingComplexSubs = positionsStrands(unbinding, :);
                unbindingComplexVals = proteins(unbinding, 1);
                this.complexBoundSites = CircularSparseMat(positionsStrands(~unbinding, :), proteins(~unbinding, 1), [this.sequenceLen this.nCompartments], 1);
            end
            positionsStrands = positionsStrands(unbinding, :);
            
            %effects on other states
            if ~suspendExternalStateUpdating
                this.updateExternalState(unbindingMonomerSubs, unbindingMonomerVals, unbindingComplexSubs, unbindingComplexVals, proteinIsDegraded);
            end
        end
        
        function updateExternalState(this, ~, ~, unbindingComplexSubs, unbindingComplexVals, proteinIsDegraded)
            posStrnds = unbindingComplexSubs(ismembc(unbindingComplexVals, this.complexIndexs_rnaPolymerase), :);
            if ~isempty(posStrnds)
                [~, footprint3Prime, footprint5Prime] = this.getDNAFootprint([], this.complexIndexs_rnaPolymerase);
                posStrnds( isodd(posStrnds(:, 2)), 1) = posStrnds( isodd(posStrnds(:, 2)), 1) + footprint5Prime;
                posStrnds(iseven(posStrnds(:, 2)), 1) = posStrnds(iseven(posStrnds(:, 2)), 1) + footprint3Prime;
                this.rnaPolymerase.releasePolymerase(posStrnds, proteinIsDegraded);
            end
        end
        
        function [positionsStrands, sideEffects] = setSiteDamaged(this, ...
                damageType, damageSubType, probDamage, maxDamages, ...
                vulnerableMotif, vulnerableMotifType)
            import edu.stanford.covert.cell.sim.SimulationStateSideEffect;
            
            %side effects
            sideEffects = SimulationStateSideEffect.empty(0, 1);
            positionsStrands = zeros(0, 2);
            
            %return if probability of damage is 0
            if probDamage == 0 || maxDamages == 0
                return;
            end
            
            %sample vulnerable sites
            if ischar(vulnerableMotif)
                positionsStrands = this.sampleAccessibleSites(probDamage, maxDamages, vulnerableMotif);
            else
                if nnz(this.(vulnerableMotifType)) == 0
                    return;
                end
                positionsStrands = find( ...
                    vulnerableMotif == this.(vulnerableMotifType) & ...
                    vulnerableMotif == this.damagedSites_nonRedundant);
                if isempty(positionsStrands)
                    return;
                end
                maxDamages = min(maxDamages, ...
                    this.randStream.stochasticRound(size(positionsStrands, 1) * probDamage));
                if maxDamages == 0
                    return;
                end
                positionsStrands = this.randStream.randomlySelectNRows(positionsStrands, maxDamages);
            end
            if isempty(positionsStrands)
                return;
            end
            
            %damage selected sites
            this.(damageType)(positionsStrands) = damageSubType;
        end
    end
    
    %private methods which modify this class' state, and possibly request
    %changes that of other parts of the simulation's state
    methods
        function mergeOwnAdjacentRegions(this)
            import edu.stanford.covert.util.CircularSparseMat;

            if nnz(this.polymerizedRegions) <= 1
                return;
            end
            
            %% add linking numbers for adjacent double-stranded regions
            
            %convert from SparseMat to coordinates
            [dsPosStrands, lengths] = find(this.doubleStrandedRegions);
            [lkPosStrands, oldLKNums] = find(this.linkingNumbers);

            %combine linking numbers for adjacent double-stranded regions
            lkNums = zeros(size(lengths));
            for i = 1:size(lkPosStrands,1)
                idx = find(...
                    lkPosStrands(i, 1) >= dsPosStrands(:, 1)           & ...
                    lkPosStrands(i, 1) <= dsPosStrands(:, 1) + lengths & ...
                    lkPosStrands(i, 2) == dsPosStrands(:, 2));
                lkNums(idx) = lkNums(idx) + oldLKNums(i);
            end
            
            %convert from coordinates to SparseMat
            this.linkingNumbers = CircularSparseMat(...
                dsPosStrands, lkNums, size(this.linkingNumbers), 1);

            %% merge adjacent regions
            this.polymerizedRegions = this.mergeAdjacentRegions(this.polymerizedRegions);            
        end
        
        function polymerizedRegions = mergeAdjacentRegions(~, polymerizedRegions)
            import edu.stanford.covert.util.CircularSparseMat;

            if nnz(polymerizedRegions) <= 1
                return;
            end

            %% merge adjacent regions
            %convert from SparseMat to coordinates
            [positionsStrands, lengths] = find(polymerizedRegions);
            
            %make sure that regions don't overlap
            if any(lengths < 0) || ...
               any(diff(positionsStrands(:,1)) < lengths(1:end-1) & diff(positionsStrands(:,2)) == 0)
                throw(MException('Chromosome:error','polymerizedRegions is corrupt'))
            end
            
            %merge adjacent regions (except over OriC)
            idxs = find(diff(positionsStrands(:,1)) == lengths(1:end-1) & diff(positionsStrands(:,2)) == 0);
            for i = numel(idxs):-1:1
                j = idxs(i);
                lengths(j, 1) = lengths(j, 1) + lengths(j+1, 1);
            end
            positionsStrands(idxs+1, :) = [];
            lengths(idxs+1, :) = [];

            %% convert from coordinates to SparseMat
            polymerizedRegions = CircularSparseMat(positionsStrands, lengths, size(polymerizedRegions), 1);
        end

        function [rgnPosStrnds, rgnLens] = excludeRegions(this, incPosStrnds, incLens, excPosStrnds, excLens)
            L = this.sequenceLen;
            
            %options
            if isscalar(excLens)
                excLens = excLens(ones(size(excPosStrnds,1), 1), 1);
            end
            
            %join included regions
            [incPosStrnds, incLens] = this.joinSplitOverOriCRegions(incPosStrnds, incLens);
            
            %exclude excluded regions
            [excPosStrnds, excLens] = this.joinSplitRegions(excPosStrnds, excLens);            
            
            excPos = [
                excPosStrnds(:, 1) - L;
                excPosStrnds(:, 1);
                excPosStrnds(:, 1) + L];
            excStrnds = [excPosStrnds(:, 2); excPosStrnds(:, 2); excPosStrnds(:, 2)];
            excLens = [excLens; excLens; excLens];
            
            rgnPos = zeros(0, 1);
            rgnEnds = zeros(0, 1);
            rgnStrnds = zeros(0, 1);
            for j = 1:size(incPosStrnds, 1)
                startCoor = incPosStrnds(j, 1);
                endCoor = startCoor + incLens(j) - 1;
                strnd = incPosStrnds(j, 2);
                
                excIdxs = find(...
                    ((excPos <= startCoor & excPos + excLens-1 >= startCoor) | ...
                    (excPos <= endCoor & excPos + excLens-1 >= endCoor) | ...
                    (excPos >= startCoor & excPos + excLens-1 <= endCoor)) & ...
                    excStrnds == strnd);
                
                if isempty(excIdxs)
                    addtlPos = startCoor;
                    addtlEnds = endCoor;
                elseif excPos(excIdxs(1)) <= startCoor
                    if excPos(excIdxs(end)) + excLens(end) - 1 >= endCoor
                        addtlPos = excPos(excIdxs(1:end-1)) + excLens(excIdxs(1:end-1));
                        addtlEnds = excPos(excIdxs(2:end))-1;
                    else
                        addtlPos = excPos(excIdxs) + excLens(excIdxs);
                        addtlEnds = [excPos(excIdxs(2:end))-1; endCoor];
                    end
                else
                    if excPos(excIdxs(end)) + excLens(end) - 1 >= endCoor
                        addtlPos = [startCoor; excPos(excIdxs(1:end-1)) + excLens(excIdxs(1:end-1))];
                        addtlEnds = excPos(excIdxs)-1;
                    else
                        addtlPos = [startCoor; excPos(excIdxs) + excLens(excIdxs)];
                        addtlEnds = [excPos(excIdxs)-1; endCoor];
                    end
                end
                
                rgnPos = [
                    rgnPos;
                    addtlPos];
                rgnEnds = [
                    rgnEnds;
                    addtlEnds];
                rgnStrnds = [
                    rgnStrnds;
                    strnd(ones(size(addtlPos)), 1)];
            end
            
            idx = find(rgnPos > L);
            rgnPos(idx) = rgnPos(idx) - L;
            rgnEnds(idx) = rgnEnds(idx) - L;
            
            idx = find(rgnPos > rgnEnds);
            rgnPos(idx,:) = [];
            rgnEnds(idx,:) = [];
            rgnStrnds(idx,:) = [];
            
            %join split regions
            [rgnPosStrnds, rgnLens] = this.joinSplitOverOriCRegions([rgnPos rgnStrnds], rgnEnds - rgnPos + 1);
            
            %format output
            rgnPosStrnds(:, 1) = mod(rgnPosStrnds(:, 1) - 1, L) + 1;
            [rgnPosStrnds, order] = edu.stanford.covert.util.SparseMat.sort_subs(rgnPosStrnds, [L this.nCompartments]);
            rgnLens = rgnLens(order, :);
        end

        %finds the intersection of two lists of regions
        function [posStrnds, lens] = intersectRegions(this, posStrndsA, lensA, posStrndsB, lensB)
            [posStrndsA lensA] = this.splitOverOriC(posStrndsA, lensA);
            posStrndsA(:, 1) = mod(posStrndsA(:, 1) - 1, this.sequenceLen) + 1;
            [posStrndsA, idxs] = edu.stanford.covert.util.SparseMat.sort_subs(posStrndsA, [this.sequenceLen this.nCompartments]);
            lensA = lensA(idxs);

            [posStrndsB lensB] = this.splitOverOriC(posStrndsB, lensB);            
            posStrndsB(:, 1) = mod(posStrndsB(:, 1) - 1, this.sequenceLen) + 1;
            [posStrndsB, idxs] = edu.stanford.covert.util.SparseMat.sort_subs(posStrndsB, [this.sequenceLen this.nCompartments]);
            lensB = lensB(idxs);

            posStrnds = zeros(0, 2);
            lens = zeros(0, 1);
            for strand = 1:this.nCompartments
                rowsA = posStrndsA(:,2) == strand;
                rowsB = posStrndsB(:,2) == strand;
                posA = posStrndsA(rowsA, 1);
                posB = posStrndsB(rowsB, 1);
                lenA = lensA(rowsA, 1);
                lenB = lensB(rowsB, 1);
                iA = 1;
                iB = 1;
                while iA <= length(posA) && iB <= length(posB)
                    if posA(iA) <= posB(iB)
                        if posA(iA) + lenA(iA) > posB(iB)
                            posStrnds(end+1,:) = [posB(iB) strand];
                            if posA(iA) + lenA(iA) < posB(iB) + lenB(iB)
                                lens(end+1,:) = posA(iA) + lenA(iA) - posB(iB);
                                iA = iA + 1;
                            else
                                lens(end+1,:) = lenB(iB);
                                iB = iB + 1;
                            end
                        else
                            iA = iA + 1;
                        end
                    else
                        if posB(iB) + lenB(iB) > posA(iA)
                            posStrnds(end+1,:) = [posA(iA) strand];
                            if posB(iB) + lenB(iB) < posA(iA) + lenA(iA)
                                lens(end+1,:) = posB(iB) + lenB(iB) - posA(iA);
                                iB = iB + 1;
                            else
                                lens(end+1,:) = lenA(iA);
                                iA = iA + 1;
                            end
                        else
                            iB = iB + 1;
                        end
                    end
                end
            end

            [posStrnds, lens] = this.joinSplitOverOriCRegions(posStrnds, lens);
        end

        %split regions into several splitLen sized pieces
        function [posStrnds, lens] = splitRegions(~, rgnPosStrnds, rgnLens, splitLen)
            nPos = sum(floor(abs(rgnLens)/splitLen));
            posStrnds = zeros(nPos, 2);
            lens = repmat(splitLen, nPos, 1);
            
            j = 0;
            for i = 1:size(rgnPosStrnds, 1)
                nPos = floor(abs(rgnLens(i))/splitLen);
                
                posStrnds(j+(1:nPos), 1) = rgnPosStrnds(i, 1) + sign(rgnLens(i))*(0:nPos-1)*splitLen;
                posStrnds(j+(1:nPos), 2) = rgnPosStrnds(i, 2);
                
                j = j + nPos;
            end
        end

        function [posStrnds, lens] = splitOverOriC(this, posStrnds, lens)
            idxs = find(posStrnds(:, 1) + lens - 1 > this.sequenceLen);
            if isempty(idxs)
                return;
            end
            posStrnds = [posStrnds; ones(numel(idxs), 1) posStrnds(idxs, 2)];
            lens = [lens; posStrnds(idxs, 1) + reshape(lens(idxs), [], 1) - this.sequenceLen - 1];
            lens(idxs) = this.sequenceLen - posStrnds(idxs,1) + 1;
        end

        %join regions which have been split over the ORI
        function [posStrnds, lens] = joinSplitOverOriCRegions(this, posStrnds, lens)
            pos = posStrnds(:, 1);
            strnds = posStrnds(:, 2);
            ends = pos + lens - 1;
            
            for i = 1:max(strnds)
                idx1 = find(pos == 1 & strnds == i, 1, 'first');
                idx2 = find(ends == this.sequenceLen & strnds == i, 1, 'first');
                
                if isempty(idx1) || isempty(idx2) || idx1 == idx2
                    continue;
                end
                
                ends(idx2) = ends(idx2) + (ends(idx1) - pos(idx1) + 1);
                pos(idx1, :) = [];
                strnds(idx1, :) = [];
                ends(idx1, :) = [];
            end
            
            posStrnds = [pos strnds];
            lens = ends - pos + 1;
        end
        
        function [posStrnds, lens] = joinSplitRegions(this, posStrnds, lens)
            %sort
            posStrnds(:, 1) = mod(posStrnds(:, 1) - 1, this.sequenceLen) + 1;
            [posStrnds, order] = edu.stanford.covert.util.SparseMat.sort_subs(posStrnds, [this.sequenceLen this.nCompartments]);
            lens = lens(order);
            
            %join
            starts = posStrnds(:, 1);
            ends = starts + lens - 1;
            strnds = posStrnds(:, 2);
            
            tfs = true(size(starts));
            for i = 1:max(strnds)
                idxs = find(strnds == i);
                for j = 1:numel(idxs)-1
                    if ends(idxs(j))+1 >= starts(idxs(j+1))
                        starts(idxs(j+1)) = starts(idxs(j));
                        ends(idxs(j+1)) = max(ends(idxs(j)), ends(idxs(j+1)));
                        tfs(idxs(j)) = false;
                    end
                end
                if numel(idxs) >= 2
                    idx = idxs(find(tfs(idxs), 1, 'first'));
                    if ends(idxs(end))+1 >= starts(idx)+this.sequenceLen
                        ends(idxs(end)) = this.sequenceLen;
                        starts(idx) = 1;
                    end
                end
            end
            
            %format output
            posStrnds = [starts(tfs) strnds(tfs)];
            lens = ends(tfs) - starts(tfs) + 1;
        end
        
        function [idxs, newIdxs] = excludeOverlappingRegions(this, ...
                idxs, newIdxs, positionsStrands, lengths, ...
                footprint, footprint3Prime, footprint5Prime, isPositionsStrandFootprintCentroid, ...
                eitherStrand)
            
            tmpPositionsStrands = positionsStrands([idxs; newIdxs], :);
            tmpLengths = lengths([idxs; newIdxs], 1);
            
            startCoors = tmpPositionsStrands(:, 1) + min(0, tmpLengths + 1);
            if isPositionsStrandFootprintCentroid
                startCoors(mod(tmpPositionsStrands(:, 2), 2) == 1) = ...
                    startCoors(mod(tmpPositionsStrands(:, 2), 2) == 1) - footprint5Prime;
                startCoors(mod(tmpPositionsStrands(:, 2), 2) == 0) = ...
                    startCoors(mod(tmpPositionsStrands(:, 2), 2) == 0) - footprint3Prime;
            end
            endCoors = startCoors + (abs(tmpLengths) - 1) + (footprint - 1);
            if eitherStrand
                strnds = ceil(tmpPositionsStrands(:, 2) / 2);
            else
                strnds = tmpPositionsStrands(:, 2);
            end
            tmpIdxs = [zeros(size(idxs)); (1:numel(newIdxs))'];
            
            tmp = find(startCoors < 0);
            startCoors = [startCoors; startCoors(tmp, :) + this.sequenceLen];
            endCoors   = [endCoors; min(this.sequenceLen, endCoors(tmp, :) + this.sequenceLen)];
            strnds       = [strnds; strnds(tmp, :)];
            tmpIdxs    = [tmpIdxs; tmpIdxs(tmp, :)];
            
            tmp = find(endCoors > this.sequenceLen);
            startCoors = [startCoors; max(1, startCoors(tmp, :) - this.sequenceLen)];
            endCoors   = [endCoors; endCoors(tmp, :) - this.sequenceLen];
            strnds       = [strnds; strnds(tmp, :)];
            tmpIdxs    = [tmpIdxs; tmpIdxs(tmp, :)];
            
            tmpTfs = true(size(newIdxs));
            for i = numel(tmpIdxs):-1:2
                if tmpIdxs(i) <= numel(idxs)
                    continue; 
                end
                tmpTfs(tmpIdxs(i)) = tmpTfs(tmpIdxs(i)) && ~any((...
                    startCoors(1:i-1) <= startCoors(i) & startCoors(i) <= endCoors(1:i-1) | ...
                    startCoors(1:i-1) <= endCoors(i)   & endCoors(  i) <= endCoors(1:i-1)) & ...
                    strnds(i) == strnds(1:i-1));
            end
            newIdxs = newIdxs(tmpTfs);
        end
    end
    
    %setters
    methods        
        %integers [positions x strands] indicating the start positions of
        %polymerized regions of strands and their lengths
        function set.polymerizedRegions(this, value)
            if isequal(this.polymerizedRegions, value)
                return;
            end
            this.polymerizedRegions = value;
            this.validated = this.validated  + 1;
            this.validated_polymerizedRegions = this.validated; %#ok<*MCSUP>
        end
        
        %integers [positions x strands] indicating the current linking number of
        %each double-stranded region
        function set.linkingNumbers(this, value)
            %NOTE: performance likely better here without checking if new value
            %is different than old
            this.linkingNumbers = value;
            this.validated = this.validated  + 1;
            this.validated_linkingNumbers = this.validated;
        end
        
        %indices [positions x strands] indicating start positions of protein
        %monomers bound to DNA bases
        function set.monomerBoundSites(this, value)
            %NOTE: performance likely better here without checking if new value
            %is different than old
            this.monomerBoundSites = value;
            this.validated = this.validated  + 1;
            this.validated_proteinBoundSites = this.validated;
        end
        
        %indices [positions x strands] indicating start positions of
        %macromolecular complexes bound to DNA bases
        function set.complexBoundSites(this, value)
            %NOTE: performance likely better here without checking if new value
            %is different than old
            this.complexBoundSites = value;
            this.validated = this.validated  + 1;
            this.validated_proteinBoundSites = this.validated;
        end
        
        %boolean [positions x strands] indicating positions of gap sites
        function set.gapSites(this, value)
            if isequal(this.gapSites, value)
                return;
            end
            this.gapSites = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_gapSites = this.validated;
        end
        
        %boolean [positions x strands] indicating positions of abasic sites
        function set.abasicSites(this, value)
            if isequal(this.abasicSites, value)
                return;
            end
            this.abasicSites = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_abasicSites = this.validated;
        end
        
        %indices [positions x strands] indicating metabolite identity of damaged
        %sugar-phosphates
        function set.damagedSugarPhosphates(this, value)
            if isequal(this.damagedSugarPhosphates, value)
                return;
            end
            this.damagedSugarPhosphates = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_damagedSugarPhosphates = this.validated;
        end
        
        %indices [positions x strands] indicating metabolite identity of damaged
        %bases
        function set.damagedBases(this, value)
            if isequal(this.damagedBases, value)
                return;
            end
            this.damagedBases = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_damagedBases = this.validated;
        end
        
        %boolean [positions x strands] indicating metabolite identity of
        %intrastrand cross links in DNA
        function set.intrastrandCrossLinks(this, value)
            if isequal(this.intrastrandCrossLinks, value)
                return;
            end
            this.intrastrandCrossLinks = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_intrastrandCrossLinks = this.validated;
        end
        
        %boolean [positions x strands] indicating positions of strand breaks in
        %strands of DNA
        function set.strandBreaks(this, value)
            if isequal(this.strandBreaks, value)
                return;
            end
            this.strandBreaks = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_strandBreaks = this.validated;
        end
        
        %boolean [positions x strands] indicating positions of holliday
        %junctions
        function set.hollidayJunctions(this, value)
            if isequal(this.hollidayJunctions, value)
                return;
            end
            this.hollidayJunctions = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_hollidayJunctions = this.validated;
        end
        
        %boolean indicating whether or not the chromsomes are segregated
        function set.segregated(this, value)
            if isequal(this.segregated, value)
                return;
            end
            this.segregated = value;
            this.validated = this.validated  + 1;
            this.validated_segregated = this.validated;
        end
    end
    
    %getters for alternative views of state
    methods
        function invalidate(this)
            this.validated = uint32(1);
            
            this.validated_polymerizedRegions     = this.validated;
            this.validated_linkingNumbers         = this.validated;
            this.validated_proteinBoundSites      = this.validated;
            this.validated_damaged                = this.validated;
            this.validated_abasicSites            = this.validated;
            this.validated_gapSites               = this.validated;
            this.validated_damagedSugarPhosphates = this.validated;
            this.validated_damagedBases           = this.validated;
            this.validated_strandBreaks           = this.validated;
            this.validated_intrastrandCrossLinks  = this.validated;
            this.validated_hollidayJunctions      = this.validated;
            this.validated_segregated             = this.validated;
            
            this.validated_unpolymerizedRegions          = uint32(0);
            this.validated_singleStrandedRegions         = uint32(0);
            this.validated_doubleStrandedRegions         = uint32(0);
            this.validated_geneCopyNumbers               = uint32(0);
            this.validated_ploidy                        = uint32(0);
            this.validated_polymerizedGenes              = uint32(0);
            this.validated_transcriptionUnitCopyNumbers  = uint32(0);
            this.validated_polymerizedTranscriptionUnits = uint32(0);
            this.validated_geneCopyNumbers_Accessible    = uint32(0);
            this.validated_transcriptionUnitCopyNumbers_Accessible = uint32(0);
            this.validated_accessibleGenes               = uint32(0);
            this.validated_accessibleTranscriptionUnits  = uint32(0);            
            this.validated_linkingNumbers_minFreeEnergy  = uint32(0);
            this.validated_supercoils                    = uint32(0);
            this.validated_superhelicalDensity           = uint32(0);
            this.validated_supercoiled                   = uint32(0);
            this.validated_damagedSites                  = uint32(0);
            this.validated_damagedSites_shifted_incm6AD  = uint32(0);
            this.validated_damagedSites_nonRedundant     = uint32(0);
            this.validated_damagedSites_excm6AD          = uint32(0);
            this.validated_gapSites3                     = uint32(0);
            this.validated_gapSites5                     = uint32(0);
            this.validated_abasicSites3                  = uint32(0);
            this.validated_abasicSites5                  = uint32(0);
            this.validated_damagedSugarPhosphates3       = uint32(0);
            this.validated_damagedSugarPhosphates5       = uint32(0);
            this.validated_damagedBases3                 = uint32(0);
            this.validated_damagedBases5                 = uint32(0);
            this.validated_strandBreaks3                 = uint32(0);
            this.validated_strandBreaks5                 = uint32(0);
            this.validated_intrastrandCrossLinks3        = uint32(0);
            this.validated_intrastrandCrossLinks5        = uint32(0);
            this.validated_hollidayJunctions3            = uint32(0);
            this.validated_hollidayJunctions5            = uint32(0);
            this.validated_singleStrandBreaks            = uint32(0);
            this.validated_doubleStrandBreaks            = uint32(0);
            this.validated_strandBreakClassification     = uint32(0);
            this.validated_munIRMSiteMethylationStatus   = uint32(0);
            this.validated_munIRMSiteRestrictionStatus   = uint32(0);
            this.validated_dryWeight                     = uint32(0);
        end
        
        %integers indicating the start positions of unpolymerized regions (ie.
        %not yet replicated) of strands and their lengths
        function value = get.unpolymerizedRegions(this)
            if this.validated_polymerizedRegions > this.validated_unpolymerizedRegions
                this.unpolymerizedRegions = this.calcUnpolymerizedRegions();
                this.validated_unpolymerizedRegions = this.validated;
            end
                
            value = this.unpolymerizedRegions;
        end
                
        function value = calcUnpolymerizedRegions(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            [polPosStrndsTimes, polLens] = find(this.polymerizedRegions);
            polPos = polPosStrndsTimes(:, 1);
            polStrnd = polPosStrndsTimes(:, 2);
            polTimes = polPosStrndsTimes(:, 3:end);
            if isempty(polTimes)
                polTimes = ones(size(polPos));
            end
            
            unpolPos = [];
            unpolStrnds = [];
            unpolTimes = [];
            unpolLens = [];
            
            for i = 1:this.nCompartments
                idxs = find(polStrnd == i);
                
                for j = 1:size(this.polymerizedRegions, 3)
                    idxs2 = idxs(polTimes(idxs) == j);
                    if isempty(idxs2)
                        unpolPos = [unpolPos; 1];
                        unpolLens = [unpolLens; this.sequenceLen];
                        unpolStrnds = [unpolStrnds; i];
                        unpolTimes = [unpolTimes; j];
                        continue;
                    end
                    
                    unpolPos = [
                        unpolPos;
                        1;
                        polPos(idxs2) + polLens(idxs2)]; %#ok<*AGROW>
                    unpolLens = [
                        unpolLens;
                        polPos(idxs2(1))-1;
                        polPos(idxs2(2:end)) - (polPos(idxs2(1:end-1)) + polLens(idxs2(1:end-1)));
                        this.sequenceLen - (polPos(idxs2(end)) + polLens(idxs2(end))) + 1;
                        ]; %#ok<*AGROW>
                    unpolStrnds = [
                        unpolStrnds;
                        i(ones(numel(idxs2) + 1, 1), 1)]; %#ok<*AGROW>
                    unpolTimes = [
                        unpolTimes;
                        j(ones(numel(idxs2)+1, 1), 1)]; %#ok<*AGROW>
                end
            end
            
            idxs = find(unpolLens > 0);
            value = CircularSparseMat([unpolPos(idxs) unpolStrnds(idxs) unpolTimes(idxs)], unpolLens(idxs), size(this.polymerizedRegions), 1);            
        end
                
        function value = get.singleStrandedRegions(this)
            if this.validated_polymerizedRegions > this.validated_singleStrandedRegions
                this.singleStrandedRegions = this.calcSingleStrandedRegions();
                this.validated_singleStrandedRegions = this.validated;
            end
                
            value = this.singleStrandedRegions;
        end
        
        function value = calcSingleStrandedRegions(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            [positionsStrandsTimes, lengths] = find(this.polymerizedRegions);
            strnd = positionsStrandsTimes(:, 2);
            oppStrnd = strnd;
            oppStrnd(mod(strnd, 2) == 0, 1) = strnd(mod(strnd, 2) == 0, 1) - 1;
            oppStrnd(mod(strnd, 2) == 1, 1) = strnd(mod(strnd, 2) == 1, 1) + 1;
            
            starts = mod([
                positionsStrandsTimes(:,1);
                positionsStrandsTimes(:,1) + lengths] ...
                - 1, this.sequenceLen) + 1;
            strandsTimes = [
                strnd positionsStrandsTimes(:, 3:end);
                oppStrnd positionsStrandsTimes(:, 3:end)];
            oppStrandsTimes = [
                oppStrnd positionsStrandsTimes(:, 3:end);
                strnd positionsStrandsTimes(:, 3:end)];
            lengths = [
                lengths;
                repmat(max(lengths), size(lengths))];
            
            [~, ~, ~, extents1] = this.isRegionPolymerized([starts strandsTimes], lengths, true);
            [~, ~, ~, extents2] = this.isRegionNotPolymerized([starts oppStrandsTimes], lengths, true);
            extents = min(abs(extents1), abs(extents2));
            idxs = find(extents > 0);
            
            if size(strandsTimes, 2) == 1
                tmp = edu.stanford.covert.util.SparseMat.unique_subs(...
                    [starts(idxs, :) strandsTimes(idxs, :) extents(idxs)], ...
                    [this.sequenceLen this.nCompartments this.sequenceLen]);
            else
                tmp = edu.stanford.covert.util.SparseMat.unique_subs(...
                    [starts(idxs, :) strandsTimes(idxs, :) extents(idxs)], ...
                    [this.sequenceLen this.nCompartments size(this.polymerizedRegions, 3) this.sequenceLen]);
            end
            
            for i = size(tmp, 1):-1:1
                if tmp(i, 1) == 1
                    idx = find(tmp(:, 1) + tmp(:, 3) - 1 > this.sequenceLen & tmp(:, 2) == tmp(i, 2));                    
                    if ~isempty(idx)
                        tmp(idx, 3) = max(tmp(idx, 3), (this.sequenceLen - tmp(idx,1) + 1) + tmp(i, 3));
                        tmp(i, :) = [];
                    end
                end
            end
            
            posStrnds = tmp(:, 1:end-1);
            lens = tmp(:, end);
            [posStrnds, lens] = this.splitOverOriC(posStrnds, lens);
            value = this.mergeAdjacentRegions(CircularSparseMat(posStrnds, lens, size(this.polymerizedRegions), 1));
        end
                
        function value = get.doubleStrandedRegions(this)
            if this.validated_polymerizedRegions > this.validated_doubleStrandedRegions
                this.doubleStrandedRegions = this.calcDoubleStrandedRegions();
                this.validated_doubleStrandedRegions = this.validated;
            end
                
            value = this.doubleStrandedRegions;
        end
        
        function value = calcDoubleStrandedRegions(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            [positionsStrandsTimes, lengths] = find(this.polymerizedRegions);
            
                        
            if size(positionsStrandsTimes, 2) == 2
                tmp = edu.stanford.covert.util.SparseMat.unique_subs([
                    positionsStrandsTimes(:, 1) ceil(positionsStrandsTimes(:, 2)/2) lengths;
                    ], [this.sequenceLen this.nCompartments this.sequenceLen]);
                
                [~, ~, ~, extents1] = this.isRegionPolymerized([tmp(:, 1) 2*tmp(:,2)-1], tmp(:, end), true);
                [~, ~, ~, extents2] = this.isRegionPolymerized([tmp(:, 1) 2*tmp(:,2)  ], tmp(:, end), true);
                extents = min(abs(extents1), abs(extents2));
                idxs = find(extents > 0);
                
                tmp = edu.stanford.covert.util.SparseMat.unique_subs([
                    tmp(idxs, 1) 2*tmp(idxs, 2)-1 extents(idxs)
                    tmp(idxs, 1) 2*tmp(idxs, 2)   extents(idxs)
                    ], [this.sequenceLen this.nCompartments this.sequenceLen]);
            else
                tmp = edu.stanford.covert.util.SparseMat.unique_subs([
                    positionsStrandsTimes(:, 1) ceil(positionsStrandsTimes(:, 2)/2) positionsStrandsTimes(:, 3:end) lengths;
                    ], [this.sequenceLen this.nCompartments size(this.polymerizedRegions, 3) this.sequenceLen]);
                
                [~, ~, ~, extents1] = this.isRegionPolymerized([tmp(:, 1) 2*tmp(:,2)-1 tmp(:,3:end-1)], tmp(:, end), true);
                [~, ~, ~, extents2] = this.isRegionPolymerized([tmp(:, 1) 2*tmp(:,2)   tmp(:,3:end-1)], tmp(:, end), true);
                extents = min(abs(extents1), abs(extents2));
                idxs = find(extents > 0);
                
                tmp = edu.stanford.covert.util.SparseMat.unique_subs([
                    tmp(idxs, 1) 2*tmp(idxs, 2)-1 tmp(idxs, 3:end-1) extents(idxs)
                    tmp(idxs, 1) 2*tmp(idxs, 2)   tmp(idxs, 3:end-1) extents(idxs)
                    ], [this.sequenceLen this.nCompartments size(this.polymerizedRegions, 3) this.sequenceLen]);
            end
            
            for i = size(tmp, 1):-1:1
                if tmp(i, 1) == 1
                    idx = find(tmp(:, 1) + tmp(:, 3) - 1 > this.sequenceLen & tmp(:, 2) == tmp(i, 2));
                    if ~isempty(idx)
                        tmp(idx, 3) = max(tmp(idx, 3), (this.sequenceLen - tmp(idx,1) + 1) + tmp(i, 3));
                        tmp(i, :) = [];
                    end
                end
            end
            
            posStrnds = tmp(:, 1:end-1);
            lens = tmp(:, end);
            [posStrnds, lens] = this.splitOverOriC(posStrnds, lens);
            value = this.mergeAdjacentRegions(CircularSparseMat(posStrnds, lens, size(this.polymerizedRegions), 1));
        end
        
        
        function value = get.geneCopyNumbers(this)
            if this.validated_polymerizedRegions > this.validated_geneCopyNumbers
                this.geneCopyNumbers = this.calcGeneCopyNumbers();
                this.validated_geneCopyNumbers = this.validated;
            end
                
            value = this.geneCopyNumbers;
        end
        
        %number of copies of each gene that have been polymerized (Nx1)
        function value = calcGeneCopyNumbers(this)
            value = sum(this.polymerizedGenes, 2);
        end
        
        function value = get.ploidy(this)
            if this.validated_polymerizedRegions > this.validated_ploidy
                this.ploidy = this.calcPloidy();
                this.validated_ploidy = this.validated;
            end
            
            value = this.ploidy;
        end
        
        function value = calcPloidy(this)
            value = collapse(this.polymerizedRegions)/(2*this.sequenceLen);
        end
        
        function value = get.polymerizedGenes(this)
            if this.validated_polymerizedRegions > this.validated_polymerizedGenes
                this.polymerizedGenes = this.calcPolymerizedGenes();
                this.validated_polymerizedGenes = this.validated;
            end
                
            value = this.polymerizedGenes;
        end

        %whether each copy of each gene has been polymerized (Nx2)
        function value = calcPolymerizedGenes(this)
            value = this.isRegionPolymerized(...
                [this.gene.startCoordinates this.gene.strands;
                 this.gene.startCoordinates this.gene.strands+2],...
                [this.gene.lengths; this.gene.lengths], false);
            value = reshape(value, [], 2);
        end
        
        
        function value = get.transcriptionUnitCopyNumbers(this)
            if this.validated_polymerizedRegions > this.validated_transcriptionUnitCopyNumbers
                this.transcriptionUnitCopyNumbers = this.getTranscriptionUnitCopyNumbers();
                this.validated_transcriptionUnitCopyNumbers = this.validated;
            end
                
            value = this.transcriptionUnitCopyNumbers;
        end

        %number of copies of each transcription unit that have been polymerized (Nx1)
        function value = getTranscriptionUnitCopyNumbers(this)
            value = sum(this.polymerizedTranscriptionUnits, 2);
        end

        
        function value = get.polymerizedTranscriptionUnits(this)
            if this.validated_polymerizedRegions > this.validated_polymerizedTranscriptionUnits
                this.polymerizedTranscriptionUnits = this.calcPolymerizedTranscriptionUnits();
                this.validated_polymerizedTranscriptionUnits = this.validated;
            end
                
            value = this.polymerizedTranscriptionUnits;
        end
        
        %whether each copy of each transcription unit has been polymerized (Nx2)
        function value = calcPolymerizedTranscriptionUnits(this)
            value = this.isRegionPolymerized(...
                [this.transcriptionUnitStartCoordinates this.transcriptionUnitStrands;
                 this.transcriptionUnitStartCoordinates this.transcriptionUnitStrands+2],...
                [this.transcriptionUnitLengths; this.transcriptionUnitLengths], false);
            value = reshape(value, [], 2);
        end
        
        
        function value = get.geneCopyNumbers_Accessible(this)
            val = max([
                this.validated_polymerizedRegions; 
                this.validated_damaged;
                this.validated_proteinBoundSites]);
            if val > this.validated_geneCopyNumbers_Accessible
                this.geneCopyNumbers_Accessible = this.calcCopyNumbers_Accessible();
                this.validated_geneCopyNumbers_Accessible = this.validated;
            end
                
            value = this.geneCopyNumbers_Accessible;
        end

        %number of copies of each gene that are accessible
        function value = calcCopyNumbers_Accessible(this)
            value = sum(this.accessibleGenes, 2);
        end
        
        
        function value = get.transcriptionUnitCopyNumbers_Accessible(this)
            val = max([
                this.validated_polymerizedRegions; 
                this.validated_damaged;
                this.validated_proteinBoundSites]);
            if val > this.validated_transcriptionUnitCopyNumbers_Accessible
                this.transcriptionUnitCopyNumbers_Accessible = this.calcTranscriptionUnitCopyNumbers_Accessible();
                this.validated_transcriptionUnitCopyNumbers_Accessible = this.validated;
            end
                
            value = this.transcriptionUnitCopyNumbers_Accessible;
        end
                
        %number of copies of each transcription unit that are accessible
        function value = calcTranscriptionUnitCopyNumbers_Accessible(this)
            value = sum(this.accessibleTranscriptionUnits, 2);
        end
        
        
        function value = get.accessibleGenes(this)
            val = max([
                this.validated_polymerizedRegions; 
                this.validated_damaged;
                this.validated_proteinBoundSites]);
            if val > this.validated_accessibleGenes
                this.accessibleGenes = this.calcAccessibleGenes();
                this.validated_accessibleGenes = this.validated;
            end
                
            value = this.accessibleGenes;
        end
        
        %boolean indicator of undamaged, unoccupied genes
        %true  ==> gene is accessible
        %false ==> gene is inaccessible
        function value = calcAccessibleGenes(this)
            value = reshape(this.isRegionAccessible([...
                this.gene.startCoordinates this.gene.strands
                this.gene.startCoordinates this.gene.strands+2], ...
                [this.gene.lengths; this.gene.lengths], [], [], true, [], false, true), [], this.nCompartments/2);
        end
        
        
        function value = get.accessibleTranscriptionUnits(this)
            val = max([
                this.validated_polymerizedRegions; 
                this.validated_damaged;
                this.validated_proteinBoundSites]);
            if val > this.validated_accessibleTranscriptionUnits
                this.accessibleTranscriptionUnits = this.calcAccessibleTranscriptionUnits();
                this.validated_accessibleTranscriptionUnits = this.validated;
            end
                
            value = this.accessibleTranscriptionUnits;
        end
        
        %boolean indicator of undamaged, unoccupied transcription units
        %true  ==> transcription unit is accessible
        %false ==> transcription unit is inaccessible
        function value = calcAccessibleTranscriptionUnits(this)
            value = reshape(this.isRegionAccessible([...
                this.transcriptionUnitStartCoordinates this.transcriptionUnitStrands;
                this.transcriptionUnitStartCoordinates this.transcriptionUnitStrands+2], ...
                [this.transcriptionUnitLengths; this.transcriptionUnitLengths], [], [], true, [], false, true), ...
                [], this.nCompartments/2);
        end
        
        function value = get.strandBreakClassification(this)
            if max(this.validated_polymerizedRegions, this.validated_damaged) > this.validated_strandBreakClassification
                this.strandBreakClassification = this.calcStrandBreakClassification();
                this.validated_strandBreakClassification = this.validated;
            end
                
            value = this.strandBreakClassification;
        end
        
        %numbers of each class of strand break (SSB, SSB+, 2SSB, DSB, DSB+,
        %DSB++) in DNA
        function value = calcStrandBreakClassification(this)
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.util.countUnique;
            
            %parameters
            segmentLength = this.strandBreakClassification_segmentLength;
            dsbSep = this.strandBreakClassification_doubleStrandBreakSeparation;
            genomeLength = this.sequenceLen;
            numStrands = this.nCompartments;
            numTime = size(this.strandBreaks, 3);
            numSegments = ceil(genomeLength / segmentLength);
            
            %damaged Sites
            [polymerizedPositionsStrands, polymerizedLengths] = find(this.polymerizedRegions);
            idxs = find(polymerizedPositionsStrands(:,1) <= segmentLength * numSegments - genomeLength);
            polymerizedPositionsStrands = [polymerizedPositionsStrands; polymerizedPositionsStrands(idxs, 1) + this.sequenceLen polymerizedPositionsStrands(idxs, 2:end)];
            polymerizedLengths = [polymerizedLengths; polymerizedLengths(idxs, :)];
            idxs = find(polymerizedPositionsStrands(:,1) + polymerizedLengths - 1 > segmentLength * numSegments);
            polymerizedLengths(idxs,:) = segmentLength * numSegments - polymerizedPositionsStrands(idxs, 1) + 1;
            [polymerizedPositionsStrands, polymerizedLengths] = find(this.mergeAdjacentRegions(...
                CircularSparseMat(polymerizedPositionsStrands, polymerizedLengths, [segmentLength * numSegments numStrands numTime], 1)));
            
            if numTime == 1
                polymerizedPositionsStrands = [polymerizedPositionsStrands ones(size(polymerizedPositionsStrands,1), 1)];
            end
            
            polymerizedStrands = [];
            for i = 1:size(polymerizedPositionsStrands, 1)
                regions = ceil((polymerizedPositionsStrands(i,1)-1) / segmentLength)+1 : floor((polymerizedPositionsStrands(i,1) + polymerizedLengths(i,1) -1) / segmentLength);
                polymerizedStrands = [polymerizedStrands;
                    regions' ...
                    repmat(ceil(polymerizedPositionsStrands(i,2)/2), numel(regions), 1) ...
                    repmat(polymerizedPositionsStrands(i,3), numel(regions), 1)];
            end
            [polymerizedRegions, ~, idxs]= unique(polymerizedStrands, 'rows'); %#ok<PROP>
            [idxs, counts] = countUnique(idxs);
            polymerizeTimes = polymerizedRegions(idxs(counts == 2), 3:end); %#ok<PROP>
            
            numPolymerizedSegments = zeros(1, 1, numTime);
            [idxs, counts] = countUnique(polymerizeTimes);
            numPolymerizedSegments(idxs) = counts;
            
            %initialize classification
            dmgSites = [this.strandBreaks; this.strandBreaks(1:segmentLength * numSegments - genomeLength, :, :)];
            
            subs = find(permute(reshape(dmgSites, [], numSegments, numStrands, numTime), [3 1 2 4])); %[strands X positions X segments X time]
            if numTime == 1; subs = [subs ones(size(subs,1), 1)]; end
            segmentTimeInds = sub2ind([numSegments numTime], subs(:, 3), subs(:, 4));
            
            value = zeros(this.strandBreakClassification_index_DSB__, 1, numTime);
            
            %classify segments
            while ~isempty(subs)
                %time
                time = subs(1, 4);
                
                %find extent of segment
                endIdx = find(...
                    segmentTimeInds(1) ~= segmentTimeInds | ...
                    ceil(subs(1,1)/2) ~= ceil(subs(:,1)/2), ...
                    1, 'first') - 1;
                if isempty(endIdx)
                    endIdx = size(subs, 1);
                end
                
                if ~this.isRegionDoubleStranded([(subs(:,3)-1)*segmentLength+1 subs(:,1) subs(:,4)], segmentLength, false);
                    continue;
                end
                
                %classify damaged segment
                if endIdx == 1
                    classification = this.strandBreakClassification_index_SSB;
                else
                    strands   = subs(1:endIdx, 1);
                    positions = subs(1:endIdx, 2);
                    if ~isempty(strands) && all(strands == strands(1))
                        classification = this.strandBreakClassification_index_SSB_;
                    else
                        dsbIdx = find(diff(positions) < dsbSep & diff(strands));
                        if isempty(dsbIdx)
                            classification = this.strandBreakClassification_index_2SSB;
                        else
                            if numel(dsbIdx)>1 && diff(positions(dsbIdx([1 end]))) >= dsbSep
                                classification = this.strandBreakClassification_index_DSB__;
                            elseif numel(dsbIdx)>1 || ...
                                    (dsbIdx>1 && diff(positions([dsbIdx - 1 dsbIdx])) < dsbSep) || ...
                                    (dsbIdx<numel(positions) - 1 && diff(positions([dsbIdx + 1 dsbIdx + 2])) < dsbSep)
                                classification = this.strandBreakClassification_index_DSB_;
                            else
                                classification = this.strandBreakClassification_index_DSB;
                            end
                        end
                    end
                end
                
                %update counts of classified segments
                value(classification, 1, time) = value(classification, 1, time) + 1;
                
                %shrink subs
                subs(1:endIdx, :) = [];
                segmentTimeInds(1:endIdx, :) = [];
            end
            
            %compute numbers of segments without damage
            value(this.strandBreakClassification_index_NB, 1, :) =  numPolymerizedSegments - sum(value,1);
        end
        
        %
        function value = get.linkingNumbers_minFreeEnergy(this)
            if this.validated_polymerizedRegions > this.validated_linkingNumbers_minFreeEnergy
                this.linkingNumbers_minFreeEnergy = this.calcLinkingNumbers_minFreeEnergy();
                this.validated_linkingNumbers_minFreeEnergy = this.validated;
            end
            
            value = this.linkingNumbers_minFreeEnergy;
        end
        
        function value = calcLinkingNumbers_minFreeEnergy(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            [posStrands, lens] = find(this.doubleStrandedRegions);
            
            value = CircularSparseMat(posStrands, lens / this.relaxedBasesPerTurn, size(this.linkingNumbers), 1);
        end
        
        %
        function value = get.supercoils(this)
            if any([this.validated_polymerizedRegions; this.validated_linkingNumbers] > this.validated_supercoils)
                this.supercoils = this.calcSupercoils();
                this.validated_supercoils = this.validated;
            end
            
            value = this.supercoils;
        end
        
        function value = calcSupercoils(this)
            value = this.linkingNumbers - this.linkingNumbers_minFreeEnergy;
        end
        
        %
        function value = get.superhelicalDensity(this)
            if any([this.validated_polymerizedRegions; this.validated_linkingNumbers] > this.validated_superhelicalDensity)
                this.superhelicalDensity = this.calcSuperhelicalDensity();
                this.validated_superhelicalDensity = this.validated;
            end
            
            value = this.superhelicalDensity;
        end
        
        function value = calcSuperhelicalDensity(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            lkNums_min = this.linkingNumbers_minFreeEnergy;
            [posStrnds, deltas] = find(this.linkingNumbers - lkNums_min);
            value = CircularSparseMat(posStrnds, deltas ./ lkNums_min(posStrnds), size(this.linkingNumbers), 1);
        end
        
        %check if superhelical density within tolerance of equilbrium value
        function value = get.supercoiled(this)
            if any([this.validated_polymerizedRegions; this.validated_linkingNumbers] > this.validated_supercoiled)
                this.supercoiled = this.calcSupercoiled();
                this.validated_supercoiled = this.validated;
            end
            
            value = this.supercoiled;
        end
        
        function value = calcSupercoiled(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            siz = [this.sequenceLen this.nCompartments];
            
            [posStrnds, lens] = find(this.doubleStrandedRegions);
            lks = this.linkingNumbers(posStrnds);
            
            lk0s = lens / this.relaxedBasesPerTurn;
            sigmas = (lks - lk0s) ./ lk0s;
            tfs = abs(sigmas - this.equilibriumSuperhelicalDensity) < this.supercoiledSuperhelicalDensityTolerance;
            
            value = CircularSparseMat(posStrnds, tfs, siz, 1);
        end
        
        function value = get.damagedSites(this)
            if this.validated_damaged > this.validated_damagedSites
                this.damagedSites = this.calcDamagedSites();
                this.validated_damagedSites = this.validated;
            end
            
            value = this.damagedSites;
        end
        
        function value = calcDamagedSites(this)
            value = this.getDamagedSites(true, true, true, true, false, true, false);
        end
        
        function value = get.damagedSites_shifted_incm6AD(this)
            if this.validated_damaged > this.validated_damagedSites_shifted_incm6AD
                this.damagedSites_shifted_incm6AD = this.calcDamagedSites_shifted_incm6AD();
                this.validated_damagedSites_shifted_incm6AD = this.validated;
            end
            
            value = this.damagedSites_shifted_incm6AD;
        end
        
        function value = calcDamagedSites_shifted_incm6AD(this)
            value = this.getDamagedSites(true, true, true, true, false, true, true);
        end
        
        function value = get.damagedSites_nonRedundant(this)
            if this.validated_damaged > this.validated_damagedSites_nonRedundant
                this.damagedSites_nonRedundant = this.calcDamagedSites_nonRedundant();
                this.validated_damagedSites_nonRedundant = this.validated;
            end
            
            value = this.damagedSites_nonRedundant;
        end
        
        function value = calcDamagedSites_nonRedundant(this)
            value = this.getDamagedSites(true, true, false, false, false, false, true);
        end
        
        function value = get.damagedSites_excm6AD(this)
            if this.validated_damaged > this.validated_damagedSites_excm6AD
                this.damagedSites_excm6AD = this.calcDamagedSites_excm6AD();
                this.validated_damagedSites_excm6AD = this.validated;
            end
            
            value = this.damagedSites_excm6AD;
        end
        
        function value = calcDamagedSites_excm6AD(this)
            value = this.getDamagedSites(true, true, false, false, false, false, false);
        end
        
        %boolean (genome length x 2) indicating positions of gap sites 3'
        %to bases
        function value = get.gapSites3(this)
            if this.validated_gapSites > this.validated_gapSites3
                this.gapSites3 = this.calcGapSites3();
                this.validated_gapSites3 = this.validated;
            end
            
            value = this.gapSites3;
        end
        
        function value = calcGapSites3(this)
            value = this.shiftCircularSparseMatBase5Prime(this.gapSites);
        end
        
        %boolean (genome length x 2) indicating positions of gap sites 5'
        %to bases
        function value = get.gapSites5(this)
            if this.validated_gapSites > this.validated_gapSites5
                this.gapSites5 = this.calcGapSites5();
                this.validated_gapSites5 = this.validated;
            end
            
            value = this.gapSites5;
        end
        
        function value = calcGapSites5(this)
            value = this.shiftCircularSparseMatBase3Prime(this.gapSites);
        end
        
        %boolean (genome length x 2) indicating positions of abasic sites 3'
        %to bases
        function value = get.abasicSites3(this)
            if this.validated_abasicSites > this.validated_abasicSites3
                this.abasicSites3 = this.calcAbasicSites3();
                this.validated_abasicSites3 = this.validated;
            end
            
            value = this.abasicSites3;
        end
        
        function value = calcAbasicSites3(this)
            value = this.shiftCircularSparseMatBase5Prime(this.abasicSites);
        end
        
        %boolean (genome length x 2) indicating positions of abasic sites 5'
        %to bases
        function value = get.abasicSites5(this)
            if this.validated_abasicSites > this.validated_abasicSites5
                this.abasicSites5 = this.calcAbasicSites5();
                this.validated_abasicSites5 = this.validated;
            end
            
            value = this.abasicSites5;
        end
        
        function value = calcAbasicSites5(this)
            value = this.shiftCircularSparseMatBase3Prime(this.abasicSites);
        end
        
        %integer (genome length x 2) indicating indices of damaged sugar
        %phosphates 3' to bases
        function value = get.damagedSugarPhosphates3(this)
            if this.validated_damagedSugarPhosphates > this.validated_damagedSugarPhosphates3
                this.damagedSugarPhosphates3 = this.calcDamagedSugarPhosphates3();
                this.validated_damagedSugarPhosphates3 = this.validated;
            end
            
            value = this.damagedSugarPhosphates3;
        end
        
        function value = calcDamagedSugarPhosphates3(this)
            value = this.shiftCircularSparseMatBase5Prime(this.damagedSugarPhosphates);
        end
        
        %integer (genome length x 2) indicating indices of damaged sugar
        %phosphates 5' to bases
        function value = get.damagedSugarPhosphates5(this)
            if this.validated_damagedSugarPhosphates > this.validated_damagedSugarPhosphates5
                this.damagedSugarPhosphates5 = this.calcDamagedSugarPhosphates5();
                this.validated_damagedSugarPhosphates5 = this.validated;
            end
            
            value = this.damagedSugarPhosphates5;
        end
        
        function value = calcDamagedSugarPhosphates5(this)
            value = this.shiftCircularSparseMatBase3Prime(this.damagedSugarPhosphates);
        end
        
        %integer (genome length x 2) indicating indices of damaged bases 3' to bases
        function value = get.damagedBases3(this)
            if this.validated_damagedBases > this.validated_damagedBases3
                this.damagedBases3 = this.calcDamagedBases3();
                this.validated_damagedBases3 = this.validated;
            end
            
            value = this.damagedBases3;
        end
        
        function value = calcDamagedBases3(this)
            value = this.shiftCircularSparseMatBase5Prime(this.damagedBases);
        end
        
        %integer (genome length x 2) indicating indices of damaged bases 5' to bases
        function value = get.damagedBases5(this)
            if this.validated_damagedBases > this.validated_damagedBases5
                this.damagedBases5 = this.calcDamagedBases5();
                this.validated_damagedBases5 = this.validated;
            end
            
            value = this.damagedBases5;
        end
        
        function value = calcDamagedBases5(this)
            value = this.shiftCircularSparseMatBase3Prime(this.damagedBases);
        end
        
        %boolean (genome length x 2) indicating positions of intrastrand cross
        %links 3' to bases
        function value = get.intrastrandCrossLinks3(this)
            if this.validated_intrastrandCrossLinks > this.validated_intrastrandCrossLinks3
                this.intrastrandCrossLinks3 = this.calcIntrastrandCrossLinks3();
                this.validated_intrastrandCrossLinks3 = this.validated;
            end
            
            value = this.intrastrandCrossLinks3;
        end
        
        function value = calcIntrastrandCrossLinks3(this)
            value = this.shiftCircularSparseMatBase5Prime(this.intrastrandCrossLinks);
        end
        
        %boolean (genome length x 2) indicating positions of intrastrand cross
        %links 5' to bases
        function value = get.intrastrandCrossLinks5(this)
            if this.validated_intrastrandCrossLinks > this.validated_intrastrandCrossLinks5
                this.intrastrandCrossLinks5 = this.calcIntrastrandCrossLinks5();
                this.validated_intrastrandCrossLinks5 = this.validated;
            end
            
            value = this.intrastrandCrossLinks5;
        end
        
        function value = calcIntrastrandCrossLinks5(this)
            value = this.shiftCircularSparseMatBase3Prime(this.intrastrandCrossLinks);
        end
        
        %boolean (genome length x 2) indicating positions of strand breaks
        %3' to bases
        function value = get.strandBreaks3(this)
            if this.validated_strandBreaks > this.validated_strandBreaks3
                this.strandBreaks3 = this.calcStrandBreaks3();
                this.validated_strandBreaks3 = this.validated;
            end
            
            value = this.strandBreaks3;
        end
        
        function value = calcStrandBreaks3(this)
            value = this.unshiftCircularSparseMatBond3Prime(this.strandBreaks);
        end
        
        %boolean (genome length x 2) indicating positions of strand breaks
        %5' to bases
        function value = get.strandBreaks5(this)
            if this.validated_strandBreaks > this.validated_strandBreaks5
                this.strandBreaks5 = this.calcStrandBreaks5();
                this.validated_strandBreaks5 = this.validated;
            end
            
            value = this.strandBreaks5;
        end
        
        function value = calcStrandBreaks5(this)
            value = this.unshiftCircularSparseMatBond5Prime(this.strandBreaks);
        end
        
        %boolean (genome length x 2) indicating positions of holliday junctions
        %3' to bases
        function value = get.hollidayJunctions3(this)
            if this.validated_hollidayJunctions > this.validated_hollidayJunctions3
                this.hollidayJunctions3 = this.calcHollidayJunctions3();
                this.validated_hollidayJunctions3 = this.validated;
            end
            
            value = this.hollidayJunctions3;
        end
        
        function value = calcHollidayJunctions3(this)
            value = this.shiftCircularSparseMatBond5Prime(this.hollidayJunctions);
        end
        
        %boolean (genome length x 2) indicating positions of holliday junctions
        %5' to bases
        function value = get.hollidayJunctions5(this)
            if this.validated_hollidayJunctions > this.validated_hollidayJunctions5
                this.hollidayJunctions5 = this.calcHollidayJunctions5();
                this.validated_hollidayJunctions5 = this.validated;
            end
            
            value = this.hollidayJunctions5;
        end
        
        function value = calcHollidayJunctions5(this)
            value = this.shiftCircularSparseMatBond3Prime(this.hollidayJunctions);
        end
        
        %boolean (genome length x 2) indicating positions of single
        %strand breaks -- strand breaks excluding
        %- double strand breaks that are part of double strand breaks
        %- strand breaks adjacent to gap sites
        function value = get.singleStrandBreaks(this)
            if this.validated_damaged > this.validated_singleStrandBreaks
                this.singleStrandBreaks = this.calcSingleStrandBreaks();
                this.validated_singleStrandBreaks = this.validated;
            end
            
            value = this.singleStrandBreaks;
        end
        
        function value = calcSingleStrandBreaks(this)           
            value = this.strandBreaks;
            
            %exclude strand breaks that are part of double strand breaks
            value(find(this.doubleStrandBreaks)) = 0; %#ok<FNDSB>
            
            %exclude strand breaks adjacent to other damage (except holliday
            %junctions)
            otherDamages = this.getDamagedSites(true, false, true, true, false, false, false);
            value(find( ...
                this.shiftCircularSparseMatBond3Prime(otherDamages) | ...
                this.shiftCircularSparseMatBond5Prime(otherDamages) ...
                )) = 0; %#ok<FNDSB>
            
            %cast to logical sparse mat
            value = valueCast(value, 'logical');
        end
        
        function value = get.doubleStrandBreaks(this)
            if this.validated_strandBreaks > this.validated_doubleStrandBreaks
                this.doubleStrandBreaks = this.calcDoubleStrandBreaks();
                this.validated_doubleStrandBreaks = this.validated;
            end
            
            value = this.doubleStrandBreaks;
        end
        
        function value = calcDoubleStrandBreaks(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            if this.doubleStrandBreakSeparation ~= 1
                throw(MException('DNARepair:error','Simulation only valid for doubleStrandBreakSeparation equal 1'));
            end
            
            value = this.strandBreaks(:, 1:2:end, :) & this.strandBreaks(:, 2:2:end, :);
            value = [value(:,1) value(:,1) value(:,2) value(:,2)];
        end
                
        function value = get.restrictableMunIRMSites(this)
            if ...
                    this.validated_damagedBases > this.validated_munIRMSiteRestrictionStatus || ...
                    this.validated_strandBreaks > this.validated_munIRMSiteRestrictionStatus
                this.restrictableMunIRMSites = this.calcRestrictableMunIRMSites();
                this.validated_munIRMSiteRestrictionStatus = this.validated;
            end
            
            value = this.restrictableMunIRMSites;
        end
        
        function value = get.hemiunmethylatedMunIRMSites(this)
            if this.validated_damagedBases > this.validated_munIRMSiteMethylationStatus
                this.hemiunmethylatedMunIRMSites = this.calcHemiunmethylatedMunIRMSites();
                this.validated_munIRMSiteMethylationStatus = this.validated;
            end
            
            value = this.hemiunmethylatedMunIRMSites;
        end
        
        function value = calcRestrictableMunIRMSites(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            dr = this.dnaRepair;
            nSites = size(dr.RM_MunI_RecognitionSites, 1);
            methylationPosStrnds = [
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))   ones(nSites, 1);
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1)) 3*ones(nSites, 1);
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2)) 2*ones(nSites, 1);
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2)) 4*ones(nSites, 1)];
            restrictionPosStrnds = [
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_RestrictionPositions(1))   ones(nSites, 1);
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_RestrictionPositions(1)) 3*ones(nSites, 1);
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_RestrictionPositions(2)) 2*ones(nSites, 1);
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_RestrictionPositions(2)) 4*ones(nSites, 1)];
            
            isPositionMethylated = reshape(this.damagedBases(methylationPosStrnds), [], 2) == this.metabolite.m6ADIndexs;
            isPositionDamaged = [
                cat(3, ...
                    reshape(this.damagedSites_shifted_incm6AD([dr.RM_MunI_RecognitionSites(:)   ones(6*nSites, 1)]), [], 6), ...
                    reshape(this.damagedSites_shifted_incm6AD([dr.RM_MunI_RecognitionSites(:) 2*ones(6*nSites, 1)]), [], 6))
                cat(3, ...
                    reshape(this.damagedSites_shifted_incm6AD([dr.RM_MunI_RecognitionSites(:) 3*ones(6*nSites, 1)]), [], 6), ...
                    reshape(this.damagedSites_shifted_incm6AD([dr.RM_MunI_RecognitionSites(:) 4*ones(6*nSites, 1)]), [], 6))];
            isPositionStrandBreaks = [
                cat(3, ...
                    this.strandBreaks([reshape(dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_RestrictionPositions(1)), [], 1),   ones(nSites, 1)]),...
                    this.strandBreaks([reshape(dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_RestrictionPositions(2)), [], 1), 2*ones(nSites, 1)])), ...
                cat(3, ...
                    this.strandBreaks([reshape(dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_RestrictionPositions(1)), [], 1), 3*ones(nSites, 1)]),...
                    this.strandBreaks([reshape(dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_RestrictionPositions(2)), [], 1), 4*ones(nSites, 1)]))];
            isPositionDamaged(isPositionStrandBreaks) = 0;
            isSiteUndamaged = ~any(any(isPositionDamaged, 2), 3);
            
            isSiteUnmethylated = ~any(isPositionMethylated, 2);
            isSitePolymerized = this.isRegionPolymerized(restrictionPosStrnds, 1, false, false, false);
            
            value = CircularSparseMat(...
                restrictionPosStrnds([isSiteUnmethylated; isSiteUnmethylated] & [isSiteUndamaged; isSiteUndamaged] & isSitePolymerized, :), ...
                true, [this.sequenceLen this.nCompartments], 1);
        end
        
        function value = calcHemiunmethylatedMunIRMSites(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            dr = this.dnaRepair;
            nSites = size(dr.RM_MunI_RecognitionSites, 1);
            methylationPosStrnds = [
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))   ones(nSites, 1);
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1)) 3*ones(nSites, 1);
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2)) 2*ones(nSites, 1);
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2)) 4*ones(nSites, 1)];
            
            isPositionMethylated = reshape(this.damagedBases(methylationPosStrnds), [], 2) == this.metabolite.m6ADIndexs;
            
            isSiteMethylated = all(isPositionMethylated, 2);
            isSiteUnmethylated = ~any(isPositionMethylated, 2);
            isSiteHemimethylated = ~isSiteMethylated & ~isSiteUnmethylated;
            
            value = CircularSparseMat(...
                methylationPosStrnds([isSiteHemimethylated; isSiteHemimethylated] & ~isPositionMethylated(:), :),  ...
                1, [this.sequenceLen this.nCompartments], 1);
        end
                
        function value = get.dryWeight(this)
            if max(this.validated_polymerizedRegions, this.validated_damaged) > this.validated_dryWeight
                this.dryWeight = this.calcDryWeight();
                this.validated_dryWeight = this.validated;
            end
           
            value = this.dryWeight;
        end
        
        function value = calcDryWeight(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            %time
            numTime = size(this.abasicSites, 3);
            
            %mass of undamaged DNA
            baseCounts = zeros(4, numTime);
            bonds = zeros(1, numTime);
            [positionsStrandTimes, lengths] = find(this.getStrandView('polymerizedRegions'));
            positionsStrandTimes = [positionsStrandTimes ones(size(positionsStrandTimes, 1), 3 - size(positionsStrandTimes, 2))];
            
            for i = 1:size(positionsStrandTimes, 1)
                if lengths(i) == this.sequenceLen
                    baseCounts(:, positionsStrandTimes(i, 3)) = ...
                        baseCounts(:, positionsStrandTimes(i, 3)) + ...
                        getBaseCounts(this.sequence, positionsStrandTimes(i, 2));
                else
                    baseCounts(:, positionsStrandTimes(i, 3)) = ...
                        baseCounts(:, positionsStrandTimes(i, 3)) + ...
                        this.sequence.subsequenceBaseCounts(positionsStrandTimes(i,1) + (0:lengths(i)-1)', positionsStrandTimes(i, 2));
                end
                bonds(1, positionsStrandTimes(i,3)) = ...
                    bonds(1, positionsStrandTimes(i,3)) + lengths(i) - 1;
                
                if ...
                        lengths(i) == this.sequenceLen || ...
                        (positionsStrandTimes(i, 1) + lengths(i) - 1 == this.sequenceLen && ...
                        ismember([1 positionsStrandTimes(i, 2:end)], positionsStrandTimes, 'rows'))
                    bonds(1, positionsStrandTimes(i,3)) = ...
                        bonds(1, positionsStrandTimes(i,3)) + 1;
                end
            end
            
            value = this.metabolite.molecularWeights(this.metabolite.dnmpIndexs)' * baseCounts ...
                - (ConstantUtil.elements.H + ConstantUtil.elements.O) * bonds;
            
            %mass represented by damage
            for k = 1:numTime
                %gap sites
                value(:,k) = value(:,k) - ...
                    (this.metabolite.molecularWeights(this.metabolite.dr5pIndexs) - ...
                    2 * ConstantUtil.elements.H) * collapse(this.gapSites(:,:,k));
                
                %abasic sites
                [position, index] = find(this.abasicSites(:,:,k));
                value(:,k) = value(:,k) - ...
                    length(index) * this.metabolite.molecularWeights(this.metabolite.waterIndexs) - ...
                    this.sequence.subsequenceBaseCounts(position)' * ...
                    this.metabolite.molecularWeights(this.metabolite.dnmpIndexs);
                
                %damaged sugar-phosphates
                [~, index] = find(this.damagedSugarPhosphates(:,:,k));
                value(:,k) = value(:,k) + ...
                    sum(this.metabolite.molecularWeights(index)) - ...
                    length(index) * this.metabolite.molecularWeights(this.metabolite.dr5pIndexs);
                
                %damaged bases
                [position, index] = find(this.damagedBases(:,:,k));
                value(:,k) = value(:,k) + ...
                    sum(this.metabolite.molecularWeights(index)) - ...
                    this.sequence.subsequenceBaseCounts(position)' * ...
                    this.metabolite.molecularWeights(this.metabolite.unmodifiedBaseIndexs);
                
                %intrastrand cross links
                [position, index] = find(this.intrastrandCrossLinks(:,:,k));
                value(:,k) = value(:,k) + ...
                    sum(this.metabolite.molecularWeights(index)) - ...
                    this.sequence.subsequenceBaseCounts([position(:,1) position(:,1) + 1], position(:,2))' * ...
                    this.metabolite.molecularWeights(this.metabolite.dnmpIndexs);
                
                %strand breaks
                value(:,k) = value(:,k) + ...
                    (ConstantUtil.elements.H + ConstantUtil.elements.O)*...
                    collapse(this.strandBreaks(:,:,k));
            end
            
            value = permute(value, [1 3 2]);
            value = value / ConstantUtil.nAvogadro;
        end
    end
    
    %helper methods
    methods (Static)
        function value = shiftCircularSparseMatBase3Prime(sparseMat, varargin)
            import edu.stanford.covert.cell.sim.state.Chromosome;
            import edu.stanford.covert.util.CircularSparseMat;
            
            if ~isa(sparseMat,'CircularSparseMat') || ~isDimCircular(sparseMat, 1)
                throw(MException('Chromosome:error','sparseMat must be a circular sparse mat with a circular first dimension'));
            end
            
            [subs, vals] = find(sparseMat);
            subs = Chromosome.shiftPositionsStrandsBase3Prime(subs, varargin{:});
            value = CircularSparseMat(subs, vals, size(sparseMat), 1);
        end
        
        function value = shiftCircularSparseMatBase5Prime(sparseMat, varargin)
            import edu.stanford.covert.cell.sim.state.Chromosome;
            import edu.stanford.covert.util.CircularSparseMat;
            
            if ~isa(sparseMat,'CircularSparseMat') || ~isDimCircular(sparseMat, 1)
                throw(MException('Chromosome:error','sparseMat must be a circular sparse mat with a circular first dimension'));
            end
            
            [subs, vals] = find(sparseMat);
            subs = Chromosome.shiftPositionsStrandsBase5Prime(subs, varargin{:});
            value = CircularSparseMat(subs, vals, size(sparseMat), 1);
        end
        
        function value = shiftCircularSparseMatBond3Prime(sparseMat, varargin)
            import edu.stanford.covert.cell.sim.state.Chromosome;
            import edu.stanford.covert.util.CircularSparseMat;
            
            if ~isa(sparseMat,'CircularSparseMat') || ~isDimCircular(sparseMat, 1)
                throw(MException('Chromosome:error','sparseMat must be a circular sparse mat with a circular first dimension'));
            end
            
            [subs, vals] = find(sparseMat);
            subs = Chromosome.shiftPositionsStrandsBond3Prime(subs, varargin{:});
            value = CircularSparseMat(subs, vals, size(sparseMat), 1);
        end
        
        function value = unshiftCircularSparseMatBond3Prime(sparseMat, varargin)
            import edu.stanford.covert.cell.sim.state.Chromosome;
            import edu.stanford.covert.util.CircularSparseMat;
            
            if ~isa(sparseMat,'CircularSparseMat') || ~isDimCircular(sparseMat, 1)
                throw(MException('Chromosome:error','sparseMat must be a circular sparse mat with a circular first dimension'));
            end
            
            [subs, vals] = find(sparseMat);
            subs = Chromosome.unshiftPositionsStrandsBond3Prime(subs, varargin{:});
            value = CircularSparseMat(subs, vals, size(sparseMat), 1);
        end
        
        function value = shiftCircularSparseMatBond5Prime(sparseMat, varargin)
            import edu.stanford.covert.cell.sim.state.Chromosome;
            import edu.stanford.covert.util.CircularSparseMat;
            
            if ~isa(sparseMat,'CircularSparseMat') || ~isDimCircular(sparseMat, 1)
                throw(MException('Chromosome:error','sparseMat must be a circular sparse mat with a circular first dimension'));
            end
            
            [subs, vals] = find(sparseMat);
            subs = Chromosome.shiftPositionsStrandsBond5Prime(subs, varargin{:});
            value = CircularSparseMat(subs, vals, size(sparseMat), 1);
        end
        
        function value = unshiftCircularSparseMatBond5Prime(sparseMat, varargin)
            import edu.stanford.covert.cell.sim.state.Chromosome;
            import edu.stanford.covert.util.CircularSparseMat;
            
            if ~isa(sparseMat,'CircularSparseMat') || ~isDimCircular(sparseMat, 1)
                throw(MException('Chromosome:error','sparseMat must be a circular sparse mat with a circular first dimension'));
            end
            
            [subs, vals] = find(sparseMat);
            subs = Chromosome.unshiftPositionsStrandsBond5Prime(subs, varargin{:});
            value = CircularSparseMat(subs, vals, size(sparseMat), 1);
        end
        
        function positionsStrands = shiftPositionsStrandsBase3Prime(positionsStrands, lengths)
            if nargin < 2
                lengths = 1; 
            end
            
            positionsStrands(mod(positionsStrands(:,2),2)==1, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==1, 1) + lengths;
            positionsStrands(mod(positionsStrands(:,2),2)==0, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==0, 1) - lengths;
        end
        
        function positionsStrands = shiftPositionsStrandsBase5Prime(positionsStrands, lengths)
            if nargin < 2
                lengths = 1; 
            end
            
            positionsStrands(mod(positionsStrands(:,2),2)==1, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==1, 1) - lengths;
            positionsStrands(mod(positionsStrands(:,2),2)==0, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==0, 1) + lengths;
        end
        
        function positionsStrands = shiftPositionsStrandsBond3Prime(positionsStrands, lengths)
            if nargin < 2
                lengths = 1; 
            end
            
            positionsStrands(mod(positionsStrands(:,2),2)==0, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==0, 1) - lengths;
        end
        
        function positionsStrands = unshiftPositionsStrandsBond3Prime(positionsStrands, lengths)
            if nargin < 2
                lengths = 1; 
            end
            
            positionsStrands(mod(positionsStrands(:,2),2)==0, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==0, 1) + lengths;
        end
        
        function positionsStrands = shiftPositionsStrandsBond5Prime(positionsStrands, lengths)
            if nargin < 2
                lengths = 1;
            end
            
            positionsStrands(mod(positionsStrands(:,2),2)==1, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==1, 1) - lengths;
        end
        
        function positionsStrands = unshiftPositionsStrandsBond5Prime(positionsStrands, lengths)
            if nargin < 2
                lengths = 1;
            end
            
            positionsStrands(mod(positionsStrands(:,2),2)==1, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==1, 1) + lengths;
        end
        
        function [footprint3Prime, footprint5Prime] = calculateFootprintOverhangs(footprint)
            footprint5Prime = ceil((footprint-1)/2);
            footprint3Prime = footprint - 1 - footprint5Prime;
        end
    end
    
    %printing
    methods
        %print state
        function disp(this)
            %superclass method
            this.disp@edu.stanford.covert.cell.sim.CellState();

            %numbers of DNA damages
            fprintf('%24s\t%4s\n','Damage Type','No.');
            fprintf('%24s\t%4s\n',repmat('=',1,24), repmat('=',1,4));
            fprintf('%24s\t%4d\n','Gap Sites',               collapse(this.gapSites));
            fprintf('%24s\t%4d\n','Abasic Sites',            collapse(this.abasicSites));
            fprintf('%24s\t%4d\n','Damaged sugar phosphates',collapse(this.damagedSugarPhosphates));
            fprintf('%24s\t%4d\n','Damaged bases',           collapse(this.damagedBases));
            fprintf('%24s\t%4d\n','Intrastrand cross links', collapse(this.intrastrandCrossLinks));
            fprintf('%24s\t%4d\n','Strand breaks',           collapse(this.strandBreaks));
            fprintf('%24s\t%4d\n','Holliday junctions',      collapse(this.hollidayJunctions));
            fprintf('\n');
            
            %strand break classification used in track structure models
            sbc = this.strandBreakClassification();
            fprintf('%24s\t%4s\n','Strand Break','No.')
            fprintf('%24s\t%4s\n',repmat('=',1,24), repmat('=',1,4))
            fprintf('%24s\t%4d\n','NB',    sbc(this.strandBreakClassification_index_NB));
            fprintf('%24s\t%4d\n','SSB',   sbc(this.strandBreakClassification_index_SSB));
            fprintf('%24s\t%4d\n','SSB+',  sbc(this.strandBreakClassification_index_SSB_));
            fprintf('%24s\t%4d\n','2SSB',  sbc(this.strandBreakClassification_index_2SSB));
            fprintf('%24s\t%4d\n','DSB',   sbc(this.strandBreakClassification_index_DSB));
            fprintf('%24s\t%4d\n','DSB+',  sbc(this.strandBreakClassification_index_DSB_));
            fprintf('%24s\t%4d\n','DSB++', sbc(this.strandBreakClassification_index_DSB__));
            fprintf('\n');
            
            %list of damaged genes
            fprintf('%13s\n','Damaged Genes');
            fprintf('%13s\n',repmat('=',1,13));
            damagedGeneWholeCellModelIDs = this.gene.wholeCellModelIDs(~all(this.accessibleGenes,2));
            damagedGeneWholeCellNames    = this.gene.names(~all(this.accessibleGenes,2));
            for i=1:length(damagedGeneWholeCellModelIDs);
                damagedGeneWholeCellModelID = damagedGeneWholeCellModelIDs{i};
                damagedGeneWholeCellName    = damagedGeneWholeCellNames{i};
                if length(damagedGeneWholeCellModelID)>12; damagedGeneWholeCellModelID=[damagedGeneWholeCellModelID(1:8) ' ...']; end;
                if length(damagedGeneWholeCellName)>32; damagedGeneWholeCellName=[damagedGeneWholeCellName(1:28) ' ...']; end;
                fprintf('%12s\t%32s\n',damagedGeneWholeCellModelID,damagedGeneWholeCellName);
            end
            fprintf('\n');
            
            %list of damaged transcription units
            fprintf('%27s\n','Damaged Transcription Units');
            fprintf('%27s\n',repmat('=',1,27));
            damagedTranscriptionUnitWholeCellModelIDs = this.transcriptionUnitWholeCellModelIDs(~all(this.accessibleTranscriptionUnits,2));
            damagedTranscriptionUnitWholeCellNames    = this.transcriptionUnitNames(~all(this.accessibleTranscriptionUnits,2));
            for i=1:length(damagedTranscriptionUnitWholeCellModelIDs);
                damagedTranscriptionUnitWholeCellModelID = damagedTranscriptionUnitWholeCellModelIDs{i};
                damagedTranscriptionUnitWholeCellName    = damagedTranscriptionUnitWholeCellNames{i};
                if length(damagedTranscriptionUnitWholeCellModelID)>12; damagedTranscriptionUnitWholeCellModelID=[damagedTranscriptionUnitWholeCellModelID(1:8) ' ...']; end;
                if length(damagedTranscriptionUnitWholeCellName)>32; damagedTranscriptionUnitWholeCellName=[damagedTranscriptionUnitWholeCellName(1:28) ' ...']; end;
                fprintf('%12s\t%32s\n',damagedTranscriptionUnitWholeCellModelID,damagedTranscriptionUnitWholeCellName);
            end
            fprintf('\n');
        end
    end
end
