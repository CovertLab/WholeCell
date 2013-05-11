%Replication process test case
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/9/2010
classdef Replication_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = Replication_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end
    
    methods
        function loadSimpleTestFixture(this)
            %% import classes
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            import edu.stanford.covert.cell.sim.state.Chromosome;
            
            %% process
            m = this.process;
            
            %% constants
            m.oriCPosition = 1;
            m.terCPosition = 290038;
            
            m.primerLength                      = 11;
            m.dnaPolymeraseElongationRate       = 100;
            m.ssbComplexSpacing                 = 30;
            m.okazakiFragmentMeanLength         = 1500;
            m.ligaseRate                        = 0.04;
            m.laggingBackupClampReloadingLength = 750;
            m.startingOkazakiLoopLength         = 50;
            
            %% metabolites
            m.substrateWholeCellModelIDs = { %whole cell model IDs of substrates
                'DATP'; 'DCTP'; 'DGTP'; 'DTTP';
                'ATP'; 'CTP'; 'GTP'; 'UTP';
                'PPI'; 'H2O'; 'H'; 'NAD'; 'NMN'; 'ADP'; 'AMP'; 'PI';
                };
            m.substrateNames = m.substrateWholeCellModelIDs;
            
            m.substrateIndexs_dntp        = (1:4)';
            m.substrateIndexs_ntp         = (5:8)';
            m.substrateIndexs_diphosphate = 9;
            m.substrateIndexs_water       = 10;
            m.substrateIndexs_hydrogen    = 11;
            m.substrateIndexs_nad         = 12;
            m.substrateIndexs_nmn         = 13;
            m.substrateIndexs_atp         = 5;
            m.substrateIndexs_adp         = 14;
            m.substrateIndexs_amp         = 15;
            m.substrateIndexs_phosphate   = 16;
            
            m.substrateMetaboliteLocalIndexs = (1:numel(m.substrateWholeCellModelIDs))';
            m.substrateMetaboliteGlobalIndexs = m.substrateMetaboliteLocalIndexs;
            
            m.substrateMolecularWeights = [
                487.1495;  463.1248;  503.1489;  478.1361;  503.1489;  479.1242;  519.1483;
                480.1090;  174.9513;   18.0152;    1.0079;  662.4162;  333.2107;  424.1769;
                345.2049;   95.9793];
            
            %% proteins
            m.enzymeWholeCellModelIDs = { %enzyme whole cell model ids
                'REPLISOME';                                             %replisome
                'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE'; %DNA-directed DNA polymerase (2) core + beta-clamp + gamma-complex
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX';          %DNA-directed DNA polymerase core + beta-clamp + gamma-complex
                'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE';                %DNA-directed DNA polymerase core + beta-clamp
                'DNA_POLYMERASE_CORE';                                   %DNA-directed DNA polymerase core
                'DNA_POLYMERASE_GAMMA_COMPLEX';                          %DNA-directed DNA polymerase gamma complex
                'MG_001_MONOMER';                                        %DNA polymerase III, beta subunit
                'MG_001_DIMER';                                          %DNA polymerase III, beta subunit
                'MG_094_HEXAMER';                                        %replicative DNA helicase
                'MG_254_MONOMER';                                        %DNA ligase, NAD-dependent
                'MG_250_MONOMER';                                        %DNA primase
                'MG_091_TETRAMER';                                       %single-strand binding protein family
                'MG_091_OCTAMER';                                        %single-strand binding protein family
                };
            m.enzymeNames = m.enzymeWholeCellModelIDs;
            
            m.enzymeIndexs_replisome                         = 1;
            m.enzymeIndexs_2coreBetaClampGammaComplexPrimase = 2;
            m.enzymeIndexs_coreBetaClampGammaComplex         = 3;
            m.enzymeIndexs_coreBetaClampPrimase              = 4;
            m.enzymeIndexs_core                              = 5;
            m.enzymeIndexs_gammaComplex                      = 6;
            m.enzymeIndexs_betaClampMonomer                  = 7;
            m.enzymeIndexs_betaClamp                         = 8;
            m.enzymeIndexs_helicase                          = 9;
            m.enzymeIndexs_ligase                            = 10;
            m.enzymeIndexs_primase                           = 11;
            m.enzymeIndexs_ssb4mer                           = 12;
            m.enzymeIndexs_ssb8mer                           = 13;
            
            m.enzymeMonomerLocalIndexs = [m.enzymeIndexs_betaClampMonomer; m.enzymeIndexs_ligase; m.enzymeIndexs_primase];
            m.enzymeComplexLocalIndexs = setdiff((1:numel(m.enzymeWholeCellModelIDs))', m.enzymeMonomerLocalIndexs);
            m.enzymeMonomerGlobalIndexs = (1:numel(m.enzymeMonomerLocalIndexs))';
            m.enzymeComplexGlobalIndexs = (1:numel(m.enzymeComplexLocalIndexs))';
            
            m.enzymeGlobalIndexs = zeros(size(m.enzymeWholeCellModelIDs));
            m.enzymeGlobalIndexs(m.enzymeMonomerLocalIndexs) = m.enzymeMonomerGlobalIndexs;
            m.enzymeGlobalIndexs(m.enzymeComplexLocalIndexs) = m.enzymeComplexGlobalIndexs;
            
            m.enzymeComposition = zeros(numel(m.enzymeWholeCellModelIDs));
            m.enzymeComposition(m.enzymeIndexs_core,                m.enzymeIndexs_replisome)                         = 2;
            m.enzymeComposition(m.enzymeIndexs_betaClampMonomer,    m.enzymeIndexs_replisome)                         = 4;
            m.enzymeComposition(m.enzymeIndexs_gammaComplex,        m.enzymeIndexs_replisome)                         = 1;
            m.enzymeComposition(m.enzymeIndexs_helicase,            m.enzymeIndexs_replisome)                         = 1;
            m.enzymeComposition(m.enzymeIndexs_primase,             m.enzymeIndexs_replisome)                         = 1;
            m.enzymeComposition(m.enzymeIndexs_core,                m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            m.enzymeComposition(m.enzymeIndexs_betaClampMonomer,    m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            m.enzymeComposition(m.enzymeIndexs_gammaComplex,        m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 1;
            m.enzymeComposition(m.enzymeIndexs_primase,             m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 1;
            m.enzymeComposition(m.enzymeIndexs_core,                m.enzymeIndexs_coreBetaClampGammaComplex)         = 1;
            m.enzymeComposition(m.enzymeIndexs_betaClampMonomer,    m.enzymeIndexs_coreBetaClampGammaComplex)         = 2;
            m.enzymeComposition(m.enzymeIndexs_gammaComplex,        m.enzymeIndexs_coreBetaClampGammaComplex)         = 1;
            m.enzymeComposition(m.enzymeIndexs_core,                m.enzymeIndexs_coreBetaClampPrimase)              = 1;
            m.enzymeComposition(m.enzymeIndexs_betaClampMonomer,    m.enzymeIndexs_coreBetaClampPrimase)              = 2;
            m.enzymeComposition(m.enzymeIndexs_primase,             m.enzymeIndexs_coreBetaClampPrimase)              = 1;
            m.enzymeComposition(m.enzymeIndexs_betaClampMonomer,    m.enzymeIndexs_betaClamp)                         = 2;
            m.enzymeComposition(m.enzymeIndexs_ssb4mer,             m.enzymeIndexs_ssb8mer)                           = 2;
            
            m.enzymeMolecularWeights = zeros(size(m.enzymeWholeCellModelIDs));
            m.enzymeMolecularWeights(m.enzymeIndexs_core) = 267939.067600;
            m.enzymeMolecularWeights(m.enzymeIndexs_gammaComplex) = 339655.731900;
            m.enzymeMolecularWeights(m.enzymeIndexs_betaClampMonomer) = 44289.824700;
            m.enzymeMolecularWeights(m.enzymeIndexs_helicase) = 323295.523800;
            m.enzymeMolecularWeights(m.enzymeIndexs_ligase) = 75426.424900;
            m.enzymeMolecularWeights(m.enzymeIndexs_primase) = 71135.076600;
            m.enzymeMolecularWeights(m.enzymeIndexs_ssb4mer) = 71854.750000;
            m.enzymeMolecularWeights(m.enzymeMolecularWeights == 0) = m.enzymeComposition(:, m.enzymeMolecularWeights == 0)' * m.enzymeMolecularWeights;
            
            %% chromosome
            c = Chromosome([], []);
            m.states = {c};
            m.chromosome = c;
            
            %constants
            c.relaxedBasesPerTurn = 10.5;
            c.equilibriumSuperhelicalDensity = -0.06;
            c.supercoiledSuperhelicalDensityTolerance = 0.1;
            
            %sequence
            bases = 'ACGT';
            seq = bases(randi(4, 580076, 1));
            c.sequence = ChromosomeSequence(seq);
            c.sequenceLen = length(seq);
            c.sequenceGCContent = getGCContent(c.sequence);
            
            %metabolites
            c.metabolite.molecularWeights = [m.substrateMolecularWeights; 329.2055; 305.1808; 345.2049; 320.1921; 212.0942; 149.1530];
            c.metabolite.dnmpIndexs = numel(m.substrateMolecularWeights) + (1:4)';
            c.metabolite.waterIndexs = m.substrateIndexs_water;
            c.metabolite.dr5pIndexs = numel(m.substrateMolecularWeights) + 5;
            c.metabolite.m6ADIndexs = numel(m.substrateMolecularWeights) + 6;
            
            %other proteins
            m.complexIndexs_DnaA_1mer_ATP = numel(m.enzymeComplexGlobalIndexs)+1;
            m.complexIndexs_DnaA_7mer_ATP = numel(m.enzymeComplexGlobalIndexs)+2;
            
            %footprints
            m.enzymeDNAFootprints = zeros(size(m.enzymeWholeCellModelIDs));
            m.enzymeDNAFootprints(m.enzymeIndexs_replisome) = 44;
            m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 49;
            m.enzymeDNAFootprints(m.enzymeIndexs_coreBetaClampGammaComplex)  = 49;
            m.enzymeDNAFootprints(m.enzymeIndexs_coreBetaClampPrimase)       = 49;
            m.enzymeDNAFootprints(m.enzymeIndexs_core)                       = 24;
            m.enzymeDNAFootprints(m.enzymeIndexs_gammaComplex)               = 26;
            m.enzymeDNAFootprints(m.enzymeIndexs_betaClamp)                  = 25;
            m.enzymeDNAFootprints(m.enzymeIndexs_betaClampMonomer)           = 25;
            m.enzymeDNAFootprints(m.enzymeIndexs_helicase)                   = 20;
            m.enzymeDNAFootprints(m.enzymeIndexs_ligase)                     = 19;
            m.enzymeDNAFootprints(m.enzymeIndexs_primase)                    = 14;
            m.enzymeDNAFootprints(m.enzymeIndexs_ssb4mer)                    = 72;
            m.enzymeDNAFootprints(m.enzymeIndexs_ssb8mer)                    = 145;
            [m.enzymeDNAFootprints3Prime, m.enzymeDNAFootprints5Prime] = m.chromosome.calculateFootprintOverhangs(m.enzymeDNAFootprints);
            
            c.monomerDNAFootprints = m.enzymeDNAFootprints(m.enzymeMonomerLocalIndexs);
            c.complexDNAFootprints = m.enzymeDNAFootprints(m.enzymeComplexLocalIndexs);
            c.monomerDNAFootprintBindingStrandedness = repmat(c.dnaStrandedness_dsDNA, size(c.monomerDNAFootprints));
            c.complexDNAFootprintBindingStrandedness = repmat(c.dnaStrandedness_dsDNA, size(c.complexDNAFootprints));
            c.monomerDNAFootprintRegionStrandedness = repmat(c.dnaStrandedness_dsDNA, size(c.monomerDNAFootprints));
            c.complexDNAFootprintRegionStrandedness = repmat(c.dnaStrandedness_dsDNA, size(c.complexDNAFootprints));
            
            c.monomerDNAFootprintBindingStrandedness(m.enzymeMonomerGlobalIndexs(m.enzymeMonomerLocalIndexs == m.enzymeIndexs_betaClampMonomer)) = c.dnaStrandedness_ssDNA;
            c.monomerDNAFootprintBindingStrandedness(m.enzymeMonomerGlobalIndexs(m.enzymeMonomerLocalIndexs == m.enzymeIndexs_ligase)) = c.dnaStrandedness_dsDNA;
            c.monomerDNAFootprintBindingStrandedness(m.enzymeMonomerGlobalIndexs(m.enzymeMonomerLocalIndexs == m.enzymeIndexs_primase)) = c.dnaStrandedness_ssDNA;
            c.complexDNAFootprintBindingStrandedness(m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase)) = c.dnaStrandedness_ssDNA;
            c.complexDNAFootprintBindingStrandedness(m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex)) = c.dnaStrandedness_ssDNA;
            c.complexDNAFootprintBindingStrandedness(m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase)) = c.dnaStrandedness_ssDNA;
            c.complexDNAFootprintBindingStrandedness(m.enzymeGlobalIndexs(m.enzymeIndexs_betaClamp)) = c.dnaStrandedness_ssDNA;
            c.complexDNAFootprintBindingStrandedness(m.enzymeGlobalIndexs(m.enzymeIndexs_helicase)) = c.dnaStrandedness_ssDNA;
            c.complexDNAFootprintBindingStrandedness(m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_ssb4mer)) = c.dnaStrandedness_ssDNA;
            c.complexDNAFootprintBindingStrandedness(m.enzymeGlobalIndexs(m.enzymeIndexs_ssb8mer)) = c.dnaStrandedness_ssDNA;
            c.complexDNAFootprintBindingStrandedness(m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_core)) = c.dnaStrandedness_ssDNA;
            c.complexDNAFootprintBindingStrandedness(m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_gammaComplex)) = c.dnaStrandedness_ssDNA;
            c.complexDNAFootprintBindingStrandedness(m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_replisome)) = c.dnaStrandedness_ssDNA;
            
            c.monomerDNAFootprintRegionStrandedness(m.enzymeMonomerGlobalIndexs(m.enzymeMonomerLocalIndexs == m.enzymeIndexs_betaClampMonomer)) = c.dnaStrandedness_xsDNA;
            c.monomerDNAFootprintRegionStrandedness(m.enzymeMonomerGlobalIndexs(m.enzymeMonomerLocalIndexs == m.enzymeIndexs_ligase)) = c.dnaStrandedness_dsDNA;
            c.monomerDNAFootprintRegionStrandedness(m.enzymeMonomerGlobalIndexs(m.enzymeMonomerLocalIndexs == m.enzymeIndexs_primase)) = c.dnaStrandedness_xsDNA;
            c.complexDNAFootprintRegionStrandedness(m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase)) = c.dnaStrandedness_xsDNA;
            c.complexDNAFootprintRegionStrandedness(m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex)) = c.dnaStrandedness_xsDNA;
            c.complexDNAFootprintRegionStrandedness(m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase)) = c.dnaStrandedness_xsDNA;
            c.complexDNAFootprintRegionStrandedness(m.enzymeGlobalIndexs(m.enzymeIndexs_betaClamp)) = c.dnaStrandedness_xsDNA;
            c.complexDNAFootprintRegionStrandedness(m.enzymeGlobalIndexs(m.enzymeIndexs_helicase)) = c.dnaStrandedness_xsDNA;
            c.complexDNAFootprintRegionStrandedness(m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_ssb4mer)) = c.dnaStrandedness_ssDNA;
            c.complexDNAFootprintRegionStrandedness(m.enzymeGlobalIndexs(m.enzymeIndexs_ssb8mer)) = c.dnaStrandedness_ssDNA;
            c.complexDNAFootprintRegionStrandedness(m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_core)) = c.dnaStrandedness_xsDNA;
            c.complexDNAFootprintRegionStrandedness(m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_gammaComplex)) = c.dnaStrandedness_xsDNA;
            c.complexDNAFootprintRegionStrandedness(m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_replisome)) = c.dnaStrandedness_xsDNA;
            
            %release reactions
            nReactions = 18;
            c.reactionBoundMonomer = zeros(nReactions, 1);
            c.reactionBoundComplex = zeros(nReactions, 1);
            c.reactionMonomerCatalysisMatrix = zeros(nReactions, numel(m.enzymeMonomerGlobalIndexs));
            c.reactionComplexCatalysisMatrix = zeros(nReactions, numel(m.enzymeComplexGlobalIndexs));
            
            c.reactionBoundComplex(1:7) = m.enzymeGlobalIndexs(m.enzymeIndexs_ssb8mer);
            c.reactionBoundMonomer(8:13) = m.enzymeGlobalIndexs(m.enzymeIndexs_ligase);
            c.reactionBoundComplex(14) = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            c.reactionBoundComplex(15) = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex);
            c.reactionBoundComplex(16) = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase);
            c.reactionBoundComplex(17) = m.enzymeGlobalIndexs(m.enzymeIndexs_betaClamp);
            c.reactionBoundComplex(18) = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            
            c.reactionMonomerCatalysisMatrix(1, m.enzymeGlobalIndexs(m.enzymeIndexs_betaClampMonomer)) = 1;
            c.reactionComplexCatalysisMatrix(2, m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase)) = 1;
            c.reactionComplexCatalysisMatrix(3, m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex)) = 1;
            c.reactionComplexCatalysisMatrix(4, m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase)) = 1;
            c.reactionComplexCatalysisMatrix(5, m.enzymeGlobalIndexs(m.enzymeIndexs_helicase)) = 1;
            c.reactionComplexCatalysisMatrix(6, m.enzymeGlobalIndexs(m.enzymeIndexs_betaClamp)) = 1;
            c.reactionMonomerCatalysisMatrix(7, m.enzymeGlobalIndexs(m.enzymeIndexs_ligase)) = 1;
            
            c.reactionComplexCatalysisMatrix(8, m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase)) = 1;
            c.reactionComplexCatalysisMatrix(9, m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex)) = 1;
            c.reactionComplexCatalysisMatrix(10, m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase)) = 1;
            c.reactionComplexCatalysisMatrix(11, m.enzymeGlobalIndexs(m.enzymeIndexs_helicase)) = 1;
            c.reactionComplexCatalysisMatrix(12, m.enzymeGlobalIndexs(m.enzymeIndexs_betaClamp)) = 1;
            c.reactionMonomerCatalysisMatrix(13, m.enzymeGlobalIndexs(m.enzymeIndexs_betaClampMonomer)) = 1;
            
            c.reactionMonomerCatalysisMatrix(14:18, m.enzymeGlobalIndexs(m.enzymeIndexs_betaClampMonomer)) = 1;
            
            c.reactionThresholds = ...
                sum(c.reactionMonomerCatalysisMatrix, 2) + ...
                sum(c.reactionComplexCatalysisMatrix, 2);
            
            %primase binding locations
            m.primaseBindingLocations = {
                size(c.sequence, 1) + 1 - m.calculatePrimaseBindingLocations(size(c.sequence, 1) - m.terCPosition)
                m.calculatePrimaseBindingLocations(m.terCPosition)
                }';
            
            %% initial state
            m.substrates = zeros(size(m.substrateWholeCellModelIDs));
            m.enzymes = zeros(size(m.enzymeWholeCellModelIDs));
            m.enzymes(m.enzymeIndexs_core) = 30;
            m.enzymes(m.enzymeIndexs_gammaComplex) = 9;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 51;
            m.enzymes(m.enzymeIndexs_helicase) = 5;
            m.enzymes(m.enzymeIndexs_ligase) = 14;
            m.enzymes(m.enzymeIndexs_primase) = 11;
            m.enzymes(m.enzymeIndexs_ssb4mer) = 4;
            m.boundEnzymes = zeros(size(m.enzymeWholeCellModelIDs));
            
            c.initialize();
        end
    end
    
    %test evolve state
    methods
        function testConstants(this)
            m = this.process;
            c = m.chromosome;
            
            enzymeComposition = m.enzymeComposition;
            enzymeMolecularWeights = m.enzymeMolecularWeights;
            enzymeDNAFootprints = m.enzymeDNAFootprints;
            monomerDNAFootprintBindingStrandedness = c.monomerDNAFootprintBindingStrandedness(m.enzymeMonomerGlobalIndexs);
            complexDNAFootprintBindingStrandedness = c.complexDNAFootprintBindingStrandedness(m.enzymeComplexGlobalIndexs);
            monomerDNAFootprintRegionStrandedness = c.monomerDNAFootprintRegionStrandedness(m.enzymeMonomerGlobalIndexs);
            complexDNAFootprintRegionStrandedness = c.complexDNAFootprintRegionStrandedness(m.enzymeComplexGlobalIndexs);
            rxnIdxs = find(...
                (ismember(c.reactionBoundMonomer, m.enzymeMonomerGlobalIndexs) | ...
                ismember(c.reactionBoundComplex, m.enzymeComplexGlobalIndexs)) & ...
                sum(c.reactionMonomerCatalysisMatrix(:, m.enzymeMonomerGlobalIndexs), 2) == sum(c.reactionMonomerCatalysisMatrix, 2) & ...
                sum(c.reactionComplexCatalysisMatrix(:, m.enzymeComplexGlobalIndexs), 2) == sum(c.reactionComplexCatalysisMatrix, 2));
            [~, reactionBoundMonomer] = ismember(c.reactionBoundMonomer(rxnIdxs, :), m.enzymeMonomerGlobalIndexs);
            [~, reactionBoundComplex] = ismember(c.reactionBoundComplex(rxnIdxs, :), m.enzymeComplexGlobalIndexs);
            reactionMonomerCatalysisMatrix = c.reactionMonomerCatalysisMatrix(rxnIdxs, m.enzymeMonomerGlobalIndexs);
            reactionComplexCatalysisMatrix = c.reactionComplexCatalysisMatrix(rxnIdxs, m.enzymeComplexGlobalIndexs);
            
            this.loadSimpleTestFixture();
            c = m.chromosome;
            
            %check that oriC is at first base
            assertEqual(1, m.oriCPosition);
            
            %assert enzyme composition, molecular weights are correct
            assertEqual(enzymeComposition, m.enzymeComposition);
            assertElementsAlmostEqual(enzymeMolecularWeights, m.enzymeMolecularWeights);
            
            %assert footprints equal for (2)core-beta-gamma, core-beta-gamma, core-beta
            assertEqual(1, numel(unique(enzymeDNAFootprints([
                m.enzymeIndexs_2coreBetaClampGammaComplexPrimase;
                m.enzymeIndexs_coreBetaClampGammaComplex;
                m.enzymeIndexs_coreBetaClampPrimase]))));
            assertEqual(1, numel(unique(enzymeDNAFootprints([
                m.enzymeIndexs_betaClampMonomer;
                m.enzymeIndexs_betaClamp]))));
            
            %assert footprints additive for beta-clamp and core
            assertEqual(sum(enzymeDNAFootprints([m.enzymeIndexs_core; m.enzymeIndexs_betaClamp])), enzymeDNAFootprints(m.enzymeIndexs_coreBetaClampPrimase));
            
            %assert footprint values and strandedness
            assertEqual(enzymeDNAFootprints, m.enzymeDNAFootprints);
            assertEqual(monomerDNAFootprintBindingStrandedness, c.monomerDNAFootprintBindingStrandedness(m.enzymeMonomerGlobalIndexs));
            assertEqual(complexDNAFootprintBindingStrandedness', c.complexDNAFootprintBindingStrandedness(m.enzymeComplexGlobalIndexs)');
            assertEqual(monomerDNAFootprintRegionStrandedness, c.monomerDNAFootprintRegionStrandedness(m.enzymeMonomerGlobalIndexs));
            assertEqual(complexDNAFootprintRegionStrandedness', c.complexDNAFootprintRegionStrandedness(m.enzymeComplexGlobalIndexs)');
            
            %assert
            assertEqual(...
                sortrows([  reactionBoundMonomer   reactionBoundComplex   reactionMonomerCatalysisMatrix   reactionComplexCatalysisMatrix], 1:2+numel(m.enzymes)), ...
                sortrows([c.reactionBoundMonomer c.reactionBoundComplex c.reactionMonomerCatalysisMatrix c.reactionComplexCatalysisMatrix], 1:2+numel(m.enzymes)));
            
            %check m.dnaAFunctionalBoxStartPositions is sorted
            assertTrue(issorted(m.dnaAFunctionalBoxStartPositions));
            
            %check sorted so that ismembc can be used in place of ismember
            assertTrue(issorted(m.enzymeGlobalIndexs([m.enzymeIndexs_2coreBetaClampGammaComplexPrimase; m.enzymeIndexs_coreBetaClampGammaComplex; m.enzymeIndexs_helicase])));
            
            %check that maximum okazaki fragment length less than twice
            %mean okazaki fragment length
            assertIn(max(abs([diff(m.primaseBindingLocations{1}); diff(m.primaseBindingLocations{2})])), [0 2*m.okazakiFragmentMeanLength-1]);
        end
        
        function testCompleteReplication_SmallChromosome(this)
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            
            m = this.process;
            c = m.chromosome;
            
            m.okazakiFragmentMeanLength         = 500;
            m.laggingBackupClampReloadingLength = 250;
            
            bases = 'ACGT';
            seq = bases(randi(4, 2*6*m.okazakiFragmentMeanLength, 1));
            c.sequence = ChromosomeSequence(seq);
            c.sequenceLen = length(seq);
            c.sequenceGCContent = getGCContent(c.sequence);
            c.terCPosition = c.sequenceLen / 2;
            m.terCPosition = c.sequenceLen / 2;
            m.primaseBindingLocations = {
                size(c.sequence, 1) + 1 - m.calculatePrimaseBindingLocations(size(c.sequence, 1) - m.terCPosition)
                m.calculatePrimaseBindingLocations(m.terCPosition)
                }';
            m.dnaAFunctionalBoxStartPositions = c.sequenceLen - (40 + 10*(1:5)');
            
            %% initial state
            c.initialize();
            c.linkingNumbers(1, 1:2) = 0;
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_dntp) = [
                sum(sum(c.sequence(:, 2:3) == 'A'));
                sum(sum(c.sequence(:, 2:3) == 'C'));
                sum(sum(c.sequence(:, 2:3) == 'G'));
                sum(sum(c.sequence(:, 2:3) == 'T'))];
            m.substrates(m.substrateIndexs_atp) = ...
                + size(c.sequence, 1) ...                              %helicase unwinding each base pair
                + sum(cellfun(@numel, m.primaseBindingLocations)) ...  %lagging strand beta-clamp formation
                + 2;                                                   %leading strand beta-clamp formation
            m.substrates(m.substrateIndexs_water) = m.substrates(m.substrateIndexs_atp);
            m.substrates(m.substrateIndexs_nad) = sum(cellfun(@numel,m.primaseBindingLocations)) + 2;
            
            m.enzymes = m.enzymeComposition(:, [m.enzymeIndexs_replisome; m.enzymeIndexs_betaClamp; m.enzymeIndexs_ssb8mer]) * [2; 2; 100];
            m.enzymes(m.enzymeIndexs_ligase) = 100;
            m.boundEnzymes(:) = 0;
            
            this.setUpDnaAOriComplex();
            
            %% evolve state
            for i = 1:60
                m.evolveState();
                if all(m.strandDuplicated)
                    break;
                end
            end
            
            %% assert final state
            chrLen = c.sequenceLen;
            assertEqual(CircularSparseMat([ones(4, 1) (1:4)'], repmat(chrLen, 4, 1), [chrLen 4], 1), c.polymerizedRegions);
            assertEqual(CircularSparseMat([], [], [chrLen 4], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [chrLen 4], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [chrLen 4], 1), c.complexBoundSites);
            
            final_substrates = zeros(size(m.substrates));
            final_substrates(m.substrateIndexs_diphosphate) = 2 * size(c.sequence, 1);
            final_substrates(m.substrateIndexs_adp) = size(c.sequence, 1) + sum(cellfun(@numel,m.primaseBindingLocations)) + 2;
            final_substrates(m.substrateIndexs_phosphate) = m.substrates(m.substrateIndexs_adp);
            final_substrates(m.substrateIndexs_nmn) =  sum(cellfun(@numel,m.primaseBindingLocations)) + 2;
            final_substrates(m.substrateIndexs_amp) = m.substrates(m.substrateIndexs_nmn);
            final_substrates(m.substrateIndexs_hydrogen) = m.substrates(m.substrateIndexs_adp) + m.substrates(m.substrateIndexs_nmn);
            
            final_enzymes = m.enzymeComposition(:, [m.enzymeIndexs_replisome; m.enzymeIndexs_betaClamp; m.enzymeIndexs_ssb8mer]) * [2; 2; 100];
            final_enzymes(m.enzymeIndexs_ligase) = 100;
            
            assertEqual(final_substrates', m.substrates');
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(zeros(size(m.enzymes)), m.boundEnzymes);
            
            assertFalse(m.isDnaAORIComplexAssembled);
            assertEqual([false false], m.leadingStrandElongating);
            assertEqual([false false], m.laggingStrandElongating);
            assertEqual([true true], m.leadingStrandPolymerized);
            assertEqual([true true], m.laggingStrandPolymerized);
            assertEqual([true true], m.strandPolymerized);
            assertEqual(cellfun(@numel,m.primaseBindingLocations) + 1, m.numLigations);
            assertEqual([true true], m.strandLigated);
            assertEqual([true true], m.strandDuplicated);
            
            assertEqual([0 0], m.helicasePosition);
            assertEqual([0 0], m.leadingPolymerasePosition);
            assertEqual([0 0], m.laggingPolymerasePosition);
            assertEqual([0 0], m.leadingPosition);
            assertEqual([0 0], m.laggingPosition);
            assertEqual([0 0], m.laggingBackupBetaClampPosition);
            
            assertEqual([0 0], m.okazakiFragmentIndex);
            assertEqual([0 0], m.okazakiFragmentPosition);
            assertEqual([0 0], m.okazakiFragmentLength);
            assertEqual([0 0], m.okazakiFragmentProgress);
            
            assertEqual(CircularSparseMat([], [], [chrLen 2], 1), m.leadingStrandBoundSSBs);
            assertEqual(CircularSparseMat([], [], [chrLen 2], 1), m.laggingStrandBoundSSBs);
            assertEqual([0 0], m.numLeadingTemplateBoundSSBs);
            assertEqual([0 0], m.numLaggingTemplateBoundSSBs);
            assertEqual([false false], m.areLaggingStrandSSBSitesBound);
        end
        
        function testCompleteReplicationInSteps(this)
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            
            m = this.process;
            c = m.chromosome;
            
            m.okazakiFragmentMeanLength         = 500;
            m.laggingBackupClampReloadingLength = 250;
            
            bases = 'ACGT';
            seq = bases(randi(4, 2*6*m.okazakiFragmentMeanLength, 1));
            c.sequence = ChromosomeSequence(seq);
            c.sequenceLen = length(seq);
            c.sequenceGCContent = getGCContent(c.sequence);
            c.terCPosition = c.sequenceLen / 2;
            m.terCPosition = c.sequenceLen / 2;
            m.primaseBindingLocations = {
                size(c.sequence, 1) + 1 - m.calculatePrimaseBindingLocations(size(c.sequence, 1) - m.terCPosition)
                m.calculatePrimaseBindingLocations(m.terCPosition)
                }';
            m.dnaAFunctionalBoxStartPositions = c.sequenceLen - (40 + 10*(1:5)');
            
            helGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            ssbGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_ssb8mer);
            
            chrLen = size(m.chromosome.sequence, 1);
            polRate = m.dnaPolymeraseElongationRate;
            
            helFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            corFtpt = m.enzymeDNAFootprints(m.enzymeIndexs_core);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            corFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_core);
            bClmpFtpt = m.enzymeDNAFootprints(m.enzymeIndexs_betaClamp);
            
            %% initial state
            c.initialize();
            c.linkingNumbers(1, 1:2) = 0;
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_dntp) = [
                sum(sum(c.sequence(:, 2:3) == 'A'));
                sum(sum(c.sequence(:, 2:3) == 'C'));
                sum(sum(c.sequence(:, 2:3) == 'G'));
                sum(sum(c.sequence(:, 2:3) == 'T'))];
            m.substrates(m.substrateIndexs_atp) = ...
                + size(c.sequence, 1) ...                              %helicase unwinding each base pair
                + sum(cellfun(@numel, m.primaseBindingLocations)) ...  %lagging strand beta-clamp formation
                + 2;                                                   %leading strand beta-clamp formation
            m.substrates(m.substrateIndexs_water) = m.substrates(m.substrateIndexs_atp);
            m.substrates(m.substrateIndexs_nad) = sum(cellfun(@numel,m.primaseBindingLocations)) + 2;
            
            m.enzymes = m.enzymeComposition(:, [m.enzymeIndexs_replisome; m.enzymeIndexs_betaClamp; m.enzymeIndexs_ssb8mer]) * [2; 2; 100];
            m.enzymes(m.enzymeIndexs_ligase) = 100;
            m.boundEnzymes(:) = 0;
            
            this.setUpDnaAOriComplex();
            
            %% initiate
            %store current state
            tmpSubstrates = m.substrates;
            tmpEnzymes = m.enzymes;
            tmpBoundEnzymes = m.boundEnzymes;
            
            %evolve state
            m.initiateReplication();
            
            %assertions
            uwdLen = helFtpt3 + corFtpt5 + 1;
            substrates = tmpSubstrates;
            substrates([m.substrateIndexs_atp;m.substrateIndexs_water]) = ...
                substrates([m.substrateIndexs_atp;m.substrateIndexs_water]) - ...
                2*(uwdLen+1);
            substrates([m.substrateIndexs_adp;m.substrateIndexs_phosphate;m.substrateIndexs_hydrogen]) = ...
                substrates([m.substrateIndexs_adp;m.substrateIndexs_phosphate;m.substrateIndexs_hydrogen]) + ...
                2*(uwdLen+1);
            assertEqual(substrates, m.substrates);
            
            enzymes = tmpEnzymes;
            enzymes(m.enzymeIndexs_helicase) = enzymes(m.enzymeIndexs_helicase) - 2;
            enzymes = enzymes - 2 * m.enzymeComposition(:, m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            boundEnzymes = tmpBoundEnzymes;
            boundEnzymes(m.enzymeIndexs_helicase) = boundEnzymes(m.enzymeIndexs_helicase) + 2;
            boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) + 2;
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            
            assertEqual(CircularSparseMat([1 1;uwdLen+1 2;1 4; chrLen-uwdLen+1 4], [chrLen; chrLen-2*uwdLen; uwdLen; uwdLen], [chrLen 4], 1), c.polymerizedRegions);
            posStrnds = [
                m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R1234) ones(4, 1);
                m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R5) 1;
                chrLen-helFtpt-corFtpt5 m.leadingStrandIndexs(1); 1+corFtpt5+1 m.leadingStrandIndexs(2);
                chrLen-corFtpt5 m.leadingStrandIndexs(1); 1-(holFtpt-corFtpt5)+1 m.leadingStrandIndexs(2)];
            prots = [
                repmat(m.complexIndexs_DnaA_7mer_ATP, 4, 1);
                m.complexIndexs_DnaA_1mer_ATP;
                helGblIdx; helGblIdx;
                polGblIdx; polGblIdx];
            assertEqual(CircularSparseMat(posStrnds, prots, [chrLen 4], 1), c.complexBoundSites);
            
            assertTrue(m.isDnaAORIComplexAssembled);
            assertEqual([chrLen 1], m.leadingPosition);
            assertEqual([0 0], m.laggingPosition);
            assertEqual([0 0], m.laggingBackupBetaClampPosition);
            assertEqual([-helFtpt holFtpt-size(c.sequence,1)], m.helicasePosition-m.leadingPolymerasePosition);
            
            %% unwind and polymerize
            for i = 1:35
                m.unwindAndPolymerizeDNA();
                m.freeAndBindSSBs();
                
                assertEqual([chrLen-helFtpt holFtpt], mod(m.helicasePosition-m.leadingPolymerasePosition-1, chrLen)+1);
            end
            
            assertFalse(m.isDnaAORIComplexAssembled);
            
            helicasePos = m.helicasePosition;
            assertTrue(chrLen - (helicasePos(1) + helFtpt) >= 2 * m.okazakiFragmentMeanLength);
            assertTrue(chrLen - (helicasePos(1) + helFtpt) < 2 * m.okazakiFragmentMeanLength + polRate);
            assertTrue(helicasePos(2) - 1 >= 2 * m.okazakiFragmentMeanLength);
            assertTrue(helicasePos(2) - 1 < 2 * m.okazakiFragmentMeanLength + polRate);
            
            assertEqual([0 0], m.laggingPosition);
            assertEqual([0 0], m.laggingBackupBetaClampPosition);
            
            %% initiate first Okazaki fragment
            %bind beta-clamp
            tmpSubstrates = m.substrates;
            tmpEnzymes = m.enzymes;
            tmpBoundEnzymes = m.boundEnzymes;
            
            m.initiateOkazakiFragment();
            
            substrates = tmpSubstrates;
            substrates([m.substrateIndexs_atp;m.substrateIndexs_water]) = ...
                substrates([m.substrateIndexs_atp;m.substrateIndexs_water]) - ...
                2;
            substrates([m.substrateIndexs_adp;m.substrateIndexs_phosphate;m.substrateIndexs_hydrogen]) = ...
                substrates([m.substrateIndexs_adp;m.substrateIndexs_phosphate;m.substrateIndexs_hydrogen]) + ...
                2;
            assertEqual(substrates, m.substrates);
            
            enzymes = tmpEnzymes;
            enzymes(m.enzymeIndexs_betaClampMonomer) = enzymes(m.enzymeIndexs_betaClampMonomer) - 4;
            boundEnzymes = tmpBoundEnzymes;
            boundEnzymes(m.enzymeIndexs_betaClamp) = boundEnzymes(m.enzymeIndexs_betaClamp) + 2;
            assertEqual(enzymes(setdiff(1:end,[m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])), m.enzymes(setdiff(1:end,[m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])));
            assertEqual(boundEnzymes(setdiff(1:end,m.enzymeIndexs_ssb8mer)), m.boundEnzymes(setdiff(1:end, m.enzymeIndexs_ssb8mer)));
            assertTrue(all(m.enzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer]) >= 0));
            assertTrue(all(m.boundEnzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer]) >= 0));
            assertEqual(m.boundEnzymes(m.enzymeIndexs_ssb8mer), nnz(c.complexBoundSites == ssbGblIdx));
            assertEqual(200, m.enzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])' * [1;2] + m.boundEnzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])' * [1;2]);
            
            assertEqual([0 0], m.laggingPosition);
            assertEqual([m.primaseBindingLocations{1}(1)-bClmpFtpt-corFtpt3 m.primaseBindingLocations{2}(1)+corFtpt3+1], m.laggingBackupBetaClampPosition);
            
            %bind polymerase core to lagging strand
            %synthesize lagging strand primer
            %synthesize more leading strand
            tmpSubstrates = m.substrates;
            tmpEnzymes = m.enzymes;
            tmpBoundEnzymes = m.boundEnzymes;
            tmpLeadingPosition = m.leadingPosition;
            
            m.unwindAndPolymerizeDNA();
            m.dissociateFreeSSBComplexes();
            m.freeAndBindSSBs();
            
            substrates = tmpSubstrates;
            substrates(m.substrateIndexs_diphosphate) = ...
                substrates(m.substrateIndexs_diphosphate) ...
                + 2 * m.primerLength ...
                + 2 * polRate;
            substrates(m.substrateIndexs_dntp) = ...
                substrates(m.substrateIndexs_dntp) ...
                - c.sequence.subsequenceBaseCounts(m.primaseBindingLocations{1}(1) + (0:m.primerLength-1), 3) ...
                - c.sequence.subsequenceBaseCounts(m.primaseBindingLocations{2}(1) - (0:m.primerLength-1), 2) ...
                - c.sequence.subsequenceBaseCounts(tmpLeadingPosition(1) - (0:polRate-1), 2) ...
                - c.sequence.subsequenceBaseCounts(tmpLeadingPosition(2) + (0:polRate-1), 3);
            substrates([m.substrateIndexs_atp;m.substrateIndexs_water]) = ...
                substrates([m.substrateIndexs_atp;m.substrateIndexs_water]) - ...
                2*polRate;
            substrates([m.substrateIndexs_adp;m.substrateIndexs_phosphate;m.substrateIndexs_hydrogen]) = ...
                substrates([m.substrateIndexs_adp;m.substrateIndexs_phosphate;m.substrateIndexs_hydrogen]) + ...
                2*polRate;
            assertEqual(substrates, m.substrates);
            
            enzymes = tmpEnzymes;
            boundEnzymes = tmpBoundEnzymes;
            boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) - 2;
            boundEnzymes(m.enzymeIndexs_betaClamp) = boundEnzymes(m.enzymeIndexs_betaClamp) - 2;
            boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) + 2;
            boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) + 2;
            assertEqual(enzymes(setdiff(1:end,[m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])), m.enzymes(setdiff(1:end,[m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])));
            assertEqual(boundEnzymes(setdiff(1:end,m.enzymeIndexs_ssb8mer)), m.boundEnzymes(setdiff(1:end, m.enzymeIndexs_ssb8mer)));
            assertTrue(all(m.enzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer]) >= 0));
            assertTrue(all(m.boundEnzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer]) >= 0));
            assertEqual(m.boundEnzymes(m.enzymeIndexs_ssb8mer), nnz(c.complexBoundSites == ssbGblIdx));
            assertEqual(200, m.enzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])' * [1;2] + m.boundEnzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])' * [1;2]);
            
            assertEqual([tmpLeadingPosition(1)-polRate tmpLeadingPosition(2)+polRate], m.leadingPosition);
            assertEqual([chrLen-helFtpt holFtpt], mod(m.helicasePosition-m.leadingPolymerasePosition-1, chrLen)+1);
            assertEqual([m.primaseBindingLocations{1}(1)+m.primerLength m.primaseBindingLocations{2}(1)-m.primerLength], m.laggingPosition);
            assertEqual([0 0], m.laggingBackupBetaClampPosition);
            
            %% unwind and polymerize leading strand and first lagging Okazaki fragment
            for i = 1:25
                m.unwindAndPolymerizeDNA();
                m.dissociateFreeSSBComplexes();
                m.freeAndBindSSBs();
                
                assertEqual([chrLen-helFtpt holFtpt], mod(m.helicasePosition-m.leadingPolymerasePosition-1, chrLen)+1);
            end
            
            helicasePos = m.helicasePosition;
            assertTrue(m.primaseBindingLocations{1}(1) - (helicasePos(1) + helFtpt) >= 2 * m.okazakiFragmentMeanLength);
            assertTrue(m.primaseBindingLocations{1}(1) - (helicasePos(1) + helFtpt) < 2 * m.okazakiFragmentMeanLength + polRate);
            assertTrue(helicasePos(2) - m.primaseBindingLocations{2}(1) >= 2 * m.okazakiFragmentMeanLength);
            assertTrue(helicasePos(2) - m.primaseBindingLocations{2}(1) < 2 * m.okazakiFragmentMeanLength + polRate);
            
            assertEqual([1 chrLen], m.laggingPosition);
            assertEqual([1 1], m.okazakiFragmentIndex);
            assertEqual(m.okazakiFragmentLength, m.okazakiFragmentProgress);
            assertEqual([0 0], m.laggingBackupBetaClampPosition);
            
            %% initiate second Okazaki fragment
            m.initiateOkazakiFragment();
            
            enzymes = zeros(size(m.enzymes));
            enzymes(m.enzymeIndexs_ligase) = 100;
            boundEnzymes = zeros(size(m.enzymes));
            boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = 2;
            boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = 2;
            boundEnzymes(m.enzymeIndexs_betaClamp) = 2;
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            assertEqual(enzymes(setdiff(1:end,[m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])), m.enzymes(setdiff(1:end,[m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])));
            assertEqual(boundEnzymes(setdiff(1:end,m.enzymeIndexs_ssb8mer)), m.boundEnzymes(setdiff(1:end, m.enzymeIndexs_ssb8mer)));
            assertTrue(all(m.enzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer]) >= 0));
            assertTrue(all(m.boundEnzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer]) >= 0));
            assertEqual(m.boundEnzymes(m.enzymeIndexs_ssb8mer), nnz(c.complexBoundSites == ssbGblIdx));
            assertEqual(200, m.enzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])' * [1;2] + m.boundEnzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])' * [1;2]);
            
            assertEqual([1 chrLen], m.laggingPosition);
            assertEqual([1 1], m.okazakiFragmentIndex);
            assertEqual(m.okazakiFragmentLength, m.okazakiFragmentProgress);
            assertEqual([m.primaseBindingLocations{1}(2)-bClmpFtpt-corFtpt3 m.primaseBindingLocations{2}(2)+corFtpt3+1], m.laggingBackupBetaClampPosition);
            
            %% terminate first Okazaki fragment
            m.terminateOkazakiFragment();
            
            enzymes = zeros(size(m.enzymes));
            enzymes(m.enzymeIndexs_ligase) = 100;
            enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            boundEnzymes = zeros(size(m.enzymes));
            boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = 2;
            boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = 2;
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            assertEqual(enzymes(setdiff(1:end, [m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])), m.enzymes(setdiff(1:end,[m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])));
            assertEqual(boundEnzymes(setdiff(1:end,m.enzymeIndexs_ssb8mer)), m.boundEnzymes(setdiff(1:end, m.enzymeIndexs_ssb8mer)));
            assertTrue(all(m.enzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer]) >= 0));
            assertTrue(all(m.boundEnzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer]) >= 0));
            assertEqual(m.boundEnzymes(m.enzymeIndexs_ssb8mer), nnz(c.complexBoundSites == ssbGblIdx));
            assertEqual(200, m.enzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])' * [1;2] + m.boundEnzymes([m.enzymeIndexs_ssb4mer; m.enzymeIndexs_ssb8mer])' * [1;2]);
            
            assertEqual([m.primaseBindingLocations{1}(2)-bClmpFtpt-corFtpt3 m.primaseBindingLocations{2}(2)+corFtpt3+1-corFtpt], m.laggingPolymerasePosition);
            assertEqual([m.primaseBindingLocations{1}(2) m.primaseBindingLocations{2}(2)], m.laggingPosition);
            assertEqual([2 2], m.okazakiFragmentIndex);
            assertEqual([0 0], m.laggingBackupBetaClampPosition);
            
            assertEqual(CircularSparseMat([chrLen 2; chrLen 3], [1;1], [chrLen 4], 1), c.strandBreaks);
            assertEqual(0, nnz(c.complexBoundSites(m.primaseBindingLocations{1}(1):chrLen,1) == ssbGblIdx))
            assertEqual(0, nnz(c.complexBoundSites(1:m.primaseBindingLocations{2}(1),4) == ssbGblIdx))
            
            %% ligate first Okazaki fragments and leading strand
            tmpSubstrates = m.substrates;
            
            m.ligateDNA();
            
            tmpSubstrates(m.substrateIndexs_nad) = tmpSubstrates(m.substrateIndexs_nad) - 2;
            tmpSubstrates(m.substrateIndexs_nmn) = tmpSubstrates(m.substrateIndexs_nmn) + 2;
            tmpSubstrates(m.substrateIndexs_amp) = tmpSubstrates(m.substrateIndexs_amp) + 2;
            tmpSubstrates(m.substrateIndexs_hydrogen) = tmpSubstrates(m.substrateIndexs_hydrogen) + 2;
            assertEqual(tmpSubstrates', m.substrates');
            
            assertEqual([1 1], m.numLigations);
            assertEqual(CircularSparseMat([], [], [chrLen 4], 1), c.strandBreaks);
            
            %% unwind and polymerize
            nOkazakiFragments = cellfun(@numel, m.primaseBindingLocations);
            for i = 1:35
                m.unwindAndPolymerizeDNA();
                m.freeAndBindSSBs();
                m.dissociateFreeSSBComplexes();
                m.initiateOkazakiFragment();
                if ~all(nOkazakiFragments == m.okazakiFragmentIndex)
                    m.terminateOkazakiFragment();
                end
                m.ligateDNA();
                
                helicasePos = m.helicasePosition;
                okazakiFragmentPosition = m.okazakiFragmentPosition;
                assertTrue(okazakiFragmentPosition(1) - (helicasePos(1) + helFtpt) < 2 * m.okazakiFragmentMeanLength + polRate);
                if okazakiFragmentPosition ~= 0
                    assertTrue(helicasePos(2) - okazakiFragmentPosition(2) < 2 * m.okazakiFragmentMeanLength + polRate);
                end
                assertEqual([chrLen-helFtpt holFtpt], mod(helicasePos - m.leadingPolymerasePosition-1, chrLen) + 1);
            end
            
            enzymes = zeros(size(m.enzymes));
            enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            enzymes(m.enzymeIndexs_ligase) = 100;
            enzymes(m.enzymeIndexs_ssb4mer) = 200;
            boundEnzymes = zeros(size(enzymes));
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = 2;
            boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = 2;
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            
            assertEqual([m.terCPosition-(corFtpt5+bClmpFtpt+corFtpt5+helFtpt3) m.terCPosition+1+(corFtpt5+bClmpFtpt+corFtpt5+helFtpt3)], m.leadingPosition);
            assertEqual([m.primaseBindingLocations{1}(end-1) m.primaseBindingLocations{2}(end-1)], m.laggingPosition);
            assertEqual(nOkazakiFragments, m.okazakiFragmentIndex);
            assertEqual(m.okazakiFragmentLength, m.okazakiFragmentProgress);
            assertEqual([0 0], m.laggingBackupBetaClampPosition);
            assertEqual([0 0], m.numLeadingTemplateBoundSSBs);
            assertEqual([0 0], m.numLaggingTemplateBoundSSBs);
            
            assertEqual(CircularSparseMat([ones(4, 1) (1:4)'], repmat(chrLen, 4, 1), [chrLen 4], 1), c.polymerizedRegions);
            assertEqual(CircularSparseMat([], [], [chrLen 4], 1), c.strandBreaks);
            
            %% terminate last Okazaki fragment
            m.terminateOkazakiFragment();
            
            enzymes = zeros(size(m.enzymes));
            enzymes(m.enzymeIndexs_betaClampMonomer) = 8;
            enzymes(m.enzymeIndexs_ligase) = 100;
            enzymes(m.enzymeIndexs_ssb4mer) = 200;
            boundEnzymes = zeros(size(enzymes));
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            
            assertEqual(CircularSparseMat([ones(4, 1) (1:4)'], repmat(chrLen, 4, 1), [chrLen 4], 1), c.polymerizedRegions);
            assertEqual(CircularSparseMat([m.primaseBindingLocations{1}(end-1)-1 3; m.primaseBindingLocations{2}(end-1) 2], [1;1],[chrLen 4],1), c.strandBreaks);
            
            assertEqual([m.terCPosition-(corFtpt5+bClmpFtpt+corFtpt5+helFtpt3) m.terCPosition+1+(corFtpt5+bClmpFtpt+corFtpt5+helFtpt3)], m.leadingPosition);
            assertEqual([0 0], m.laggingPosition);
            assertEqual([0 0], m.okazakiFragmentIndex);
            assertEqual([0 0], m.laggingBackupBetaClampPosition);
            assertEqual([0 0], m.numLaggingTemplateBoundSSBs);
            
            %% ligate last Okazaki fragment
            m.ligateDNA();
            
            assertEqual(CircularSparseMat([], [], [chrLen 4], 1), c.strandBreaks);
            
            %% terminate replication
            m.terminateReplication();
            
            assertEqual(CircularSparseMat([m.terCPosition 2; m.terCPosition 3], [1;1],[chrLen 4],1), c.strandBreaks);
            
            %% final ligation
            m.ligateDNA();
            
            %% assert final state
            assertEqual(CircularSparseMat([ones(4, 1) (1:4)'], repmat(chrLen, 4, 1), [chrLen 4], 1), c.polymerizedRegions);
            assertEqual(CircularSparseMat([], [], [chrLen 4], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [chrLen 4], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [chrLen 4], 1), c.complexBoundSites);
            
            final_substrates = zeros(size(m.substrates));
            final_substrates(m.substrateIndexs_diphosphate) = 2 * size(c.sequence, 1);
            final_substrates(m.substrateIndexs_adp) = size(c.sequence, 1) + sum(cellfun(@numel,m.primaseBindingLocations)) + 2;
            final_substrates(m.substrateIndexs_phosphate) = m.substrates(m.substrateIndexs_adp);
            final_substrates(m.substrateIndexs_nmn) =  sum(cellfun(@numel,m.primaseBindingLocations)) + 2;
            final_substrates(m.substrateIndexs_amp) = m.substrates(m.substrateIndexs_nmn);
            final_substrates(m.substrateIndexs_hydrogen) = m.substrates(m.substrateIndexs_adp) + m.substrates(m.substrateIndexs_nmn);
            
            final_enzymes = m.enzymeComposition(:, [m.enzymeIndexs_replisome; m.enzymeIndexs_betaClamp; m.enzymeIndexs_ssb8mer]) * [2; 2; 100];
            final_enzymes(m.enzymeIndexs_ligase) = 100;
            
            assertEqual(final_substrates', m.substrates');
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(zeros(size(m.enzymes)), m.boundEnzymes);
            
            assertEqual(false, m.isDnaAORIComplexAssembled);
            assertEqual([false false], m.leadingStrandElongating);
            assertEqual([false false], m.laggingStrandElongating);
            assertEqual([true true], m.leadingStrandPolymerized);
            assertEqual([true true], m.laggingStrandPolymerized);
            assertEqual([true true], m.strandPolymerized);
            assertEqual(cellfun(@numel,m.primaseBindingLocations) + 1, m.numLigations);
            assertEqual([true true], m.strandLigated);
            assertEqual([true true], m.strandDuplicated);
            
            assertEqual([0 0], m.helicasePosition);
            assertEqual([0 0], m.leadingPolymerasePosition);
            assertEqual([0 0], m.laggingPolymerasePosition);
            assertEqual([0 0], m.leadingPosition);
            assertEqual([0 0], m.laggingPosition);
            assertEqual([0 0], m.laggingBackupBetaClampPosition);
            
            assertEqual([0 0], m.okazakiFragmentIndex);
            assertEqual([0 0], m.okazakiFragmentPosition);
            assertEqual([0 0], m.okazakiFragmentLength);
            assertEqual([0 0], m.okazakiFragmentProgress);
            
            assertEqual(CircularSparseMat([], [], [chrLen 2], 1), m.leadingStrandBoundSSBs);
            assertEqual(CircularSparseMat([], [], [chrLen 2], 1), m.laggingStrandBoundSSBs);
            assertEqual([0 0], m.numLeadingTemplateBoundSSBs);
            assertEqual([0 0], m.numLaggingTemplateBoundSSBs);
            assertEqual([false false], m.areLaggingStrandSSBSitesBound);
        end
        
        function testReplicationDuration(this)
            m = this.process;
            c = m.chromosome;
            s = edu.stanford.covert.cell.sim.SimulationFixture.load();
            t = s.state('Time');
            
            assertIn((c.sequenceLen / 2) / m.okazakiFragmentMeanLength * (1 + 1 + (m.okazakiFragmentMeanLength-m.primerLength)/m.dnaPolymeraseElongationRate), ...
                [0 t.replicationDuration]);
            
            assertEqual(...
                t.cellCycleLength - t.replicationInitiationDuration - t.cytokinesisDuration, ...
                t.replicationDuration);
        end
        
        function testReplicationWithFewSSBs(this)
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            
            m = this.process;
            c = m.chromosome;
            
            m.okazakiFragmentMeanLength         = 500;
            m.laggingBackupClampReloadingLength = 250;
            
            bases = 'ACGT';
            seq = bases(randi(4, 2*6*m.okazakiFragmentMeanLength, 1));
            c.sequence = ChromosomeSequence(seq);
            c.sequenceLen = length(seq);
            c.sequenceGCContent = getGCContent(c.sequence);
            c.terCPosition = c.sequenceLen / 2;
            m.terCPosition = c.sequenceLen / 2;
            m.primaseBindingLocations = {
                [5578 5102 4576 4050 3523 3001]', ...
                [ 532 1054 1562 2040 2520 3000]'};
            m.dnaAFunctionalBoxStartPositions = c.sequenceLen - (40 + 10*(1:5)');
            
            %% initial state such that more sum of SSBs needed to proceed at both replisomes greater that number of SSBs,
            %sufficient SSBs to proceed either replisome individually
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_dntp) = [
                sum(sum(c.sequence(:, 2:3) == 'A'));
                sum(sum(c.sequence(:, 2:3) == 'C'));
                sum(sum(c.sequence(:, 2:3) == 'G'));
                sum(sum(c.sequence(:, 2:3) == 'T'))];
            m.substrates(m.substrateIndexs_atp) = ...
                + size(c.sequence, 1) ...                              %helicase unwinding each base pair
                + sum(cellfun(@numel, m.primaseBindingLocations)) ...  %lagging strand beta-clamp formation
                + 2;                                                   %leading strand beta-clamp formation
            m.substrates(m.substrateIndexs_water) = m.substrates(m.substrateIndexs_atp);
            m.substrates(m.substrateIndexs_nad) = sum(cellfun(@numel,m.primaseBindingLocations)) + 2;
            
            chrLen = c.sequenceLen;
            c.initialize();
            initPolymerizedRegions = CircularSparseMat([
                1 1 2634 3590 1 4576 1 3568
                1 2 2 2 3 3 4 4
                ]', [6000 1562 934 2411 2611 1425 2633 2433]', [chrLen 4], 1);
            c.polymerizedRegions = initPolymerizedRegions;
            c.linkingNumbers = CircularSparseMat([
                1 3590 1 3590 1 4576 1 4576
                1 1 2 2 3 3 4 4
                ]', [
                148.7630 229.6208 148.7630 229.6208 248.6686 135.7153 248.6686 135.7153
                ]', [chrLen 4], 1);
            c.monomerBoundSites = CircularSparseMat([], [], [chrLen 4], 1);
            c.complexBoundSites = CircularSparseMat([
                2263 3557 3577 1042 2052 4014 5066 2576 2625 3731
                1 1 1 2 2 3 3 4 4 4
                ]', [
                48 50  4  5 10 10  5  4 50 48
                ]', [chrLen 4], 1);
            c.strandBreaks = CircularSparseMat([], [], [chrLen 4], 1);
            c.damagedBases = CircularSparseMat([], [], [chrLen 4], 1);
            m.enzymes = [0 0 0 0 0 0 0 0 0 100 0 0 0]';
            m.enzymes(m.enzymeIndexs_ssb4mer) = 2;
            m.boundEnzymes = [0 0 2 2 0 0 0 2 2 0 0 0 2]';
            
            %% evolve state with no SSB dissociation and check that replication doesn't progress
            m.ssbDissociationRate = 0;
            for i = 1:600
                m.evolveState();
                if all(m.strandDuplicated)
                    break;
                end
            end
            assertEqual(initPolymerizedRegions, c.polymerizedRegions);
            
            %% increase SSB dissociation rate, evolve state and check that replication finishes
            m.ssbDissociationRate = 2;
            for i = 1:100
                m.evolveState();
                if all(m.strandDuplicated)
                    break;
                end
            end
            
            assertEqual(CircularSparseMat([ones(4, 1) (1:4)'], repmat(chrLen, 4, 1), [chrLen 4], 1), c.polymerizedRegions);
            assertEqual([true true], m.strandLigated);
            assertEqual([true true], m.strandDuplicated);
        end
    end
    
    %test initiateReplication
    methods
        function testInitiateReplication(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            helGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            chrLen = size(m.chromosome.sequence, 1);
            helFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            uwdLen = helFtpt3 + corFtpt5 + 1;
            
            c.initialize();
            this.setUpDnaAOriComplex();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2*(uwdLen+1);
            m.substrates(m.substrateIndexs_water) = 2*(uwdLen+1);
            m.enzymes(:) = 2 * m.enzymeComposition(:, m.enzymeIndexs_replisome);
            m.boundEnzymes(:) = 0;
            
            m.initiateReplication();
            
            substrates = zeros(size(m.substrates));
            substrates([m.substrateIndexs_adp;m.substrateIndexs_phosphate;m.substrateIndexs_hydrogen]) = 2*(uwdLen+1);
            assertEqual(substrates, m.substrates);
            
            enzymes = 2 * m.enzymeComposition(:, m.enzymeIndexs_replisome);
            enzymes(m.enzymeIndexs_helicase) = 0;
            enzymes = enzymes - 2 * m.enzymeComposition(:, m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            boundEnzymes = zeros(size(m.enzymes));
            boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            
            posStrnds = [
                m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R1234) ones(4, 1);
                m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R5) 1;
                chrLen-helFtpt-corFtpt5 1; 1+corFtpt5+1 4;
                chrLen-corFtpt5 1; 1-(holFtpt-corFtpt5)+1 4];
            prots = [
                repmat(m.complexIndexs_DnaA_7mer_ATP, 4, 1);
                m.complexIndexs_DnaA_1mer_ATP;
                helGblIdx; helGblIdx;
                polGblIdx; polGblIdx];
            assertEqual(CircularSparseMat(posStrnds, prots, [chrLen 4], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([1 1;uwdLen+1 2;1 4; chrLen-uwdLen+1 4], [chrLen; chrLen-2*uwdLen; uwdLen; uwdLen], [chrLen 4], 1), c.polymerizedRegions);
            
            assertTrue(m.isDnaAORIComplexAssembled);
            assertEqual([chrLen 1], m.leadingPosition);
            assertEqual([-helFtpt holFtpt-size(c.sequence,1)], m.helicasePosition-m.leadingPolymerasePosition);
        end
        
        function testInitiateReplication_NoATP(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            chrLen = size(m.chromosome.sequence, 1);
            helFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            uwdLen = helFtpt3 + corFtpt5 + 1;
            
            c.initialize();
            this.setUpDnaAOriComplex();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 0;
            m.substrates(m.substrateIndexs_water) = 2*(uwdLen+1);
            m.enzymes(:) = 2 * m.enzymeComposition(:, m.enzymeIndexs_replisome);
            m.boundEnzymes(:) = 0;
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            
            m.initiateReplication();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            
            posStrnds = [
                m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R1234) ones(4, 1);
                m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R5) 1];
            prots = [
                repmat(m.complexIndexs_DnaA_7mer_ATP, 4, 1);
                m.complexIndexs_DnaA_1mer_ATP];
            assertEqual(CircularSparseMat(posStrnds, prots, [chrLen 4], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([1 1;1 2], [chrLen; chrLen], [chrLen 4], 1), c.polymerizedRegions);
            
            assertTrue(m.isDnaAORIComplexAssembled);
            assertEqual([0 0], m.leadingPosition);
            assertEqual([0 0], m.helicasePosition-m.leadingPolymerasePosition);
        end
        
        function testInitiateReplication_NoWater(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            chrLen = size(m.chromosome.sequence, 1);
            helFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            uwdLen = helFtpt3 + corFtpt5 + 1;
            
            c.initialize();
            this.setUpDnaAOriComplex();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2*(uwdLen+1);
            m.substrates(m.substrateIndexs_water) = 0;
            m.enzymes(:) = 2 * m.enzymeComposition(:, m.enzymeIndexs_replisome);
            m.boundEnzymes(:) = 0;
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            
            m.initiateReplication();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            
            posStrnds = [
                m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R1234) ones(4, 1);
                m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R5) 1];
            prots = [
                repmat(m.complexIndexs_DnaA_7mer_ATP, 4, 1);
                m.complexIndexs_DnaA_1mer_ATP];
            assertEqual(CircularSparseMat(posStrnds, prots, [chrLen 4], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([1 1;1 2], [chrLen; chrLen], [chrLen 4], 1), c.polymerizedRegions);
            
            assertTrue(m.isDnaAORIComplexAssembled);
            assertEqual([0 0], m.leadingPosition);
            assertEqual([0 0], m.helicasePosition-m.leadingPolymerasePosition);
        end
        
        function testInitiateReplication_NoHelicase(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            chrLen = size(m.chromosome.sequence, 1);
            helFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            uwdLen = helFtpt3 + corFtpt5 + 1;
            
            c.initialize();
            this.setUpDnaAOriComplex();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2*(uwdLen+1);
            m.substrates(m.substrateIndexs_water) = 2*(uwdLen+1);
            m.enzymes(:) = 2 * m.enzymeComposition(:, m.enzymeIndexs_replisome);
            m.enzymes(m.enzymeIndexs_helicase) = 0;
            m.boundEnzymes(:) = 0;
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            
            m.initiateReplication();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            
            posStrnds = [
                m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R1234) ones(4, 1);
                m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R5) 1];
            prots = [
                repmat(m.complexIndexs_DnaA_7mer_ATP, 4, 1);
                m.complexIndexs_DnaA_1mer_ATP];
            assertEqual(CircularSparseMat(posStrnds, prots, [chrLen 4], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([1 1;1 2], [chrLen; chrLen], [chrLen 4], 1), c.polymerizedRegions);
            
            assertTrue(m.isDnaAORIComplexAssembled);
            assertEqual([0 0], m.leadingPosition);
            assertEqual([0 0], m.helicasePosition-m.leadingPolymerasePosition);
        end
    end
    
    %test unwindAndPolymerizeDNA
    methods
        function testUnwindAndPolymerizeDNA(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            pol2GblIdx  = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helGblIdx   = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            bcplxGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_betaClamp);
            
            chrLen = size(c.sequence, 1);
            polRate = m.dnaPolymeraseElongationRate;
            helFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            helFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_helicase);
            corFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_core);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            bClmpFtpt= m.enzymeDNAFootprints(m.enzymeIndexs_betaClamp);
            
            %% Ex 1. No lagging polymerase bound yet
            c.initialize();
            this.setUpDnaAOriComplex();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_dntp) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 8;
            m.enzymes(m.enzymeIndexs_core) = 4;
            m.enzymes(m.enzymeIndexs_gammaComplex) = 2;
            m.enzymes(m.enzymeIndexs_helicase) = 2;
            m.enzymes(m.enzymeIndexs_primase) = 2;
            m.boundEnzymes(:) = 0;
            
            m.initiateReplication();
            assertEqual([chrLen 1], m.leadingPosition);
            assertTrue(m.isDnaAORIComplexAssembled);
            assertEqual([-helFtpt holFtpt-size(c.sequence,1)], m.helicasePosition-m.leadingPolymerasePosition);
            
            m.unwindAndPolymerizeDNA();
            
            uwdLen = helFtpt3 + corFtpt5 + 1;
            polProg = m.primerLength;
            substrates = zeros(size(m.substrates));
            substrates(m.substrateIndexs_dntp) = ...
                + repmat(1e6, 4, 1) ...
                - c.sequence.subsequenceBaseCounts(1:polProg, 3) ...
                - c.sequence.subsequenceBaseCounts(chrLen-polProg+1:chrLen, 2);
            substrates(m.substrateIndexs_diphosphate) = 2*polProg;
            substrates(m.substrateIndexs_atp)         = 1e6 - 2*polProg - 2*(uwdLen+1);
            substrates(m.substrateIndexs_water)       = 1e6 - 2*polProg - 2*(uwdLen+1);
            substrates(m.substrateIndexs_adp)         = 2*polProg + 2*(uwdLen+1);
            substrates(m.substrateIndexs_hydrogen)    = 2*polProg + 2*(uwdLen+1);
            substrates(m.substrateIndexs_phosphate)   = 2*polProg + 2*(uwdLen+1);
            assertEqual(substrates', m.substrates');
            
            enzymes = zeros(size(m.enzymes));
            enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            boundEnzymes = zeros(size(m.boundEnzymes));
            boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            
            assertEqual([chrLen-polProg 1+polProg], m.leadingPosition);
            assertEqual([0 0], m.laggingPosition);
            assertEqual([-helFtpt holFtpt-size(c.sequence,1)], m.helicasePosition-m.leadingPolymerasePosition);
            
            helPos = m.helicasePosition;
            assertEqual(CircularSparseMat([
                1 1;
                helPos(2)+helFtpt3 2;
                chrLen-polProg+1 2;
                1 3;
                1 4;
                helPos(1)+helFtpt5+1 4
                ],[
                chrLen;
                helPos(1)+helFtpt5 - (helPos(2)+helFtpt3)+1;
                polProg;
                polProg;
                helPos(2)+helFtpt3-1;
                chrLen - (helPos(1)+helFtpt5);
                ], [chrLen 4], 1), c.polymerizedRegions);
            
            %% Ex 2. Lagging polymerase also bound
            c.initialize();
            
            %initiate Okazaki fragment
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4]) = pol2GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100         1; m.primaseBindingLocations{2}(1)+100         4]) = helGblIdx;
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
            
            assertEqual([0 0], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            m.initiateOkazakiFragment();
            
            substrates = zeros(size(m.substrates));
            substrates(m.substrateIndexs_adp) = 2;
            substrates(m.substrateIndexs_phosphate) = 2;
            substrates(m.substrateIndexs_hydrogen) = 2;
            assertEqual(substrates, m.substrates);
            
            enzymes = zeros(size(m.enzymes));
            boundEnzymes = zeros(size(m.enzymes));
            boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            boundEnzymes(m.enzymeIndexs_betaClamp) = 2;
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            
            posStrnds  = [
                m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4;
                m.primaseBindingLocations{1}(1)-100 1; m.primaseBindingLocations{2}(1)+100 4;
                m.primaseBindingLocations{1}(1)-corFtpt3-bClmpFtpt 3; m.primaseBindingLocations{2}(1)+corFtpt3+1 2
                ];
            prots = [
                repmat(pol2GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1);
                repmat(bcplxGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            %polymerize
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_dntp) = 1e6;
            leadingPosition = m.leadingPosition;
            laggingPosition = [m.primaseBindingLocations{1}(1) m.primaseBindingLocations{2}(1)];
            
            m.unwindAndPolymerizeDNA();
            
            leadPolProg = polRate;
            lagPolProg = m.primerLength;
            uwdLen = polRate;
            substrates = zeros(size(m.substrates));
            substrates(m.substrateIndexs_dntp) = ...
                + repmat(1e6, 4, 1) ...
                - c.sequence.subsequenceBaseCounts(leadingPosition(1) - (0:leadPolProg-1), 2) ...
                - c.sequence.subsequenceBaseCounts(leadingPosition(2) + (0:leadPolProg-1), 3) ...
                - c.sequence.subsequenceBaseCounts(laggingPosition(1) + (0:lagPolProg-1),  3) ...
                - c.sequence.subsequenceBaseCounts(laggingPosition(2) - (0:lagPolProg-1),  2);
            substrates(m.substrateIndexs_diphosphate) = 2*(leadPolProg+lagPolProg);
            substrates(m.substrateIndexs_atp)         = 1e6 - 2*uwdLen;
            substrates(m.substrateIndexs_water)       = 1e6 - 2*uwdLen;
            substrates(m.substrateIndexs_adp)         = 2*uwdLen;
            substrates(m.substrateIndexs_hydrogen)    = 2*uwdLen;
            substrates(m.substrateIndexs_phosphate)   = 2*uwdLen;
            assertEqual(substrates', m.substrates');
            
            enzymes = zeros(size(m.enzymes));
            boundEnzymes = zeros(size(m.enzymes));
            boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = 2;
            boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = 2;
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            
            assertEqual(leadingPosition + polRate*[-1 1], m.leadingPosition);
            assertEqual([m.primaseBindingLocations{1}(1)+m.primerLength m.primaseBindingLocations{2}(1)-m.primerLength], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            helPos = m.helicasePosition;
            assertEqual(CircularSparseMat([
                1 1;
                helPos(2)+helFtpt3 2;
                leadingPosition(1)-polRate+1 2;
                1 3;
                1 4;
                helPos(1)+helFtpt5+1 4;
                m.primaseBindingLocations{1}(1) 3;
                m.primaseBindingLocations{2}(1)-m.primerLength+1 2;
                ],[
                chrLen;
                helPos(1)+helFtpt5 - (helPos(2)+helFtpt3)+1;
                chrLen-(leadingPosition(1) - polRate+1)+1;
                leadingPosition(2) + polRate-1;
                helPos(2)+helFtpt3-1;
                chrLen-(helPos(1)+helFtpt5);
                m.primerLength;
                m.primerLength;
                ], [chrLen 4], 1), c.polymerizedRegions);
        end
        
        function testUnwindAndPolymerizeDNA_NoDNTPs(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            pol1GblIdx  = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex);
            pol2GblIdx  = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            pol3GblIdx  = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase);
            helGblIdx   = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            bcplxGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_betaClamp);
            
            chrLen = size(c.sequence, 1);
            helFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            corFtpt = m.enzymeDNAFootprints(m.enzymeIndexs_core);
            corFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_core);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            bClmpFtpt= m.enzymeDNAFootprints(m.enzymeIndexs_betaClamp);
            
            %% Ex 1. No lagging polymerase bound yet
            c.initialize();
            this.setUpDnaAOriComplex();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_dntp) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 8;
            m.enzymes(m.enzymeIndexs_core) = 4;
            m.enzymes(m.enzymeIndexs_gammaComplex) = 2;
            m.enzymes(m.enzymeIndexs_helicase) = 2;
            m.enzymes(m.enzymeIndexs_primase) = 2;
            m.boundEnzymes(:) = 0;
            
            m.initiateReplication();
            assertEqual([chrLen 1], m.leadingPosition);
            assertTrue(m.isDnaAORIComplexAssembled);
            assertEqual([-helFtpt holFtpt-size(c.sequence,1)], m.helicasePosition-m.leadingPolymerasePosition);
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            initial_polymerizedRegions = c.polymerizedRegions;
            initial_complexBoundSites = c.complexBoundSites;
            
            m.unwindAndPolymerizeDNA();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(initial_polymerizedRegions, c.polymerizedRegions);
            assertEqual(initial_complexBoundSites, c.complexBoundSites);
            
            %% Ex 2. Lagging polymerase also bound
            c.initialize();
            
            %initiate Okazaki fragment
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4]) = pol2GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100         1; m.primaseBindingLocations{2}(1)+100         4]) = helGblIdx;
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
            
            assertEqual([0 0], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            m.initiateOkazakiFragment();
            
            substrates = zeros(size(m.substrates));
            substrates(m.substrateIndexs_adp) = 2;
            substrates(m.substrateIndexs_phosphate) = 2;
            substrates(m.substrateIndexs_hydrogen) = 2;
            assertEqual(substrates, m.substrates);
            
            enzymes = zeros(size(m.enzymes));
            boundEnzymes = zeros(size(m.enzymes));
            boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            boundEnzymes(m.enzymeIndexs_betaClamp) = 2;
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            
            posStrnds  = [
                m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4;
                m.primaseBindingLocations{1}(1)-100 1; m.primaseBindingLocations{2}(1)+100 4;
                m.primaseBindingLocations{1}(1)-corFtpt3-bClmpFtpt 3; m.primaseBindingLocations{2}(1)+corFtpt3+1 2
                ];
            prots = [
                repmat(pol2GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1);
                repmat(bcplxGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            %polymerize
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_dntp) = 0;
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            initial_polymerizedRegions = c.polymerizedRegions;
            initial_complexBoundSites = c.complexBoundSites;
            
            m.unwindAndPolymerizeDNA();
            
            final_boundEnzymes = initial_boundEnzymes;
            final_boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 0;
            final_boundEnzymes(m.enzymeIndexs_betaClamp) = 0;
            final_boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = 2;
            final_boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = 2;
            
            final_complexBoundSites = initial_complexBoundSites;
            final_complexBoundSites([m.leadingPolymerasePosition; 1 4]') = pol1GblIdx;
            final_complexBoundSites([m.primaseBindingLocations{1}(1)-corFtpt3-bClmpFtpt 3; m.primaseBindingLocations{2}(1)+corFtpt3+1 2]) = 0;
            final_complexBoundSites([m.primaseBindingLocations{1}(1)-corFtpt3-bClmpFtpt 3; m.primaseBindingLocations{2}(1)+corFtpt3+1-corFtpt 2]) = pol3GblIdx;
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(initial_polymerizedRegions, c.polymerizedRegions);
            assertEqual(final_complexBoundSites, c.complexBoundSites);
        end
        
        function testUnwindAndPolymerizeDNA_NoWater(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            pol2GblIdx  = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helGblIdx   = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            bcplxGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_betaClamp);
            
            chrLen = size(c.sequence, 1);
            helFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            corFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_core);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            bClmpFtpt= m.enzymeDNAFootprints(m.enzymeIndexs_betaClamp);            
            
            %% Ex 1. Lagging polymerase also bound
            c.initialize();
            
            %initiate Okazaki fragment
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4]) = pol2GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100         1; m.primaseBindingLocations{2}(1)+100         4]) = helGblIdx;
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
            
            assertEqual([0 0], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            m.initiateOkazakiFragment();
            
            substrates = zeros(size(m.substrates));
            substrates(m.substrateIndexs_adp) = 2;
            substrates(m.substrateIndexs_phosphate) = 2;
            substrates(m.substrateIndexs_hydrogen) = 2;
            assertEqual(substrates, m.substrates);
            
            enzymes = zeros(size(m.enzymes));
            boundEnzymes = zeros(size(m.enzymes));
            boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            boundEnzymes(m.enzymeIndexs_betaClamp) = 2;
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            
            posStrnds  = [
                m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4;
                m.primaseBindingLocations{1}(1)-100 1; m.primaseBindingLocations{2}(1)+100 4;
                m.primaseBindingLocations{1}(1)-corFtpt3-bClmpFtpt 3; m.primaseBindingLocations{2}(1)+corFtpt3+1 2
                ];
            prots = [
                repmat(pol2GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1);
                repmat(bcplxGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            %polymerize
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            m.substrates(m.substrateIndexs_water) = 0;
            m.substrates(m.substrateIndexs_dntp) = 1e6;
            
            initial_helicasePosition = m.helicasePosition;
            initial_leadingPolPosition = m.leadingPosition;
            
            m.unwindAndPolymerizeDNA();
            
            assertEqual(initial_helicasePosition, m.helicasePosition);
            assertEqual(initial_leadingPolPosition, m.leadingPosition);
        end
        
        function testUnwindAndPolymerizeDNA_ProteinDisplacement(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            pol2GblIdx  = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helGblIdx   = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            bcplxGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_betaClamp);
            
            chrLen = size(c.sequence, 1);
            polRate = m.dnaPolymeraseElongationRate;
            helFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            helFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_helicase);
            corFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_core);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            bClmpFtpt = m.enzymeDNAFootprints(m.enzymeIndexs_betaClamp);
            
            %% Ex 1. Lagging polymerase also bound
            c.initialize();
            
            %initiate Okazaki fragment
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
                        
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4]) = pol2GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100         1; m.primaseBindingLocations{2}(1)+100         4]) = helGblIdx;
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
            
            assertEqual([0 0], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            m.initiateOkazakiFragment();
            
            posStrnds  = [
                m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4;
                m.primaseBindingLocations{1}(1)-100 1; m.primaseBindingLocations{2}(1)+100 4;
                m.primaseBindingLocations{1}(1)-corFtpt3-bClmpFtpt 3; m.primaseBindingLocations{2}(1)+corFtpt3+1 2
                ];
            prots = [
                repmat(pol2GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1);
                repmat(bcplxGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            %polymerize
            assertFalse(any(m.enzymeComplexGlobalIndexs == 1));
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100-c.complexDNAFootprints(1) 4]) = 1; %add protein in path of helicase
            c.complexBoundSites([m.primaseBindingLocations{2}(1)+100+helFtpt 1]) = 1; %add protein in path of helicase
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_dntp) = 1e6;
            
            leadingPosition = m.leadingPosition;
            
            m.unwindAndPolymerizeDNA();
            
            assertEqual(0, nnz(c.complexBoundSites == 1));
            assertEqual(6, nnz(c.complexBoundSites));
            
            assertEqual(leadingPosition + polRate*[-1 1], m.leadingPosition);
            assertEqual([m.primaseBindingLocations{1}(1)+m.primerLength m.primaseBindingLocations{2}(1)-m.primerLength], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            helPos = m.helicasePosition;
            assertEqual(CircularSparseMat([
                1 1;
                helPos(2)+helFtpt3 2;
                leadingPosition(1)-polRate+1 2;
                1 3;
                1 4;
                helPos(1)+helFtpt5+1 4;
                m.primaseBindingLocations{1}(1) 3;
                m.primaseBindingLocations{2}(1)-m.primerLength+1 2;
                ],[
                chrLen;
                helPos(1)+helFtpt5 - (helPos(2)+helFtpt3)+1;
                chrLen-(leadingPosition(1) - polRate+1)+1;
                leadingPosition(2) + polRate-1;
                helPos(2)+helFtpt3-1;
                chrLen-(helPos(1)+helFtpt5);
                m.primerLength;
                m.primerLength;
                ], [chrLen 4], 1), c.polymerizedRegions);
        end
    end
    
    %test dissociateFreeSSBComplexes
    methods
        function testDissociateFreeSSBComplexes(this)
            m = this.process;
            
            m.enzymes(:) = (1:numel(m.enzymes))';
            m.boundEnzymes(:) = numel(m.enzymes)+(1:numel(m.enzymes))';
            
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.dissociateFreeSSBComplexes();
            
            finalEnzymes = initialEnzymes;
            finalEnzymes(m.enzymeIndexs_ssb8mer) = 0;
            finalEnzymes(m.enzymeIndexs_ssb4mer) = initialEnzymes(m.enzymeIndexs_ssb4mer) + 2*initialEnzymes(m.enzymeIndexs_ssb8mer);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(initialBoundEnzymes, m.boundEnzymes);
        end
    end
    
    %test freeAndBindSSBs
    methods
        function testFreeAndBindSSBs(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            chrLen = size(c.sequence, 1);
            ftpt = m.enzymeDNAFootprints(m.enzymeIndexs_ssb8mer);
            spcg = m.ssbComplexSpacing;
            ssbGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_ssb8mer);
            
            %Ex 1: region right of oriC
            c.initialize();
            c.setRegionUnwound(1, 2000);
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ssb4mer) = 1e4;
            
            m.freeAndBindSSBs();
            
            finalEnzymes = zeros(size(m.enzymes));
            finalBoundEnzymes = zeros(size(m.enzymes));
            n = floor(2000 / (ftpt+spcg));
            finalEnzymes(m.enzymeIndexs_ssb4mer) = 1e4 - 4*n;
            finalBoundEnzymes(m.enzymeIndexs_ssb8mer) = 2*n;
            
            assertEqual(CircularSparseMat([1+(ftpt+spcg)*(0:n-1)' ones(n, 1); 1+(ftpt+spcg)*(0:n-1)' repmat(4, n, 1)], repmat(ssbGblIdx, 2*n, 1), [size(c.sequence, 1) 4], 1), c.complexBoundSites);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            
            %Ex 2: region left of oriC
            c.initialize();
            c.setRegionUnwound(size(c.sequence, 1), -2000);
            c.setRegionPolymerized([size(c.sequence, 1) 1], -2000);
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ssb4mer) = 1e4;
            
            m.freeAndBindSSBs();
            
            finalEnzymes = zeros(size(m.enzymes));
            finalBoundEnzymes = zeros(size(m.enzymes));
            n = floor(2000 / (ftpt+spcg));
            finalEnzymes(m.enzymeIndexs_ssb4mer) = 1e4 - 2*n;
            finalBoundEnzymes(m.enzymeIndexs_ssb8mer) = n;
            
            assertEqual(CircularSparseMat([size(c.sequence, 1)-ftpt+1-(ftpt+spcg)*(0:n-1)' repmat(4, n, 1)], repmat(ssbGblIdx, n, 1), [size(c.sequence, 1) 4], 1), c.complexBoundSites);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            
            %Ex 3: region about terC
            c.initialize();
            c.linkingNumbers(:, :) = 0;
            c.setRegionUnwound(1, m.terCPosition);
            c.setRegionUnwound(chrLen, -(chrLen-c.terCPosition));
            c.setRegionPolymerized([1 1], c.terCPosition-2000);
            c.setRegionPolymerized([1 2], c.terCPosition-2000);
            c.setRegionPolymerized([chrLen 1], -(chrLen-c.terCPosition-2000));
            c.setRegionPolymerized([chrLen 2], -(chrLen-c.terCPosition-2000));
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ssb4mer) = 1e4;
            
            m.freeAndBindSSBs();
            
            finalEnzymes = zeros(size(m.enzymes));
            finalBoundEnzymes = zeros(size(m.enzymes));
            n = floor(2000 / (ftpt+spcg));
            finalEnzymes(m.enzymeIndexs_ssb4mer) = 1e4 - 4*2*n;
            finalBoundEnzymes(m.enzymeIndexs_ssb8mer) = 4*n;
            
            posStrnds = [
                m.terCPosition-2000+1+(ftpt+spcg)*(0:n-1)' ones(n, 1);
                m.terCPosition+2000-ftpt+1-(ftpt+spcg)*(0:n-1)' ones(n, 1);
                m.terCPosition-2000+1+(ftpt+spcg)*(0:n-1)' repmat(4, n, 1);
                m.terCPosition+2000-ftpt+1-(ftpt+spcg)*(0:n-1)' repmat(4, n, 1)];
            assertEqual(CircularSparseMat(posStrnds, repmat(ssbGblIdx, 4*n, 1), [size(c.sequence, 1) 4], 1), c.complexBoundSites);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
        end
        
        function testFreeAndBindSSBs_NoSSBs(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            chrLen = size(c.sequence, 1);
            
            %Ex 1: region about terC, no SSB
            c.initialize();
            c.linkingNumbers(:, :) = 0;
            c.setRegionUnwound(1, m.terCPosition);
            c.setRegionUnwound(chrLen, -(chrLen-c.terCPosition));
            c.setRegionPolymerized([1 1], c.terCPosition-2000);
            c.setRegionPolymerized([1 2], c.terCPosition-2000);
            c.setRegionPolymerized([chrLen 1], -(chrLen-c.terCPosition-2000));
            c.setRegionPolymerized([chrLen 2], -(chrLen-c.terCPosition-2000));
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ssb4mer) = 0;
            
            m.freeAndBindSSBs();
            
            finalEnzymes = zeros(size(m.enzymes));
            finalBoundEnzymes = zeros(size(m.enzymes));
            finalEnzymes(m.enzymeIndexs_ssb4mer) = 0;
            
            assertEqual(CircularSparseMat([], [], [size(c.sequence, 1) 4], 1), c.complexBoundSites);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            
            %Ex 2: region about terC, limited SSB
            c.initialize();
            c.linkingNumbers(:, :) = 0;
            c.setRegionUnwound(1, m.terCPosition);
            c.setRegionUnwound(chrLen, -(chrLen-c.terCPosition));
            c.setRegionPolymerized([1 1], c.terCPosition-2000);
            c.setRegionPolymerized([1 2], c.terCPosition-2000);
            c.setRegionPolymerized([chrLen 1], -(chrLen-c.terCPosition-2000));
            c.setRegionPolymerized([chrLen 2], -(chrLen-c.terCPosition-2000));
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ssb4mer) = 10;
            
            m.freeAndBindSSBs();
            
            finalEnzymes = zeros(size(m.enzymes));
            finalBoundEnzymes = zeros(size(m.enzymes));
            finalEnzymes(m.enzymeIndexs_ssb4mer) = 0;
            finalBoundEnzymes(m.enzymeIndexs_ssb8mer) = 5;
            
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
        end
    end
    
    %test initiateOkazakiFragment
    methods
        function testInitiateOkazakiFragment(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            pol2GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            pol1GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex);
            pol3GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase);
            
            helGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            bcplxGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_betaClamp);
            chrLen = size(c.sequence, 1);
            helFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            helFtpt3  = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            corFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_core);
            bClmpFtpt= m.enzymeDNAFootprints(m.enzymeIndexs_betaClamp);
            
            %Ex 1. first Okazaki fragment
            c.initialize();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4]) = pol2GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100 1; m.primaseBindingLocations{2}(1)+100 4]) = helGblIdx;
            
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
            
            assertEqual([0 0], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            m.initiateOkazakiFragment();
            
            substrates = zeros(size(m.substrates));
            substrates(m.substrateIndexs_adp) = 2;
            substrates(m.substrateIndexs_phosphate) = 2;
            substrates(m.substrateIndexs_hydrogen) = 2;
            assertEqual(substrates, m.substrates);
            
            enzymes = zeros(size(m.enzymes));
            boundEnzymes = zeros(size(m.enzymes));
            boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            boundEnzymes(m.enzymeIndexs_betaClamp) = 2;
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            
            posStrnds  = [
                m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4;
                m.primaseBindingLocations{1}(1)-100 1; m.primaseBindingLocations{2}(1)+100 4;
                m.primaseBindingLocations{1}(1)-corFtpt3-bClmpFtpt 3; m.primaseBindingLocations{2}(1)+corFtpt3+1 2;
                ];
            prots = [
                repmat(pol2GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1);
                repmat(bcplxGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            m.initiateOkazakiFragment();
            assertEqual(substrates, m.substrates);
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            %Ex 2. other Okazaki fragment
            c.initialize();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = 2;
            m.boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.primaseBindingLocations{1}(3)-100+helFtpt 1; m.primaseBindingLocations{2}(3)+100-holFtpt 4]) = pol1GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(3)-100 1; m.primaseBindingLocations{2}(3)+100 4]) = helGblIdx;
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.complexBoundSites([m.primaseBindingLocations{1}(1)+1000-(holFtpt-corFtpt5)+1 3; m.primaseBindingLocations{2}(1)-1000-corFtpt5 2]) = pol3GblIdx;
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
            c.setRegionPolymerized([m.laggingPosition(1)+1 2], chrLen - (m.laggingPosition(1)+1) +1);
            c.setRegionPolymerized([1 1], m.laggingPosition(2) - 1);
            
            assertEqual([m.primaseBindingLocations{1}(1)+1000 m.primaseBindingLocations{2}(1)-1000], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            m.initiateOkazakiFragment();
            
            substrates = zeros(size(m.substrates));
            substrates(m.substrateIndexs_adp) = 2;
            substrates(m.substrateIndexs_phosphate) = 2;
            substrates(m.substrateIndexs_hydrogen) = 2;
            assertEqual(substrates, m.substrates);
            
            enzymes = zeros(size(m.enzymes));
            boundEnzymes = zeros(size(m.enzymes));
            boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = 2;
            boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = 2;
            boundEnzymes(m.enzymeIndexs_betaClamp) = 2;
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            
            posStrnds  = [
                m.primaseBindingLocations{1}(3)-100+helFtpt 1; m.primaseBindingLocations{2}(3)+100-holFtpt 4;
                m.primaseBindingLocations{1}(1)+1000-(holFtpt-corFtpt5)+1 3; m.primaseBindingLocations{2}(1)-1000-corFtpt5 2;
                m.primaseBindingLocations{1}(3)-100 1; m.primaseBindingLocations{2}(3)+100 4;
                m.primaseBindingLocations{1}(2)-corFtpt3-bClmpFtpt 3; m.primaseBindingLocations{2}(2)+corFtpt3+1 2;
                ];
            prots = [
                repmat(pol1GblIdx, 2, 1);
                repmat(pol3GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1);
                repmat(bcplxGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            m.initiateOkazakiFragment();
            assertEqual(substrates, m.substrates);
            assertEqual(enzymes, m.enzymes);
            assertEqual(boundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
        end
        
        function testInitiateOkazakiFragment_NoATP(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            pol2GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            pol1GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex);
            pol3GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase);
            
            helGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            chrLen = size(c.sequence, 1);
            helFtpt   = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            helFtpt3  = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            holFtpt   = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            corFtpt5  = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            
            %Ex 1. first Okazaki fragment
            c.initialize();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 0;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4]) = pol2GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100 1; m.primaseBindingLocations{2}(1)+100 4]) = helGblIdx;
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
            
            assertEqual([0 0], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            
            m.initiateOkazakiFragment();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            
            posStrnds  = [
                m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4;
                m.primaseBindingLocations{1}(1)-100 1; m.primaseBindingLocations{2}(1)+100 4;
                ];
            prots = [
                repmat(pol2GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            m.initiateOkazakiFragment();
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            %Ex 2. other Okazaki fragment
            c.initialize();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 0;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = 2;
            m.boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.primaseBindingLocations{1}(3)-100+helFtpt 1; m.primaseBindingLocations{2}(3)+100-holFtpt 4]) = pol1GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(3)-100 1; m.primaseBindingLocations{2}(3)+100 4]) = helGblIdx;
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.complexBoundSites([m.primaseBindingLocations{1}(1)+1000-(holFtpt-corFtpt5)+1 3; m.primaseBindingLocations{2}(1)-1000-corFtpt5 2]) = pol3GblIdx;
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
            c.setRegionPolymerized([m.laggingPosition(1)+1 2], chrLen - (m.laggingPosition(1)+1) +1);
            c.setRegionPolymerized([1 1], m.laggingPosition(2) - 1);
            
            assertEqual([m.primaseBindingLocations{1}(1)+1000 m.primaseBindingLocations{2}(1)-1000], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            
            m.initiateOkazakiFragment();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            
            posStrnds  = [
                m.primaseBindingLocations{1}(3)-100+helFtpt 1; m.primaseBindingLocations{2}(3)+100-holFtpt 4;
                m.primaseBindingLocations{1}(1)+1000-(holFtpt-corFtpt5)+1 3; m.primaseBindingLocations{2}(1)-1000-corFtpt5 2;
                m.primaseBindingLocations{1}(3)-100 1; m.primaseBindingLocations{2}(3)+100 4;
                ];
            prots = [
                repmat(pol1GblIdx, 2, 1);
                repmat(pol3GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            m.initiateOkazakiFragment();
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
        end
        
        function testInitiateOkazakiFragment_NoWater(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            pol2GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            
            helGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            chrLen = size(c.sequence, 1);
            helFtpt   = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            helFtpt3  = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            holFtpt   = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            corFtpt5  = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            
            %Ex 1. first Okazaki fragment
            c.initialize();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_water) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4]) = pol2GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100 1; m.primaseBindingLocations{2}(1)+100 4]) = helGblIdx;
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
            
            assertEqual([0 0], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            
            m.initiateOkazakiFragment();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            
            posStrnds  = [
                m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4;
                m.primaseBindingLocations{1}(1)-100 1; m.primaseBindingLocations{2}(1)+100 4;
                ];
            prots = [
                repmat(pol2GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            m.initiateOkazakiFragment();
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);            
        end
        
        function testInitiateOkazakiFragment_NoBetaClampMonomer(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            pol2GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            pol1GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex);
            pol3GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase);
            
            helGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            chrLen = size(c.sequence, 1);
            helFtpt   = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            helFtpt3  = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            holFtpt   = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            corFtpt5  = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            
            %Ex 1. first Okazaki fragment
            c.initialize();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 0;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4]) = pol2GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100 1; m.primaseBindingLocations{2}(1)+100 4]) = helGblIdx;
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
            
            assertEqual([0 0], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            
            m.initiateOkazakiFragment();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            
            posStrnds  = [
                m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4;
                m.primaseBindingLocations{1}(1)-100 1; m.primaseBindingLocations{2}(1)+100 4;
                ];
            prots = [
                repmat(pol2GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            m.initiateOkazakiFragment();
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            %Ex 2. other Okazaki fragment
            c.initialize();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 0;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = 2;
            m.boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.primaseBindingLocations{1}(3)-100+helFtpt 1; m.primaseBindingLocations{2}(3)+100-holFtpt 4]) = pol1GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(3)-100 1; m.primaseBindingLocations{2}(3)+100 4]) = helGblIdx;
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.complexBoundSites([m.primaseBindingLocations{1}(1)+1000-(holFtpt-corFtpt5)+1 3; m.primaseBindingLocations{2}(1)-1000-corFtpt5 2]) = pol3GblIdx;
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
            c.setRegionPolymerized([m.laggingPosition(1)+1 2], chrLen - (m.laggingPosition(1)+1) +1);
            c.setRegionPolymerized([1 1], m.laggingPosition(2) - 1);
            
            assertEqual([m.primaseBindingLocations{1}(1)+1000 m.primaseBindingLocations{2}(1)-1000], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            
            m.initiateOkazakiFragment();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            
            posStrnds  = [
                m.primaseBindingLocations{1}(3)-100+helFtpt 1; m.primaseBindingLocations{2}(3)+100-holFtpt 4;
                m.primaseBindingLocations{1}(1)+1000-(holFtpt-corFtpt5)+1 3; m.primaseBindingLocations{2}(1)-1000-corFtpt5 2;
                m.primaseBindingLocations{1}(3)-100 1; m.primaseBindingLocations{2}(3)+100 4;
                ];
            prots = [
                repmat(pol1GblIdx, 2, 1);
                repmat(pol3GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            m.initiateOkazakiFragment();
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
        end
    end
    
    %test terminateOkazakiFragment
    methods
        function testTerminateOkazakiFragment(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            ssbGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_ssb8mer);
            pol2GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            pol1GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex);
            pol3GblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase);
            
            helGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            bcplxGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_betaClamp);
            chrLen = size(c.sequence, 1);
            helFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            helFtpt3  = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            corFtpt = m.enzymeDNAFootprints(m.enzymeIndexs_core);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            
            %Ex 1. first Okazaki fragment
            c.initialize();
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ssb4mer) = 200;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = 2;
            m.boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_betaClamp) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.primaseBindingLocations{1}(4)-100+helFtpt 1; m.primaseBindingLocations{2}(4)+100-holFtpt 4]) = pol1GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(4)-100 1; m.primaseBindingLocations{2}(4)+100 4]) = helGblIdx;
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
            c.setRegionPolymerized([m.primaseBindingLocations{1}(1) 2], chrLen - m.primaseBindingLocations{1}(1) +1);
            c.setRegionPolymerized([1 1], m.primaseBindingLocations{2}(1));
            c.complexBoundSites([chrLen+1-(holFtpt-corFtpt5)+1 3; 0-corFtpt5 2]) = pol3GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(2)-(holFtpt-corFtpt5)+1 3; m.primaseBindingLocations{2}(2)-corFtpt5+corFtpt 2]) = bcplxGblIdx;
            
            assertEqual([1 chrLen], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            assertEqual([1 1], m.okazakiFragmentIndex);
            assertEqual([chrLen-m.primaseBindingLocations{1}(1)+1 m.primaseBindingLocations{2}(1)], m.okazakiFragmentLength);
            assertEqual(m.okazakiFragmentLength, m.okazakiFragmentProgress);
            
            m.freeAndBindSSBs();
            assertEqual([true true], m.areLaggingStrandSSBSitesBound);
            
            m.terminateOkazakiFragment();
            m.dissociateFreeSSBComplexes();
            
            assertEqual([m.primaseBindingLocations{1}(2) m.primaseBindingLocations{2}(2)], m.laggingPosition);
            assertEqual([2 2], m.okazakiFragmentIndex);
            assertEqual([0 0], m.okazakiFragmentProgress);
            
            substrates = zeros(size(m.substrates));
            assertEqual(substrates, m.substrates);
            
            ssbPosStrnds = find(c.complexBoundSites == ssbGblIdx);
            enzymes = zeros(size(m.enzymes));
            enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            enzymes(m.enzymeIndexs_ssb4mer) = 200-2*size(ssbPosStrnds, 1);
            boundEnzymes = zeros(size(m.enzymes));
            boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = 2;
            boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = 2;
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            boundEnzymes(m.enzymeIndexs_ssb8mer) = size(ssbPosStrnds, 1);
            assertEqual(enzymes', m.enzymes');
            assertEqual(boundEnzymes', m.boundEnzymes');
            
            posStrnds  = [
                ssbPosStrnds;
                m.primaseBindingLocations{1}(4)-100+helFtpt 1; m.primaseBindingLocations{2}(4)+100-holFtpt 4;
                m.primaseBindingLocations{1}(2)-(holFtpt-corFtpt5)+1 3; m.primaseBindingLocations{2}(2)-corFtpt5 2;
                m.primaseBindingLocations{1}(4)-100 1; m.primaseBindingLocations{2}(4)+100 4;
                ];
            prots = [
                repmat(ssbGblIdx, size(ssbPosStrnds, 1), 1);
                repmat(pol1GblIdx, 2, 1);
                repmat(pol3GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            assertEqual(CircularSparseMat([chrLen 2; chrLen 3], [1; 1], [size(c.sequence, 1), 4], 1), c.strandBreaks);
            
            %Ex 2. last Okazaki fragment
            c.initialize();
            c.polymerizedRegions(1, :) = chrLen;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_coreBetaClampGammaComplex) = 2;
            m.boundEnzymes(m.enzymeIndexs_coreBetaClampPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.terCPosition-corFtpt5 1; m.terCPosition+1-holFtpt+corFtpt5+1 4]) = pol1GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(end-1)-(holFtpt-corFtpt5)+1 3; m.primaseBindingLocations{2}(end-1)-corFtpt5 2]) = pol3GblIdx;
            c.complexBoundSites([m.terCPosition-corFtpt5-helFtpt 1; m.terCPosition+1+corFtpt5+1 4]) = helGblIdx;
            
            assertEqual([numel(m.primaseBindingLocations{1}) numel(m.primaseBindingLocations{2})], m.okazakiFragmentIndex);
            assertEqual(m.okazakiFragmentLength, m.okazakiFragmentProgress);
            assertEqual([m.terCPosition m.terCPosition+1], m.leadingPosition);
            assertEqual([m.primaseBindingLocations{1}(end-1) m.primaseBindingLocations{2}(end-1)], m.laggingPosition);
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            
            m.freeAndBindSSBs();
            assertEqual([true true], m.areLaggingStrandSSBSitesBound);
            
            m.terminateOkazakiFragment();
            m.dissociateFreeSSBComplexes();
            
            assertEqual([0 0], m.laggingPosition);
            assertEqual([0 0], m.okazakiFragmentIndex);
            assertEqual([0 0], m.okazakiFragmentProgress);
            
            substrates = zeros(size(m.substrates));
            assertEqual(substrates, m.substrates);
            
            enzymes = zeros(size(m.enzymes));
            enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            boundEnzymes = zeros(size(m.enzymes));
            boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            boundEnzymes(m.enzymeIndexs_helicase) = 2;
            assertEqual(enzymes', m.enzymes');
            assertEqual(boundEnzymes', m.boundEnzymes');
            
            posStrnds  = [
                m.terCPosition-corFtpt5 1; m.terCPosition+1-holFtpt+corFtpt5+1 4;
                m.terCPosition-corFtpt5-helFtpt 1; m.terCPosition+1+corFtpt5+1 4;
                ];
            prots = [
                repmat(pol2GblIdx, 2, 1);
                repmat(helGblIdx, 2, 1)];
            assertEqual(CircularSparseMat(posStrnds, prots, [size(c.sequence, 1), 4], 1), c.complexBoundSites);
            
            assertEqual(CircularSparseMat([m.primaseBindingLocations{1}(end-1)-1 3; m.primaseBindingLocations{2}(end-1) 2], [1; 1], [size(c.sequence, 1), 4], 1), c.strandBreaks);
        end
    end
    
    %test ligateDNA
    methods
        function testLigateDNA(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            %example 1: no damages
            c.initialize();
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.ligateDNA();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 2: single strand break
            c.initialize();
            c.strandBreaks(100,1)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_ligase)=1e6;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalSubstrates(m.substrateIndexs_nad) = finalSubstrates(m.substrateIndexs_nad)-1;
            finalSubstrates(m.substrateIndexs_nmn) = finalSubstrates(m.substrateIndexs_nmn)+1;
            finalSubstrates(m.substrateIndexs_amp) = finalSubstrates(m.substrateIndexs_amp)+1;
            finalSubstrates(m.substrateIndexs_hydrogen) = finalSubstrates(m.substrateIndexs_hydrogen)+1;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.ligateDNA();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 3: single strand break + nearby gap
            c.initialize();
            c.strandBreaks(100,1)=1;
            c.gapSites(100,1)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_ligase)=1e6;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.ligateDNA();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
        end
        
        function testLigateDNA_InsufficientMetabolites(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
            c.strandBreaks(100,1)=1;
            m.substrates(:)=1e6;
            m.substrates(m.substrateIndexs_nad)=0;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_ligase)=1e6;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.ligateDNA();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
        end
        
        function testLigateDNA_InsufficientEnzymes(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
            c.strandBreaks(100,1)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_ligase)=0;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.ligateDNA();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
        end
    end
    
    %test terminateReplication
    methods
        function testTerminateReplication(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            helFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            
            %Ex 1
            c.initialize();
            c.polymerizedRegions(1, :) = c.polymerizedRegions(1, 1);
            c.linkingNumbers(1, :) = c.linkingNumbers(1, 1);
            c.complexBoundSites(m.terCPosition+1-corFtpt5-1, 1) = polGblIdx;
            c.complexBoundSites(m.terCPosition+1-(holFtpt-corFtpt5)+1, 4) = polGblIdx;
            c.complexBoundSites(m.terCPosition+1-helFtpt-corFtpt5-1, 1) = helGblIdx;
            c.complexBoundSites(m.terCPosition+1+corFtpt5+1, 4) = helGblIdx;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            m.boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            
            assertEqual([-helFtpt holFtpt], m.helicasePosition-m.leadingPolymerasePosition);
            assertEqual([m.terCPosition, m.terCPosition+1], m.leadingPosition);
            assertEqual([0 0], m.laggingPosition);
            
            m.terminateReplication();
            
            assertEqual(0, unique(m.substrates));
            enzymes = 2 * m.enzymeComposition(:, m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            enzymes(m.enzymeIndexs_helicase) = 2;
            assertEqual(enzymes, m.enzymes);
            assertEqual(0, unique(m.boundEnzymes));
            assertEqual(0, nnz(c.complexBoundSites));
            
            assertEqual(CircularSparseMat([m.terCPosition 2; m.terCPosition 3], [1;1], [size(c.sequence, 1) 4], 1), c.strandBreaks);
        end
        
        function testRNAPolymeraseCollision(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            comp = m.compartment;
            c = m.chromosome;
            pc = m.complex;
            r = c.rnaPolymerase;
            t = r.transcripts;
            
            pol2GblIdx  = m.enzymeGlobalIndexs(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helGblIdx   = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            
            chrLen = size(c.sequence, 1);
            polRate = m.dnaPolymeraseElongationRate;
            helFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_helicase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            
            c.initialize();
            
            %initiate Okazaki fragment
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_betaClampMonomer) = 4;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = 2;
            m.boundEnzymes(m.enzymeIndexs_helicase) = 2;
            
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100+helFtpt 1; m.primaseBindingLocations{2}(1)+100-holFtpt 4]) = pol2GblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-100         1; m.primaseBindingLocations{2}(1)+100         4]) = helGblIdx;
            c.setRegionUnwound(chrLen, -(chrLen - (m.leadingPosition(1)+1) +1 + corFtpt5+helFtpt3+1));
            c.setRegionUnwound(1, m.leadingPosition(2) - 1 + corFtpt5+helFtpt3+1);
            c.setRegionPolymerized([m.leadingPosition(1)+1 1], chrLen - (m.leadingPosition(1)+1) +1);
            c.setRegionPolymerized([1 2], m.leadingPosition(2) - 1);
                       
            m.initiateOkazakiFragment();
            
            %polymerize
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_dntp) = 1e6;
            
            m.unwindAndPolymerizeDNA();
            
            helicasePosition = m.leadingPosition;
            laggingPosition = m.laggingPosition;
            
            r.states = repmat(r.activelyTranscribingValue + 1, 4, 1);
            r.positionStrands = [
                helicasePosition(1)-polRate 1
                helicasePosition(2)+polRate 2
                laggingPosition(1)+polRate 2
                laggingPosition(2)-polRate 1
                ];
            t.boundTranscriptionUnits = ones(4, 1);
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptChromosome = ones(4, 1);
            t.abortedTranscripts = zeros(0, 2);
            t.transcriptionUnitFivePrimeCoordinates(1:4) = r.positionStrands(:, 1)-1;
            t.transcriptionUnitDirections(1:4) = [1; 0; 0; 1];
            [~, rnaPolFtpt3, rnaPolFtpt5] = c.getDNAFootprint([], pc.rnaPolymeraseIndexs(1));
            c.complexBoundSites(r.positionStrands(1, :) - [rnaPolFtpt5 0]) = pc.rnaPolymeraseIndexs(1);
            c.complexBoundSites(r.positionStrands(2, :) - [rnaPolFtpt3 0]) = pc.rnaPolymeraseIndexs(1);
            c.complexBoundSites(r.positionStrands(3, :) - [rnaPolFtpt3 0]) = pc.rnaPolymeraseIndexs(1);
            c.complexBoundSites(r.positionStrands(4, :) - [rnaPolFtpt5 0]) = pc.rnaPolymeraseIndexs(1);
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs) = 4;           
            
            m.rnaPolymeraseCollisionMeanDwellTime = Inf;
            m.unwindAndPolymerizeDNA();
            assertEqual(4, nnz(c.complexBoundSites == pc.rnaPolymeraseIndexs(1)));
            assertEqual(0, numel(t.abortedSequences));
            
            m.rnaPolymeraseCollisionMeanDwellTime = 0;
            m.unwindAndPolymerizeDNA();
            assertEqual(0, nnz(c.complexBoundSites == pc.rnaPolymeraseIndexs(1)));
            assertEqual(4, numel(t.abortedSequences));
        end
        
        function testReplicationTerminationInteractionWithSupercoils(this)
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            
            m = this.process;
            c = m.chromosome;
            
            m.okazakiFragmentMeanLength = 500;
            m.laggingBackupClampReloadingLength = 250;
            
            bases = 'ACGT';
            seq = bases(randi(4, 2*6*m.okazakiFragmentMeanLength, 1));
            c.sequence = ChromosomeSequence(seq);
            c.sequenceLen = length(seq);
            c.sequenceGCContent = getGCContent(c.sequence);
            c.terCPosition = c.sequenceLen / 2;
            m.terCPosition = c.sequenceLen / 2;
            m.primaseBindingLocations = {
                size(c.sequence, 1) + 1 - m.calculatePrimaseBindingLocations(size(c.sequence, 1) - m.terCPosition)
                m.calculatePrimaseBindingLocations(m.terCPosition)
                }';
            m.dnaAFunctionalBoxStartPositions = c.sequenceLen - (40 + 10*(1:5)');
            
            %% initial state
            c.initialize();
            LK_ss = collapse(c.linkingNumbers)/2;
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_dntp) = [
                sum(sum(c.sequence(:, 2:3) == 'A'));
                sum(sum(c.sequence(:, 2:3) == 'C'));
                sum(sum(c.sequence(:, 2:3) == 'G'));
                sum(sum(c.sequence(:, 2:3) == 'T'))];
            m.substrates(m.substrateIndexs_atp) = ...
                + size(c.sequence, 1) ...                              %helicase unwinding each base pair
                + sum(cellfun(@numel, m.primaseBindingLocations)) ...  %lagging strand beta-clamp formation
                + 2;                                                   %leading strand beta-clamp formation
            m.substrates(m.substrateIndexs_water) = m.substrates(m.substrateIndexs_atp);
            m.substrates(m.substrateIndexs_nad) = sum(cellfun(@numel,m.primaseBindingLocations)) + 2;
            
            m.enzymes = m.enzymeComposition(:, [m.enzymeIndexs_replisome; m.enzymeIndexs_betaClamp; m.enzymeIndexs_ssb8mer]) * [2; 2; 100];
            m.enzymes(m.enzymeIndexs_ligase) = 100;
            m.boundEnzymes(:) = 0;
            
            this.setUpDnaAOriComplex();
            
            %% evolve state -- sufficient time to replicate, see testCompleteReplication_SmallChromosome
            for i = 1:60
                m.evolveState();
            end
            
            %% assertions 
            
            %replication finishes, and isn't stuck behind unresolved
            %supercoils
            assertFalse(any(m.leadingPosition));
            assertFalse(any(m.laggingPosition));
            
            %supercoils added behind replication forks, and push toward
            %terC ahead of replication forks
            helFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_helicase);
            unPolLen = c.polymerizedRegions([min(c.terCPosition+1, m.helicasePosition(2)+ helFtpt5-1) 2]);
            assertElementsAlmostEqual(LK_ss + (collapse(c.doubleStrandedRegions) / 2 - unPolLen) / c.relaxedBasesPerTurn, ...
                collapse(c.linkingNumbers)/2, ...
                'relative', 1e-6);
        end
        
        function testInitiationWithSlowATPProduction(this)
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            
            m = this.process;
            c = m.chromosome;
            
            m.okazakiFragmentMeanLength = 500;
            m.laggingBackupClampReloadingLength = 250;
            
            bases = 'ACGT';
            seq = bases(randi(4, 2*6*m.okazakiFragmentMeanLength, 1));
            c.sequence = ChromosomeSequence(seq);
            c.sequenceLen = length(seq);
            c.sequenceGCContent = getGCContent(c.sequence);
            c.terCPosition = c.sequenceLen / 2;
            m.terCPosition = c.sequenceLen / 2;
            m.primaseBindingLocations = {
                size(c.sequence, 1) + 1 - m.calculatePrimaseBindingLocations(size(c.sequence, 1) - m.terCPosition)
                m.calculatePrimaseBindingLocations(m.terCPosition)
                }';
            m.dnaAFunctionalBoxStartPositions = c.sequenceLen - (40 + 10*(1:5)');
            
            %% initial state
            c.initialize();
            LK_ss = collapse(c.linkingNumbers)/2;
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_dntp) = [
                sum(sum(c.sequence(:, 2:3) == 'A'));
                sum(sum(c.sequence(:, 2:3) == 'C'));
                sum(sum(c.sequence(:, 2:3) == 'G'));
                sum(sum(c.sequence(:, 2:3) == 'T'))];
            m.substrates(m.substrateIndexs_nad) = sum(cellfun(@numel,m.primaseBindingLocations)) + 2;
            
            m.enzymes = m.enzymeComposition(:, [m.enzymeIndexs_replisome; m.enzymeIndexs_betaClamp; m.enzymeIndexs_ssb8mer]) * [2; 2; 100];
            m.enzymes(m.enzymeIndexs_ligase) = 100;
            m.boundEnzymes(:) = 0;
            
            this.setUpDnaAOriComplex();
            
            %% evolve state -- sufficient time to replicate, see testCompleteReplication_SmallChromosome
            atpProd = 10;
            helFtpt3 = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_helicase);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            uwdLen = helFtpt3 + corFtpt5 + 1;
            assertTrue(atpProd < uwdLen);
            for i = 1:200
                if ~any(m.leadingPosition)
                    %not so much ATP that initiation can complete in one
                    %step
                    m.substrates(m.substrateIndexs_atp) = atpProd;
                    m.substrates(m.substrateIndexs_water) = atpProd;
                else
                    %ample ATP
                    m.substrates(m.substrateIndexs_atp) = 1000;
                    m.substrates(m.substrateIndexs_water) = 1000;
                end
                
                m.evolveState();
                
                if collapse(c.polymerizedRegions) == 4 * c.sequenceLen && ~any(m.helicasePosition)
                    break;
                end
            end
            
            %% assertions
            
            %replication finishes, and isn't stuck behind unresolved
            %supercoils
            assertFalse(any(m.helicasePosition));
            assertFalse(any(m.leadingPosition));
            assertFalse(any(m.laggingPosition));
            
            %supercoils added behind replication forks, and push toward
            %terC ahead of replication forks
            helFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_helicase);
            unPolLen = c.polymerizedRegions([min(c.terCPosition+1, m.helicasePosition(2)+ helFtpt5-1) 2]);
            assertElementsAlmostEqual(LK_ss + (collapse(c.doubleStrandedRegions) / 2 - unPolLen) / c.relaxedBasesPerTurn, ...
                collapse(c.linkingNumbers)/2, ...
                'relative', 1e-6);
        end
    end
    
    %test helper functions
    methods
        function testCalculatePrimaseBindingLocations(this)
            m = this.process;
            c = m.chromosome;
            seq = c.sequence;
            
            primaseBindingLocations = {
                size(seq, 1) + 1 - m.calculatePrimaseBindingLocations(size(seq, 1) - m.terCPosition)
                m.calculatePrimaseBindingLocations(m.terCPosition)
                }';
            
            assertEqual(primaseBindingLocations{1}(end), m.terCPosition+1);
            assertEqual(primaseBindingLocations{2}(end), m.terCPosition);
            
            %Check that Okazaki fragment lengths are poission distributed with
            %lambda = okazakiFragmentMeanLength (that is the distribution has
            %mean and variance okazakiFragmentMeanLength)
            assertTrue(numel(primaseBindingLocations{1}) < 1.1 * m.terCPosition / m.okazakiFragmentMeanLength);
            assertTrue(numel(primaseBindingLocations{1}) > 0.9 * m.terCPosition / m.okazakiFragmentMeanLength);
            assertTrue(numel(primaseBindingLocations{2}) < 1.1 * m.terCPosition / m.okazakiFragmentMeanLength);
            assertTrue(numel(primaseBindingLocations{2}) > 0.9 * m.terCPosition / m.okazakiFragmentMeanLength);
            
            assertTrue(mean(abs(diff(primaseBindingLocations{1}))) > 0.9 * m.okazakiFragmentMeanLength);
            assertTrue(mean(abs(diff(primaseBindingLocations{1}))) < 1.1 * m.okazakiFragmentMeanLength);
            assertTrue(mean(abs(diff(primaseBindingLocations{2}))) > 0.9 * m.okazakiFragmentMeanLength);
            assertTrue(mean(abs(diff(primaseBindingLocations{2}))) < 1.1 * m.okazakiFragmentMeanLength);
            
            assertTrue(std(abs(diff(primaseBindingLocations{1}))) > 0.9 * sqrt(m.okazakiFragmentMeanLength));
            assertTrue(std(abs(diff(primaseBindingLocations{1}))) < 1.1 * sqrt(m.okazakiFragmentMeanLength));
            assertTrue(std(abs(diff(primaseBindingLocations{2}))) > 0.9 * sqrt(m.okazakiFragmentMeanLength));
            assertTrue(std(abs(diff(primaseBindingLocations{2}))) < 1.1 * sqrt(m.okazakiFragmentMeanLength));
            
            assertTrue(min(abs(diff(primaseBindingLocations{1}))) > m.laggingBackupClampReloadingLength);
            assertTrue(min(abs(diff(primaseBindingLocations{2}))) > m.laggingBackupClampReloadingLength);
        end
    end
    
    %test getters
    methods
        function testLeadingStrandElongating(this)
            m = this.process;
            c = m.chromosome;
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex);
            helGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            
            %Ex 1: nothing bound
            c.complexBoundSites(:, :) = 0;
            assertEqual([false false], m.leadingStrandElongating);
            
            %Ex 2
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([100 2; m.terCPosition+200 1; 300 3]) = polGblIdx;
            c.complexBoundSites([m.terCPosition+200-24 1]) = helGblIdx;
            assertEqual([true false], m.leadingStrandElongating);
            
            %Ex 3
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+200 1; 300 4]) = polGblIdx;
            c.complexBoundSites([300-24 4]) = helGblIdx;
            assertEqual([false true], m.leadingStrandElongating);
            
            %Ex 4
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+200 1; 300 4; 1305 4]) = polGblIdx;
            assertExceptionThrown(@() m.leadingStrandElongating, 'Replication:error');
        end
        
        function testLaggingStrandElongating(this)
            m = this.process;
            c = m.chromosome;
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase);
            
            %Ex 1: nothing bound
            c.complexBoundSites(:, :) = 0;
            assertEqual([false false], m.laggingStrandElongating);
            
            %Ex 2
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+300 3; m.terCPosition+1305 3; 400 1]) = polGblIdx;
            assertEqual([true false], m.laggingStrandElongating);
            
            %Ex 3
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([200 2; 400 2;m.terCPosition+300 3; m.terCPosition+1305 3; 400 1]) = polGblIdx;
            assertEqual([true true], m.laggingStrandElongating);
            
            %Ex 4
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([200 2; 400 2; 500 2;m.terCPosition+300 3; m.terCPosition+1305 3; 400 1]) = polGblIdx;
            assertExceptionThrown(@() m.laggingStrandElongating, 'Replication:error');
        end
        
        function testLeadingStrandPolymerized(this)
            m = this.process;
            c = m.chromosome;
            
            %Ex 1
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, 1:2) = size(c.sequence, 1);
            assertEqual([false false], m.leadingStrandPolymerized);
            
            %Ex 2
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, 1:2) = size(c.sequence, 1);
            c.setRegionUnwound(1, size(c.sequence, 1));
            c.setRegionPolymerized([1 2], m.terCPosition-1);
            assertEqual([false false], m.leadingStrandPolymerized);
            
            %Ex 3
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, 1:2) = size(c.sequence, 1);
            c.setRegionUnwound(1, size(c.sequence, 1));
            c.setRegionPolymerized([1 2], m.terCPosition);
            c.setRegionPolymerized([m.terCPosition+2 1], m.terCPosition-1);
            assertEqual([false true], m.leadingStrandPolymerized);
            
            %Ex 3
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, 1:2) = size(c.sequence, 1);
            c.setRegionUnwound(1, size(c.sequence, 1));
            c.setRegionPolymerized([2 2], m.terCPosition-1);
            c.setRegionPolymerized([m.terCPosition+1 1], m.terCPosition);
            assertEqual([true false], m.leadingStrandPolymerized);
        end
        
        function testLaggingStrandPolymerized(this)
            m = this.process;
            c = m.chromosome;
            
            %Ex 1
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, 1:2) = size(c.sequence, 1);
            assertEqual([false false], m.laggingStrandPolymerized);
            
            %Ex 2
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, 1:2) = size(c.sequence, 1);
            c.setRegionUnwound(1, size(c.sequence, 1));
            c.setRegionPolymerized([1 1], m.terCPosition-1);
            assertEqual([false false], m.laggingStrandPolymerized);
            
            %Ex 3
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, 1:2) = size(c.sequence, 1);
            c.setRegionUnwound(1, size(c.sequence, 1));
            c.setRegionPolymerized([1 1], m.terCPosition);
            c.setRegionPolymerized([m.terCPosition+2 2], m.terCPosition-1);
            assertEqual([false true], m.laggingStrandPolymerized);
            
            %Ex 3
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, 1:2) = size(c.sequence, 1);
            c.setRegionUnwound(1, size(c.sequence, 1));
            c.setRegionPolymerized([2 1], m.terCPosition-1);
            c.setRegionPolymerized([m.terCPosition+1 2], m.terCPosition);
            assertEqual([true false], m.laggingStrandPolymerized);
        end
        
        function testNumLigations(this)
            m = this.process;
            c = m.chromosome;
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_coreBetaClampGammaComplex);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            
            %Ex 1
            c.polymerizedRegions(:, :) = 0;
            c.strandBreaks(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            assertEqual([0 0], m.numLigations);
            assertEqual([false false], m.strandLigated);
            
            %Ex 2
            c.polymerizedRegions(:, :) = 0;
            c.strandBreaks(:, :) = 0;
            c.strandBreaks(m.primaseBindingLocations{1}(7)-1, 3) = 1;
            c.strandBreaks(m.primaseBindingLocations{2}(8),   2) = 1;
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(3)-holFtpt+corFtpt5+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(4)-corFtpt5           2]) = polGblIdx;
            assertEqual([2 3], m.numLigations);
            assertEqual([false false], m.strandLigated);
            
            %Ex 3
            c.polymerizedRegions(:, :) = 0;
            c.strandBreaks(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-holFtpt+corFtpt5+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(1)-corFtpt5 2]) = polGblIdx;
            assertEqual([0 0], m.numLigations);
            assertEqual([false false], m.strandLigated);
            
            %Ex 4
            c.polymerizedRegions(:, :) = 0;
            c.strandBreaks(:, :) = 0;
            c.strandBreaks(m.primaseBindingLocations{1}(7)-1, 3) = 1;
            c.strandBreaks(m.primaseBindingLocations{2}(8), 2) = 1;
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(end)-holFtpt+corFtpt5+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(end)-corFtpt5 2]) = polGblIdx;
            assertEqual([numel(m.primaseBindingLocations{1})-2 numel(m.primaseBindingLocations{2})-2], m.numLigations);
            assertEqual([false false], m.strandLigated);
            
            %Ex 5
            c.polymerizedRegions(:, :) = 0;
            c.strandBreaks(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(end)-holFtpt+corFtpt5+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(end)-corFtpt5 2]) = polGblIdx;
            assertEqual([numel(m.primaseBindingLocations{1})-1 numel(m.primaseBindingLocations{2})-1], m.numLigations);
            assertEqual([false false], m.strandLigated);
            
            %Ex 6
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, :) = size(c.sequence, 1);
            c.strandBreaks(:, :) = 0;
            c.strandBreaks(m.primaseBindingLocations{2}(3), 2) = 1;
            c.complexBoundSites(:, :) = 0;
            assertEqual([numel(m.primaseBindingLocations{1})+1 numel(m.primaseBindingLocations{2})], m.numLigations);
            assertEqual([true false], m.strandLigated);
        end
        
        function testStrandDuplicated(this)
            m = this.process;
            c = m.chromosome;
            
            %Ex 1
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, :) = size(c.sequence, 1);
            c.strandBreaks(:, :) = 0;
            c.strandBreaks(m.primaseBindingLocations{2}(3), 2) = 1;
            c.complexBoundSites(:, :) = 0;
            assertEqual([true false], m.strandDuplicated);
            
            %Ex 2
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, :) = size(c.sequence, 1);
            c.strandBreaks(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            assertEqual([true true], m.strandDuplicated);
            
            %Ex 3
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, 1:2) = size(c.sequence, 1);
            c.strandBreaks(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            assertEqual([false false], m.strandDuplicated);
            
            %Ex 4
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, 1:3) = size(c.sequence, 1);
            c.strandBreaks(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            assertEqual([false false], m.strandDuplicated);
            
            %Ex 5
            c.polymerizedRegions(:, :) = 0;
            c.polymerizedRegions(1, [1 3 4]) = size(c.sequence, 1);
            c.strandBreaks(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            assertEqual([false false], m.strandDuplicated);
        end
        
        function testHelicasePosition(this)
            m = this.process;
            c = m.chromosome;
            helGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            
            %Ex 1: nothing bound
            c.complexBoundSites(:, :) = 0;
            assertEqual([0 0], m.helicasePosition);
            
            %Ex 2
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([100 4; m.terCPosition+200 1; 300 3]) = helGblIdx;
            assertEqual([m.terCPosition+200 100], m.helicasePosition);
            
            %Ex 3
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([100 4; 300 3]) = helGblIdx;
            assertEqual([0 100], m.helicasePosition);
            
            %Ex 4
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([100 4; m.terCPosition+200 1; 300 4]) = helGblIdx;
            assertExceptionThrown(@() m.helicasePosition, 'Replication:error');
        end
        
        function testLeadingPolymerasePosition(this)
            m = this.process;
            c = m.chromosome;
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex);
            
            %Ex 1: nothing bound
            c.complexBoundSites(:, :) = 0;
            assertEqual([0 0], m.leadingPolymerasePosition);
            
            %Ex 2
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([100 2; m.terCPosition+200 1; 300 3]) = polGblIdx;
            assertEqual([m.terCPosition+200 0], m.leadingPolymerasePosition);
            
            %Ex 3
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+200 1; 300 4]) = polGblIdx;
            assertEqual([m.terCPosition+200 300], m.leadingPolymerasePosition);
            
            %Ex 4
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+200 1; 300 4; 1305 4]) = polGblIdx;
            assertExceptionThrown(@() m.leadingPolymerasePosition, 'Replication:error');
        end
        
        function testLaggingPolymerasePosition(this)
            m = this.process;
            c = m.chromosome;
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase);
            
            %Ex 1: nothing bound
            c.complexBoundSites(:, :) = 0;
            assertEqual([0 0], m.laggingPolymerasePosition);
            
            %Ex 2
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+300 3; m.terCPosition+1305 3; 400 1]) = polGblIdx;
            assertEqual([m.terCPosition+1305 0], m.laggingPolymerasePosition);
            
            %Ex 3
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([200 2; 400 2;m.terCPosition+300 3; m.terCPosition+1305 3; 400 1]) = polGblIdx;
            assertEqual([m.terCPosition+1305 200], m.laggingPolymerasePosition);
            
            %Ex 4
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([200 2; 400 2; 500 2;m.terCPosition+300 3; m.terCPosition+1305 3; 400 1]) = polGblIdx;
            assertExceptionThrown(@() m.laggingPolymerasePosition, 'Replication:error');
        end
        
        function testLeadingPosition(this)
            m = this.process;
            c = m.chromosome;
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_coreBetaClampGammaComplex);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            
            %Ex 1: nothing bound
            c.complexBoundSites(:, :) = 0;
            assertEqual([0 0], m.leadingPosition);
            
            %Ex 2
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([100 2; m.terCPosition+200 1; 300 3]) = polGblIdx;
            assertEqual([m.terCPosition+200+corFtpt5 0], m.leadingPosition);
            
            %Ex 3
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+200 1; 300 4]) = polGblIdx;
            assertEqual([m.terCPosition+200+corFtpt5 300+holFtpt-corFtpt5-1], m.leadingPosition);
            
            %Ex 4
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+200 1; 300 4; 1305 4]) = polGblIdx;
            assertExceptionThrown(@() m.leadingPosition, 'Replication:error');
        end
        
        function testLaggingPosition(this)
            m = this.process;
            c = m.chromosome;
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_coreBetaClampGammaComplex);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            
            %Ex 1: nothing bound
            c.complexBoundSites(:, :) = 0;
            assertEqual([0 0], m.laggingPosition);
            
            %Ex 2
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+300 3; m.terCPosition+1305 3; 400 1]) = polGblIdx;
            assertEqual([m.terCPosition+1305+holFtpt-corFtpt5-1 0], m.laggingPosition);
            
            %Ex 3
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([200 2; 400 2;m.terCPosition+300 3; m.terCPosition+1305 3; 400 1]) = polGblIdx;
            assertEqual([m.terCPosition+1305+holFtpt-corFtpt5-1 200+corFtpt5], m.laggingPosition);
            
            %Ex 4
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([200 2; 400 2; 500 2;m.terCPosition+300 3; m.terCPosition+1305 3; 400 1]) = polGblIdx;
            assertExceptionThrown(@() m.laggingPosition, 'Replication:error');
        end
        
        function testOkazakiFragment(this)
            m = this.process;
            c = m.chromosome;
            c.polymerizedRegions(:, :) = 0;
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_coreBetaClampGammaComplex);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            
            %Ex 1
            c.complexBoundSites(:, :) = 0;
            assertEqual([0 0], m.okazakiFragmentIndex);
            assertEqual([0 0], m.okazakiFragmentPosition);
            assertEqual([0 0], m.okazakiFragmentLength);
            assertEqual([0 0], m.okazakiFragmentProgress);
            
            %Ex 2
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-holFtpt+corFtpt5+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(3)-corFtpt5 2]) = polGblIdx;
            assertEqual([m.primaseBindingLocations{1}(1) m.primaseBindingLocations{2}(3)], m.laggingPosition);
            assertEqual([1 3], m.okazakiFragmentIndex);
            assertEqual([m.primaseBindingLocations{1}(1) m.primaseBindingLocations{2}(3)], m.okazakiFragmentPosition);
            assertEqual([size(c.sequence, 1) - m.primaseBindingLocations{1}(1)+1 diff(m.primaseBindingLocations{2}([2 3]))], m.okazakiFragmentLength);
            assertEqual([0 0], m.okazakiFragmentProgress);
            
            %Ex 3
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-holFtpt+corFtpt5+1-1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(3)-corFtpt5-1 2]) = polGblIdx;
            assertEqual([2 3], m.okazakiFragmentIndex);
            assertEqual([m.primaseBindingLocations{1}(2) m.primaseBindingLocations{2}(3)], m.okazakiFragmentPosition);
            assertEqual([diff(m.primaseBindingLocations{1}([2 1])) diff(m.primaseBindingLocations{2}([2 3]))], m.okazakiFragmentLength);
            assertEqual([diff(m.primaseBindingLocations{1}([2 1]))-1 1], m.okazakiFragmentProgress);
            
            %Ex 4
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-holFtpt+corFtpt5+1+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(3)-corFtpt5+1 2]) = polGblIdx;
            assertEqual([1 4], m.okazakiFragmentIndex);
            assertEqual([m.primaseBindingLocations{1}(1) m.primaseBindingLocations{2}(4)], m.okazakiFragmentPosition);
            assertEqual([size(c.sequence, 1) - m.primaseBindingLocations{1}(1)+1 diff(m.primaseBindingLocations{2}([3 4]))], m.okazakiFragmentLength);
            assertEqual([1 diff(m.primaseBindingLocations{2}([3 4]))-1], m.okazakiFragmentProgress);
            
            %Ex 5
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(3)-holFtpt+corFtpt5+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(1)-corFtpt5 2]) = polGblIdx;
            assertEqual([3 1], m.okazakiFragmentIndex);
            assertEqual([m.primaseBindingLocations{1}(3) m.primaseBindingLocations{2}(1)], m.okazakiFragmentPosition);
            assertEqual([diff(m.primaseBindingLocations{1}([3 2])) m.primaseBindingLocations{2}(1)], m.okazakiFragmentLength);
            assertEqual([0 0], m.okazakiFragmentProgress);
            
            %Ex 6
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(3)-holFtpt+corFtpt5+1-1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(1)-corFtpt5-1 2]) = polGblIdx;
            assertEqual([4 1], m.okazakiFragmentIndex);
            assertEqual([m.primaseBindingLocations{1}(4) m.primaseBindingLocations{2}(1)], m.okazakiFragmentPosition);
            assertEqual([diff(m.primaseBindingLocations{1}([4 3])) m.primaseBindingLocations{2}(1)], m.okazakiFragmentLength);
            assertEqual([diff(m.primaseBindingLocations{1}([4 3]))-1 1], m.okazakiFragmentProgress);
            
            %Ex 7
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(3)-holFtpt+corFtpt5+1+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(1)-corFtpt5+1 2]) = polGblIdx;
            assertEqual([3 2], m.okazakiFragmentIndex);
            assertEqual([m.primaseBindingLocations{1}(3) m.primaseBindingLocations{2}(2)], m.okazakiFragmentPosition);
            assertEqual([diff(m.primaseBindingLocations{1}([3 2])) diff(m.primaseBindingLocations{2}([1 2]))], m.okazakiFragmentLength);
            assertEqual([1 diff(m.primaseBindingLocations{2}([1 2]))-1], m.okazakiFragmentProgress);
            
            %Ex 8
            n1 = numel(m.primaseBindingLocations{1});
            n2 = numel(m.primaseBindingLocations{2});
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(end)-holFtpt+corFtpt5+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(end)-corFtpt5 2]) = polGblIdx;
            assertEqual([n1 n2], m.okazakiFragmentIndex);
            assertEqual([m.primaseBindingLocations{1}(n1) m.primaseBindingLocations{2}(n2)], m.okazakiFragmentPosition);
            assertEqual([diff(m.primaseBindingLocations{1}([n1 n1-1])) diff(m.primaseBindingLocations{2}([n2-1 n2]))], m.okazakiFragmentLength);
            assertEqual([0 0], m.okazakiFragmentProgress);
            
            %Ex 9
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(end)-holFtpt+corFtpt5+1+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(end)-corFtpt5+1 2]) = polGblIdx;
            assertEqual([n1 0], m.okazakiFragmentIndex);
            assertEqual([m.primaseBindingLocations{1}(n1) 0], m.okazakiFragmentPosition);
            assertEqual([diff(m.primaseBindingLocations{1}([n1 n1-1])) 0], m.okazakiFragmentLength);
            assertEqual([1 0], m.okazakiFragmentProgress);
            
            %Ex 10
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(end)-holFtpt+corFtpt5+1-1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(end)-corFtpt5-1 2]) = polGblIdx;
            assertEqual([0 n2], m.okazakiFragmentIndex);
            assertEqual([0 m.primaseBindingLocations{2}(n2)], m.okazakiFragmentPosition);
            assertEqual([0 diff(m.primaseBindingLocations{2}([n2-1 n2]))], m.okazakiFragmentLength);
            assertEqual([0 1], m.okazakiFragmentProgress);
            
            %Ex 11
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(3)-holFtpt+corFtpt5+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(6)-holFtpt+corFtpt5+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(4)-corFtpt5 2]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(7)-corFtpt5 2]) = polGblIdx;
            assertEqual([3 4], m.okazakiFragmentIndex);
            assertEqual([m.primaseBindingLocations{1}(3) m.primaseBindingLocations{2}(4)], m.okazakiFragmentPosition);
            assertEqual([diff(m.primaseBindingLocations{1}([3 2])) diff(m.primaseBindingLocations{2}([3 4]))], m.okazakiFragmentLength);
            assertEqual([0 0], m.okazakiFragmentProgress);
        end
        
        function testLeadingStrandBoundSSBs(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            ssbGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_ssb8mer);
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampGammaComplex);
            helGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            
            %Ex 1
            c.complexBoundSites(:, :) = 0;
            assertEqual(CircularSparseMat([],[],[size(c.sequence, 1) 2], 1), m.leadingStrandBoundSSBs);
            assertEqual([0 0], m.numLeadingTemplateBoundSSBs);
            
            %Ex 2
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+100 1; 200 4]) = helGblIdx;
            c.complexBoundSites([m.terCPosition+200 1; 100 4]) = polGblIdx;
            assertEqual(CircularSparseMat([],[],[size(c.sequence, 1) 2], 1), m.leadingStrandBoundSSBs);
            assertEqual([0 0], m.numLeadingTemplateBoundSSBs);
            
            %Ex 3
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+100 1; 200 4]) = helGblIdx;
            c.complexBoundSites([m.terCPosition+200 1; 100 4]) = polGblIdx;
            c.complexBoundSites([m.terCPosition+150 1; 150 4; 200 1; m.terCPosition+500 1; 50 4; 400 4; 130 1; 140 2; 160 3; m.terCPosition+130 2;m.terCPosition+140 3;m.terCPosition+160 4]) = ssbGblIdx;
            assertEqual(CircularSparseMat([m.terCPosition+150 1; 150 2],[1; 1],[size(c.sequence, 1) 2], 1), m.leadingStrandBoundSSBs);
            assertEqual([1 1], m.numLeadingTemplateBoundSSBs);
        end
        
        function testLaggingStrandBoundSSBs(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            c.polymerizedRegions(:, :) = 0;
            ssbGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_ssb8mer);
            polGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_coreBetaClampPrimase);
            helGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_helicase);
            holFtpt  = m.enzymeDNAFootprints(m.enzymeIndexs_coreBetaClampGammaComplex);
            corFtpt5 = m.enzymeDNAFootprints5Prime(m.enzymeIndexs_core);
            ssbSpcg = m.ssbComplexSpacing;
            ssbFtpt = m.enzymeDNAFootprints3Prime(m.enzymeIndexs_ssb8mer);
            
            %Ex 1
            c.complexBoundSites(:, :) = 0;
            assertEqual(CircularSparseMat([],[],[size(c.sequence, 1) 2], 1), m.laggingStrandBoundSSBs);
            assertEqual([0 0], m.numLaggingTemplateBoundSSBs);
            assertEqual([false false], m.areLaggingStrandSSBSitesBound);
            
            %Ex 2
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+10000 1; 5000 4]) = helGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-holFtpt+corFtpt5+1 4]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(1)-corFtpt5 1]) = polGblIdx;
            assertEqual(CircularSparseMat([],[],[size(c.sequence, 1) 2], 1), m.laggingStrandBoundSSBs);
            assertEqual([0 0], m.numLaggingTemplateBoundSSBs);
            assertEqual([false false], m.areLaggingStrandSSBSitesBound);
            
            %Ex 3
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.terCPosition+10000 1; 20000 4]) = helGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-holFtpt+corFtpt5+1 4]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(1)-corFtpt5 1]) = polGblIdx;
            c.complexBoundSites([
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-1 4;
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-2 4;
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-3 4;
                m.primaseBindingLocations{2}(1) + (ssbSpcg+ssbFtpt)*1 1;
                m.primaseBindingLocations{2}(1) + (ssbSpcg+ssbFtpt)*2 1;
                m.primaseBindingLocations{2}(1) + (ssbSpcg+ssbFtpt)*3 1;
                ]) = ssbGblIdx;
            assertEqual(CircularSparseMat([
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-1 1;
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-2 1;
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-3 1;
                m.primaseBindingLocations{2}(1) + (ssbSpcg+ssbFtpt)*1 2;
                m.primaseBindingLocations{2}(1) + (ssbSpcg+ssbFtpt)*2 2;
                m.primaseBindingLocations{2}(1) + (ssbSpcg+ssbFtpt)*3 2;
                ],[1; 1; 1; 1; 1; 1],[size(c.sequence, 1) 2], 1), m.laggingStrandBoundSSBs);
            assertEqual([3 3], m.numLaggingTemplateBoundSSBs);
            assertEqual([false false], m.areLaggingStrandSSBSitesBound);
            
            %Ex 4
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-1000 1; m.primaseBindingLocations{2}(1)+2000 4]) = helGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{1}(1)-holFtpt+corFtpt5+1 3]) = polGblIdx;
            c.complexBoundSites([m.primaseBindingLocations{2}(1)-corFtpt5 2]) = polGblIdx;
            c.complexBoundSites([
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-1 4;
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-2 4;
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-3 4;
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-4 4;
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-5 4;
                m.primaseBindingLocations{2}(1) + (ssbSpcg+ssbFtpt)*1 1;
                m.primaseBindingLocations{2}(1) + (ssbSpcg+ssbFtpt)*2 1;
                m.primaseBindingLocations{2}(1) + (ssbSpcg+ssbFtpt)*3 1;
                ]) = ssbGblIdx;
            assertEqual(CircularSparseMat([
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-1 1;
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-2 1;
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-3 1;
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-4 1;
                m.primaseBindingLocations{1}(1) + (ssbSpcg+ssbFtpt)*-5 1;
                m.primaseBindingLocations{2}(1) + (ssbSpcg+ssbFtpt)*1 2;
                m.primaseBindingLocations{2}(1) + (ssbSpcg+ssbFtpt)*2 2;
                m.primaseBindingLocations{2}(1) + (ssbSpcg+ssbFtpt)*3 2;
                ],[1; 1; 1; 1; 1; 1; 1; 1],[size(c.sequence, 1) 2], 1), m.laggingStrandBoundSSBs);
            assertEqual([5 3], m.numLaggingTemplateBoundSSBs);
            assertEqual([true false], m.areLaggingStrandSSBSitesBound);
        end
    end
    
    methods
        function testGeneEssentiality(this)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            
            %process
            m = this.process;
            c = m.chromosome;
            
            % initialize constants and state
            m.okazakiFragmentMeanLength         = 200;
            m.laggingBackupClampReloadingLength = 100;
            m.dnaPolymeraseElongationRate       = 100;
            m.primerLength                      = 11;
            m.ssbComplexSpacing                 = 30;
            m.startingOkazakiLoopLength         = 50;
            
            m.enzymeDNAFootprints(m.enzymeIndexs_ssb8mer) = 25;
            c.monomerDNAFootprints(m.enzymeMonomerGlobalIndexs) = m.enzymeDNAFootprints(m.enzymeMonomerLocalIndexs);
            c.complexDNAFootprints(m.enzymeComplexGlobalIndexs) = m.enzymeDNAFootprints(m.enzymeComplexLocalIndexs);
            
            bases = 'ACGT';
            seq = bases(randi(4, 2*4*m.okazakiFragmentMeanLength, 1));
            c.sequence = ChromosomeSequence(seq);
            c.sequenceLen = length(seq);
            c.sequenceGCContent = getGCContent(c.sequence);
            c.terCPosition = c.sequenceLen / 2;
            m.terCPosition = c.sequenceLen / 2;
            m.primaseBindingLocations = {
                size(c.sequence, 1) + 1 - m.calculatePrimaseBindingLocations(size(c.sequence, 1) - m.terCPosition)
                m.calculatePrimaseBindingLocations(m.terCPosition)
                }';
            m.dnaAFunctionalBoxStartPositions = c.sequenceLen - (40 + 10*(1:5)');
            
            c.initialize();
            c.linkingNumbers(1, :) = 0;
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_dntp) = [
                sum(sum(c.sequence(:, 2:3) == 'A'));
                sum(sum(c.sequence(:, 2:3) == 'C'));
                sum(sum(c.sequence(:, 2:3) == 'G'));
                sum(sum(c.sequence(:, 2:3) == 'T'))];
            m.substrates(m.substrateIndexs_atp) = ...
                + size(c.sequence, 1) ...                              %helicase unwinding each base pair
                + sum(cellfun(@numel, m.primaseBindingLocations)) ...  %lagging strand beta-clamp formation
                + 2;                                                   %leading strand beta-clamp formation
            m.substrates(m.substrateIndexs_water) = m.substrates(m.substrateIndexs_atp);
            m.substrates(m.substrateIndexs_nad) = sum(cellfun(@numel,m.primaseBindingLocations)) + 2;
            
            m.enzymes = m.enzymeComposition(:, [m.enzymeIndexs_replisome; m.enzymeIndexs_betaClamp; m.enzymeIndexs_ssb8mer]) * [2; 2; 100];
            m.enzymes(m.enzymeIndexs_ligase) = 100;
            m.boundEnzymes(:) = 0;
            
            this.setUpDnaAOriComplex();
            
            % evolve state
            this.helpTestGeneEssentiality(...
                {'MG_001'; %DNA polymerase III, beta subunit
                'MG_007'; %DNA polymerase III delta prime subunit, putative
                'MG_031'; %DNA polymerase III, alpha subunit
                'MG_091'; %single-strand binding protein family
                'MG_094'; %replicative DNA helicase
                'MG_254'; %DNA ligase, NAD-dependent
                'MG_250'; %DNA primase
                'MG_261'; %DNA polymerase III, alpha subunit
                'MG_315'; %DNA polymerase III, delta subunit, putative
                'MG_419'},... %DNA polymerase III, subunit gamma and tau
                @(m,i) all(m.strandDuplicated),...
                struct('lengthSec', 40));
        end
    end
    
    %helper methods
    methods
        function setUpDnaAOriComplex(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R1234), 1) = m.complexIndexs_DnaA_7mer_ATP;
            c.complexBoundSites(m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R5), 1) = m.complexIndexs_DnaA_1mer_ATP;
        end
    end
end
