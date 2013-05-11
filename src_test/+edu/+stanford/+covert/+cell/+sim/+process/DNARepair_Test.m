%DNA Repair process test case
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef DNARepair_Test < edu.stanford.covert.cell.sim.ReactionProcessTestCase
    %constants
    properties (Constant = true)
        expected_essentialGenes = {};
    end
    
    %constructor
    methods
        function this = DNARepair_Test(name)
            this = this@edu.stanford.covert.cell.sim.ReactionProcessTestCase(name);
        end
    end
    
    %simple test fixture
    methods
        function loadSimpleTestFixture(this)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            import edu.stanford.covert.cell.sim.state.Chromosome;
            
            %% process
            m = this.process;
            m.stepSizeSec = 1;
            
            %% chromosome state
            c = Chromosome([], []);
            m.states = {c};
            m.chromosome = c;
            c.dnaRepair = m;
            
            c.doubleStrandBreakSeparation = 1;
            c.relaxedBasesPerTurn = 10.5;
            c.equilibriumSuperhelicalDensity = -0.06;
            c.supercoiledSuperhelicalDensityTolerance = 0.10;
            
            seq = [
                'AGCAATTGCAGGCACGGGGATACGGTTCAATTGGCTCCCGGAGACGACAA' ...
                'TCCAACCCGCTTGCGCCTTAACGTGCTGTCCCGAGCGAAAAACAAAAGCC' ...
                'ATAGTGTCTCTGCGCATCACTCGCACCAATGAATTAATGCTACGTCTTCC' ...
                'TGCTTTAATGCTCAGTCCATGTACAGTCCGATCACCCAACACCGGCAAGC' ...
                'GTTGCGGCACTAAGCTGGCACGATATACCGGTCGGGCGCCTGGCAGCACC' ...
                'AACTGTCGATAGGATCAGCGCTGTATCTAGAATGTGAGCTTGGCGGTCAG' ...
                'CTATTCTTGAATGAATATTCTTAGGTCGAAGCCTGCATTCAGCACCGCGG' ...
                'CGCTGCTAATTCTTTAGATTCATGTCCCAGCGAATGTGCGCCCTTCCGTA' ...
                'CTGGTCAATTGGCGGAGAAGCAATTGAATGACCGGAACATCCATTACGAC' ...
                'CCTAATCGAAGCTGACAGTTACTAGCGCTAGACAATTACGTATAAGTCTG'];
            c.sequence = ChromosomeSequence(seq);
            c.sequenceLen = length(seq);
            c.sequenceGCContent = getGCContent(c.sequence);
            
            %% parameters
            m.HR_PolA_ResectionLength       = 8;
            m.HR_RecA_Spacing               = 3;
            m.HR_RuvAB_JunctionMigrationHop = 2.2;
            m.HR_RecU_CleavageSequence      = 'WTTS';
            m.HR_RecU_CleavagePosition      = 3;
            m.NER_UvrABC_IncisionMargin3    = 4;
            m.NER_UvrABC_IncisionMargin5    = 7;
            m.NER_PcrA_StepSize             = 1;
            m.RM_EcoD_RecognitionSequence   = 'TTANNNNNNNGTCY';
            m.RM_EcoD_MethylationPosition   = 3;
            m.RM_EcoD_RestrictionPosition   = [];
            m.RM_MunI_RecognitionSequence   = 'CAATTG';
            m.RM_MunI_MethylationPosition   = 3;
            m.RM_MunI_RestrictionPosition   = 1;
            
            %% substrates
            m.substrateWholeCellModelIDs = {
                'AD'; 'ADP';'AHCYS';'AMET';'AMP';'ATP';'CSN';'DAMP';'DATP';'DCMP';'DCTP';'DGMP';'DGTP';
                'DR5P';'DTMP';'DTTP';'DUMP';'FAPyAD';'FAPyGN';'FAPydA';'FAPydG';'FAPydAMP';
                'FAPydGMP';'GN';'H';'H2O';'NAD';'NMN';'PI';'PPI';'THY';'THY64CSN';'URA';'URI';'dApdAp';
                'dCpdCp';'dGpdGp';'dRibose5P_dRibose5P';'dT64dC';'dTMP64dCMP';'dTpdTp';
                'ho5URA';'ho5dU';'ho5dUMP';'m6AD';'m6dA';'m6dAMP';'oxo8AD';'oxo8GN';'oxo8dA';'oxo8dG';
                'oxo8dAMP';'oxo8dGMP'};
            m.substrateNames = m.substrateWholeCellModelIDs;
            
            m.substrateMolecularWeights = [
                135.1265; 424.1769; 384.4100; 399.4444; 345.2049; 503.1489; 111.1018; 329.2055;
                487.1495; 305.1808; 463.1248; 345.2049; 503.1489; 212.0942; 320.1921; 478.1361;
                306.1656; 153.1417; 169.1411; 269.2566; 285.2560; 347.2207; 363.2201; 151.1259;
                1.0079;  18.0152; 662.4162; 333.2107;  95.9793; 174.9513; 126.1131; 237.2149;
                112.0866; 244.2009; 641.4037; 593.3543; 673.4025; 407.1811; 391.1817; 469.4447;
                625.3729; 623.3769; 196.0948; 128.0860; 244.2009; 322.1650; 149.1530; 265.2679;
                343.2320; 151.1259; 167.1253; 267.2408; 283.2402; 345.2049; 361.2043];
            
            %indices
            m.substrateIndexs_water                = m.substrateIndexs({'H2O'});
            m.substrateIndexs_hydrogen             = m.substrateIndexs({'H'});
            m.substrateIndexs_DR5P                 = m.substrateIndexs({'DR5P'});
            m.substrateIndexs_dRibose5P_dRibose5P  = m.substrateIndexs({'dRibose5P_dRibose5P'});
            m.substrateIndexs_NAD                  = m.substrateIndexs({'NAD'});
            m.substrateIndexs_AMP                  = m.substrateIndexs({'AMP'});
            m.substrateIndexs_NMN                  = m.substrateIndexs({'NMN'});
            m.substrateIndexs_PPI                  = m.substrateIndexs({'PPI'});
            m.substrateIndexs_m6AD                 = m.substrateIndexs({'m6AD'});
            m.substrateIndexs_undamagedNucleobases = m.substrateIndexs({'AD';'CSN';'GN';'THY'});
            m.substrateIndexs_dNMPs                = m.substrateIndexs({'DAMP';'DCMP';'DGMP';'DTMP'});
            m.substrateIndexs_dNTPs                = m.substrateIndexs({'DATP';'DCTP';'DGTP';'DTTP'});
            m.substrateIndexs_dNpdNps              = m.substrateIndexs({'dApdAp';'dCpdCp';'dGpdGp';'dTpdTp'});
            m.substrateIndexs_damagedNucleobases   = m.substrateIndexs({
                'FAPyAD'; 'FAPyGN'; 'THY64CSN'; 'URA'; 'ho5URA'; 'm6AD'; 'oxo8AD'; 'oxo8GN'});
            m.substrateIndexs_damagedNucleosides = m.substrateIndexs({
                'FAPydA'; 'FAPydG'; 'URI'; 'dT64dC'; 'ho5dU'; 'm6dA'; 'oxo8dA'; 'oxo8dG'});
            m.substrateIndexs_damagedDNMPs = m.substrateIndexs({
                'DUMP'; 'FAPydAMP'; 'FAPydGMP'; 'dTMP64dCMP';'ho5dUMP'; 'm6dAMP'; 'oxo8dAMP'; 'oxo8dGMP'});
            
            m.substrateMetaboliteLocalIndexs = (1:numel(m.substrateWholeCellModelIDs))';
            m.substrateMetaboliteGlobalIndexs = m.substrateMetaboliteLocalIndexs;
            
            c.metabolite.dr5pIndexs = m.substrateIndexs({'DR5P'});
            c.metabolite.waterIndexs = m.substrateIndexs({'H2O'});
            c.metabolite.dnmpIndexs = m.substrateIndexs({'AD'; 'CSN'; 'GN'; 'THY'});
            c.metabolite.m6ADIndexs = m.substrateIndexs({'m6AD'});
            c.metabolite.molecularWeights = m.substrateMolecularWeights;
            
            %% enzymes
            m.enzymeWholeCellModelIDs = {
                'MG_001_MONOMER'; 'MG_097_MONOMER'; 'MG_190_MONOMER'; 'MG_235_MONOMER'; 'MG_254_MONOMER';
                'MG_262_MONOMER'; 'MG_498_MONOMER'; 'MG_339_MONOMER'; 'MG_438_MONOMER'; 'MG_073_206_421_TETRAMER';
                'MG_105_OCTAMER'; 'MG_184_DIMER'; 'MG_244_DIMER'; 'MG_352_DIMER'; 'MG_358_359_10MER'};
            m.enzymeNames = m.enzymeWholeCellModelIDs;
            m.enzymeMolecularWeights = 1e5 * [
                0.4429; 0.2819; 0.3642; 0.3243; 0.7540; 0.3314; 0.3265; 0.3480;
                0.4481; 3.5693; 1.8130; 0.7414; 1.6393; 0.3878; 3.0496];
            
            %indices
            m.enzymeIndexs_BER_baseExcision               = m.enzymeIndexs({'MG_097_MONOMER';'MG_498_MONOMER'});
            m.enzymeIndexs_formamidopyrimidineGlycosylase = m.enzymeIndexs({'MG_498_MONOMER'});          %BER  fpg          formamidopyrimidine-DNA glycosylase
            m.enzymeIndexs_uracilGlycosylase              = m.enzymeIndexs({'MG_097_MONOMER'});          %BER  ung          uracil-DNA glycosylase, putative
            m.enzymeIndexs_apurinicEndonuclease           = m.enzymeIndexs({'MG_235_MONOMER'});          %BER  nfo          apurinic endonuclease
            m.enzymeIndexs_incisionComplex                = m.enzymeIndexs({'MG_073_206_421_TETRAMER'}); %NER  uvrABC       DNA incision complex
            m.enzymeIndexs_helicase35                     = m.enzymeIndexs({'MG_244_DIMER'});            %NER  pcrA         3'-5' helicase
            m.enzymeIndexs_recombinationStrandExchange    = m.enzymeIndexs({'MG_339_MONOMER'});          %HR   recA         recombination protein, strand exchange
            m.enzymeIndexs_hollidayJunctionEndonuclease   = m.enzymeIndexs({'MG_352_DIMER'});            %HR   recU         Holliday junction endonuclease
            m.enzymeIndexs_hollidayJunctionHelicase       = m.enzymeIndexs({'MG_358_359_10MER'});        %HR   ruvAB        Holliday junction DNA helicase
            m.enzymeIndexs_exonuclease53                  = m.enzymeIndexs({'MG_262_MONOMER'});          %HR   pol I-like   5'-3' exonuclease, putative
            m.enzymeIndexs_polymerase                     = m.enzymeIndexs({'MG_001_MONOMER'});          %     dnaN         DNA polymerase III, beta subunit
            m.enzymeIndexs_ligase                         = m.enzymeIndexs({'MG_254_MONOMER'});          %     ligA         DNA ligase, NAD-dependent
            m.enzymeIndexs_phosphoesterase                = m.enzymeIndexs({'MG_190_MONOMER'});          %     mgpA         phosphoesterase
            m.enzymeIndexs_DisA                           = m.enzymeIndexs({'MG_105_OCTAMER'});          %     disA         DNA integrity scanning protein
            m.enzymeIndexs_RM_typeI                       = m.enzymeIndexs({'MG_438_MONOMER'});          %                  type I restriction modification DNA specificity domain protein
            m.enzymeIndexs_RM_typeII                      = m.enzymeIndexs({'MG_184_DIMER'});            %                  adenine-specific DNA modification methylase
            
            m.enzymeMonomerLocalIndexs =  (1:9)';
            m.enzymeComplexLocalIndexs  = (10:15)';
            m.enzymeMonomerGlobalIndexs = (1:9)';
            m.enzymeComplexGlobalIndexs = (1:6)';
            
            %protein release reactions
            c.monomerDNAFootprints = ones(size(m.enzymeMonomerLocalIndexs));
            c.complexDNAFootprints = ones(size(m.enzymeComplexLocalIndexs));
            c.monomerDNAFootprints(m.enzymeMonomerGlobalIndexs(m.enzymeMonomerLocalIndexs == m.enzymeIndexs_recombinationStrandExchange)) = m.HR_RecA_Spacing;
            
            c.monomerDNAFootprintBindingStrandedness = repmat(c.dnaStrandedness_dsDNA, size(m.enzymeMonomerLocalIndexs));
            c.complexDNAFootprintBindingStrandedness = repmat(c.dnaStrandedness_dsDNA, size(m.enzymeComplexLocalIndexs));
            c.monomerDNAFootprintRegionStrandedness = repmat(c.dnaStrandedness_dsDNA, size(m.enzymeMonomerLocalIndexs));
            c.complexDNAFootprintRegionStrandedness = repmat(c.dnaStrandedness_dsDNA, size(m.enzymeComplexLocalIndexs));
            c.monomerDNAFootprintBindingStrandedness(m.enzymeMonomerGlobalIndexs(m.enzymeMonomerLocalIndexs == m.enzymeIndexs_polymerase)) = c.dnaStrandedness_ssDNA;
            c.monomerDNAFootprintRegionStrandedness(m.enzymeMonomerGlobalIndexs(m.enzymeMonomerLocalIndexs == m.enzymeIndexs_polymerase)) = c.dnaStrandedness_xsDNA;
            
            nReactions = numel(m.enzymeMonomerGlobalIndexs) + numel(m.enzymeComplexGlobalIndexs);
            c.reactionBoundMonomer = zeros(nReactions, 1);
            c.reactionBoundComplex = zeros(nReactions, 1);
            c.reactionMonomerCatalysisMatrix = zeros(nReactions, numel(m.enzymeMonomerGlobalIndexs));
            c.reactionComplexCatalysisMatrix = zeros(nReactions, numel(m.enzymeComplexGlobalIndexs));
            
            c.reactionBoundComplex(1:nReactions-1) = ...
                m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_DisA);
            c.reactionMonomerCatalysisMatrix(1:numel(m.enzymeMonomerGlobalIndexs), :) = eye(numel(m.enzymeMonomerGlobalIndexs));
            c.reactionComplexCatalysisMatrix(numel(m.enzymeMonomerGlobalIndexs) + (1:numel(m.enzymeComplexGlobalIndexs)-1), ...
                m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs ~= m.enzymeIndexs_DisA)) = eye(numel(m.enzymeComplexGlobalIndexs)-1);
            
            c.reactionBoundMonomer(nReactions) = ...
                m.enzymeMonomerGlobalIndexs(m.enzymeMonomerLocalIndexs == m.enzymeIndexs_recombinationStrandExchange);
            c.reactionMonomerCatalysisMatrix(nReactions, m.enzymeMonomerGlobalIndexs(m.enzymeMonomerLocalIndexs == m.enzymeIndexs_polymerase)) = 1;
            
            c.reactionThresholds = ...
                sum(c.reactionMonomerCatalysisMatrix, 2) + ...
                sum(c.reactionComplexCatalysisMatrix, 2);
            
            %% reactions
            nSubstrates = length(m.substrateWholeCellModelIDs);
            nEnzymes    = length(m.enzymeWholeCellModelIDs);
            nReactions  = 32;
            m.reactionWholeCellModelIDs                = cell(nReactions,1);
            m.reactionTypes                            = cell(nReactions,1);
            m.reactionSmallMoleculeStoichiometryMatrix = zeros(nSubstrates, nReactions);
            m.reactionDNAStoichiometryMatrix           = zeros(nSubstrates, nReactions);
            m.reactionStoichiometryMatrix              = zeros(nSubstrates, nReactions);
            m.reactionDNABase                          = zeros(nReactions, 1);
            m.reactionCatalysisMatrix                  = zeros(nReactions, nEnzymes);
            m.enzymeBounds                             = repmat([-Inf Inf], nReactions, 1);
            
            %% base excision repair
            m.reactionWholeCellModelIDs{1}  = 'AP_endonuclease';
            m.reactionTypes{1} = 'base excision repair';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'H'}), 1) = [-1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dRibose5P_dRibose5P';'DR5P'}), 1) = [-1 2];
            m.reactionCatalysisMatrix(1, m.enzymeIndexs({'MG_235_MONOMER'})) = 1;
            m.enzymeBounds(1, 2) = 3.0800;
            
            m.reactionWholeCellModelIDs{2}  = 'AP_lyase';
            m.reactionTypes{2} = 'base excision repair';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'H'}), 2) = [-1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dRibose5P_dRibose5P';'DR5P'}), 2) = [-1 2];
            m.reactionCatalysisMatrix(2, m.enzymeIndexs({'MG_498_MONOMER'})) = 1;
            m.enzymeBounds(2, 2) = 8.9000e-004;
            
            m.reactionWholeCellModelIDs{3}  = 'deoxyribosephosphodiesterase';
            m.reactionTypes{3} = 'base excision repair';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O'; 'H'; 'DR5P'}), 3) = [-1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dRibose5P_dRibose5P'; 'DR5P'}), 3) = [-1 1];
            m.reactionCatalysisMatrix(3, m.enzymeIndexs({'MG_235_MONOMER'})) = 1;
            m.enzymeBounds(3, 2) = 0.0500;
                        
            m.reactionWholeCellModelIDs{4} = 'dRiboseP_lyase';
            m.reactionTypes{4} = 'base excision repair';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O'; 'H'; 'DR5P'}), 4) = [-1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dRibose5P_dRibose5P'; 'DR5P'}), 4) = [-1 1];
            m.reactionCatalysisMatrix(4, m.enzymeIndexs({'MG_001_MONOMER'})) = 1;
            m.enzymeBounds(4, 2) = 0.0750;
                       
            m.reactionWholeCellModelIDs{5} = 'MutM_FAPyAD';
            m.reactionTypes{5} = 'base excision repair, base excision';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'FAPyAD'}), 5) = [-1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'FAPydAMP';'DR5P'}), 5) = [-1 1];
            m.reactionDNABase(5,1) = m.substrateIndexs({'FAPydAMP'});
            m.reactionCatalysisMatrix(5, m.enzymeIndexs({'MG_498_MONOMER'})) = 1;
            m.enzymeBounds(5, 2) = 0.0014;
            
            m.reactionWholeCellModelIDs{6} = 'MutM_FAPyGN';
            m.reactionTypes{6} = 'base excision repair, base excision';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'FAPyGN'}), 6) = [-1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'FAPydGMP';'DR5P'}), 6) = [-1 1];
            m.reactionDNABase(6,1) = m.substrateIndexs({'FAPydGMP'});
            m.reactionCatalysisMatrix(6, m.enzymeIndexs({'MG_498_MONOMER'})) = 1;
            m.enzymeBounds(6, 2) = 0.0014;
            
            m.reactionWholeCellModelIDs{7} = 'MutM_oxo8AD';
            m.reactionTypes{7} = 'base excision repair, base excision';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'oxo8AD'}), 7) = [-1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'oxo8dAMP';'DR5P'}), 7) = [-1 1];
            m.reactionDNABase(7,1) = m.substrateIndexs({'oxo8dAMP'});
            m.reactionCatalysisMatrix(7, m.enzymeIndexs({'MG_498_MONOMER'})) = 1;
            m.enzymeBounds(7, 2) = 0.0014;
            
            m.reactionWholeCellModelIDs{8} = 'MutM_oxo8GN';
            m.reactionTypes{8} = 'base excision repair, base excision';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'oxo8GN'}), 8) = [-1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'oxo8dGMP';'DR5P'}), 8) = [-1 1];
            m.reactionDNABase(8,1) = m.substrateIndexs({'oxo8dGMP'});
            m.reactionCatalysisMatrix(8, m.enzymeIndexs({'MG_498_MONOMER'})) = 1;
            m.enzymeBounds(8, 2) = 0.0014;
            
            m.reactionWholeCellModelIDs{9} = 'Ung_ho5URA';
            m.reactionTypes{9} = 'base excision repair, base excision';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'ho5URA'}), 9) =  [-1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'ho5dUMP';'DR5P'}), 9) =  [-1 1];
            m.reactionDNABase(9,1) = m.substrateIndexs({'ho5dUMP'});
            m.reactionCatalysisMatrix(9, m.enzymeIndexs({'MG_097_MONOMER'})) = 1;
            m.enzymeBounds(9, 2) = 42;
            
            m.reactionWholeCellModelIDs{10} = 'Ung_URA';
            m.reactionTypes{10} = 'base excision repair, base excision';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'URA'}), 10) =  [-1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'DUMP';'DR5P'}), 10) =  [-1 1];
            m.reactionDNABase(10,1) = m.substrateIndexs({'DUMP'});
            m.reactionCatalysisMatrix(10, m.enzymeIndexs({'MG_097_MONOMER'})) = 1;
            m.enzymeBounds(10, 2) = 42;
            
            %% nucleotide excision repair
            m.reactionWholeCellModelIDs{11} = 'NER_DNAExcision';
            m.reactionTypes{11} = 'nucleotide excision repair';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'ATP';'H2O';'ADP';'H';'PI'}), 11) = [-1 -1 1 1 1];
            m.reactionCatalysisMatrix(11, m.enzymeIndexs({'MG_244_DIMER'})) = 1;
            m.enzymeBounds(11, 2) = 193;
            
            m.reactionWholeCellModelIDs{12} = 'NER_DNAIncision';
            m.reactionTypes{12} = 'nucleotide excision repair';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'ATP';'H2O';'ADP';'H';'PI'}), 12) = [-1 -3 1 3 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dRibose5P_dRibose5P';'DR5P'}), 12) = [-2 4];
            m.reactionCatalysisMatrix(12, m.enzymeIndexs({'MG_073_206_421_TETRAMER'})) = 1;
            m.enzymeBounds(12, 2) = 0.0083;
            
            %% cleavage
            m.reactionWholeCellModelIDs{13}  = 'DNAProcessiveCleavage_dAMP';
            m.reactionTypes{13} = 'DNA cleavage';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'H';'DAMP'}), 13) = [-1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dApdAp';'DAMP'}), 13) = [-1 1];
            m.reactionCatalysisMatrix(13, m.enzymeIndexs({'MG_190_MONOMER'})) = 1;
            m.enzymeBounds(13, 2) = 0.0740;
            
            m.reactionWholeCellModelIDs{14}  = 'DNAProcessiveCleavage_dCMP';
            m.reactionTypes{14} = 'DNA cleavage';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'H';'DCMP'}), 14) = [-1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dCpdCp';'DCMP'}), 14) = [-1 1];
            m.reactionCatalysisMatrix(14, m.enzymeIndexs({'MG_190_MONOMER'})) = 1;
            m.enzymeBounds(14, 2) = 0.0740;
            
            m.reactionWholeCellModelIDs{15}  = 'DNAProcessiveCleavage_dGMP';
            m.reactionTypes{15} = 'DNA cleavage';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'H';'DGMP'}), 15) = [-1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dGpdGp';'DGMP'}), 15) = [-1 1];
            m.reactionCatalysisMatrix(15, m.enzymeIndexs({'MG_190_MONOMER'})) = 1;
            m.enzymeBounds(15, 2) = 0.0740;
            
            m.reactionWholeCellModelIDs{16}  = 'DNAProcessiveCleavage_dTMP';
            m.reactionTypes{16} = 'DNA cleavage';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'H';'DTMP'}), 16) = [-1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dTpdTp';'DTMP'}), 16) = [-1 1];
            m.reactionCatalysisMatrix(16, m.enzymeIndexs({'MG_190_MONOMER'})) = 1;
            m.enzymeBounds(16, 2) = 0.0740;
            
            %% homologous recombination
            m.reactionWholeCellModelIDs{17} = 'HR_DNAResection_dAMP';
            m.reactionTypes{17} = 'homologous recombination';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'DAMP';'H'}), 17) = [-1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dApdAp';'DAMP'}), 17) = [-1 1];
            m.reactionCatalysisMatrix(17, m.enzymeIndexs({'MG_262_MONOMER'})) = 1;
            m.enzymeBounds(17, 2) = 0.1100;
            
            m.reactionWholeCellModelIDs{18} = 'HR_DNAResection_dCMP';
            m.reactionTypes{18} = 'homologous recombination';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'DCMP';'H'}), 18) = [-1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dCpdCp';'DCMP'}), 18) = [-1 1];
            m.reactionCatalysisMatrix(18, m.enzymeIndexs({'MG_262_MONOMER'})) = 1;
            m.enzymeBounds(18, 2) = 0.1100;
            
            m.reactionWholeCellModelIDs{19} = 'HR_DNAResection_dGMP';
            m.reactionTypes{19} = 'homologous recombination';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'DGMP';'H'}), 19) = [-1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dGpdGp';'DGMP'}), 19) = [-1 1];
            m.reactionCatalysisMatrix(19, m.enzymeIndexs({'MG_262_MONOMER'})) = 1;
            m.enzymeBounds(19, 2) = 0.1100;
            
            m.reactionWholeCellModelIDs{20} = 'HR_DNAResection_dTMP';
            m.reactionTypes{20} = 'homologous recombination';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'DTMP';'H'}), 20) = [-1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dTpdTp';'DTMP'}), 20) = [-1 1];
            m.reactionCatalysisMatrix(20, m.enzymeIndexs({'MG_262_MONOMER'})) = 1;
            m.enzymeBounds(20, 2) = 0.1100;
            
            m.reactionWholeCellModelIDs{21} = 'HR_junctionMigration';
            m.reactionTypes{21} = 'homologous recombination';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'ATP';'H2O';'ADP';'H';'PI'}), 21) = [-1 -1 1 1 1];
            m.reactionCatalysisMatrix(21, m.enzymeIndexs({'MG_358_359_10MER'})) = 1;
            m.enzymeBounds(21, 2) = 4.5454;
            
            m.reactionWholeCellModelIDs{22} = 'HR_junctionResolution';
            m.reactionTypes{22} = 'homologous recombination';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'H'}), 22) = [-2 2];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dRibose5P_dRibose5P';'DR5P'}), 22) = [-2 4];
            m.reactionCatalysisMatrix(22, m.enzymeIndexs({'MG_352_DIMER'})) = 1;
            m.enzymeBounds(22, 2) = 0.0035;
            
            m.reactionWholeCellModelIDs{23} = 'HR_strandExchange';
            m.reactionTypes{23} = 'homologous recombination';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'ATP';'H2O';'ADP';'H';'PI'}), 23) = [-1 -1 1 1 1];
            m.reactionCatalysisMatrix(23, m.enzymeIndexs({'MG_339_MONOMER'})) = 1;
            m.enzymeBounds(23, 2) = 0.0600;
            
            %% polymerization
            m.reactionWholeCellModelIDs{24} = 'DNA_polymerization_dATP_repair';
            m.reactionTypes{24} = 'DNA polymerization';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'DATP';'PPI'}), 24) = [-1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'DAMP';'dApdAp'}), 24) = [-1 1];
            m.reactionCatalysisMatrix(24, m.enzymeIndexs({'MG_001_MONOMER'})) = 1;
            m.enzymeBounds(24, 2) = 0.4500;
            
            m.reactionWholeCellModelIDs{25} = 'DNA_polymerization_dCTP_repair';
            m.reactionTypes{25} = 'DNA polymerization';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'DCTP';'PPI'}), 25) =[-1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'DCMP';'dCpdCp'}), 25) = [-1 1];
            m.reactionCatalysisMatrix(25, m.enzymeIndexs({'MG_001_MONOMER'})) = 1;
            m.enzymeBounds(25, 2) = 0.4500;
            
            m.reactionWholeCellModelIDs{26} = 'DNA_polymerization_dGTP_repair';
            m.reactionTypes{26} = 'DNA polymerization';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'DGTP';'PPI'}), 26) = [-1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'DGMP';'dGpdGp'}), 26) = [-1 1];
            m.reactionCatalysisMatrix(26, m.enzymeIndexs({'MG_001_MONOMER'})) = 1;
            m.enzymeBounds(26, 2) = 0.4500;
            
            m.reactionWholeCellModelIDs{27} = 'DNA_polymerization_dTTP_repair';
            m.reactionTypes{27} = 'DNA polymerization';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'DTTP';'PPI'}), 27) = [-1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'DTMP';'dTpdTp'}), 27) = [-1 1];
            m.reactionCatalysisMatrix(27, m.enzymeIndexs({'MG_001_MONOMER'})) = 1;
            m.enzymeBounds(27, 2) = 0.4500;
            
            %% ligation
            m.reactionWholeCellModelIDs{28}  = 'DNA_ligation_repair';
            m.reactionTypes{28} = 'DNA ligation';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'NAD';'AMP';'H';'NMN'}), 28) = [-1 1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'DR5P';'dRibose5P_dRibose5P'}), 28) = [-2 1];
            m.reactionCatalysisMatrix(28, m.enzymeIndexs({'MG_254_MONOMER'})) = 1;
            m.enzymeBounds(28, 2) = 0.0400;
            
            %% restriction/modification
            m.reactionWholeCellModelIDs{29} = 'DNA_RM_EcoD_Methylation';
            m.reactionTypes{29} = 'DNA restriction/modification';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'AMET';'AHCYS';'H'}), 29) = [-1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'AD';'m6AD'}), 29) = [-1 1];
            m.reactionCatalysisMatrix(29, m.enzymeIndexs({'MG_438_MONOMER'})) = 1;
            m.enzymeBounds(29, 2) = Inf;
            
            m.reactionWholeCellModelIDs{30} = 'DNA_RM_EcoD_Restriction';
            m.reactionTypes{30} = 'DNA restriction/modification';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'ATP';'H2O';'ADP';'H';'PI'}), 30) = [-1 -3 1 3 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dRibose5P_dRibose5P';'DR5P'}), 30) = [-2 4];
            m.reactionCatalysisMatrix(30, m.enzymeIndexs({'MG_438_MONOMER'})) = 1;
            m.enzymeBounds(30, 2) = Inf;
            
            m.reactionWholeCellModelIDs{31} = 'DNA_RM_MunI_Methylation';
            m.reactionTypes{31} = 'DNA restriction/modification';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'AMET';'AHCYS';'H'}), 31) = [-1 1 1];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'AD';'m6AD'}), 31) = [-1 1];
            m.reactionCatalysisMatrix(31, m.enzymeIndexs({'MG_184_DIMER'})) = 1;
            m.enzymeBounds(31, 2) =  0.0011;
            
            m.reactionWholeCellModelIDs{32} = 'DNA_RM_MunI_Restriction';
            m.reactionTypes{32} = 'DNA restriction/modification';
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O';'H'}), 32) = [-2 2];
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dRibose5P_dRibose5P';'DR5P'}), 32) = [-2 4];
            
            %names
            m.reactionNames  = m.reactionWholeCellModelIDs;
            
            %indices
            m.reactionIndexs_BER                   = find(strcmp(m.reactionTypes,'base excision repair') | strcmp(m.reactionTypes,'base excision repair, base excision'));
            m.reactionIndexs_BER_baseexcision      = find(strcmp(m.reactionTypes,'base excision repair, base excision'));
            m.reactionIndexs_NER                   = find(strcmp(m.reactionTypes,'nucleotide excision repair'));
            m.reactionIndexs_DNACleavage           = find(strcmp(m.reactionTypes,'DNA cleavage'));
            m.reactionIndexs_HR_dsbr               = find(strcmp(m.reactionTypes,'homologous recombination'));
            m.reactionIndexs_polymerization        = find(strcmp(m.reactionTypes,'DNA polymerization'));
            m.reactionIndexs_RM                    = find(strcmp(m.reactionTypes, 'DNA restriction/modification'));
            m.reactionIndexs_BER_aplyase           = m.reactionIndexs({'AP_lyase'});
            m.reactionIndexs_BER_drpesterase       = m.reactionIndexs({'deoxyribosephosphodiesterase'});
            m.reactionIndexs_BER_apendonuclease    = m.reactionIndexs({'AP_endonuclease'});
            m.reactionIndexs_BER_dr5plyase         = m.reactionIndexs({'dRiboseP_lyase'});
            m.reactionIndexs_NER_incision          = m.reactionIndexs({'NER_DNAIncision'});
            m.reactionIndexs_NER_excision          = m.reactionIndexs({'NER_DNAExcision'});
            m.reactionIndexs_HR_DNAResection       = m.reactionIndexs({'HR_DNAResection_dAMP';'HR_DNAResection_dCMP';'HR_DNAResection_dGMP';'HR_DNAResection_dTMP'});
            m.reactionIndexs_HR_strandExchange     = m.reactionIndexs({'HR_strandExchange'});
            m.reactionIndexs_HR_junctionMigration  = m.reactionIndexs({'HR_junctionMigration'});
            m.reactionIndexs_HR_junctionResolution = m.reactionIndexs({'HR_junctionResolution'});
            m.reactionIndexs_ligation              = m.reactionIndexs({'DNA_ligation_repair'});
            m.reactionIndexs_RM_EcoD_Methylation   = m.reactionIndexs({'DNA_RM_EcoD_Methylation'});
            m.reactionIndexs_RM_EcoD_Restriction   = m.reactionIndexs({'DNA_RM_EcoD_Restriction'});
            m.reactionIndexs_RM_MunI_Methylation   = m.reactionIndexs({'DNA_RM_MunI_Methylation'});
            m.reactionIndexs_RM_MunI_Restriction   = m.reactionIndexs({'DNA_RM_MunI_Restriction'});
            
            %% DNA
            m.hollidayJunctionResolutionSites = c.sequence.findSubsequence(m.HR_RecU_CleavageSequence, m.HR_RecU_CleavagePosition, true);
            assertFalse(size(m.hollidayJunctionResolutionSites,1) < 2);
            
            m.RM_MunI_MethylatedPositions = [
                m.RM_MunI_MethylationPosition;
                length(m.RM_MunI_RecognitionSequence)-m.RM_MunI_MethylationPosition+1]';
            m.RM_MunI_RestrictionPositions = [
                m.RM_MunI_RestrictionPosition;
                length(m.RM_MunI_RecognitionSequence)-m.RM_MunI_RestrictionPosition]';
            
            m.RM_EcoD_RecognitionSites = c.sequence.findSubsequence(m.RM_EcoD_RecognitionSequence, 1);
            m.RM_EcoD_RecognitionSites = mod(...
                repmat(m.RM_EcoD_RecognitionSites,1,length(m.RM_EcoD_RecognitionSequence)) + ...
                repmat(0:length(m.RM_EcoD_RecognitionSequence)-1,size(m.RM_EcoD_RecognitionSites,1),1) - ...
                1, size(c.sequence,1))+1;
            assertFalse(size(m.RM_EcoD_RecognitionSites,1) < 4);
            
            m.RM_MunI_RecognitionSites = c.sequence.findSubsequence(m.RM_MunI_RecognitionSequence, 1);
            m.RM_MunI_RecognitionSites = mod(...
                repmat(m.RM_MunI_RecognitionSites,1,length(m.RM_MunI_RecognitionSequence)) + ...
                repmat(0:length(m.RM_MunI_RecognitionSequence)-1,size(m.RM_MunI_RecognitionSites,1),1) - ...
                1, size(c.sequence,1))+1;
            assertFalse(size(m.RM_MunI_RecognitionSites,1) < 4);
            
            %% initial state
            m.substrates   = zeros(length(m.substrateWholeCellModelIDs), 1);
            m.enzymes      = zeros(length(m.enzymeWholeCellModelIDs),    1);
            m.boundEnzymes = zeros(length(m.enzymeWholeCellModelIDs),    1);
            
            c.initialize();
        end
    end
    
    %tests
    methods
        function testConstants(this)
            m = this.process;
            
            assertEqual((1:numel(m.substrateMetaboliteLocalIndexs))', m.substrateMetaboliteLocalIndexs);
            
            %test sorted so ismembc, ismembc2 can be used instead of ismember
            assertTrue(issorted(m.substrateGlobalIndexs));
        end
        
        function testReactionPropertiesLoadedCorrectly(this)
            m = this.process;
            
            substrateWholeCellModelIDs               = m.substrateWholeCellModelIDs;
            enzymeWholeCellModelIDs                  = m.enzymeWholeCellModelIDs;
            reactionWholeCellModelIDs                = m.reactionWholeCellModelIDs;
            reactionTypes                            = m.reactionTypes;
            reactionDNABase                          = m.reactionDNABase;
            reactionSmallMoleculeStoichiometryMatrix = m.reactionSmallMoleculeStoichiometryMatrix;
            reactionDNAStoichiometryMatrix           = m.reactionDNAStoichiometryMatrix;
            reactionCatalysisMatrix                  = m.reactionCatalysisMatrix;
            enzymeBounds                             = m.enzymeBounds;
            
            this.loadSimpleTestFixture();
            
            [tf, sbIdxs] = ismember(m.substrateWholeCellModelIDs, substrateWholeCellModelIDs);
            assertTrue(all(tf));
            
            [tf, ezIdxs] = ismember(m.enzymeWholeCellModelIDs, enzymeWholeCellModelIDs);
            assertTrue(all(tf));
            
            [tf, rxIdxs] = ismember(m.reactionWholeCellModelIDs, reactionWholeCellModelIDs);
            assertTrue(all(tf));
            
            for i = 1:length(rxIdxs)
                assertEqual(reactionTypes{rxIdxs(i)},                                    m.reactionTypes{i},                               sprintf('Error in definition of reaction %s',m.reactionWholeCellModelIDs{i}));
                assertEqual(reactionSmallMoleculeStoichiometryMatrix(sbIdxs,rxIdxs(i))', m.reactionSmallMoleculeStoichiometryMatrix(:,i)', sprintf('Error in definition of reaction %s',m.reactionWholeCellModelIDs{i}));
                assertEqual(reactionDNAStoichiometryMatrix(sbIdxs,rxIdxs(i))',           m.reactionDNAStoichiometryMatrix(:,i)',           sprintf('Error in definition of reaction %s',m.reactionWholeCellModelIDs{i}));
                assertEqual(reactionCatalysisMatrix(rxIdxs(i),ezIdxs)',                  m.reactionCatalysisMatrix(i,:)',                  sprintf('Error in definition of reaction %s',m.reactionWholeCellModelIDs{i}));
                
                assertEqual(enzymeBounds(rxIdxs(i), 1), m.enzymeBounds(i, 1), ...
                    sprintf('Error in definition of reaction %s',m.reactionWholeCellModelIDs{i}));
                assertElementsAlmostEqual(enzymeBounds(rxIdxs(i),2), m.enzymeBounds(i,2), ...
                    'relative', 5e-2, ...
                    sprintf('Error in definition of reaction %s', m.reactionWholeCellModelIDs{i}));
                
                if reactionDNABase(rxIdxs(i)) ~= 0 || m.reactionDNABase(i) ~= 0
                    assertTrue(reactionDNABase(rxIdxs(i))~=0 && reactionDNABase(rxIdxs(i))~=0,  ...
                        sprintf('Error in definition of reaction %s', m.reactionWholeCellModelIDs{i}));
                    assertEqual(...
                        substrateWholeCellModelIDs{reactionDNABase(rxIdxs(i))}, ...
                        m.substrateWholeCellModelIDs{m.reactionDNABase(i)},  ...
                        sprintf('Error in definition of reaction %s', m.reactionWholeCellModelIDs{i}));
                end
            end
        end
        
        %tests base excision repair (BER)
        function testBER(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            %example 1: nothing starts damaged, no additional damages
            this.clearDNADamages();
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.evolveState_BER();
            assertEqual(repmat(1e6, size(m.substrates)), m.substrates);
            assertEqual(repmat(1e6, size(m.enzymes)), m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 2: oxo8dAMP base excision
            this.clearDNADamages();
            c.damagedBases(100,1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_498_MONOMER'})) = 1e2;
            m.enzymeBounds(m.reactionIndexs_BER_baseexcision, 2) = 1e3;
            initialEnzymes = m.enzymes;
            m.evolveState_BER();
            assertEqual(max(0, m.reactionSmallMoleculeStoichiometryMatrix(:, m.reactionIndexs({'MutM_oxo8AD'}))), m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 3: no oxo8dAMP base excision b/c no enzyme
            this.clearDNADamages();
            c.damagedBases(100,1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:)=0;
            m.substrates(m.substrateIndexs_water)=1;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_097_MONOMER'}))=1e2;
            m.enzymeBounds(m.reactionIndexs_BER_aplyase,2)=1e3;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_BER();
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'})),[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 4: no oxo8dAMP base excision b/c no substrates
            this.clearDNADamages();
            c.damagedBases(100, 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:)=0;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_498_MONOMER'}))=1e2;
            m.enzymeBounds(m.reactionIndexs_BER_aplyase,2)=1e3;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_BER();
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'})),[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 5: ho5dUMP base excision
            this.clearDNADamages();
            c.damagedBases(100,1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'ho5dUMP'}));
            m.substrates(:)=0;
            m.substrates(m.substrateIndexs_water)=1;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_097_MONOMER'}))=1e6;
            initialEnzymes = m.enzymes;
            m.evolveState_BER();
            assertEqual(max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,m.reactionIndexs({'Ung_ho5URA'}))), m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 6: oxo8dAMP base excision, 3' AP lyase
            this.clearDNADamages();
            rxIdxs = m.reactionIndexs({'MutM_oxo8AD';'AP_lyase';});
            c.damagedBases(100,1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:)=0;
            m.substrates = sum(max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)),2);
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_498_MONOMER'}))=1e2;
            m.enzymeBounds(rxIdxs,2)=1;
            initialEnzymes = m.enzymes;
            m.evolveState_BER();
            assertEqual(sum(max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)),2), m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments],1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 7: oxo8dAMP base excision, 3' AP lyase, 5'-deoxyribosephosphodiesterase
            this.clearDNADamages();
            rxIdxs = m.reactionIndexs({'MutM_oxo8AD'; 'AP_lyase'; 'deoxyribosephosphodiesterase'});
            c.damagedBases(100,1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates = sum(max(0, -m.reactionSmallMoleculeStoichiometryMatrix(:, rxIdxs)), 2);
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_498_MONOMER';'MG_235_MONOMER'})) = 1e2;
            m.enzymeBounds(rxIdxs, 2) = 1e3;
            m.enzymeBounds(m.reactionIndexs('AP_endonuclease'), 2) = 0;
            initialEnzymes = m.enzymes;
            m.evolveState_BER();
            assertEqual(sum(max(0, m.reactionSmallMoleculeStoichiometryMatrix(:, rxIdxs)), 2), m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([99 1; 100 1],[1;1],[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 8: oxo8dAMP base excision, 5'-deoxyribosephosphodiesterase, no lyase
            this.clearDNADamages();
            rxIdxs = m.reactionIndexs({'MutM_oxo8AD'; 'AP_lyase'; 'deoxyribosephosphodiesterase'});
            c.damagedBases(100, 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates = sum(max(0, -m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)),2);
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_498_MONOMER'; 'MG_235_MONOMER'})) = 1e2;
            m.enzymeBounds(rxIdxs, 2) = 1;
            m.enzymeBounds(m.reactionIndexs({'AP_lyase'}), 2) = 0;
            m.enzymeBounds(m.reactionIndexs({'AP_endonuclease'}), 2) = 0;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_BER();
            assertEqual(initialSubstrates + m.reactionSmallMoleculeStoichiometryMatrix(:, rxIdxs) * [1; 0; 0], m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 9: oxo8dAMP base excision, 5'-AP endonuclease
            this.clearDNADamages();
            rxIdxs = m.reactionIndexs({'Ung_ho5URA';'AP_endonuclease'});
            c.damagedBases(100,1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'ho5dUMP'}));
            m.substrates = sum(max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)),2);
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_097_MONOMER';'MG_235_MONOMER'}))=1e2;
            m.enzymeBounds(rxIdxs,2)=1;
            initialEnzymes = m.enzymes;
            m.evolveState_BER();
            assertEqual(sum(max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)),2), m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([99 1],1,[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 10: oxo8dAMP base excision, 5'-AP endonuclease, 3'-(deoxyribose-5'-phosphate) lyase
            this.clearDNADamages();
            rxIdxs = m.reactionIndexs({'Ung_ho5URA'; 'AP_endonuclease'; 'dRiboseP_lyase'});
            c.damagedBases(100, 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'ho5dUMP'}));
            m.substrates(:) = 0;
            m.substrates = sum(max(0, -m.reactionSmallMoleculeStoichiometryMatrix(:, rxIdxs)), 2);
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_097_MONOMER'; 'MG_235_MONOMER'; 'MG_001_MONOMER'})) = 1e2;
            m.enzymeBounds(rxIdxs, 2) = 1;
            initialEnzymes = m.enzymes;
            m.evolveState_BER();
            assertEqual(sum(max(0, m.reactionSmallMoleculeStoichiometryMatrix(:, rxIdxs)), 2), m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([99 1;100 1],[1;1],[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 11: oxo8dAMP base excision,  3'-(deoxyribose-5'-phosphate) lyase, no endonuclease
            this.clearDNADamages();
            rxIdxs = m.reactionIndexs({'Ung_ho5URA';'AP_endonuclease';'dRiboseP_lyase'});
            c.damagedBases(100,1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'ho5dUMP'}));
            m.substrates(:)=0;
            m.substrates = sum(max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)),2);
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_097_MONOMER';'MG_001_MONOMER'}))=1e2;
            m.enzymeBounds(rxIdxs,2)=1;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_BER();
            assertEqual(initialSubstrates+m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)*[1;0;0], m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 12: 2nd strand oxo8dAMP base excision, 5'-AP endonuclease, 3'-(deoxyribose-5'-phosphate) lyase
            this.clearDNADamages();
            rxIdxs = m.reactionIndexs({'Ung_ho5URA';'AP_endonuclease';'dRiboseP_lyase'});
            c.damagedBases(100,2)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'ho5dUMP'}));
            m.substrates(:)=0;
            m.substrates = sum(max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)),2);
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_097_MONOMER';'MG_235_MONOMER';'MG_001_MONOMER'}))=1e2;
            m.enzymeBounds(rxIdxs,2)=1;
            initialEnzymes = m.enzymes;
            m.evolveState_BER();
            assertEqual(sum(max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)),2), m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([100 2],1,[size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([100 2],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([99 2;100 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 13: fairness
            count = zeros(2,1);
            rxIdx = m.reactionIndexs({'Ung_ho5URA'});
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_097_MONOMER'}))=1e2;
            m.enzymeBounds(rxIdx,2)=1;
            for i = 1:150
                this.clearDNADamages();
                c.damagedBases([100 1;200 2]) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'ho5dUMP'}));
                m.substrates = max(0, -m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
                m.evolveState_BER();
                
                assertEqual(1, nnz(c.damagedBases));
                
                positionStrands = find(c.damagedBases);
                count(positionStrands(:,2)) = count(positionStrands(:, 2)) + 1;
            end
            assertTrue(range(count) < 0.3 * max(count));
        end
        
        %tests nucleotide excision repair (NER)
        function testNER(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            %example 1: no damage
            this.clearDNADamages();
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 2: damaged base
            this.clearDNADamages();
            c.damagedBases(100, 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER';'MG_244_DIMER';'MG_190_MONOMER'})) = 1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            
            rxIdxs = [m.reactionIndexs_NER_incision; m.reactionIndexs_NER_excision; m.reactionIndexs_DNACleavage];
            rxCounts = [
                1;
                ceil((m.NER_UvrABC_IncisionMargin5 + m.NER_UvrABC_IncisionMargin3 + 1) / m.NER_PcrA_StepSize);
                c.sequence.subsequenceBaseCounts(100 + (-m.NER_UvrABC_IncisionMargin5:m.NER_UvrABC_IncisionMargin3), 1)];
            finalSubstrates = initialSubstrates + ...
                m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs) * rxCounts;
            finalSubstrates(m.substrateIndexs({'oxo8dAMP'})) = ...
                finalSubstrates(m.substrateIndexs({'oxo8dAMP'})) + 1;
            finalSubstrates(m.substrateIndexs_dNMPs) = ...
                finalSubstrates(m.substrateIndexs_dNMPs) - ...
                c.sequence.subsequenceBaseCounts(100, 1);
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([100+(-7:4)' ones(12,1)], ones(12,1), [size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([100+(-7:4)' ones(12,1)], ones(12,1), [size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([100+([-8 4])' ones(2,1)], ones(2,1), [size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 3: damaged base + damaged sugar phosphate outside excision region
            this.clearDNADamages();
            c.damagedBases(100,1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            c.damagedSugarPhosphates(200,1)=m.substrateMetaboliteGlobalIndexs(1);
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER';'MG_244_DIMER';'MG_190_MONOMER'}))=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            
            rxIdxs = [m.reactionIndexs_NER_incision; m.reactionIndexs_NER_excision; m.reactionIndexs_DNACleavage];
            rxCounts = [
                1;
                ceil((m.NER_UvrABC_IncisionMargin5 + m.NER_UvrABC_IncisionMargin3 + 1)/m.NER_PcrA_StepSize);
                c.sequence.subsequenceBaseCounts(100 + (-m.NER_UvrABC_IncisionMargin5:m.NER_UvrABC_IncisionMargin3), 1)];
            finalSubstrates = initialSubstrates + ...
                m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs) * rxCounts;
            finalSubstrates(m.substrateIndexs({'oxo8dAMP'})) = ...
                finalSubstrates(m.substrateIndexs({'oxo8dAMP'})) + 1;
            finalSubstrates(m.substrateIndexs_dNMPs) = ...
                finalSubstrates(m.substrateIndexs_dNMPs) - ...
                c.sequence.subsequenceBaseCounts(100,1);
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([100+(-7:4)' ones(12,1)], ones(12,1), [size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([100+(-7:4)' ones(12,1)], ones(12,1), [size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([200 1],m.substrateMetaboliteGlobalIndexs(1),[size(c.sequence,1) c.nCompartments],1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([100+([-8 4])' ones(2,1)], ones(2,1), [size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 4: damaged base + damaged sugar phosphate inside excision region
            this.clearDNADamages();
            c.damagedBases(100,1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            c.damagedSugarPhosphates(101,1)=m.substrateMetaboliteGlobalIndexs(1);
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER'; 'MG_244_DIMER'; 'MG_190_MONOMER'}))=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([101 1],m.substrateMetaboliteGlobalIndexs(1),[size(c.sequence,1) c.nCompartments],1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'})),[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 5: damaged base, insufficient metabolites
            this.clearDNADamages();
            c.damagedBases(100,1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:)=1e6;
            m.substrates(m.substrateIndexs_water)=0;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER';'MG_244_DIMER';'MG_190_MONOMER'}))=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'})),[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 6: damaged base, insufficient incision enzyme
            this.clearDNADamages();
            c.damagedBases(100,1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_244_DIMER';'MG_190_MONOMER'}))=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'})),[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 7: damaged base, insufficient excision enzyme
            this.clearDNADamages();
            c.damagedBases(100,1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER';'MG_190_MONOMER'}))=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'})),[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 8: damaged base, insufficient cleavage enzyme
            this.clearDNADamages();
            c.damagedBases(100,1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER';'MG_244_DIMER'}))=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'})),[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 9: damaged base, near origin, + strand
            this.clearDNADamages();
            c.damagedBases(2,1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER';'MG_244_DIMER';'MG_190_MONOMER'}))=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            
            rxIdxs = [m.reactionIndexs_NER_incision; m.reactionIndexs_NER_excision; m.reactionIndexs_DNACleavage];
            rxCounts = [
                1;
                ceil((m.NER_UvrABC_IncisionMargin5 + m.NER_UvrABC_IncisionMargin3 + 1)/m.NER_PcrA_StepSize);
                c.sequence.subsequenceBaseCounts(2 + (-m.NER_UvrABC_IncisionMargin5:m.NER_UvrABC_IncisionMargin3), 1)];
            finalSubstrates = initialSubstrates + ...
                m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs) * rxCounts;
            finalSubstrates(m.substrateIndexs({'oxo8dAMP'})) = ...
                finalSubstrates(m.substrateIndexs({'oxo8dAMP'})) + 1;
            finalSubstrates(m.substrateIndexs_dNMPs) = ...
                finalSubstrates(m.substrateIndexs_dNMPs) - ...
                c.sequence.subsequenceBaseCounts(2,1);
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([2+(-7:4)' ones(12,1)], ones(12,1), [size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([2+(-7:4)' ones(12,1)], ones(12,1), [size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([2+([-8 4])' ones(2,1)], ones(2,1), [size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 10: damaged base, near origin, - strand
            this.clearDNADamages();
            c.damagedBases(2,2)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER';'MG_244_DIMER';'MG_190_MONOMER'}))=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            
            rxIdxs = [m.reactionIndexs_NER_incision; m.reactionIndexs_NER_excision; m.reactionIndexs_DNACleavage];
            rxCounts = [
                1;
                ceil((m.NER_UvrABC_IncisionMargin5 + m.NER_UvrABC_IncisionMargin3 + 1)/m.NER_PcrA_StepSize);
                c.sequence.subsequenceBaseCounts(2 + (-m.NER_UvrABC_IncisionMargin3:m.NER_UvrABC_IncisionMargin5), 2)];
            finalSubstrates = initialSubstrates + ...
                m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs) * rxCounts;
            finalSubstrates(m.substrateIndexs({'oxo8dAMP'})) = ...
                finalSubstrates(m.substrateIndexs({'oxo8dAMP'})) + 1;
            finalSubstrates(m.substrateIndexs_dNMPs) = ...
                finalSubstrates(m.substrateIndexs_dNMPs) - ...
                c.sequence.subsequenceBaseCounts(2,2);
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([2+(-4:7)' repmat(2, 12,1)], ones(12,1), [size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([2+(-4:7)' repmat(2, 12,1)], ones(12,1), [size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([2+([-5 7])' repmat(2, 2,1)], ones(2,1), [size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 11: two damaged bases, insufficient resource to repair both
            counts = zeros(2,1);
            
            rxIdxs = [m.reactionIndexs_NER_incision; m.reactionIndexs_NER_excision; m.reactionIndexs_DNACleavage(1)];
            rxCounts = [
                1;
                ceil((m.NER_UvrABC_IncisionMargin5 + m.NER_UvrABC_IncisionMargin3 + 1)/m.NER_PcrA_StepSize);
                m.NER_UvrABC_IncisionMargin3 + m.NER_UvrABC_IncisionMargin5 + 1];
            initialSubstrates = max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)) * rxCounts;
            
            for i = 1:200
                this.clearDNADamages();
                c.damagedBases([100 1; 200 2]) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
                m.substrates = initialSubstrates;
                m.enzymes(:) = 0;
                m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER'; 'MG_244_DIMER'; 'MG_190_MONOMER'})) = 1e6;
                m.evolveState_NER();
                
                positionStrands = find(c.damagedBases);
                counts(positionStrands(:, 2)) = counts(positionStrands(:,2)) + 1;
            end
            assertTrue(range(counts) < 0.30 * max(counts));
            
            %example 12: intrastrand cross link
            this.clearDNADamages();
            c.intrastrandCrossLinks(100,1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'dTMP64dCMP'}));
            m.substrates(:) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER';'MG_244_DIMER';'MG_190_MONOMER'}))=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            
            rxIdxs = [m.reactionIndexs_NER_incision; m.reactionIndexs_NER_excision; m.reactionIndexs_DNACleavage];
            rxCounts = [
                1;
                ceil((m.NER_UvrABC_IncisionMargin5 + m.NER_UvrABC_IncisionMargin3 + 1)/m.NER_PcrA_StepSize);
                c.sequence.subsequenceBaseCounts(100 + (-m.NER_UvrABC_IncisionMargin5:m.NER_UvrABC_IncisionMargin3), 1)];
            finalSubstrates = initialSubstrates + ...
                m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs) * rxCounts;
            finalSubstrates(m.substrateIndexs({'dTMP64dCMP'})) = ...
                finalSubstrates(m.substrateIndexs({'dTMP64dCMP'})) + 1;
            finalSubstrates(m.substrateIndexs_dNMPs) = ...
                finalSubstrates(m.substrateIndexs_dNMPs) - ...
                c.sequence.subsequenceBaseCounts(100:101,1);
            finalSubstrates(m.substrateIndexs_water)=finalSubstrates(m.substrateIndexs_water)+1;
            finalSubstrates(m.substrateIndexs_hydrogen)=finalSubstrates(m.substrateIndexs_hydrogen)-1;
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([100+(-7:4)' ones(12,1)], ones(12,1), [size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([100+(-7:4)' ones(12,1)], ones(12,1), [size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([100+([-8 4])' ones(2,1)], ones(2,1), [size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 13: intrastrand cross link + gap site
            this.clearDNADamages();
            c.intrastrandCrossLinks(100,1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'dTMP64dCMP'}));
            c.gapSites(99,1)=1; %TODO: implement more precise damage filter; change to 100,1 and still pass
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER';'MG_244_DIMER';'MG_190_MONOMER'}))=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([99 1],1,[size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'dTMP64dCMP'})),[size(c.sequence,1) c.nCompartments],1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 14: intrastrand cross link + damaged base
            this.clearDNADamages();
            c.intrastrandCrossLinks(100,1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'dTMP64dCMP'}));
            c.damagedBases(99,1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'})); %TODO: implement more precise damage filter; change to 100,1 and still pass
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER';'MG_244_DIMER';'MG_190_MONOMER'}))=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_NER();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([99 1],m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'})),[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'dTMP64dCMP'})),[size(c.sequence,1) c.nCompartments],1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 15: several damaged bases and intrastrand cross links
            this.clearDNADamages();
            c.damagedBases((50:40:450)',1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:) = 1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs({'MG_073_206_421_TETRAMER';'MG_244_DIMER';'MG_190_MONOMER'}))=1e6;
            m.evolveState_NER();
            
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(12*11, nnz(c.gapSites));
            assertEqual(12*11, nnz(c.abasicSites));
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(2*11, nnz(c.strandBreaks));
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
        end
        
        %tests homologous recombination double strand break repair (HR-DSBR)
        function testHR(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            %Simulation only valid for doubleStrandBreakSeparation=1
            assertEqual(1, c.doubleStrandBreakSeparation, 'Simulation only valid for doubleStrandBreakSeparation equal 1');
            
            %example 1: no damage
            this.clearDNADamages();
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            m.evolveState_HR();
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 2: double strand break
            this.clearDNADamages();
            c.strandBreaks(100, 1:2) = 1;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.boundEnzymes(:) = 0;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_HR();
            
            rxIdxs = [m.reactionIndexs_HR_strandExchange; m.reactionIndexs_HR_DNAResection];
            rxCounts = [2; ...
                c.sequence.subsequenceBaseCounts(100+(1:8),1) + c.sequence.subsequenceBaseCounts(100+(-7:0),2)];
            finalSubstrates = initialSubstrates + m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)*rxCounts;
            
            finalEnzymes = initialEnzymes;
            finalBoundEnzymes = initialBoundEnzymes;
            finalEnzymes(m.enzymeIndexs_recombinationStrandExchange) = finalEnzymes(m.enzymeIndexs_recombinationStrandExchange)-2*ceil(8/3)+1;
            finalBoundEnzymes(m.enzymeIndexs_recombinationStrandExchange) = finalBoundEnzymes(m.enzymeIndexs_recombinationStrandExchange)+2*ceil(8/3)-1;
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(2*ceil(8/3)-1, nnz(c.monomerBoundSites) + nnz(c.complexBoundSites));
            assertEqual(ceil(8/3)-1, nnz(c.monomerBoundSites(100+(1:8)',1)) + nnz(c.complexBoundSites(100+(1:8)',1)));
            assertEqual(ceil(8/3)-1, nnz(c.monomerBoundSites(100+(-7:0)',2)) + nnz(c.complexBoundSites(100+(-7:0)',2)));
            assertEqual(CircularSparseMat([100+(1:8)' ones(8,1); 100+(-7:0)' repmat(2,8,1)],ones(16,1),[size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([100+(1:8)' ones(8,1); 100+(-7:0)' repmat(2,8,1)],ones(16,1),[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([100 1; 100 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([100+8 1; 100-8 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.hollidayJunctions);
            
            %example 3: double strand break, with other damage nearby
            this.clearDNADamages();
            c.strandBreaks(100, 1:2)=1;
            c.strandBreaks(102, 1:2)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_HR();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(initialBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([100 1; 100 2; 102 1;102 2],[1;1;1;1],[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 4: double strand break, insufficient metabolites
            this.clearDNADamages();
            c.strandBreaks(100, 1:2)=1;
            m.substrates(:)=1e6;
            m.substrates(m.substrateIndexs_water)=0;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_HR();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(initialBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([100 1; 100 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 5: double strand break, insufficient enzymes
            this.clearDNADamages();
            c.strandBreaks(100, 1:2)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.enzymes(m.enzymeIndexs_recombinationStrandExchange)=0;
            m.boundEnzymes(:)=0;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_HR();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(initialBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([100 1; 100 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 6: double strand break, near ORI
            this.clearDNADamages();
            c.strandBreaks(2, 1:2)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_HR();
            
            rxIdxs = [m.reactionIndexs_HR_strandExchange; m.reactionIndexs_HR_DNAResection];
            rxCounts = [2; ...
                c.sequence.subsequenceBaseCounts(2+(1:8),1) + c.sequence.subsequenceBaseCounts(2+(-7:0),2)];
            finalSubstrates = initialSubstrates + m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)*rxCounts;
            
            finalEnzymes = initialEnzymes;
            finalBoundEnzymes = initialBoundEnzymes;
            finalEnzymes(m.enzymeIndexs_recombinationStrandExchange) = finalEnzymes(m.enzymeIndexs_recombinationStrandExchange)-2*ceil(8/3)+1;
            finalBoundEnzymes(m.enzymeIndexs_recombinationStrandExchange) = finalBoundEnzymes(m.enzymeIndexs_recombinationStrandExchange)+2*ceil(8/3)-1;
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(finalBoundEnzymes(m.enzymeIndexs_recombinationStrandExchange), nnz(c.monomerBoundSites) + nnz(c.complexBoundSites));
            assertEqual(ceil(8/3)-1, nnz(c.monomerBoundSites(2+(1:8)',1)) + nnz(c.complexBoundSites(2+(1:8)',1)));
            assertEqual(ceil(8/3)-1, nnz(c.monomerBoundSites(2+(-7:0)',2)) + nnz(c.complexBoundSites(2+(-7:0)',2)));
            assertEqual(CircularSparseMat([2+(1:8)' ones(8,1); 2+(-7:0)' repmat(2,8,1)],ones(16,1),[size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([2+(1:8)' ones(8,1); 2+(-7:0)' repmat(2,8,1)],ones(16,1),[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([2 1;2 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([2+8 1; 2-8 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.hollidayJunctions);
            
            %example 7: holliday junctions
            this.clearDNADamages();
            c.hollidayJunctions([100+8 1; 100-8 2])=1;
            m.hollidayJunctionResolutionSites = [
                (10:20:190)' repmat(1, 10,1);
                (20:20:200)' repmat(2, 10,1)]; %#ok<RPMT1>
            m.HR_RuvAB_JunctionMigrationHop=2;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_HR();
            
            junctionMigrationHops = ...
                110-(100+8) + ...
                (100-8)-80;
            rxIdxs = [m.reactionIndexs_HR_strandExchange; m.reactionIndexs_HR_junctionMigration;  m.reactionIndexs_HR_junctionResolution];
            rxCounts = [2; junctionMigrationHops/m.HR_RuvAB_JunctionMigrationHop; 1];
            finalSubstrates = initialSubstrates + m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)*rxCounts;
            
            finalEnzymes = initialEnzymes;
            finalBoundEnzymes = initialBoundEnzymes;
            
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
            assertEqual(CircularSparseMat([110 1;80 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 8: holliday junctions, with incorrect spacing to be repaired
            this.clearDNADamages();
            c.hollidayJunctions([100+7 1; 100-8 2])=1;
            m.hollidayJunctionResolutionSites = [
                (10:20:190)' repmat(1, 10,1);
                (20:20:200)' repmat(2, 10,1)]; %#ok<RPMT1>
            m.HR_RuvAB_JunctionMigrationHop=2;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_HR();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(initialBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([107 1; 92 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.hollidayJunctions);
            
            %example 9: holliday junctions, with other damage nearby
            this.clearDNADamages();
            c.hollidayJunctions([100+8 1; 100-8 2])=1;
            c.damagedBases([103 1])=m.substrateMetaboliteGlobalIndexs(1);
            m.hollidayJunctionResolutionSites = [
                (10:20:190)' repmat(1, 10,1);
                (20:20:200)' repmat(2, 10,1)]; %#ok<RPMT1>
            m.HR_RuvAB_JunctionMigrationHop=2;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_HR();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(initialBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([103 1],m.substrateMetaboliteGlobalIndexs(1),[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([108 1; 92 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.hollidayJunctions);
            
            %example 10: holliday junctions, insufficient metabolites
            this.clearDNADamages();
            c.hollidayJunctions([100+8 1; 100-8 2])=1;
            m.hollidayJunctionResolutionSites = [
                (10:20:190)' repmat(1, 10,1);
                (20:20:200)' repmat(2, 10,1)]; %#ok<RPMT1>
            m.HR_RuvAB_JunctionMigrationHop=2;
            m.substrates(:)=1e6;
            m.substrates(m.substrateIndexs_water)=0;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_HR();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(initialBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([108 1; 92 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.hollidayJunctions);
            
            %example 12: holliday junctions, insufficient enzymes
            this.clearDNADamages();
            c.hollidayJunctions([100+8 1; 100-8 2])=1;
            m.hollidayJunctionResolutionSites = [
                (10:20:190)' repmat(1, 10,1);
                (20:20:200)' repmat(2, 10,1)]; %#ok<RPMT1>
            m.HR_RuvAB_JunctionMigrationHop=2;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.enzymes(m.enzymeIndexs_hollidayJunctionHelicase)=0;
            m.boundEnzymes(:)=0;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_HR();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(initialBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([108 1; 92 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.hollidayJunctions);
            
            %example 12: holliday junctions, near ORI
            this.clearDNADamages();
            dnaLength=size(c.sequence,1);
            site=dnaLength-2;
            c.hollidayJunctions([site+8 1; site-8 2])=1;
            m.hollidayJunctionResolutionSites = [
                (10:20:190)' repmat(1, 10,1);
                (20:20:200)' repmat(2, 10,1)]; %#ok<RPMT1>
            m.HR_RuvAB_JunctionMigrationHop=2;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_HR();
            
            junctionMigrationHops = ...
                10-(mod(site+8-1,dnaLength)+1) + ...
                (site-8)-200;
            rxIdxs = [m.reactionIndexs_HR_strandExchange; m.reactionIndexs_HR_junctionMigration;  m.reactionIndexs_HR_junctionResolution];
            rxCounts = [2; junctionMigrationHops/m.HR_RuvAB_JunctionMigrationHop; 1];
            finalSubstrates = initialSubstrates + m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)*rxCounts;
            
            finalEnzymes = initialEnzymes;
            finalBoundEnzymes = initialBoundEnzymes;
            
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
            assertEqual(CircularSparseMat([10 1;200 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
        end
        
        %tests polymerization and ligation
        function testPolymerizeLigate(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            %example 1: no damages
            this.clearDNADamages();
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
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
            
            %example 2: gap produced by BER
            this.clearDNADamages();
            c.abasicSites(100, 1) = 1;
            c.gapSites(100, 1) = 1;
            c.strandBreaks(99:100, 1) = 1;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_polymerase) = 1e6;
            m.boundEnzymes(:) = 0;
            
            rxIdxs = m.reactionIndexs_polymerization;
            rxCounts = c.sequence.subsequenceBaseCounts(100, 1);
            finalSubstrates = m.substrates + m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs) * rxCounts;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
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
            
            %example 3: gap produced by NER
            this.clearDNADamages();
            c.abasicSites(100:110,1)=1;
            c.gapSites(100:110,1)=1;
            c.strandBreaks([99;110],1)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_polymerase)=1e6;
            m.enzymeBounds(m.reactionIndexs_polymerization,2)=1e2;
            m.boundEnzymes(:)=0;
            
            rxIdxs = m.reactionIndexs_polymerization;
            rxCounts = c.sequence.subsequenceBaseCounts(100:110,1);
            finalSubstrates = m.substrates + m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)*rxCounts;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
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
            assertEqual(CircularSparseMat([110 1],1,[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 4: gap produced by HR
            this.clearDNADamages();
            
            c.abasicSites(100:110,1)=1;
            c.gapSites(100:110,1)=1;
            c.strandBreaks([99;110],1)=1;
            c.hollidayJunctions(110,1)=1;
            
            c.abasicSites(90:99,2)=1;
            c.gapSites(90:99,2)=1;
            c.strandBreaks([89;99],2)=1;
            c.hollidayJunctions(89,2)=1;
            
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_polymerase) = 1e6;
            m.enzymeBounds(m.reactionIndexs_polymerization, 2) = 1e2;
            m.boundEnzymes(:)=0;
            
            rxIdxs = m.reactionIndexs_polymerization;
            rxCounts =  ...
                c.sequence.subsequenceBaseCounts(100:110, 1) + ...
                c.sequence.subsequenceBaseCounts(90:99, 2);
            finalSubstrates = m.substrates + m.reactionSmallMoleculeStoichiometryMatrix(:, rxIdxs) * rxCounts;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
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
            assertEqual(CircularSparseMat([110 1;89 2],1,[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([110 1; 89 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.hollidayJunctions);
            
            %example 5: gap + nearby damage
            this.clearDNADamages();
            c.abasicSites(100,1)=1;
            c.gapSites(100,1)=1;
            c.strandBreaks(99:100,1)=1;
            c.damagedSugarPhosphates(99,1)=m.substrateMetaboliteGlobalIndexs(1);
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_polymerase)=1e6;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([99 1],m.substrateMetaboliteGlobalIndexs(1),[size(c.sequence,1) c.nCompartments],1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([99 1; 100 1],1,[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 6: gap, insufficient metabolites
            this.clearDNADamages();
            c.abasicSites(100,1)=1;
            c.gapSites(100,1)=1;
            c.strandBreaks(99:100,1)=1;
            m.substrates(:)=1e6;
            m.substrates(m.substrateIndexs_dNTPs)=0;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_polymerase)=1e6;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([99 1; 100 1],1,[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 7: gap, insufficient enzymes
            this.clearDNADamages();
            c.abasicSites(100,1)=1;
            c.gapSites(100,1)=1;
            c.strandBreaks(99:100,1)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_polymerase)=0;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([99 1; 100 1],1,[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 8: no single strand breaks
            this.clearDNADamages();
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
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
            
            %example 9: double strand break
            this.clearDNADamages();
            c.strandBreaks(100,1:2)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
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
            assertEqual(CircularSparseMat([100 1;100 2],[1;1],[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 10: single strand break
            this.clearDNADamages();
            c.strandBreaks(100,1)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_ligase)=1e6;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates + m.reactionSmallMoleculeStoichiometryMatrix(:,m.reactionIndexs_ligation);
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
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
            
            %example 11: single strand break + nearby gap
            this.clearDNADamages();
            c.strandBreaks(100,1)=1;
            c.gapSites(100,1)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_ligase)=1e6;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
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
            
            %example 12: single strand break, insufficient metabolites
            this.clearDNADamages();
            c.strandBreaks(100,1)=1;
            m.substrates(:)=1e6;
            m.substrates(m.substrateIndexs_NAD)=0;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_ligase)=1e6;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
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
            
            %example 13: single strand break, insufficient enzymes
            this.clearDNADamages();
            c.strandBreaks(100,1)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_ligase)=0;
            m.boundEnzymes(:)=0;
            
            finalSubstrates = m.substrates;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
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
            
            %example 14: multiple gaps
            this.clearDNADamages();
            
            c.abasicSites(100,1)=1;
            c.gapSites(100,1)=1;
            c.strandBreaks(99:100,1)=1;
            
            c.abasicSites(200:210,1)=1;
            c.gapSites(200:210,1)=1;
            c.strandBreaks([199;210],1)=1;
            
            c.abasicSites(300:310,1)=1;
            c.gapSites(300:310,1)=1;
            c.strandBreaks([299;310],1)=1;
            c.hollidayJunctions(310,1)=1;
            
            c.abasicSites(290:299,2)=1;
            c.gapSites(290:299,2)=1;
            c.strandBreaks([289;299],2)=1;
            c.hollidayJunctions(289,2)=1;
            
            m.substrates(:)=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_polymerase)=1e6;
            m.enzymes(m.enzymeIndexs_ligase)=1e6;
            m.boundEnzymes(:)=0;
            
            rxIdxs = [m.reactionIndexs_polymerization;m.reactionIndexs_ligation];
            rxCounts = [
                c.sequence.subsequenceBaseCounts(100,1) + ...
                c.sequence.subsequenceBaseCounts(200:210,1) + ...
                c.sequence.subsequenceBaseCounts(300:310,1) + ...
                c.sequence.subsequenceBaseCounts(290:299,2);
                4];
            finalSubstrates = m.substrates + m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)*rxCounts;
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
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
            assertEqual(CircularSparseMat([310 1; 289 2], [1; 1], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 15: two gaps, insufficient resources to repair both
            counts = zeros(2,1);
            m.enzymeBounds(m.reactionIndexs_ligation, 2) = 1;
            for i = 1:100
                this.clearDNADamages();
                
                c.abasicSites(100, 1) = 1;
                c.gapSites(100, 1) = 1;
                c.strandBreaks(99:100, 1) = 1;
                
                c.abasicSites(200, 2) = 1;
                c.gapSites(200, 2) = 1;
                c.strandBreaks(199:200, 2) = 1;
                
                m.substrates(:) = 1e6;
                m.enzymes(:) = 0;
                m.enzymes(m.enzymeIndexs_polymerase) = 1e6;
                m.enzymes(m.enzymeIndexs_ligase) = 1;
                m.boundEnzymes(:) = 0;
                
                m.evolveState_Polymerize();
                m.evolveState_Ligate();
                
                positionStrands = find(c.strandBreaks);
                assertEqual([1 2], size(positionStrands));
                
                counts(positionStrands(:, 2)) = counts(positionStrands(:, 2)) + 1;
            end
            assertTrue(range(counts) < 0.5 * max(counts));
        end
        
        function testRMSiteStrandBreakRepair(this)
            r = this.process;
            c = r.chromosome;
            
            c.initialize();
            r.initializeState();
            c.strandBreaks([282585 2]) = 1;
            c.invalidate();
            
            r.substrates(:) = 1e6;
            r.enzymes(:) = 1e2;
            r.enzymeBounds = 1e9 * r.enzymeBounds;
            for i = 1:10
                r.evolveState();
            end
            
            assertEqual(0, nnz(c.strandBreaks));
            assertEqual(0, nnz(c.abasicSites));
        end
        
        %tests DisA regulation
        function testDisA(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;            
            disAGlobalIdx = m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_DisA);
            [~, ~, footprint5Prime] = c.getDNAFootprint([], disAGlobalIdx);
            
            %example 1: no damage
            this.clearDNADamages();
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            m.evolveState_DisA();
            
            assertEqual(initialSubstrates, m.substrates);
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(initialBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 2: damaged sugar phosphate
            this.clearDNADamages();
            c.damagedSugarPhosphates(100,1)=m.substrateMetaboliteGlobalIndexs(1);
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.enzymes(m.enzymeIndexs_DisA)=1;
            m.boundEnzymes(:)=0;
            
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            finalSubstrates = initialSubstrates;
            finalEnzymes = initialEnzymes;
            finalBoundEnzymes = initialBoundEnzymes;
            
            finalEnzymes(m.enzymeIndexs_DisA) = 0;
            finalBoundEnzymes(m.enzymeIndexs_DisA) = 1;
            
            m.evolveState_DisA();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([100-footprint5Prime 1], disAGlobalIdx, [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(1), [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 3: damage is gap site, disA already bound
            this.clearDNADamages();
            c.gapSites(100,1)=1;
            c.complexBoundSites(100-footprint5Prime,1) = disAGlobalIdx;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.enzymes(m.enzymeIndexs_DisA)=1;
            m.boundEnzymes(:)=0;
            m.boundEnzymes(m.enzymeIndexs_DisA)=1;
            
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            finalSubstrates = initialSubstrates;
            finalEnzymes = initialEnzymes;
            finalBoundEnzymes = initialBoundEnzymes;
            
            m.evolveState_DisA();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([100-footprint5Prime 1],disAGlobalIdx,[size(c.sequence,1) c.nCompartments],1), c.complexBoundSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 4: damaged sugar phosphate, no DisA
            this.clearDNADamages();
            c.damagedSugarPhosphates(100,1)=m.substrateMetaboliteGlobalIndexs(1);
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.enzymes(m.enzymeIndexs_DisA)=0;
            m.boundEnzymes(:)=0;
            
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            finalSubstrates = initialSubstrates;
            finalEnzymes = initialEnzymes;
            finalBoundEnzymes = initialBoundEnzymes;
            
            finalEnzymes(m.enzymeIndexs_DisA)=0;
            finalBoundEnzymes(m.enzymeIndexs_DisA)=0;
            
            m.evolveState_DisA();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(1),[size(c.sequence,1) c.nCompartments],1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 5: damaged sugar phosphate, no metabolites
            this.clearDNADamages();
            c.damagedSugarPhosphates(100,1)=m.substrateMetaboliteGlobalIndexs(1);
            m.substrates(:)=0;
            m.enzymes(:)=1e6;
            m.enzymes(m.enzymeIndexs_DisA)=1;
            m.boundEnzymes(:)=0;
            
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            finalSubstrates = initialSubstrates;
            finalEnzymes = initialEnzymes;
            finalBoundEnzymes = initialBoundEnzymes;
            
            finalEnzymes(m.enzymeIndexs_DisA)=0;
            finalBoundEnzymes(m.enzymeIndexs_DisA)=1;
            
            m.evolveState_DisA();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([100-footprint5Prime 1],disAGlobalIdx,[size(c.sequence,1) c.nCompartments],1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(1),[size(c.sequence,1) c.nCompartments],1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 6: damage is methylation of restriction/modification site
            %DisA ignores m6AD
            this.clearDNADamages();
            m.RM_MunI_RecognitionSites(1, m.RM_MunI_MethylatedPositions(1)) = 100;
            c.damagedBases(100, 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD);
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.enzymes(m.enzymeIndexs_DisA)=1;
            m.boundEnzymes(:)=0;
            
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            finalSubstrates = initialSubstrates;
            finalEnzymes = initialEnzymes;
            finalBoundEnzymes = initialBoundEnzymes;
            
            finalEnzymes(m.enzymeIndexs_DisA) = 1;
            finalBoundEnzymes(m.enzymeIndexs_DisA) = 0;
            
            m.evolveState_DisA();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], disAGlobalIdx,[size(c.sequence,1) c.nCompartments],1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([100 1],m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD),[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 7: damage is to restriction/modification site
            this.clearDNADamages();
            m.RM_MunI_RecognitionSites(1,m.RM_MunI_MethylatedPositions(1))=100;
            c.damagedSugarPhosphates(100,1)=m.substrateIndexs_m6AD;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.enzymes(m.enzymeIndexs_DisA)=1;
            m.boundEnzymes(:)=0;
            
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            finalSubstrates = initialSubstrates;
            finalEnzymes = initialEnzymes;
            finalBoundEnzymes = initialBoundEnzymes;
            
            finalEnzymes(m.enzymeIndexs_DisA)=0;
            finalBoundEnzymes(m.enzymeIndexs_DisA)=1;
            
            m.evolveState_DisA();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([100-footprint5Prime 1],disAGlobalIdx,[size(c.sequence,1) c.nCompartments],1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([100 1],m.substrateIndexs_m6AD,[size(c.sequence,1) c.nCompartments],1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 8: multiple damages
            this.clearDNADamages();
            c.complexDNAFootprints(disAGlobalIdx)=1;
            footprint5Prime = 0;
            c.strandBreaks((99:100)',1)=1;
            c.gapSites(100,1)=1;
            c.abasicSites(100,1)=1;
            c.damagedSugarPhosphates([120;150;200],1)=m.substrateMetaboliteGlobalIndexs(1);
            c.damagedBases([50; 70; 190],2)=m.substrateMetaboliteGlobalIndexs(1);
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.enzymes(m.enzymeIndexs_DisA)=8;
            m.boundEnzymes(:)=0;
            
            initialSubstrates = m.substrates;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            finalSubstrates = initialSubstrates;
            finalEnzymes = initialEnzymes;
            finalBoundEnzymes = initialBoundEnzymes;
            
            finalEnzymes(m.enzymeIndexs_DisA)=0;
            finalBoundEnzymes(m.enzymeIndexs_DisA)=8;
            
            m.evolveState_DisA();
            
            assertEqual(finalSubstrates, m.substrates);
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            positionsStrands = [99 1;100 1; 120 1; 150 1; 200 1;50 2; 70 2;190 2];
            positionsStrands(:,1) = positionsStrands(:,1) - footprint5Prime;
            assertEqual(CircularSparseMat(positionsStrands, disAGlobalIdx,[size(c.sequence,1) c.nCompartments],1), c.complexBoundSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.gapSites);
            assertEqual(CircularSparseMat([100 1],1,[size(c.sequence,1) c.nCompartments],1), c.abasicSites);
            assertEqual(CircularSparseMat([120 1; 150 1; 200 1],m.substrateMetaboliteGlobalIndexs([1;1;1]),[size(c.sequence,1) c.nCompartments],1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([50 2; 70 2; 190 2],m.substrateMetaboliteGlobalIndexs([1;1;1]),[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([99 1;100 1],[1;1],[size(c.sequence,1) c.nCompartments],1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %example 9: multiple damages, insufficient DisA to bind all
            counts = zeros(2,1);
            for i = 1:100
                this.clearDNADamages();
                c.damagedSugarPhosphates([100 1; 200 2]) = m.substrateMetaboliteGlobalIndexs(1);
                m.substrates(:) = 1e6;
                m.enzymes(:) = 1e6;
                m.enzymes(m.enzymeIndexs_DisA) = 1;
                m.boundEnzymes(:) = 0;
                
                m.evolveState_DisA();
                
                positionStrands = find(c.complexBoundSites == disAGlobalIdx);
                counts(positionStrands(:,2)) = counts(positionStrands(:,2)) + 1;
            end
            assertTrue(range(counts) < 0.3*max(counts));
        end
        
        %tests restriction/modification systems
        function testRestrictionModification(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            %% setup
            rxIdxs = [m.reactionIndexs_RM_MunI_Methylation; m.reactionIndexs_RM_MunI_Restriction];
            m.RM_MunI_RecognitionSites = repmat(30*(1:15)',1,6)+repmat(1:6,15,1);
            m.RM_MunI_MethylatedPositions = [3 4];
            m.RM_MunI_RestrictionPositions =[1 5];
            m.enzymeBounds(rxIdxs,2)=1;
            
            this.clearDNADamages();
            
            %example 1
            c.damagedBases([m.RM_MunI_RecognitionSites(1, m.RM_MunI_MethylatedPositions)' [1;2]]) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %methylated
            
            %example 2
            c.damagedBases(m.RM_MunI_RecognitionSites(2, m.RM_MunI_MethylatedPositions(1)), 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %hemimethylated at first strand
            
            %example 3
            c.damagedBases(m.RM_MunI_RecognitionSites(3, m.RM_MunI_MethylatedPositions(2)), 2) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %hemimethylated at second strand
            
            %example 4
            c.damagedBases(m.RM_MunI_RecognitionSites(4, m.RM_MunI_MethylatedPositions(1)), 2) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %damaged at second strand
            
            %example 5
            c.damagedBases(m.RM_MunI_RecognitionSites(5, m.RM_MunI_MethylatedPositions(2)), 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %damaged at first strand
            
            %example 6
            c.damagedBases(m.RM_MunI_RecognitionSites(6, 1), 2) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %damaged at second strand
            
            %example 7
            c.damagedBases(m.RM_MunI_RecognitionSites(7, 1), 1) = m.substrateMetaboliteGlobalIndexs(1); %damaged at first strand
            
            %example 8
            c.gapSites(m.RM_MunI_RecognitionSites(8, 1), 1) = 1; %damaged at first strand
            
            %example 9
            c.abasicSites(m.RM_MunI_RecognitionSites(9, 1), 1) = 1; %damaged at first strand
            
            %example 10
            c.strandBreaks(m.RM_MunI_RecognitionSites(10, 1), 2) = 1; %damaged at second strand
            
            %example 11
            
            %example 12
            c.strandBreaks([m.RM_MunI_RecognitionSites(12,m.RM_MunI_RestrictionPositions)' [1;2]]) = 1; %cleaved
            
            %example 13
            c.strandBreaks([m.RM_MunI_RecognitionSites(13,m.RM_MunI_RestrictionPositions(1)) 1]) = 1; %cleaved at first strand
            
            %example 14
            c.strandBreaks([m.RM_MunI_RecognitionSites(14,m.RM_MunI_RestrictionPositions(2)) 2]) = 1; %cleaved at second strand
            
            %example 15
            c.strandBreaks([m.RM_MunI_RecognitionSites(15,m.RM_MunI_RestrictionPositions(2)) 1]) = 1; %cleaved at first strand
            
            initial_monomerBoundSites = c.monomerBoundSites;
            initial_complexBoundSites = c.complexBoundSites;
            initial_abasicSites = c.abasicSites;
            initial_gapSites = c.gapSites;
            initial_damagedSugarPhosphates = c.damagedSugarPhosphates;
            initial_damagedBases = c.damagedBases;
            initial_intrastrandCrossLinks = c.intrastrandCrossLinks;
            initial_strandBreaks = c.strandBreaks;
            initial_hollidayJunctions = c.hollidayJunctions;
            
            %% All resources
            c.monomerBoundSites = initial_monomerBoundSites;
            c.complexBoundSites = initial_complexBoundSites;
            c.abasicSites = initial_abasicSites;
            c.gapSites = initial_gapSites;
            c.damagedSugarPhosphates = initial_damagedSugarPhosphates;
            c.damagedBases = initial_damagedBases;
            c.intrastrandCrossLinks = initial_intrastrandCrossLinks;
            c.strandBreaks = initial_strandBreaks;
            c.hollidayJunctions = initial_hollidayJunctions;
            
            final_damagedBases = c.damagedBases;
            final_damagedBases(m.RM_MunI_RecognitionSites(2,m.RM_MunI_MethylatedPositions(2)), 2) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD);
            final_damagedBases(m.RM_MunI_RecognitionSites(3,m.RM_MunI_MethylatedPositions(1)), 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD);
            
            final_strandBreaks = c.strandBreaks;
            final_strandBreaks([m.RM_MunI_RecognitionSites(11,m.RM_MunI_RestrictionPositions)' [1;2]]) = 1;
            
            m.substrates(:)=0;
            m.substrates(m.substrateIndexs_water)=1e6;
            m.substrates(m.substrateIndexs({'AMET'}))=1e6;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_RM_typeII)=1e6;
            m.boundEnzymes(:)=0;
            
            final_substrates = m.substrates + m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdxs)*[2;2];
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            
            m.evolveState_Modification();
            m.evolveState_Restriction();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(initial_monomerBoundSites, c.monomerBoundSites);
            assertEqual(initial_complexBoundSites, c.complexBoundSites);
            assertEqual(initial_gapSites, c.gapSites);
            assertEqual(initial_abasicSites, c.abasicSites);
            assertEqual(initial_damagedSugarPhosphates, c.damagedSugarPhosphates);
            assertEqual(final_damagedBases, c.damagedBases);
            assertEqual(initial_intrastrandCrossLinks, c.intrastrandCrossLinks);
            assertEqual(final_strandBreaks, c.strandBreaks);
            assertEqual(initial_hollidayJunctions, c.hollidayJunctions);
            
            %% insufficient substrates
            c.monomerBoundSites = initial_monomerBoundSites;
            c.complexBoundSites = initial_complexBoundSites;
            c.abasicSites = initial_abasicSites;
            c.gapSites = initial_gapSites;
            c.damagedSugarPhosphates = initial_damagedSugarPhosphates;
            c.damagedBases = initial_damagedBases;
            c.intrastrandCrossLinks = initial_intrastrandCrossLinks;
            c.strandBreaks = initial_strandBreaks;
            c.hollidayJunctions = initial_hollidayJunctions;
            
            m.substrates(:)=0;
            m.enzymes(:)=0;
            m.enzymes(m.enzymeIndexs_RM_typeII)=1e6;
            m.boundEnzymes(:)=0;
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            
            m.evolveState_Modification();
            m.evolveState_Restriction();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(initial_monomerBoundSites, c.monomerBoundSites);
            assertEqual(initial_complexBoundSites, c.complexBoundSites);
            assertEqual(initial_gapSites, c.gapSites);
            assertEqual(initial_abasicSites, c.abasicSites);
            assertEqual(initial_damagedSugarPhosphates, c.damagedSugarPhosphates);
            assertEqual(initial_damagedBases, c.damagedBases);
            assertEqual(initial_intrastrandCrossLinks, c.intrastrandCrossLinks);
            assertEqual(initial_strandBreaks, c.strandBreaks);
            assertEqual(initial_hollidayJunctions, c.hollidayJunctions);
            
            %% insufficient enzymes
            c.monomerBoundSites = initial_monomerBoundSites;
            c.complexBoundSites = initial_complexBoundSites;
            c.abasicSites = initial_abasicSites;
            c.gapSites = initial_gapSites;
            c.damagedSugarPhosphates = initial_damagedSugarPhosphates;
            c.damagedBases = initial_damagedBases;
            c.intrastrandCrossLinks = initial_intrastrandCrossLinks;
            c.strandBreaks = initial_strandBreaks;
            c.hollidayJunctions = initial_hollidayJunctions;
            
            m.substrates(:)=0;
            m.substrates(m.substrateIndexs_water)=1e6;
            m.substrates(m.substrateIndexs({'AMET'}))=1e6;
            m.enzymes(:)=0;
            m.boundEnzymes(:)=0;
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            
            m.evolveState_Modification();
            m.evolveState_Restriction();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(initial_monomerBoundSites, c.monomerBoundSites);
            assertEqual(initial_complexBoundSites, c.complexBoundSites);
            assertEqual(initial_gapSites, c.gapSites);
            assertEqual(initial_abasicSites, c.abasicSites);
            assertEqual(initial_damagedSugarPhosphates, c.damagedSugarPhosphates);
            assertEqual(initial_damagedBases, c.damagedBases);
            assertEqual(initial_intrastrandCrossLinks, c.intrastrandCrossLinks);
            assertEqual(initial_strandBreaks, c.strandBreaks);
            assertEqual(initial_hollidayJunctions, c.hollidayJunctions);
        end
        
        %test NER + polymerization/ligation
        function testAllRepair(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            disAGlobalIdx = m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_DisA);
            
            %% example 1: multiple double strand breaks, BER +
            %polymerization/ligation
            %example 2: oxo8dAMP base excision
            this.clearDNADamages();
            c.damagedBases((50:50:400)', 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.enzymeBounds(m.reactionIndexs_polymerization, 2) = 1e2;
            
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_BER();
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            m.evolveState_BER();
            m.evolveState_Polymerize();
            m.evolveState_Ligate();            
            
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
            
            %% example 2: multiple double strand breaks, NER +
            %polymerization/ligation
            this.clearDNADamages();
            c.damagedBases((50:50:400)',1)=m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.enzymeBounds(m.reactionIndexs_polymerization,2)=1e2;
            
            finalEnzymes = m.enzymes;
            finalBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_NER();
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
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
            
            %% example 3: multiple double strand breaks, HR + polymerization/ligation + HR + polymerization/ligation
            this.clearDNADamages();
            m.hollidayJunctionResolutionSites = [
                (60:50:410)' repmat(1, 8,1);
                (40:50:390)' repmat(2, 8,1)]; %#ok<RPMT1>
            m.HR_RuvAB_JunctionMigrationHop=2;
            m.enzymeBounds(m.reactionIndexs_polymerization, 2) = 1e2;
            c.strandBreaks((50:50:400)', 1:2)=1;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.boundEnzymes(:)=0;
            initialEnzymes = m.enzymes;
            initialBoundEnzymes = m.boundEnzymes;
            
            m.evolveState_HR();
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            m.evolveState_HR();
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
            assertEqual(initialEnzymes, m.enzymes);
            assertEqual(initialBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            % example 4: two double strand breaks, HR + polymerization/ligation + HR + polymerization/ligation, insufficient resources to repair both
            counts = zeros(2,1);
            m.hollidayJunctionResolutionSites = [
                110 1;
                210 1;
                80 2;
                180 2];
            m.HR_RuvAB_JunctionMigrationHop = 2;
            m.enzymeBounds(m.reactionIndexs_ligation, 2) = 1;
            for i = 1:50
                this.clearDNADamages();
                c.strandBreaks([100; 200], 1:2)=1;
                m.substrates(:) = 1e6;
                m.enzymes(:) =1e6;
                m.enzymes(m.enzymeIndexs_ligase) = 3;
                m.boundEnzymes(:) = 0;
                initialEnzymes = m.enzymes;
                initialBoundEnzymes = m.boundEnzymes;
                
                m.evolveState_HR();
                m.evolveState_Polymerize();
                m.evolveState_Ligate();
                m.evolveState_HR();
                m.evolveState_Polymerize();
                m.evolveState_Ligate();
                
                assertEqual(initialEnzymes, m.enzymes);
                assertEqual(initialBoundEnzymes, m.boundEnzymes);
                
                damagedSites = c.getDamagedSites(true, true, false, false, false, false, true);
                damagedSites = collapse(reshape(damagedSites([75:115 175:215], 1:2), 41, 2, 2), -2) ~= 0;
                if ~any(damagedSites)
                    throw(MException('Chromosome:error','At least one double strand break shouldn''t be completely repaired'));
                end
                
                counts = counts + damagedSites;
            end
            
            assertTrue(range(counts) < 0.3 * max(counts));
            
            %% example 5: unmethylated restriction/modification site -> cleaved -> repaired
            m = this.process;
            
            c.complexDNAFootprints(:) = 1;
            c.monomerDNAFootprints(:) = 1;
            
            rxIdxs = [m.reactionIndexs_RM_MunI_Methylation; m.reactionIndexs_RM_MunI_Restriction];
            m.RM_MunI_RecognitionSites = repmat(30*(1:6)',1,6)+repmat(1:6,6,1);
            m.RM_MunI_MethylatedPositions = [3 4];
            m.RM_MunI_RestrictionPositions =[1 5];
            m.enzymeBounds(rxIdxs,2)=1;
            m.enzymeBounds(m.reactionIndexs_polymerization, 2) = 1e2;
            
            this.clearDNADamages();
            c.damagedBases([m.RM_MunI_RecognitionSites(1,m.RM_MunI_MethylatedPositions)' [1;2]]) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %methylated
            c.damagedBases(m.RM_MunI_RecognitionSites(2,m.RM_MunI_MethylatedPositions(1)), 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %hemimethylated at first strand
            c.damagedBases(m.RM_MunI_RecognitionSites(3,m.RM_MunI_MethylatedPositions(2)), 2) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %hemimethylated at second strand
            c.strandBreaks([m.RM_MunI_RecognitionSites(4,m.RM_MunI_RestrictionPositions)' [1;2]]) = 1; %cleaved
            c.damagedBases(m.RM_MunI_RecognitionSites(6,1), 1) = m.substrateMetaboliteGlobalIndexs(1); %damaged
            
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.boundEnzymes(:) = 0;            
            
            final_enzymes = m.enzymes;
            final_enzymes(m.enzymeIndexs_DisA) = m.enzymes(m.enzymeIndexs_DisA);
            final_boundEnzymes = m.boundEnzymes;
            final_boundEnzymes(m.enzymeIndexs_DisA) = m.boundEnzymes(m.enzymeIndexs_DisA);
            
            m.evolveState_Modification();
            m.evolveState_Restriction();
            m.evolveState_DisA();
            assertEqual(CircularSparseMat([33 1; 63 1; 93 1;34 2;64 2;94 2;181 1],[repmat(m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD),6,1); 9],[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([121 1;151 1;181 1;125 2;155 2],repmat(disAGlobalIdx,5,1),[size(c.sequence,1) c.nCompartments],1), c.complexBoundSites);
            m.evolveState_NER();
            m.evolveState_Polymerize();
            m.evolveState_Ligate();
            
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([33 1; 63 1; 93 1;34 2;64 2;94 2],repmat(m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD),6,1),[size(c.sequence,1) c.nCompartments],1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
            
            %% example 6: all kinds of damage
            m = this.process;
            m.seedRandStream();
            
            c.complexDNAFootprints(:) = 1;
            c.monomerDNAFootprints(:) = 1;
            
            m.hollidayJunctionResolutionSites = [
                (45:50:445)' repmat(1, 9, 1);
                (25:50:425)' repmat(2, 9, 1)]; %#ok<RPMT1>
            m.HR_RuvAB_JunctionMigrationHop = 2;
            rxIdxs = [m.reactionIndexs_RM_MunI_Methylation; m.reactionIndexs_RM_MunI_Restriction];
            m.RM_MunI_RecognitionSites = repmat(20*(1:6)',1,6)+repmat(1:6,6,1);
            m.RM_MunI_MethylatedPositions = [3 4];
            m.RM_MunI_RestrictionPositions =[1 5];
            
            this.clearDNADamages();
            c.damagedBases((200:50:400)', 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            c.damagedBases((220:50:470)', 1) = m.substrateMetaboliteGlobalIndexs(1);
            c.strandBreaks((235:50:435)', 1:2) = 1;
            c.damagedBases([m.RM_MunI_RecognitionSites(1,m.RM_MunI_MethylatedPositions)' [1;2]]) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %methylated
            c.damagedBases(m.RM_MunI_RecognitionSites(2,m.RM_MunI_MethylatedPositions(1)), 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %hemimethylated at first strand
            c.damagedBases(m.RM_MunI_RecognitionSites(3,m.RM_MunI_MethylatedPositions(2)), 2) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %hemimethylated at second strand
            c.strandBreaks([m.RM_MunI_RecognitionSites(4,m.RM_MunI_RestrictionPositions)' [1; 2]]) = 1; %cleaved
            c.damagedBases(m.RM_MunI_RecognitionSites(6,1), 1) = m.substrateMetaboliteGlobalIndexs(1); %damaged
            
            m.enzymeBounds(rxIdxs, 2) = 1;
            m.enzymeBounds(m.reactionIndexs_polymerization, 2) = 1e2;
            
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.boundEnzymes(:) = 0;
            m.enzymeBounds(m.reactionIndexs_polymerization, 2) = 1e2;
            
            finalEnzymes = m.enzymes;
            finalEnzymes(m.enzymeIndexs_DisA) = m.enzymes(m.enzymeIndexs_DisA);
            finalBoundEnzymes = m.boundEnzymes;
            finalBoundEnzymes(m.enzymeIndexs_DisA) = m.boundEnzymes(m.enzymeIndexs_DisA);
            
            for i = 1:6
                m.evolveState();
            end
            
            assertEqual(finalEnzymes, m.enzymes);
            assertEqual(finalBoundEnzymes, m.boundEnzymes);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.complexBoundSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.gapSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.abasicSites);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([23 1; 43 1; 24 2; 44 2; 63 1; 64 2], ...
                repmat(m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD), 6, 1), [size(c.sequence, 1) c.nCompartments], 1), c.damagedBases);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.strandBreaks);
            assertEqual(CircularSparseMat([], [], [size(c.sequence,1) c.nCompartments], 1), c.hollidayJunctions);
        end
        
        function testGeneEssentiality(this)
            %% process
            m = this.process;
            c = m.chromosome;
            
            c.monomerDNAFootprints(:) = 1;
            c.complexDNAFootprints(:) = 1;
            
            m.hollidayJunctionResolutionSites = [
                (160:50:560)' repmat(1, 9, 1);
                (140:50:540)' repmat(2, 9, 1)]; %#ok<RPMT1>
            m.HR_RuvAB_JunctionMigrationHop = 2;
            
            m.RM_MunI_RecognitionSites = repmat(20*(1:6)',1,6)+repmat(1:6,6,1);
            m.RM_MunI_MethylatedPositions = [3 4];
            m.RM_MunI_RestrictionPositions =[1 5];
            rxIdxs = [m.reactionIndexs_RM_MunI_Methylation; m.reactionIndexs_RM_MunI_Restriction];
            m.enzymeBounds(rxIdxs,2) = 1;
            
            m.enzymeBounds(m.reactionIndexs_polymerization,2) = 1e2;
            
            %% damage DNA                
            this.clearDNADamages();
            
            %BER
            c.damagedBases(310, 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'oxo8dAMP'}));
            
            c.abasicSites(325, 1) = 1;
            
            c.abasicSites(335, 1)  = 1;
            c.strandBreaks(335, 1) = 1;
            
            c.abasicSites(345, 1) = 1;
            c.strandBreaks(344, 1) = 1;
            
            %NER
            c.intrastrandCrossLinks(200, 2) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs({'dTMP64dCMP'}));
            
            %HR
            c.strandBreaks(150, 1:2) = 1;
            
            %restriction/modification
            c.damagedBases(m.RM_MunI_RecognitionSites(1, m.RM_MunI_MethylatedPositions(1)), 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %methylated
            c.damagedBases(m.RM_MunI_RecognitionSites(2, m.RM_MunI_MethylatedPositions(1)), 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %methylated
            c.damagedBases(m.RM_MunI_RecognitionSites(3, m.RM_MunI_MethylatedPositions(1)), 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %methylated
            c.damagedBases(m.RM_MunI_RecognitionSites(4, m.RM_MunI_MethylatedPositions(1)), 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %hemimethylated at first strand
            c.damagedBases(m.RM_MunI_RecognitionSites(5, m.RM_MunI_MethylatedPositions(1)), 1) = m.substrateMetaboliteGlobalIndexs(m.substrateIndexs_m6AD); %hemimethylated at second strand
            c.strandBreaks([m.RM_MunI_RecognitionSites(6, m.RM_MunI_RestrictionPositions)' [1;2]]) = 1; %cleaved
            
            %polymerization, ligation
            c.gapSites(450,1) = 1;
            c.abasicSites(450,1) = 1;
            c.strandBreaks([449;450],1) = 1;
            
            c.gapSites((451:462)',1) = 1;
            c.abasicSites((451:462)',1) = 1;
            c.strandBreaks([450 462],1) = 1;
            
            %% resources for repair
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.boundEnzymes(:) = 1e6;
            
            %% super class method
            this.helpTestGeneEssentiality({
                'MG_001'; %DNA polymerase III, beta subunit
                'MG_073'; %DNA incision complex
                %'MG_097'; %uracil-DNA glycosylase, putative
                %'MG_105'; %DNA integrity scanning protein
                'MG_184'; %adenine-specific DNA modification methylase
                'MG_190'; %phosphoesterase
                'MG_206'; %DNA incision complex
                'MG_235'; %apurinic endonuclease
                'MG_244'; %3'-5' helicase
                'MG_254'; %DNA ligase, NAD-dependent
                'MG_262'; %5'-3' exonuclease, putative
                'MG_339'; %recombination protein, strand exchange
                'MG_352'; %Holliday junction endonuclease
                'MG_358'; %Holliday junction DNA helicase
                'MG_359'; %Holliday junction DNA helicase
                'MG_421'; %DNA incision complex
                %'MG_438'; %type I restriction modification DNA specificity domain protein
                %'MG_498'; %formamidopyrimidine-DNA glycosylase
                }, @this.isProperlyFunctioning, ...
                struct('lengthSec', 5));
        end
    end
    
    %helper methods
    methods
        function result = isProperlyFunctioning(~, m, ~)
            result = false;
            
            c = m.chromosome;
            
            excludedPositions = reshape(m.RM_MunI_RecognitionSites, [], 1);
            
            if nnz(c.gapSites)               > 0 && nnz(c.gapSites)               > nnz(c.gapSites(              excludedPositions, :)); return; end;
            if nnz(c.abasicSites)            > 0 && nnz(c.abasicSites)            > nnz(c.abasicSites(           excludedPositions, :)); return; end;
            if nnz(c.damagedSugarPhosphates) > 0 && nnz(c.damagedSugarPhosphates) > nnz(c.damagedSugarPhosphates(excludedPositions, :)); return; end;
            if nnz(c.damagedBases)           > 0 && nnz(c.damagedBases)           > nnz(c.damagedBases(          excludedPositions, :)); return; end;
            if nnz(c.strandBreaks)           > 0 && nnz(c.strandBreaks)           > nnz(c.strandBreaks(          excludedPositions, :)); return; end;
            if nnz(c.intrastrandCrossLinks)  > 0 && nnz(c.intrastrandCrossLinks)  > nnz(c.intrastrandCrossLinks( excludedPositions, :)); return; end;
            if nnz(c.hollidayJunctions)      > 0 && nnz(c.hollidayJunctions)      > nnz(c.hollidayJunctions(     excludedPositions, :)); return; end;

            [~, ~, methylatedSites] = c.rmStatus(...
                m.RM_MunI_RecognitionSites, m.RM_MunI_MethylatedPositions, m.RM_MunI_RestrictionPositions, ...
                [], m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_RM_typeII));
            methylatedSites = find(methylatedSites);
            methylatedSites = methylatedSites(:,1);
            if isempty(methylatedSites); return; end;
            
            result = true;
        end
        
        function clearDNADamages(this)
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
        end
    end
end