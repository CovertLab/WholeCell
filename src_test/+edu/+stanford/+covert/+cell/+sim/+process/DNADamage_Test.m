%DNA damage process test case
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef DNADamage_Test < edu.stanford.covert.cell.sim.ReactionProcessTestCase
    %constants
    properties (Constant=true)
        expected_essentialGenes = {};
    end
    
    %constructor
    methods
        function this = DNADamage_Test(name)
            this = this@edu.stanford.covert.cell.sim.ReactionProcessTestCase(name);
        end
    end
    
    %simple test fixture
    methods
        function loadSimpleTestFixture(this)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            import edu.stanford.covert.cell.sim.state.Chromosome;
            
            %process
            m = this.process;
            
            %% chromosome state
            c = Chromosome([], []);
            c.monomerDNAFootprints = zeros(0, 1);
            c.complexDNAFootprints = zeros(0, 1);
            c.monomerDNAFootprintBindingStrandedness = repmat(c.dnaStrandedness_dsDNA, 0, 1);
            c.complexDNAFootprintBindingStrandedness = repmat(c.dnaStrandedness_dsDNA, 0, 1);
            c.monomerDNAFootprintRegionStrandedness = repmat(c.dnaStrandedness_dsDNA, 0, 1);
            c.complexDNAFootprintRegionStrandedness = repmat(c.dnaStrandedness_dsDNA, 0, 1);
            c.reactionBoundMonomer = zeros(0, 1);
            c.reactionBoundComplex = zeros(0, 1);
            c.reactionMonomerCatalysisMatrix = zeros(0, 0);
            c.reactionComplexCatalysisMatrix = zeros(0, 0);
            c.reactionThresholds = ...
                sum(c.reactionMonomerCatalysisMatrix, 2) + ...
                sum(c.reactionComplexCatalysisMatrix, 2);
            m.states = {c};
            m.chromosome = c;
            
            %parameters
            c.relaxedBasesPerTurn = 10.5;
            c.equilibriumSuperhelicalDensity = -0.06;
            c.supercoiledSuperhelicalDensityTolerance = 0.10;
            
            %genome
            seq = reshape([
                'AGCATTTTCAGGCACGGGGATACGGTTCATGGCGCTCCCGGAGACGACAA'
                'TCCAACCCGCTTGCGCCTTATAATCATTTCTCGAGCGAAAAACAAAAGCC'
                'ATAGTGTCTCTGCGCATCACTCGCACCAATGAATTCCGCGTCTGATAGCC'
                'TGCTTTGAGTCCCGGCTAATGTACAGTCCGATCACCCAACACCGGCAAGC'
                'GTTGCGGCACTAAGCTGGCACGATATACCGGTCGGGCGCCTGGCAGCACC'
                'AACTGTCGATAGGATCAGCGCTGTATCTAGAATGTGAGCTTGGCGGTCAG'
                'CTATTCTTGAATGAATATTCTTAGGTCGAAGCCTGCATTCAGCACCGCGG'
                'CGCTGCTAATTCTTTCGTTGGTTGGGCCAGCGAATGTGCGCCCTTCCGTA'
                'CTGGTCATGGTGCGGAGAAGCAACGTAATGACCGGAACATCCATTACGAC'
                'CCTAATCGAAGCTGACAGTTACTAGCGCTAGACGAGACGTACCAGGGAAG'],...
                1,[]);
            c.sequence = ChromosomeSequence(seq);
            c.sequenceLen = length(seq);
            c.sequenceGCContent = getGCContent(c.sequence);
            
            %% enzymes
            m.enzymeWholeCellModelIDs = cell(0,1);
            m.enzymeNames = cell(0,1);
            m.enzymeMolecularWeights = zeros(0,1);
            
            m.enzymes       = zeros(0, 1);
            m.boundEnzymes  = zeros(0, 1);
            
            %% substrates
            %whole cell model IDs, molecular weights
            m.substrateWholeCellModelIDs = {
                'AD'; 'CSN'; 'DA'; 'DAMP'; 'DC'; 'DCMP'; 'DG'; 'DGMP'; 'DHTHY'; 'DR5P';
                'DT'; 'DTMP'; 'ENU'; 'GLC'; 'GN'; 'H'; 'H2O'; 'MNU'; 'NH3'; 'NU'; 'O4eTHY';
                'O4mTHY'; 'PTRC'; 'THY'; 'THY64THY'; 'URA'; 'UVB_radiation'; 'butylene_glycol';
                'cyclodG'; 'dRibose5P_dRibose5P'; 'dhpURA'; 'gamma_radiation'; 'glchm5URA';
                'hm5URA'; 'ho5URA'; 'hydrogen_radical'; 'hydroxyl_radical'; 'm5CSN'; 'putTHY'};
            m.substrateMolecularWeights = [
                135.1265; 111.1018; 251.2414; 329.2055; 227.2167; 305.1808; 267.2408;
                345.2049; 128.1289; 212.0942; 242.2280; 320.1921; 117.1063; 180.1554;
                151.1259;   1.0079;  18.0152; 103.0798;  17.0304;  89.0533; 154.1661;
                140.1396;  90.1668; 126.1131; 252.2262; 112.0866;        0;  90.1206;
                265.2250; 407.1811; 214.2179;        0; 304.2527; 142.1125; 128.0860;
                1.0079;  17.0073; 125.1283; 212.2483];
            
            %names
            m.substrateNames = m.substrateWholeCellModelIDs;
            
            %indices
            m.substrateIndexs_radiation             = m.substrateIndexs({'UVB_radiation';'gamma_radiation'});
            m.substrateIndexs_hydrogen              = m.substrateIndexs({'H'});
            m.substrateIndexs_hydroxylRadical       = m.substrateIndexs({'hydroxyl_radical'});
            m.substrateIndexs_water                 = m.substrateIndexs({'H2O'});
            m.substrateIndexs_DR5P                  = m.substrateIndexs({'DR5P'});
            m.substrateIndexs_dRibose5P_dRibose5P   = m.substrateIndexs({'dRibose5P_dRibose5P'});
            m.substrateIndexs_unmodifiedNucleobases = m.substrateIndexs({'AD';'CSN';'GN';'THY'});
            m.substrateIndexs_unmodifiedNucleosides = m.substrateIndexs({'DA';'DC';'DG';'DT'});
            m.substrateIndexs_unmodifiedDNMPs       = m.substrateIndexs({'DAMP';'DCMP';'DGMP';'DTMP'});
            m.substrateIndexs_modifiedNucleobases   = m.substrateIndexs({
                'DHTHY'; 'O4eTHY'; 'O4mTHY'; 'THY64THY'; 'URA'; 'dhpURA';
                'glchm5URA'; 'hm5URA'; 'ho5URA'; 'm5CSN'; 'putTHY'});
            m.substrateIndexs_modifiedNucleosides   = m.substrateIndexs({'cyclodG'});
            m.substrateIndexs_modifiedDNMPs         = m.substrateIndexs({});
            
            m.substrateMetaboliteLocalIndexs = (1:numel(m.substrateWholeCellModelIDs))';
            m.substrateMetaboliteGlobalIndexs = m.substrateMetaboliteLocalIndexs;
            
            c.metabolite.dr5pIndexs = m.substrateIndexs({'DR5P'}) ;
            c.metabolite.waterIndexs = m.substrateIndexs({'H2O'});
            c.metabolite.dnmpIndexs = m.substrateIndexs({'AD'; 'CSN'; 'GN'; 'THY'});
            c.metabolite.m6ADIndexs = numel(m.substrateWholeCellModelIDs)+1;
            
            c.metabolite.molecularWeights = [m.substrateMolecularWeights; 0];
            
            %% reactions -- 1 reaction of each type
            %whole cell model IDs
            m.reactionWholeCellModelIDs = {
                'DNADamage_THYTHY_THY64THY_UVB_radiation';
                'DNADamage_hm5U_dhpURA_butylene_glycol';
                'DNADamage_hm5U_putTHY_PTRC';
                'DNADamage_THY_O4eTHY_ENU';
                'DNADamage_hm5U_glchm5URA_GLC';
                'DNADamage_THY_O4mTHY_MNU';
                'DNADamage_THY_DHTHY_hydrogen_radical';
                'DNADamage_DG_cyclodG_UVB_radiation';
                'DNADamage_URA_ho5URA_hydroxyl_radical';
                'DNADamage_SpontaneousBaseDeamination_m5C';
                'DNADamage_SpontaneousBaseLoss_thymine';
                'DNADamage_SpontaneousAbasicSiteStrandBreak'};
            
            %names
            m.reactionNames = m.reactionWholeCellModelIDs;
            
            %types
            m.reactionTypes = {
                'UV-B photodimerization';
                'base alkylation';
                'base amination';
                'base ethylation';
                'base glucosyl transfer';
                'base methylation';
                'base reduction';
                'photooxidation';
                'radiation (gamma-ray) induced base oxidation';
                'spontaneous base deamination';
                'spontaneous base loss';
                'strand break'};
            
            %indices
            m.reactionIndexs_spontaneousBaseLoss        = find(strcmp(m.reactionTypes,'spontaneous base loss'));
            m.reactionIndexs_spontaneousBaseDeamination = find(strcmp(m.reactionTypes,'spontaneous base deamination'));
            m.reactionIndexs_gammaRayBaseOxidation      = find(strcmp(m.reactionTypes,'radiation (gamma-ray) induced base oxidation'));
            m.reactionIndexs_uvbDimerization            = find(strcmp(m.reactionTypes,'UV-B photodimerization'));
            m.reactionIndexs_strandBreak                = find(strcmp(m.reactionTypes,'strand break'));
            m.reactionIndexs_baseReduction              = find(strcmp(m.reactionTypes,'base reduction'));
            m.reactionIndexs_baseAmination              = find(strcmp(m.reactionTypes,'base amination'));
            m.reactionIndexs_baseGlucosylTransfer       = find(strcmp(m.reactionTypes,'base glucosyl transfer'));
            m.reactionIndexs_baseAlkylation             = find(strcmp(m.reactionTypes,'base alkylation'));
            m.reactionIndexs_baseEthylation             = find(strcmp(m.reactionTypes,'base ethylation'));
            m.reactionIndexs_baseMethylation            = find(strcmp(m.reactionTypes,'base methylation'));
            m.reactionIndexs_photooxidation             = find(strcmp(m.reactionTypes,'photooxidation'));
            
            %ReactionProcess inherited reaction properties
            m.reactionCatalysisMatrix = zeros(length(m.reactionWholeCellModelIDs), 0);
            m.enzymeBounds = repmat([-Inf Inf], length(m.reactionWholeCellModelIDs), 1);
            
            %% reaction stoichiometry matrices, reaction sequences
            %allocate memory
            m.reactionStoichiometryMatrix              = zeros(length(m.substrateWholeCellModelIDs), length(m.reactionWholeCellModelIDs));
            m.reactionRadiationStoichiometryMatrix     = zeros(length(m.substrateWholeCellModelIDs), length(m.reactionWholeCellModelIDs));
            m.reactionSmallMoleculeStoichiometryMatrix = zeros(length(m.substrateWholeCellModelIDs), length(m.reactionWholeCellModelIDs));
            m.reactionDNAStoichiometryMatrix           = zeros(length(m.substrateWholeCellModelIDs), length(m.reactionWholeCellModelIDs));
            m.reactionVulnerableMotifs                 = cell(length(m.reactionWholeCellModelIDs), 1);
            m.reactionVulnerableMotifTypes             = cell(length(m.reactionWholeCellModelIDs), 1);
            m.reactionDamageTypes                      = cell(length(m.reactionWholeCellModelIDs), 1);
            m.reactionBounds                           = zeros(length(m.reactionWholeCellModelIDs), 2);
            
            %DNADamage_THYTHY_THY64THY_UVB_radiation
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'THY'}), m.reactionIndexs_uvbDimerization) = -2;
            m.reactionRadiationStoichiometryMatrix(m.substrateIndexs({'UVB_radiation'}), m.reactionIndexs_uvbDimerization) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'THY64THY'}), m.reactionIndexs_uvbDimerization) = 1;
            m.reactionVulnerableMotifs{m.reactionIndexs_uvbDimerization} = 'TT';
            m.reactionVulnerableMotifTypes{m.reactionIndexs_uvbDimerization} = 'sequence';
            m.reactionDamageTypes{m.reactionIndexs_uvbDimerization} = 'intrastrandCrossLinks';
            m.reactionBounds(m.reactionIndexs_uvbDimerization, 2) = 1.8300e-008;
            
            %DNADamage_hm5U_dhpURA_butylene_glycol
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'butylene_glycol'}), m.reactionIndexs_baseAlkylation) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'hm5URA'}), m.reactionIndexs_baseAlkylation) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dhpURA'}), m.reactionIndexs_baseAlkylation) = 1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O'}), m.reactionIndexs_baseAlkylation) = 1;
            m.reactionVulnerableMotifs{m.reactionIndexs_baseAlkylation}  = m.substrateIndexs({'hm5URA'});
            m.reactionVulnerableMotifTypes{m.reactionIndexs_baseAlkylation} = 'damagedBases';
            m.reactionDamageTypes{m.reactionIndexs_baseAlkylation} = 'damagedBases';
            m.reactionBounds(m.reactionIndexs_baseAlkylation, 2) = 0;
            
            %DNADamage_hm5U_putTHY_PTRC
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'hm5URA'}), m.reactionIndexs_baseAmination) = -1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'PTRC'}), m.reactionIndexs_baseAmination) = -1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H'}), m.reactionIndexs_baseAmination) = 2;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O'}), m.reactionIndexs_baseAmination) = 1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'putTHY'}), m.reactionIndexs_baseAmination) = 1;
            m.reactionVulnerableMotifs{m.reactionIndexs_baseAmination}  = m.substrateIndexs({'hm5URA'});
            m.reactionVulnerableMotifTypes{m.reactionIndexs_baseAmination} = 'damagedBases';
            m.reactionDamageTypes{m.reactionIndexs_baseAmination} = 'damagedBases';
            m.reactionBounds(m.reactionIndexs_baseAmination, 2) = 0;
            
            %DNADamage_THY_O4eTHY_ENU
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'ENU'}), m.reactionIndexs_baseEthylation) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'THY'}), m.reactionIndexs_baseEthylation) = -1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'NU'}), m.reactionIndexs_baseEthylation) = 1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'O4eTHY'}), m.reactionIndexs_baseEthylation) = 1;
            m.reactionVulnerableMotifs{m.reactionIndexs_baseEthylation}  = 'T';
            m.reactionVulnerableMotifTypes{m.reactionIndexs_baseEthylation} = 'sequence';
            m.reactionDamageTypes{m.reactionIndexs_baseEthylation} = 'damagedBases';
            m.reactionBounds(m.reactionIndexs_baseEthylation, 2) = 0;
            
            %DNADamage_hm5U_glchm5URA_GLC
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'GLC'}), m.reactionIndexs_baseGlucosylTransfer) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'hm5URA'}), m.reactionIndexs_baseGlucosylTransfer) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'glchm5URA'}), m.reactionIndexs_baseGlucosylTransfer) = 1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O'}), m.reactionIndexs_baseGlucosylTransfer) = 1;
            m.reactionVulnerableMotifs{m.reactionIndexs_baseGlucosylTransfer}  = m.substrateIndexs({'hm5URA'});
            m.reactionVulnerableMotifTypes{m.reactionIndexs_baseGlucosylTransfer} = 'damagedBases';
            m.reactionDamageTypes{m.reactionIndexs_baseGlucosylTransfer} = 'damagedBases';
            m.reactionBounds(m.reactionIndexs_baseGlucosylTransfer, 2) = 0;
            
            %DNADamage_THY_O4mTHY_MNU
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'MNU'}), m.reactionIndexs_baseMethylation) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'THY'}), m.reactionIndexs_baseMethylation) = -1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'NU'}), m.reactionIndexs_baseMethylation) = 1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'O4mTHY'}), m.reactionIndexs_baseMethylation) = 1;
            m.reactionVulnerableMotifs{m.reactionIndexs_baseMethylation}  = 'T';
            m.reactionVulnerableMotifTypes{m.reactionIndexs_baseMethylation} = 'sequence';
            m.reactionDamageTypes{m.reactionIndexs_baseMethylation} = 'damagedBases';
            m.reactionBounds(m.reactionIndexs_baseMethylation, 2) = 0;
            
            %DNADamage_THY_DHTHY_hydrogen_radical
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'hydrogen_radical'}), m.reactionIndexs_baseReduction) = -2;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'THY'}), m.reactionIndexs_baseReduction) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'DHTHY'}), m.reactionIndexs_baseReduction) = 1;
            m.reactionVulnerableMotifs{m.reactionIndexs_baseReduction}  = 'T';
            m.reactionVulnerableMotifTypes{m.reactionIndexs_baseReduction} = 'sequence';
            m.reactionDamageTypes{m.reactionIndexs_baseReduction} = 'damagedBases';
            m.reactionBounds(m.reactionIndexs_baseReduction, 2) = 0;
            
            %DNADamage_DG_cyclodG_UVB_radiation
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'DG'}), m.reactionIndexs_photooxidation) = -1;
            m.reactionRadiationStoichiometryMatrix(m.substrateIndexs({'UVB_radiation'}), m.reactionIndexs_photooxidation) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'cyclodG'}), m.reactionIndexs_photooxidation) = 1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H'}), m.reactionIndexs_photooxidation) = 2;
            m.reactionVulnerableMotifs{m.reactionIndexs_photooxidation}  = 'G';
            m.reactionVulnerableMotifTypes{m.reactionIndexs_photooxidation} = 'sequence';
            m.reactionDamageTypes{m.reactionIndexs_photooxidation} = 'damagedBases';
            m.reactionBounds(m.reactionIndexs_photooxidation, 2) = 0;
            
            %DNADamage_URA_ho5URA_hydroxyl_radical
            m.reactionRadiationStoichiometryMatrix(m.substrateIndexs({'gamma_radiation'}), m.reactionIndexs_gammaRayBaseOxidation) = -1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O'}), m.reactionIndexs_gammaRayBaseOxidation) = -1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H'}), m.reactionIndexs_gammaRayBaseOxidation) = 2;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'URA'}), m.reactionIndexs_gammaRayBaseOxidation) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'ho5URA'}), m.reactionIndexs_gammaRayBaseOxidation) = 1;
            m.reactionVulnerableMotifs{m.reactionIndexs_gammaRayBaseOxidation}  = m.substrateIndexs({'URA'});
            m.reactionVulnerableMotifTypes{m.reactionIndexs_gammaRayBaseOxidation} = 'damagedBases';
            m.reactionDamageTypes{m.reactionIndexs_gammaRayBaseOxidation} = 'damagedBases';
            m.reactionBounds(m.reactionIndexs_gammaRayBaseOxidation, 2) = 5.8999e-008;
            
            %DNADamage_SpontaneousBaseDeamination_m5C
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O'}), m.reactionIndexs_spontaneousBaseDeamination) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'m5CSN'}), m.reactionIndexs_spontaneousBaseDeamination) = -1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'NH3'}), m.reactionIndexs_spontaneousBaseDeamination) = 1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'THY'}), m.reactionIndexs_spontaneousBaseDeamination) = 1;
            m.reactionVulnerableMotifs{m.reactionIndexs_spontaneousBaseDeamination}  = m.substrateIndexs({'m5CSN'});
            m.reactionVulnerableMotifTypes{m.reactionIndexs_spontaneousBaseDeamination} = 'damagedBases';
            m.reactionDamageTypes{m.reactionIndexs_spontaneousBaseDeamination} = 'damagedBases';
            m.reactionBounds(m.reactionIndexs_spontaneousBaseDeamination, 2) = 2.2000e-014;
            
            %DNADamage_SpontaneousBaseLoss_thymine
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'DTMP'}), m.reactionIndexs_spontaneousBaseLoss) = -1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O'}), m.reactionIndexs_spontaneousBaseLoss) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'DR5P'}), m.reactionIndexs_spontaneousBaseLoss) = 1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'THY'}), m.reactionIndexs_spontaneousBaseLoss) = 1;
            m.reactionVulnerableMotifs{m.reactionIndexs_spontaneousBaseLoss}  = 'T';
            m.reactionVulnerableMotifTypes{m.reactionIndexs_spontaneousBaseLoss} = 'sequence';
            m.reactionDamageTypes{m.reactionIndexs_spontaneousBaseLoss} = 'abasicSites';
            m.reactionBounds(m.reactionIndexs_spontaneousBaseLoss, 2) = 6.0000e-012;
            
            %DNADamage_SpontaneousAbasicSiteStrandBreak
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'dRibose5P_dRibose5P'}), m.reactionIndexs_strandBreak) = -1;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H2O'}), m.reactionIndexs_strandBreak) = -1;
            m.reactionDNAStoichiometryMatrix(m.substrateIndexs({'DR5P'}), m.reactionIndexs_strandBreak) = 2;
            m.reactionSmallMoleculeStoichiometryMatrix(m.substrateIndexs({'H'}), m.reactionIndexs_strandBreak) = 1;
            m.reactionVulnerableMotifs{m.reactionIndexs_strandBreak}  = true;
            m.reactionVulnerableMotifTypes{m.reactionIndexs_strandBreak}  = 'abasicSites';
            m.reactionDamageTypes{m.reactionIndexs_strandBreak} = 'strandBreaks';
            m.reactionBounds(m.reactionIndexs_strandBreak, 2) = 1.9254e-006;
            
            %reaction DNA product
            [i, j] = find(m.reactionDNAStoichiometryMatrix > 0);
            m.reactionDNAProduct = zeros(numel(m.reactionWholeCellModelIDs), 1);
            m.reactionDNAProduct(j) = m.substrateMetaboliteGlobalIndexs(ismembc2(i, m.substrateMetaboliteLocalIndexs));
            
            %reaction radiation
            [i, j] = find(m.reactionRadiationStoichiometryMatrix);
            m.reactionRadiation = zeros(numel(m.reactionWholeCellModelIDs), 1);
            m.reactionRadiation(j) = i;
            
            %reactionStoichiometryMatrix
            m.reactionStoichiometryMatrix = ...
                m.reactionRadiationStoichiometryMatrix + ...
                m.reactionSmallMoleculeStoichiometryMatrix + ...
                m.reactionDNAStoichiometryMatrix;            
            
            %% initial state
            m.substrates = zeros(length(m.substrateWholeCellModelIDs), 1);
            
            m.chromosome.initialize();
        end
    end
    
    %tests
    methods
        function testConstants(this)
            m = this.process;
            c = m.chromosome;
            
            %assert that are sorted
            assertTrue(issorted(m.substrateMetaboliteLocalIndexs));
            
            %assert that this process doesn't create or modify m6AD            
            assertFalse(any(m.substrateMetaboliteGlobalIndexs == c.metabolite.m6ADIndexs));
            for i = 1:numel(m.reactionVulnerableMotifs)
                assertFalse(isequal(m.reactionVulnerableMotifs{i}, c.metabolite.m6ADIndexs));
            end
            
            %check that reaction properties were loaded correctly from KB
            substrateWholeCellModelIDs               = m.substrateWholeCellModelIDs;            
            substrateMetaboliteGlobalIndexs          = m.substrateMetaboliteGlobalIndexs;
            substrateMetaboliteLocalIndexs           = m.substrateMetaboliteLocalIndexs;
            reactionWholeCellModelIDs                = m.reactionWholeCellModelIDs;
            reactionTypes                            = m.reactionTypes;
            reactionDamageTypes                      = m.reactionDamageTypes;
            reactionVulnerableMotifs                 = m.reactionVulnerableMotifs;
            reactionVulnerableMotifTypes             = m.reactionVulnerableMotifTypes;
            reactionSmallMoleculeStoichiometryMatrix = m.reactionSmallMoleculeStoichiometryMatrix;
            reactionDNAStoichiometryMatrix           = m.reactionDNAStoichiometryMatrix;
            reactionRadiationStoichiometryMatrix     = m.reactionRadiationStoichiometryMatrix;
            reactionBounds                           = m.reactionBounds;
            
            this.loadSimpleTestFixture();
            
            [~, sbIdxs1, sbIdxs2] = intersect(substrateWholeCellModelIDs, m.substrateWholeCellModelIDs);
            [~, rxIdxs1, rxIdxs2] = intersect(reactionWholeCellModelIDs, m.reactionWholeCellModelIDs);
            
            for i=1:length(rxIdxs1)
                motif = reactionVulnerableMotifs{rxIdxs1(i)};
                mMotif = m.reactionVulnerableMotifs{rxIdxs2(i)};
                if isnumeric(motif) || isnumeric(mMotif)
                    motif = substrateWholeCellModelIDs{substrateMetaboliteLocalIndexs(motif == substrateMetaboliteGlobalIndexs)};
                    mMotif = m.substrateWholeCellModelIDs{m.substrateMetaboliteLocalIndexs(mMotif == m.substrateMetaboliteGlobalIndexs)};
                end
                
                assertEqual(reactionTypes{rxIdxs1(i)},                                    m.reactionTypes{rxIdxs2(i)},                                    sprintf('Error in definition of reaction %s', m.reactionWholeCellModelIDs{rxIdxs2(i)}));
                assertEqual(reactionDamageTypes{rxIdxs1(i)},                              m.reactionDamageTypes{rxIdxs2(i)},                              sprintf('Error in definition of reaction %s', m.reactionWholeCellModelIDs{rxIdxs2(i)}));
                assertEqual(motif,                                                        mMotif,                                                         sprintf('Error in definition of reaction %s', m.reactionWholeCellModelIDs{rxIdxs2(i)}));
                assertEqual(reactionVulnerableMotifTypes{rxIdxs1(i)},                     m.reactionVulnerableMotifTypes{rxIdxs2(i)},                     sprintf('Error in definition of reaction %s', m.reactionWholeCellModelIDs{rxIdxs2(i)}));
                assertEqual(reactionSmallMoleculeStoichiometryMatrix(sbIdxs1,rxIdxs1(i)), m.reactionSmallMoleculeStoichiometryMatrix(sbIdxs2,rxIdxs2(i)), sprintf('Error in definition of reaction %s', m.reactionWholeCellModelIDs{rxIdxs2(i)}));
                assertEqual(reactionDNAStoichiometryMatrix(sbIdxs1,rxIdxs1(i)),           m.reactionDNAStoichiometryMatrix(sbIdxs2,rxIdxs2(i)),           sprintf('Error in definition of reaction %s', m.reactionWholeCellModelIDs{rxIdxs2(i)}));
                assertEqual(reactionRadiationStoichiometryMatrix(sbIdxs1,rxIdxs1(i)),     m.reactionRadiationStoichiometryMatrix(sbIdxs2,rxIdxs2(i)),     sprintf('Error in definition of reaction %s', m.reactionWholeCellModelIDs{rxIdxs2(i)}));
                assertEqual(reactionBounds(rxIdxs1(i),1),                                 m.reactionBounds(rxIdxs2(i),1),                                 sprintf('Error in definition of reaction %s', m.reactionWholeCellModelIDs{rxIdxs2(i)}));
                assertElementsAlmostEqual(reactionBounds(rxIdxs1(i),2), m.reactionBounds(rxIdxs2(i),2), 'relative','1e-4', sprintf('Error in definition of reaction %s',m.reactionWholeCellModelIDs{rxIdxs2(i)}));
            end           
        end
        
        function testReactionStoichiometryMatrixFactorization(this)
            m = this.process;
            
            dMW1 = m.reactionStoichiometryMatrix'*m.substrateMolecularWeights;
            dMW2 = (...
                m.reactionSmallMoleculeStoichiometryMatrix+...
                m.reactionDNAStoichiometryMatrix+...
                m.reactionRadiationStoichiometryMatrix)'*m.substrateMolecularWeights;
            
            assertElementsAlmostEqual(dMW1, dMW2);
        end
        
        function testBaseAlkylation(this)
            this.loadSimpleTestFixture();
            
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_hm5U_dhpURA_butylene_glycol'});
            
            %set amount of substrates
            m.substrates = ...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 1;
            
            %set DNA to undamaged
            this.clearDNADamages();
            m.chromosome.damagedBases(1,1) = m.substrateMetaboliteGlobalIndexs(m.substrateMetaboliteLocalIndexs == find(m.reactionDNAStoichiometryMatrix(:,rxIdx)<0));
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.damagedBases);
            [~, idxs] = ismember(vals, m.substrateMetaboliteGlobalIndexs);
            vals = m.substrateMetaboliteLocalIndexs(idxs, :);
            
            %correct sequence damaged
            assertEqual([1 1], subs);
            
            %correct damage recorded
            assertEqual(find(m.reactionDNAStoichiometryMatrix(:,rxIdx)>0), unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.damagedBases, m.chromosome.getDamagedSites(true, true, true, true, false, false, true));
            
            %substrates updated correctly
            nReactions = size(subs,1);
            assertEqual(m.substrates, ...
                nReactions * max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)));
        end
        
        function testBaseAmination(this)
            this.loadSimpleTestFixture();
            
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_hm5U_putTHY_PTRC'});
            
            %set amount of substrates
            m.substrates = ...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 1;
            
            %set DNA to undamaged
            this.clearDNADamages();
            m.chromosome.damagedBases(1,1) = m.substrateMetaboliteGlobalIndexs(m.substrateMetaboliteLocalIndexs == find(m.reactionDNAStoichiometryMatrix(:,rxIdx)<0));
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.damagedBases);
            [~, idxs] = ismember(vals, m.substrateMetaboliteGlobalIndexs);
            vals = m.substrateMetaboliteLocalIndexs(idxs, :);
            
            %correct sequence damaged
            assertEqual([1 1], subs);
            
            %correct damage recorded
            assertEqual(find(m.reactionDNAStoichiometryMatrix(:,rxIdx)>0), unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.damagedBases, m.chromosome.getDamagedSites(true, true, true, true, false, false, true));
            
            %substrates updated correctly
            nReactions = size(subs,1);
            assertEqual(m.substrates, ...
                nReactions * max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)));
        end
        
        function testBaseEthylation(this)
            this.loadSimpleTestFixture();
            
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_THY_O4eTHY_ENU'});
            
            %set amount of substrates
            m.substrates = ...
                max(0, -m.reactionSmallMoleculeStoichiometryMatrix(:, rxIdx));
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 10 / size(m.chromosome.sequence,1);
            
            %set DNA to undamaged
            this.clearDNADamages();
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.damagedBases);
            [~, idxs] = ismember(vals, m.substrateMetaboliteGlobalIndexs);
            vals = m.substrateMetaboliteLocalIndexs(idxs, :);
            
            %correct sequence damaged
            assertTrue(~isempty(subs));
            assertEqual(m.reactionVulnerableMotifs{rxIdx}, this.reactedSequences(subs, rxIdx));
            
            %correct damage recorded
            assertEqual(find(m.reactionDNAStoichiometryMatrix(:,rxIdx)>0), unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.damagedBases, m.chromosome.getDamagedSites(true, true, true, true, false, false, true));
            
            %substrates updated correctly
            nReactions = size(subs,1);
            assertEqual(m.substrates, ...
                nReactions * max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)));
        end
        
        function testBaseGlucosylTransfer(this)
            this.loadSimpleTestFixture();
            
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_hm5U_glchm5URA_GLC'});
            
            %set amount of substrates
            m.substrates = ...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 1;
            
            %set DNA to undamaged
            this.clearDNADamages();
            m.chromosome.damagedBases(1,1) = m.substrateMetaboliteGlobalIndexs(m.substrateMetaboliteLocalIndexs == find(m.reactionDNAStoichiometryMatrix(:,rxIdx)<0));
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.damagedBases);
            [~, idxs] = ismember(vals, m.substrateMetaboliteGlobalIndexs);
            vals = m.substrateMetaboliteLocalIndexs(idxs, :);
            
            %correct sequence damaged
            assertEqual([1 1], subs);
            
            %correct damage recorded
            assertEqual(find(m.reactionDNAStoichiometryMatrix(:,rxIdx)>0), unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.damagedBases, m.chromosome.getDamagedSites(true, true, true, true, false, false, true));
            
            %substrates updated correctly
            nReactions = size(subs,1);
            assertEqual(m.substrates, ...
                nReactions * max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)));
        end
        
        function testBaseMethylation(this)
            this.loadSimpleTestFixture();
            
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_THY_O4mTHY_MNU'});
            
            %set amount of substrates
            m.substrates = ...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 10 / size(m.chromosome.sequence,1);
            
            %set DNA to undamaged
            this.clearDNADamages();
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.damagedBases);
            [~, idxs] = ismember(vals, m.substrateMetaboliteGlobalIndexs);
            vals = m.substrateMetaboliteLocalIndexs(idxs, :);
            
            %correct sequence damaged
            assertEqual(m.reactionVulnerableMotifs{rxIdx}, this.reactedSequences(subs, rxIdx));
            
            %correct damage recorded
            assertEqual(find(m.reactionDNAStoichiometryMatrix(:,rxIdx)>0), unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.damagedBases, m.chromosome.getDamagedSites(true, true, true, true, false, false, true));
            
            %substrates updated correctly
            nReactions = size(subs,1);
            assertEqual(m.substrates, ...
                nReactions * max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)));
        end
        
        function testBaseReduction(this)
            this.loadSimpleTestFixture();
            
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_THY_DHTHY_hydrogen_radical'});
            
            %set amount of substrates
            m.substrates = ...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 10 / size(m.chromosome.sequence,1);
            
            %set DNA to undamaged
            this.clearDNADamages();
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.damagedBases);
            [~, idxs] = ismember(vals, m.substrateMetaboliteGlobalIndexs);
            vals = m.substrateMetaboliteLocalIndexs(idxs, :);
            
            %correct sequence damaged
            assertEqual(m.reactionVulnerableMotifs{rxIdx}, this.reactedSequences(subs, rxIdx));
            
            %correct damage recorded
            assertEqual(find(m.reactionDNAStoichiometryMatrix(:,rxIdx)>0), unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.damagedBases, m.chromosome.getDamagedSites(true, true, true, true, false, false, true));
            
            %substrates updated correctly
            nReactions = size(subs,1);
            assertEqual(m.substrates, ...
                nReactions * max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)));
        end
        
        function testPhotooxidation(this)
            this.loadSimpleTestFixture();
            
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_DG_cyclodG_UVB_radiation'});
            
            %set amount of substrates
            m.substrates = ...
                max(0,-m.reactionRadiationStoichiometryMatrix(:,rxIdx)) + ...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 10 / size(m.chromosome.sequence,1);
            
            %set DNA to undamaged
            this.clearDNADamages();
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.damagedBases);
            [~, idxs] = ismember(vals, m.substrateMetaboliteGlobalIndexs);
            vals = m.substrateMetaboliteLocalIndexs(idxs, :);
            
            %correct sequence damaged
            assertEqual(m.reactionVulnerableMotifs{rxIdx}, this.reactedSequences(subs, rxIdx));
            
            %correct damage recorded
            assertEqual(find(m.reactionDNAStoichiometryMatrix(:,rxIdx)>0), unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.damagedBases, m.chromosome.getDamagedSites(true, true, true, true, false, false, true));
            
            %substrates updated correctly
            nReactions = size(subs,1);
            assertEqual(m.substrates, ...
                max(0,-m.reactionRadiationStoichiometryMatrix(:,rxIdx)) + ...
                nReactions * max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)));
        end
        
        function testGammaRayBaseOxidation(this)
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_URA_ho5URA_hydroxyl_radical'});
            
            %set amount of substrates
            m.substrates = ...
                max(0,-m.reactionRadiationStoichiometryMatrix(:,rxIdx)) + ...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 1;
            
            %set DNA to undamaged
            this.clearDNADamages();
            m.chromosome.damagedBases(1,1) = m.substrateMetaboliteGlobalIndexs(m.substrateMetaboliteLocalIndexs == find(m.reactionDNAStoichiometryMatrix(:,rxIdx)<0));
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.damagedBases);
            [~, idxs] = ismember(vals, m.substrateMetaboliteGlobalIndexs);
            vals = m.substrateMetaboliteLocalIndexs(idxs, :);
            
            %correct sequence damagedthis.reactedSequences(subs, rxIdx)
            assertEqual([1 1], subs);
            
            %correct damage recorded
            assertEqual(find(m.reactionDNAStoichiometryMatrix(:,rxIdx)>0), unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.damagedBases, m.chromosome.getDamagedSites(true, true, true, true, false, false, true));
            
            %substrates updated correctly
            nReactions = size(subs,1);
            assertEqual(m.substrates, ...
                max(0,-m.reactionRadiationStoichiometryMatrix(:,rxIdx)) + ...
                nReactions * max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)));
        end
        
        function testSpontaneousBaseDeamination(this)
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_SpontaneousBaseDeamination_m5C'});
            
            %set amount of substrates
            m.substrates = ...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 1;
            
            %set DNA to undamaged
            this.clearDNADamages();
            m.chromosome.damagedBases(1,1) = m.substrateMetaboliteGlobalIndexs(m.substrateMetaboliteLocalIndexs == find(m.reactionDNAStoichiometryMatrix(:,rxIdx)<0));
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.damagedBases);
            [~, idxs] = ismember(vals, m.substrateMetaboliteGlobalIndexs);
            vals = m.substrateMetaboliteLocalIndexs(idxs, :);
            
            %correct sequence damaged
            assertEqual([1 1], subs);
            
            %correct damage recorded
            assertEqual(find(m.reactionDNAStoichiometryMatrix(:,rxIdx)>0), unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.damagedBases, m.chromosome.getDamagedSites(true, true, true, true, false, false, true));
            
            %substrates updated correctly
            nReactions = size(subs,1);
            assertEqual(m.substrates, ...
                nReactions * max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)));
        end
        
        function testSpontaneousBaseLoss(this)
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_SpontaneousBaseLoss_thymine'});
            
            %set amount of substrates
            m.substrates = ...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 10 / size(m.chromosome.sequence,1);
            
            %set DNA to undamaged
            this.clearDNADamages();
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.abasicSites);
            
            %correct sequence damaged
            assertEqual(m.reactionVulnerableMotifs{rxIdx}, this.reactedSequences(subs, rxIdx));
            
            %correct damage recorded
            assertEqual(true, unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.abasicSites, m.chromosome.getDamagedSites(true, true, true, true, false, false, true));
            
            %substrates updated correctly
            nReactions = size(subs,1);
            assertEqual(m.substrates, ...
                nReactions * max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)));
        end
        
        function testStrandBreak(this)
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_SpontaneousAbasicSiteStrandBreak'});
            
            %set amount of substrates
            m.substrates = ...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 1;
            
            %set DNA to undamaged
            this.clearDNADamages();
            m.chromosome.abasicSites(1, 1) = true;
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.strandBreaks);
            
            %correct sequence damaged
            assertEqual([1 1], subs);
            
            %correct damage recorded
            assertEqual(true, unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.strandBreaks | m.chromosome.strandBreaks5, m.chromosome.getDamagedSites(true, true, true, true, false, false, true) ~= 0);
            
            %substrates updated correctly
            nReactions = size(subs,1);
            assertEqual(m.substrates, ...
                nReactions * max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)));
        end
        
        function testUVBDimerization(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_THYTHY_THY64THY_UVB_radiation'});
            
            %set amount of substrates
            m.substrates = ...
                max(0,-m.reactionRadiationStoichiometryMatrix(:,rxIdx)) + ...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 10 / size(m.chromosome.sequence,1);
            
            %set DNA to undamaged
            this.clearDNADamages();
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.intrastrandCrossLinks);
            [~, idxs] = ismember(vals, m.substrateMetaboliteGlobalIndexs);
            vals = m.substrateMetaboliteLocalIndexs(idxs, :);
            
            %correct sequence damaged
            assertEqual(m.reactionVulnerableMotifs{rxIdx}, this.reactedSequences(subs, rxIdx));
            
            %correct damage recorded
            assertEqual(find(m.reactionDNAStoichiometryMatrix(:,rxIdx)>0), unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.intrastrandCrossLinks + m.chromosome.intrastrandCrossLinks5, m.chromosome.getDamagedSites(true, true, true, true, false, false, true));
            
            %substrates updated correctly
            nReactions = size(subs,1);
            assertEqual(m.substrates, ...
                max(0,-m.reactionRadiationStoichiometryMatrix(:,rxIdx)) + ...
                nReactions * max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)));
        end
        
        function testMultipleReactions(this)
            this.loadSimpleTestFixture();
            
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_hm5U_dhpURA_butylene_glycol'; 'DNADamage_THY_O4eTHY_ENU'});
            
            %set amount of substrates
            m.substrates = sum(...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)),...
                2);
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = [1; 10 / size(m.chromosome.sequence,1)];
            
            %set DNA to undamaged
            this.clearDNADamages();
            m.chromosome.damagedBases(1,1) = m.substrateMetaboliteGlobalIndexs(m.substrateMetaboliteLocalIndexs == find(m.reactionDNAStoichiometryMatrix(:,rxIdx(1))<0));
            
            %evolve state
            m.evolveState();
            
            %find locations of damages
            [subs, vals] = find(m.chromosome.damagedBases);
            [~, idxs] = ismember(vals, m.substrateMetaboliteGlobalIndexs);
            vals = m.substrateMetaboliteLocalIndexs(idxs, :);
            
            %correct sequence damaged
            assertEqual([1 1], subs(vals == find(m.reactionDNAStoichiometryMatrix(:,rxIdx(1))>0),:));
            assertEqual(m.reactionVulnerableMotifs{rxIdx(2)}, this.reactedSequences(subs(vals == find(m.reactionDNAStoichiometryMatrix(:,rxIdx(2))>0),:), rxIdx(2)));
            
            %correct damage recorded
            [~, vals] = find(m.chromosome.damagedBases);
            [~, idxs] = ismember(vals, m.substrateMetaboliteGlobalIndexs);
            vals = m.substrateMetaboliteLocalIndexs(idxs, :);
            assertEqual(find(any(m.reactionDNAStoichiometryMatrix(:,rxIdx)>0,2)), unique(vals));
            
            %no other damages occured
            assertEqual(m.chromosome.damagedBases, m.chromosome.getDamagedSites(true, true, true, true, false, false, true));
            
            %substrates updated correctly
            nReactions = sum(m.reactionDNAStoichiometryMatrix(vals,rxIdx)>0,1)';
            assertEqual(m.substrates, ...
                max(0,m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx)) * nReactions);
        end
        
        function testFairness(this)
            this.loadSimpleTestFixture();
            
            m = this.process;
            rxIdx = m.reactionIndexs({'DNADamage_hm5U_dhpURA_butylene_glycol'});
            dnaSubstrate = find(m.reactionDNAStoichiometryMatrix(:,rxIdx)<0);
            dnaProduct = find(m.reactionDNAStoichiometryMatrix(:,rxIdx)>0);
            
            dnaSubstrate = m.substrateMetaboliteGlobalIndexs(m.substrateMetaboliteLocalIndexs == dnaSubstrate);
            dnaProduct = m.substrateMetaboliteGlobalIndexs(m.substrateMetaboliteLocalIndexs == dnaProduct);
            
            %set amount of substrates
            m.substrates = ...
                max(0,-m.reactionSmallMoleculeStoichiometryMatrix(:,rxIdx));
            initial_substrates = m.substrates;
            
            %increase rate of damage
            m.reactionBounds(:,2) = 0;
            m.reactionBounds(rxIdx, 2) = 1;
            
            %set DNA to undamaged
            this.clearDNADamages();
            m.chromosome.damagedBases(1:2,1) = dnaSubstrate;
            damagedBases = m.chromosome.damagedBases;
            
            counts = zeros(2,1);
            for i=1:200
                %reset state
                m.substrates = initial_substrates;
                m.chromosome.damagedBases = damagedBases;
                
                %evolve state
                m.evolveState();
                
                %find locations of damages
                [subs, vals] = find(m.chromosome.damagedBases);
                subs = subs(vals == dnaProduct,1);
                
                %update counts
                counts(subs) = counts(subs) + 1;
            end
            
            %check positions modified fairly
            assertTrue(max(counts)-min(counts) < 0.2 * max(counts));
        end
        
        function testGeneEssentiality(this)
            m = this.process;
            c = m.chromosome;
            c.polymerizedRegions(1, 1:2) = size(c.sequence, 1);
            
            m.substrates(:) = 1e3;
            m.substrates(m.substrateIndexs_radiation) = 1;
            
            %increase rates of damage
            m.reactionBounds(:,2) = 10 / size(m.chromosome.sequence,1);
            
            damagedSites = m.chromosome.getDamagedSites(true, true, true, true, false, false, true);
            this.helpTestGeneEssentiality(...
                {},...
                @(m,i) nnz(m.chromosome.getDamagedSites(true, true, true, true, false, false, true)) > nnz(damagedSites));
        end
    end
    
    %helper methods
    methods
        function clearDNADamages(this)
            m = this.process;
            
            m.chromosome.initialize();
        end
        
        function sequence = reactedSequences(this, subs, rxIdx)
            m = this.process;
            positions = mod(...
                repmat(subs(:,1), 1, length(m.reactionVulnerableMotifs{rxIdx})) + ...
                repmat(2*mod(subs(:,2),2)-1, 1, length(m.reactionVulnerableMotifs{rxIdx})) .* ...
                repmat(0:length(m.reactionVulnerableMotifs{rxIdx})-1, size(subs,1), 1) ...
                -1, size(m.chromosome.sequence,1))+1;
            strands = repmat(subs(:,2), 1, length(m.reactionVulnerableMotifs{rxIdx}));
            
            sequence = reshape(unique(m.chromosome.sequence.subsequence(positions, strands),'rows'), size(m.reactionVulnerableMotifs{rxIdx}));
        end
    end
end
