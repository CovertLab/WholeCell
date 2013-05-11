%tRNA aminoacylation process test case
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef tRNAAminoacylation_Test < edu.stanford.covert.cell.sim.ReactionProcessTestCase
    methods
        function this = tRNAAminoacylation_Test(name)
            this = this@edu.stanford.covert.cell.sim.ReactionProcessTestCase(name);
        end
    end

    %fixtures
    methods       
        function loadSimpleTestFixture(this)
            %process
            m = this.process;

            %whole cell model IDs
            m.compartmentWholeCellModelIDs={'c';'e';'m';'tc';'tm'};
            m.reactionWholeCellModelIDs = {'MG475_Aminoacylation';
                'MG485_Aminoacylation'; 'MG488_Aminoacylation';
                'MG488_Formyltransferase'; 'MG502_Amidotransferase';
                'MG502_Aminoacylation';'MG_0004_Aminoacylation'};
            m.substrateWholeCellModelIDs = {'ADP'; 'ALA'; 'AMP'; 'ARG'; 'ASN';
                'ASP'; 'ATP'; 'CYS'; 'FMET'; 'FTHF10'; 'GLU'; 'GLY'; 'H'; 'H2O';
                'HIS'; 'ILE'; 'LEU'; 'LYS'; 'MET'; 'PHE'; 'PI'; 'PPI'; 'PRO';
                'SER'; 'THF'; 'THR'; 'TRP'; 'TYR'; 'VAL'; 'GLN'};
            m.enzymeWholeCellModelIDs = {'MG_253_MONOMER'; 'MG_266_MONOMER'; 'MG_334_MONOMER';
                'MG_345_MONOMER'; 'MG_365_MONOMER'; 'MG_378_MONOMER'; 'MG_462_MONOMER';
                'MG_005_DIMER'; 'MG_021_DIMER'; 'MG_035_DIMER'; 'MG_036_DIMER';
                'MG_098_099_100_TRIMER'; 'MG_113_DIMER'; 'MG_126_DIMER'; 'MG_136_DIMER';
                'MG_194_195_TETRAMER'; 'MG_251_DIMER'; 'MG_283_DIMER'; 'MG_292_TETRAMER';
                'MG_375_DIMER'; 'MG_455_DIMER'; 'MG475'; 'MG485'; 'MG488'; 'MG502';'MG_0004'};
            m.freeRNAWholeCellModelIDs = {'MG475';'MG485';'MG488';'MG502';'MG_0004'};
            m.aminoacylatedRNAWholeCellModelIDs = m.freeRNAWholeCellModelIDs;

            %names
            m.reactionNames  = m.reactionWholeCellModelIDs;
            m.substrateNames = m.substrateWholeCellModelIDs;
            m.enzymeNames    = m.enzymeWholeCellModelIDs;

            %types
            m.reactionTypes = {'aminoacylation';'aminoacylation';'aminoacylation';
                'transfer';'transfer';'aminoacylation';'aminoacylation'};

            %indices
            m.compartmentIndexs_cytosol = find(strcmp(m.compartmentWholeCellModelIDs, 'c'));

            m.substrateIndexs_aminoAcids  = m.substrateIndexs({'ALA'; 'ARG'; 'ASN'; 'ASP'; 'CYS'; 'GLN'; 'GLU'; 'GLY'; 'HIS'; 'ILE'; 'LEU'; 'LYS'; 'MET'; 'PHE'; 'PRO'; 'SER'; 'THR'; 'TRP'; 'TYR'; 'VAL'});
            m.substrateIndexs_glutamate   = m.substrateIndexs({'GLU'});
            m.substrateIndexs_glutamine   = m.substrateIndexs({'GLN'});
            m.substrateIndexs_methionine  = m.substrateIndexs({'MET'});
            m.substrateIndexs_fmethionine = m.substrateIndexs({'FMET'});
            m.substrateIndexs_atp         = m.substrateIndexs({'ATP'});
            m.substrateIndexs_amp         = m.substrateIndexs({'AMP'});
            m.substrateIndexs_adp         = m.substrateIndexs({'ADP'});
            m.substrateIndexs_diphosphate = m.substrateIndexs({'PPI'});
            m.substrateIndexs_phosphate   = m.substrateIndexs({'PI'});
            m.substrateIndexs_water       = m.substrateIndexs({'H2O'});
            m.substrateIndexs_hydrogen    = m.substrateIndexs({'H'});
            m.substrateIndexs_fthf10      = m.substrateIndexs({'FTHF10'});
            m.substrateIndexs_thf         = m.substrateIndexs({'THF'});

            m.enzymeIndexs_tRNASynthetases = m.enzymeIndexs({...
                'MG_292_TETRAMER';...       %alanyl-tRNA synthetase
                'MG_378_MONOMER';...        %arginyl-tRNA synthetase
                'MG_036_DIMER';...          %aspartyl-tRNA synthetase
                'MG_113_DIMER';...          %asparaginyl-tRNA synthetase
                'MG_253_MONOMER';...        %cysteinyl-tRNA synthetase
                'MG_462_MONOMER';...        %glutamyl-tRNA synthetase
                'MG_251_DIMER';...          %glycyl-tRNA synthetase
                'MG_035_DIMER';...          %histidyl-tRNA synthetase
                'MG_345_MONOMER';...        %isoleucyl-tRNA synthetase
                'MG_266_MONOMER';...        %leucyl-tRNA synthetase
                'MG_136_DIMER';...          %lysyl-tRNA synthetase
                'MG_021_DIMER';...          %methionyl-tRNA synthetase
                'MG_194_195_TETRAMER';...   %phenylalanyl-tRNA synthetase
                'MG_283_DIMER';...          %prolyl-tRNA synthetase
                'MG_005_DIMER';...          %seryl-tRNA synthetase
                'MG_375_DIMER';...          %threonyl-tRNA synthetase
                'MG_126_DIMER';...          %tryptophanyl-tRNA synthetase
                'MG_455_DIMER';...          %tyrosyl-tRNA synthetase
                'MG_334_MONOMER'});         %valyl-tRNA synthetase
            m.enzymeIndexs_tRNATransferases = m.enzymeIndexs({...
                'MG_098_099_100_TRIMER';... %glutamyl-tRNA(Gln) and/or aspartyl-tRNA(Asn) amidotransferase
                'MG_365_MONOMER'});         %methionyl-tRNA formyltransferase
            m.enzymeIndexs_tRNAGlutamylSynthetase         = m.enzymeIndexs({'MG_462_MONOMER'});        %glutamyl-tRNA synthetase
            m.enzymeIndexs_tRNAMethionylSynthetase        = m.enzymeIndexs({'MG_021_DIMER'});          %methionyl-tRNA synthetase
            m.enzymeIndexs_tRNAGlutamylAmidotransferase   = m.enzymeIndexs({'MG_098_099_100_TRIMER'}); %glutamyl-tRNA(Gln) and/or aspartyl-tRNA(Asn) amidotransferase
            m.enzymeIndexs_tRNAMethionylFormyltransferase = m.enzymeIndexs({'MG_365_MONOMER'});        %methionyl-tRNA formyltransferase
            m.enzymeIndexs_tRNAs = m.enzymeIndexs(m.freeRNAWholeCellModelIDs);
            m.enzymeIndexs_tmRNA = m.enzymeIndexs({'MG_0004'});

            m.reactionIndexs_aminoacylation          = find(strcmp(m.reactionTypes,'aminoacylation'));
            m.reactionIndexs_transfer                = find(strcmp(m.reactionTypes,'transfer'));
            m.reactionIndexs_glutamylamidotransfer   = m.reactionIndexs({'MG502_Amidotransferase'});
            m.reactionIndexs_methionylformyltransfer = m.reactionIndexs({'MG488_Formyltransferase'});

            %reactions
            %MG475_Aminoacylation: ATP + SER ==> AMP + PPI
            %MG485_Aminoacylation: ATP + MET ==> AMP + PPI
            %MG488_Aminoacylation: ATP + MET ==> AMP + PPI
            %MG488_Formyltransferase: FTHF10 + MET ==> FMET + THF
            %MG502_Amidotransferase: ATP + GLN + GLU + H2O ==> ADP + GLN + GLU + H + PI
            %MG502_Aminoacylation: ATP + GLU ==> AMP + PPI
            %MG_0004_Aminoacylation: ALA + ATP ==> AMP + PPI
            m.reactionStoichiometryMatrix = zeros(length(m.substrateWholeCellModelIDs), length(m.reactionWholeCellModelIDs));

            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'ATP'), strcmp(m.reactionWholeCellModelIDs, 'MG475_Aminoacylation'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'SER'), strcmp(m.reactionWholeCellModelIDs, 'MG475_Aminoacylation'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'AMP'), strcmp(m.reactionWholeCellModelIDs, 'MG475_Aminoacylation'))=1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'PPI'), strcmp(m.reactionWholeCellModelIDs, 'MG475_Aminoacylation'))=1;

            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'ATP'), strcmp(m.reactionWholeCellModelIDs, 'MG485_Aminoacylation'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'MET'), strcmp(m.reactionWholeCellModelIDs, 'MG485_Aminoacylation'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'AMP'), strcmp(m.reactionWholeCellModelIDs, 'MG485_Aminoacylation'))=1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'PPI'), strcmp(m.reactionWholeCellModelIDs, 'MG485_Aminoacylation'))=1;

            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'ATP'), strcmp(m.reactionWholeCellModelIDs, 'MG488_Aminoacylation'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'MET'), strcmp(m.reactionWholeCellModelIDs, 'MG488_Aminoacylation'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'AMP'), strcmp(m.reactionWholeCellModelIDs, 'MG488_Aminoacylation'))=1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'PPI'), strcmp(m.reactionWholeCellModelIDs, 'MG488_Aminoacylation'))=1;

            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'FTHF10'), strcmp(m.reactionWholeCellModelIDs, 'MG488_Formyltransferase'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'THF'), strcmp(m.reactionWholeCellModelIDs, 'MG488_Formyltransferase'))=1;

            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'ATP'), strcmp(m.reactionWholeCellModelIDs, 'MG502_Amidotransferase'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'GLN'), strcmp(m.reactionWholeCellModelIDs, 'MG502_Amidotransferase'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'H2O'), strcmp(m.reactionWholeCellModelIDs, 'MG502_Amidotransferase'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'GLU'), strcmp(m.reactionWholeCellModelIDs, 'MG502_Amidotransferase'))=1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'ADP'), strcmp(m.reactionWholeCellModelIDs, 'MG502_Amidotransferase'))=1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'H'), strcmp(m.reactionWholeCellModelIDs, 'MG502_Amidotransferase'))=1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'PI'), strcmp(m.reactionWholeCellModelIDs, 'MG502_Amidotransferase'))=1;

            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'ATP'), strcmp(m.reactionWholeCellModelIDs, 'MG502_Aminoacylation'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'GLU'), strcmp(m.reactionWholeCellModelIDs, 'MG502_Aminoacylation'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'AMP'), strcmp(m.reactionWholeCellModelIDs, 'MG502_Aminoacylation'))=1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'PPI'), strcmp(m.reactionWholeCellModelIDs, 'MG502_Aminoacylation'))=1;

            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'ATP'), strcmp(m.reactionWholeCellModelIDs, 'MG_0004_Aminoacylation'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'ALA'), strcmp(m.reactionWholeCellModelIDs, 'MG_0004_Aminoacylation'))=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'AMP'), strcmp(m.reactionWholeCellModelIDs, 'MG_0004_Aminoacylation'))=1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'PPI'), strcmp(m.reactionWholeCellModelIDs, 'MG_0004_Aminoacylation'))=1;

            m.reactionModificationMatrix = zeros(length(m.reactionWholeCellModelIDs), length(m.freeRNAWholeCellModelIDs));
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG475_Aminoacylation'), strcmp(m.freeRNAWholeCellModelIDs,    'MG475'))   = 1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG485_Aminoacylation'), strcmp(m.freeRNAWholeCellModelIDs,    'MG485'))   = 1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG488_Aminoacylation'), strcmp(m.freeRNAWholeCellModelIDs,    'MG488'))   = 1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG488_Formyltransferase'), strcmp(m.freeRNAWholeCellModelIDs, 'MG488'))   = 1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG502_Amidotransferase'), strcmp(m.freeRNAWholeCellModelIDs,  'MG502'))   = 1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG502_Aminoacylation'), strcmp(m.freeRNAWholeCellModelIDs,    'MG502'))   = 1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG_0004_Aminoacylation'), strcmp(m.freeRNAWholeCellModelIDs,  'MG_0004')) = 1;

            m.reactionCatalysisMatrix = zeros(length(m.reactionWholeCellModelIDs), length(m.enzymeWholeCellModelIDs));
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG475_Aminoacylation'),    strcmp(m.enzymeWholeCellModelIDs, 'MG_005_DIMER'))=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG485_Aminoacylation'),    strcmp(m.enzymeWholeCellModelIDs, 'MG_021_DIMER'))=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG488_Aminoacylation'),    strcmp(m.enzymeWholeCellModelIDs, 'MG_021_DIMER'))=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG488_Formyltransferase'), strcmp(m.enzymeWholeCellModelIDs, 'MG_365_MONOMER'))=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG502_Amidotransferase'),  strcmp(m.enzymeWholeCellModelIDs, 'MG_098_099_100_TRIMER'))=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG502_Aminoacylation'),    strcmp(m.enzymeWholeCellModelIDs, 'MG_462_MONOMER'))=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs, 'MG_0004_Aminoacylation'),  strcmp(m.enzymeWholeCellModelIDs, 'MG_292_TETRAMER'))=1;

            m.enzymeBounds = repmat([-Inf Inf], length(m.reactionWholeCellModelIDs), 1);
            m.enzymeBounds(strcmp(m.reactionWholeCellModelIDs, 'MG475_Aminoacylation'),    2) = 7.6;
            m.enzymeBounds(strcmp(m.reactionWholeCellModelIDs, 'MG485_Aminoacylation'),    2) = 4.2627;
            m.enzymeBounds(strcmp(m.reactionWholeCellModelIDs, 'MG488_Aminoacylation'),    2) = 7.6;
            m.enzymeBounds(strcmp(m.reactionWholeCellModelIDs, 'MG488_Formyltransferase'), 2) = 90;
            m.enzymeBounds(strcmp(m.reactionWholeCellModelIDs, 'MG502_Amidotransferase'),  2) = 141.9;
            m.enzymeBounds(strcmp(m.reactionWholeCellModelIDs, 'MG502_Aminoacylation'),    2) = 14;
            m.enzymeBounds(strcmp(m.reactionWholeCellModelIDs, 'MG_0004_Aminoacylation'),  2) = 3.2;
            
            m.initializeSpeciesNetwork();

            %molecular weights
            m.substrateMolecularWeights = [
                424.1769;  89.0929; 345.2049; 175.2083; 132.1176; 132.0945; 503.1489;
                121.1579; 177.2210; 471.4226; 146.1210;  75.0664;   1.0079;  18.0152;
                155.1542; 131.1724; 131.1724; 147.1949; 149.2109; 165.1887;  95.9793;
                174.9513; 115.1301; 105.0923; 443.4125; 119.1188; 204.2247; 181.1881;
                117.1459; 146.1441];
            m.enzymeMolecularWeights = ones(size(m.enzymeWholeCellModelIDs));
            m.freeRNAMolecularWeights = 1e5 * (1:length(m.freeRNAWholeCellModelIDs))';
            m.aminoacylatedRNAMolecularWeights = m.freeRNAMolecularWeights - ...
                m.reactionModificationMatrix'*m.reactionStoichiometryMatrix'*m.substrateMolecularWeights;

            %initial state            
            m.substrates        = zeros(length(m.substrateWholeCellModelIDs),        1);
            m.enzymes           = zeros(length(m.enzymeWholeCellModelIDs),           1);
            m.boundEnzymes      = zeros(length(m.enzymeWholeCellModelIDs),           1);
            m.freeRNAs          = zeros(length(m.freeRNAWholeCellModelIDs),          1);
            m.aminoacylatedRNAs = zeros(length(m.aminoacylatedRNAWholeCellModelIDs), 1);
        end
    end

    %tests
    methods
        function testOneAminoacylation(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.substrates(m.substrateIndexs({'SER'})) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_005_DIMER'})) = 1e6;
            m.freeRNAs(:) = 0;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG475')) = 1;
            m.aminoacylatedRNAs(:) = 0;

            m.evolveState();
            assertEqual(0,   m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG475')));
            assertEqual(1,   m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG475')));
            assertEqual(0,   m.substrates(m.substrateIndexs_atp));
            assertEqual(0,   m.substrates(m.substrateIndexs({'SER'})));
            assertEqual(1,   m.substrates(m.substrateIndexs_amp));
            assertEqual(1,   m.substrates(m.substrateIndexs_diphosphate));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_005_DIMER'})));
        end

        function testOneAminoacylation_noEnyzme(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.substrates(m.substrateIndexs({'SER'})) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_005_DIMER'})) = 0;
            m.freeRNAs(:) = 0;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG475')) = 1;
            m.aminoacylatedRNAs(:) = 0;

            m.evolveState();
            assertEqual(1, m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG475')));
            assertEqual(0, m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG475')));
            assertEqual(1, m.substrates(m.substrateIndexs_atp));
            assertEqual(1, m.substrates(m.substrateIndexs({'SER'})));
            assertEqual(0, m.substrates(m.substrateIndexs_amp));
            assertEqual(0, m.substrates(m.substrateIndexs_diphosphate));
            assertEqual(0, m.enzymes(m.enzymeIndexs({'MG_005_DIMER'})));
        end

        function testOneAminoacylation_noATP(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 0;
            m.substrates(m.substrateIndexs({'SER'})) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_005_DIMER'})) = 1e6;
            m.freeRNAs(:) = 0;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG475')) = 1;
            m.aminoacylatedRNAs(:) = 0;

            m.evolveState();
            assertEqual(1,   m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG475')));
            assertEqual(0,   m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG475')));
            assertEqual(0,   m.substrates(m.substrateIndexs_atp));
            assertEqual(1,   m.substrates(m.substrateIndexs({'SER'})));
            assertEqual(0,   m.substrates(m.substrateIndexs_amp));
            assertEqual(0,   m.substrates(m.substrateIndexs_diphosphate));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_005_DIMER'})));
        end

        function testTRNARequiringFormylTransfer(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.substrates(m.substrateIndexs_methionine) = 1;
            m.substrates(m.substrateIndexs_fthf10)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_021_DIMER'})) = 1e6;
            m.enzymes(m.enzymeIndexs({'MG_365_MONOMER'})) = 1e6;
            m.freeRNAs(:) = 0;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG488')) = 1;
            m.aminoacylatedRNAs(:) = 0;

            m.evolveState();
            assertEqual(0,   m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG488')));
            assertEqual(1,   m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG488')));
            assertEqual(0,   m.substrates(m.substrateIndexs_atp));
            assertEqual(0,   m.substrates(m.substrateIndexs_methionine));
            assertEqual(0,   m.substrates(m.substrateIndexs_fthf10));
            assertEqual(1,   m.substrates(m.substrateIndexs_amp));
            assertEqual(1,   m.substrates(m.substrateIndexs_diphosphate));
            assertEqual(1,   m.substrates(m.substrateIndexs_thf));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_021_DIMER'})));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_365_MONOMER'})));
        end

        function testTRNARequiringFormylTransfer_insufficientEnzyme(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.substrates(m.substrateIndexs_methionine) = 1;
            m.substrates(m.substrateIndexs_fthf10)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_021_DIMER'})) = 0;
            m.enzymes(m.enzymeIndexs({'MG_365_MONOMER'})) = 1e6;
            m.freeRNAs(:) = 0;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG488')) = 1;
            m.aminoacylatedRNAs(:) = 0;

            m.evolveState();
            assertEqual(1,   m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG488')));
            assertEqual(0,   m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG488')));
            assertEqual(1,   m.substrates(m.substrateIndexs_atp));
            assertEqual(1,   m.substrates(m.substrateIndexs_methionine));
            assertEqual(1,   m.substrates(m.substrateIndexs_fthf10));
            assertEqual(0,   m.substrates(m.substrateIndexs_amp));
            assertEqual(0,   m.substrates(m.substrateIndexs_diphosphate));
            assertEqual(0,   m.substrates(m.substrateIndexs_thf));
            assertEqual(0, m.enzymes(m.enzymeIndexs({'MG_021_DIMER'})));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_365_MONOMER'})));
        end

        function testTRNARequiringFormylTransfer_oneTooFewATP(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 0;
            m.substrates(m.substrateIndexs_methionine) = 1;
            m.substrates(m.substrateIndexs_fthf10)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_021_DIMER'})) = 1e6;
            m.enzymes(m.enzymeIndexs({'MG_365_MONOMER'})) = 1e6;
            m.freeRNAs(:) = 0;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG488')) = 1;
            m.aminoacylatedRNAs(:) = 0;

            m.evolveState();
            assertEqual(1,   m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG488')));
            assertEqual(0,   m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG488')));
            assertEqual(0,   m.substrates(m.substrateIndexs_atp));
            assertEqual(1,   m.substrates(m.substrateIndexs_methionine));
            assertEqual(1,   m.substrates(m.substrateIndexs_fthf10));
            assertEqual(0,   m.substrates(m.substrateIndexs_amp));
            assertEqual(0,   m.substrates(m.substrateIndexs_diphosphate));
            assertEqual(0,   m.substrates(m.substrateIndexs_thf));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_021_DIMER'})));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_365_MONOMER'})));
        end

        function testTRNARequiringAmidoTransfer(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_glutamate) = 1;
            m.substrates(m.substrateIndexs_glutamine) = 1;
            m.substrates(m.substrateIndexs_water)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_098_099_100_TRIMER'})) = 1e6;
            m.enzymes(m.enzymeIndexs({'MG_462_MONOMER'})) = 1e6;
            m.freeRNAs(:) = 0;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG502')) = 1;
            m.aminoacylatedRNAs(:) = 0;

            m.evolveState();
            assertEqual(0,   m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG502')));
            assertEqual(1,   m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG502')));
            assertEqual(0,   m.substrates(m.substrateIndexs_atp));
            assertEqual(1,   m.substrates(m.substrateIndexs_glutamate));
            assertEqual(0,   m.substrates(m.substrateIndexs_glutamine));
            assertEqual(0,   m.substrates(m.substrateIndexs_water));
            assertEqual(1,   m.substrates(m.substrateIndexs_adp));
            assertEqual(1,   m.substrates(m.substrateIndexs_diphosphate));
            assertEqual(1,   m.substrates(m.substrateIndexs_phosphate));
            assertEqual(1,   m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_098_099_100_TRIMER'})));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_462_MONOMER'})));
        end

        function testTRNARequiringAmidoTransfer_oneTooATP(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.substrates(m.substrateIndexs_glutamate) = 1;
            m.substrates(m.substrateIndexs_glutamine) = 1;
            m.substrates(m.substrateIndexs_water) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_098_099_100_TRIMER'})) = 1e6;
            m.enzymes(m.enzymeIndexs({'MG_462_MONOMER'})) = 1e6;
            m.freeRNAs(:) = 0;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs, 'MG502')) = 1;
            m.aminoacylatedRNAs(:) = 0;

            m.evolveState();
            
            assertEqual(1,   m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs, 'MG502')));
            assertEqual(0,   m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs, 'MG502')));
            assertEqual(1,   m.substrates(m.substrateIndexs_atp));
            assertEqual(1,   m.substrates(m.substrateIndexs_glutamate));
            assertEqual(1,   m.substrates(m.substrateIndexs_glutamine));
            assertEqual(1,   m.substrates(m.substrateIndexs_water));
            assertEqual(0,   m.substrates(m.substrateIndexs_adp));
            assertEqual(0,   m.substrates(m.substrateIndexs_diphosphate));
            assertEqual(0,   m.substrates(m.substrateIndexs_phosphate));
            assertEqual(0,   m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_098_099_100_TRIMER'})));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_462_MONOMER'})));
        end

        function testTRNARequiringAmidoTransfer_oneTooFewGlutamate(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_glutamate) = 0;
            m.substrates(m.substrateIndexs_glutamine) = 1;
            m.substrates(m.substrateIndexs_water)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_098_099_100_TRIMER'})) = 1e6;
            m.enzymes(m.enzymeIndexs({'MG_462_MONOMER'})) = 1e6;
            m.freeRNAs(:) = 0;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG502')) = 1;
            m.aminoacylatedRNAs(:) = 0;

            m.evolveState();
            assertEqual(1,   m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG502')));
            assertEqual(0,   m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG502')));
            assertEqual(2,   m.substrates(m.substrateIndexs_atp));
            assertEqual(0,   m.substrates(m.substrateIndexs_glutamate));
            assertEqual(1,   m.substrates(m.substrateIndexs_glutamine));
            assertEqual(1,   m.substrates(m.substrateIndexs_water));
            assertEqual(0,   m.substrates(m.substrateIndexs_adp));
            assertEqual(0,   m.substrates(m.substrateIndexs_diphosphate));
            assertEqual(0,   m.substrates(m.substrateIndexs_phosphate));
            assertEqual(0,   m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_098_099_100_TRIMER'})));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_462_MONOMER'})));
        end

        function testTRNARequiringAmidoTransfer_insufficientEnzyme(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 2;
            m.substrates(m.substrateIndexs_glutamate) = 1;
            m.substrates(m.substrateIndexs_glutamine) = 1;
            m.substrates(m.substrateIndexs_water)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_098_099_100_TRIMER'})) = 1e6;
            m.enzymes(m.enzymeIndexs({'MG_462_MONOMER'})) = 0;
            m.freeRNAs(:) = 0;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG502')) = 1;
            m.aminoacylatedRNAs(:) = 0;

            m.evolveState();
            assertEqual(1,   m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG502')));
            assertEqual(0,   m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG502')));
            assertEqual(2,   m.substrates(m.substrateIndexs_atp));
            assertEqual(1,   m.substrates(m.substrateIndexs_glutamate));
            assertEqual(1,   m.substrates(m.substrateIndexs_glutamine));
            assertEqual(1,   m.substrates(m.substrateIndexs_water));
            assertEqual(0,   m.substrates(m.substrateIndexs_adp));
            assertEqual(0,   m.substrates(m.substrateIndexs_diphosphate));
            assertEqual(0,   m.substrates(m.substrateIndexs_phosphate));
            assertEqual(0,   m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs({'MG_098_099_100_TRIMER'})));
            assertEqual(0, m.enzymes(m.enzymeIndexs({'MG_462_MONOMER'})));
        end

        function testOneOfEveryRNA(this)
            m = this.process;

            m.substrates = sum(max(0,-m.reactionStoichiometryMatrix),2);
            m.enzymes(:) = 1e6;
            m.freeRNAs(:) = 1;
            m.aminoacylatedRNAs(:) = 0;

            m.evolveState();
            assertEqual(zeros(size(m.freeRNAs)), m.freeRNAs);
            assertEqual(ones(size(m.aminoacylatedRNAs)), m.aminoacylatedRNAs);
            assertEqual(sum(max(0,m.reactionStoichiometryMatrix),2), m.substrates);
            assertEqual(repmat(1e6, size(m.enzymes)), m.enzymes);
        end

        % Start with a lot of each of two kinds of RNAs requiring
        % aminoacylation and transfer, and insufficient ATP to modify all
        % of both kinds, then verify that we modify them "fairly"
        function testFairAllocationOfResources(this)
            numATP = 1000;
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = numATP;
            m.substrates(m.substrateIndexs({'SER'})) = numATP;
            m.substrates(m.substrateIndexs({'MET'})) = numATP;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs({'MG_005_DIMER'})) = 1e6;
            m.enzymes(m.enzymeIndexs({'MG_021_DIMER'})) = 1e6;
            m.freeRNAs(:) = 0;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG475')) = numATP;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG485')) = numATP;
            m.aminoacylatedRNAs(:) = 0;

            m.evolveState();
            aminoacylatedRNAs = [m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG475'));m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs,'MG485'))];
            assertEqual(0, m.substrates(m.substrateIndexs_atp));
            assertEqual(numATP, sum(aminoacylatedRNAs));
            assertTrue(...
                abs(diff(aminoacylatedRNAs)) < 0.15 * max(aminoacylatedRNAs), ...
                sprintf('significant monomer imbalance: %d, %d', ...
                aminoacylatedRNAs(1), aminoacylatedRNAs(2)));
        end

        function testGeneEssentiality(this)
            m = this.process;
            m.substrates = sum(max(0,-m.reactionStoichiometryMatrix),2);
            m.enzymes(:) = 1e6;
            m.freeRNAs(:) = 1;
            m.aminoacylatedRNAs(:) = 0;

            this.helpTestGeneEssentiality({
                'MG_005'; %seryl-tRNA synthetase
                'MG_021'; %methionyl-tRNA synthetase
                'MG_035'; %histidyl-tRNA synthetase
                'MG_036'; %aspartyl-tRNA synthetase
                'MG_098'; %glutamyl-tRNA(Gln) and/or aspartyl-tRNA(Asn) amidotransferase, C subunit
                'MG_099'; %glutamyl-tRNA(Gln) and/or aspartyl-tRNA(Asn) amidotransferase, A subunit
                'MG_100'; %glutamyl-tRNA(Gln) and/or aspartyl-tRNA(Asn) amidotransferase, B subunit
                'MG_113'; %asparaginyl-tRNA synthetase
                'MG_126'; %tryptophanyl-tRNA synthetase
                'MG_136'; %lysyl-tRNA synthetase
                'MG_194'; %phenylalanyl-tRNA synthetase, alpha subunit
                'MG_195'; %phenylalanyl-tRNA synthetase, beta subunit
                'MG_251'; %glycyl-tRNA synthetase
                'MG_253'; %cysteinyl-tRNA synthetase
                'MG_266'; %leucyl-tRNA synthetase
                'MG_283'; %prolyl-tRNA synthetase
                'MG_292'; %alanyl-tRNA synthetase
                'MG_334'; %valyl-tRNA synthetase
                'MG_345'; %isoleucyl-tRNA synthetase
                'MG_365'; %methionyl-tRNA formyltransferase
                'MG_375'; %threonyl-tRNA synthetase
                'MG_378'; %arginyl-tRNA synthetase
                'MG_455'; %tyrosyl-tRNA synthetase
                'MG_462'; %glutamyl-tRNA synthetase
                'MG471';  %GCA;GCC;GCG;GCT
                'MG472';  %ATC;ATT
                'MG475';  %AGC;AGT
                'MG479';  %ACC;ACT
                'MG483';  %TGC;TGT
                'MG484';  %CCA;CCC;CCG;CCT
                'MG485';  %ATG
                'MG486';  %ATA
                'MG487';  %TCA
                'MG488';  %ATA;ATC;ATG;ATT;TTA;TTG;GTG
                'MG489';  %GAC;GAT
                'MG490';  %TTC;TTT
                'MG492';  %CGC;CGT
                'MG493';  %GGC;GGT
                'MG495';  %AGG
                'MG496';  %TGG
                'MG497';  %CGA;CGG
                'MG499';  %GGA;GGG
                'MG500';  %TTA
                'MG501';  %AAA
                'MG502';  %CAA;CAG
                'MG503';  %TAC;TAT
                'MG504';  %TGA
                'MG506';  %TCC;TCT
                'MG507';  %TCG
                'MG508';  %CTA;CTG
                'MG509';  %AAG
                'MG510';  %ACA
                'MG511';  %GTA;GTC;GTG;GTT
                'MG512';  %ACG
                'MG513';  %GAA;GAG
                'MG514';  %AAC;AAT
                'MG518';  %CAC;CAT
                'MG519';  %CTC;CTT
                'MG520';  %TTG
                'MG523';  %AGA
                }, @this.isProperlyFunctioning);
        end
    end

    %helper methods
    methods
        function result = isProperlyFunctioning(~, m, i)
            indexs_requiredTRNAs = setdiff(1:length(m.freeRNAWholeCellModelIDs), find(strcmp(m.freeRNAWholeCellModelIDs, 'MG_0004')));
            result = all(i.aminoacylatedRNAs(indexs_requiredTRNAs) < m.aminoacylatedRNAs(indexs_requiredTRNAs));
        end

        %knockout tRNAs
        function helpKnockoutGene(this, geneWholeCellModelID)
            m = this.process;
            m.freeRNAs(strcmp(m.freeRNAWholeCellModelIDs,geneWholeCellModelID))=0;
            m.aminoacylatedRNAs(strcmp(m.freeRNAWholeCellModelIDs,geneWholeCellModelID))=0;
        end
    end
end
