%Protein decay test cases
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/04/2010
classdef ProteinDecay_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    %constructor
    methods
        function this = ProteinDecay_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end
    
    %simple test fixtures
    methods        
        function loadSimpleFixture(this)
            m = this.process;
            
            %% constants
            m.lonProteaseSpecificRate         = 1.667;  %lon protease kinetic rate [PUB_0029]
            m.lonProteaseEnergyCost    	      = 6;      %lon protease energy requirement [PUB_0029]
            m.lonProteaseFragmentLength       = 20;     %length of peptide fragments after cleavage by lon protease [PUB_0029]
            m.ftsHProteaseSpecificRate        = 0.03;   %ftsH protease kinetic rate [PUB_0031]
            m.ftsHProteaseEnergyCost          = 8.3;    %ftsH protease energy requirement [PUB_0031]
            m.ftsHProteaseFragmentLength      = 15;     %length of peptide fragments after cleavage by ftsH protease [PUB_0030]
            m.oligoendopeptidaseFSpecificRate = 27.474; %oligoendopeptidase F kinetic rate [PUB_0035]
            m.proteinMisfoldingRate           = 1e-6;   %frequency of protein misfolds
            
            %% compartments
            m.compartment.wholeCellModelIDs = {'c'; 'e'; 'm'; 'tc'; 'tm'};
            m.compartment.cytosolIndexs  = find(strcmp(m.compartment.wholeCellModelIDs, 'c'));
            m.compartment.membraneIndexs = find(strcmp(m.compartment.wholeCellModelIDs, 'm'));
            m.compartment.terminalOrganelleCytosolIndexs  = find(strcmp(m.compartment.wholeCellModelIDs, 'tc'));
            m.compartment.terminalOrganelleMembraneIndexs = find(strcmp(m.compartment.wholeCellModelIDs, 'tm'));
            
            %% metabolites
            m.substrateWholeCellModelIDs = {
                'ADP'; 'ALA'; 'ARG'; 'ASN'; 'ASP'; 'ATP'; 'COA'; 'CYS'; 'FMET'; 'FOR';
                'GLN'; 'GLU'; 'GLY'; 'H'; 'H2O'; 'HIS'; 'ILE'; 'LEU'; 'LIPOYLLYS'; 'LYS';
                'MET'; 'MG'; 'NH3'; 'PAP'; 'PHE'; 'PI'; 'PRO'; 'SER'; 'THR'; 'TRP'; 'TYR';
                'VAL'; 'ZN'; 'diacylglycerolCys'; 'pSER'; 'pTHR'; 'pTYR'};
            
            m.substrateNames = m.substrateWholeCellModelIDs;
            
            m.substrateMolecularWeights = [
                424.1769;  89.0929; 175.2083; 132.1176; 132.0945; 503.1489; 763.5012; 121.1579;
                177.2210;  45.0174; 146.1441; 146.1210;  75.0664;   1.0079;  18.0152; 155.1542;
                131.1724; 131.1724; 334.4968; 147.1949; 149.2109;   24.305;  17.0304; 423.1690;
                165.1887;  95.9793; 115.1301; 105.0923; 119.1188; 204.2247; 181.1881; 117.1459;
                65.39; 672.0514; 183.0564; 197.0829; 259.1522];
            
            m.substrateIndexs_atp         = m.substrateIndexs({'ATP'});
            m.substrateIndexs_adp         = m.substrateIndexs({'ADP'});
            m.substrateIndexs_phosphate   = m.substrateIndexs({'PI'});
            m.substrateIndexs_hydrogen    = m.substrateIndexs({'H'});
            m.substrateIndexs_water       = m.substrateIndexs({'H2O'});
            m.substrateIndexs_aminoAcids  = m.substrateIndexs({'ALA'; 'ARG'; 'ASN'; 'ASP'; 'CYS'; 'GLN'; 'GLU'; 'GLY'; 'HIS'; 'ILE'; 'LEU'; 'LYS'; 'MET'; 'PHE'; 'PRO'; 'SER'; 'THR'; 'TRP'; 'TYR'; 'VAL'; 'FMET'});
            m.substrateIndexs_methionine  = m.substrateIndexs({'MET'});
            m.substrateIndexs_fmethionine = m.substrateIndexs({'FMET'});
            m.substrateIndexs_glutamate   = m.substrateIndexs({'GLU'});
            m.substrateIndexs_glutamine   = m.substrateIndexs({'GLN'});
            m.substrateIndexs_ammonia     = m.substrateIndexs({'NH3'});
            m.substrateIndexs_formate     = m.substrateIndexs({'FOR'});
            
            %% monomers
            m.monomer.wholeCellModelIDs = repmat({...
                'MG_001_MONOMER'; %complex subunit
                'MG_003_MONOMER'; %N-terminal methionine cleavage
                'MG_007_MONOMER'; %complex subunit
                'MG_031_MONOMER'; %complex subunit
                'MG_038_MONOMER'; %complex subunit
                'MG_046_MONOMER'; %prosthetic metal ion
                'MG_048_MONOMER'; %complex subunit
                'MG_063_MONOMER'; %complex subunit
                'MG_067_MONOMER'; %lipoprotein
                'MG_070_MONOMER'; %phosphorlyated
                'MG_074_MONOMER'; %extracellular
                'MG_120_MONOMER'; %membrane
                'MG_124_MONOMER'; %cytosolic
                'MG_160_MONOMER'; %alpha glutamate ligation
                'MG_191_MONOMER'; %terminal organelle, cytosol
                'MG_261_MONOMER'; %complex subunit
                'MG_271_MONOMER'; %complex subunit
                'MG_272_MONOMER'; %complex subunit
                'MG_273_MONOMER'; %complex subunit
                'MG_274_MONOMER'; %complex subunit
                'MG_315_MONOMER'; %complex subunit
                'MG_318_MONOMER'; %terminal organelle, membrane
                'MG_419_MONOMER'; %complex subunit
                }, 10, 1);
            
            m.monomer.signalSequenceIndexs = (1:23)' + 3*23;
            m.monomer.matureIndexs         = (1:23)' + 5*23;
            m.monomer.inactivatedIndexs    = (1:23)' + 6*23;
            m.monomer.misfoldedIndexs      = (1:23)' + 8*23;
            m.monomer.damagedIndexs        = (1:23)' + 9*23;
            
            m.monomerMolecularWeights = [
                repmat(1e5 * [0.4429; 0.7346; 0.2939; 1.6750; 0.5690; 0.3472; 0.5021; 0.3393; 0.5709; 0.3315; 0.1386; 0.5886; 0.1150; 0.0999; 1.5966; 1.0044; 0.4999; 0.4164; 0.3597; 0.4059; 0.3407; 0.3216; 0.69056], 10, 1);
                repmat(1e5 * 0.0998, 5, 1)];
            
            m.monomer.decayRates = repmat(0.96271e-5, size(m.monomer.wholeCellModelIDs));
            m.monomer.decayRates(m.monomer.signalSequenceIndexs) = Inf;
            m.monomer.decayRates(m.monomer.misfoldedIndexs)      = Inf;
            m.monomer.decayRates(m.monomer.damagedIndexs)        = Inf;
            
            monomerLengths = [
                repmat([380; 649; 254; 1451; 508; 315; 446; 303; 494; 284; 114; 520; 102;  85; 1444; 874; 456; 384; 325; 357; 297; 280; 597], 10, 1);
                repmat(84,5,1)];
            m.monomerLonProteaseCleavages = ceil(monomerLengths/m.lonProteaseFragmentLength) - 1;
            
            m.monomerDecayReactions = zeros(numel(m.substrateWholeCellModelIDs), numel(m.monomer.wholeCellModelIDs));
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_001_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([ 6  8 14 46  3 38 12 10  4 43 33 40  5 31  9 37 13  0 11 16  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_003_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([34 34 39 41  4 53 23 44 13 43 49 49 12 36 22 43 33  5 19 53  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_007_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([ 9  9  8 27  6 16  9  5  2 28 29 24  0 16  7 22 15  3  7 11  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_031_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([84 42 87 96 15 106 50 65 38 131 138 152 25 89 48 68 74 16 54 72  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_038_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([45 10 25 30  7 28 23 27  9 41 49 47  9 21 17 32 29 16 10 32  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_046_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([25  7 19 18  7 11  8 21 10 30 37 24  3 15 11 27 13  2 12 14  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_048_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([26 17 20 21  3 31 26 26  4 35 56 49 20 15 12 25 28  3  4 24  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_063_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([16  2 19 22  2 11  9 11  2 26 44 32  4 15  9 18 25  2  9 24  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_067_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([23 13 29 49  2 19 28 30 10 30 47 40  3 32 19 38 44  7 31 21  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_070_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([14 13 12 32  2 19 17  9  3 22 34 27  4  9 12  9 17  2  8 18  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_074_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([ 5  6  5 11  3  7  6  1  0  5 20 12  0 19  0 12  4  3  7 10  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_120_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([36 14  9 33  2 18 19 30  7 48 57 35 12 51 12 44 24 10 16 42  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_124_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([10  2  7  6  2  7  3  3  1  9  5 10  1  8  4  5  6  2  1  9  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_160_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([ 5  9  2  4  0  3  3  3  5  7  5 12  0  2  2  4  6  1  2  9  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_191_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([74 35 74 130 0 53 67 103 18 63 130 106 15 77 94 129 121 28 41 85  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_261_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([45 20 53 79  9 35 59 27 11 80 116 75  9 51 21 47 51 10 27 48  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_271_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([36  8 23 28  8 25 22 34  7 44 37 39  8 15 12 20 29  1 18 42  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_272_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([37  6 18 17  1 29 12 23  5 36 21 37  5 14 17 19 34  1  7 44  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_273_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([27 15 13 13  1 23 13 26  3 27 27 19  7 17 16 19 18  2 12 27  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_274_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([27 13 19 23  3 27 19 24  6 22 31 26 14 12 14 17 17  4 19 20  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_315_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([12  3 14 26  6 13 16  9  4 28 41 31  3 18  8 18 19  3  9 15  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_318_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([11 16  3 30  0 12 45 15  4 13 24  9  9 21 34 10  8  1  1 13  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_419_MONOMER', m.monomer.wholeCellModelIDs)) = repmat([32 14 33 52  5 42 26 16  9 55 61 69  4 46 10 41 38  4 16 23  1]',1,10);
            m.monomerDecayReactions(m.substrateIndexs_aminoAcids, strcmp('MG_287_MONOMER_ACP', m.monomer.wholeCellModelIDs)) = repmat([3 1 6 8 0 8 2 2 1 7 15 11 3 4 1 5 2 0 0 5 0]',1,5);
            
            m.monomerDecayReactions(m.substrateIndexs_fmethionine, 24:230) = m.monomerDecayReactions(m.substrateIndexs_fmethionine, 24:230) -1;
            
            [~, idxs] = ismember({
                'MG_007_MONOMER'; 'MG_031_MONOMER'; 'MG_038_MONOMER'; 'MG_046_MONOMER'; 'MG_048_MONOMER';
                'MG_063_MONOMER'; 'MG_067_MONOMER'; 'MG_070_MONOMER'; 'MG_074_MONOMER'; 'MG_120_MONOMER';
                'MG_124_MONOMER'; 'MG_160_MONOMER'; 'MG_191_MONOMER'; 'MG_261_MONOMER'; 'MG_272_MONOMER';
                'MG_315_MONOMER'; 'MG_318_MONOMER'; 'MG_419_MONOMER'}, m.monomer.wholeCellModelIDs(24:230));
            m.monomerDecayReactions(m.substrateIndexs_methionine, 23 + idxs) = ...
                m.monomerDecayReactions(m.substrateIndexs_methionine, 23 + idxs) + 1;
            
            m.monomerDecayReactions(strcmp('MG', m.substrateWholeCellModelIDs), find(strcmp('MG_038_MONOMER', m.monomer.wholeCellModelIDs),6,'last')) = 1;
            m.monomerDecayReactions(strcmp('MG', m.substrateWholeCellModelIDs), find(strcmp('MG_063_MONOMER', m.monomer.wholeCellModelIDs),6,'last')) = 1;
            m.monomerDecayReactions(strcmp('MG', m.substrateWholeCellModelIDs), find(strcmp('MG_273_MONOMER', m.monomer.wholeCellModelIDs),6,'last')) = 1;
            m.monomerDecayReactions(strcmp('ZN', m.substrateWholeCellModelIDs), find(strcmp('MG_046_MONOMER', m.monomer.wholeCellModelIDs),6,'last')) = 1;
            
            m.monomerDecayReactions(strcmp('LIPOYLLYS', m.substrateWholeCellModelIDs), find(strcmp('MG_272_MONOMER', m.monomer.wholeCellModelIDs), 5, 'last')) = 1;
            m.monomerDecayReactions(strcmp('diacylglycerolCys', m.substrateWholeCellModelIDs), find(strcmp('MG_067_MONOMER', m.monomer.wholeCellModelIDs), 5, 'last')) = 1;
            m.monomerDecayReactions(strcmp('pSER', m.substrateWholeCellModelIDs), find(strcmp('MG_070_MONOMER', m.monomer.wholeCellModelIDs), 5, 'last')) = 2;
            m.monomerDecayReactions(strcmp('pSER', m.substrateWholeCellModelIDs), find(strcmp('MG_273_MONOMER', m.monomer.wholeCellModelIDs), 5, 'last')) = 1;
            m.monomerDecayReactions(strcmp('pTHR', m.substrateWholeCellModelIDs), find(strcmp('MG_070_MONOMER', m.monomer.wholeCellModelIDs), 5, 'last')) = 3;
            m.monomerDecayReactions(strcmp('pTHR', m.substrateWholeCellModelIDs), find(strcmp('MG_274_MONOMER', m.monomer.wholeCellModelIDs), 5, 'last')) = 1;
            
            m.monomerDecayReactions(strcmp('COA', m.substrateWholeCellModelIDs), strcmp('MG_287_MONOMER_ACP', m.monomer.wholeCellModelIDs)) = 1;
            m.monomerDecayReactions(strcmp('PAP', m.substrateWholeCellModelIDs), strcmp('MG_287_MONOMER_ACP', m.monomer.wholeCellModelIDs)) = -1;
            m.monomerDecayReactions(strcmp('H', m.substrateWholeCellModelIDs), strcmp('MG_287_MONOMER_ACP', m.monomer.wholeCellModelIDs)) = -1;
            
            m.monomerDecayReactions(m.substrateIndexs_water, :) = ...
                m.monomerDecayReactions(m.substrateIndexs_water, :) - (monomerLengths - 1)';
            
            %% RNAs
            m.rnaWholeCellModelIDs = {'MG_0001'};
            m.rnaMolecularWeights = 1e5 * 2.9624;
            
            %% complexes
            m.complex.wholeCellModelIDs = repmat({
                'DNA_POLYMERASE_HOLOENZYME'; %complex with complex subunits
                'MG_0001_048'; %complex with RNA subunits
                'MG_038_TETRAMER'; %prosthetic metal ion
                'MG_063_DIMER'; %complex
                'MG_271_272_273_274_192MER'; %complex including modified monomers
                },6,1);
            
            m.complex.matureIndexs                      = (1:5)' + 1*5;
            m.complex.inactivatedIndexs                 = (1:5)' + 2*5;
            m.complex.misfoldedIndexs                   = (1:5)' + 4*5;
            m.complex.damagedIndexs                     = (1:5)' + 5*5;
            
            m.complexMolecularWeights = [
                repmat(1e6 * [0.9641; 0.0766; 0.2276; 0.0679; 7.6916],6,1);
                repmat(1e6 * 0.0348,5,1)];
            
            m.complex.decayRates = repmat(0.96271e-5, size(m.complex.wholeCellModelIDs));
            m.complex.decayRates(m.complex.misfoldedIndexs) = Inf;
            m.complex.decayRates(m.complex.damagedIndexs)   = Inf;
            
            m.complexDecayReactions = zeros(numel(m.substrateWholeCellModelIDs), numel(m.complex.wholeCellModelIDs));
            m.complexDecayReactions(strcmp('H',m.substrateWholeCellModelIDs), find(strcmp('MG_454_DIMER_ox', m.complex.wholeCellModelIDs),5,'last')) = -4;
            m.complexDecayReactions(strcmp('ZN',m.substrateWholeCellModelIDs), strcmp('MG_038_TETRAMER', m.complex.wholeCellModelIDs)) = 2;
            
            m.proteinComplexRNAComposition = zeros(numel(m.rnaWholeCellModelIDs), numel(m.complex.wholeCellModelIDs));
            m.proteinComplexMonomerComposition = zeros(numel(m.monomer.damagedIndexs), numel(m.complex.wholeCellModelIDs));
            m.proteinComplexMonomerComposition(strcmp('MG_001_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('DNA_POLYMERASE_HOLOENZYME',m.complex.wholeCellModelIDs)) = 2;
            m.proteinComplexMonomerComposition(strcmp('MG_031_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('DNA_POLYMERASE_HOLOENZYME',m.complex.wholeCellModelIDs)) = 2;
            m.proteinComplexMonomerComposition(strcmp('MG_261_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('DNA_POLYMERASE_HOLOENZYME',m.complex.wholeCellModelIDs)) = 2;
            m.proteinComplexMonomerComposition(strcmp('MG_007_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('DNA_POLYMERASE_HOLOENZYME',m.complex.wholeCellModelIDs)) = 1;
            m.proteinComplexMonomerComposition(strcmp('MG_315_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('DNA_POLYMERASE_HOLOENZYME',m.complex.wholeCellModelIDs)) = 1;
            m.proteinComplexMonomerComposition(strcmp('MG_419_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('DNA_POLYMERASE_HOLOENZYME',m.complex.wholeCellModelIDs)) = 4;
            m.proteinComplexRNAComposition(    strcmp('MG_0001',       m.rnaWholeCellModelIDs),                             strcmp('MG_0001_048',              m.complex.wholeCellModelIDs)) = 1;
            m.proteinComplexMonomerComposition(strcmp('MG_048_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('MG_0001_048',              m.complex.wholeCellModelIDs)) = 1;
            m.proteinComplexMonomerComposition(strcmp('MG_038_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('MG_038_TETRAMER',          m.complex.wholeCellModelIDs)) = 4;
            m.proteinComplexMonomerComposition(strcmp('MG_063_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('MG_063_DIMER',             m.complex.wholeCellModelIDs)) = 2;
            m.proteinComplexMonomerComposition(strcmp('MG_271_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('MG_271_272_273_274_192MER',m.complex.wholeCellModelIDs)) = 12;
            m.proteinComplexMonomerComposition(strcmp('MG_272_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('MG_271_272_273_274_192MER',m.complex.wholeCellModelIDs)) = 60;
            m.proteinComplexMonomerComposition(strcmp('MG_273_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('MG_271_272_273_274_192MER',m.complex.wholeCellModelIDs)) = 60;
            m.proteinComplexMonomerComposition(strcmp('MG_274_MONOMER',m.monomer.wholeCellModelIDs(m.monomer.damagedIndexs)), strcmp('MG_271_272_273_274_192MER',m.complex.wholeCellModelIDs)) = 60;
            
            %% enzymes
            m.enzymeWholeCellModelIDs = {
                'MG_355_HEXAMER';    %protease clpB            cytoplasmic
                'MG_239_HEXAMER';    %protease La              cytoplasmic
                'MG_457_HEXAMER';    %metalloprotease FtsH     integral membrane
                'MG_183_MONOMER';    %oligoendopeptidase F     cytoplasmic
                'MG_324_MONOMER';    %aminopeptidase           cytoplasmic
                'MG_391_HEXAMER';    %cytosol aminopeptidase   cytoplasmic
                'MG_208_DIMER';      %glycoprotease            cytoplasmic
                'MG_046_DIMER';      %metalloendopeptidase     cytoplasmic
                'MG_020_MONOMER'};   %proline iminopeptidase   cytoplasmic
            m.enzymeNames = m.enzymeWholeCellModelIDs;
            
            m.enzymeIndexs_clpBProtease          = 1;
            m.enzymeIndexs_lonProtease           = 2;
            m.enzymeIndexs_ftsHProtease          = 3;
            m.enzymeIndexs_oligoendopeptidaseF   = 4;
            m.enzymeIndexs_aminopeptidase        = 5;
            m.enzymeIndexs_cytosolAminopeptidase = 6;
            m.enzymeIndexs_glycoprotease         = 7;
            m.enzymeIndexs_metalloendopeptidase  = 8;
            m.enzymeIndexs_prolineIminopeptidase = 9;
            m.enzymeIndexs_peptidases            = (4:9)';
            
            m.enzymeMolecularWeights = 1e5 * [4.8626; 5.3992; 4.6052; 0.7094; 0.3985; 2.9486; 0.4507; 0.6945; 0.3530];
            
            %% initial state
            numCompartments  = length(m.compartment.wholeCellModelIDs);
            m.substrates     = zeros(length(m.substrateWholeCellModelIDs), 1);
            m.enzymes        = zeros(length(m.enzymeWholeCellModelIDs),    1);
            m.boundEnzymes   = zeros(length(m.enzymeWholeCellModelIDs),    1);
            
            m.RNAs           = zeros(length(m.rnaWholeCellModelIDs),      numCompartments);
            m.monomers       = zeros(length(m.monomer.wholeCellModelIDs), numCompartments);
            m.complexs       = zeros(length(m.complex.wholeCellModelIDs), numCompartments);
            
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
        end
    end
    
    %tests
    methods
        function testConstants(this)
            m = this.process;
            comp = m.compartment;
            pm = m.monomer;
            pc = m.complex;
            
            %check metabolite weights
            assertTrue(all(m.substrateMolecularWeights > 0));
            assertEqual(24.305, m.substrateMolecularWeights(m.substrateIndexs({'MG'})));
            assertEqual(65.39,  m.substrateMolecularWeights(m.substrateIndexs({'ZN'})));
            
            %check composition of complexes with complex subunits
            assertEqual(repmat(2,1,6), m.proteinComplexMonomerComposition(strcmp('MG_001_MONOMER', pm.wholeCellModelIDs(pm.damagedIndexs)), strcmp('DNA_POLYMERASE_HOLOENZYME', pc.wholeCellModelIDs)));
            assertEqual(repmat(2,1,6), m.proteinComplexMonomerComposition(strcmp('MG_031_MONOMER', pm.wholeCellModelIDs(pm.damagedIndexs)), strcmp('DNA_POLYMERASE_HOLOENZYME', pc.wholeCellModelIDs)));
            assertEqual(repmat(2,1,6), m.proteinComplexMonomerComposition(strcmp('MG_261_MONOMER', pm.wholeCellModelIDs(pm.damagedIndexs)), strcmp('DNA_POLYMERASE_HOLOENZYME', pc.wholeCellModelIDs)));
            assertEqual(repmat(1,1,6), m.proteinComplexMonomerComposition(strcmp('MG_007_MONOMER', pm.wholeCellModelIDs(pm.damagedIndexs)), strcmp('DNA_POLYMERASE_HOLOENZYME', pc.wholeCellModelIDs))); %#ok<RPMT1>
            assertEqual(repmat(1,1,6), m.proteinComplexMonomerComposition(strcmp('MG_315_MONOMER', pm.wholeCellModelIDs(pm.damagedIndexs)), strcmp('DNA_POLYMERASE_HOLOENZYME', pc.wholeCellModelIDs))); %#ok<RPMT1>
            assertEqual(repmat(4,1,6), m.proteinComplexMonomerComposition(strcmp('MG_419_MONOMER', pm.wholeCellModelIDs(pm.damagedIndexs)), strcmp('DNA_POLYMERASE_HOLOENZYME', pc.wholeCellModelIDs)));
            
            %check monomer decay reactions >0 except water, meaning monomers only require energy and water to be degraded
            assertEqual({'H2O'},m.substrateWholeCellModelIDs(any(m.monomerDecayReactions < 0,2)));
            
            %check decay rates
            infTfs = false(size(pm.decayRates));
            infTfs([pm.signalSequenceIndexs; pm.misfoldedIndexs; pm.damagedIndexs]) = true;
            infTfs(~ismember(pm.compartments, [comp.cytosolIndexs; comp.terminalOrganelleCytosolIndexs])) = false;
            assertAllEqual(0, pm.decayRates(~ismember(pm.compartments, [comp.cytosolIndexs; comp.terminalOrganelleCytosolIndexs])));
            assertAllEqual(true, isfinite(pm.decayRates(~infTfs)));
            assertAllEqual(Inf, pm.decayRates(infTfs));
            
            infTfs = false(size(pc.decayRates));
            infTfs([pc.misfoldedIndexs; pc.damagedIndexs]) = true;
            infTfs(~ismember(pc.compartments, [comp.cytosolIndexs; comp.terminalOrganelleCytosolIndexs])) = false;
            assertAllEqual(0, pc.decayRates(~ismember(pc.compartments, [comp.cytosolIndexs; comp.terminalOrganelleCytosolIndexs])));
            assertAllEqual(true, isfinite(pc.decayRates(~infTfs)));
            assertAllEqual(Inf, pc.decayRates(infTfs));
        end
    end
    
    %misfolding/refolding tests
    methods
        function testMisfoldingOfOneInactivatedComplex(this)
            m = this.process;
            m.proteinMisfoldingRate = 1e6;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.complexs(m.complex.inactivatedIndexs(1), m.compartment.cytosolIndexs) = 1;
            m.RNAs(:) = 0;
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_complexs(m.complex.inactivatedIndexs(1), m.compartment.cytosolIndexs) = 0;
            final_complexs(m.complex.misfoldedIndexs(1), m.compartment.cytosolIndexs) = 1;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            m.evolveState_MisfoldProteins();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testMisfoldingAllProteins(this)
            m = this.process;
            m.proteinMisfoldingRate = 1e6;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 1;
            m.complexs(:) = 1;
            m.monomer.chromosome.initialize();
            m.monomer.chromosome.monomerBoundSites([2000*(1:numel(m.monomer.boundIndexs))' ones(size(m.monomer.boundIndexs))]) = 1:numel(m.monomer.boundIndexs);
            m.monomer.chromosome.complexBoundSites([1000+2000*(1:numel(m.complex.boundIndexs))' ones(size(m.complex.boundIndexs))]) = 1:numel(m.complex.boundIndexs);
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            m.complex.ftsZRing.numEdgesOneStraight = 1;
            m.complex.ftsZRing.numResidualBent = 1;
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers([m.monomer.matureIndexs; m.monomer.inactivatedIndexs; m.monomer.boundIndexs], setdiff(1:end, [m.compartment.membraneIndexs; m.compartment.terminalOrganelleMembraneIndexs])) = 0;
            final_monomers(m.monomer.misfoldedIndexs, setdiff(1:end, [m.compartment.membraneIndexs; m.compartment.terminalOrganelleMembraneIndexs]))=4;
            final_complexs = m.complexs;
            idxs = [m.complex.replisomeIndexs; m.complex.ftsZGTPIndexs; m.complex.ftsZGDPIndexs; m.complex.dnaAPolymerIndexs];
            final_complexs([m.complex.matureIndexs; m.complex.inactivatedIndexs; m.complex.boundIndexs], setdiff(1:end, [m.compartment.membraneIndexs; m.compartment.terminalOrganelleMembraneIndexs])) = 0;            
            final_complexs(m.complex.misfoldedIndexs, setdiff(1:end, [m.compartment.membraneIndexs; m.compartment.terminalOrganelleMembraneIndexs]))=4;
            final_complexs(m.complex.boundIndexs(idxs), m.compartment.cytosolIndexs) = 1;
            final_complexs(m.complex.misfoldedIndexs(idxs), m.compartment.cytosolIndexs) = 3;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = 1 + numel(m.polypeptide.abortedSequences);

            warnState = warning('query', 'WholeCell:warning');
            warning('off', 'WholeCell:warning');
            m.evolveState_MisfoldProteins();
            warning(warnState.state, 'WholeCell:warning');
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, numel(m.polypeptide.abortedSequences));
        end
        
        function testMisfoldingAndRepairOfOneMatureMonomer(this)
            m = this.process;
            m.proteinMisfoldingRate = 1e6;
            m.substrates(:)=0;
            m.substrates(m.substrateIndexs_atp)=1;
            m.substrates(m.substrateIndexs_water)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_clpBProtease) = 1;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(m.monomer.matureIndexs(1), m.compartment.cytosolIndexs) = 1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_substrates(m.substrateIndexs_atp)=0;
            final_substrates(m.substrateIndexs_water)=0;
            final_substrates(m.substrateIndexs_adp)=1;
            final_substrates(m.substrateIndexs_phosphate)=1;
            final_substrates(m.substrateIndexs_hydrogen)=1;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % misfold
            m.evolveState_MisfoldProteins();
            assertEqual(0, m.monomers(m.monomer.matureIndexs(1), m.compartment.cytosolIndexs));
            assertEqual(1, m.monomers(m.monomer.misfoldedIndexs(1), m.compartment.cytosolIndexs));
            
            % repair
            m.evolveState_RefoldProteins();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testMisfoldingAndRepairOfOneMatureMonomer_NoATP(this)
            m = this.process;
            m.proteinMisfoldingRate = 1e6;
            m.substrates(:)=0;
            m.substrates(m.substrateIndexs_atp)=0;
            m.substrates(m.substrateIndexs_water)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_clpBProtease) = 1;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(m.monomer.matureIndexs(1), m.compartment.cytosolIndexs) = 1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(m.monomer.matureIndexs(1), m.compartment.cytosolIndexs)=0;
            final_monomers(m.monomer.misfoldedIndexs(1), m.compartment.cytosolIndexs)=1;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % misfold
            m.evolveState_MisfoldProteins();
            assertEqual(0, m.monomers(m.monomer.matureIndexs(1), m.compartment.cytosolIndexs));
            assertEqual(1, m.monomers(m.monomer.misfoldedIndexs(1), m.compartment.cytosolIndexs));
            
            % repair
            m.evolveState_RefoldProteins();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testMisfoldingAndRepairOfOneMatureMonomer_NoEnzyme(this)
            m = this.process;
            m.proteinMisfoldingRate = 1e6;
            m.substrates(:)=0;
            m.substrates(m.substrateIndexs_atp)=1;
            m.substrates(m.substrateIndexs_water)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_clpBProtease) = 0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(m.monomer.matureIndexs(1), m.compartment.cytosolIndexs) = 1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(m.monomer.matureIndexs(1), m.compartment.cytosolIndexs)=0;
            final_monomers(m.monomer.misfoldedIndexs(1), m.compartment.cytosolIndexs)=1;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % misfold
            m.evolveState_MisfoldProteins();
            assertEqual(0, m.monomers(m.monomer.matureIndexs(1), m.compartment.cytosolIndexs));
            assertEqual(1, m.monomers(m.monomer.misfoldedIndexs(1), m.compartment.cytosolIndexs));
            
            % repair
            m.evolveState_RefoldProteins();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testMisfoldingAndRepairOfOneMatureMonomer_InMembrane(this)
            m = this.process;
            m.proteinMisfoldingRate = 1e6;
            m.substrates(:)=0;
            m.substrates(m.substrateIndexs_atp)=1;
            m.substrates(m.substrateIndexs_water)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_clpBProtease) = 1;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(m.monomer.matureIndexs(1), m.compartment.membraneIndexs) = 1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(m.monomer.matureIndexs(1), m.compartment.membraneIndexs)=1;
            final_monomers(m.monomer.misfoldedIndexs(1), m.compartment.membraneIndexs)=0;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % misfold
            m.evolveState_MisfoldProteins();
            assertEqual(1, m.monomers(m.monomer.matureIndexs(1), m.compartment.membraneIndexs));
            assertEqual(0, m.monomers(m.monomer.misfoldedIndexs(1), m.compartment.membraneIndexs));
            
            % repair
            m.evolveState_RefoldProteins();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testMisfoldingAndRepairOfOneMatureMonomer_InTermOrgCytosol(this)
            m = this.process;
            m.proteinMisfoldingRate = 1e6;
            m.substrates(:)=0;
            m.substrates(m.substrateIndexs_atp)=1;
            m.substrates(m.substrateIndexs_water)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_clpBProtease) = 1;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(m.monomer.matureIndexs(1), m.compartment.terminalOrganelleCytosolIndexs) = 1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_substrates(m.substrateIndexs_atp)=0;
            final_substrates(m.substrateIndexs_water)=0;
            final_substrates(m.substrateIndexs_adp)=1;
            final_substrates(m.substrateIndexs_phosphate)=1;
            final_substrates(m.substrateIndexs_hydrogen)=1;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % misfold
            m.evolveState_MisfoldProteins();
            assertEqual(0, m.monomers(m.monomer.matureIndexs(1), m.compartment.terminalOrganelleCytosolIndexs));
            assertEqual(1, m.monomers(m.monomer.misfoldedIndexs(1), m.compartment.terminalOrganelleCytosolIndexs));
            
            % repair
            m.evolveState_RefoldProteins();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testMisfoldingAndRepairOfOneMatureMonomer_InTermOrgMembrane(this)
            m = this.process;
            m.proteinMisfoldingRate = 1e6;
            m.substrates(:)=0;
            m.substrates(m.substrateIndexs_atp)=1;
            m.substrates(m.substrateIndexs_water)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_clpBProtease) = 1;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(m.monomer.matureIndexs(1), m.compartment.terminalOrganelleMembraneIndexs) = 1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(m.monomer.matureIndexs(1), m.compartment.terminalOrganelleMembraneIndexs)=1;
            final_monomers(m.monomer.misfoldedIndexs(1), m.compartment.terminalOrganelleMembraneIndexs)=0;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % misfold
            m.evolveState_MisfoldProteins();
            assertEqual(1, m.monomers(m.monomer.matureIndexs(1), m.compartment.terminalOrganelleMembraneIndexs));
            assertEqual(0, m.monomers(m.monomer.misfoldedIndexs(1), m.compartment.terminalOrganelleMembraneIndexs));
            
            % repair
            m.evolveState_RefoldProteins();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testMisfoldingAndRepairOfOneInactiveComplexMonomer(this)
            m = this.process;
            m.proteinMisfoldingRate = 1e6;
            m.substrates(:)=0;
            m.substrates(m.substrateIndexs_atp)=1;
            m.substrates(m.substrateIndexs_water)=1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_clpBProtease) = 1;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.complexs(m.complex.inactivatedIndexs(1), m.compartment.cytosolIndexs) = 1;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_substrates(m.substrateIndexs_atp)=0;
            final_substrates(m.substrateIndexs_water)=0;
            final_substrates(m.substrateIndexs_adp)=1;
            final_substrates(m.substrateIndexs_phosphate)=1;
            final_substrates(m.substrateIndexs_hydrogen)=1;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_complexs(m.complex.inactivatedIndexs(1), m.compartment.cytosolIndexs) = 0;
            final_complexs(m.complex.matureIndexs(1), m.compartment.cytosolIndexs) = 1;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % misfold
            m.evolveState_MisfoldProteins();
            assertEqual(0, m.complexs(m.complex.inactivatedIndexs(1), m.compartment.cytosolIndexs));
            assertEqual(0, m.complexs(m.complex.matureIndexs(1), m.compartment.cytosolIndexs));
            assertEqual(1, m.complexs(m.complex.misfoldedIndexs(1), m.compartment.cytosolIndexs));
            
            % repair
            m.evolveState_RefoldProteins();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testMisfoldingAndRepairOfAllProteins(this)
            m = this.process;
            idxs = [m.complex.replisomeIndexs; m.complex.ftsZGTPIndexs; m.complex.ftsZGDPIndexs; m.complex.dnaAPolymerIndexs];
            nProteins = 2 * 4 * (numel(m.monomer.matureIndexs) + numel(m.complex.matureIndexs)) - numel(idxs);
            m.proteinMisfoldingRate = 1e6;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = nProteins;
            m.substrates(m.substrateIndexs_water) = nProteins;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_clpBProtease) = 1;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 1;
            m.complexs(:) = 1;
            m.monomer.chromosome.initialize();
            m.monomer.chromosome.monomerBoundSites([2000*(1:numel(m.monomer.boundIndexs))' ones(size(m.monomer.boundIndexs))]) = 1:numel(m.monomer.boundIndexs);
            m.monomer.chromosome.complexBoundSites([1000+2000*(1:numel(m.complex.boundIndexs))' ones(size(m.complex.boundIndexs))]) = 1:numel(m.complex.boundIndexs);
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            m.complex.ftsZRing.numEdgesOneStraight = 1;
            m.complex.ftsZRing.numResidualBent = 1;
            
            final_substrates = m.substrates;
            final_substrates(m.substrateIndexs_atp) = 0;
            final_substrates(m.substrateIndexs_water) = 0;
            final_substrates(m.substrateIndexs_adp) = nProteins;
            final_substrates(m.substrateIndexs_phosphate) = nProteins;
            final_substrates(m.substrateIndexs_hydrogen) = nProteins;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            
            final_monomers = m.monomers;
            final_monomers(m.monomer.misfoldedIndexs, :) = 4;
            final_monomers(m.monomer.matureIndexs, :) = 0;
            final_monomers(m.monomer.boundIndexs, :) = 0;
            final_monomers(m.monomer.inactivatedIndexs, :) = 0;
            final_monomers(m.monomer.matureIndexs, [m.compartment.cytosolIndexs; m.compartment.terminalOrganelleCytosolIndexs]) = 4;
            final_monomers(m.monomer.misfoldedIndexs, [m.compartment.cytosolIndexs; m.compartment.terminalOrganelleCytosolIndexs]) = 0;
            final_monomers(m.monomer.matureIndexs, [m.compartment.membraneIndexs; m.compartment.terminalOrganelleMembraneIndexs]) = 1;
            final_monomers(m.monomer.boundIndexs, [m.compartment.membraneIndexs; m.compartment.terminalOrganelleMembraneIndexs]) = 1;
            final_monomers(m.monomer.misfoldedIndexs, [m.compartment.membraneIndexs; m.compartment.terminalOrganelleMembraneIndexs]) = 1;
            final_monomers(m.monomer.inactivatedIndexs, [m.compartment.membraneIndexs; m.compartment.terminalOrganelleMembraneIndexs]) = 1;
                        
            final_complexs = m.complexs;
            final_complexs(m.complex.misfoldedIndexs, :) = 4;
            final_complexs(m.complex.matureIndexs, :) = 0;
            final_complexs(m.complex.boundIndexs, :) = 0;
            final_complexs(m.complex.inactivatedIndexs, :) = 0;
            final_complexs(m.complex.matureIndexs, [m.compartment.cytosolIndexs; m.compartment.terminalOrganelleCytosolIndexs]) = 4;
            final_complexs(m.complex.misfoldedIndexs, [m.compartment.cytosolIndexs; m.compartment.terminalOrganelleCytosolIndexs]) = 0;
            final_complexs(m.complex.matureIndexs, [m.compartment.membraneIndexs; m.compartment.terminalOrganelleMembraneIndexs]) = 1;
            final_complexs(m.complex.boundIndexs, [m.compartment.membraneIndexs; m.compartment.terminalOrganelleMembraneIndexs]) = 1;
            final_complexs(m.complex.misfoldedIndexs, [m.compartment.membraneIndexs; m.compartment.terminalOrganelleMembraneIndexs]) = 1;
            final_complexs(m.complex.inactivatedIndexs, [m.compartment.membraneIndexs; m.compartment.terminalOrganelleMembraneIndexs]) = 1;
            final_complexs(m.complex.boundIndexs(idxs), m.compartment.cytosolIndexs) = 1;
            final_complexs(m.complex.matureIndexs(idxs), m.compartment.cytosolIndexs) = 3;
            
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = ...
                numel(m.polypeptide.abortedSequences) + ...
                1;
            
            % misfold
            warnState = warning('query', 'WholeCell:warning');
            warning('off', 'WholeCell:warning');
            m.evolveState_MisfoldProteins();
            warning(warnState.state, 'WholeCell:warning');
            
            % repair
            m.evolveState_RefoldProteins();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, numel(m.polypeptide.abortedSequences));
        end
        
        function testMisfoldingAndRepairFairness(this)
            m = this.process;
            m.proteinMisfoldingRate = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_clpBProtease) = 1;
            m.boundEnzymes(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            counts = zeros(2,1);
            for i=1:50
                m.substrates(:)=0;
                m.substrates(m.substrateIndexs_atp)=1;
                m.substrates(m.substrateIndexs_water)=1;
                m.monomers(:) = 0;
                m.complexs(:) = 0;
                m.monomers(m.monomer.inactivatedIndexs(1), m.compartment.cytosolIndexs)=1;
                m.complexs(m.complex.misfoldedIndexs(1), m.compartment.cytosolIndexs)=1;
                
                m.evolveState_MisfoldProteins();
                m.evolveState_RefoldProteins();
                
                if m.monomers(m.monomer.matureIndexs(1), m.compartment.cytosolIndexs)
                    counts(1)=counts(1)+1;
                end
                if m.complexs(m.complex.matureIndexs(1), m.compartment.cytosolIndexs)
                    counts(2)=counts(2)+1;
                end
                assertEqual(1, m.monomers(m.monomer.matureIndexs(1), m.compartment.cytosolIndexs)+ m.complexs(m.complex.matureIndexs(1), m.compartment.cytosolIndexs));
            end
            assertTrue(range(counts) < 0.5 * max(counts));
        end
    end
    
    %complex degradation tests
    methods
        function testDecayofMatureComplex(this)
            m = this.process;
            
            idx = m.complex.matureIndexs(1);
            
            m.complex.decayRates(:)=1e6;
            m.substrates(:)=0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.complexs(idx,m.compartment.cytosolIndexs) =1;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.complexDecayReactions(:, idx);
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(m.monomer.damagedIndexs, m.compartment.cytosolIndexs) = m.proteinComplexMonomerComposition(:,idx);
            final_complexs = m.complexs;
            final_complexs(idx,m.compartment.cytosolIndexs)=0;
            final_RNAs = m.RNAs;
            final_RNAs(:, m.compartment.cytosolIndexs) = m.proteinComplexRNAComposition(:,idx);
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade complexes
            m.evolveState_DegradeComplexes();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayofMatureComplexWithRNASubunits(this)
            m = this.process;
            
            idx = m.complex.matureIndexs(strcmp('MG_0001_048',m.complex.wholeCellModelIDs(m.complex.matureIndexs)));
            
            m.complex.decayRates(:)=1e6;
            m.substrates(:)=0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.complexs(idx,m.compartment.cytosolIndexs) =1;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.complexDecayReactions(:, idx);
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(m.monomer.damagedIndexs, m.compartment.cytosolIndexs) = m.proteinComplexMonomerComposition(:,idx);
            final_complexs = m.complexs;
            final_complexs(idx,m.compartment.cytosolIndexs)=0;
            final_RNAs = m.RNAs;
            final_RNAs(:, m.compartment.cytosolIndexs) = m.proteinComplexRNAComposition(:,idx);
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade complexes
            m.evolveState_DegradeComplexes();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayofComplexWithRNASubunits(this)
            m = this.process;
            
            idx = m.complex.matureIndexs(strcmp('MG_454_DIMER_ox',m.complex.wholeCellModelIDs(m.complex.matureIndexs)));
            
            m.complex.decayRates(:)=1e6;
            m.substrates(:)=0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.complexs(idx,m.compartment.cytosolIndexs) =1;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.complexDecayReactions(:, idx);
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(m.monomer.damagedIndexs, m.compartment.cytosolIndexs) = m.proteinComplexMonomerComposition(:,idx);
            final_complexs = m.complexs;
            final_complexs(idx,m.compartment.cytosolIndexs)=0;
            final_RNAs = m.RNAs;
            final_RNAs(:, m.compartment.cytosolIndexs) = m.proteinComplexRNAComposition(:,idx);
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade complexes
            m.evolveState_DegradeComplexes();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayofMisfoldedComplex_InMembrane(this)
            m = this.process;
            
			idx = m.complex.misfoldedIndexs(find(m.complex.compartments(m.complex.misfoldedIndexs) == m.compartment.membraneIndexs, 1, 'first'));
            
            m.complex.decayRates(m.complex.compartments == m.compartment.cytosolIndexs) = 1e6; 
			m.complex.decayRates(m.complex.compartments == m.compartment.terminalOrganelleCytosolIndexs) = 1e6; 
            
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.complexs(idx, m.compartment.membraneIndexs) = 1;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_RNAs(:, m.compartment.cytosolIndexs) = m.proteinComplexRNAComposition(:, idx);
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade complexes
            m.evolveState_DegradeComplexes();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayofMatureComplex_InTerminalOrganelleCytosol(this)
            m = this.process;
            
            idx = m.complex.matureIndexs(1);
            
            m.complex.decayRates(:)=1e6;
            m.substrates(:)=0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.complexs(idx,m.compartment.terminalOrganelleCytosolIndexs) =1;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.complexDecayReactions(:, idx);
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(m.monomer.damagedIndexs, m.compartment.cytosolIndexs) = m.proteinComplexMonomerComposition(:,idx);
            final_complexs = m.complexs;
            final_complexs(idx,m.compartment.terminalOrganelleCytosolIndexs)=0;
            final_RNAs = m.RNAs;
            final_RNAs(:, m.compartment.cytosolIndexs) = m.proteinComplexRNAComposition(:,idx);
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade complexes
            m.evolveState_DegradeComplexes();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayofMatureComplex_InTerminalOrganelleMembrane(this)
            m = this.process;
            
			idx = m.complex.misfoldedIndexs(find(m.complex.compartments(m.complex.misfoldedIndexs) == m.compartment.membraneIndexs, 1, 'first'));
            
            m.complex.decayRates(m.complex.compartments == m.compartment.cytosolIndexs) = 1e6; 
			m.complex.decayRates(m.complex.compartments == m.compartment.terminalOrganelleCytosolIndexs) = 1e6; 

            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.complexs(idx, m.compartment.terminalOrganelleMembraneIndexs) = 1;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_RNAs(:, m.compartment.cytosolIndexs) = m.proteinComplexRNAComposition(:, idx);
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade complexes
            m.evolveState_DegradeComplexes();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayofInactivatedComplex_InMembrane(this)
            m = this.process;
            
			idx = m.complex.inactivatedIndexs(find(m.complex.compartments(m.complex.inactivatedIndexs) == m.compartment.membraneIndexs, 1, 'first'));
            
            m.complex.decayRates(m.complex.compartments == m.compartment.cytosolIndexs) = 1e6; 
			m.complex.decayRates(m.complex.compartments == m.compartment.terminalOrganelleCytosolIndexs) = 1e6; 
            
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.complexs(idx, m.compartment.membraneIndexs) = 1;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_RNAs(:, m.compartment.cytosolIndexs) = m.proteinComplexRNAComposition(:, idx);
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade complexes
            m.evolveState_DegradeComplexes();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayofMatureComplex_InMembrane(this)
            m = this.process;
            
            idx = m.complex.matureIndexs(find(m.complex.compartments(m.complex.matureIndexs) == m.compartment.membraneIndexs, 1, 'first'));
            
            m.complex.decayRates(m.complex.compartments == m.compartment.cytosolIndexs) = 1e6; 
			m.complex.decayRates(m.complex.compartments == m.compartment.terminalOrganelleCytosolIndexs) = 1e6; 
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.complexs(idx, m.compartment.membraneIndexs) = 1;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_RNAs(:, m.compartment.cytosolIndexs) = m.proteinComplexRNAComposition(:, idx);
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade complexes
            m.evolveState_DegradeComplexes();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
    end
    
    %proteolysis tagged monomer degradation tests
    methods
        function testDecayOfTaggedMonomer(this)
            m = this.process;
            
            m.monomer.decayRates(:) = 1e6;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ftsHProtease) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidases) = 1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = [1 30 0];
            
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = cell(0, 1);
            
            % degrade monomers
            m.evolveState_DegradeAbortedPolypeptides();
            
            % assert
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayOfMultipleTaggedMonomers(this)
            m = this.process;
            
            m.monomer.decayRates(:) = 1e6;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ftsHProtease) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidases) = 1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = repmat([1 30 0], 6, 1);
            
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = cell(0, 1);
            
            % degrade monomers
            m.evolveState_DegradeAbortedPolypeptides();
            
            % assert
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayOfTaggedMonomer_LimitedWater(this)
            m = this.process;
            
            m.monomer.decayRates(:) = 1e6;
            m.substrates(:) = 1e6;            
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ftsHProtease) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidases) = 1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = [1 30 0];
            m.substrates(m.substrateIndexs_water) = ...
                + floor(m.ftsHProteaseEnergyCost * ceil(length(m.polypeptide.abortedSequences{1}) / m.ftsHProteaseFragmentLength - 1)) ...
                + length(m.polypeptide.abortedSequences{1}) - 1 ...
                - 1;
            
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade monomers
            m.evolveState_DegradeAbortedPolypeptides();
            
            % assert
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayOfTaggedMonomer_NoATP(this)
            m = this.process;
            
            m.monomer.decayRates(:) = 1e6;
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_atp) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ftsHProtease) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidases) = 1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = [1 30 0];
            
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade monomers
            m.evolveState_DegradeAbortedPolypeptides();
            
            % assert
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayOfTaggedMonomer_NoProteaseFtsH(this)
            m = this.process;
            
            m.monomer.decayRates(:) = 1e6;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ftsHProtease) = 0;
            m.enzymes(m.enzymeIndexs_peptidases) = 1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = [1 30 0];
            
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade monomers
            m.evolveState_DegradeAbortedPolypeptides();
            
            % assert
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayOfTaggedMonomer_NoPeptidase(this)
            m = this.process;
            
            m.monomer.decayRates(:) = 1e6;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ftsHProtease) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidases) = 0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = [1 30 0];
            
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade monomers
            m.evolveState_DegradeAbortedPolypeptides();
            
            % assert
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
    end
    
    %monomer degradation tests
    methods
        function testDecayOfMonomer_InCytosol(this)
            m = this.process;
            
            idx = m.monomer.matureIndexs(1);
            energyIdxs = [m.substrateIndexs_atp;m.substrateIndexs_adp;m.substrateIndexs_phosphate;m.substrateIndexs_water;m.substrateIndexs_hydrogen];
            
            m.monomer.decayRates(:)=1e6;
            m.substrates = max(0,-m.monomerDecayReactions(:,idx));
            m.substrates(energyIdxs) = m.substrates(energyIdxs) + ...
                max([1 -1 -1 1 -1]',0) * m.monomerLonProteaseCleavages(idx) * m.lonProteaseEnergyCost;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lonProtease)=1e6;
            m.enzymes(m.enzymeIndexs_peptidases)=1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(idx,m.compartment.cytosolIndexs) =1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_substrates = final_substrates + ...
                m.monomerDecayReactions(:,idx);
            final_substrates(energyIdxs) = ...
                final_substrates(energyIdxs) - ...
                [1 -1 -1 1 -1]' * m.monomerLonProteaseCleavages(idx) * m.lonProteaseEnergyCost;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(idx,m.compartment.cytosolIndexs) = 0;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade monomers
            m.evolveState_DegradeMonomers();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayOfMonomer_NoATP(this)
            m = this.process;
            
            idx = m.monomer.matureIndexs(1);
            energyIdxs = [m.substrateIndexs_atp;m.substrateIndexs_adp;m.substrateIndexs_phosphate;m.substrateIndexs_water;m.substrateIndexs_hydrogen];
            
            m.monomer.decayRates(:)=1e6;
            m.substrates = max(0,-m.monomerDecayReactions(:,idx));
            m.substrates(energyIdxs) = m.substrates(energyIdxs) + ...
                max([0 -1 -1 1 -1]',0) * m.monomerLonProteaseCleavages(idx) * m.lonProteaseEnergyCost;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lonProtease)=1e6;
            m.enzymes(m.enzymeIndexs_peptidases)=1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(idx,m.compartment.cytosolIndexs) =1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade monomers
            m.evolveState_DegradeMonomers();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayOfMonomer_NoLonProtease(this)
            m = this.process;
            
            idx = m.monomer.matureIndexs(1);
            energyIdxs = [m.substrateIndexs_atp;m.substrateIndexs_adp;m.substrateIndexs_phosphate;m.substrateIndexs_water;m.substrateIndexs_hydrogen];
            
            m.monomer.decayRates(:)=1e6;
            m.substrates = max(0,-m.monomerDecayReactions(:,idx));
            m.substrates(energyIdxs) = ...
                m.substrates(energyIdxs) + ...
                max([1 -1 -1 1 -1]',0) * m.monomerLonProteaseCleavages(idx) * m.lonProteaseEnergyCost;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lonProtease)=0;
            m.enzymes(m.enzymeIndexs_peptidases)=1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(idx,m.compartment.cytosolIndexs) =1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade monomers
            m.evolveState_DegradeMonomers();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayOfMonomer_NoPeptidase(this)
            m = this.process;
            
            idx = m.monomer.matureIndexs(1);
            energyIdxs = [m.substrateIndexs_atp;m.substrateIndexs_adp;m.substrateIndexs_phosphate;m.substrateIndexs_water;m.substrateIndexs_hydrogen];
            
            m.monomer.decayRates(:)=1e6;
            m.substrates = max(0,-m.monomerDecayReactions(:,idx));
            m.substrates(energyIdxs) = ...
                m.substrates(energyIdxs) + ...
                max([1 -1 -1 1 -1]',0) * m.monomerLonProteaseCleavages(idx) * m.lonProteaseEnergyCost;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lonProtease)=1e6;
            m.enzymes(m.enzymeIndexs_peptidases)=0;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(idx,m.compartment.cytosolIndexs) =1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade monomers
            m.evolveState_DegradeMonomers();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayOfMonomer_InTerminalOrganelleCytosol(this)
            m = this.process;
            
            idx = m.monomer.matureIndexs(1);
            energyIdxs = [m.substrateIndexs_atp;m.substrateIndexs_adp;m.substrateIndexs_phosphate;m.substrateIndexs_water;m.substrateIndexs_hydrogen];
            
            m.monomer.decayRates(:)=1e6;
            m.substrates = max(0,-m.monomerDecayReactions(:,idx));
            m.substrates(energyIdxs) = ...
                m.substrates(energyIdxs) + ...
                max([1 -1 -1 1 -1]',0) * m.monomerLonProteaseCleavages(idx) * m.lonProteaseEnergyCost;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lonProtease)=1e6;
            m.enzymes(m.enzymeIndexs_peptidases)=1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(idx,strcmp(m.compartment.wholeCellModelIDs,'tc')) =1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates + ...
                m.monomerDecayReactions(:,idx);
            final_substrates(energyIdxs) = ...
                final_substrates(energyIdxs) - ...
                [1 -1 -1 1 -1]' * m.monomerLonProteaseCleavages(idx) * m.lonProteaseEnergyCost;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(idx,strcmp(m.compartment.wholeCellModelIDs,'tc')) = 0;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade monomers
            m.evolveState_DegradeMonomers();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayOfMonomer_InMembrane(this)
            m = this.process;
                       
			idx = m.monomer.matureIndexs(find(m.monomer.compartments(m.monomer.matureIndexs) == m.compartment.membraneIndexs, 1, 'first'));
            
            energyIdxs = [m.substrateIndexs_atp; m.substrateIndexs_adp; m.substrateIndexs_phosphate; m.substrateIndexs_water; m.substrateIndexs_hydrogen];
            
            m.monomer.decayRates(m.monomer.compartments == m.compartment.cytosolIndexs) = 1e6; 
			m.monomer.decayRates(m.monomer.compartments == m.compartment.terminalOrganelleCytosolIndexs) = 1e6; 			
            m.substrates = max(0, -m.monomerDecayReactions(:, idx));
            m.substrates(energyIdxs) = ...
                m.substrates(energyIdxs) + ...
                max([1 -1 -1 1 -1]', 0) * m.monomerLonProteaseCleavages(idx) * m.lonProteaseEnergyCost;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lonProtease) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidases) = 1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(idx, m.compartment.membraneIndexs) = 1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade monomers
            m.evolveState_DegradeMonomers();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayOfMonomer_InTerminalOrganelleMembrane(this)
            m = this.process;
            
			idx = m.monomer.matureIndexs(find(m.monomer.compartments(m.monomer.matureIndexs) == m.compartment.membraneIndexs, 1, 'first'));
            
            energyIdxs = [m.substrateIndexs_atp; m.substrateIndexs_adp; m.substrateIndexs_phosphate; m.substrateIndexs_water; m.substrateIndexs_hydrogen];
            
            m.monomer.decayRates(m.monomer.compartments == m.compartment.cytosolIndexs) = 1e6; 
			m.monomer.decayRates(m.monomer.compartments == m.compartment.terminalOrganelleCytosolIndexs) = 1e6; 
            m.substrates = max(0, -m.monomerDecayReactions(:, idx));
            m.substrates(energyIdxs) = ...
                m.substrates(energyIdxs) + ...
                max([1 -1 -1 1 -1]', 0) * m.monomerLonProteaseCleavages(idx) * m.lonProteaseEnergyCost;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lonProtease) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidases) = 1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(idx, strcmp(m.compartment.wholeCellModelIDs, 'tm')) = 1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % degrade monomers
            m.evolveState_DegradeMonomers();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayMonomerAndTaggedMonomer(this)
            m = this.process;
            
            m.monomer.decayRates(:) = 1e6;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lonProtease) = 1e6;
            m.enzymes(m.enzymeIndexs_ftsHProtease) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidases) = 1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(m.monomer.matureIndexs(1), m.compartment.cytosolIndexs) = 1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = [1 30 0];
            
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(m.monomer.matureIndexs,m.compartment.cytosolIndexs) = 0;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = cell(0, 1);
            
            % degrade monomers
            m.evolveState_DegradeAbortedPolypeptides();
            m.evolveState_DegradeMonomers();
            
            % assert
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testDecayAllMonomersAndATaggedMonomer(this)
            m = this.process;
            
            m.monomer.decayRates(m.monomer.compartments == m.compartment.cytosolIndexs) = 1e6; 
			m.monomer.decayRates(m.monomer.compartments == m.compartment.terminalOrganelleCytosolIndexs) = 1e6; 
			
            m.substrates(:) = 1e10;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lonProtease) = 1e10;
            m.enzymes(m.enzymeIndexs_ftsHProtease) = 1e10;
            m.enzymes(m.enzymeIndexs_peptidases) = 1e10;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(m.monomer.compartments == m.compartment.cytosolIndexs, [m.compartment.cytosolIndexs; m.compartment.terminalOrganelleCytosolIndexs]) = 1;
			m.monomers(m.monomer.compartments == m.compartment.terminalOrganelleCytosolIndexs, [m.compartment.cytosolIndexs; m.compartment.terminalOrganelleCytosolIndexs]) = 1;
            m.complexs(:) = 0;
            m.monomer.chromosome.initialize();
            m.monomer.chromosome.monomerBoundSites([2000 * (1:numel(m.monomer.boundIndexs))' ones(size(m.monomer.boundIndexs))]) = 1:numel(m.monomer.boundIndexs);
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = [1 30 0];
            
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(:, [m.compartment.cytosolIndexs; m.compartment.terminalOrganelleCytosolIndexs]) = 0;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = cell(0, 1);
            
            % degrade monomers
            m.evolveState_DegradeAbortedPolypeptides();
            m.evolveState_DegradeMonomers();
            
            % assert
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testMonomerDecayFairness(this)
            m = this.process;
            
            m.monomer.decayRates(:) = 1e6;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lonProtease) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidases) = 1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            % degrade monomers
            counts = zeros(1,2);
            for i = 1:200
                m.substrates(:) = 1e6;
                m.substrates(m.substrateIndexs_atp) = ...
                    2 * m.lonProteaseEnergyCost * m.monomerLonProteaseCleavages(1) - 1;
                m.monomers(1, [m.compartment.cytosolIndexs; m.compartment.terminalOrganelleCytosolIndexs]) = 1;
                
                m.evolveState_DegradeMonomers();
                
                dCounts = m.monomers(1, [m.compartment.cytosolIndexs; m.compartment.terminalOrganelleCytosolIndexs]);
                assertEqual(1, sum(dCounts));
                counts = counts + dCounts;
            end
            assertTrue(range(counts) < 0.25 * max(counts));
        end
    end
    
    %comprehensive evolve state tests
    methods
        function testMonomerMisfoldingAndDegradation(this)
            m = this.process;
            
            idx = m.monomer.matureIndexs(1);
            energyIdxs = [
                m.substrateIndexs_atp;
                m.substrateIndexs_adp;
                m.substrateIndexs_phosphate;
                m.substrateIndexs_water;
                m.substrateIndexs_hydrogen];
            
            m.proteinMisfoldingRate = 1e6;
            m.monomer.decayRates(:) = Inf;
            m.monomer.decayRates(idx) = 0;
            m.substrates = max(0, -m.monomerDecayReactions(:, idx));
            m.substrates(energyIdxs) = ...
                m.substrates(energyIdxs) + ...
                max([1 -1 -1 1 -1]', 0) * m.monomerLonProteaseCleavages(idx) * m.lonProteaseEnergyCost; %#ok<*UDIM>
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lonProtease) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidases) = 1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.monomers(idx, m.compartment.cytosolIndexs) = 1;
            m.complexs(:) = 0;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_substrates = m.substrates + ...
                m.monomerDecayReactions(:, idx);
            final_substrates(energyIdxs) = ...
                final_substrates(energyIdxs) - ...
                [1 -1 -1 1 -1]' * m.monomerLonProteaseCleavages(idx) * m.lonProteaseEnergyCost;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_monomers(idx, m.compartment.cytosolIndexs) = 0;
            final_complexs = m.complexs;
            final_RNAs = m.RNAs;
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % misfold and degrade monomer
            m.evolveState();
            
            % assert
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testComplexMisfoldingAndCompleteDegradation(this)
            m = this.process;
            
            idx = m.complex.matureIndexs(10);
            
            m.proteinMisfoldingRate = 1e6;
            m.monomer.decayRates(:) = 1e6;
            m.complex.decayRates(:) = Inf;
            m.complex.decayRates(idx) = 0;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lonProtease) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidases) = 1e6;
            m.boundEnzymes(:) = 0;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.complexs(idx, m.compartment.cytosolIndexs) = 1;
            m.RNAs(:) = 0;
            m.polypeptide.abortedPolypeptides = zeros(0, 3);
            
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            final_monomers = m.monomers;
            final_complexs = m.complexs;
            final_complexs(idx, m.compartment.cytosolIndexs) = 0;
            final_RNAs = m.RNAs;
            final_RNAs(:, m.compartment.cytosolIndexs) = m.proteinComplexRNAComposition(:, idx);
            final_polypeptide.abortedSequences = m.polypeptide.abortedSequences;
            
            % misfold and degrade complexes
            m.evolveState();
            
            % assert
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(final_monomers, m.monomers);
            assertEqual(final_complexs, m.complexs);
            assertEqual(final_RNAs, m.RNAs);
            assertEqual(final_polypeptide.abortedSequences, m.polypeptide.abortedSequences);
        end
        
        function testExpectations(this)
            m = this.process;
            
            m.substrates(:) = 1e12;
            m.monomers(:) = 0;
            m.complexs(:) = 0;
            m.monomers(1) = 1e6;
            m.complexs(1) = 1e6;
            m.enzymes(:) = 1e12;
            
            monomers0 = m.monomers;
            complexs0 = m.complexs;
            
            iterMax = 10000;
            for i = 1:iterMax
                m.evolveState();
            end
            
            assertElementsAlmostEqual(monomers0(setdiff(1:end, m.monomer.damagedIndexs), :) .* exp(-m.monomer.decayRates(setdiff(1:end, m.monomer.damagedIndexs), ones(size(complexs0, 2), 1)) * iterMax), ...
                m.monomers(setdiff(1:end, m.monomer.damagedIndexs), :), 'relative', 1e-3);
            assertElementsAlmostEqual(complexs0 .* exp(-m.complex.decayRates(:, ones(size(complexs0, 2), 1)) * iterMax), ...
                m.complexs, 'relative', 1e-3);
        end
        
        function testSynthesisDecayBalance(this)
            m = this.process;
            pm = m.monomer;
            pc = m.complex;
            
            cellCycleLength = 9.0 * 3600;
            
            initMonomers = m.monomers;
            initMonomers(m.monomer.matureIndexs, :) = initMonomers(m.monomer.matureIndexs, :) + initMonomers(m.monomer.boundIndexs, :);
            initMonomers(m.monomer.boundIndexs, :) = 0;
            initComplexs = m.complexs;
            initComplexs(m.complex.matureIndexs, :) = initComplexs(m.complex.matureIndexs, :) + initComplexs(m.complex.boundIndexs, :);
            initComplexs(m.complex.boundIndexs, :) = 0;
            initEnzymes = m.enzymes;
            
            m.substrates(:) = 1e10;
            newMonomers = zeros(size(initMonomers));
            newComplexs = zeros(size(initComplexs));
            for i = 1:10000
                tmp = m.randStream.stochasticRound(...
                    + initMonomers * log(2) / cellCycleLength * exp(i * log(2) / cellCycleLength) ...
                    + initMonomers .* pm.decayRates(:, ones(1, size(initMonomers, 2))) * exp(i * log(2) / cellCycleLength) ...
                    );
                tmp(isnan(tmp)) = 0;
                newMonomers = newMonomers + tmp;
                m.monomers = m.monomers + tmp;
                
                tmp = m.randStream.stochasticRound(...
                    + initComplexs * log(2) / cellCycleLength * exp(i * log(2) / cellCycleLength) ...
                    + initComplexs .* pc.decayRates(:, ones(1, size(initComplexs, 2))) * exp(i * log(2) / cellCycleLength) ...
                    );
                tmp(isnan(tmp)) = 0;
                newComplexs = newComplexs + tmp;
                m.complexs = m.complexs + tmp;
                
                m.enzymes = initEnzymes * exp(i * log(2) / cellCycleLength);
                
                m.evolveState();
            end
            
            assertIn(max(max(m.monomers(pm.damagedIndexs, :))), [0 5]);
            assertIn(max(max(m.complexs(pc.damagedIndexs, :))), [0 0]);
            
            expNewMonomers = ...
                + initMonomers * (exp(i * log(2) / cellCycleLength) - 1) ...
                + initMonomers .* min(1, pm.decayRates(:, ones(1, size(initMonomers, 2)))) * cellCycleLength / log(2) * (exp(i * log(2) / cellCycleLength) - 1);
            expNewComplexs = ...
                + initComplexs * (exp(i * log(2) / cellCycleLength) - 1) ...
                + initComplexs .* min(1, pc.decayRates(:, ones(1, size(initComplexs, 2)))) * cellCycleLength / log(2) * (exp(i * log(2) / cellCycleLength) - 1);
            assertElementsAlmostEqual(sum(newMonomers(:)), sum(expNewMonomers(:)), 'relative', 10e-2);
            assertElementsAlmostEqual(sum(newComplexs(:)), sum(expNewComplexs(:)), 'relative', 10e-2);
            
            assertElementsAlmostEqual(sum(m.monomers' * pm.molecularWeights), ...
                sum((initMonomers * exp(i * log(2) / cellCycleLength))' * pm.molecularWeights), ...
                'relative', 5e-2);
            assertElementsAlmostEqual(sum(m.complexs' * pc.molecularWeights), ...
                sum((initComplexs * exp(i * log(2) / cellCycleLength))' * pc.molecularWeights), ...
                'relative', 5e-2);
        end
        
        function testGeneEssentiality(this)
            m = this.process;
            
            m.polypeptide.abortedPolypeptides = [1 30 0];
            
            m.proteinMisfoldingRate = 1e6;
            
            m.monomer.decayRates(m.monomer.matureIndexs(1))    = 1e6;
            m.monomer.decayRates(m.monomer.misfoldedIndexs(1)) = 1e6;
            m.monomer.decayRates(m.monomer.matureIndexs(2))    = 0;
            m.monomer.decayRates(m.monomer.misfoldedIndexs(2)) = 0;
            
            m.complex.decayRates(:)  = 1e6;
            m.complex.decayRates(m.complex.matureIndexs(1))    = 1e6;
            m.complex.decayRates(m.complex.misfoldedIndexs(1)) = 1e6;
            m.complex.decayRates(m.complex.matureIndexs(2))    = 0;
            m.complex.decayRates(m.complex.misfoldedIndexs(2)) = 0;
            
            m.enzymes(:) = 1e6;
            m.substrates(:) = 1e6;
            
            m.monomers(:)=0;
            m.complexs(:)=0;
            m.monomers(m.monomer.matureIndexs(1),    m.compartment.cytosolIndexs) = 1;
            m.monomers(m.monomer.misfoldedIndexs(2), m.compartment.cytosolIndexs) = 1;
            m.complexs(m.complex.matureIndexs(1),    m.compartment.cytosolIndexs) = 1;
            m.complexs(m.complex.misfoldedIndexs(2), m.compartment.cytosolIndexs) = 1;
            
            this.helpTestGeneEssentiality({
                'MG_020'; %proline iminopeptidase
                'MG_046'; %metalloendopeptidase
                'MG_183'; %oligoendopeptidase F
                'MG_208'; %glycoprotease
                'MG_239'; %protease La
                'MG_324'; %aminopeptidase
                'MG_355'; %protease clpB
                'MG_391'; %cytosol aminopeptidase
                'MG_457'; %metalloprotease FtsH
                }, @(m,~) ...
                    m.monomers(m.monomer.matureIndexs(2),    m.compartment.cytosolIndexs) == 1 && ...
                    m.monomers(m.monomer.matureIndexs(1),    m.compartment.cytosolIndexs) == 0 && ...
                    m.monomers(m.monomer.misfoldedIndexs(1), m.compartment.cytosolIndexs) == 0 && ...
                    m.complexs(m.complex.matureIndexs(2),    m.compartment.cytosolIndexs) == 1 && ...
                    m.complexs(m.complex.matureIndexs(1),    m.compartment.cytosolIndexs) == 0 && ...
                    m.complexs(m.complex.misfoldedIndexs(1), m.compartment.cytosolIndexs) == 0 && ...
                    isempty(m.polypeptide.abortedSequences));
        end
    end
end
