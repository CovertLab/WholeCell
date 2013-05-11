%RNA processing process test case
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef RNAProcessing_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    %constructor
    methods
        function this = RNAProcessing_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end
    
    %hard coded fixture
    methods
        function loadSimpleTestFixture(this)
            %% process
            m = this.process;
            
            %% constants
            m.enzymeSpecificRate_DeaD     = 1.48;
            m.enzymeSpecificRate_RsgA     = 0.2917;
            m.enzymeSpecificRate_RNAseP   = 6;
            m.enzymeSpecificRate_RNAseIII = 7.7;
            m.enzymeSpecificRate_RNAseJ   = 0.37;
            
            m.enzymeEnergyCost_DeaD       = 2.37;
            m.enzymeEnergyCost_RsgA       = 1;
            
            %% substrates
            m.substrateWholeCellModelIDs = { %whole cell model IDs of substrates
                'ATP';'GTP';
                'ADP';'GDP';
                'PI';'H2O';'H'};
            m.substrateNames = m.substrateWholeCellModelIDs;
            
            m.substrateIndexs_NTPs  = (1:2)';
            m.substrateIndexs_ATP   = 1;
            m.substrateIndexs_GTP   = 2;
            m.substrateIndexs_NDPs  = (3:4)';
            m.substrateIndexs_ADP   = 3;
            m.substrateIndexs_GDP   = 4;
            m.substrateIndexs_PI    = 5;
            m.substrateIndexs_H2O   = 6;
            m.substrateIndexs_H     = 7;
            
            m.substrateMolecularWeights = [503.1000; 519.1000; 424.2000; 440.2000; 96.0000; 18.0000; 1.0000];
            
            %% enzymes
            m.enzymeWholeCellModelIDs={
                'MG_367_DIMER';        %ribonuclease III
                'MG_139_DIMER';        %ribonuclease J
                'MG_0003_465';         %ribonuclease P
                'MG_110_MONOMER';      %ribosome small subunit-dependent GTPase A
                'MG_425_DIMER'};       %ATP-dependent RNA helicase
            m.enzymeNames    = m.enzymeWholeCellModelIDs;
            
            m.enzymeIndexs_RNAseIII = 1; %index within enzymes of ribonuclease III
            m.enzymeIndexs_RNAseJ   = 2; %index within enzymes of ribonuclease J
            m.enzymeIndexs_RNAseP   = 3; %index within enzymes of ribonuclease P
            m.enzymeIndexs_RsgA     = 4; %index within enzymes of ribosome small subunit-dependent GTPase A
            m.enzymeIndexs_DeaD     = 5; %index within enzymes of ATP-dependent RNA helicase
            
            m.enzymeMolecularWeights = 1e5 * [0.6067; 1.2846; 1.3499; 0.3218; 1.0256];
            
            %% RNAs
            %unprocessed
            m.unprocessedRNAWholeCellModelIDs = {...
                'TU_mRNA_001'; 'TU_mRNA_002'; 'TU_rRNA_001'; 'TU_sRNA_001'; 'TU_sRNA_002'; 'TU_tRNA_001'; 'TU_tRNA_002'};
            
            m.unprocessedRNAIndexs_mRNA = (1:2)';
            m.unprocessedRNAIndexs_rRNA = 3;
            m.unprocessedRNAIndexs_sRNA = (4:5);
            m.unprocessedRNAIndexs_tRNA = (6:7)';
            m.unprocessedRNAIndexs_scRNA = 4;
            m.unprocessedRNAIndexs_tmRNA = 5;
            
            %processed
            m.processedRNAWholeCellModelIDs = {...
                'TU_mRNA_001';
                'TU_mRNA_002';
                'MGrrnA16S'; 'MGrrnA23S'; 'MGrrnA5S';
                'MG_0001'; 'MG_0004';
                'MG508'; 'MG509'; 'MG510'; 'MG511'; 'MG512'; 'MG513'; 'MG514';
                'MG492'};
            
            m.processedRNAIndexs_mRNA = (1:2)';
            m.processedRNAIndexs_rRNA = (3:5)';
            m.processedRNAIndexs_sRNA = (6:7);
            m.processedRNAIndexs_tRNA = (8:15)';
            m.processedRNAIndexs_scRNA = 6;
            m.processedRNAIndexs_tmRNA = 7;
            
            %intergenic
            m.intergenicRNAWholeCellModelIDs = {...
                'TU_rRNA_001-01';
                'TU_rRNA_001-02';
                'TU_tRNA_001-01';
                'TU_tRNA_001-02';
                'TU_tRNA_001-03';
                'TU_tRNA_001-04';
                'TU_tRNA_001-05';
                'TU_tRNA_001-06'};
            
            %% adjacency matrices
            
            %processed RNAs X unprocessed RNAs
            m.rna.nascentRNAMatureRNAComposition = zeros(numel(m.processedRNAWholeCellModelIDs), numel(m.unprocessedRNAWholeCellModelIDs));
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'TU_mRNA_001'), strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_mRNA_001')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'TU_mRNA_002'), strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_mRNA_002')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'MGrrnA16S'),   strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_rRNA_001')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'MGrrnA23S'),   strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_rRNA_001')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'MGrrnA5S'),    strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_rRNA_001')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'MG_0001'),     strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_sRNA_001')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'MG508'),       strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'MG509'),       strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'MG510'),       strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'MG511'),       strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'MG512'),       strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'MG513'),       strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'MG514'),       strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs, 'MG492'),       strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_002')) = 1;
            
            %intergenic RNAs X unprocessed RNAs
            m.rna.intergenicRNAMatrix = zeros(numel(m.intergenicRNAWholeCellModelIDs), numel(m.unprocessedRNAWholeCellModelIDs));
            m.rna.intergenicRNAMatrix(strcmp(m.intergenicRNAWholeCellModelIDs, 'TU_rRNA_001-01'), strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_rRNA_001')) = 1;
            m.rna.intergenicRNAMatrix(strcmp(m.intergenicRNAWholeCellModelIDs, 'TU_rRNA_001-02'), strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_rRNA_001')) = 1;
            m.rna.intergenicRNAMatrix(strcmp(m.intergenicRNAWholeCellModelIDs, 'TU_tRNA_001-01'), strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            m.rna.intergenicRNAMatrix(strcmp(m.intergenicRNAWholeCellModelIDs, 'TU_tRNA_001-02'), strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            m.rna.intergenicRNAMatrix(strcmp(m.intergenicRNAWholeCellModelIDs, 'TU_tRNA_001-03'), strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            m.rna.intergenicRNAMatrix(strcmp(m.intergenicRNAWholeCellModelIDs, 'TU_tRNA_001-04'), strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            m.rna.intergenicRNAMatrix(strcmp(m.intergenicRNAWholeCellModelIDs, 'TU_tRNA_001-05'), strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            m.rna.intergenicRNAMatrix(strcmp(m.intergenicRNAWholeCellModelIDs, 'TU_tRNA_001-06'), strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            
            %reactants and byproducts X unprocessed RNAs
            m.reactantByproductMatrix = zeros(numel(m.substrateWholeCellModelIDs), numel(m.unprocessedRNAWholeCellModelIDs));
            
            m.reactantByproductMatrix(m.substrateIndexs_H2O, :) = -(sum(m.rna.nascentRNAMatureRNAComposition,1) + sum(m.rna.intergenicRNAMatrix,1) - 1);
            m.reactantByproductMatrix(m.substrateIndexs_H, :)   =  (sum(m.rna.nascentRNAMatureRNAComposition,1) + sum(m.rna.intergenicRNAMatrix,1) - 1);
            
            energyIdxs_DeaD = [m.substrateIndexs_ATP; m.substrateIndexs_ADP; m.substrateIndexs_PI; m.substrateIndexs_H2O; m.substrateIndexs_H];
            m.reactantByproductMatrix(energyIdxs_DeaD, m.unprocessedRNAIndexs_rRNA) = ...
                m.reactantByproductMatrix(energyIdxs_DeaD, m.unprocessedRNAIndexs_rRNA) + ...
                round(m.enzymeEnergyCost_RsgA) * [-1; 1; 1; -1; 1];
            
            energyIdxs_RsgA = [m.substrateIndexs_GTP; m.substrateIndexs_GDP; m.substrateIndexs_PI; m.substrateIndexs_H2O; m.substrateIndexs_H];
            m.reactantByproductMatrix(energyIdxs_RsgA, m.unprocessedRNAIndexs_rRNA) = ...
                m.reactantByproductMatrix(energyIdxs_RsgA, m.unprocessedRNAIndexs_rRNA) + ...
                round(m.enzymeEnergyCost_DeaD) * [-1; 1; 1; -1; 1];
            
            %enzymes X unprocessed RNAs
            catalysisMatrix = zeros(numel(m.enzymeWholeCellModelIDs), numel(m.unprocessedRNAWholeCellModelIDs));
            catalysisMatrix(m.enzymeIndexs_RNAseIII, strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_rRNA_001')) = 2;
            catalysisMatrix(m.enzymeIndexs_RNAseJ,   strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_rRNA_001')) = 2;
            catalysisMatrix(m.enzymeIndexs_RsgA,     strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_rRNA_001')) = 1;
            catalysisMatrix(m.enzymeIndexs_DeaD,     strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_rRNA_001')) = 1;
            catalysisMatrix(m.enzymeIndexs_RNAseIII, strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_sRNA_001')) = 1;
            catalysisMatrix(m.enzymeIndexs_RNAseIII, strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            catalysisMatrix(m.enzymeIndexs_RNAseP,   strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_001')) = 1;
            catalysisMatrix(m.enzymeIndexs_RNAseIII, strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_002')) = 1;
            catalysisMatrix(m.enzymeIndexs_RNAseP,   strcmp(m.unprocessedRNAWholeCellModelIDs, 'TU_tRNA_002')) = 1;
            
            enzymeSpecificRates = zeros(5,1);
            enzymeSpecificRates(m.enzymeIndexs_DeaD)     = m.enzymeSpecificRate_DeaD;
            enzymeSpecificRates(m.enzymeIndexs_RsgA)     = m.enzymeSpecificRate_RsgA;
            enzymeSpecificRates(m.enzymeIndexs_RNAseP)   = m.enzymeSpecificRate_RNAseP;
            enzymeSpecificRates(m.enzymeIndexs_RNAseIII) = m.enzymeSpecificRate_RNAseIII;
            enzymeSpecificRates(m.enzymeIndexs_RNAseJ)   = m.enzymeSpecificRate_RNAseJ;
            
            m.catalysisMatrix = catalysisMatrix ./ repmat(enzymeSpecificRates, 1, numel(m.unprocessedRNAWholeCellModelIDs));
            
            %% RNA molecular weights
            m.rna.molecularWeights(m.rna.processedIndexs) = 1e6 + (1:numel(m.processedRNAWholeCellModelIDs))';
            m.rna.molecularWeights(m.rna.intergenicIndexs) = 1e5 + (1:numel(m.intergenicRNAWholeCellModelIDs))';
            m.rna.molecularWeights(m.rna.nascentIndexs) = ...
                m.rna.nascentRNAMatureRNAComposition' * m.rna.molecularWeights(m.rna.processedIndexs)  + ...
                m.rna.intergenicRNAMatrix' * m.rna.molecularWeights(m.rna.intergenicIndexs) + ...
                - (sum(m.rna.nascentRNAMatureRNAComposition,1) + sum(m.rna.intergenicRNAMatrix,1) - 1)' * ...
                (m.substrateMolecularWeights(m.substrateIndexs_H2O) - m.substrateMolecularWeights(m.substrateIndexs_H));
            
            %% initial state
            m.substrates      = zeros(length(m.substrateWholeCellModelIDs),      1);
            m.enzymes         = zeros(length(m.enzymeWholeCellModelIDs),         1);
            m.boundEnzymes    = zeros(length(m.enzymeWholeCellModelIDs),         1);
            m.unprocessedRNAs = zeros(length(m.unprocessedRNAWholeCellModelIDs), 1);
            m.processedRNAs   = zeros(length(m.processedRNAWholeCellModelIDs),   1);
            m.intergenicRNAs  = zeros(length(m.intergenicRNAWholeCellModelIDs),  1);
        end
    end
    
    %tests
    methods
        function testConstants(this)
            %% process
            m = this.process;
            
            %% mRNA
            assertEqual([0;1], unique(m.rna.nascentRNAMatureRNAComposition(:,m.unprocessedRNAIndexs_mRNA)));
            assertEqual(numel(m.unprocessedRNAIndexs_mRNA), sum(sum(m.rna.nascentRNAMatureRNAComposition(:,m.unprocessedRNAIndexs_mRNA))));
            assertEqual(0, unique(m.rna.intergenicRNAMatrix(:,m.unprocessedRNAIndexs_mRNA)));
            
            assertEqual(0, unique(m.reactantByproductMatrix(:,m.unprocessedRNAIndexs_mRNA)));
            
            assertEqual(0, unique(m.catalysisMatrix(:, m.unprocessedRNAIndexs_mRNA)));
            
            %% rRNA transcription unit
            rRNATranscriptionUnitIdx = find(m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs,'MGrrnA16S'),:));
            
            assertEqual({'MGrrnA16S';'MGrrnA23S';'MGrrnA5S'}, m.processedRNAWholeCellModelIDs(m.rna.nascentRNAMatureRNAComposition(:,rRNATranscriptionUnitIdx)>0));
            assertEqual(2, sum(m.rna.intergenicRNAMatrix(:,rRNATranscriptionUnitIdx)));
            
            assertEqual(23, sum(abs(m.reactantByproductMatrix(:,rRNATranscriptionUnitIdx))));
            assertEqual(-1, m.reactantByproductMatrix(m.substrateIndexs_ATP, rRNATranscriptionUnitIdx));
            assertEqual(-2, m.reactantByproductMatrix(m.substrateIndexs_GTP, rRNATranscriptionUnitIdx));
            assertEqual(1,  m.reactantByproductMatrix(m.substrateIndexs_ADP, rRNATranscriptionUnitIdx));
            assertEqual(2,  m.reactantByproductMatrix(m.substrateIndexs_GDP, rRNATranscriptionUnitIdx));
            assertEqual(3,  m.reactantByproductMatrix(m.substrateIndexs_PI,  rRNATranscriptionUnitIdx));
            assertEqual(7,  m.reactantByproductMatrix(m.substrateIndexs_H,   rRNATranscriptionUnitIdx));
            assertEqual(-7, m.reactantByproductMatrix(m.substrateIndexs_H2O, rRNATranscriptionUnitIdx));
            
            assertEqual(4, numel(find(m.catalysisMatrix(:, rRNATranscriptionUnitIdx))));
            assertEqual(2 / m.enzymeSpecificRate_RNAseIII, m.catalysisMatrix(m.enzymeIndexs_RNAseIII, rRNATranscriptionUnitIdx));
            assertEqual(2 / m.enzymeSpecificRate_RNAseJ,   m.catalysisMatrix(m.enzymeIndexs_RNAseJ,   rRNATranscriptionUnitIdx));
            assertEqual(1 / m.enzymeSpecificRate_RsgA,     m.catalysisMatrix(m.enzymeIndexs_RsgA,     rRNATranscriptionUnitIdx));
            assertEqual(1 / m.enzymeSpecificRate_DeaD,     m.catalysisMatrix(m.enzymeIndexs_DeaD,     rRNATranscriptionUnitIdx));
            
            %% sRNA
            assertEqual([0;1], unique(m.rna.nascentRNAMatureRNAComposition(:,m.unprocessedRNAIndexs_sRNA)));
            assertEqual(numel(m.unprocessedRNAIndexs_sRNA), sum(sum(m.rna.nascentRNAMatureRNAComposition(:,m.unprocessedRNAIndexs_sRNA))));
            assertEqual(0, unique(m.rna.intergenicRNAMatrix(:,m.unprocessedRNAIndexs_sRNA)));
            
            assertEqual(0, unique(m.reactantByproductMatrix(:,m.unprocessedRNAIndexs_sRNA)));

            assertEqual(0, unique(m.catalysisMatrix(:, setdiff(m.unprocessedRNAIndexs_sRNA, [m.unprocessedRNAIndexs_scRNA; m.unprocessedRNAIndexs_tmRNA]))));
            
            assertEqual(0, unique(m.catalysisMatrix(setdiff(1:end, m.enzymeIndexs_RNAseIII), m.unprocessedRNAIndexs_scRNA)));
            assertEqual(2 / m.enzymeSpecificRate_RNAseIII, m.catalysisMatrix(m.enzymeIndexs_RNAseIII, m.unprocessedRNAIndexs_scRNA));
            
            assertEqual(0, unique(m.catalysisMatrix(setdiff(1:end, m.enzymeIndexs_RNAseP), m.unprocessedRNAIndexs_tmRNA)));
            assertEqual(1 / m.enzymeSpecificRate_RNAseP, m.catalysisMatrix(m.enzymeIndexs_RNAseP, m.unprocessedRNAIndexs_tmRNA));
            
            %% tRNA
            tRNATranscriptionUnitIdx = find(m.rna.nascentRNAMatureRNAComposition(strcmp(m.processedRNAWholeCellModelIDs,'MG508'),:));
            
            assertEqual({'MG508'; 'MG509'; 'MG510'; 'MG511'; 'MG512'; 'MG513'; 'MG514'; }, m.processedRNAWholeCellModelIDs(m.rna.nascentRNAMatureRNAComposition(:,tRNATranscriptionUnitIdx)>0));
            assertEqual(6, sum(m.rna.intergenicRNAMatrix(:,tRNATranscriptionUnitIdx)));
            
            assertEqual(24, sum(abs(m.reactantByproductMatrix(:,tRNATranscriptionUnitIdx))));
            assertEqual(12,  m.reactantByproductMatrix(m.substrateIndexs_H,   tRNATranscriptionUnitIdx));
            assertEqual(-12, m.reactantByproductMatrix(m.substrateIndexs_H2O, tRNATranscriptionUnitIdx));
            
            assertEqual(2, numel(find(m.catalysisMatrix(:, tRNATranscriptionUnitIdx))));
            assertEqual(1 / m.enzymeSpecificRate_RNAseIII, m.catalysisMatrix(m.enzymeIndexs_RNAseIII, tRNATranscriptionUnitIdx));
            assertEqual(1 / m.enzymeSpecificRate_RNAseP,   m.catalysisMatrix(m.enzymeIndexs_RNAseP,   tRNATranscriptionUnitIdx));
            
            %% molecular weights
            assertTrue(all(m.rna.molecularWeights > 0));
            
            assertTrue(all(m.rna.nascentRNAMatureRNAComposition  * m.rna.molecularWeights(m.rna.nascentIndexs) >= m.rna.molecularWeights(m.rna.processedIndexs)));
            assertTrue(all(m.rna.intergenicRNAMatrix * m.rna.molecularWeights(m.rna.nascentIndexs) >= m.rna.molecularWeights(m.rna.intergenicIndexs)));
            
            assertElementsAlmostEqual(m.rna.molecularWeights(m.rna.nascentIndexs), ...
                m.rna.nascentRNAMatureRNAComposition' * m.rna.molecularWeights(m.rna.processedIndexs)  + ...
                m.rna.intergenicRNAMatrix' * m.rna.molecularWeights(m.rna.intergenicIndexs) + ...
                - (sum(m.rna.nascentRNAMatureRNAComposition,1) + sum(m.rna.intergenicRNAMatrix,1) - 1)' * ...
                (m.substrateMolecularWeights(m.substrateIndexs_H2O) - m.substrateMolecularWeights(m.substrateIndexs_H)));
            
             assertElementsAlmostEqual(m.rna.molecularWeights(m.rna.nascentIndexs), ...
                m.rna.nascentRNAMatureRNAComposition' * m.rna.molecularWeights(m.rna.processedIndexs)  + ...
                m.rna.intergenicRNAMatrix' * m.rna.molecularWeights(m.rna.intergenicIndexs) + ...
                m.reactantByproductMatrix' * m.substrateMolecularWeights);
            
            %% adjacency matrices
            assertTrue(all(sum(m.rna.nascentRNAMatureRNAComposition,1)));
            assertEqual([0;1], unique(m.rna.nascentRNAMatureRNAComposition));
            assertEqual([0;1], unique(m.rna.intergenicRNAMatrix));
            assertEqual(ones(size(m.processedRNAs)), sum(m.rna.nascentRNAMatureRNAComposition,2));
            assertEqual(ones(size(m.intergenicRNAs)), sum(m.rna.intergenicRNAMatrix,2));
            
            notWaterIdxs = setdiff(1:numel(m.substrateWholeCellModelIDs), [m.substrateIndexs_H2O; m.substrateIndexs_H]);
            assertEqual(1, size(unique(m.reactantByproductMatrix(notWaterIdxs, m.unprocessedRNAIndexs_mRNA)','rows'),1));
            assertEqual(1, size(unique(m.reactantByproductMatrix(notWaterIdxs, m.unprocessedRNAIndexs_rRNA)','rows'),1));
            assertEqual(1, size(unique(m.reactantByproductMatrix(notWaterIdxs, m.unprocessedRNAIndexs_sRNA)','rows'),1));
            assertEqual(1, size(unique(m.reactantByproductMatrix(notWaterIdxs, m.unprocessedRNAIndexs_tRNA)','rows'),1));
            assertEqual(1, size(unique(m.catalysisMatrix(:,m.unprocessedRNAIndexs_mRNA)','rows'),1));
            assertEqual(1, size(unique(m.catalysisMatrix(:,m.unprocessedRNAIndexs_rRNA)','rows'),1));
            assertEqual(1, size(unique(m.catalysisMatrix(:,m.unprocessedRNAIndexs_sRNA)','rows'),3));
            assertEqual(1, size(unique(m.catalysisMatrix(:,m.unprocessedRNAIndexs_tRNA)','rows'),1));
        end
        
        function testMRNAProcessing(this)
            m = this.process;
            
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.unprocessedRNAs(m.unprocessedRNAIndexs_mRNA) = 1;
            m.processedRNAs(:)   = 0;
            m.intergenicRNAs(:)   = 0;
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            
            m.evolveState();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(0, unique(m.unprocessedRNAs));
            assertEqual(numel(m.processedRNAIndexs_mRNA), sum(abs(m.processedRNAs)));
            assertEqual(1, unique(m.processedRNAs(m.processedRNAIndexs_mRNA)));
            assertEqual(0, sum(abs([])));
        end
        
        function testRRNAProcessing(this)
            m = this.process;
            
            %with water, hydrogen, and enzymes
            m.substrates = sum(max(0, -m.reactantByproductMatrix(:, m.unprocessedRNAIndexs_rRNA)),2);
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_RNAseIII) = 1e6;
            m.enzymes(m.enzymeIndexs_RNAseJ) = 1e6;
            m.enzymes(m.enzymeIndexs_RsgA) = 1e6;
            m.enzymes(m.enzymeIndexs_DeaD) = 1e6;
            m.boundEnzymes(:) = 0;
            m.unprocessedRNAs(m.unprocessedRNAIndexs_rRNA) = 1;
            m.processedRNAs(:)   = 0;
            m.intergenicRNAs(:)   = 0;
            
            final_substrates = sum(max(0, m.reactantByproductMatrix(:, m.unprocessedRNAIndexs_rRNA)),2);
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            
            m.evolveState();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(0, unique(m.unprocessedRNAs));
            assertEqual(numel(m.processedRNAIndexs_rRNA), sum(abs(m.processedRNAs)));
            assertEqual(1, unique(m.processedRNAs(m.processedRNAIndexs_rRNA)));
            if ~isempty(m.intergenicRNAs)
                assertEqual([0;1], unique(m.intergenicRNAs));
                assertEqual(2, sum(abs(m.intergenicRNAs)));
                assertEqual(sum(m.rna.intergenicRNAMatrix(:, m.unprocessedRNAIndexs_rRNA),2), m.intergenicRNAs);
            end
            
            %no RsgA
            m.substrates = sum(max(0, -m.reactantByproductMatrix(:, m.unprocessedRNAIndexs_rRNA)),2);
            m.enzymes(m.enzymeIndexs_RsgA) = 0;
            m.unprocessedRNAs(m.unprocessedRNAIndexs_rRNA) = 1;
            m.processedRNAs(:)   = 0;
            m.intergenicRNAs(:)   = 0;
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(numel(m.unprocessedRNAIndexs_rRNA), sum(abs(m.unprocessedRNAs)));
            assertEqual(1, unique(m.unprocessedRNAs(m.unprocessedRNAIndexs_rRNA)));
            assertEqual(0, unique(m.processedRNAs));
            assertEqual(0, sum(abs((m.intergenicRNAs))));
            
            %no ATP
            m.substrates = sum(max(0, -m.reactantByproductMatrix(:, m.unprocessedRNAIndexs_rRNA)),2);
            m.substrates(m.substrateIndexs_ATP) = 0;
            m.enzymes(m.enzymeIndexs_RsgA) = 1e6;
            m.unprocessedRNAs(m.unprocessedRNAIndexs_rRNA) = 1;
            m.processedRNAs(:)   = 0;
            m.intergenicRNAs(:)   = 0;
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(numel(m.unprocessedRNAIndexs_rRNA), sum(abs(m.unprocessedRNAs)));
            assertEqual(1, unique(m.unprocessedRNAs(m.unprocessedRNAIndexs_rRNA)));
            assertEqual(0, unique(m.processedRNAs));
            assertEqual(0, sum(abs((m.intergenicRNAs))));
        end
        
        function testSRNAProcessing(this)
            m = this.process;
            
            %with RNAseIII
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_RNAseIII) = 1e6;
            m.boundEnzymes(:) = 0;
            m.unprocessedRNAs(m.unprocessedRNAIndexs_scRNA) = 1;
            m.processedRNAs(:)   = 0;
            m.intergenicRNAs(:)   = 0;
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            
            m.evolveState();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(0, unique(m.unprocessedRNAs));
            assertEqual(numel(m.processedRNAIndexs_scRNA), sum(abs(m.processedRNAs)));
            assertEqual(1, unique(m.processedRNAs(m.processedRNAIndexs_scRNA)));
            assertEqual(0, sum(abs((m.intergenicRNAs))));
            
            %without RNAseIII
            m.enzymes(m.enzymeIndexs_RNAseIII) = 0;
            m.unprocessedRNAs(m.unprocessedRNAIndexs_scRNA) = 1;
            m.processedRNAs(:)   = 0;
            
            final_enzymes = m.enzymes;
            
            m.evolveState();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(numel(m.processedRNAIndexs_scRNA), sum(abs(m.unprocessedRNAs)));
            assertEqual(1, unique(m.unprocessedRNAs(m.unprocessedRNAIndexs_scRNA)));
            assertEqual(0, unique(m.processedRNAs));
            assertEqual(0, sum(abs((m.intergenicRNAs))));
        end
        
        function testTRNAProcessing(this)
            m = this.process;
            
            %with water, hydrogen, RNAseIII, and RNAseP
            m.substrates = sum(max(0, -m.reactantByproductMatrix(:, m.unprocessedRNAIndexs_tRNA)),2);
            m.substrates(m.substrateIndexs_H2O)=1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_RNAseIII) = 1e6;
            m.enzymes(m.enzymeIndexs_RNAseP) = 1e6;
            m.boundEnzymes(:) = 0;
            m.unprocessedRNAs(m.unprocessedRNAIndexs_tRNA) = 1;
            m.processedRNAs(:)   = 0;
            m.intergenicRNAs(:)   = 0;
            
            final_substrates = m.substrates + sum(m.reactantByproductMatrix(:, m.unprocessedRNAIndexs_tRNA),2);
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            
            m.evolveState();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(0, unique(m.unprocessedRNAs));
            assertEqual(numel(m.processedRNAIndexs_tRNA), sum(abs(m.processedRNAs)));
            assertEqual(1, unique(m.processedRNAs(m.processedRNAIndexs_tRNA)));
            assertEqual(sum(m.rna.intergenicRNAMatrix(:, m.unprocessedRNAIndexs_tRNA),2), m.intergenicRNAs);
            
            %without RNAseP
            m.substrates = sum(max(0, -m.reactantByproductMatrix(:, m.unprocessedRNAIndexs_tRNA)),2);
            m.substrates(m.substrateIndexs_H2O)=1e6;
            m.enzymes(m.enzymeIndexs_RNAseP) = 0;
            m.unprocessedRNAs(m.unprocessedRNAIndexs_tRNA) = 1;
            m.processedRNAs(:)   = 0;
            m.intergenicRNAs(:)   = 0;
            
            final_substrates = m.substrates;
            final_enzymes = m.enzymes;
            
            m.evolveState();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(numel(m.unprocessedRNAIndexs_tRNA), sum(abs(m.unprocessedRNAs)));
            assertEqual(1, unique(m.unprocessedRNAs(m.unprocessedRNAIndexs_tRNA)));
            assertEqual(0, unique(m.processedRNAs));
            assertEqual(0, sum(abs((m.intergenicRNAs))));           
        end
        
        function testTRNAProcessingWaterLimits(this)
            m = this.process;
            
            %without water
            if any(m.reactantByproductMatrix(m.substrateIndexs_H2O,  m.unprocessedRNAIndexs_tRNA))
                m.substrates = sum(max(0, -m.reactantByproductMatrix(:, m.unprocessedRNAIndexs_tRNA)),2);
                m.substrates(m.substrateIndexs_H2O)=1e6;
                m.substrates(m.substrateIndexs_H2O) = 0;
                m.enzymes(m.enzymeIndexs_RNAseP) = 1e6;
                m.unprocessedRNAs(m.unprocessedRNAIndexs_tRNA) = 1;
                m.processedRNAs(:) = 0;
                m.intergenicRNAs(:) = 0;
                
                final_substrates = m.substrates;
                final_enzymes = m.enzymes;
                final_boundEnzymes = m.boundEnzymes;
                
                m.evolveState();
                
                assertEqual(final_substrates, m.substrates);
                assertEqual(final_enzymes, m.enzymes);
                assertEqual(final_boundEnzymes, m.boundEnzymes);
                assertEqual(numel(m.unprocessedRNAIndexs_tRNA), sum(abs(m.unprocessedRNAs)));
                assertEqual(1, unique(m.unprocessedRNAs(m.unprocessedRNAIndexs_tRNA)));
                assertEqual(0, unique(m.processedRNAs));
                assertEqual(0, sum(abs((m.intergenicRNAs))));
            end
        end
        
        function testAllRNAProcessing(this)
            m = this.process;
            
            m.substrates = sum(max(0, -m.reactantByproductMatrix),2);
            m.substrates(m.substrateIndexs_H2O)=1e6;
            m.enzymes(:) = 1e6;
            m.boundEnzymes(:) = 0;
            m.unprocessedRNAs(:) = 1;
            m.processedRNAs(:)   = 0;
            m.intergenicRNAs(:)   = 0;
            
            final_substrates =  m.substrates + sum(m.reactantByproductMatrix,2);
            final_enzymes = m.enzymes;
            final_boundEnzymes = m.boundEnzymes;
            
            m.evolveState();
            
            assertEqual(final_substrates, m.substrates);
            assertEqual(final_enzymes, m.enzymes);
            assertEqual(final_boundEnzymes, m.boundEnzymes);
            assertEqual(0, unique(m.unprocessedRNAs));
            assertEqual(1, unique(m.processedRNAs));
            assertEqual((numel(m.intergenicRNAs)>0) + 0, numel(unique(m.intergenicRNAs)));
            assertEqual(numel(m.intergenicRNAs), sum(abs(m.intergenicRNAs)));
        end
        
        function testFairTRNAProcessing(this)
            m = this.process;
            
            idx1 = m.unprocessedRNAIndexs_tRNA(1);
            idx2 = m.unprocessedRNAIndexs_tRNA(2);
            
            m.rna.nascentRNAMatureRNAComposition(:,m.unprocessedRNAIndexs_tRNA)      = repmat(m.rna.nascentRNAMatureRNAComposition(:,idx1),      1, numel(m.unprocessedRNAIndexs_tRNA));
            m.rna.intergenicRNAMatrix(:,m.unprocessedRNAIndexs_tRNA)     = repmat(m.rna.intergenicRNAMatrix(:,idx1),     1, numel(m.unprocessedRNAIndexs_tRNA));
            m.reactantByproductMatrix(:,m.unprocessedRNAIndexs_tRNA) = repmat(m.reactantByproductMatrix(:,idx1), 1, numel(m.unprocessedRNAIndexs_tRNA));
            m.catalysisMatrix(:,m.unprocessedRNAIndexs_tRNA)         = repmat(m.catalysisMatrix(:,idx1),         1, numel(m.unprocessedRNAIndexs_tRNA));
            
            counts = repmat(200, 2,1);
            for i = 1:50
                m.unprocessedRNAs(:)=0;
                m.unprocessedRNAs([idx1 idx2]) = 1;
                m.processedRNAs(:) = 0;
                m.intergenicRNAs(:) = 0;
                m.substrates = max(0, -m.reactantByproductMatrix(:,idx1));
                m.enzymes(:)=1e6;
                m.enzymes(m.enzymeIndexs_RNAseIII) = 1/m.enzymeSpecificRate_RNAseIII;
                
                m.evolveState();
                
                assertEqual(1, sum(m.unprocessedRNAs([idx1 idx2])));
                counts = counts - m.unprocessedRNAs([idx1 idx2]);
            end
            assertTrue(range(counts) < 0.1 * max(counts));
        end
        
        function testGeneEssentiality(this)
            m = this.process;
            
            m.enzymes(:) = 1e6;
            m.boundEnzymes(:) = 0;
            m.unprocessedRNAs(:) = 1;
            m.processedRNAs(:)   = 0;
            m.intergenicRNAs(:)   = 0;
            m.substrates(:) = 1e6;
            
            this.helpTestGeneEssentiality({
                'MG_0003'; %ribonuclease P
                'MG_110';  %ribosome small subunit-dependent GTPase A
                'MG_139';  %ribonuclease J
                'MG_367';  %ribonuclease III
                'MG_425';  %ATP-dependent RNA helicase});
                'MG_465';  %ribonuclease P
                }, @this.isProperlyFunctioning);
        end
    end
    
    %helper methods
    methods
        function result = isProperlyFunctioning(~, m, i)
            result = ...
                all(m.unprocessedRNAs < i.unprocessedRNAs | i.unprocessedRNAs==0) & ...
                all(m.processedRNAs > i.processedRNAs | (m.rna.nascentRNAMatureRNAComposition * i.unprocessedRNAs)==0);
        end
    end
end
