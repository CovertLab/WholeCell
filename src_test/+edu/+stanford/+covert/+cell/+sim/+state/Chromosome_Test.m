%Chromosome_Test
% Chromosome test class.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef Chromosome_Test < edu.stanford.covert.cell.sim.CellStateTestCase
    methods
        function this = Chromosome_Test(name)
            this = this@edu.stanford.covert.cell.sim.CellStateTestCase(name);
        end
    end
    
    %fixtures
    methods
        function loadSimpleTestFixture(this)
            s = this.state;
            s.invalidate();
            
            %% constants
            s.strandBreakClassification_doubleStrandBreakSeparation = 10;
            s.strandBreakClassification_segmentLength               = 216;
            
            s.relaxedBasesPerTurn                                   = 10.5;
            s.equilibriumSuperhelicalDensity                        = -0.06;
            s.supercoiledSuperhelicalDensityTolerance               = 0.1;
            
            %% compartments
            s.compartmentIndexs_cytosol = 1;
            
            %% genome
            s.sequence = edu.stanford.covert.cell.sim.constant.ChromosomeSequence(reshape([
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
                1,[]));
            s.sequenceLen = size(s.sequence, 1);
            s.sequenceGCContent = getGCContent(s.sequence);
            
            %% genes
            s.gene.wholeCellModelIDs = {'Gene_1';'Gene_2';'Gene_3';'Gene_4';'Gene_5';'Gene_6'};
            s.gene.startCoordinates = 6+15*(0:length(s.gene.wholeCellModelIDs)-1)';
            s.gene.lengths = repmat(9, numel(s.gene.wholeCellModelIDs), 1);
            s.gene.strands = [2; 1; 1; 2; 2; 1];
            s.gene.names = s.gene.wholeCellModelIDs;
            
            %% transcription units
            s.transcriptionUnitWholeCellModelIDs = {'TU_1';'TU_2';'TU_3';'TU_4';};
            s.transcriptionUnitStartCoordinates = s.gene.startCoordinates([1; 2; 4; 6]);
            s.transcriptionUnitLengths = (s.gene.startCoordinates([1; 3; 5; 6]) + s.gene.lengths([1; 3; 5; 6])) - s.transcriptionUnitStartCoordinates + 1;
            s.transcriptionUnitStrands  = s.gene.strands([1; 2; 4; 6]);
            s.transcriptionUnitNames = s.transcriptionUnitWholeCellModelIDs;
            
            %% metabolites
            s.metabolite.dr5pIndexs = 1;
            s.metabolite.waterIndexs  = 2;
            s.metabolite.dnmpIndexs = (3:6)';
            s.metabolite.m6ADIndexs = 7;
            
            this.metabolite.molecularWeights = zeros(7,1);
            
            %% proteins
            monomerWholeCellModelIDs = {
                'MG_001_MONOMER'; %DnaN
                'MG_097_MONOMER'; %Fpg
                'MG_190_MONOMER'; %MgpA
                'MG_235_MONOMER'; %Nfo
                'MG_254_MONOMER'; %LigA
                'MG_262_MONOMER'; %Pol I-like
                'MG_339_MONOMER'; %RecA
                'MG_438_MONOMER'; %EcoD
                'MG_498_MONOMER'; %Ung
                };
            
            complexWholeCellModelIDs = {
                'DNA_POLYMERASE_CORE';     %
                'MG_073_206_421_TETRAMER'; %UvrABC
                'MG_105_OCTAMER';          %DisA
                'MG_184_DIMER';            %MunI
                'MG_244_DIMER';            %PcrA
                'MG_352_DIMER';            %RecU
                'MG_358_359_10MER';        %RuvAB
                };
            
            [~,s.monomerIndexs_ligase]        = ismember('MG_254_MONOMER', monomerWholeCellModelIDs);
            [~,s.complexIndexs_dnaPolymerase] = ismember('DNA_POLYMERASE_CORE', complexWholeCellModelIDs);
            [~,s.complexIndexs_DisA]          = ismember('MG_105_OCTAMER', complexWholeCellModelIDs);
            
            s.monomerDNAFootprints = ones(size(monomerWholeCellModelIDs));
            s.complexDNAFootprints = ones(size(complexWholeCellModelIDs));
            s.monomerDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, size(monomerWholeCellModelIDs));
            s.complexDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, size(complexWholeCellModelIDs));
            s.monomerDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, size(monomerWholeCellModelIDs));
            s.complexDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, size(complexWholeCellModelIDs));
            
            s.reactionMonomerCatalysisMatrix = zeros(15, numel(monomerWholeCellModelIDs));
            s.reactionComplexCatalysisMatrix = zeros(15, numel(complexWholeCellModelIDs));
            s.reactionBoundMonomer = zeros(15, 1);
            s.reactionBoundComplex = zeros(15, 1);
            
            rxIdx = 0;
            
            %DisA can be released by DNA repair enzymes
            [~, dnaRepairEnzymes_monomers] = ismember({
                'MG_001_MONOMER'; %DnaN
                'MG_097_MONOMER'; %Fpg
                'MG_190_MONOMER'; %MgpA
                'MG_235_MONOMER'; %Nfo
                'MG_254_MONOMER'; %LigA
                'MG_262_MONOMER'; %Pol I-like
                'MG_339_MONOMER'; %RecA
                'MG_438_MONOMER'; %EcoD
                'MG_498_MONOMER'; %Ung
                }, monomerWholeCellModelIDs);
            s.reactionMonomerCatalysisMatrix(sub2ind(...
                size(s.reactionMonomerCatalysisMatrix),...
                rxIdx + (1:numel(dnaRepairEnzymes_monomers))',...
                dnaRepairEnzymes_monomers)) = 1;
            s.reactionBoundComplex(rxIdx+(1:numel(dnaRepairEnzymes_monomers))) = s.complexIndexs_DisA;
            rxIdx = rxIdx + numel(dnaRepairEnzymes_monomers);
            
            [~, dnaRepairEnzymes_complexs] = ismember({
                'MG_073_206_421_TETRAMER'; %UvrABC
                'MG_184_DIMER';            %MunI
                'MG_244_DIMER';            %PcrA
                'MG_352_DIMER';            %RecU
                'MG_358_359_10MER';        %RuvAB
                }, complexWholeCellModelIDs);
            s.reactionComplexCatalysisMatrix(sub2ind(...
                size(s.reactionComplexCatalysisMatrix),...
                rxIdx + (1:numel(dnaRepairEnzymes_complexs))',...
                dnaRepairEnzymes_complexs)) = 1;
            s.reactionBoundComplex(rxIdx+(1:numel(dnaRepairEnzymes_complexs))) = s.complexIndexs_DisA;
            rxIdx = rxIdx + numel(dnaRepairEnzymes_complexs);
            
            %RecA can be released by DnaN
            s.reactionMonomerCatalysisMatrix(rxIdx + 1,...
                ismember({'MG_001_MONOMER'}, monomerWholeCellModelIDs)) = 1;
            [~, s.reactionBoundMonomer(rxIdx+1)] = ismember({'MG_339_MONOMER'}, monomerWholeCellModelIDs);
            
            s.reactionThresholds = ...
                sum(s.reactionMonomerCatalysisMatrix, 2) + ...
                sum(s.reactionComplexCatalysisMatrix, 2);
            
            %% state
            s.allocateMemory(1);
        end
    end
    
    %public methods to query state
    methods
        function testSampleAccessibleSites(this)
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            
            s = this.state;
            
            seq = 'AG';
            seqLen = length(seq);
            s.sequence = ChromosomeSequence(seq);
            s.sequenceLen = size(s.sequence, 1);
            s.sequenceGCContent = getGCContent(s.sequence);
            assertEqual((sum(seq == 'G') + sum(seq == 'C')) / seqLen, s.sequence.getGCContent());
            
            %example 0
            subs = [1 1; 1 2; 2 2];
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            s.damagedBases(subs) = ones(size(subs, 1), 1);
            assertEqual(zeros(0, 2), s.sampleAccessibleSites([], 1, char([])));
            
            subs = [1 1; 1 2];
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            s.damagedBases(subs) = ones(size(subs, 1), 1);
            assertEqual([2 1; 2 2], s.sampleAccessibleSites([], 2, char([])));
            
            subs = [1 1; 2 1];
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            s.damagedBases(subs) = ones(size(subs, 1), 1);
            assertEqual(zeros(0,2), s.sampleAccessibleSites([], 2, char([])));
            
            %example 1
            subs = zeros(0,2);
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            s.damagedBases(subs) = ones(size(subs, 1), 1);
            [positions, strands]=find(s.sequence(:, 1:2) == 'A');
            assertEqual(reshape([positions strands], [], 2), ...
                s.sampleAccessibleSites([], 1, 'A'));
            
            %example 2
            subs = [1 1; 1 2];
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            s.damagedBases(subs) = ones(size(subs, 1), 1);
            assertEqual(zeros(0, 2), s.sampleAccessibleSites([], 1, 'A'));
            
            %example 3
            subs = [2 1; 2 2];
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            s.damagedBases(subs) = ones(size(subs, 1), 1);
            assertEqual([1 1], s.sampleAccessibleSites([], 1, 'A'));
            
            %example 4
            subs = [1 1; 2 1];
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            s.damagedBases(subs) = ones(size(subs, 1), 1);
            assertEqual(zeros(0, 2), s.sampleAccessibleSites([], 1, 'A'));
            
            %example 5
            subs = [1 2; 2 2]';
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            s.damagedBases(subs) = ones(size(subs, 1), 1);
            assertEqual(zeros(0,2), s.sampleAccessibleSites([], 1, 'A'));
            
            %example 6
            subs = zeros(0,2);
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            s.damagedBases(subs) = ones(size(subs, 1), 1);
            positions1 = strfind([s.sequence(:,1); s.sequence(1,1)]', 'GA');
            positions2 = strfind([s.sequence(:,2); s.sequence(1,2)]', seqreverse('GA'));
            assertEqual(seqLen, positions1);
            assertEqual([], positions2);
            assertEqual(reshape([positions1 ones(size(positions1)); positions2 2*ones(size(positions2))], [], 2), ...
                s.sampleAccessibleSites([], 1, 'GA'));
            
            %example 7
            subs = zeros(0,2);
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            s.damagedBases(subs) = ones(size(subs, 1), 1);
            positions1=strfind([s.sequence(:,1); s.sequence(1,1)]', 'AG');
            positions2=strfind([s.sequence(:,2); s.sequence(1,2)]', seqreverse('AG'));
            assertEqual(1,positions1);
            assertEqual([],positions2);
            assertEqual(reshape([positions1 ones(size(positions1)); positions2 2*ones(size(positions2))],[],2), ...
                s.sampleAccessibleSites([], 1, 'AG'));
            
            %example 8
            subs = zeros(0,2);
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            s.damagedBases(subs) = ones(size(subs,1),1);
            positions1=strfind([s.sequence(:,1); s.sequence(1,1)]', 'TC');
            positions2=strfind([s.sequence(end,2); s.sequence(:,2)]', seqreverse('TC'));
            assertEqual([],positions1);
            assertEqual(1,positions2);
            assertEqual(reshape([positions1 ones(size(positions1)); positions2 2*ones(size(positions2))],[],2), ...
                s.sampleAccessibleSites([], 1, 'TC'));
            
            %example 9
            subs = zeros(0,2);
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            s.damagedBases(subs) = ones(size(subs,1),1);
            positions1=strfind([s.sequence(:,1); s.sequence(1,1)]', 'CT');
            positions2=strfind([s.sequence(end,2); s.sequence(:,2)]', seqreverse('CT'));
            assertEqual([],positions1);
            assertEqual(2,positions2);
            assertEqual(reshape([positions1 ones(size(positions1)); positions2 2*ones(size(positions2))],[],2), ...
                s.sampleAccessibleSites([], 1, 'CT'));
            
            %example 10
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = 2;
            seq = 'AT';
            s.sequence = ChromosomeSequence(seq);
            s.sequenceLen = size(s.sequence, 1);
            s.sequenceGCContent = getGCContent(s.sequence);
            assertEqual((sum(seq == 'G') + sum(seq == 'C')) / seqLen, s.sequence.getGCContent());
            
            counts = zeros(2,1);
            for i = 1:250
                positionStrands = s.sampleAccessibleSites(0.5, 10, 'A');
                counts(positionStrands(:,1)) = counts(positionStrands(:,1)) + 1;
            end
            assertTrue(range(counts) < 0.20 * max(counts));
        end
        
        function testIsSiteAccessible(this)
            s = this.state;
            
            monomerWholeCellModelIDs = {'m1';'m2'};
            complexWholeCellModelIDs = {'c1';'c2'};
            s.monomerDNAFootprints = ones(size(monomerWholeCellModelIDs));
            s.complexDNAFootprints = ones(size(complexWholeCellModelIDs));
            s.monomerDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, size(monomerWholeCellModelIDs));
            s.complexDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_ssDNA, size(complexWholeCellModelIDs));
            s.monomerDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, size(monomerWholeCellModelIDs));
            s.complexDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_ssDNA, size(complexWholeCellModelIDs));
            
            s.reactionBoundMonomer = [1; 2; 0; 0];
            s.reactionBoundComplex = [0; 0; 1; 2];
            s.reactionMonomerCatalysisMatrix = zeros(4, 2);
            s.reactionComplexCatalysisMatrix = zeros(4, 2);
            s.reactionMonomerCatalysisMatrix(3, 1) = 1;
            s.reactionMonomerCatalysisMatrix(4, 2) = 2;
            s.reactionComplexCatalysisMatrix(1, 1) = 1;
            s.reactionComplexCatalysisMatrix(2, 2) = 2;
            
            s.reactionThresholds = ...
                sum(s.reactionMonomerCatalysisMatrix, 2) + ...
                sum(s.reactionComplexCatalysisMatrix, 2);
            
            %ex 1: chromosomes all accessible
            s.allocateMemory(1);
            s.polymerizedRegions(1, :) = size(s.sequence, 1);
            result = cell(3,1);
            [result{:}] = s.isRegionAccessible(zeros(0,2), 1, [], [], true, [], false, true);
            
            assertEqual({false(0,1); zeros(0,1); zeros(0, 2)}, result);
            
            [tfs, idxs, positionStrands] = s.isRegionAccessible([100 1; 100 2], 1, [], [], true, [], false, true);
            assertEqual(true(2,1), tfs);
            assertEqual([1;2], idxs);
            assertEqual([100 1; 100 2], positionStrands);
            
            %ex 2: some areas not polymerized
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.polymerizedRegions(1, 3) = 10;

            result = cell(3,1);
            [result{:}] = s.isRegionAccessible([1 1; 1 2; 1 3; 1 4], 1, [], [], true, [], false, true);
            assertEqual({logical([1;1;0;0]); [1; 2]; [1 1; 1 2]}, result);

            result = cell(3,1);
            [result{:}] = s.isRegionAccessible([1 1; 1 2; 1 3; 1 4], 1, 1, [], true, [], false, true);
            assertEqual({logical([1;1;0;0]); [1; 2]; [1 1; 1 2]}, result);

            result = cell(3,1);
            [result{:}] = s.isRegionAccessible([1 1; 1 2; 1 3; 1 4], 1, [], 1, true, [], false, true);
            assertEqual({logical([0;0;1;0]); 3; [1 3]}, result);

            %ex 3: protein bound to chromosome
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.monomerBoundSites(100, 1:2:end) = 1;
            result = cell(3,1);
            [result{:}] = s.isRegionAccessible([100 1; 100 2], 1, [], [], true, [], false, true);
            assertEqual({[false;false]; zeros(0, 1); zeros(0, 2)}, result);
            
            [tfs, idxs, positionStrands] = s.isRegionAccessible([100 1; 100 2], 1, [], [], true, [], false, true);
            assertEqual(zeros(0,2), positionStrands);
            assertEqual(zeros(0,1), idxs);
            assertEqual(false(2,1), tfs);
            
            %ex 4: DNA is damaged
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.damagedBases(100, :) = 1;
            assertEqual([false;false], s.isRegionAccessible([100 1; 100 2], 1, [], [], true, [], false, true));
            
            %ex 5: bound proteins
            s.monomerDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, size(monomerWholeCellModelIDs));
            s.complexDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, size(complexWholeCellModelIDs));
            s.monomerDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, size(monomerWholeCellModelIDs));
            s.complexDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, size(complexWholeCellModelIDs));
            
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.monomerBoundSites([100 1]) = 1;
            s.complexBoundSites([101 2]) = 1;
            assertEqual(logical([0;1;1;0;1;1]), s.isRegionAccessible([100 1; 101 1; 102 1; 100 2; 101 2; 102 2], 1, 1, [], true, [], false, true));
            assertEqual(logical([0;0;1;0;0;1]), s.isRegionAccessible([100 1; 101 1; 102 1; 100 2; 101 2; 102 2], 1, 2, [], true, [], false, true));
            assertEqual(logical([1;0;1;1;0;1]), s.isRegionAccessible([100 1; 101 1; 102 1; 100 2; 101 2; 102 2], 1, [], 1, true, [], false, true));
            assertEqual(logical([0;0;1;0;0;1]), s.isRegionAccessible([100 1; 101 1; 102 1; 100 2; 101 2; 102 2], 1, [], 2, true, [], false, true));
            assertEqual(logical([1;1;1;1;1;1]), s.isRegionAccessible([100 1; 101 1; 102 1; 100 2; 101 2; 102 2], 1, 1, 1, true, [], false, true));
            assertEqual(logical([0;1;1;0;1;1]), s.isRegionAccessible([100 1; 101 1; 102 1; 100 2; 101 2; 102 2], 1, 1, 2, true, [], false, true));
            assertEqual(logical([1;0;1;1;0;1]), s.isRegionAccessible([100 1; 101 1; 102 1; 100 2; 101 2; 102 2], 1, 2, 1, true, [], false, true));
            assertEqual(logical([0;0;1;0;0;1]), s.isRegionAccessible([100 1; 101 1; 102 1; 100 2; 101 2; 102 2], 1, 2, 2, true, [], false, true));
            assertEqual(logical([0;1;1;0;1;1]), s.isRegionAccessible([100 1; 101 1; 102 1; 100 2; 101 2; 102 2], 1, [1 2], 2, true, [], false, true));
            assertEqual(logical([1;1;1;1;1;1]), s.isRegionAccessible([100 1; 101 1; 102 1; 100 2; 101 2; 102 2], 1, [1 2], [2 1], true, [], false, true));
            
            [tfs, idxs, positionStrands] = s.isRegionAccessible([100 1; 101 1; 102 1; 100 2; 101 2; 102 2], 1, [1 2], 2, true, [], false, true);
            assertEqual(logical([0;1;1;0;1;1]), tfs);
            assertEqual([101 1; 102 1; 101 2; 102 2], positionStrands);
            assertEqual([2;3;5;6], idxs);
        end
        
        function testIsRegionAccessible(this)
            s = this.state;
            
            monomerWholeCellModelIDs = {'m1';'m2'};
            complexWholeCellModelIDs = {'c1';'c2'};
            s.monomerDNAFootprints = [1; 1];
            s.complexDNAFootprints = [3; 5];
            s.monomerDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, size(monomerWholeCellModelIDs));
            s.complexDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, size(complexWholeCellModelIDs));
            s.monomerDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, size(monomerWholeCellModelIDs));
            s.complexDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, size(complexWholeCellModelIDs));
            
            s.reactionBoundMonomer = [1; 2; 0; 0];
            s.reactionBoundComplex = [0; 0; 1; 2];
            s.reactionMonomerCatalysisMatrix = zeros(4, 2);
            s.reactionComplexCatalysisMatrix = zeros(4, 2);
            s.reactionMonomerCatalysisMatrix(3, 1) = 1;
            s.reactionMonomerCatalysisMatrix(4, 2) = 2;
            s.reactionComplexCatalysisMatrix(1, 1) = 1;
            s.reactionComplexCatalysisMatrix(2, 2) = 2;
            
            s.reactionThresholds = ...
                sum(s.reactionMonomerCatalysisMatrix, 2) + ...
                sum(s.reactionComplexCatalysisMatrix, 2);
            
            %ex 1: chromosomes all accessible
            s.allocateMemory(1);
            s.polymerizedRegions(1, :) = size(s.sequence, 1);
            
            result = cell(4,1);
            [result{:}] = s.isRegionAccessible(zeros(0,2), 1, [], [], true, [], false, true);
            assertEqual({false(0,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            result = cell(4,1);
            [result{:}] = s.isRegionAccessible(zeros(0,2), zeros(0,1), [], [], true, [], false, true);
            assertEqual({false(0,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            %ex 2: some damage on strand 1
            s.allocateMemory(1);
            s.polymerizedRegions(1, :) = size(s.sequence, 1);
            s.damagedBases(5, 1) = 1;
            result = cell(4,1);
            
            [result{:}] = s.isRegionAccessible([5 1], 2, [], 2, true, 1, false, true);
            assertEqual({true(1,1); 1; [5 1]; 2}, result);
            
            [result{:}] = s.isRegionAccessible([5 2], 2, [], 2, true, 1, false, true);
            assertEqual({true(1,1); 1; [5 2]; 2}, result);
            
            [result{:}] = s.isRegionAccessible([5 1], 2, [], 2, true, [-1 1 2], false, true);
            assertEqual({true(1,1); 1; [5 1]; 2}, result);
            
            [result{:}] = s.isRegionAccessible([5 2], 2, [], 2, true, [0 1 6 7], false, true);
            assertEqual({true(1,1); 1; [5 2]; 2}, result);
            
            [result{:}] = s.isRegionAccessible([5 1], 2, [], 2, true, 2, false, true);
            assertEqual({false(1,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            [result{:}] = s.isRegionAccessible([5 2], 2, [], 2, true, 2, false, true);
            assertEqual({false(1,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            [result{:}] = s.isRegionAccessible([5 1], 2, [], 2, true, [-1 2], false, true);
            assertEqual({false(1,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            [result{:}] = s.isRegionAccessible([5 2], 2, [], 2, true, [2 6], false, true);
            assertEqual({false(1,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            [result{:}] = s.isRegionAccessible([5 1], 2, [], 2, true, [], false, true);
            assertEqual({false(1,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            [result{:}] = s.isRegionAccessible([5 2], 2, [], 2, true, [], false, true);
            assertEqual({false(1,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            %ex 3: some damage on strand 2
            s.allocateMemory(1);
            s.polymerizedRegions(1, :) = size(s.sequence, 1);
            s.damagedBases(5, 2) = 1;
            result = cell(4,1);
            
            [result{:}] = s.isRegionAccessible([5 1], 2, [], 2, true, 1, false, true);
            assertEqual({true(1,1); 1; [5 1]; 2}, result);
            
            [result{:}] = s.isRegionAccessible([5 2], 2, [], 2, true, 1, false, true);
            assertEqual({true(1,1); 1; [5 2]; 2}, result);
            
            [result{:}] = s.isRegionAccessible([5 1], 2, [], 2, true, [-1 1 2], false, true);
            assertEqual({true(1,1); 1; [5 1]; 2}, result);
            
            [result{:}] = s.isRegionAccessible([5 2], 2, [], 2, true, [0 1 6 7], false, true);
            assertEqual({true(1,1); 1; [5 2]; 2}, result);
            
            [result{:}] = s.isRegionAccessible([5 1], 2, [], 2, true, 2, false, true);
            assertEqual({false(1,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            [result{:}] = s.isRegionAccessible([5 2], 2, [], 2, true, 2, false, true);
            assertEqual({false(1,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            [result{:}] = s.isRegionAccessible([5 1], 2, [], 2, true, [-1 2], false, true);
            assertEqual({false(1,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            [result{:}] = s.isRegionAccessible([5 2], 2, [], 2, true, [2 6], false, true);
            assertEqual({false(1,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            [result{:}] = s.isRegionAccessible([5 1], 2, [], 2, true, [], false, true);
            assertEqual({false(1,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
            
            [result{:}] = s.isRegionAccessible([5 2], 2, [], 2, true, [], false, true);
            assertEqual({false(1,1); zeros(0,1); zeros(0, 2); zeros(0,1)}, result);
        end
    end
    
    %private methods for query state by region
    methods
        %tests isSiteSingleStranded, isSiteDoubleStranded
        function testIsRegionSingleDoubleStranded(this)
            s = this.state;
            
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1) = 50;
            s.polymerizedRegions(1, 2) = 70;
            s.polymerizedRegions(91, 1) = 10;
            s.polymerizedRegions(91, 3) = 10;
            s.polymerizedRegions(91, 4) = 20;
            
            %test sites
            [tfs, idxs, positionsStrands, lengths] = s.isRegionSingleStranded(...
                [1 1; 1 2; 50 1; 51 1; 51 2; 90 1; 90 2; 91 1; 92 1; 100 1; 100 2; 101 2; 91 3; 91 4; 100 3; 100 4; 101 3; 101 4; 110 3; 110 4; 111 3; 111 4], 1, false);
            assertEqual(logical([0;0;0;0;1;0;0;1;1;1;0;0;0;0;0;0;0;1;0;1;0;0]), tfs);
            assertEqual([5; 8; 9; 10; 18; 20], idxs);
            assertEqual([51 2; 91 1; 92 1; 100 1; 101 4; 110 4], positionsStrands);
            assertEqual([1;1;1;1;1;1], lengths);

            [tfs, idxs, positionsStrands, lengths] = s.isRegionDoubleStranded(...
                [1 1; 1 2; 50 1; 51 1; 51 2; 90 1; 90 2; 91 1; 92 1; 100 1; 100 2; 101 2; 91 3; 91 4; 100 3; 100 4; 101 3; 101 4; 110 3; 110 4; 111 3; 111 4], 1, false);
            assertEqual(logical([1;1;1;0;0;0;0;0;0;0;0;0;1;1;1;1;0;0;0;0;0;0]), tfs);
            assertEqual([1; 2; 3; 13; 14; 15; 16], idxs);
            assertEqual([1 1; 1 2; 50 1; 91 3; 91 4; 100 3; 100 4], positionsStrands);
            assertEqual([1;1;1;1;1;1;1], lengths);
            
            %test regions
            [tfs, idxs, positionsStrands, lengths] = s.isRegionSingleStranded(...
                [1 1; 1 1; 50 2; 51 2; 51 2; 90 2; 91 1; 92 2; 91 3; 91 4; 101 3; 101 4],...
                [50; 70; 20; 20; 21; 10; 10; 10; 10; 10; 10; 10], false);
            assertEqual(logical([0;0;0;1;0;0;1;0;0;0;0;1]), tfs);
            assertEqual([4; 7; 12], idxs);
            assertEqual([51 2; 91 1; 101 4], positionsStrands);
            assertEqual([20; 10; 10], lengths);
            
            [tfs, idxs, positionsStrands, lengths] = s.isRegionDoubleStranded(...
                [1 1; 1 1; 50 2; 51 2; 51 2; 90 2; 91 1; 92 2; 91 3; 91 4; 101 3; 101 4],...
                [50; 70; 20; 20; 21; 10; 10; 10; 10; 10; 10; 10], false);
            assertEqual(logical([1;0;0;0;0;0;0;0;1;1;0;0]), tfs);
            assertEqual([1; 9; 10], idxs);
            assertEqual([1 1; 91 3; 91 4], positionsStrands);
            assertEqual([50;10;10], lengths);
        end
        
        %tests isRegionPolymerized
        function testIsRegionPolymerized(this)
            s = this.state;
            
            %% test sites
            %ex 1: nothing polymerized, 1 query site
            s.allocateMemory(1);
            result = cell(3,1);
            [result{:}] = s.isRegionPolymerized([1 1], 1, false);
            assertEqual({false; zeros(0,1); zeros(0,2)}, result);
            
            %ex 2: nothing polymerized, multiple query sites
            s.allocateMemory(1);
            result = cell(3,1);
            [result{:}] = s.isRegionPolymerized([1 1; 2 2; 3 3], 1, false);
            assertEqual({false(3,1); zeros(0,1); zeros(0,2)}, result);
            
            %ex 3: something polymerized, 1 query sites that is polymerized
            s.allocateMemory(1);
            s.polymerizedRegions(1,1) = 1;
            result = cell(3,1);
            [result{:}] = s.isRegionPolymerized([1 1], 1, false);
            assertEqual({true(1,1); 1; [1 1]}, result);
            
            %ex 4: something polymerized, 1 query sites that is not polymerized
            s.allocateMemory(1);
            s.polymerizedRegions(1,1) = 1;
            result = cell(3,1);
            [result{:}] = s.isRegionPolymerized([2 1], 1, false);
            assertEqual({false(1,1); zeros(0,1); zeros(0,2)}, result);
            
            s.allocateMemory(1);
            s.polymerizedRegions(1,1) = 1;
            result = cell(3,1);
            [result{:}] = s.isRegionPolymerized([1 2], 1, false);
            assertEqual({false(1,1); zeros(0,1); zeros(0,2)}, result);
            
            s.allocateMemory(1);
            s.polymerizedRegions(1,1) = 1;
            result = cell(3,1);
            [result{:}] = s.isRegionPolymerized([1 3], 1, false);
            assertEqual({false(1,1); zeros(0,1); zeros(0,2)}, result);
            
            s.allocateMemory(1);
            s.polymerizedRegions(2,1) = 1;
            result = cell(3,1);
            [result{:}] = s.isRegionPolymerized([1 1], 1, false);
            assertEqual({false(1,1); zeros(0,1); zeros(0,2)}, result);
            
            %ex 5: something polymerized, multiple query sites that are not polymerized
            s.allocateMemory(1);
            s.polymerizedRegions(2,1) = 1;
            result = cell(3,1);
            [result{:}] = s.isRegionPolymerized([1 1;2 1; 2 2; 3 3], 1, false);
            assertEqual({[false;true;false;false]; 2; [2 1]}, result);
            
            s.allocateMemory(1);
            s.polymerizedRegions(2,3) = 5;
            result = cell(3,1);
            [result{:}] = s.isRegionPolymerized([1 1; 2 1; 2 2; 3 3; 6 3; 7 3; 4 4], 1, false);
            assertEqual({[false(3,1); true(2,1); false(2, 1)]; [4; 5]; [3 3; 6 3]}, result);
            
            %% test regions
            %ex 1: nothing polymerized, 1 query site
            s.allocateMemory(1);
            result = cell(4,1);
            [result{:}] = s.isRegionPolymerized([1 1], 1, false);
            assertEqual({false; zeros(0,1); zeros(0,2); zeros(0,1)}, result);
            
            %ex 2: nothing polymerized, 1 query region
            s.allocateMemory(1);
            result = cell(4,1);
            [result{:}] = s.isRegionPolymerized([1 1], 10, false);
            assertEqual({false; zeros(0,1); zeros(0,2); zeros(0,1)}, result);
            
            %ex 3: site polymerized, 1 query site
            s.allocateMemory(1);
            s.polymerizedRegions(1,1) = 1;
            result = cell(4,1);
            [result{:}] = s.isRegionPolymerized([1 1], 1, false);
            assertEqual({true; 1; [1 1]; 1}, result);
            
            %ex 4: site polymerized, 1 query region
            s.allocateMemory(1);
            s.polymerizedRegions(1,1) = 1;
            result = cell(4,1);
            [result{:}] = s.isRegionPolymerized([1 1], 2, false);
            assertEqual({false; zeros(0,1); zeros(0,2); zeros(0,1)}, result);
            
            %ex 5: region polymerized, 1 query region
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1) = 10;
            result = cell(4,1);
            [result{:}] = s.isRegionPolymerized([1 1], 10, false);
            assertEqual({true; 1; [1 1]; 10}, result);
            
            %ex 6: region partially polymerized, 1 query region
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1) = 10;
            result = cell(4,1);
            [result{:}] = s.isRegionPolymerized([3 1], 10, false);
            assertEqual({false; zeros(0,1); zeros(0,2); zeros(0,1)}, result);
            
            %ex 7: various polymerization, multiple query regions including two
            %that wrap ORI
            s.allocateMemory(1);
            s.polymerizedRegions(2, 1) = 15;
            s.polymerizedRegions(19, 1) = 3;
            s.polymerizedRegions(9, 2) = 2;
            s.polymerizedRegions(end, 3) = 1;
            s.polymerizedRegions(1, 3) = 8;
            s.polymerizedRegions(end, 4) = 1;
            s.polymerizedRegions(1, 4) = 9;
            result = cell(4,1);
            [result{:}] = s.isRegionPolymerized([3 1; 20 1; 8 2; 20 3; 40 3; size(s.sequence,1) 3; size(s.sequence,1) 4], [10; 3; 5; 1; 2; 10; 10], false);
            assertEqual({[true; false(5,1); true]; [1; 7]; [3 1; size(s.sequence,1) 4]; [10; 10]}, result);
        end
        
        %tests isRegionProteinFree
        function testIsRegionProteinFree(this)
            s = this.state;
            
            %% site
            s.allocateMemory(1);
            s.monomerDNAFootprints = [1; 2; 4];
            s.complexDNAFootprints = [3; 5; 7];
            s.monomerBoundSites(end, 1) = 1;
            s.complexBoundSites(end, 3) = 1;
            s.complexBoundSites(4, 3) = 2;
            
            result = cell(6,1);
            positionsStrands = [
                size(s.sequence, 1) 1;
                size(s.sequence, 1) 2;
                1 1;
                1 2;
                size(s.sequence, 1) 3;
                size(s.sequence, 1) 4;
                1 3;
                1 4;
                4 3;
                4 4;
                8 3;
                8 4;
                9 3;
                9 4];
            [result{:}] = s.isRegionProteinFree(positionsStrands, 1, false, [], [], true);
            
            tfs = logical([0 0 1 1 0 0 0 0 0 0 0 0 1 1]');
            assertEqual(tfs, result{1});
            assertEqual(reshape(find(tfs), [], 1), result{2});
            assertEqual(positionsStrands(tfs, :), result{3});
            assertEqual(ones(sum(tfs), 1), result{4});
            assertEqual([1 1 0 0 0 0 0 0 0 0 0 0 0 0], full(result{5})');
            assertEqual([0 0 0 0 1 1 1 1 2 2 2 2 0 0], full(result{6})');

            %% region
            s.allocateMemory(1);
            s.monomerDNAFootprints = [1; 2; 4];
            s.complexDNAFootprints = [3; 5; 7];
            s.monomerBoundSites(end, 1) = 1;
            s.complexBoundSites(end, 3) = 1;
            s.complexBoundSites(4, 3) = 2;
            
            result = cell(6,1);
            positionsStrands = [
                size(s.sequence, 1) 1;
                size(s.sequence, 1) 2;
                1 1;
                1 2;
                size(s.sequence, 1) 3;
                size(s.sequence, 1) 4;
                1 3;
                1 4;
                4 3;
                4 4];
            lengths = [2; 2; 10; 10; 5; 5; 5; 5; 6; 6];
            [result{:}] = s.isRegionProteinFree(positionsStrands, lengths, false, [], [], true);
            
            tfs = ~[true(2,1); false(2,1); true(6,1)];
            assertEqual(tfs, result{1});
            assertEqual(reshape(find(tfs), [], 1), result{2});
            assertEqual(positionsStrands(tfs, :), result{3});
            assertEqual(lengths(tfs, :), result{4});
            mat = zeros(10,10);
            mat(1:2, 1) = 1;
            assertEqual(mat, full(result{5}));
            mat = zeros(10,10);
            mat(5:8, 1) = 1;
            mat(9:10, 1) = 2;
            mat(5:8,4) = 2;
            assertEqual(mat, full(result{6}));
        end

        %tests isRegionDamaged
        function testIsRegionUndamaged(this)
            s = this.state;
            
            %% site
            s.allocateMemory(1);
            s.damagedBases(1,1) = 1;
            s.gapSites(3,2) = 1;
            s.intrastrandCrossLinks(7,1) = 1;
            s.strandBreaks(11,1) = 1;
            s.hollidayJunctions(15,2) = 1;
            
            posStrnds = find(s.damagedSites);
            assertEqual([
                1 1
                7 1
                8 1
                11 1
                12 1
                3 2
                14 2
                15 2
                ], posStrnds);
            
            result = cell(3,1);
            positionsStrands = [1 1; 1 2; 3 1; 3 2; 6 1; 7 1; 8 1; 10 1; 11 1; 12 1; 14 2; 15 2; 16 2];
            [result{:}] = s.isRegionUndamaged(positionsStrands, 1, false, [], false);
            tfs = ~[true; false; false; true; false; true; true; false; true; true; true; true; false];
            assertEqual(tfs', result{1}');
            assertEqual(reshape(find(tfs),[],1), result{2});
            assertEqual(positionsStrands(tfs, :), result{3});
            
            %% region
            s.allocateMemory(1);
            s.damagedBases(1,1) = 1;
            s.gapSites(10,2) = 1;
            s.intrastrandCrossLinks(20,1) = 1;
            s.strandBreaks(30,1) = 1;
            s.hollidayJunctions(40,2) = 1;
            
            posStrnds = find(s.damagedSites);
            assertEqual([
                1 1
                20 1
                21 1
                30 1
                31 1
                10 2
                39 2
                40 2                
                ], posStrnds);
            
            result = cell(4,1);
            posStrnds = [1 1; 1 2; 1 2; 21 1; 22 1; 39 2; 40 2; 41 2; 42 2];
            lens = [9; 9; 10; 5; 4; 10; 9; 8; 7];
            [result{:}] = s.isRegionUndamaged(posStrnds, lens, false, [], false);
            tfs = ~[true; false; true; true; false; true; true; false; false];
            assertEqual(tfs', result{1}');
            assertEqual(reshape(find(tfs),[],1), result{2});
            assertEqual(posStrnds(tfs, :), result{3});
            assertEqual(lens(tfs, :), result{4});
        end
    end
    
    %additional private methods used to query state
    methods
        function testGetStrandView(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            s = this.state;
            
            %ex 1
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.setRegionUnwound(1, 10);
            s.setRegionPolymerized([1 1], 10);
            s.damagedBases(3:5, 4) = [1; 2; 3];
            s.damagedBases(2:4, 2) = [4; 5; 6];
            
            polymerizedRegions = CircularSparseMat([1 1; 1 2; 1 4], [size(s.sequence, 1); size(s.sequence, 1); 10], [size(s.sequence, 1) s.nCompartments], 1);
            assertEqual(polymerizedRegions, s.getStrandView());
            assertEqual(polymerizedRegions, s.getStrandView('polymerizedRegions'));
            assertEqual(polymerizedRegions, s.getStrandView({'polymerizedRegions'}));
            
            damagedBases = CircularSparseMat([3 2; 4 2; 5 2; 2 4; 3 4; 4 4], [1; 2; 3; 4; 5; 6], [size(s.sequence, 1) s.nCompartments], 1);
            assertEqual(damagedBases, s.getStrandView('damagedBases'));
            
            %ex 2
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.damagedBases(3:5, 3) = [1; 2; 3];
            s.damagedBases(2:4, 1) = [4; 5; 6];
            s.damagedBases(5:7, 2) = [7; 8; 9];
            s.setRegionUnwound(1, 10);
            s.setRegionPolymerized([1 2], 10);
            
            polymerizedRegions = CircularSparseMat([1 1; 1 2; 1 3], [size(s.sequence, 1); size(s.sequence, 1); 10], [size(s.sequence, 1) s.nCompartments], 1);
            assertEqual(polymerizedRegions, s.getStrandView('polymerizedRegions'));
            
            damagedBases = CircularSparseMat([3 3; 4 3; 5 3; 2 1; 3 1; 4 1; 5 2; 6 2; 7 2], [1;2;3;4;5;6;7;8;9], [size(s.sequence, 1) s.nCompartments], 1);
            assertEqual(damagedBases, s.getStrandView('damagedBases'));
        end
        
        function testGetDNAFootprint(this)
            s = this.state;
                        
            s.monomerDNAFootprints = [1; 2; 4];
            s.complexDNAFootprints = [3; 5; 7];
            s.monomerDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, 3, 1);
            s.complexDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, 3, 1);
            s.monomerDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, 3, 1);
            s.complexDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, 3, 1);
            
            result = cell(3,1);
            
            %ex 1
            [result{:}] = s.getDNAFootprint(1, []);
            assertEqual({1; 0; 0}, result);
            
            %ex 2
            [result{:}] = s.getDNAFootprint(2, []);
            assertEqual({2; 0; 1}, result);
            
            %ex 3
            [result{:}] = s.getDNAFootprint([], 3);
            assertEqual({7; 3; 3}, result);
            
            %ex 4
            [result{:}] = s.getDNAFootprint([1 3], []);
            assertEqual({4; 1; 2}, result);
            
            %ex 5
            [result{:}] = s.getDNAFootprint(3, 1);
            assertEqual({4; 1; 2}, result);
            
            %ex 6
            [result{:}] = s.getDNAFootprint(3, [1 2]);
            assertEqual({5; 2; 2}, result);
        end
        
        function testGetAccessibleRegions(this)
            c = this.state;
            
            L = length(c.sequence);
            protDNAFootprint = 630;
            c.monomerDNAFootprints = protDNAFootprint;
            c.complexDNAFootprints = protDNAFootprint;
            c.monomerDNAFootprintBindingStrandedness = c.dnaStrandedness_dsDNA;
            c.complexDNAFootprintBindingStrandedness = c.dnaStrandedness_dsDNA;
            c.monomerDNAFootprintRegionStrandedness = c.dnaStrandedness_dsDNA;
            c.complexDNAFootprintRegionStrandedness = c.dnaStrandedness_dsDNA;
            c.reactionMonomerCatalysisMatrix = zeros(0, 1);
            c.reactionComplexCatalysisMatrix = zeros(0, 1);
            c.reactionBoundMonomer = zeros(0, 1);
            c.reactionBoundComplex = zeros(0, 1);
            c.reactionThresholds = ...
                sum(c.reactionMonomerCatalysisMatrix, 2) + ...
                sum(c.reactionComplexCatalysisMatrix, 2);
            
            %% ex 1: first chromosome completely synthesized
            c.initialize();
            
            %ex 1.1
            c.monomerBoundSites(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions([], 1);
            assertEqual(1, rgnPosStrnds(:, 1));
            assertEqual(L, rgnPosStrnds(:, 1) + rgnLens - 1);
            assertEqual(1, 2*ceil(rgnPosStrnds(:,2)/2)-1);
            
            %ex 1.2
            c.monomerBoundSites(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            c.monomerBoundSites([1 1; 1000 1; 5000 2]) = 1;
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions(1, []);
            assertEqual([1000 + protDNAFootprint; 5000 + protDNAFootprint], rgnPosStrnds(:, 1));
            assertEqual([5000-1; L], rgnPosStrnds(:, 1) + rgnLens - 1);
            assertEqual([1; 1], 2*ceil(rgnPosStrnds(:,2)/2)-1);
            
            %ex 1.3
            c.monomerBoundSites(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([50 1; 1000 2; 5000 1]) = 1;
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions([], 1);
            assertEqual([1000 + protDNAFootprint; 5000 + protDNAFootprint], rgnPosStrnds(:, 1));
            assertEqual([5000-1; L + 50-1], rgnPosStrnds(:, 1) + rgnLens - 1);
            assertEqual([1; 1], 2*ceil(rgnPosStrnds(:,2)/2)-1);
            
            %ex 1.4
            c.monomerBoundSites(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([1000 end], 1) = 1;
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions([], 1);
            assertEqual(1000 + protDNAFootprint, rgnPosStrnds(:, 1));
            assertEqual(L-1, rgnPosStrnds(:, 1) + rgnLens - 1);
            assertEqual(1, 2*ceil(rgnPosStrnds(:,2)/2)-1);
            
            %ex 1.5
            c.monomerBoundSites(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([1000 end-50], 2) = 1;
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions([], 1);
            assertEqual(1000 + protDNAFootprint, rgnPosStrnds(:, 1));
            assertEqual(L-50-1, rgnPosStrnds(:, 1) + rgnLens - 1);
            assertEqual(1, 2*ceil(rgnPosStrnds(:,2)/2)-1);
            
            %% ex 2:
            c.initialize();
            c.setRegionUnwound(1, 1000);
            c.setRegionPolymerized([1 2], 1000);
            
            %ex 2.1
            c.monomerBoundSites(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([5000 1; 200 3]) = 1;
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions([], 1);
            assertEqual([
                1001                    1;
                5000 + protDNAFootprint 1],...
                rgnPosStrnds);
            assertEqual([5000-1; size(c.sequence,1)], rgnPosStrnds(:, 1) + rgnLens - 1);
            
            %ex 2.2
            c.monomerBoundSites(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites([2 4]) = 1;
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions([], 1);
            assertEqual([1001 1], rgnPosStrnds);
            assertEqual(L, rgnPosStrnds(:, 1) + rgnLens - 1);
            
            %% ex 3:
            c.initialize();
            c.setRegionUnwound(L, -1000);
            c.setRegionPolymerized([L 2], -1000);
            
            %ex 3.1
            c.monomerBoundSites(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(end-50,4) = 1;
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions([], 1);
            assertEqual([1; L-1000+1], rgnPosStrnds(:, 1));
            assertEqual([L-1000; L-50-1], rgnPosStrnds(:, 1) + rgnLens - 1);
            assertEqual([1; 3], 2*ceil(rgnPosStrnds(:,2)/2)-1);
            
            %ex 3.2
            c.monomerBoundSites(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(end-800,3) = 1;
            c.complexBoundSites(50,4) = 1;
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions([], 1);
            assertEqual(1, rgnPosStrnds(:, 1));
            assertEqual(L-1000, rgnPosStrnds(:, 1) + rgnLens - 1);
            assertEqual(1, 2*ceil(rgnPosStrnds(:,2)/2)-1);
        end
    end
    
    %accessors for particular processes
    methods
        function testPolymerizedGenes(this)
            s = this.state;
            A = s.polymerizedGenes;
            assertEqual([525 2], size(A));
            assertAllEqual(true, A(:,1)');
            assertAllEqual(false, A(:,2)');

            s.polymerizedRegions(1, s.strandIndexs_ch2) = size(s.sequence, 1);
            A = s.polymerizedGenes;
            assertEqual([525 2], size(A));
            assertAllEqual(true, A(:,1)');
            assertAllEqual(true, A(:,2)');
        end
        
        function testGetCopyNumberGenes(this)
            s = this.state;
            assertEqual(ones(525, 1), s.geneCopyNumbers);

            s.polymerizedRegions(1, s.strandIndexs_ch2) = size(s.sequence, 1);
            assertEqual(repmat(2, [525 1]), s.geneCopyNumbers);
        end
        
        function testAreTranscriptionUnitsPolymerized(this)
            s = this.state;            
            
            A = s.polymerizedTranscriptionUnits;
            assertEqual([335 2], size(A));
            assertAllEqual(true, A(:,1)');
            assertAllEqual(false, A(:,2)');

            s.polymerizedRegions(1, s.strandIndexs_ch2) = size(s.sequence, 1);
            A = s.polymerizedTranscriptionUnits;
            assertEqual([335 2], size(A));
            assertAllEqual(true, A(:,1)');
            assertAllEqual(true, A(:,2)');
        end

        function testGetCopyNumberTranscriptionUnits(this)
            s = this.state;
            assertEqual(ones(335, 1), s.transcriptionUnitCopyNumbers);

            s.polymerizedRegions(1, s.strandIndexs_ch2) = size(s.sequence, 1);
            assertEqual(repmat(2, [335 1]), s.transcriptionUnitCopyNumbers);
        end
        
        %tests accessibleGenes, geneCopyNumbers_Accessible, and implicitly tests isRegionAccessible
        function testAccessibleGenes(this)
            s = this.state;
            
            nGenes = size(s.gene.wholeCellModelIDs,1);
            geneStartCoordinates = s.gene.startCoordinates;
            geneEndCoordinates = geneStartCoordinates + s.gene.lengths - 1;
            
            geneIdx = find(...
                geneStartCoordinates(2:end-1) > geneEndCoordinates(1:end-2) + 1 & ...
                geneEndCoordinates(2:end-1)   < geneStartCoordinates(3:end) - 1, 1, 'first') + 1;
            startCoor = geneStartCoordinates(geneIdx);
            endCoor = geneEndCoordinates(geneIdx);
            dir = mod(s.gene.strands(geneIdx), 2);
            
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence,1);
            accessibleGenes = [true(nGenes, 1) false(nGenes,1)];
            assertEqual(accessibleGenes, s.accessibleGenes)
            assertEqual(sum(accessibleGenes,2), s.geneCopyNumbers_Accessible);
            
            outsideCoors = [startCoor-1;endCoor+1];
            s.damagedBases(outsideCoors, :) = 1;
            accessibleGenes(geneIdx)=1;
            assertEqual(accessibleGenes, s.accessibleGenes);
            
            s.damagedBases(startCoor, :) = 1;
            accessibleGenes(geneIdx)=0;
            assertEqual(accessibleGenes, s.accessibleGenes);
            
            s.damagedBases(startCoor, ~dir+1) = 0;
            accessibleGenes(geneIdx)=0;
            assertEqual(accessibleGenes, s.accessibleGenes);
            
            s.damagedBases(startCoor, dir+1) = 0;
            accessibleGenes(geneIdx)=1;
            assertEqual(accessibleGenes, s.accessibleGenes);
            
            s.damagedBases(endCoor, :) = 1;
            accessibleGenes(geneIdx)=0;
            assertEqual(accessibleGenes, s.accessibleGenes);
            
            s.damagedBases(endCoor, ~dir+1) = 0;
            accessibleGenes(geneIdx)=0;
            assertEqual(accessibleGenes, s.accessibleGenes);
            
            s.damagedBases(endCoor, dir+1) = 0;
            accessibleGenes(geneIdx)=1;
            assertEqual(accessibleGenes, s.accessibleGenes);
        end
        
        %tests accessibleTranscriptionUnits,
        %getCopyNumberTranscriptionUnits_Accessible, and implicitly tests
        %isRegionAccessible
        function testAccessibleTranscriptionUnits(this)
            s = this.state;
            
            nTranscriptionUnits = size(s.transcriptionUnitWholeCellModelIDs,1);
            transcriptionUnitStartCoordinates = s.transcriptionUnitStartCoordinates;
            transcriptionUnitEndCoordinates = transcriptionUnitStartCoordinates + s.transcriptionUnitLengths - 1;
            transcriptionUnitIdx = find(...
                transcriptionUnitStartCoordinates(2:end-1) > transcriptionUnitEndCoordinates(1:end-2)+1 & ...
                transcriptionUnitEndCoordinates(2:end-1) < transcriptionUnitStartCoordinates(3:end)-1,1,'first')+1;
            startCoor = transcriptionUnitStartCoordinates(transcriptionUnitIdx);
            endCoor = transcriptionUnitEndCoordinates(transcriptionUnitIdx);
            dir = mod(s.transcriptionUnitStrands(transcriptionUnitIdx), 2);
            
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence,1);
            accessibleTranscriptionUnits = [true(nTranscriptionUnits, 1) false(nTranscriptionUnits, 1)];
            assertEqual(accessibleTranscriptionUnits, s.accessibleTranscriptionUnits)
            assertEqual(sum(accessibleTranscriptionUnits,2), s.transcriptionUnitCopyNumbers_Accessible);
            
            outsideCoors = [startCoor-1;endCoor+1];
            s.damagedBases(outsideCoors, :) = 1;
            accessibleTranscriptionUnits(transcriptionUnitIdx)=1;
            assertEqual(accessibleTranscriptionUnits, s.accessibleTranscriptionUnits);
            
            s.damagedBases(startCoor, :) = 1;
            accessibleTranscriptionUnits(transcriptionUnitIdx)=0;
            assertEqual(accessibleTranscriptionUnits, s.accessibleTranscriptionUnits);
            
            s.damagedBases(startCoor, ~dir+1) = 0;
            accessibleTranscriptionUnits(transcriptionUnitIdx)=0;
            assertEqual(accessibleTranscriptionUnits, s.accessibleTranscriptionUnits);
            
            s.damagedBases(startCoor, dir+1) = 0;
            accessibleTranscriptionUnits(transcriptionUnitIdx)=1;
            assertEqual(accessibleTranscriptionUnits, s.accessibleTranscriptionUnits);
            
            s.damagedBases(endCoor, :) = 1;
            accessibleTranscriptionUnits(transcriptionUnitIdx)=0;
            assertEqual(accessibleTranscriptionUnits, s.accessibleTranscriptionUnits);
            
            s.damagedBases(endCoor, ~dir+1) = 0;
            accessibleTranscriptionUnits(transcriptionUnitIdx)=0;
            assertEqual(accessibleTranscriptionUnits, s.accessibleTranscriptionUnits);
            
            s.damagedBases(endCoor, dir+1) = 0;
            accessibleTranscriptionUnits(transcriptionUnitIdx)=1;
            assertEqual(accessibleTranscriptionUnits, s.accessibleTranscriptionUnits);
        end
        
        function testStrandBreakClassification(this)
            s = this.state;
            s.strandBreakClassification_segmentLength=51;
            
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence,1);
            initial_strandBreaks = s.strandBreaks;
            N = ceil(size(s.sequence,1)/s.strandBreakClassification_segmentLength);
            
            %example 1
            s.strandBreaks  = initial_strandBreaks;
            strandBreakClassification = s.strandBreakClassification;
            assertEqual(N, strandBreakClassification(s.strandBreakClassification_index_NB));
            assertEqual([0;N], unique(strandBreakClassification));
            
            %example 2
            s.strandBreaks  = initial_strandBreaks;
            s.strandBreaks(100,1)=1;
            strandBreakClassification = s.strandBreakClassification;
            assertEqual(N-1, strandBreakClassification(s.strandBreakClassification_index_NB));
            assertEqual(1, strandBreakClassification(s.strandBreakClassification_index_SSB));
            assertEqual([0;1;N-1], unique(strandBreakClassification));
            
            %example 3
            s.strandBreaks  = initial_strandBreaks;
            s.strandBreaks(s.strandBreakClassification_segmentLength+[1;1+s.strandBreakClassification_segmentLength],1)=1;
            strandBreakClassification = s.strandBreakClassification;
            assertEqual(N-2, strandBreakClassification(s.strandBreakClassification_index_NB));
            assertEqual(2, strandBreakClassification(s.strandBreakClassification_index_SSB));
            assertEqual([0;2;N-2], unique(strandBreakClassification));
            
            %example 4
            s.strandBreaks  = initial_strandBreaks;
            s.strandBreaks(s.strandBreakClassification_segmentLength+[1;s.strandBreakClassification_segmentLength],1)=1;
            strandBreakClassification = s.strandBreakClassification;
            assertEqual(N-1, strandBreakClassification(s.strandBreakClassification_index_NB));
            assertEqual(1, strandBreakClassification(s.strandBreakClassification_index_SSB_));
            assertEqual([0;1;N-1], unique(strandBreakClassification));
            
            %example 5
            s.strandBreaks  = initial_strandBreaks;
            s.strandBreaks([s.strandBreakClassification_segmentLength+1 1;s.strandBreakClassification_segmentLength+1+s.strandBreakClassification_doubleStrandBreakSeparation 2])=1;
            strandBreakClassification = s.strandBreakClassification;
            assertEqual(N-1, strandBreakClassification(s.strandBreakClassification_index_NB));
            assertEqual(1, strandBreakClassification(s.strandBreakClassification_index_2SSB));
            assertEqual([0;1;N-1], unique(strandBreakClassification));
            
            %example 6
            s.strandBreaks  = initial_strandBreaks;
            s.strandBreaks([s.strandBreakClassification_segmentLength+1 1;s.strandBreakClassification_segmentLength+s.strandBreakClassification_doubleStrandBreakSeparation 2])=1;
            strandBreakClassification = s.strandBreakClassification;
            assertEqual(N-1, strandBreakClassification(s.strandBreakClassification_index_NB));
            assertEqual(1, strandBreakClassification(s.strandBreakClassification_index_DSB));
            assertEqual([0;1;N-1], unique(strandBreakClassification));
            
            %example 7
            s.strandBreaks  = initial_strandBreaks;
            s.strandBreaks(1,1)=1;
            s.strandBreaks(end,2)=1;
            strandBreakClassification = s.strandBreakClassification;
            assertEqual(N-2, strandBreakClassification(s.strandBreakClassification_index_NB));
            assertEqual(1, strandBreakClassification(s.strandBreakClassification_index_SSB));
            assertEqual(1, strandBreakClassification(s.strandBreakClassification_index_DSB));
            assertEqual([0;1;N-2], unique(strandBreakClassification));
            
            %example 8
            s.strandBreaks  = initial_strandBreaks;
            s.strandBreaks(s.strandBreakClassification_segmentLength+[1 s.strandBreakClassification_doubleStrandBreakSeparation+1],1)=1;
            s.strandBreaks(s.strandBreakClassification_segmentLength+1,2)=1;
            strandBreakClassification = s.strandBreakClassification;
            assertEqual(N-1, strandBreakClassification(s.strandBreakClassification_index_NB));
            assertEqual(1, strandBreakClassification(s.strandBreakClassification_index_DSB));
            assertEqual([0;1;N-1], unique(strandBreakClassification));
            
            %example 9
            s.strandBreaks  = initial_strandBreaks;
            s.strandBreaks(s.strandBreakClassification_segmentLength+[1 s.strandBreakClassification_doubleStrandBreakSeparation],1)=1;
            s.strandBreaks(s.strandBreakClassification_segmentLength+1,2)=1;
            strandBreakClassification = s.strandBreakClassification;
            assertEqual(N-1, strandBreakClassification(s.strandBreakClassification_index_NB));
            assertEqual(1, strandBreakClassification(s.strandBreakClassification_index_DSB_));
            assertEqual([0;1;N-1], unique(strandBreakClassification));
            
            %example 10
            s.strandBreaks  = initial_strandBreaks;
            s.strandBreaks(s.strandBreakClassification_segmentLength+1,1:2)=1;
            s.strandBreaks(s.strandBreakClassification_segmentLength+1+s.strandBreakClassification_doubleStrandBreakSeparation,1:2)=1;
            strandBreakClassification = s.strandBreakClassification;
            assertEqual(N-1, strandBreakClassification(s.strandBreakClassification_index_NB));
            assertEqual(1, strandBreakClassification(s.strandBreakClassification_index_DSB__));
            assertEqual([0;1;N-1], unique(strandBreakClassification));
            
            %example 11
            s.strandBreaks  = initial_strandBreaks;
            s.strandBreaks(s.strandBreakClassification_segmentLength+1,1:2)=1;
            s.strandBreaks(s.strandBreakClassification_segmentLength+s.strandBreakClassification_doubleStrandBreakSeparation,1:2)=1;
            strandBreakClassification = s.strandBreakClassification;
            assertEqual(N-1, strandBreakClassification(s.strandBreakClassification_index_NB));
            assertEqual(1, strandBreakClassification(s.strandBreakClassification_index_DSB_));
            assertEqual([0;1;N-1], unique(strandBreakClassification));
        end
        
        function testRM_MunI_Status(this)
            import edu.stanford.covert.util.SparseMat;
            
            s = this.state;
            
            m = struct;
            m.RM_MunI_RecognitionSites     = repmat((10:20:400)', 1, 6) + repmat(1:6, 20, 1);
            m.RM_MunI_MethylatedPositions  = [3 4];
            m.RM_MunI_RestrictionPositions = [1 5];
            
            nExamples = 15;
            nSites = size(m.RM_MunI_RecognitionSites,1);
            nChromosomes = s.nCompartments / 2;
            
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            exp_unmethylatedSites   = SparseMat([(nExamples+1:nSites)' ones(nSites-nExamples,1)],ones(nSites-nExamples,1),[nSites nChromosomes]);
            exp_hemimethylatedSites = SparseMat([],[],[nSites nChromosomes]);
            exp_methylatedSites     = SparseMat([],[],[nSites nChromosomes]);
            exp_cleavedSites        = SparseMat([],[],[nSites nChromosomes]);
            exp_damagedRegions      = SparseMat([(1:nSites)' repmat(2, nSites, 1)],true(nSites, 1),[nSites nChromosomes]);
            exp_sum                 = SparseMat([(1:nSites)' ones(nSites,1); (1:nSites)' repmat(2, nSites,1)], true(2*nSites, 1), [nSites nChromosomes]);
            
            %example 1
            exp_unmethylatedSites(1,1) = 1;
            
            %example 2
            s.damagedBases(m.RM_MunI_RecognitionSites(2,m.RM_MunI_MethylatedPositions(1)),1) = s.metabolite.m6ADIndexs;
            exp_hemimethylatedSites(2,1) = 1;
            
            %example 3
            s.damagedBases(m.RM_MunI_RecognitionSites(3,m.RM_MunI_MethylatedPositions(1)),2) = s.metabolite.m6ADIndexs;
            exp_damagedRegions(3,1) = 1;
            
            %example 4
            s.damagedBases(m.RM_MunI_RecognitionSites(4,m.RM_MunI_MethylatedPositions(2)),1) = s.metabolite.m6ADIndexs;
            exp_damagedRegions(4,1) = 1;
            
            %example 5
            s.damagedBases(m.RM_MunI_RecognitionSites(5,m.RM_MunI_MethylatedPositions(2)),2) = s.metabolite.m6ADIndexs;
            exp_hemimethylatedSites(5,1) = 1;
            
            %example 6
            s.damagedBases(m.RM_MunI_RecognitionSites(6,m.RM_MunI_MethylatedPositions(1)),1) = s.metabolite.m6ADIndexs;
            s.damagedBases(m.RM_MunI_RecognitionSites(6,m.RM_MunI_MethylatedPositions(2)),2) = s.metabolite.m6ADIndexs;
            exp_methylatedSites(6,1)=1;
            
            %example 7
            s.strandBreaks(m.RM_MunI_RecognitionSites(7,m.RM_MunI_RestrictionPositions(1)),1) = 1;
            s.strandBreaks(m.RM_MunI_RecognitionSites(7,m.RM_MunI_RestrictionPositions(2)),2) = 1;
            exp_cleavedSites(7,1)=1;
            
            %example 8
            s.strandBreaks(m.RM_MunI_RecognitionSites(8,m.RM_MunI_RestrictionPositions(1)),1) = 1;
            exp_damagedRegions(8,1)=1;
            
            %example 9
            s.strandBreaks(m.RM_MunI_RecognitionSites(9,m.RM_MunI_RestrictionPositions(2)),2) = 1;
            exp_damagedRegions(9,1)=1;
            
            %example 10
            s.damagedBases(m.RM_MunI_RecognitionSites(10,m.RM_MunI_MethylatedPositions(1)),1) = 1;
            exp_damagedRegions(10,1)=1;
            
            %example 11
            s.damagedBases(m.RM_MunI_RecognitionSites(11,m.RM_MunI_MethylatedPositions(1)-1),1) = s.metabolite.waterIndexs;
            exp_damagedRegions(11,1)=1;
            
            %example 12
            s.damagedBases(m.RM_MunI_RecognitionSites(12,m.RM_MunI_MethylatedPositions(1)),1) = 1;
            s.strandBreaks(m.RM_MunI_RecognitionSites(12,m.RM_MunI_RestrictionPositions(1)),1) = 1;
            s.strandBreaks(m.RM_MunI_RecognitionSites(12,m.RM_MunI_RestrictionPositions(2)),2) = 1;
            exp_damagedRegions(12,1)=1;
            
            %example 13
            s.gapSites(m.RM_MunI_RecognitionSites(13,1),1) = 1;
            exp_damagedRegions(13,1)=1;
            
            %example 14
            s.damagedBases(m.RM_MunI_RecognitionSites(14,m.RM_MunI_MethylatedPositions(1)),1) = 1;
            s.abasicSites(m.RM_MunI_RecognitionSites(14,1),1) = 1;
            exp_damagedRegions(14,1)=1;
            
            %example 15
            s.strandBreaks(m.RM_MunI_RecognitionSites(15,m.RM_MunI_RestrictionPositions(1)),1) = 1;
            s.strandBreaks(m.RM_MunI_RecognitionSites(15,m.RM_MunI_RestrictionPositions(2)),2) = 1;
            s.intrastrandCrossLinks(m.RM_MunI_RecognitionSites(15,1), 1) = 1;
            exp_damagedRegions(15,1)=1;
            
            %calculate status
            [unmethylatedSites, hemimethylatedSites, methylatedSites, cleavedSites, damagedRegions] = ...
                s.rmStatus(m.RM_MunI_RecognitionSites, m.RM_MunI_MethylatedPositions, m.RM_MunI_RestrictionPositions);
            
            %assertions
            assertEqual(exp_unmethylatedSites, unmethylatedSites);
            assertEqual(exp_hemimethylatedSites, hemimethylatedSites);
            assertEqual(exp_methylatedSites, methylatedSites);
            assertEqual(exp_cleavedSites, cleavedSites);
            assertEqual(exp_damagedRegions, damagedRegions);
            assertEqual(exp_sum, ...
                unmethylatedSites + ...
                hemimethylatedSites + ...
                methylatedSites + ...
                cleavedSites + ...
                damagedRegions);
        end
    end
    
    %mutators
    methods
        function testSetRegionUnwound(this)            
            import edu.stanford.covert.util.CircularSparseMat;
            
            s = this.state;
            L = size(s.sequence,1);
            
            %ex1: starting in forward direction from oriC
            s.initialize();
            k = s.linkingNumbers([1 1]);
            s.setRegionUnwound(1, 100);
            
            assertEqual([k; k], s.linkingNumbers([101 1; 101 2]));
            assertEqual(2, nnz(s.linkingNumbers));
            
            assertEqual(L, s.polymerizedRegions([1 1]));
            assertEqual(L-100, s.polymerizedRegions([101 2]));
            assertEqual(100, s.polymerizedRegions([1 4]));
            assertEqual(3, nnz(s.polymerizedRegions));
            
            assertEqual(2*L, collapse(s.polymerizedRegions));
            assertEqual(2*k, collapse(s.linkingNumbers));
            
            %ex2: continuing in forward direction
            s.setRegionUnwound(101, 100);
            
            assertEqual([k; k], s.linkingNumbers([201 1; 201 2]));
            assertEqual(2, nnz(s.linkingNumbers));
            
            assertEqual(L, s.polymerizedRegions([1 1]));
            assertEqual(L-200, s.polymerizedRegions([201 2]));
            assertEqual(200, s.polymerizedRegions([1 4]));
            assertEqual(3, nnz(s.polymerizedRegions));
            
            assertEqual(2*L, collapse(s.polymerizedRegions));
            assertEqual(2*k, collapse(s.linkingNumbers));
            
            %ex3: starting in backward direction from oriC
            s.initialize();
            k = s.linkingNumbers([1 1]);
            s.setRegionUnwound(L, -100);
            
            assertEqual([k; k], s.linkingNumbers([1 1; 1 2]));
            assertEqual(2, nnz(s.linkingNumbers));
            
            assertEqual(L, s.polymerizedRegions([1 1]));
            assertEqual(L-100, s.polymerizedRegions([1 2]));
            assertEqual(100, s.polymerizedRegions([L-100+1 4]));
            assertEqual(3, nnz(s.polymerizedRegions));
            
            assertEqual(2*L, collapse(s.polymerizedRegions));
            assertEqual(2*k, collapse(s.linkingNumbers));
            
            %ex4: continuing in backward direction
            s.setRegionUnwound(L-100, -100);
            
            assertEqual([k; k], s.linkingNumbers([1 1; 1 2]));
            assertEqual(2, nnz(s.linkingNumbers));
            
            assertEqual(L, s.polymerizedRegions([1 1]));
            assertEqual(L-200, s.polymerizedRegions([1 2]));
            assertEqual(200, s.polymerizedRegions([L-200+1 4]));
            assertEqual(3, nnz(s.polymerizedRegions));
            
            assertEqual(2*L, collapse(s.polymerizedRegions));
            assertEqual(2*k, collapse(s.linkingNumbers));

            %ex5: starting somewhere in the middle
            s.initialize();
            assertExceptionThrown(@() s.setRegionUnwound(34567, 100), 'Chromosome:invalidInput');

            %ex6: unwinding a single-stranded region in forward direction
            s.initialize();
            s.setRegionUnwound(1, 100);
            assertExceptionThrown(@() s.setRegionUnwound(50, 100), 'Chromosome:invalidInput');

            %ex7: unwinding a single-stranded region in reverse direction
            s.initialize();
            s.setRegionUnwound(L, -100);
            assertExceptionThrown(@() s.setRegionUnwound(L-50, -100), 'Chromosome:invalidInput');
        end
        
        function testSetRegionPolymerized(this)
            import edu.stanford.covert.cell.sim.SimulationStateSideEffect;
            import edu.stanford.covert.util.CircularSparseMat;

            s = this.state;
            
            %ex1: polymerize strand 1 and then strand 2 in forward direction
            s.initialize();
            L = size(s.sequence,1);
            R = s.linkingNumbers([1 1]);
            r = s.relaxedBasesPerTurn;
                        
            assertExceptionThrown(@() s.setRegionPolymerized([1 0], 10), 'Chromosome:invalidInput');
            assertExceptionThrown(@() s.setRegionPolymerized([1 3], 10), 'Chromosome:invalidInput');
            assertExceptionThrown(@() s.setRegionPolymerized([1 4], 10), 'Chromosome:invalidInput');
            assertExceptionThrown(@() s.setRegionPolymerized([1 5], 10), 'Chromosome:invalidInput');
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionUnwound(1, 0));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([1 1], 0));
            assertEqual(CircularSparseMat([1 1; 1 2], [L;L], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 1 2], [R;R], [L s.nCompartments], 1), s.linkingNumbers);
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionUnwound(1, 10));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([1 1], 5));
            assertEqual(CircularSparseMat([1 1; 1 2; 11 2; 1 4], [L; 5; L-10; 10], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 11 1; 1 2; 11 2], [5/r; R; 5/r; R], [L s.nCompartments], 1), s.linkingNumbers);
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([1 2], 5));
            assertEqual(CircularSparseMat([1 1; 1 2; 11 2; 1 3; 1 4], [L; 5; L-10; 5; 10], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 11 1; 1 2; 11 2; 1 3; 1 4], [5/r; R; 5/r; R; 5/r; 5/r], [L s.nCompartments], 1), s.linkingNumbers);
            
            assertExceptionThrown(@() s.setRegionUnwound(1, 10), 'Chromosome:invalidInput');
            assertExceptionThrown(@() s.setRegionUnwound(9, 10), 'Chromosome:invalidInput');
            assertExceptionThrown(@() s.setRegionPolymerized([1 1], 10), 'Chromosome:error');
            assertExceptionThrown(@() s.setRegionPolymerized([1 2], 10), 'Chromosome:error');
            assertExceptionThrown(@() s.setRegionPolymerized([9 1], 10), 'Chromosome:error');
            assertExceptionThrown(@() s.setRegionPolymerized([9 2], 10), 'Chromosome:error');
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionUnwound(11, 10));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([6 2], 10));
            assertEqual(CircularSparseMat([1 1; 1 2; 21 2; 1 3; 1 4], [L; 5; L-20; 15; 20], [L s.nCompartments], 1), s.polymerizedRegions);
            assertElementsAlmostEqual(CircularSparseMat([1 1; 21 1; 1 2; 21 2; 1 3; 1 4], [5/r; R; 5/r; R; 15/r; 15/r], [L s.nCompartments], 1), s.linkingNumbers);
            
            %ex2: polymerize strand 2 and then strand 1 in forward direction
            s.initialize();
            assertEqual(CircularSparseMat([1 1; 1 2], [L;L], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 1 2], [R;R], [L s.nCompartments], 1), s.linkingNumbers);
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionUnwound(1, 10));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([1 2], 5));
            assertEqual(CircularSparseMat([1 1; 1 3; 11 2; 1 4], [L; 5; L-10; 10], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 3; 11 1; 1 4; 11 2], [5/r; R; 5/r; R], [L s.nCompartments], 1), s.linkingNumbers);
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([1 1], 5));
            assertEqual(CircularSparseMat([1 1; 1 2; 11 2; 1 3; 1 4], [L; 5; L-10; 5; 10], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 11 1; 1 2; 11 2; 1 3; 1 4], [5/r; R; 5/r; R; 5/r; 5/r], [L s.nCompartments], 1), s.linkingNumbers);
            
            %ex3: polymerize strand 1 and then strand 2 in reverse direction
            s.initialize();

            assertExceptionThrown(@() s.setRegionUnwound(L-10, 20), 'Chromosome:invalidInput');
            assertExceptionThrown(@() s.setRegionPolymerized([L-10 1], 20), 'Chromosome:invalidInput');
            assertExceptionThrown(@() s.setRegionPolymerized([L-10 2], 20), 'Chromosome:invalidInput');

            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionUnwound(L, -11));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([L 2], -6));
            assertEqual(CircularSparseMat([1 1; 1 2; L-5 3; L-10 4], [L; L-11; 6; 11], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 1 2; L-5 3; L-5 4], [R; R; 6/r; 6/r], [L s.nCompartments], 1), s.linkingNumbers);

            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionUnwound(L-11, -20));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([L-6 2], -20));
            assertEqual(CircularSparseMat([1 1; 1 2; L-25 3; L-30 4], [L; L-31; 26; 31], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 1 2; L-25 3; L-25 4], [R; R; 26/r; 26/r], [L s.nCompartments], 1), s.linkingNumbers);
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([L 1], -26));
            assertEqual(CircularSparseMat([1 1; 1 2; L-25 2; L-25 3; L-30 4], [L; L-31; 26; 26; 31], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 1 2; L-25 1; L-25 2; L-25 3; L-25 4], [R; R; 26/r; 26/r; 26/r; 26/r], [L s.nCompartments], 1), s.linkingNumbers);
            
            %ex4: polymerize strand 2 and then strand 1 in reverse direction
            s.initialize();
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionUnwound(L, -11));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([L 1], -6));
            assertEqual(CircularSparseMat([1 1; 1 2; L-5 2; L-10 4], [L; L-11; 6; 11], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 1 2; L-5 1; L-5 2], [R; R; 6/r; 6/r], [L s.nCompartments], 1), s.linkingNumbers);
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionUnwound(L-11, -20));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([L-6 1], -20));
            assertEqual(CircularSparseMat([1 1; 1 2; L-25 2; L-30 4], [L; L-31; 26; 31], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 1 2; L-25 1; L-25 2], [R; R; 26/r; 26/r], [L s.nCompartments], 1), s.linkingNumbers);
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([L 2], -26));
            assertEqual(CircularSparseMat([1 1; 1 2; L-25 2; L-25 3; L-30 4], [L; L-31; 26; 26; 31], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 1 2; L-25 1; L-25 2; L-25 3; L-25 4], [R; R; 26/r; 26/r; 26/r; 26/r], [L s.nCompartments], 1), s.linkingNumbers);            
            
            %ex5: polymerize both strand 1 and strand 2 in one method call (strand 1 then strand 2)
            s.initialize();
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionUnwound(1, 10));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([1 1; 1 2], [5; 5]));
            assertEqual(CircularSparseMat([1 1; 1 2; 11 2; 1 3; 1 4], [L; 5; L-10; 5; 10], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 11 1; 1 2; 11 2; 1 3; 1 4], [5/r; R; 5/r; R; 5/r; 5/r], [L s.nCompartments], 1), s.linkingNumbers);
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionUnwound(11, 10));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([6 1; 6 2], [10; 10]));
            assertEqual(CircularSparseMat([1 1; 1 2; 21 2; 1 3; 1 4], [L; 15; L-20; 15; 20], [L s.nCompartments], 1), s.polymerizedRegions);
            assertElementsAlmostEqual(CircularSparseMat([1 1; 21 1; 1 2; 21 2; 1 3; 1 4], [15/r; R; 15/r; R; 15/r; 15/r], [L s.nCompartments], 1), s.linkingNumbers);
            
            %ex6: polymerize both strand 1 and strand 2 in one method call (strand 2 then strand 1)
            s.initialize();
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionUnwound(1, 10));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([1 2; 1 1], [5; 5]));
            assertEqual(CircularSparseMat([1 1; 1 2; 11 2; 1 3; 1 4], [L; 5; L-10; 5; 10], [L s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 11 1; 1 2; 11 2; 1 3; 1 4], [5/r; R; 5/r; R; 5/r; 5/r], [L s.nCompartments], 1), s.linkingNumbers);
            
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionUnwound(11, 10));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), s.setRegionPolymerized([6 2; 6 1], [10; 10]));
            assertEqual(CircularSparseMat([1 1; 1 2; 21 2; 1 3; 1 4], [L; 15; L-20; 15; 20], [L s.nCompartments], 1), s.polymerizedRegions);
            assertElementsAlmostEqual(CircularSparseMat([1 1; 21 1; 1 2; 21 2; 1 3; 1 4], [15/r; R; 15/r; R; 15/r; 15/r], [L s.nCompartments], 1), s.linkingNumbers);            
        end
        
        function testSetSiteProteinBound(this)
            import edu.stanford.covert.cell.sim.SimulationStateSideEffect;
            import edu.stanford.covert.cell.sim.SimulationStateSideEffectItem;
            import edu.stanford.covert.util.CircularSparseMat;
            
            s = this.state;
            
            %% ex 1
            s.monomerDNAFootprints = [1; 2; 3; 1; 2; 2];
            s.complexDNAFootprints = repmat(3, 6, 1);
            s.monomerDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, 6, 1);
            s.complexDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, 6, 1);
            s.monomerDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, 6, 1);
            s.complexDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, 6, 1);
            
            s.reactionMonomerCatalysisMatrix = zeros(0, 6);
            s.reactionComplexCatalysisMatrix = zeros(0, 6);
            s.reactionBoundMonomer = zeros(0, 1);
            s.reactionBoundComplex = zeros(0, 1);
            s.reactionThresholds = ...
                sum(s.reactionMonomerCatalysisMatrix, 2) + ...
                sum(s.reactionComplexCatalysisMatrix, 2);
            
            %ex 1.1
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.monomerBoundSites(100, 1) = 2;
            s.monomerBoundSites(103, 2) = 2;
            s.monomerBoundSites(89, 1) = 2;
            s.monomerBoundSites(91, 2) = 3;
            s.monomerBoundSites(105, 2) = 5;
            s.monomerBoundSites(107, 2) = 5;
            s.monomerBoundSites(109, 1) = 6;
            monomerBoundSites = s.monomerBoundSites;
            monomerBoundSites(1, 1) = 3;
            monomerBoundSites(120, 1) = 3;
            monomerBoundSites(130, 2) = 3;
            
            [releasedMonomers, releasedComplexs, sideEffects, tfs, idxs, positionsStrands] = ...
                s.setSiteProteinBound([1 1; 91 1; 101 1; 104 1; 107 2; 120 1; 130 2], 3, [], 3, [], (1:4)', (1:4)', true, false, 1, false, []);

            assertEqual(monomerBoundSites, s.monomerBoundSites);
            assertEqual([0; 0; -3; 0], releasedMonomers);
            assertEqual([0; 0; 0; 0], releasedComplexs);
            assertEqual(SimulationStateSideEffect.empty(0, 1), sideEffects);
            assertEqual([true; false(4,1); true(2,1)], tfs);
            assertEqual([1; 6; 7], idxs);
            assertEqual([1 1; 120 1; 130 2], positionsStrands);
            
            %ex 1.2
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.monomerBoundSites(100, 1) = 2;
            s.monomerBoundSites(103, 2) = 2;
            s.monomerBoundSites(89, 1) = 2;
            s.monomerBoundSites(91, 2) = 3;
            s.monomerBoundSites(105, 2) = 5;
            s.monomerBoundSites(107, 2) = 5;
            s.monomerBoundSites(109, 1) = 6;
            monomerBoundSites = s.monomerBoundSites;
            monomerBoundSites(1, 1) = 3;
            monomerBoundSites(120, 1) = 3;
            monomerBoundSites(130, 2) = 3;
            
            [releasedMonomers, releasedComplexs, sideEffects, tfs, idxs, positionsStrands] = ...
                s.setSiteProteinBound([1 1; 91 1; 101 1; 104 1; 107 2; 120 1; 130 2], 10, [], 3, [], (1:4)', (1:4)', true, false, 1, false, []);  
            
            assertEqual(monomerBoundSites, s.monomerBoundSites);
            assertEqual([0; 0; -3; 0], releasedMonomers);
            assertEqual([0; 0; 0; 0], releasedComplexs);
            assertEqual(SimulationStateSideEffect.empty(0, 1), sideEffects);
            assertEqual([true; false(4,1); true(2,1)], tfs);
            assertEqual([1; 6; 7], idxs);
            assertEqual([1 1; 120 1; 130 2], positionsStrands);
            
            %ex 1.3
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.monomerBoundSites(100, 1) = 2;
            s.monomerBoundSites(103, 2) = 2;
            s.monomerBoundSites(89, 1) = 2;
            s.monomerBoundSites(91, 2) = 3;
            s.monomerBoundSites(105, 2) = 5;
            s.monomerBoundSites(107, 2) = 5;
            s.monomerBoundSites(109, 1) = 6;
            monomerBoundSites = s.monomerBoundSites;
            monomerBoundSites(1-1, 1) = 3;
            monomerBoundSites(120-1, 1) = 3;
            monomerBoundSites(130-1, 2) = 3;
            
            [releasedMonomers, releasedComplexs, sideEffects, tfs, idxs, positionsStrands] = ...
                s.setSiteProteinBound([1 1; 91 1; 101 1; 104 1; 107 2; 120 1; 130 2], 3, [], 3, [], (1:4)', (1:4)', true, true, 1, false, []);
            
            assertEqual(monomerBoundSites, s.monomerBoundSites);
            assertEqual([0; 0; -3; 0], releasedMonomers);
            assertEqual([0; 0; 0; 0], releasedComplexs);
            assertEqual(SimulationStateSideEffect.empty(0, 1), sideEffects);
            assertEqual([true; false(4,1); true(2,1)], tfs);
            assertEqual([1; 6; 7], idxs);
            assertEqual([1 1; 120 1; 130 2], positionsStrands);
            
            %ex 1.4
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.monomerBoundSites(100, 1) = 2;
            s.monomerBoundSites(103, 2) = 2;
            s.monomerBoundSites(89, 1) = 2;
            s.monomerBoundSites(91, 2) = 3;
            s.monomerBoundSites(105, 2) = 5;
            s.monomerBoundSites(107, 2) = 5;
            s.monomerBoundSites(109, 1) = 6;
            
            [releasedMonomers, releasedComplexs, sideEffects, tfs, idxs, positionsStrands] = ...
                s.setSiteProteinBound([1 1; 91 1; 101 1; 104 1; 107 2; 120 1; 130 1], 2, [], 3, [], (1:4)', (1:4)', true, false, 1, false, []);
            
            assertEqual([0; 0; -2; 0], releasedMonomers);
            assertEqual([0; 0; 0; 0], releasedComplexs);
            assertEqual(SimulationStateSideEffect.empty(0, 1), sideEffects);
            assertEqual(2, sum(tfs));
            assertEqual(2, numel(idxs));
            assertEqual(2, size(positionsStrands, 2));
            
            %% ex 2
            s.monomerDNAFootprints = [1; 2; 3; 1; 3; 2];
            s.complexDNAFootprints = repmat(3, 6, 1);
            s.monomerDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, 6, 1);
            s.complexDNAFootprintBindingStrandedness = repmat(s.dnaStrandedness_dsDNA, 6, 1);
            s.monomerDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, 6, 1);
            s.complexDNAFootprintRegionStrandedness = repmat(s.dnaStrandedness_dsDNA, 6, 1);
            
            s.reactionMonomerCatalysisMatrix = zeros(2, 6);
            s.reactionComplexCatalysisMatrix = zeros(2, 6);
            s.reactionComplexCatalysisMatrix(:, 3) = 1;
            s.reactionBoundMonomer = [3; 5];
            s.reactionBoundComplex = [0; 0];
            s.reactionThresholds = ...
                sum(s.reactionMonomerCatalysisMatrix, 2) + ...
                sum(s.reactionComplexCatalysisMatrix, 2);
            
            %ex 2.1
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.monomerBoundSites(100, 1) = 2;
            s.monomerBoundSites(102, 2) = 2;
            s.monomerBoundSites(89, 1) = 2;
            s.monomerBoundSites(91, 2) = 3;
            s.monomerBoundSites(104, 2) = 5;
            s.monomerBoundSites(107, 2) = 5;
            s.monomerBoundSites(110, 1) = 6;
            
            [releasedMonomers, releasedComplexs, sideEffects, tfs, idxs, positionsStrands] = ...
                s.setSiteProteinBound([1 1; 91 1; 101 1; 105 1; 108 2; 120 1; 130 2], 10, [], [], 3, (1:4)', (1:4)', true, false, 1, false, []);
            
            assertEqual(CircularSparseMat([89 1; 100 1; 110 1; 102 2], [2; 2; 6; 2], [size(s.sequence, 1) s.nCompartments], 1), s.monomerBoundSites);
            assertEqual(CircularSparseMat([1 1; 91 1; 105 1; 120 1; 130 2], [3; 3; 3; 3; 3], [size(s.sequence, 1) s.nCompartments], 1), s.complexBoundSites);
            assertEqual([0; 0; 1; 0], releasedMonomers);
            assertEqual([0; 0; -5; 0], releasedComplexs);
            assertEqual(SimulationStateSideEffect([...
                SimulationStateSideEffectItem('ProteinMonomer','counts','matureIndexs', 5, 1,  2); ...
                SimulationStateSideEffectItem('ProteinMonomer','counts','boundIndexs',  5, 1, -2)]),...
                sideEffects);
            assertEqual([true; true; false; true; false; true(2,1)], tfs);
            assertEqual([1; 2; 4; 6; 7], idxs);
            assertEqual([1 1; 91 1; 105 1; 120 1; 130 2], positionsStrands);
            
            %ex 2.2
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.monomerBoundSites(100, 1) = 2;
            s.monomerBoundSites(102, 2) = 2;
            s.monomerBoundSites(89, 1) = 2;
            s.monomerBoundSites(91, 2) = 3;
            s.monomerBoundSites(105, 2) = 5;
            s.monomerBoundSites(107, 2) = 5;
            s.monomerBoundSites(110, 1) = 6;
            
            [releasedMonomers, releasedComplexs, sideEffects, tfs, idxs, positionsStrands] = ...
                s.setSiteProteinBound([1 1; 91 1; 101 1; 104 1; 108 2; 120 1; 130 2], 10, [], [], 3, (1:4)', (1:4)', true, false, 1, false, []);
            
            assertEqual(CircularSparseMat([89 1; 100 1; 110 1; 102 2; 107 2], [2; 2; 6; 2; 5], [size(s.sequence, 1) s.nCompartments], 1), s.monomerBoundSites);
            assertEqual(CircularSparseMat([1 1; 91 1; 104 1; 120 1; 130 2], [3; 3; 3; 3; 3], [size(s.sequence, 1) s.nCompartments], 1), s.complexBoundSites);
            assertEqual([0; 0; 1; 0], releasedMonomers);
            assertEqual([0; 0; -5; 0], releasedComplexs);
            assertEqual(SimulationStateSideEffect([...
                SimulationStateSideEffectItem('ProteinMonomer','counts','matureIndexs', 5, 1,  1); ...
                SimulationStateSideEffectItem('ProteinMonomer','counts','boundIndexs',  5, 1, -1)]),...
                sideEffects);
            assertEqual([true; true; false; true; false; true(2,1)], tfs);
            assertEqual([1; 2; 4; 6; 7], idxs);
            assertEqual([1 1; 91 1; 104 1; 120 1; 130 2], positionsStrands);
            
            %ex 2.3
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.monomerBoundSites(100, 1) = 2;
            s.monomerBoundSites(102, 2) = 2;
            s.monomerBoundSites(89, 1) = 2;
            s.monomerBoundSites(91, 2) = 3;
            s.monomerBoundSites(105, 2) = 5;
            s.monomerBoundSites(107, 2) = 5;
            s.monomerBoundSites(110, 1) = 6;
            
            [releasedMonomers, releasedComplexs, sideEffects, tfs, idxs, positionsStrands] = ...
                s.setSiteProteinBound([1 1; 90 1; 101 1; 104 1; 108 2; 120 1; 130 2], 10, [], [], 3, (1:4)', (1:4)', true, false, 1, false, []);
            
            assertEqual(CircularSparseMat([89 1; 100 1; 110 1; 91 2; 102 2; 107 2], [2; 2; 6; 3; 2; 5], [size(s.sequence, 1) s.nCompartments], 1), s.monomerBoundSites);
            assertEqual(CircularSparseMat([1 1; 104 1; 120 1; 130 2], [3; 3; 3; 3], [size(s.sequence, 1) s.nCompartments], 1), s.complexBoundSites);
            assertEqual([0; 0; 0; 0], releasedMonomers);
            assertEqual([0; 0; -4; 0], releasedComplexs);
            assertEqual(SimulationStateSideEffect([...
                SimulationStateSideEffectItem('ProteinMonomer','counts','matureIndexs', 5, 1,  1); ...
                SimulationStateSideEffectItem('ProteinMonomer','counts','boundIndexs',  5, 1, -1)]),...
                sideEffects);
            assertEqual([true; false; false; true; false; true(2,1)], tfs);
            assertEqual([1; 4; 6; 7], idxs);
            assertEqual([1 1; 104 1; 120 1; 130 2], positionsStrands);
        end
        
        function testSetRegionProteinUnbound(this)
            import edu.stanford.covert.cell.sim.SimulationStateSideEffect;
            import edu.stanford.covert.cell.sim.SimulationStateSideEffectItem;
            import edu.stanford.covert.util.CircularSparseMat;
            
            s = this.state;
            
            %example 1: monomers bound
            s.monomerDNAFootprints = ones(6, 1);
            s.complexDNAFootprints = ones(6, 1);
            
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.monomerBoundSites(100, 1) = 2;
            s.monomerBoundSites(103, 2) = 2;
            s.monomerBoundSites(90, 1) = 2;
            s.monomerBoundSites(91, 2) = 3;
            s.monomerBoundSites(104, 2) = 5;
            s.monomerBoundSites(105, 2) = 5;
            s.monomerBoundSites(106, 1) = 6;
            
            [releasedMonomers, releasedComplexs, sideEffects] = s.setRegionProteinUnbound([91 2; 100 1], 5, [1; 2; 3; 4], [], true, true, false, false);
            
            assertEqual([0; 2; 1; 0], releasedMonomers);
            assertEqual(zeros(0, 1), releasedComplexs);
            assertEqual(...
                SimulationStateSideEffect([...
                SimulationStateSideEffectItem('ProteinMonomer','counts','matureIndexs', 5, 1,  1); ...
                SimulationStateSideEffectItem('ProteinMonomer','counts','boundIndexs',  5, 1, -1)]),...
                sideEffects);
            assertEqual(CircularSparseMat([1 1; 1 2], [size(s.sequence, 1); size(s.sequence, 1)], [size(s.sequence,1) s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([90 1; 106 1; 105 2], [2; 6; 5], [size(s.sequence,1) s.nCompartments], 1), s.monomerBoundSites);
            assertEqual(CircularSparseMat([], [], [size(s.sequence,1) s.nCompartments], 1), s.complexBoundSites);
            assertEqual(CircularSparseMat([], false(0, 1), [size(s.sequence,1) s.nCompartments], 1), s.gapSites);
            assertEqual(CircularSparseMat([], false(0, 1), [size(s.sequence,1) s.nCompartments], 1), s.abasicSites);
            assertEqual(CircularSparseMat([], [],          [size(s.sequence,1) s.nCompartments], 1), s.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [],          [size(s.sequence,1) s.nCompartments], 1), s.damagedBases);
            assertEqual(CircularSparseMat([], [],          [size(s.sequence,1) s.nCompartments], 1), s.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], false(0, 1), [size(s.sequence,1) s.nCompartments], 1), s.strandBreaks);
            assertEqual(CircularSparseMat([], false(0, 1), [size(s.sequence,1) s.nCompartments], 1), s.hollidayJunctions);
            
            %example 2: complexes bound
            s.monomerDNAFootprints = ones(6, 1);
            s.complexDNAFootprints = repmat(2, 6, 1);
            
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            s.complexBoundSites(100, 1) = 2;
            s.complexBoundSites(103, 2) = 2;
            s.complexBoundSites(90, 1) = 2;
            s.complexBoundSites(91, 2) = 3;
            s.complexBoundSites(104, 2) = 5;
            s.complexBoundSites(105, 2) = 5;
            s.complexBoundSites(106, 1) = 6;
            
            [releasedMonomers, releasedComplexs, sideEffects] = s.setRegionProteinUnbound([91 2; 100 1], 5, [], [1; 2; 3; 4], true, true, false, false);
            
            assertEqual([0; 3; 1; 0], releasedComplexs);
            assertEqual(zeros(0, 1), releasedMonomers);
            assertEqual(...
                SimulationStateSideEffect([...
                SimulationStateSideEffectItem('ProteinComplex','counts','matureIndexs', 5, 1,  1); ...
                SimulationStateSideEffectItem('ProteinComplex','counts','boundIndexs',  5, 1, -1)]),...
                sideEffects);
            assertEqual(CircularSparseMat([1 1; 1 2], [size(s.sequence, 1); size(s.sequence, 1)], [size(s.sequence,1) s.nCompartments], 1), s.polymerizedRegions);
            assertEqual(CircularSparseMat([], [], [size(s.sequence,1) s.nCompartments], 1), s.monomerBoundSites);
            assertEqual(CircularSparseMat([106 1; 105 2], [6; 5], [size(s.sequence,1) s.nCompartments], 1), s.complexBoundSites);
            assertEqual(CircularSparseMat([], false(0, 1), [size(s.sequence,1) s.nCompartments], 1), s.gapSites);
            assertEqual(CircularSparseMat([], false(0, 1), [size(s.sequence,1) s.nCompartments], 1), s.abasicSites);
            assertEqual(CircularSparseMat([], [],          [size(s.sequence,1) s.nCompartments], 1), s.damagedSugarPhosphates);
            assertEqual(CircularSparseMat([], [],          [size(s.sequence,1) s.nCompartments], 1), s.damagedBases);
            assertEqual(CircularSparseMat([], [],          [size(s.sequence,1) s.nCompartments], 1), s.intrastrandCrossLinks);
            assertEqual(CircularSparseMat([], false(0, 1), [size(s.sequence,1) s.nCompartments], 1), s.strandBreaks);
            assertEqual(CircularSparseMat([], false(0, 1), [size(s.sequence,1) s.nCompartments], 1), s.hollidayJunctions);
        end
        
        function testSetSiteDamaged(this)
            import edu.stanford.covert.cell.sim.SimulationStateSideEffect;
            import edu.stanford.covert.util.CircularSparseMat;
            
            s = this.state;
            
            %ex 1: any site can be damaged
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            
            [positionsStrands, sideEffects] = s.setSiteDamaged('damagedSugarPhosphates', 1, 0.5, 10, '', '');
            assertEqual([10 2], size(positionsStrands));
            assertEqual(CircularSparseMat(positionsStrands, 1, [size(s.sequence, 1) s.nCompartments], 1), s.damagedSugarPhosphates);
            assertEqual(s.damagedSugarPhosphates, s.getDamagedSites(true, true, true, true, false, false, true));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), sideEffects);
            
            %ex 2: any site with particular sequence can be damaged
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            
            [positionsStrands, sideEffects] = s.setSiteDamaged('damagedBases', 1, 0.5, 10, 'A', '');
            assertEqual([10 2], size(positionsStrands));
            assertEqual('A', unique(s.sequence.subsequence(positionsStrands)));
            assertEqual(CircularSparseMat(positionsStrands, 1, [size(s.sequence, 1) s.nCompartments], 1), s.damagedBases);
            assertEqual(s.damagedBases, s.getDamagedSites(true, true, true, true, false, false, true));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), sideEffects);
            
            %ex 3: only damaged sites can be damaged
            [positionsStrands2, sideEffects] = s.setSiteDamaged('damagedBases', 2, 1, 20, 1, 'damagedBases');
            assertEqual(positionsStrands, positionsStrands2);
            assertEqual(CircularSparseMat(positionsStrands, 2, [size(s.sequence, 1) s.nCompartments], 1), s.damagedBases);
            assertEqual(s.damagedBases, s.getDamagedSites(true, true, true, true, false, false, true));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), sideEffects);
            
            %ex 4: zero probability damage
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            
            [positionsStrands, sideEffects] = s.setSiteDamaged('damagedBases', 2, 0, 20, 1, '');
            assertEqual(zeros(0, 2), positionsStrands);
            assertEqual(CircularSparseMat([], [], [size(s.sequence, 1) s.nCompartments], 1), s.getDamagedSites(true, true, true, true, false, false, true));
            assertEqual(repmat(SimulationStateSideEffect, 0, 1), sideEffects);
        end
    end
    
    %private mutators
    methods
        function testMergeAdjacentRegions(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            s = this.state;
            
            %ex 1: nothing polymerized
            s.allocateMemory(1);
            assertEqual(CircularSparseMat([], [], [size(s.sequence,1) s.nCompartments], 1), s.mergeAdjacentRegions(s.polymerizedRegions));
            
            %ex 2: one thing polymerized
            s.allocateMemory(1);
            s.polymerizedRegions(1,1) = 50;
            assertEqual(CircularSparseMat([1 1], 50, [size(s.sequence,1) s.nCompartments], 1), s.mergeAdjacentRegions(s.polymerizedRegions));
            
            %ex 3: two things polymerized, same strand, not contiguous
            s.allocateMemory(1);
            s.polymerizedRegions(1,1) = 50;
            s.polymerizedRegions(100, 1) = 50;
            assertEqual(CircularSparseMat([1 1; 100 1], [50; 50], [size(s.sequence,1) s.nCompartments], 1), s.mergeAdjacentRegions(s.polymerizedRegions));
            
            %ex 4: two things polymerized, different strands, contiguous
            s.allocateMemory(1);
            s.polymerizedRegions(1,1) = 50;
            s.polymerizedRegions(51,2) = 50;
            assertEqual(CircularSparseMat([1 1; 51 2], [50;50], [size(s.sequence,1) s.nCompartments], 1), s.mergeAdjacentRegions(s.polymerizedRegions));
            
            %ex 5: two things polymerized, contiguous
            s.allocateMemory(1);
            s.polymerizedRegions(1,1) = 50;
            s.polymerizedRegions(51,1) = 50;
            assertEqual(CircularSparseMat([1 1], 100, [size(s.sequence,1) s.nCompartments], 1), s.mergeAdjacentRegions(s.polymerizedRegions));
            
            %ex 6: multiple polymerized regions, some contiguous, some not;
            %some wrapping over ORI
            s.allocateMemory(1);
            s.polymerizedRegions(end, 1) = 1;
            s.polymerizedRegions(end, 3) = 1;
            s.polymerizedRegions(1, 1) = 50;
            s.polymerizedRegions(51, 1) = 50;
            s.polymerizedRegions(1, 2) = 50;
            s.polymerizedRegions(51, 2) = 75;
            s.polymerizedRegions(102, 1) = 40;
            s.polymerizedRegions(142, 1) = 30;
            s.polymerizedRegions(172, 1) = 100;
            s.polymerizedRegions(272, 1) = 50;
            s.polymerizedRegions(1, 3) = 104;
            s.polymerizedRegions(150, 3) = 100;
            s.polymerizedRegions(150, 4) = 100;
            assertEqual(CircularSparseMat([1 1; 102 1; size(s.sequence,1) 1; 1 2; 1 3; size(s.sequence,1) 3; 150 3; 150 4], ...
                [100; 220; 1; 125; 104; 1; 100; 100], ...
                [size(s.sequence,1) s.nCompartments], 1), s.mergeAdjacentRegions(s.polymerizedRegions));
        end

        function testExcludeRegions(this)
            s = this.state;
            L = length(s.sequence);

            [ps, len] = s.excludeRegions([1 1], L, [400000 1], 140);
            assertEqual([400140 1], ps);
            assertEqual(L - 140, len);
        end
        
        function testIntersectRegions(this)
            s = this.state;

            [psC, lenC] = s.intersectRegions(...  %  1234567890123
                [1 1; 9 1], [5; 4],...            %A -----   ---- 
                [4 1; 12 1], [6; 2]);             %B    ------  --
            assertEqual([4 1; 9 1; 12 1], psC);   %C    --   -  - 
            assertEqual([  2;   1;    1], lenC);

            L = length(s.sequence);
            [psC, lenC] = s.intersectRegions(...
                [L-1 3; 7 2], [6; 4],...
                [L-2 3; 1 3; 5 2; 1 1], [3; 3; 4; 99]);
            assertEqual([7 2; L-1 3], psC);
            assertEqual([  2;     5], lenC);
        end

        function testJoinSplitRegions(this)
            c = this.state;
            
            %ex 1: first chromosome completely synthesized
            c.initialize();
                        
            [polRgnPosChrs, polRgnLens] = find(c.polymerizedRegions(:, 1:2:end));
            [polRgnPosChrs, polRgnLens] = c.joinSplitRegions(polRgnPosChrs, polRgnLens);
            assertEqual(1, polRgnPosChrs(:, 1));
            assertEqual(size(c.sequence, 1), polRgnPosChrs(:, 1) + polRgnLens - 1);
            assertEqual(1, polRgnPosChrs(:, 2));
            
            %ex 2: 
            c.initialize();
            c.setRegionUnwound(1, 1000);
            c.setRegionPolymerized([1 1], 1000);
            
            [polRgnPosChrs, polRgnLens] = find(c.polymerizedRegions(:, 1:2:end));
            [polRgnPosChrs, polRgnLens] = c.joinSplitRegions(polRgnPosChrs, polRgnLens);
            assertEqual(1, polRgnPosChrs(:, 1));
            assertEqual(size(c.sequence, 1), polRgnPosChrs(:, 1) + polRgnLens - 1);
            assertEqual(1, polRgnPosChrs(:, 2));
            
            %ex 3: 
            c.initialize();
            c.setRegionUnwound(1, 1000);
            c.setRegionPolymerized([1 2], 1000);
            
            [polRgnPosChrs, polRgnLens] = find(c.polymerizedRegions(:, 1:2:end));
            [polRgnPosChrs, polRgnLens] = c.joinSplitRegions(polRgnPosChrs, polRgnLens);
            assertEqual([1; 1], polRgnPosChrs(:, 1));
            assertEqual([size(c.sequence, 1); 1000], polRgnPosChrs(:, 1) + polRgnLens - 1);
            assertEqual([1; 2], polRgnPosChrs(:, 2));
        end
        
        function testJoinSplitOverOriCRegions(this)
            c = this.state;
            
            %ex 1: first chromosome completely synthesized
            c.initialize();
                        
            [polRgnPosChrs, polRgnLens] = find(c.polymerizedRegions(:, 1:2:end));
            [polRgnPosChrs, polRgnLens] = c.joinSplitOverOriCRegions(polRgnPosChrs, polRgnLens);
            assertEqual(1, polRgnPosChrs(:, 1));
            assertEqual(size(c.sequence, 1), polRgnPosChrs(:, 1) + polRgnLens - 1);
            assertEqual(1, polRgnPosChrs(:, 2));
            
            %ex 2: 
            c.initialize();
            c.setRegionUnwound(1, 1000);
            c.setRegionPolymerized([1 1], 1000);
            
            [polRgnPosChrs, polRgnLens] = find(c.polymerizedRegions(:, 1:2:end));
            [polRgnPosChrs, polRgnLens] = c.joinSplitOverOriCRegions(polRgnPosChrs, polRgnLens);
            assertEqual(1, polRgnPosChrs(:, 1));
            assertEqual(size(c.sequence, 1), polRgnPosChrs(:, 1) + polRgnLens - 1);
            assertEqual(1, polRgnPosChrs(:, 2));
            
            %ex 3: 
            c.initialize();
            c.setRegionUnwound(1, 1000);
            c.setRegionPolymerized([1 2], 1000);
            
            [polRgnPosChrs, polRgnLens] = find(c.polymerizedRegions(:, 1:2:end));
            [polRgnPosChrs, polRgnLens] = c.joinSplitOverOriCRegions(polRgnPosChrs, polRgnLens);
            assertEqual([1; 1], polRgnPosChrs(:, 1));
            assertEqual([size(c.sequence, 1); 1000], polRgnPosChrs(:, 1) + polRgnLens - 1);
            assertEqual([1; 2], polRgnPosChrs(:, 2));
            
            %ex 4: 
            c.initialize();
            c.setRegionUnwound(1, 1000);
            c.setRegionUnwound(size(c.sequence, 1), -1000);
            c.setRegionPolymerized([1 2], 1000);            
            c.setRegionPolymerized([size(c.sequence, 1) 2], -1000);
            
            [polRgnPosChrs, polRgnLens] = find(c.polymerizedRegions(:, 1:2:end));
            [polRgnPosChrs, polRgnLens] = c.joinSplitOverOriCRegions(polRgnPosChrs, polRgnLens);
            assertEqual([1; size(c.sequence, 1) - 1000 + 1], polRgnPosChrs(:, 1));
            assertEqual([size(c.sequence, 1); size(c.sequence, 1)+1000], polRgnPosChrs(:, 1) + polRgnLens - 1);
            assertEqual([1; 2], polRgnPosChrs(:, 2));
        end
    end
    
    %getters
    methods
        function testDependentProperties(this)
            s = this.state;
            s.invalidate();
            
            this.testDependentProperties@edu.stanford.covert.cell.sim.CellStateTestCase();            
        end
        
        function testShift35Getters(this)
            s = this.state;
            
            s.allocateMemory(1);
            s.intrastrandCrossLinks(1,1:2) = 1;
            [subs, vals] = find(s.intrastrandCrossLinks3);
            assertEqual([size(s.sequence,1) 1;2 2], subs);
            assertEqual([1; 1], vals);
            
            s.allocateMemory(1);
            s.intrastrandCrossLinks(2,1:2)=1;
            [subs,vals]=find(s.intrastrandCrossLinks3);
            assertEqual([1 1;3 2], subs);
            assertEqual([1; 1], vals);
            
            s.allocateMemory(1);
            s.intrastrandCrossLinks(end,1:2)=1;
            [subs,vals]=find(s.intrastrandCrossLinks3);
            assertEqual([size(s.sequence,1)-1 1;1 2], subs);
            assertEqual([1; 1], vals);
            
            s.allocateMemory(1);
            s.strandBreaks(1,1:2) = 1;
            [subs, vals] = find(s.strandBreaks3);
            assertEqual([1 1;2 2], subs);
            assertEqual([true; true], vals);
            
            s.allocateMemory(1);
            s.strandBreaks(1,1:2) = 1;
            [subs, vals] = find(s.strandBreaks5);
            assertEqual([2 1;1 2], subs);
            assertEqual([true; true], vals);
        end
        
        function testSingleStrandBreaks(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            s = this.state;
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            
            s.strandBreaks(:, :) = 0;
            s.strandBreaks(1:2, :) = 1;
            s.strandBreaks(3:5, 1) = 1;
            s.strandBreaks(5:7, 2) = 1;
            
            %example 1
            s.gapSites(:, :) = 0;
            
            singleStrandBreaks = CircularSparseMat([], false(0,1), size(s.strandBreaks), 1);
            singleStrandBreaks(3:4, 1) = true;
            singleStrandBreaks(6:7, 2) = true;
            assertEqual(singleStrandBreaks, s.singleStrandBreaks);
            
            %example 2
            s.gapSites(:, :) = 0;
            s.gapSites(2, 1) = 1;
            
            singleStrandBreaks = CircularSparseMat([], false(0,1), size(s.strandBreaks), 1);
            singleStrandBreaks(3:4,1) = true;
            singleStrandBreaks(6:7,2) = true;
            assertEqual(singleStrandBreaks, s.singleStrandBreaks);
            
            %example 3
            s.gapSites(:,:) = 0;
            s.gapSites(3,1) = 1;
            
            singleStrandBreaks = CircularSparseMat([], false(0,1), size(s.strandBreaks), 1);
            singleStrandBreaks(4, 1) = true;
            singleStrandBreaks(6:7, 2) = true;
            assertEqual(singleStrandBreaks, s.singleStrandBreaks);
            
            %example 4
            s.gapSites(:,:) = 0;
            s.gapSites(8,2) = 1;
            
            singleStrandBreaks = CircularSparseMat([], false(0,1), size(s.strandBreaks), 1);
            singleStrandBreaks(3:4, 1) = true;
            singleStrandBreaks(6, 2) = true;
            assertEqual(singleStrandBreaks, s.singleStrandBreaks);
            
            %example 5
            s.gapSites(:,:) = 0;
            s.gapSites(7,2) = 1;
            
            singleStrandBreaks = CircularSparseMat([], false(0,1), size(s.strandBreaks), 1);
            singleStrandBreaks(3:4, 1) = true;
            assertEqual(singleStrandBreaks, s.singleStrandBreaks);
        end
        
        function testDoubleStrandBreaks(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            s = this.state;
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = size(s.sequence, 1);
            
            s.strandBreaks(1:2, 1:2) = 1;
            s.strandBreaks(3:5, 1) = 1;
            s.strandBreaks(5:7, 2) = 1;
            
            doubleStrandBreaks = CircularSparseMat(zeros(0,2), false(0, 1), [size(s.sequence,1) s.nCompartments], 1);
            doubleStrandBreaks([1;2;5], 1:2) = true;
            
            assertEqual(doubleStrandBreaks, s.doubleStrandBreaks);
        end
        
        function testPloidy(this)
 
            s = this.state;
            
            s.allocateMemory(1);
            s.polymerizedRegions(1, 1:2) = s.sequenceLen;
            
            assertEqual(1, s.ploidy);
            
            s.polymerizedRegions(1, 1:3) = s.sequenceLen;
            
            assertEqual(1.5, s.ploidy);
            
            s.polymerizedRegions(1, 1:4) = s.sequenceLen;
            
            assertEqual(2, s.ploidy);
        end
    end
    
    %static methods
    methods
        function testShift35(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            s = this.state;
            
            assertEqual(CircularSparseMat([2 1;500 2],[1;2],[500 2],1), s.shiftCircularSparseMatBase3Prime(CircularSparseMat([1 1; 1 2],[1;2],[500 2],1)))
            assertEqual(CircularSparseMat([500 1;2 2],[1;2],[500 2],1), s.shiftCircularSparseMatBase5Prime(CircularSparseMat([1 1; 1 2],[1;2],[500 2],1)))
            assertEqual(CircularSparseMat([1 1;500 2],[1;2],[500 2],1), s.shiftCircularSparseMatBond3Prime(CircularSparseMat([1 1; 1 2],[1;2],[500 2],1)))
            assertEqual(CircularSparseMat([500 1;1 2],[1;2],[500 2],1), s.shiftCircularSparseMatBond5Prime(CircularSparseMat([1 1; 1 2],[1;2],[500 2],1)))
            assertEqual(CircularSparseMat([1 1;2 2],[1;2],[500 2],1), s.unshiftCircularSparseMatBond3Prime(CircularSparseMat([1 1; 1 2],[1;2],[500 2],1)))
            assertEqual(CircularSparseMat([2 1;1 2],[1;2],[500 2],1), s.unshiftCircularSparseMatBond5Prime(CircularSparseMat([1 1; 1 2],[1;2],[500 2],1)))
        end
    end
end