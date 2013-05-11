%Transcriptional regulation process test case
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/19/2010
classdef TranscriptionalRegulation_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = TranscriptionalRegulation_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end
    
    %tests
    methods
        function testAssumptions(this)
            m = this.process;
            
            %assert that enzymeMonomerLocalIndexs, enzymeComplexLocalIndexs are
            %sorted so isDnaBound in ChromosomeProcessAspect can be used safely
            assertTrue(issorted(m.enzymeMonomerLocalIndexs));
            assertTrue(issorted(m.enzymeComplexLocalIndexs));
        end
        
        %Tests if binding sites overlap. Overlap doesn't pose any problem to the
        %model, but the model could be made more computationally efficient for
        %the case of overlapping binding sites.
        function testNoTFBindingSiteOverlap(this)
            m = this.process;
            c = m.chromosome;
            footprints = zeros(size(m.enzymes));
            footprints(m.enzymeMonomerLocalIndexs) = c.monomerDNAFootprints(m.enzymeMonomerGlobalIndexs);
            footprints(m.enzymeComplexLocalIndexs) = c.complexDNAFootprints(m.enzymeComplexGlobalIndexs);
            
            starts = m.tfPositionStrands(1:end/2, 1);
            ends = starts + footprints(m.tfIndexs(:,1)) - 1;
            idxs =  (1:numel(starts))';
            
            tmp = find(ends > c.sequenceLen);
            starts = [starts; max(1, starts(tmp) - c.sequenceLen)];
            ends = [ends; ends(tmp) - c.sequenceLen];
            idxs = [idxs; idxs(tmp)];
            
            tmp = find(starts < 0);
            starts = [starts; starts(tmp) + c.sequenceLen];
            ends = [ends; min(c.sequenceLen, ends(tmp) + c.sequenceLen)];
            idxs = [idxs; idxs(tmp)];
            
            overlaps = zeros(size(m.tfPositionStrands, 1));
            for i = 1:size(m.tfPositionStrands, 1) / 2-1
                overlap = find(...
                    (starts(i) >= starts(i+1:end) & starts(i) <= ends(i+1:end)) | ...
                    (ends(i)   >= starts(i+1:end) & ends(i)   <= ends(i+1:end)) | ...
                    (starts(i) <= starts(i+1:end) & ends(i)   >= ends(i+1:end)));
                overlaps(i, idxs(overlap)) = 1;
                overlaps(idxs(overlap), i) = 1;
            end
            
            assertAllEqual(0, overlaps);
        end
        
        function testBindingProbabilityFoldChange(this)
            this.releaseAllTranscriptionFactors();
            m = this.process;
            r = m.rnaPolymerase;
            
            m.tfIndexs = [1 1 1 1 3 3;
                1 1 1 1 3 3]';
            m.tuIndexs = [1 2 3 4 1 4;
                1 2 3 4 1 4]';
            m.tfActivities = [0.5 0.6 0.8 2.0 3.4 4.4;
                0.5 0.6 0.8 2.0 3.4 4.4]';
            m.tfAffinities = ones(6, 2);
            m.tfPositionStrands = [
                (1000:1000:6000)' ones(6,1);
                (1000:1000:6000)' ones(6,1)*3];
            m.otherActivities = ones([length(m.transcriptionUnitWholeCellModelIDs) 5]);
            m.otherActivities(1, 2) = 3;
            m.otherActivities(2, 4) = 1.5;
            m.enzymes = [4 1 3 2 0]';
            m.boundEnzymes = [0 0 0 0 0]';
            
            m.evolveState();
            assertEqual([0 1 1 2 0], m.enzymes');
            assertEqual([4 0 2 0 0], m.boundEnzymes');
            assertEqual(logical([1 1 1 1 1 1; 0 0 0 0 0 0]), m.tfBoundPromoters');
            
            fc = r.transcriptionFactorBindingProbFoldChange;
            assertElementsAlmostEqual([5.1 0.9 0.8 8.8; 3.0 1.5 1 1], fc(1:4,:)');
            assertAllEqual(1, fc(5:end,:));
        end
        
        function testTranscriptionFactorBinding(this)
            m = this.process;
            c = m.chromosome;
            r = m.rnaPolymerase;
            c.initialize();
            
            A = [
                1  15  0.200   5.000;
                1  31  0.780   1.282;
                1  22  0.537   1.862;
                1  40  0.753   1.328;
                1  97  0.500   2.000;
                1  39  0.884   1.131;
                1  28  1.140   1.140;
                1  67  1.248   1.248;
                3  17  0.343   2.916;
                3  84  0.478   2.092;
                3  33  0.175   5.714;
                3  10  0.382   2.618;
                3  77  0.078  12.821;
                3  17  0.901   1.110;
                4  91  1.000   1.000;
                4  67  1.000   1.000;
                5  15  6.253   6.253;
                5  27 10.012  10.012;
                5  45  3.598   3.598;
                5  70  0.322   3.106;
                5  71  0.309   3.236;
                5  90  9.257   9.257;
                5   8  0.203   4.926;
                5  11 14.847  14.847;
                5  70  5.437   5.437;
                5  71  5.715   5.715;
                5  89  3.377   3.377;
                5  99  0.253   3.953;
                5  39  0.297   3.337;
                5  47  3.812   3.812];
            m.tfIndexs = repmat(A(:,1), [1 2]);
            m.tuIndexs = repmat(A(:,2), [1 2]);
            m.tfActivities = repmat(A(:,3), [1 2]);
            m.tfAffinities = repmat(A(:,4), [1 2]);
            m.tfPositionStrands = [
                (1:length(A))' * 1000  ones(length(A),1);
                (1:length(A))' * 1000  ones(length(A),1)*3];
            m.otherActivities = ones([length(m.transcriptionUnitWholeCellModelIDs) 5]);
            
            this.releaseAllTranscriptionFactors();
            
            m.enzymes(:) = 10; %saturating for some transcription factors, not for others
            m.boundEnzymes(:) = 0;
            assertAllEqual(1, r.transcriptionFactorBindingProbFoldChange);
            assertAllEqual(false, m.tfBoundPromoters);
            
            m.evolveState();
            
            %assert enzymes bound, and to which promoters
            assertEqual([2 10 4 8 0], m.enzymes');
            assertEqual([8 0 6 2 10], m.boundEnzymes');
            for i = 1:numel(m.enzymes)
                assertEqual(m.boundEnzymes(i), sum(m.tfBoundPromoters(m.tfIndexs == i)));
            end
            
            %assert fold changes have changed as expected (based on activities)
            fc = r.transcriptionFactorBindingProbFoldChange;
            assertAllEqual(1, fc(:,2));
            assertElementsAlmostEqual(...
                m.calcBindingProbabilityFoldChange(m.tfBoundPromoters),...
                fc, 'absolute', 1e-8);
            
            unregulatedTUs = true(size(m.transcriptionUnitWholeCellModelIDs));
            unregulatedTUs(m.tuIndexs(:,1)) = false;
            assertAllEqual(1, fc(unregulatedTUs, :));
        end
        
        function testTranscriptionFactorBindingUnsaturated(this)
            m = this.process;
            c = m.chromosome;
            c.initialize();
            
            this.releaseAllTranscriptionFactors();
            
            m.enzymes(:) = 100;  %ample
            m.boundEnzymes(:) = 0;
            
            m.evolveState();
            assertAllEqual(true, m.tfBoundPromoters(:,1)'); %all promoters bound
            assertEqual(m.boundTFs, m.boundEnzymes);
            assertEqual(100 - m.boundEnzymes, m.enzymes);
        end
        
        function testTranscriptionFactorBindingSaturated(this)
            m = this.process;
            
            this.bindAllTranscriptionFactors();
            assertAllEqual(true, m.tfBoundPromoters(:,1)'); %all still bound
            
            m.enzymes(:) = 10; %excess
            m.evolveState();
            
            assertAllEqual(true, m.tfBoundPromoters(:,1)'); %all still bound
            assertEqual(m.boundTFs, m.boundEnzymes);
            assertAllEqual(10, m.enzymes);
        end
        
        function testAllTranscriptionUnitPromotersStayBound(this)
            this.bindAllTranscriptionFactors();
            m = this.process;
            m.enzymes(:) = 0;  %no excess
            m.boundEnzymes = m.boundTFs;
            
            m.evolveState();
            assertAllEqual(true, m.tfBoundPromoters(:,1)'); %all still bound
            assertEqual(m.boundTFs, m.boundEnzymes);
            assertAllEqual(0, m.enzymes);
        end
        
        function testFairnessOfTranscriptionFactorBinding(this)
            this.releaseAllTranscriptionFactors();
            m = this.process;
            m.tfIndexs = [1 1 1 1 1 1; 1 1 1 1 1 1]';
            m.tuIndexs = m.tuIndexs(1:6,:);
            m.tfAffinities = [1 1 1 1 1 1; 1 1 1 1 0 0]';  %ten binding sites
            m.tfActivities = zeros(6,2);
            m.tfPositionStrands = [
                1000 1; 2000 1; 3000 1; 4000 1; 5000 1; 6000 1;
                1000 3; 2000 3; 3000 3; 4000 3; 5000 3; 6000 3];
            m.chromosome.polymerizedRegions(1, :) = m.chromosome.sequenceLen;
            m.chromosome.linkingNumbers(1, :) = m.chromosome.linkingNumbers([1 1]);
            
            totalBound = zeros(6,2);
            for i = 1:50
                this.releaseAllTranscriptionFactors();
                m.enzymes(:) = 0;
                m.enzymes(1) = 5; %only one kind of transcription factor
                m.boundEnzymes(:) = 0;
                
                m.evolveState();
                b = m.tfBoundPromoters;
                assertEqual(5, sum(sum(b(1:10))));  %all bound
                assertEqual([false false], b(11:end));
                assertEqual([false true ], unique(b)');
                totalBound = totalBound + b;
            end
            
            assertEqual(25, mean(totalBound(1:10)));
            assertIn(std(totalBound(1:10)), [0 4]);
        end
        
        function testChromosomeFairness(this)
             m = this.process;
             c = m.chromosome;
             rnaPol = m.rnaPolymerase;
             
             fc = zeros(size(rnaPol.transcriptionFactorBindingProbFoldChange));
             tfBoundPromoters = zeros(size(m.tfBoundPromoters));
             for i = 1:1000
                 m.enzymes(:) = 1;
                 m.boundEnzymes(:) = 0;
                 
                 c.initialize();
                 c.polymerizedRegions(1, :) = c.polymerizedRegions([1 1]);
                 c.linkingNumbers(1, :) = c.linkingNumbers([1 1]);
                 
                 m.evolveState();
                 tfBoundPromoters = tfBoundPromoters + m.tfBoundPromoters;
                 fc = fc + rnaPol.transcriptionFactorBindingProbFoldChange;
             end
             
             assertElementsAlmostEqual(tfBoundPromoters(:, 1)', tfBoundPromoters(:, 2)', 'relative', 0.20, 10);
             assertElementsAlmostEqual(fc(:, 1), fc(:, 2), 'relative', 0.10);
        end
        
        function testGeneEssentiality(this)
            this.releaseAllTranscriptionFactors();
            m = this.process;
            c = m.chromosome;
            c.initialize();
            m.enzymes(:) = 1;
            this.helpTestGeneEssentiality({
                'MG_101';     %gntR, uncharacterized HTH-type transcriptional regulator
                'MG_127';     %spx, Spx subfamily protein
                'MG_205';     %heat-inducible transcription repressor HrcA, putative
                'MG_236';     %fur, ferric uptake repressor
                'MG_428'},... %gerE / luxR, LuxR bacterial regulatory protein, putative
                @(m,~) all(m.boundTFs == [1 0 1 1 1]') && m.enzymes(m.enzymeIndexs_spx));
        end
    end
    
    methods (Access = private)
        function bindAllTranscriptionFactors(this)
            m = this.process;
            c = m.chromosome;
            r = m.rnaPolymerase;
            
            c.monomerBoundSites(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            
            m.enzymes(:) = 100; %ample
            m.boundEnzymes(:) = 0;
            for i = 1:length(m.enzymes)
                m.bindProteinToChromosome(...
                    m.tfPositionStrands(m.tfIndexs(:) == i, :), i, [], [], [], false);
            end
            r.transcriptionFactorBindingProbFoldChange = m.calcBindingProbabilityFoldChange(m.tfBoundPromoters);
        end
        
        function releaseAllTranscriptionFactors(this)
            m = this.process;
            c = m.chromosome;
            r = m.rnaPolymerase;
            
            c.monomerBoundSites(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            c.gapSites(:, :) = 0;
            c.abasicSites(:, :) = 0;
            c.damagedSugarPhosphates(:, :) = 0;
            c.damagedBases(:, :) = 0;
            c.strandBreaks(:, :) = 0;
            c.intrastrandCrossLinks(:, :) = 0;
            c.hollidayJunctions(:, :) = 0;
            
            for i = 1:length(m.enzymes)
                m.releaseProteinFromChromosome(i, inf, [], []);
            end
            r.transcriptionFactorBindingProbFoldChange = m.calcBindingProbabilityFoldChange(m.tfBoundPromoters);
        end
    end
end
