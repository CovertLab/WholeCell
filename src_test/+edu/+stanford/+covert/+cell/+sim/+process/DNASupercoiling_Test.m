% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef DNASupercoiling_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = DNASupercoiling_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
        
        function loadTestFixture(this)
            this.loadTestFixture@edu.stanford.covert.cell.sim.ProcessTestCase();
            
            m = this.process;
            c = m.chromosome;
            m.enzymeProperties = m.buildEnzymeProperties(); %rebuild function handles, required in MATLAB < 2010a
            
            c.initialize();
        end
        
        %test enzymeProperties built correctly
        function testEnzymeProperties(this)
            m = this.process;
            
            assertEqual(1:numel(m.enzymeProperties), [m.enzymeProperties.idx]);
        end
        
        function testSigmaEqualsZero(this)
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
            c.linkingNumbers(1, 1:2) = size(c.sequence, 1) / c.relaxedBasesPerTurn;
            m.enzymes(:) = 1;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            initialLK = double(c.linkingNumbers(1, 1));

            m.evolveState();
            nAtpUsed = 1e6 - m.substrates(m.substrateIndexs_atp);
            assertNumBound(1, m, m.enzymeIndexs_gyrase);
            assertNumBound(0, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            assertIn(nAtpUsed, [2 4]);
            assertEqual(...
                initialLK + m.gyraseDeltaLK * nAtpUsed / m.gyraseATPCost,...
                double(c.linkingNumbers(1, 1)));
        end

        function testSigmaLessThanGyraseLimit(this)
            m = this.process;
            c = m.chromosome;
            m.enzymes(:) = 1;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            relaxedLK = length(c.sequence) / c.relaxedBasesPerTurn;
            c.linkingNumbers(1,1:2) = relaxedLK + relaxedLK*(-0.2); %set sigma to -0.2
            initialLK = double(c.linkingNumbers(1, 1));
            
            m.evolveState();
            assertEqual(1e6, m.substrates(m.substrateIndexs_atp));
            assertNumBound(0, m, m.enzymeIndexs_gyrase);
            assertNumBound(0, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            assertEqual(...
                initialLK + m.topoIDeltaLK, double(c.linkingNumbers(1, 1)));
        end
        
        function testSigmaLessThanTopoIVLimit(this)
            m = this.process;
            c = m.chromosome;
            m.enzymes(:) = 5;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            relaxedLK = length(c.sequence) / c.relaxedBasesPerTurn;
            c.linkingNumbers(1,1:2) = relaxedLK + relaxedLK*(-0.05); %set sigma to -0.05
            initialLK = double(c.linkingNumbers(1, 1));
            
            m.evolveState();
            %1 or 2 topoI activities; 4 or 5 gyrase activities
            nAtpUsed = 1e6 - m.substrates(m.substrateIndexs_atp);
            assertIn(nAtpUsed, [8 10]); %only gyrase can use ATP
            assertNumBound(5, m, m.enzymeIndexs_gyrase);
            assertNumBound(0, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            assertTrue(...
                initialLK + m.topoIDeltaLK + m.gyraseDeltaLK * nAtpUsed / m.gyraseATPCost == ...
                double(c.linkingNumbers(1, 1)) || ...
                initialLK + m.topoIDeltaLK*2 + m.gyraseDeltaLK * nAtpUsed / m.gyraseATPCost == ...
                double(c.linkingNumbers(1, 1)));
        end
        
        function testSigmaGreaterThanTopoILimit(this)
            m = this.process;
            c = m.chromosome;
            m.enzymes(:) = 1;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            relaxedLK = length(c.sequence) / c.relaxedBasesPerTurn;
            c.linkingNumbers(1,1:2) = relaxedLK * 1.1; %set sigma to 0.1
            initialLK = double(c.linkingNumbers(1, 1));
            
            m.evolveState();
            nAtpUsed = 1e6 - m.substrates(m.substrateIndexs_atp);
            assertNumBound(1, m, m.enzymeIndexs_gyrase);
            assertNumBound(1, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            assertIn(nAtpUsed, [6 10]);  %gyrase and topoIV use ATP
            assertEqual(...
                initialLK + m.gyraseDeltaLK * nAtpUsed / m.gyraseATPCost,...
                double(c.linkingNumbers(1, 1))); %gyrase and topoIV delta LK and ATP usage are the same. 
        end
        
        function testNoWater(this)
            m = this.process;
            c = m.chromosome;
            m.enzymes(:) = 1;
            m.substrates(m.substrateIndexs_water) = 0;
            relaxedLK = length(c.sequence) / c.relaxedBasesPerTurn;
            c.linkingNumbers(1, 1:2) = relaxedLK + relaxedLK * (0.1); %set sigma to 0.1
            initialLK = double(c.linkingNumbers(1, 1));
            
            m.evolveState();
            assertEqual(0, m.substrates(m.substrateIndexs_water));
            assertNumBound(1, m, m.enzymeIndexs_gyrase);
            assertNumBound(1, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            assertEqual(...
                initialLK, double(c.linkingNumbers(1, 1))); %gyrase and topoIV delta LK and ATP usage are the same.
        end
        
        function testNoATP(this)
            m = this.process;
            c = m.chromosome;
            m.enzymes(:) = 1;
            m.substrates(m.substrateIndexs_atp) = 0;
            relaxedLK = length(c.sequence) / c.relaxedBasesPerTurn;
            c.linkingNumbers(1, 1:2) = relaxedLK + relaxedLK * (0.1); %set sigma to 0.1
            initialLK = double(c.linkingNumbers(1, 1));
            
            m.evolveState();
            assertEqual(0, m.substrates(m.substrateIndexs_atp));
            assertNumBound(1, m, m.enzymeIndexs_gyrase);
            assertNumBound(1, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            assertEqual(...
                initialLK, double(c.linkingNumbers(1, 1))); %gyrase and topoIV delta LK and ATP usage are the same.
        end
        
        function testLimitingATP1(this)
            m = this.process;
            c = m.chromosome;
            m.enzymes(:) = 5;
            m.substrates(m.substrateIndexs_atp) = 1;
            relaxedLK = length(c.sequence) / c.relaxedBasesPerTurn;
            c.linkingNumbers(1,1:2) = relaxedLK + relaxedLK*(-0.05); %set sigma to -0.05
            initialLK = double(c.linkingNumbers(1, 1));

            m.evolveState();
            assertNumBound(5, m, m.enzymeIndexs_gyrase);
            assertNumBound(0, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            assertEqual(1, m.substrates(m.substrateIndexs_atp));
            %1 or 2 TopoI activities
            assertTrue(...
                initialLK + m.topoIDeltaLK == double(c.linkingNumbers(1, 1)) ||...
                initialLK + m.topoIDeltaLK*2 == double(c.linkingNumbers(1, 1)));
        end
        
        function testLimitingATP2(this)
            m = this.process;
            c = m.chromosome;
            m.enzymes(:) = 1;
            m.substrates(m.substrateIndexs_atp) = 2;
            relaxedLK = length(c.sequence) / c.relaxedBasesPerTurn;
            c.linkingNumbers(1,1:2) = relaxedLK + relaxedLK*(0.1); %set sigma to 0.1
            initialLK = double(c.linkingNumbers(1, 1));
            
            m.evolveState();
            %both topoIV and gyrase should bind, only one should act. 
            assertEqual(0, m.substrates(m.substrateIndexs_atp));
            assertNumBound(1, m, m.enzymeIndexs_gyrase);
            assertNumBound(1, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            assertEqual(...
                initialLK + m.topoIVDeltaLK, double(c.linkingNumbers(1, 1)));
        end
        
        function testNoEnzyme1(this)
            m = this.process;
            c = m.chromosome;
            m.enzymes(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            relaxedLK = length(c.sequence) / c.relaxedBasesPerTurn;
            c.linkingNumbers(1,1:2) = relaxedLK + relaxedLK*(0.1); %set sigma to 0.1
            initialLK = double(c.linkingNumbers(1, 1));
            
            m.evolveState();
            assertEqual(1e6, m.substrates(m.substrateIndexs_atp));
            assertNumBound(0, m, m.enzymeIndexs_gyrase);
            assertNumBound(0, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            assertEqual(...
                initialLK, double(c.linkingNumbers(1, 1)));
        end
        
        function testNoEnzyme2(this)
            m = this.process;
            c = m.chromosome;
            m.enzymes(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            relaxedLK = length(c.sequence) / c.relaxedBasesPerTurn;
            c.linkingNumbers(1,1:2) = relaxedLK + relaxedLK*(-0.05); %set sigma to -0.05
            initialLK = double(c.linkingNumbers(1, 1));
            
            m.evolveState();
            assertEqual(1e6, m.substrates(m.substrateIndexs_atp));
            assertNumBound(0, m, m.enzymeIndexs_gyrase);
            assertNumBound(0, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            assertEqual(...
                initialLK, double(c.linkingNumbers(1, 1)));
        end

        function testTopoIVUnbinding(this)
            m = this.process;
            c = m.chromosome;
            m.enzymes(:) = 10;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            relaxedLK = length(c.sequence) / c.relaxedBasesPerTurn;
            c.linkingNumbers(1, 1:2) = relaxedLK + relaxedLK * 0.1; %set sigma to 0.1
            
            %bind topoIVs
            m.evolveState();
            assertIn(...
                nnz(c.complexBoundSites == m.enzymeToComplex(m.enzymeIndexs_topoIV)),...
                [1 inf]);
            
            %set sigma to be negative. 
            c.linkingNumbers(1, 1:2) = relaxedLK + relaxedLK * -0.2; %set sigma to -0.2
            m.evolveState();
            assertNumBound(0, m, m.enzymeIndexs_topoIV);
        end
        
        function testGyraseProcessivity(this)
            m = this.process;
            c = m.chromosome;
            relaxedLK = length(c.sequence) / c.relaxedBasesPerTurn;
            gyraseComplexIdx = m.enzymeToComplex(m.enzymeIndexs_gyrase);
            m.gyraseActivityRate = 1;
            m.gyraseMeanDwellTime = 3;
            
            nTrials = 30;
            totalSteps = 0;
            for i = 1:nTrials;
                m.enzymes(:) = 0;
                m.enzymes(m.enzymeIndexs_gyrase) = 1;
                m.boundEnzymes(:) = 0;
                m.substrates(:) = 0;
                m.substrates(m.substrateIndexs_atp) = 1e6;

                c.initialize();

                %bind gyrase
                m.evolveState();
                assertIn(...
                    nnz(c.complexBoundSites == gyraseComplexIdx),...
                    [1,inf]);
                
                %set linking number to prevent gyrase from binding again after
                %it falls off
                c.linkingNumbers(1, 1:2) = -2 * relaxedLK;

                %count how many steps it takes for all gyrase to fall off
                for j = 1:1000
                   m.evolveState();
                   totalSteps = totalSteps + 1;
                   if nnz(c.complexBoundSites == gyraseComplexIdx) == 0
                       break;
                   end
                end
                assertNumBound(0, m, m.enzymeIndexs_gyrase);
            end

            assertElementsAlmostEqual(...
                m.gyraseMeanDwellTime,...
                totalSteps / nTrials, 'relative', 0.2);
        end
        
        function testThreeRegionsThreeSigmas(this)
            import edu.stanford.covert.util.CircularSparseMat;
            m = this.process;
            c = m.chromosome;
            %assemble region such that 20000 bases have been polymerized by
            %each replication loop. Add 11 bases of a gap where helicase
            %sits. 
            pos = [1,1; 
                   20001,1; 
                   size(c.sequence,1)-19999+11,1; 
                   1,2; 
                   20001,2; 
                   size(c.sequence,1)-19999+11,2; 
                   1,3; 
                   size(c.sequence,1)-19999+11,3; 
                   1,4; 
                   size(c.sequence,1)-19999+11,4];
            lengths = [20000-11; 
                   size(c.sequence,1)-2*20000;
                   20000-11;
                   20000-11; 
                   size(c.sequence,1)-2*20000;
                   20000-11;
                   20000-11;
                   20000-11;
                   20000-11;
                   20000-11];
            sigmas = [-0.05;
                0.2;
                -0.05;
                -0.05;
                0.2;
                -0.05;
                -0.2;
                -0.2;
                -0.2;
                -0.2];
            relaxedLK = lengths / c.relaxedBasesPerTurn;
            LKs = relaxedLK + relaxedLK.*sigmas;
            c.polymerizedRegions =...
                CircularSparseMat(pos, lengths, [size(c.sequence,1),4], 1);
            c.linkingNumbers =...
                CircularSparseMat(pos, LKs, [size(c.sequence,1),4], 1);
            m.enzymes(:) = 1;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            initialLK = c.linkingNumbers;

            m.evolveState();
            nAtpUsed = 1e6 - m.substrates(m.substrateIndexs_atp);  
            %gyrase and topoIV should have both bound
            assertTrue(ismember(nAtpUsed, [6 8 10]));
            assertNumBound(1, m, m.enzymeIndexs_gyrase);
            assertNumBound(1, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            %regions 2 and 5 may have a lower linking number
            assertTrue(initialLK(pos(2,:))>=c.linkingNumbers(pos(2,:)));
            assertTrue(initialLK(pos(5,:))>=c.linkingNumbers(pos(5,:)));
            %regions 7-10 may have a higher linking number
            assertTrue(initialLK(pos(7,:))<=c.linkingNumbers(pos(7,:)));
            assertTrue(initialLK(pos(8,:))<=c.linkingNumbers(pos(8,:)));
            assertTrue(initialLK(pos(9,:))<=c.linkingNumbers(pos(9,:)));
            assertTrue(initialLK(pos(10,:))<=c.linkingNumbers(pos(10,:)));
            %based paired regions should have the same linking number
            assertEqual((c.linkingNumbers(pos(1,:))), (c.linkingNumbers(pos(4,:))));
            assertEqual((c.linkingNumbers(pos(2,:))), (c.linkingNumbers(pos(5,:))));
            assertEqual((c.linkingNumbers(pos(3,:))), (c.linkingNumbers(pos(6,:))));
            assertEqual((c.linkingNumbers(pos(7,:))), (c.linkingNumbers(pos(9,:))));
            assertEqual((c.linkingNumbers(pos(8,:))), (c.linkingNumbers(pos(10,:))));
        end
        
        function testThreeRegionsTwoSigmas(this)
            import edu.stanford.covert.util.CircularSparseMat;
            m = this.process;
            c = m.chromosome;
            %assemble region such that 20000 bases have been polymerized by
            %each replication loop. Add 11 bases of a gap where helicase
            %sits. 
            pos = [1,1; 
                   20001,1; 
                   size(c.sequence,1)-19999+11,1; 
                   1,2; 
                   20001,2; 
                   size(c.sequence,1)-19999+11,2; 
                   1,3; 
                   size(c.sequence,1)-19999+11,3; 
                   1,4; 
                   size(c.sequence,1)-19999+11,4];
            lengths = [20000-11; 
                   size(c.sequence,1)-2*20000;
                   20000-11;
                   20000-11; 
                   size(c.sequence,1)-2*20000;
                   20000-11;
                   20000-11;
                   20000-11;
                   20000-11;
                   20000-11];
            sigmas = [-0.2;
                0.2;
                -0.2;
                -0.2;
                0.2;
                -0.2;
                -0.2;
                -0.2;
                -0.2;
                -0.2];
            LKs = lengths/c.relaxedBasesPerTurn + (lengths/c.relaxedBasesPerTurn).*(sigmas);
            c.polymerizedRegions =...
                CircularSparseMat(pos, lengths, [size(c.sequence,1),4], 1);
            c.linkingNumbers =...
                CircularSparseMat(pos, LKs, [size(c.sequence,1),4], 1);
            m.enzymes(:) = 1;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            initialLK = c.linkingNumbers;

            m.evolveState();
            nAtpUsed = 1e6 - m.substrates(m.substrateIndexs_atp);  
            %gyrase and topoIV should have both bound
            assertTrue(ismember(nAtpUsed, [6 8 105]));
            assertNumBound(1, m, m.enzymeIndexs_gyrase);
            assertNumBound(1, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            %regions 2 and 5 will have a lower linking number
            assertEqual(...
                (initialLK(pos(2,:))) + (m.gyraseDeltaLK * nAtpUsed) / m.gyraseATPCost,...
                (c.linkingNumbers(pos(2,:))));
            assertEqual(...
                (initialLK(pos(5,:))) + (m.gyraseDeltaLK * nAtpUsed) / m.gyraseATPCost,...
                (c.linkingNumbers(pos(5,:))));
            %regions 1,3-4,6-10 may have a higher linking number where the
            %deltaLK is m.topoIDeltaLK
            for i = [1,3,4,6,7,8,9,10]
                assertTrue((initialLK(pos(i,:)))<=(c.linkingNumbers(pos(i,:)))...
                    && (initialLK(pos(i,:)))+m.topoIDeltaLK>=(c.linkingNumbers(pos(i,:))));
            end
            %based paired regions should have the same linking number
            assertEqual((c.linkingNumbers(pos(1,:))), (c.linkingNumbers(pos(4,:))));
            assertEqual((c.linkingNumbers(pos(2,:))), (c.linkingNumbers(pos(5,:))));
            assertEqual((c.linkingNumbers(pos(3,:))), (c.linkingNumbers(pos(6,:))));
            assertEqual((c.linkingNumbers(pos(7,:))), (c.linkingNumbers(pos(9,:))));
            assertEqual((c.linkingNumbers(pos(8,:))), (c.linkingNumbers(pos(10,:))));
        end
        
        function testTwoRegionsTwoSigmas(this)
            import edu.stanford.covert.util.CircularSparseMat;
            m = this.process;
            c = m.chromosome;
            %assemble 2 chromosomes
            pos = [1,1; 
                   1,2; 
                   1,3; 
                   1,4];
            lengths = [size(c.sequence,1);
                   size(c.sequence,1);
                   size(c.sequence,1);
                   size(c.sequence,1)];
            sigmas = [-0.2;
                -0.2;
                0.2;
                0.2];

            LKs = lengths/c.relaxedBasesPerTurn + (lengths/c.relaxedBasesPerTurn).*(sigmas);
            c.polymerizedRegions =...
                CircularSparseMat(pos, lengths, [size(c.sequence,1),4], 1);
            c.linkingNumbers =...
                CircularSparseMat(pos, LKs, [size(c.sequence,1),4], 1);
            m.enzymes(:) = 1;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            initialLK = c.linkingNumbers;

            m.evolveState();
            nAtpUsed = 1e6 - m.substrates(m.substrateIndexs_atp);  
            assertTrue(initialLK(pos(1,:)) < c.linkingNumbers(pos(1,:)));
            assertTrue(initialLK(pos(2,:)) < c.linkingNumbers(pos(2,:)));
            assertTrue(initialLK(pos(3,:)) > c.linkingNumbers(pos(3,:)));
            assertTrue(initialLK(pos(4,:)) > c.linkingNumbers(pos(4,:)));
            %gyrase and topoIV should have both bound
            assertTrue(ismember(nAtpUsed, [6 8 10]));
            assertNumBound(1, m, m.enzymeIndexs_gyrase);
            assertNumBound(1, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            %based paired regions should have the same linking number
            assertEqual((c.linkingNumbers(pos(1,:))), (c.linkingNumbers(pos(2,:))));
            assertEqual((c.linkingNumbers(pos(3,:))), (c.linkingNumbers(pos(4,:))));
        end
        
        function testTwoRegionsOneSigma(this)
            import edu.stanford.covert.util.CircularSparseMat;
            m = this.process;
            c = m.chromosome;
            %assemble 2 chromosomes
            pos = [1,1; 
                   1,2; 
                   1,3; 
                   1,4];
            lengths = [size(c.sequence,1);
                   size(c.sequence,1);
                   size(c.sequence,1);
                   size(c.sequence,1)];
            sigmas = [-0.2;
                -0.2;
                -0.2;
                -0.2];

            LKs = lengths/c.relaxedBasesPerTurn + (lengths/c.relaxedBasesPerTurn).*(sigmas);
            c.polymerizedRegions =...
                CircularSparseMat(pos, lengths, [size(c.sequence,1),4], 1);
            c.linkingNumbers =...
                CircularSparseMat(pos, LKs, [size(c.sequence,1),4], 1);
            m.enzymes(:) = 1;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            initialLK = c.linkingNumbers;

            m.evolveState();
            nAtpUsed = 1e6 - m.substrates(m.substrateIndexs_atp);  
            assertEqual(1e6, m.substrates(m.substrateIndexs_atp));  
            %gyrase and topoIV should not have bound
            assertEqual(nAtpUsed, 0);
            assertNumBound(0, m, m.enzymeIndexs_gyrase);
            assertNumBound(0, m, m.enzymeIndexs_topoIV);
            assertNumBound(0, m, m.enzymeIndexs_topoI);
            %assert linking numbers have increased or stayed the same
            assertTrue(initialLK(pos(1,:)) <= c.linkingNumbers(pos(1,:)));
            assertTrue(initialLK(pos(2,:)) <= c.linkingNumbers(pos(2,:)));
            assertTrue(initialLK(pos(3,:)) <= c.linkingNumbers(pos(3,:)));
            assertTrue(initialLK(pos(4,:)) <= c.linkingNumbers(pos(4,:)));
        end
        
        function testSingleStrandedRegions(this)
            import edu.stanford.covert.util.CircularSparseMat;
            m = this.process;
            c = m.chromosome;
            %assemble 2 chromosomes with sigma = 0
            pos = [1,1; 
                   1,2; 
                   1,3; 
                   1,4];
            lengths = [size(c.sequence,1);
                   size(c.sequence,1);
                   size(c.sequence,1);
                   size(c.sequence,1)];
            sigmas = [0;
                0;
                0;
                0];

            LKs = lengths/c.relaxedBasesPerTurn + (lengths/c.relaxedBasesPerTurn).*(sigmas);
            c.polymerizedRegions =...
                CircularSparseMat(pos, lengths, [size(c.sequence,1),4], 1);
            c.linkingNumbers =...
                CircularSparseMat(pos, LKs, [size(c.sequence,1),4], 1);
            m.enzymes(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;

            m.evolveState();
            
            preBreakLK = c.linkingNumbers;
            
            %add a break of length 20000
            pos = [1,1; 
                   1,2; 
                   1,3; 
                   1,4
                   560077, 2];
            lengths = [size(c.sequence,1);
                   540076;
                   size(c.sequence,1);
                   size(c.sequence,1);
                   20000];

            c.polymerizedRegions =...
                CircularSparseMat(pos, lengths, [size(c.sequence,1),4], 1);
            m.enzymes(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;

            c.calcLinkingNumbers_minFreeEnergy();
            
            assertEqual(preBreakLK(pos(3,:)),c.linkingNumbers_minFreeEnergy(pos(3,:)));
            assertEqual(preBreakLK(pos(4,:)),c.linkingNumbers_minFreeEnergy(pos(4,:)));
            turnsWithBreak = (size(c.sequence,1)-20000)/c.relaxedBasesPerTurn;
            assertElementsAlmostEqual(turnsWithBreak, c.linkingNumbers_minFreeEnergy(pos(2,:))+...
                c.linkingNumbers_minFreeEnergy(pos(5,:)), ...
                'relative', 0.000000001);
        end
        
        function testOneRegionFoldChange(this)
            m = this.process;
            c = m.chromosome;
            r = m.rnaPolymerase;
            
            c.initialize();
            c.linkingNumbers(1, 1:2) = size(c.sequence, 1) / c.relaxedBasesPerTurn;
            m.enzymes(:) = 0;

            m.evolveState();
            %no fold change in non-existant chromosome 2
            assertEqual(ones(m.numTranscriptionUnits,1), r.supercoilingBindingProbFoldChange(:,2));
            %assert correct fold changes in chromosome 1
            assertEqual(m.foldChangeIntercepts(1),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(1),1));
            assertEqual(m.foldChangeIntercepts(2),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(2),1));
            assertEqual(m.foldChangeIntercepts(3),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(3),1));
        end
               
        function testFoldChange3Sigmas(this)
            import edu.stanford.covert.util.CircularSparseMat;
            m = this.process;
            c = m.chromosome;
            r = m.rnaPolymerase;
            
            %assemble region such that 20000 bases have been polymerized by
            %each replication loop. Add 11 bases of a gap where helicase
            %sits. 
            pos = [1,1; 
                   150001,1; 
                   size(c.sequence,1)-149999+11,1; 
                   1,2; 
                   150001,2; 
                   size(c.sequence,1)-149999+11,2; 
                   1,3; 
                   size(c.sequence,1)-149999+11,3; 
                   1,4; 
                   size(c.sequence,1)-149999+11,4];
            lengths = [150000-11; 
                   size(c.sequence,1)-2*150000;
                   150000-11;
                   150000-11; 
                   size(c.sequence,1)-2*150000;
                   150000-11;
                   150000-11;
                   150000-11;
                   150000-11;
                   150000-11];
            sigmas = [-0.05;
                0.02;
                -0.05;
                -0.05;
                0.02;
                -0.05;
                0.08;
                0.08;
                0.08;
                0.08];
            relaxedLK = lengths / c.relaxedBasesPerTurn;
            LKs = relaxedLK + relaxedLK.*sigmas;
            c.polymerizedRegions =...
                CircularSparseMat(pos, lengths, [size(c.sequence,1),4], 1);
            c.linkingNumbers =...
                CircularSparseMat(pos, LKs, [size(c.sequence,1),4], 1);
            m.enzymes(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;

            m.evolveState();
            assertElementsAlmostEqual(m.foldChangeSlopes(1)*-0.05+m.foldChangeIntercepts(1),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(1),1), 'relative', 0.000000001);
            assertElementsAlmostEqual(m.foldChangeSlopes(2)*0.02+m.foldChangeIntercepts(2),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(2),1), 'relative', 0.000000001);
            assertElementsAlmostEqual(m.foldChangeSlopes(3)*-0.05+m.foldChangeIntercepts(3),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(3),1), 'relative', 0.000000001);
            assertElementsAlmostEqual(m.foldChangeSlopes(1)*m.foldChangeUpperSigmaLimit+m.foldChangeIntercepts(1),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(1),2), 'relative', 0.000000001);
            assertEqual(1,...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(2),2));
            assertElementsAlmostEqual(m.foldChangeSlopes(3)*m.foldChangeUpperSigmaLimit+m.foldChangeIntercepts(3),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(3),2), 'relative', 0.000000001);            
        end
        
        function testFoldChange3SigmasOutOfLimits(this)
            import edu.stanford.covert.util.CircularSparseMat;
            m = this.process;
            c = m.chromosome;
            r = m.rnaPolymerase;
            
            %assemble region such that 20000 bases have been polymerized by
            %each replication loop. Add 11 bases of a gap where helicase
            %sits. 
            pos = [1,1; 
                   150001,1; 
                   size(c.sequence,1)-149999+11,1; 
                   1,2; 
                   150001,2; 
                   size(c.sequence,1)-149999+11,2; 
                   1,3; 
                   size(c.sequence,1)-149999+11,3; 
                   1,4; 
                   size(c.sequence,1)-149999+11,4];
            lengths = [150000-11; 
                   size(c.sequence,1)-2*150000;
                   150000-11;
                   150000-11; 
                   size(c.sequence,1)-2*150000;
                   150000-11;
                   150000-11;
                   150000-11;
                   150000-11;
                   150000-11];
            sigmas = [-0.1;
                -0.08;
                -0.1;
                -0.1;
                -0.08;
                -0.1;
                0.1;
                0.1;
                0.1;
                0.1];
            relaxedLK = lengths / c.relaxedBasesPerTurn;
            LKs = relaxedLK + relaxedLK.*sigmas;
            c.polymerizedRegions =...
                CircularSparseMat(pos, lengths, [size(c.sequence,1),4], 1);
            c.linkingNumbers =...
                CircularSparseMat(pos, LKs, [size(c.sequence,1),4], 1);
            m.enzymes(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;

            m.evolveState();
            assertElementsAlmostEqual(m.foldChangeSlopes(1)*m.foldChangeLowerSigmaLimit+m.foldChangeIntercepts(1),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(1),1), 'relative', 0.000000001);
            assertElementsAlmostEqual(m.foldChangeSlopes(2)*m.foldChangeLowerSigmaLimit+m.foldChangeIntercepts(2),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(2),1), 'relative', 0.000000001);
            assertElementsAlmostEqual(m.foldChangeSlopes(3)*m.foldChangeLowerSigmaLimit+m.foldChangeIntercepts(3),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(3),1), 'relative', 0.000000001);
            assertElementsAlmostEqual(m.foldChangeSlopes(1)*m.foldChangeUpperSigmaLimit+m.foldChangeIntercepts(1),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(1),2), 'relative', 0.000000001);
            assertEqual(1,...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(2),2));
            assertElementsAlmostEqual(m.foldChangeSlopes(3)*m.foldChangeUpperSigmaLimit+m.foldChangeIntercepts(3),...
                r.supercoilingBindingProbFoldChange(m.tuIndexs(3),2), 'relative', 0.000000001);            
        end
        
        function testGeneEssentiality_PositivelySupercoiled(this)
            m = this.process;
            c = m.chromosome;
            r = m.rnaPolymerase;
            
            c.allocateMemory(1);
            c.polymerizedRegions(1, 1:2) = c.sequenceLen;
            c.linkingNumbers(1, 1:2) = (1 + (m.topoISigmaLimit - 1e-2)) * c.sequenceLen / c.relaxedBasesPerTurn;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 10;
            m.topoIVActivityRate = 1e3;
            m.gyraseActivityRate = 1e6;
            m.buildEnzymeProperties();
            m.boundEnzymes(:) = 0;
            r.supercoilingBindingProbFoldChange = m.calcRNAPolymeraseBindingProbFoldChange();
            
            this.helpTestGeneEssentiality({
                'MG_003';
                'MG_004';
                }, @isFunctioningProperly);
            
            function val = isFunctioningProperly(~, i)
                val = logical(c.linkingNumbers(1, 1) < i.State_Chromosome.linkingNumbers(1, 1));
            end
        end
        
        function testGeneEssentiality_NegativelySupercoiled(this)
            m = this.process;
            c = m.chromosome;
            r = m.rnaPolymerase;
            
            c.allocateMemory(1);
            c.polymerizedRegions(1, 1:2) = c.sequenceLen;
            c.linkingNumbers(1, 1:2) = (1 + (m.gyraseSigmaLimit - 1e-2)) * c.sequenceLen / c.relaxedBasesPerTurn;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 10;
            m.topoIVActivityRate = 1e3;
            m.gyraseActivityRate = 1e6;
            m.buildEnzymeProperties();
            m.boundEnzymes(:) = 0;
            r.supercoilingBindingProbFoldChange = m.calcRNAPolymeraseBindingProbFoldChange();
            
            this.helpTestGeneEssentiality({
                'MG_122';
                }, @isFunctioningProperly);
            
            function val = isFunctioningProperly(~, i)
                val = logical(c.linkingNumbers(1, 1) > i.State_Chromosome.linkingNumbers(1, 1));
            end
        end
    end
end

function assertNumBound(n, m, enzymeIndex)
    c = m.chromosome;
    mIdx = m.enzymeToMonomer(enzymeIndex);
    if ~isempty(mIdx)
        assertEqual(n, nnz(c.monomerBoundSites == mIdx));
    else
        cIdx = m.enzymeToComplex(enzymeIndex);
        assertEqual(n, nnz(c.complexBoundSites == cIdx));
    end
end
