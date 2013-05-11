%Protein processing II process test case
%
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/9/2010
classdef ProteinProcessingII_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = ProteinProcessingII_Test(methodName)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(methodName);
        end
        
        function testOneMonomerRequiringNoProcessing(this)
            m = this.process;
            m.lipoproteinMonomerIndexs = [];
            m.secretedMonomerIndexs = [];
            m.unprocessedMonomerIndexs = 1;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.unprocessedMonomers = 1;
            m.processedMonomers = 0;
            m.signalSequenceMonomers = 0;

            m.evolveState();
            assertEqual(0, m.unprocessedMonomers);
            assertEqual(1, m.processedMonomers);
            assertEqual(0, m.signalSequenceMonomers);
            assertEqual(zeros(size(m.substrates)), m.substrates);
        end

        function testOneSecretedMonomer(this)
            m = this.process;
            m.lipoproteinMonomerIndexs = [];
            m.secretedMonomerIndexs = 1;
            m.unprocessedMonomerIndexs = [];
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalPeptidase) = 1;
            m.unprocessedMonomers = 1;
            m.processedMonomers = 0;
            m.signalSequenceMonomers = 0;

            m.evolveState();
            assertEqual(0, m.unprocessedMonomers);
            assertEqual(1, m.processedMonomers);
            assertEqual(1, m.signalSequenceMonomers);
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(1, m.enzymes(m.enzymeIndexs_signalPeptidase));
        end

        function testOneLipoprotein(this)
            m = this.process;
            m.lipoproteinMonomerIndexs = 1;
            m.secretedMonomerIndexs = [];
            m.unprocessedMonomerIndexs = [];
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1;
            m.substrates(m.substrateIndexs_PG160) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalPeptidase) = 1;
            m.enzymes(m.enzymeIndexs_diacylglycerylTransferase) = 99;
            m.unprocessedMonomers = 1;
            m.processedMonomers = 0;
            m.signalSequenceMonomers = 0;

            m.evolveState();
            assertEqual(0, m.unprocessedMonomers);
            assertEqual(1, m.processedMonomers);
            assertEqual(1, m.signalSequenceMonomers);
            assertEqual(0, m.substrates(m.substrateIndexs_water));
            assertEqual(0, m.substrates(m.substrateIndexs_PG160));
            assertEqual(1, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1, m.substrates(m.substrateIndexs_SNGLYP));
            assertEqual(1, m.enzymes(m.enzymeIndexs_signalPeptidase));
            assertEqual(99, m.enzymes(m.enzymeIndexs_diacylglycerylTransferase));
        end
        
        function testNoProcessingWithoutWater(this)
            m = this.process;
            m.lipoproteinMonomerIndexs = 1;
            m.secretedMonomerIndexs = 2;
            m.unprocessedMonomerIndexs = [];
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_PG160) = 1e3;
            m.enzymes(:) = 1e3;
            m.unprocessedMonomers = [1; 1];
            m.processedMonomers = [0; 0];
            m.signalSequenceMonomers = [0; 0];
            
            m.evolveState();
            assertEqual([1 1], m.unprocessedMonomers');
            assertEqual([0 0], m.processedMonomers');
            assertEqual([0 0], m.signalSequenceMonomers');
            assertEqual(0, m.substrates(m.substrateIndexs_water));
            assertEqual(1e3, m.substrates(m.substrateIndexs_PG160));
        end

        function testNoProcessingWithoutSignalPeptidase(this)
            m = this.process;
            m.lipoproteinMonomerIndexs = 1;
            m.secretedMonomerIndexs = 2;
            m.unprocessedMonomerIndexs = [];
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_PG160) = 1e3;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_diacylglycerylTransferase) = 1e3;
            m.unprocessedMonomers = [1;1];
            m.processedMonomers = [0;0];
            m.signalSequenceMonomers = [0;0];
            
            m.evolveState();
            assertEqual([1 1], m.unprocessedMonomers');
            assertEqual([0 0], m.processedMonomers');
            assertEqual([0 0], m.signalSequenceMonomers');
            assertEqual(1e6, m.substrates(m.substrateIndexs_water));
            assertEqual(1e3, m.substrates(m.substrateIndexs_PG160));
        end

        function testNoLipoproteinProcessingWithoutDiacylglycerylTransferase(this)
            m = this.process;
            m.lipoproteinMonomerIndexs = 1;
            m.secretedMonomerIndexs = [];
            m.unprocessedMonomerIndexs = [];
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_PG160) = 1e3;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalPeptidase) = 1e3;
            m.unprocessedMonomers = 1;
            m.processedMonomers = 0;
            m.signalSequenceMonomers = 0;

            m.evolveState();
            assertEqual(1, m.unprocessedMonomers);
            assertEqual(0, m.processedMonomers);
            assertEqual(0, m.signalSequenceMonomers);
            assertEqual(1e6, m.substrates(m.substrateIndexs_water));
            assertEqual(1e3, m.substrates(m.substrateIndexs_PG160));
        end

        % Verifies that signal peptidase processes roughly as many monomers as
        % its rate will allow when it's the limiting factor, and that the
        % monomers chosen for processing are chosen without egregious bias.
        function testLimitedSignalPeptidase_secretedProteinsOnly(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalPeptidase) = 3;
            m.unprocessedMonomers(:) = 0;
            m.unprocessedMonomers(m.secretedMonomerIndexs) = 10;
            m.processedMonomers(:) = 0;
            m.signalSequenceMonomers(:) = 0;

            m.evolveState();
            n = m.enzymes(m.enzymeIndexs_signalPeptidase) * ...
                m.lipoproteinSignalPeptidaseSpecificRate * m.stepSizeSec;
            i = m.secretedMonomerIndexs;
            assertVectorsAlmostEqual(...
                n, sum(m.processedMonomers(i)), 'relative', 0.10);
            assertTrue(10 > max(m.processedMonomers(i)));
            assertTrue(0 < min(m.processedMonomers(i)));
            assertEqual(...
                10 * ones(size(i)), ...
                m.processedMonomers(i) + m.unprocessedMonomers(i));
            assertEqual(m.processedMonomers, m.signalSequenceMonomers);
        end

        % Verifies that signal peptidase processes roughly as many monomers as
        % its rate will allow when it's the limiting factor and there is a mix
        % of lipoproteins and secreted proteins, and that the monomers chosen
        % for processing are chosen without egregious bias.
        function testLimitedSignalPeptidase_lipoproteinsAndSecretedProteins(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_PG160) = 1e3;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalPeptidase) = 20;
            m.enzymes(m.enzymeIndexs_diacylglycerylTransferase) = 1e4;
            m.unprocessedMonomers(:) = 0;
            m.unprocessedMonomers(m.secretedMonomerIndexs) = 10;
            m.unprocessedMonomers(m.lipoproteinMonomerIndexs) = 10;
            m.processedMonomers(:) = 0;
            m.signalSequenceMonomers(:) = 0;

            m.evolveState();
            n = m.enzymes(m.enzymeIndexs_signalPeptidase) * ...
                m.lipoproteinSignalPeptidaseSpecificRate * m.stepSizeSec;
            i = [m.secretedMonomerIndexs;m.lipoproteinMonomerIndexs];
            assertVectorsAlmostEqual(...
                n, sum(m.processedMonomers(i)), 'relative', 0.05);
            assertTrue(10 > max(m.processedMonomers(i)));
            assertTrue(0 < min(m.processedMonomers(i)));
            assertEqual(...
                10 * ones(size(i)), ...
                m.processedMonomers(i) + m.unprocessedMonomers(i));
            assertEqual(m.processedMonomers, m.signalSequenceMonomers);
        end

        % Verifies that diacylglyceryl transferase processes roughly as many
        % lipoproteins as its rate will allow when it's the limiting factor.
        function testLimitedDiacylglycerylTransferase(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_PG160) = 1e3;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalPeptidase) = 1e6;
            m.enzymes(m.enzymeIndexs_diacylglycerylTransferase) = 1e3;
            m.unprocessedMonomers(:) = 0;
            m.unprocessedMonomers(m.lipoproteinMonomerIndexs) = 10;
            m.processedMonomers(:) = 0;
            m.signalSequenceMonomers(:) = 0;

            m.evolveState();
            n = m.enzymes(m.enzymeIndexs_diacylglycerylTransferase) * ...
                m.lipoproteinDiacylglycerylTransferaseSpecificRate * m.stepSizeSec;
            i = m.lipoproteinMonomerIndexs;
            assertElementsAlmostEqual(n, sum(m.processedMonomers(i)), 'absolute', 3);
            assertTrue(10 > max(m.processedMonomers(i)));
            assertTrue(0 < min(m.processedMonomers(i)));
        end

        function testLimitedPG160(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_PG160) = 100;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalPeptidase) = 1e5;
            m.enzymes(m.enzymeIndexs_diacylglycerylTransferase) = 1e5;
            m.unprocessedMonomers(:) = 0;
            m.unprocessedMonomers(m.lipoproteinMonomerIndexs) = 10;
            m.processedMonomers(:) = 0;
            m.signalSequenceMonomers(:) = 0;
            assertTrue(sum(m.unprocessedMonomers(m.lipoproteinMonomerIndexs)) > m.substrates(m.substrateIndexs_PG160)); 

            m.evolveState();
            i = m.lipoproteinMonomerIndexs;
            assertIn(nnz(m.unprocessedMonomers(m.lipoproteinMonomerIndexs)), [1 Inf]);
            assertEqual(100, sum(m.processedMonomers(i)));
        end

        function testLotsOfEverything(this)
            m = this.process;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_PG160) = 1e4;
            m.enzymes(m.enzymeIndexs_signalPeptidase) = 1e3;
            m.enzymes(m.enzymeIndexs_diacylglycerylTransferase) = 1e2;
            i = [m.secretedMonomerIndexs;m.lipoproteinMonomerIndexs];
            m.unprocessedMonomers(i) = randi(100, size(i));
            m.processedMonomers(:) = 0;
            m.signalSequenceMonomers(:) = 0;

            m.evolveState();
            bounds = sort(...
                [m.enzymes(m.enzymeIndexs_signalPeptidase) * ...
                    m.lipoproteinSignalPeptidaseSpecificRate * ...
                    m.stepSizeSec;
                 m.enzymes(m.enzymeIndexs_diacylglycerylTransferase) * ...
                    m.lipoproteinDiacylglycerylTransferaseSpecificRate * ...
                    m.stepSizeSec]);
            assertTrue(0.99 * bounds(1) < sum(m.processedMonomers(i)));
            assertTrue(1.01 * bounds(2) > sum(m.processedMonomers(i)));
        end

        function testGeneEssentiality(this)
            m = this.process;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.substrates(m.substrateIndexs_PG160) = 1e3;
            m.enzymes(m.enzymeIndexs_signalPeptidase) = 1e3;
            m.enzymes(m.enzymeIndexs_diacylglycerylTransferase) = 1e3;
            m.unprocessedMonomers(:) = 1;
            m.processedMonomers(:) = 0;
            this.helpTestGeneEssentiality({
                'MG_086';     %prolipoprotein diacylglyceryl transferase
                'MG_210'},... %prolipoprotein signal peptidase, signal peptidase II
                @(m, i) any(i.processedMonomers(m.unprocessedMonomerIndexs) < ...
                            m.processedMonomers(m.unprocessedMonomerIndexs)) && ...
                        any(i.processedMonomers(m.lipoproteinMonomerIndexs) < ...
                            m.processedMonomers(m.lipoproteinMonomerIndexs)) && ...
                        any(i.processedMonomers(m.secretedMonomerIndexs) < ...
                            m.processedMonomers(m.secretedMonomerIndexs)));  
        end
    end
end
