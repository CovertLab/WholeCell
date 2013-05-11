%Protein processing I process test case
%
% Author: Markus Covert, mcovert@stanford.edu
% Author: Jared Jacobs, jmjacobs@cs.stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef ProteinProcessingI_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = ProteinProcessingI_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end

    methods
        function testOneDeformylation(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_deformylase) = 1;
            m.nascentMonomerNTerminalMethionineCleavages = false;
            m.unprocessedMonomers = 1;
            m.processedMonomers = 0;

            m.evolveState();
            assertEqual(0, m.unprocessedMonomers);
            assertEqual(1, m.processedMonomers);
            assertEqual(0, m.substrates(m.substrateIndexs_water));
            assertEqual(1, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1, m.substrates(m.substrateIndexs_formate));
            assertEqual(0, m.substrates(m.substrateIndexs_methionine));
            assertEqual(1, m.enzymes(m.enzymeIndexs_deformylase));
        end

        function testOneDeformylationAndCleavage(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_deformylase) = 1;
            m.enzymes(m.enzymeIndexs_methionineAminoPeptidase) = 1;
            m.nascentMonomerNTerminalMethionineCleavages = true;
            m.unprocessedMonomers = 1;
            m.processedMonomers = 0;

            m.evolveState();
            assertEqual(0, m.unprocessedMonomers);
            assertEqual(1, m.processedMonomers);
            assertEqual(0, m.substrates(m.substrateIndexs_water));
            assertEqual(1, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1, m.substrates(m.substrateIndexs_formate));
            assertEqual(1, m.substrates(m.substrateIndexs_methionine));
            assertEqual(1, m.enzymes(m.enzymeIndexs_deformylase));
        end

        function testNoDeformylationWithoutDeformylase(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1;
            m.enzymes(:) = 0;
            m.nascentMonomerNTerminalMethionineCleavages = false;
            m.unprocessedMonomers = 1;
            m.processedMonomers = 0;

            m.evolveState();
            assertEqual(1, m.unprocessedMonomers);
            assertEqual(0, m.processedMonomers);
            assertEqual(1, m.substrates(m.substrateIndexs_water));
            assertEqual(0, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(0, m.substrates(m.substrateIndexs_formate));
        end

        function testNoCleavageWithoutMethionineAminoPeptidase(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 2;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_deformylase) = 1;
            m.nascentMonomerNTerminalMethionineCleavages = true;
            m.unprocessedMonomers = 1;
            m.processedMonomers = 0;

            m.evolveState();
            assertEqual(1, m.unprocessedMonomers);
            assertEqual(0, m.processedMonomers);
            assertEqual(2, m.substrates(m.substrateIndexs_water));
            assertEqual(0, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(0, m.substrates(m.substrateIndexs_formate));
            assertEqual(0, m.substrates(m.substrateIndexs_methionine));
        end

        % Verifies that deformylase processes as many monomers as its rate will
        % allow when it's the limiting factor (i.e. water is plentiful).
        function testLimitedDeformylase(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_deformylase) = 1;

            numKinds = m.deformylaseSpecificRate * m.stepSizeSec;
            m.nascentMonomerNTerminalMethionineCleavages = false(numKinds,1);
            m.unprocessedMonomers = ones(numKinds,1);
            m.processedMonomers = zeros(numKinds,1);

            m.evolveState();
            assertEqual(0, sum(m.unprocessedMonomers));
            assertEqual(numKinds, sum(m.processedMonomers));
            assertEqual(1e6 - numKinds, m.substrates(m.substrateIndexs_water));
            assertEqual(numKinds, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(numKinds, m.substrates(m.substrateIndexs_formate));
            assertEqual(0, m.substrates(m.substrateIndexs_methionine));
        end

        % A regression test that covers a weakness of a past implementation.
        % This test is the same as the preceding one except that there are
        % slightly more kinds of monomers. In the past, if a deformylase could
        % not service all kinds of unprocessed monomers, then it would process
        % none.
        function testLimitedDeformylase_moreMonomerKindsThanEnzymeBandwidth(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_deformylase) = 1;

            numKinds = round(1.1 * m.deformylaseSpecificRate * m.stepSizeSec);
            m.nascentMonomerNTerminalMethionineCleavages = false(numKinds,1);
            m.unprocessedMonomers = ones(numKinds,1);
            m.processedMonomers = zeros(numKinds,1);

            m.evolveState();
            expectedNumProcessed = round(m.deformylaseSpecificRate);
            assertElementsAlmostEqual(expectedNumProcessed, sum(m.processedMonomers), 'relative', 0.05);
            assertElementsAlmostEqual(numKinds - expectedNumProcessed, sum(m.unprocessedMonomers), 'absolute', 3);
        end

        % Verifies that deformylase processes as many monomers as its rate will
        % allow when it's the limiting factor and cleaving is required for some
        % kinds of monomers, without favoring any kind of monomer.
        function testLimitedDeformylase_withCleaving(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_deformylase) = 1;
            m.enzymes(m.enzymeIndexs_methionineAminoPeptidase) = 1e6;

            m.deformylaseSpecificRate = 5;
            m.nascentMonomerNTerminalMethionineCleavages = [true;true;false;false];
            m.unprocessedMonomers = [3;1;2;1];
            m.processedMonomers = [0;0;0;0];

            m.evolveState();
            assertEqual([1 0 0 0], m.unprocessedMonomers');
            assertEqual([2 1 2 1], m.processedMonomers');
            assertEqual(1e6-9, m.substrates(m.substrateIndexs_water));
            assertEqual(6, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(6, m.substrates(m.substrateIndexs_formate));
            assertEqual(3, m.substrates(m.substrateIndexs_methionine));
        end

        % Verifies that methionine aminopeptidase processes as many monomers as
        % its rate will allow when it's the limiting factor.
        function testLimitedMethionineAminoPeptidase(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_deformylase) = 1e6;
            m.enzymes(m.enzymeIndexs_methionineAminoPeptidase) = 1;

            limit = m.methionineAminoPeptidaseSpecificRate * m.stepSizeSec;
            m.nascentMonomerNTerminalMethionineCleavages = true(limit,1);
            m.unprocessedMonomers = ones(limit,1);
            m.processedMonomers = zeros(limit,1);

            m.evolveState();
            assertEqual(0, sum(m.unprocessedMonomers));
            assertEqual(limit, sum(m.processedMonomers));
            assertEqual(1e6 - 2*limit, m.substrates(m.substrateIndexs_water));
            assertEqual(limit, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(limit, m.substrates(m.substrateIndexs_formate));
            assertEqual(limit, m.substrates(m.substrateIndexs_methionine));
        end

        function testLotsOfEverything(this)
            m = this.process;
            m.substrates(m.substrateIndexs_water) = 1e6;
            n = 9;
            m.enzymes(m.enzymeIndexs_deformylase) = n;
            m.enzymes(m.enzymeIndexs_methionineAminoPeptidase) = n;
            m.unprocessedMonomers(:) = m.randStream.randi(10, size(m.processedMonomers, 1), 1);
            m.processedMonomers(:) = 0;

            m.evolveState();
            bounds = sort(...
                [n * m.deformylaseSpecificRate * m.stepSizeSec
                 n * m.methionineAminoPeptidaseSpecificRate * m.stepSizeSec]);
            assertIn(sum(m.processedMonomers), [0.9 * bounds(1) Inf]);
            assertIn(sum(m.processedMonomers), [0 1.1 * bounds(2)]);
        end
        
        function testLittleWater(this)
            %case 1
            m = this.process;
            m.substrates(m.substrateIndexs_water) = 3;
            n = 9;
            m.enzymes(m.enzymeIndexs_deformylase) = n;
            m.enzymes(m.enzymeIndexs_methionineAminoPeptidase) = n;
            m.unprocessedMonomers(:) = 1;
            m.processedMonomers(:) = 0;

            m.evolveState();
            
            assertIn(sum(m.processedMonomers), [1 3]);
            
            %case 2
            m = this.process;
            m.substrates(m.substrateIndexs_water) = 3;
            n = 9;
            m.enzymes(m.enzymeIndexs_deformylase) = n;
            m.enzymes(m.enzymeIndexs_methionineAminoPeptidase) = 0;
            m.unprocessedMonomers(:) = 1;
            m.processedMonomers(:) = 0;

            m.evolveState();
            
            assertIn(sum(m.processedMonomers), [1 3]);
        end

        function testGeneEssentiality(this)
            m = this.process;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.enzymes(m.enzymeIndexs_deformylase) = 1;
            m.enzymes(m.enzymeIndexs_methionineAminoPeptidase) = 1;
            m.unprocessedMonomers(:) = 1;
            m.processedMonomers(:) = 0;
            this.helpTestGeneEssentiality({
                'MG_106';     %peptide deformylase
                'MG_172'},... %methionine aminopeptidase, type I
                @(m, i) ...
                    any(i.processedMonomers(~m.nascentMonomerNTerminalMethionineCleavages) < ...
                        m.processedMonomers(~m.nascentMonomerNTerminalMethionineCleavages)) && ...
                    any(i.processedMonomers(m.nascentMonomerNTerminalMethionineCleavages) < ...
                        m.processedMonomers(m.nascentMonomerNTerminalMethionineCleavages)));
        end
    end
end
