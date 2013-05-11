%ChromosomeSegregation process test case
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/29/2010
classdef ChromosomeSegregation_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = ChromosomeSegregation_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end

        function loadTestFixture(this)
            this.loadTestFixture@edu.stanford.covert.cell.sim.ProcessTestCase();

            c = this.process.chromosome;
            c.initialize();
            c.polymerizedRegions(1, :) = c.polymerizedRegions(1, 1);
            c.linkingNumbers(1, :) = c.linkingNumbers(1, 1);
            c.segregated = false;
        end
        
        function testSuccessfulSegregration(this)
            m = this.process();
            m.enzymes(:) = 1;
            m.substrates(:) = 1;

            m.evolveState();
            assertTrue(m.chromosome.segregated);
            assertEqual(0, m.substrates(m.substrateIndexs_gtp));
            assertEqual(2, m.substrates(m.substrateIndexs_gdp));
            assertEqual(0, m.substrates(m.substrateIndexs_water));
            assertEqual(2, m.substrates(m.substrateIndexs_phosphate));
            assertEqual(2, m.substrates(m.substrateIndexs_hydrogen));
        end

        function testNoCobQ(this)
            m = this.process();
            m.enzymes(:) = 1;
            m.enzymes(m.enzymeIndexs_cobQ) = 0;
            m.substrates(:) = 1;

            m.evolveState();
            assertFalse(m.chromosome.segregated);
        end

        function testNoMraZ(this)
            m = this.process();
            m.enzymes(:) = 1;
            m.enzymes(m.enzymeIndexs_mraZ) = 0;
            m.substrates(:) = 1;

            m.evolveState();
            assertFalse(m.chromosome.segregated);
        end

        function testNoObg(this)
            m = this.process();
            m.enzymes(:) = 1;
            m.enzymes(m.enzymeIndexs_obg) = 0;
            m.substrates(:) = 1;

            m.evolveState();
            assertFalse(m.chromosome.segregated);
        end

        function testNoEra(this)
            m = this.process();
            m.enzymes(:) = 1;
            m.enzymes(m.enzymeIndexs_era) = 0;
            m.substrates(:) = 1;

            m.evolveState();
            assertFalse(m.chromosome.segregated);
        end

        function testNoWater(this)
            m = this.process();
            m.enzymes(:) = 1;
            m.substrates(:) = 1;
            m.substrates(m.substrateIndexs_water) = 0;

            m.evolveState();
            assertFalse(m.chromosome.segregated);
        end
        
        function testNoGtp(this)
            m = this.process();
            m.enzymes(:) = 1;
            m.substrates(:) = 1;
            m.substrates(m.substrateIndexs_gtp) = 0;

            m.evolveState();
            assertFalse(m.chromosome.segregated);
        end

        function testGeneEssentiality(this)
            this.process.chromosome.segregated = false;
            
            this.process.enzymes(:) = 1;
            this.process.substrates(:) = 1;            
            
            this.helpTestGeneEssentiality({
                'MG_470'      %CobQ/CobB/MinD/ParA nucleotide binding domain
                'MG_221'      %mraZ protein
                'MG_387'      %GTP-binding protein Era
                'MG_384'      %GTPase1 Obg
                'MG_203'      %DNA topoisomerase IV
                'MG_204'      %DNA topoisomerase IV
                }, ...
                @(m,~) m.chromosome.segregated,...
                struct('lengthSec', 1));
        end
    end
end
