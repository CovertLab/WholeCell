%Ribosome assembly test cases
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/20/2010
classdef RibosomeAssembly_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = RibosomeAssembly_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end

        function testConstants(this)
            m = this.process;
            
            %each subunit participates in just one complex
            assertAllEqual(1, sum(m.proteinComplexRNAComposition, 2));
            assertAllEqual(1, sum(m.proteinComplexMonomerComposition, 2));
        end
        
        function testSuccessfulAssembly(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.RNAs(:) = 100;
            m.monomers(:) = 100;
            m.complexs(:) = 0;

            m.evolveState();
            assertEqual([100 100], m.complexs');
            assertAllEqual(0, m.RNAs');
            assertAllEqual(0, m.monomers');
            assertAllEqual(1e6, m.enzymes');
            assertEqual(1e6-600, m.substrates(m.substrateIndexs_water));
            assertEqual(1e6-600, m.substrates(m.substrateIndexs_gtp));
            assertEqual(1e6+600, m.substrates(m.substrateIndexs_gdp));
            assertEqual(1e6+600, m.substrates(m.substrateIndexs_phosphate));
            assertEqual(1e6+600, m.substrates(m.substrateIndexs_hydrogen));
        end
        
        function testNoWater(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_water) = 0;
            m.enzymes(:) = 1e6;
            m.RNAs(:) = 100;
            m.monomers(:) = 100;
            m.complexs(:) = 0;

            m.evolveState();
            assertEqual([0 0], m.complexs');
        end

        function testNoEnergy(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_gtp) = 0;
            m.enzymes(:) = 1e6;
            m.RNAs(:) = 100;
            m.monomers(:) = 100;
            m.complexs(:) = 0;

            m.evolveState();
            assertEqual([0 0], m.complexs');
        end

        function testNoEnzymes(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.RNAs(:) = 1e6;
            m.monomers(:) = 1e6;
            m.enzymes(:) = 0;
            m.complexs(:) = 0;

            m.evolveState();
            assertEqual([0 0], m.complexs');
        end

        function testNoRNAs(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.RNAs(:) = 0;
            m.monomers(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.complexs(:) = 0;

            m.evolveState();
            assertEqual([0 0], m.complexs');
        end

        function testNoMonomers(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.RNAs(:) = 1e6;
            m.monomers(:) = 0;
            m.enzymes(:) = 1e6;
            m.complexs(:) = 0;

            m.evolveState();
            assertEqual([0 0], m.complexs');
        end

        function testLimiting30SRNA(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.RNAs(:) = 1e3;
            m.RNAs(strcmp(m.rnaWholeCellModelIDs,'MGrrnA16S'),:) = 10;
            m.monomers(:) = 1e3;
            m.enzymes(:) = 1e6;
            m.complexs(:) = 0;

            m.evolveState();
            assertEqual(10, m.complexs(m.complexIndexs_30S_ribosome));
            assertTrue(0 < m.complexs(m.complexIndexs_50S_ribosome));
        end

        function testLimiting50SRNA(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.RNAs(:) = 1e6;
            m.RNAs(strcmp(m.rnaWholeCellModelIDs,'MGrrnA5S'),:) = 10;
            m.monomers(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.complexs(:) = 0;

            m.evolveState();
            assertTrue(0 < m.complexs(m.complexIndexs_30S_ribosome));
            assertEqual(10, m.complexs(m.complexIndexs_50S_ribosome));
        end

        function testLimiting30SMonomer(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.RNAs(:) = 1e3;
            m.monomers(:) = 1e3;
            m.monomers(strcmp(m.monomerWholeCellModelIDs,'MG_070_MONOMER'),:) = 10;
            m.enzymes(:) = 1e6;
            m.complexs(:) = 0;

            m.evolveState();
            assertEqual(10, m.complexs(m.complexIndexs_30S_ribosome));
            assertTrue(0 < m.complexs(m.complexIndexs_50S_ribosome));
        end

        function testLimiting50SMonomer(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.RNAs(:) = 1e6;
            m.monomers(:) = 1e6;
            m.monomers(strcmp(m.monomerWholeCellModelIDs,'MG_081_MONOMER'),:) = 10;
            m.enzymes(:) = 1e6;
            m.complexs(:) = 0;

            m.evolveState();
            assertTrue(0 < m.complexs(m.complexIndexs_30S_ribosome));
            assertEqual(10, m.complexs(m.complexIndexs_50S_ribosome));
        end

        function testNo30SEnzyme(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.RNAs(:) = 1e6;
            m.monomers(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.enzymes(m.enzymeIndexs_30S_assembly_gtpase(1)) = 0;
            m.complexs(:) = 0;

            m.evolveState();
            assertEqual(0, m.complexs(m.complexIndexs_30S_ribosome));
            assertTrue(0 < m.complexs(m.complexIndexs_50S_ribosome));
        end

        function testNo50SEnzyme(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.RNAs(:) = 1e6;
            m.monomers(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.enzymes(m.enzymeIndexs_50S_assembly_gtpase(1)) = 0;
            m.complexs(:) = 0;

            m.evolveState();
            assertTrue(0 < m.complexs(m.complexIndexs_30S_ribosome));
            assertEqual(0, m.complexs(m.complexIndexs_50S_ribosome));
        end

        function testLimitingEnergy(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_gtp) = ...
                length(m.enzymeIndexs_30S_assembly_gtpase) + ...
                length(m.enzymeIndexs_50S_assembly_gtpase);
            m.RNAs(:) = 1e6;
            m.monomers(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.complexs(:) = 0;

            m.evolveState();
            assertEqual([1 1], m.complexs');
            assertEqual(0, m.substrates(m.substrateIndexs_gtp));
        end

        function testGeneEssentiality(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.RNAs(:) = 1e3;
            m.monomers(:) = 1e3;
            m.enzymes(:) = 1e6;
            this.helpTestGeneEssentiality({
                'MG_329';     %GTP-binding protein engA
                'MG_335';     %GTP-binding protein engB, putative
                'MG_387';     %GTP-binding protein Era
                'MG_384';     %GTPase1 Obg
                'MG_143';     %ribosome-binding factor A
                'MG_442'},... %ribosomal biogenesis GTPase
                @(m,i) all(i.complexs < m.complexs));
        end
    end
end
