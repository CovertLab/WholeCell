%FtsZPolymerization process test case
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/1/2010
classdef FtsZPolymerization_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = FtsZPolymerization_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
        
        function setUp(this)
            this.setUp@edu.stanford.covert.cell.sim.ProcessTestCase();
            warning('off', 'WholeCell:warning');
        end
        
        function loadTestFixture(this)
            this.loadTestFixture@edu.stanford.covert.cell.sim.ProcessTestCase();
            m = this.process;
            m.geometry.volume = 1e-18;
        end
    end

    %tests
    methods
        function testInitializedToSteadyState(this)
            m = this.process;
            
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_FtsZ) = 1000;
            m.boundEnzymes(:);
            this.setGTP(1e5);
            this.setGDP(1e4);

            m.initializeState();
            initialEnzymes = m.enzymes;

            for i = 1:5
                m.evolveState();
            end
            assertNonnegativeCounts(m.enzymes);
            assertEqual(countFtsZ(initialEnzymes), countFtsZ(m.enzymes));
            assertElementsAlmostEqual(initialEnzymes', m.enzymes', 'relative', 0.2, 3);
        end

        function testFtsZAndGtpBalance(this)
            m = this.process;
            this.setGTP(400);
            this.setGDP(0);
            m.enzymes(m.enzymeIndexs_FtsZ) = 300;
            m.enzymes(m.enzymeIndexs_FtsZ_GDP) = 100;
            m.enzymes(m.enzymeIndexs_FtsZ_activated) = [200 101 0 0 0 0 0 0 22]';

            for i = 1:20
                m.evolveState();
                assertNonnegativeCounts(m.enzymes);
                assertEqual(1000, countFtsZ(m.enzymes));
                assertEqual(...
                    m.substrates(m.substrateIndexs_gtp),...
                    sum(m.enzymes([m.enzymeIndexs_FtsZ m.enzymeIndexs_FtsZ_GDP])));
            end
        end

        function testNoFtsZ(this)
            m = this.process;
            m.enzymes(:) = 0;
            this.setGTP(1e5);

            m.evolveState();
            assertEqual(zeros(size(m.enzymes)), m.enzymes);
            assertEqual(1e5, m.substrates(m.substrateIndexs_gtp));
        end

        function testNoGtp(this)
            m = this.process;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_FtsZ) = 1000;
            m.substrates(:) = 1e6;
            this.setGTP(0);
            this.setGDP(1e5);
            
            m.evolveState();
            
            assertTrue(all(m.substrates >= 0));
            assertTrue(all(m.enzymes >= 0));
            assertEqual(1000, [1 1 1:9] * m.enzymes);
            assertAllEqual(0, m.enzymes(m.enzymeIndexs_FtsZ_activated));
            assertEqual(0, m.substrates(m.substrateIndexs_gtp));
        end
        
        function testNoGdp(this)
            m = this.process;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_FtsZ) = 1000;
            m.substrates(:) = 1e6;
            this.setGTP(1e5);
            this.setGDP(0);
            
            m.evolveState();
            
            assertTrue(all(m.substrates >= 0));
            assertTrue(all(m.enzymes >= 0));
            assertEqual(1000, [1 1 1:9] * m.enzymes);
            assertAllEqual(0, m.enzymes(m.enzymeIndexs_FtsZ_GDP));
            assertEqual(0, m.substrates(m.substrateIndexs_gdp));
        end
        
        function testNoWater(this)
            m = this.process;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_FtsZ) = 1000;
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_water) = 0;
            this.setGDP(0);
            
            m.evolveState();
            
            assertTrue(all(m.substrates >= 0));
            assertTrue(all(m.enzymes >= 0));            
            assertEqual(1000, [1 1 1:9] * m.enzymes);
            assertAllEqual(0, m.enzymes(m.enzymeIndexs_FtsZ_GDP));
        end

        function testForwardActivationOnly(this)
            m = this.process;
            m.enzymes(:) = [1000 0 0 0 0 0 0 0 0 0 0]';  %inactive monomer only
            this.setGTP(1000);

            this.setAllReactionRatesToZero();
            m.activationFwd = 10;

            m.evolveState();
            assertEqual(0, m.substrates(m.substrateIndexs_gtp));
            assertEqual(0, m.enzymes(m.enzymeIndexs_FtsZ));
            assertEqual(0, m.enzymes(m.enzymeIndexs_FtsZ_GDP));
            assertEqual([1000 0 0 0 0 0 0 0 0], m.enzymes(m.enzymeIndexs_FtsZ_activated)');
        end

        function testReverseActivationOnly(this)
            m = this.process;
            m.enzymes(:) = [0 0 1000 0 0 0 0 0 0 0 0]';  %activated monomer only
            this.setGTP(0);

            this.setAllReactionRatesToZero();
            m.activationRev = 10;

            m.evolveState();
            assertEqual(1000, m.substrates(m.substrateIndexs_gtp));
            assertEqual(1000, m.enzymes(m.enzymeIndexs_FtsZ));
            assertEqual(0, m.enzymes(m.enzymeIndexs_FtsZ_GDP));
            assertEqual([0 0 0 0 0 0 0 0 0], m.enzymes(m.enzymeIndexs_FtsZ_activated)');
        end

        function testForwardExchangeOnly(this)
            m = this.process;
            m.enzymes(:) = [0 1000 0 0 0 0 0 0 0 0 0]';  %deactivated monomer only
            this.setGTP(1000);
            this.setGDP(0);

            this.setAllReactionRatesToZero();
            m.exchangeFwd = 1e6;

            m.evolveState();
            assertEqual(0, m.substrates(m.substrateIndexs_gtp));
            assertEqual(1000, m.substrates(m.substrateIndexs_gdp));
            assertEqual(0, m.enzymes(m.enzymeIndexs_FtsZ));
            assertEqual(0, m.enzymes(m.enzymeIndexs_FtsZ_GDP));
            assertEqual([1000 0 0 0 0 0 0 0 0], m.enzymes(m.enzymeIndexs_FtsZ_activated)');
        end

        function testReverseExchangeOnly(this)
            m = this.process;
            m.enzymes(:) = [0 0 1000 0 0 0 0 0 0 0 0]';  %activated monomer only
            this.setGTP(0);
            this.setGDP(1000);

            this.setAllReactionRatesToZero();
            m.exchangeRev = 1e6;

            m.evolveState();
            assertEqual(1000, m.substrates(m.substrateIndexs_gtp));
            assertEqual(0, m.substrates(m.substrateIndexs_gdp));
            assertEqual(0, m.enzymes(m.enzymeIndexs_FtsZ));
            assertEqual(1000, m.enzymes(m.enzymeIndexs_FtsZ_GDP));
            assertEqual([0 0 0 0 0 0 0 0 0], m.enzymes(m.enzymeIndexs_FtsZ_activated)');
        end
        
        function testForwardNucleationOnly(this)
            m = this.process;
            m.enzymes = [0 0 1000 0 0 0 0 0 0 0 0]';  %activated monomer only

            this.setAllReactionRatesToZero();
            m.nucleationFwd = 1e6;

            m.evolveState();
            assertElementsAlmostEqual(...
                [0 0 0 500 0 0 0 0 0 0 0], m.enzymes', 'absolute', 2);
            assertEqual(1000, countFtsZ(m.enzymes));
        end

        function testReverseNucleationOnly(this)
            m = this.process;
            m.enzymes = [0 0 0 500 0 0 0 0 0 0 0]';  %activated dimer only

            this.setAllReactionRatesToZero();
            m.nucleationRev = 10;

            m.evolveState();
            assertEqual([0 0 1000 0 0 0 0 0 0 0 0], m.enzymes');
        end

        function testForwardElongationOnly(this)
            m = this.process;
            m.enzymes = [0 0 700 100 0 0 0 0 0 0 0]';  %activated monomers and dimers

            this.setAllReactionRatesToZero();
            m.elongationFwd = 1e7;

            m.evolveState();
            assertElementsAlmostEqual(...
                [0 0 0 0 0 0 0 0 0 0 100], m.enzymes', 'absolute', 9);
            assertEqual(900, countFtsZ(m.enzymes));
        end

        function testReverseElongationOnly(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.enzymes = [0 0 0 0 0 0 0 0 0 0 100]';  %9mers only

            this.setAllReactionRatesToZero();
            m.elongationRev = 100;

            m.evolveState();
            assertElementsAlmostEqual(...
                [0 0 700 100 0 0 0 0 0 0 0], m.enzymes', 'absolute', 2);
            assertEqual(900, countFtsZ(m.enzymes));
        end
        
        function testApplySubstrateLimits(this)
            m = this.process;
            
            %% Case 1: ODEs trying to make FtsZ-GTP polymers from FtsZ-GDP, but unable to because no GTP
            %state before running ODEs
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_FtsZ_GDP) = sum(1:9);
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_gtp) = 0;
            
            %state after running ODEs
            substrates = m.substrates;
            enzymes = zeros(size(m.enzymes));
            enzymes(m.enzymeIndexs_FtsZ_activated) = 1;
            
            %final state
            [finEnzymes, finSubstrates] = m.applySubstrateLimits(enzymes, substrates);
            
            %assertions
            assertTrue(all(finSubstrates >= 0)); %counts are non-negative
            assertTrue(all(finEnzymes >= 0)); %counts are non-negative
            assertEqual([1 1 1:9] * m.enzymes, [1 1 1:9] * finEnzymes); %number of FtsZ monomer is constant
            assertEqual([0 1 1:9] * m.enzymes + [1 1] * m.substrates([m.substrateIndexs_gdp; m.substrateIndexs_gtp]), ... %number of guanosine (free + bound) is constant
                [0 1 1:9] * finEnzymes + [1 1] * finSubstrates([m.substrateIndexs_gdp; m.substrateIndexs_gtp]));
            assertIn([0 0 1:9] * finEnzymes + [0 1] * finSubstrates([m.substrateIndexs_gdp; m.substrateIndexs_gtp]), ... %number of GTP molecules (free + bound) doesn't increase / energy is not created
                [0 [0 0 1:9] * m.enzymes + [0 1] * m.substrates([m.substrateIndexs_gdp; m.substrateIndexs_gtp])]);
            assertEqual(m.substrates' * m.substrateMolecularWeights + m.enzymes' * m.enzymeMolecularWeights, ... %mass is conserved
                finSubstrates' * m.substrateMolecularWeights + finEnzymes' * m.enzymeMolecularWeights);
            
            %% Case 2: Making FtsZ-GDP from free FtsZ and GTP
            %state before running ODEs
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_FtsZ) = 100;
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_gdp) = 0;
            
            %state after running ODEs
            substrates = m.substrates;
            enzymes = zeros(size(m.enzymes));
            enzymes(m.enzymeIndexs_FtsZ_GDP) = 100;
            
            %final state
            [finEnzymes, finSubstrates] = m.applySubstrateLimits(enzymes, substrates);
            
            %assertions
            assertTrue(all(finSubstrates >= 0)); %counts are non-negative
            assertTrue(all(finEnzymes >= 0)); %counts are non-negative
            assertEqual([1 1 1:9] * m.enzymes, [1 1 1:9] * finEnzymes); %number of FtsZ monomer is constant
            assertEqual([0 1 1:9] * m.enzymes + [1 1] * m.substrates([m.substrateIndexs_gdp; m.substrateIndexs_gtp]), ... %number of guanosine (free + bound) is constant
                [0 1 1:9] * finEnzymes + [1 1] * finSubstrates([m.substrateIndexs_gdp; m.substrateIndexs_gtp]));
            assertIn([0 0 1:9] * finEnzymes + [0 1] * finSubstrates([m.substrateIndexs_gdp; m.substrateIndexs_gtp]), ... %number of GTP molecules (free + bound) doesn't increase / energy is not created
                [0 [0 0 1:9] * m.enzymes + [0 1] * m.substrates([m.substrateIndexs_gdp; m.substrateIndexs_gtp])]);
            assertElementsAlmostEqual(m.substrates' * m.substrateMolecularWeights + m.enzymes' * m.enzymeMolecularWeights, ... %mass is conserved
                finSubstrates' * m.substrateMolecularWeights + finEnzymes' * m.enzymeMolecularWeights, ...
                'relative', 1e-8, 0);
        end

        function testGeneEssentiality(this)
            m = this.process;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_FtsZ) = 100;
            m.boundEnzymes(:) = 0;
            this.setGTP(1e5);
            this.setGDP(1e4);
            this.helpTestGeneEssentiality({
                'MG_224'},...  %ftsZ
                @(m,i) m.enzymes(m.enzymeIndexs_FtsZ_9mer) > ...
                       i.enzymes(m.enzymeIndexs_FtsZ_9mer), ...
                       struct('lengthSec', 20));
        end
    end

    %helper methods
    methods
        function setAllReactionRatesToZero(this)
            m = this.process;
            m.activationFwd = 0;
            m.activationRev = 0;
            m.exchangeFwd = 0;
            m.exchangeRev = 0;
            m.nucleationFwd = 0;
            m.nucleationRev = 0;
            m.elongationFwd = 0;
            m.elongationRev = 0;
        end
        
        function setGTP(this, numMolecules)
            m = this.process;
            m.substrates(m.substrateIndexs_gtp) = numMolecules;
        end

        function setGDP(this, numMolecules)
            m = this.process;
            m.substrates(m.substrateIndexs_gdp) = numMolecules;
        end
    end
end

function assertNonnegativeCounts(values)
    assertTrue(all(values >= 0));
    assertEqual(round(values), values);
end

function result = countFtsZ(counts)
    result = [1 1 1:length(counts)-2] * counts;
end