% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/30/2010
classdef RNADecay_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = RNADecay_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end
    
    %tests
    methods
        function testConstants(this)
            m = this.process;
            r = m.rna;
            assertAllEqual(true, isfinite(r.decayRates(r.aminoacylatedIndexs)));
            assertEqual(m.substrateIndexs_water, find(any(m.decayReactions < 0, 1))); %water is only substrate required to decay RNA
        end
              
        function testFMethionineDecay(this)
            m = this.process;
            MG488 = 2223;  %aminoacylated FMET tRNA
            assertEqual('MG488', m.rna.wholeCellModelIDs{MG488});
            assertEqual(1, m.decayReactions(MG488, m.substrateIndexs_fmethionine));

            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ribonucleaseR) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidylTRNAHydrolase) = 1e6;
            m.RNAs(:) = 0;
            m.RNAs(MG488) = 1;
            m.rna.decayRates(:) = 1e6;

            m.evolveState();
            assertEqual(1, m.substrates(m.substrateIndexs_fmethionine));
            assertEqual(0, m.substrates(m.substrateIndexs_methionine));
            assertEqual(0, m.substrates(m.substrateIndexs_formate));
            assertEqual(77, sum(m.substrates(m.substrateIndexs({
                'AMP';'CMP';'GMP';'UMP';'s4UMP';'m7GMP';'PSIURIMP'}))));
            assertEqual(77, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1e6-77, m.substrates(m.substrateIndexs_water));
        end

        function testGlutamineDecay(this)
            m = this.process;
            MG502 = 2275;  %aminoacylated GLN tRNA
            assertEqual('MG502', m.rna.wholeCellModelIDs{MG502});
            assertEqual(1, m.decayReactions(MG502, m.substrateIndexs_glutamine));

            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ribonucleaseR) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidylTRNAHydrolase) = 1e6;
            m.RNAs(:) = 0;
            m.RNAs(MG502) = 1;
            m.rna.decayRates(:) = 1e6;

            m.evolveState();
            assertEqual(1, m.substrates(m.substrateIndexs_glutamine));
            assertEqual(0, m.substrates(m.substrateIndexs_glutamate));
            assertEqual(0, m.substrates(m.substrateIndexs_ammonia));
            assertEqual(75, sum(m.substrates(m.substrateIndexs({
                'AMP';'CMP';'GMP';'UMP';'s4UMP';'cmnm5s2UMP';'PSIURIMP'}))));
            assertEqual(75, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1e6-75, m.substrates(m.substrateIndexs_water));
        end

        function testImpartialityRegardingKindsOfRNA(this)
            m = this.process;
            notAminoacylatedRNAIndexs = true(size(m.RNAs));
            notAminoacylatedRNAIndexs(m.rna.aminoacylatedIndexs) = false;
            notAminoacylatedRNAIndex = find(notAminoacylatedRNAIndexs, 1, 'first');
            aminoacylatedRNAIndex = m.rna.aminoacylatedIndexs(1);
            m.stepSizeSec = 1;
            m.peptidylTRNAHydrolaseSpecificRate = 1;
            m.rna.decayRates(:) = .5;
            m.enzymes(:) = 1000;
            m.RNAs(:) = 0;
            m.RNAs([notAminoacylatedRNAIndex aminoacylatedRNAIndex]) = 1000;
            m.substrates(m.substrateIndexs_water) = 1e7;

            m.evolveState();
            assertVectorsAlmostEqual(500, m.RNAs(aminoacylatedRNAIndex), 'absolute', 30);
            assertVectorsAlmostEqual(500, m.RNAs(notAminoacylatedRNAIndex), 'absolute', 30);
        end
        
        function testNoWater(this)
            m = this.process;
            m.substrates(:) = 1e7;
            m.substrates(m.substrateIndexs_water) = 0;
            m.enzymes(:) = 1e6;
            m.RNAs(:) = 1;
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_RNAs = m.RNAs;
            
            m.evolveState();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_RNAs, m.RNAs);
        end
        
        function testNoRibonuclease(this)
            m = this.process;
            m.enzymes(:) = 1e6;
            m.enzymes(m.enzymeIndexs_ribonucleaseR) = 0;
            m.RNAs(:) = 1;
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_RNAs = m.RNAs;

            m.evolveState();
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_RNAs, m.RNAs);
        end

        function testNoHydrolase(this)
            m = this.process;
            m.enzymes(:) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidylTRNAHydrolase) = 0;
            m.RNAs(:) = 1;
            initial_enzymes = m.enzymes;
            initial_RNAs = m.RNAs;

            m.evolveState();
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_RNAs(m.rna.aminoacylatedIndexs), m.RNAs(m.rna.aminoacylatedIndexs));
        end

        function testHydrolaseRateZero(this)
            m = this.process;
            m.peptidylTRNAHydrolaseSpecificRate = 0;
            m.enzymes(:) = 1e6;
            m.RNAs(:) = 1;
            initial_enzymes = m.enzymes;
            initial_RNAs = m.RNAs;

            m.evolveState();
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_RNAs(m.rna.aminoacylatedIndexs), m.RNAs(m.rna.aminoacylatedIndexs));
        end

        function testNoRNAs(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.RNAs(:) = 0;
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_RNAs = m.RNAs;

            m.evolveState();
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_RNAs, m.RNAs);
        end

        function testDecayRatesZero(this)
            m = this.process;
            m.rna.decayRates(:) = 0;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.RNAs(:) = 1e6;
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_RNAs = m.RNAs;

            m.evolveState();
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_RNAs, m.RNAs);
        end

        function testDecayRatesHigh(this)
            m = this.process;
            m.rna.decayRates(:) = 1e6;
            m.substrates(:) = 1e7;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ribonucleaseR) = 1;
            m.enzymes(m.enzymeIndexs_peptidylTRNAHydrolase) = 800;
            m.RNAs(:) = 1;
            m.RNAs(m.rna.boundIndexs(m.rna.matureTMRNAIndexs), :) = 0;

            m.evolveState();
            assertAllEqual(0, m.RNAs);
            assertEqual(1, m.enzymes(m.enzymeIndexs_ribonucleaseR));
            assertEqual(800, m.enzymes(m.enzymeIndexs_peptidylTRNAHydrolase));
        end
        
        function testDecayAminoacylatedRNAs(this)
            m = this.process;
            m.rna.decayRates(:) = 1e6;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e7;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ribonucleaseR) = 1;
            m.enzymes(m.enzymeIndexs_peptidylTRNAHydrolase) = 1;
            m.peptidylTRNAHydrolaseSpecificRate = 1;
            
            counts = zeros(2, 1);
            for i = 1:200
                m.RNAs(:) = 0;
                m.RNAs(m.rna.aminoacylatedIndexs(1:2)) = [3; 1];
                m.evolveState();
                assertEqual(3, sum(m.RNAs(m.rna.aminoacylatedIndexs(1:2))));
                counts = counts + [1; 3] .* m.RNAs(m.rna.aminoacylatedIndexs(1:2));
            end
            
            assertIn(range(counts) / max(counts), [0 0.15]);
        end
        
        function testDecayAbortedTranscripts(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.transcripts.transcriptionUnitSequences{1}(1:6) = 'AGUCGU';
            m.transcripts.abortedTranscripts = [1 6];
            m.enzymes(:) = 1;
            m.RNAs(:) = 0;
            
            m.evolveState();
            
            assertEqual(cell(0, 1), m.transcripts.abortedSequences);
            assertEqual(1e6 + [1; 1; 2; 2], m.substrates(m.substrateIndexs_nmps));
            assertEqual(1e6-5, m.substrates(m.substrateIndexs_water));
            assertEqual(1e6+5, m.substrates(m.substrateIndexs_hydrogen));
        end
        
        function testExpectations(this)
            m = this.process;
            
            m.substrates(:) = 1e12;
            
            m.RNAs(:) = 0;
            m.RNAs(1) = 1e6;
            RNAs0 = m.RNAs;
            
            iterMax = 100;
            for i = 1:iterMax
                m.evolveState();
            end
            
            assertElementsAlmostEqual(RNAs0 .* exp(-m.rna.decayRates * iterMax), ...
                m.RNAs, 'relative', 1e-3);
        end

        function testGeneEssentiality(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e7;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ribonucleaseR) = 1e6;
            m.enzymes(m.enzymeIndexs_peptidylTRNAHydrolase) = 1e6;
            m.RNAs(:) = 1;
            m.RNAs(m.rna.boundIndexs(m.rna.matureTMRNAIndexs), :) = 0;
            m.transcripts.transcriptionUnitSequences{1}(1:6) = 'AGUCGU';
            m.transcripts.abortedTranscripts = [1 6];
            m.rna.decayRates(:) = 1e6;
            this.helpTestGeneEssentiality({
                'MG_104';     %ribonuclease R
                'MG_083'},... %peptidyl-tRNA hydrolase
                @this.isProperlyFunctioning);
        end
    end

    %helper functions
    methods
        function result = isProperlyFunctioning(~, m, i)
            notAminoacylatedRNAIndexs = true(size(m.RNAs));
            notAminoacylatedRNAIndexs(m.rna.aminoacylatedIndexs) = false;
            result = ...
                isempty(m.transcripts.abortedSequences) && ...
                any(i.RNAs(notAminoacylatedRNAIndexs) > ...
                    m.RNAs(notAminoacylatedRNAIndexs)) && ...
                any(i.RNAs(m.rna.aminoacylatedIndexs) > ...
                    m.RNAs(m.rna.aminoacylatedIndexs));
        end
    end
end
