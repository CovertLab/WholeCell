classdef Translation_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = Translation_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end

        function testSimpleFixture(this)
            m = this.process;

            %% initial conditions: one mRNA, components of one ribosome
            m.stimuli(:) = 0;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_gtp) = 1e6;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.boundEnzymes(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_translationFactors) = 1;
            m.enzymes(m.enzymeIndexs_initiationFactor3) = 0;
            m.enzymes(m.enzymeIndexs_elongationGFactor) = 2;
            m.enzymes(m.enzymeIndexs_ribosome30SIF3) = 1;
            m.enzymes(m.enzymeIndexs_ribosome50S) = 1;
            m.ribosome.states = [];
            m.ribosome.boundMRNAs = [];
            m.ribosome.mRNAPositions = [];
            m.ribosome.tmRNAPositions = [];
            m.polypeptide.boundMRNAs = [];
            m.polypeptide.nascentMonomerLengths = [];
            m.polypeptide.proteolysisTagLengths = [];

            m.polypeptide.monomerAASequences = {'MKKKKKKKKKKKKKKKKK'};
            m.polypeptide.monomerTRNASequences = {[10 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20]'};
            m.polypeptide.monomerLengths = 18;
            m.monomers = 0;

            m.mRNAs = 1;
            m.aminoacylatedTRNAs = repmat(20, size(m.rna.matureTRNAIndexs));
            m.aminoacylatedTMRNA = 0;

            %% initiation
            m.evolveState();
            assertEqual(m.ribosome.activeValue, m.ribosome.states(1));
            assertEqual(0, m.polypeptide.nascentMonomerLengths(1));
            assertEqual(0, m.polypeptide.proteolysisTagLengths(1));
            assertEqual(1, m.polypeptide.boundMRNAs(1));
            assertEqual(1, m.ribosome.boundMRNAs(1));
            assertEqual(0, m.monomers(1));
            assertEqual(zeros(20,1), m.substrates(m.substrateIndexs_aminoAcids));
            assertEqual(1e6-1, m.substrates(m.substrateIndexs_gtp));
            assertEqual(1, m.substrates(m.substrateIndexs_gdp));
            assertEqual(1, m.substrates(m.substrateIndexs_phosphate));
            assertEqual(1e6-1, m.substrates(m.substrateIndexs_water));
            assertEqual(1, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1, m.enzymes(m.enzymeIndexs_initiationFactor1));
            assertEqual(1, m.enzymes(m.enzymeIndexs_initiationFactor2));
            assertEqual(1, m.enzymes(m.enzymeIndexs_initiationFactor3));
            assertEqual([2;1;1;1], m.enzymes(m.enzymeIndexs_elongationFactors));
            assertEqual(1, m.enzymes(m.enzymeIndexs_terminationFactor));
            assertEqual(1, m.enzymes(m.enzymeIndexs_recyclingFactor));
            assertEqual(repmat(20, size(m.rna.matureTRNAIndexs)), m.aminoacylatedTRNAs);
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome30S));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome30SIF3));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome50S));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome70S));
            assertEqual(1, m.boundEnzymes(m.enzymeIndexs_ribosome70S));            

            %% elongation
            m.evolveState();
            assertEqual(m.ribosome.activeValue, m.ribosome.states(1));
            assertEqual(16, m.polypeptide.nascentMonomerLengths(1));
            assertEqual(0, m.polypeptide.proteolysisTagLengths(1));
            assertEqual(1, m.polypeptide.boundMRNAs(1));
            assertEqual(1, m.ribosome.boundMRNAs(1));
            assertEqual(0, m.monomers(1));
            assertEqual(zeros(20,1), m.substrates(m.substrateIndexs_aminoAcids));
            assertEqual(1e6-33, m.substrates(m.substrateIndexs_gtp));
            assertEqual(33, m.substrates(m.substrateIndexs_gdp));
            assertEqual(33, m.substrates(m.substrateIndexs_phosphate));
            assertEqual(1e6-34, m.substrates(m.substrateIndexs_water));
            assertEqual(49, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1, m.enzymes(m.enzymeIndexs_initiationFactor1));
            assertEqual(1, m.enzymes(m.enzymeIndexs_initiationFactor2));
            assertEqual(1, m.enzymes(m.enzymeIndexs_initiationFactor3));
            assertEqual([1;0;0;0], m.enzymes(m.enzymeIndexs_elongationFactors));
            assertEqual([1;1;1;1], m.boundEnzymes(m.enzymeIndexs_elongationFactors));
            assertEqual(1, m.enzymes(m.enzymeIndexs_terminationFactor));
            assertEqual(1, m.enzymes(m.enzymeIndexs_recyclingFactor));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome30S));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome30SIF3));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome50S));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome70S));
            assertEqual(1, m.boundEnzymes(m.enzymeIndexs_ribosome70S));
            assertEqual([19 20 20 20 20 20 20 20 20 20 5]', m.aminoacylatedTRNAs(10:20));

            %% elongation + termination
            m.evolveState();
            assertEqual(m.ribosome.notExistValue, m.ribosome.states(1));
            assertEqual(0, m.polypeptide.nascentMonomerLengths(1));
            assertEqual(0, m.polypeptide.proteolysisTagLengths(1));
            assertEqual(0, m.polypeptide.boundMRNAs(1));
            assertEqual(0, m.ribosome.boundMRNAs(1));
            assertEqual(1, m.monomers(1));
            assertEqual(zeros(20,1), m.substrates(m.substrateIndexs_aminoAcids));
            assertEqual(1e6-39, m.substrates(m.substrateIndexs_gtp));
            assertEqual(39, m.substrates(m.substrateIndexs_gdp));
            assertEqual(39, m.substrates(m.substrateIndexs_phosphate));
            assertEqual(1e6-40, m.substrates(m.substrateIndexs_water));
            assertEqual(57, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1, m.enzymes(m.enzymeIndexs_initiationFactor1));
            assertEqual(1, m.enzymes(m.enzymeIndexs_initiationFactor2));
            assertEqual(1, m.enzymes(m.enzymeIndexs_initiationFactor3));
            assertEqual([1;0;0;0], m.enzymes(m.enzymeIndexs_elongationFactors));
            assertEqual([1;1;1;1], m.boundEnzymes(m.enzymeIndexs_elongationFactors));
            assertEqual(1, m.enzymes(m.enzymeIndexs_ribosome30S));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome30SIF3));
            assertEqual(1, m.enzymes(m.enzymeIndexs_ribosome50S));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome70S));
            assertEqual(0, m.boundEnzymes(m.enzymeIndexs_ribosome70S));
            assertEqual([19 20 20 20 20 20 20 20 20 20 3]', m.aminoacylatedTRNAs(10:20));

            %% freeing the elongation factors
            m.evolveState();
            assertEqual([2;1;1;1], m.enzymes(m.enzymeIndexs_elongationFactors));
            assertEqual([0;0;0;0], m.boundEnzymes(m.enzymeIndexs_elongationFactors));
        end

        function testInitiationFromUncomplexedRibosome30SAndIF3(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_gtp) = 1;
            m.substrates(m.substrateIndexs_water) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_initiationFactors) = 1;
            m.enzymes(m.enzymeIndexs_ribosome30S) = 1;
            m.enzymes(m.enzymeIndexs_ribosome50S) = 1;
            m.boundEnzymes(:) = 0;
            m.mRNAs(:) = 1;
            m.ribosome.states(:) = m.ribosome.notExistValue;
            m.evolveState();
            assertTrue(any(m.ribosome.states==m.ribosome.activeValue),...
                'initiation failed');
        end

        function testInitiation(this)
            this.initializeProcessForSuccessfulInitiation();
            m = this.process;

            m.evolveState();
            assertTrue(any(m.ribosome.states==m.ribosome.activeValue),...
                'initiation failed');
        end

        function testNoInitiationWithoutAvailableEnergy(this)
            this.initializeProcessForSuccessfulInitiation();
            m = this.process;
            m.substrates(m.substrateIndexs_gtp) = 0;

            m.evolveState();
            assertTrue(~any(m.ribosome.states == m.ribosome.activeValue),...
                'initiation happened');
        end

        function testNoInitiationWithoutMRNAs(this)
            this.initializeProcessForSuccessfulInitiation();
            m = this.process;
            m.mRNAs(:) = 0;

            m.evolveState();
            assertTrue(~any(m.ribosome.states==m.ribosome.activeValue),...
                'initiation happened');
        end

        function testNoInitiationWithOneEnzymeMissing(this)
            m = this.process;
            essential_enzyme_indexs = [...
                m.enzymeIndexs_initiationFactor1,...
                m.enzymeIndexs_initiationFactor2,...
                m.enzymeIndexs_ribosome50S];
            for i=1:numel(essential_enzyme_indexs)
                this.initializeProcessForSuccessfulInitiation();
                m.enzymes(essential_enzyme_indexs(i)) = 0;

                m.evolveState();
                assertTrue(...
                    ~any(m.ribosome.states==m.ribosome.activeValue),...
                    ['initiation happened, i=' num2str(i)]);
            end

            % Ribosome 30S has two forms
            this.initializeProcessForSuccessfulInitiation();
            m.enzymes(m.enzymeIndexs_ribosome30S) = 0;
            m.enzymes(m.enzymeIndexs_ribosome30SIF3) = 0;

            m.evolveState();
            assertTrue(...
                ~any(m.ribosome.states == m.ribosome.activeValue),...
                'initiation happened');
        end

        function testElongation(this)
            this.initializeProcessForSuccessfulElongation(3);
            m = this.process;

            m.evolveState();
            assertEqual(3, m.polypeptide.nascentMonomerLengths);
            assertEqual(1, m.polypeptide.boundMRNAs);
            assertEqual(1, m.ribosome.boundMRNAs);
            assertEqual(0, m.monomers(1));
            assertTrue(~any(m.enzymes(m.enzymeIndexs_elongationFactors)));
            assertTrue(all(m.boundEnzymes(m.enzymeIndexs_elongationFactors)));
        end

        function testNoElongationWithoutEnoughEnergy(this)
            this.initializeProcessForSuccessfulElongation(3);
            m = this.process;
            m.substrates(m.substrateIndexs_gtp) = 1;  % 2 GTP required per AA

            m.evolveState();
            assertEqual(0, m.polypeptide.nascentMonomerLengths);
            assertEqual(1, m.polypeptide.boundMRNAs);
            assertEqual(1, m.ribosome.boundMRNAs);
            assertEqual(0, m.monomers(1));
            assertTrue(all(m.enzymes(m.enzymeIndexs_elongationFactors)));
            assertTrue(~any(m.boundEnzymes(m.enzymeIndexs_elongationFactors)));
        end
        
        function testNoElongationWithoutWater(this)
            this.initializeProcessForSuccessfulElongation(3);
            m = this.process;
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_water) = 0;

            m.evolveState();
            assertEqual(0, m.polypeptide.nascentMonomerLengths);
            assertEqual(1, m.polypeptide.boundMRNAs);
            assertEqual(1, m.ribosome.boundMRNAs);
            assertEqual(0, m.monomers(1));
            assertTrue(all(m.enzymes(m.enzymeIndexs_elongationFactors)));
            assertTrue(~any(m.boundEnzymes(m.enzymeIndexs_elongationFactors)));
            assertEqual(0, m.substrates(m.substrateIndexs_water));
        end

        function testNoSubstrates(this)
            this.initializeProcessForSuccessfulElongation(3);
            
            m = this.process;
            m.substrates(:) = 0;

            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            initial_aminoacylatedTRNAs = m.aminoacylatedTRNAs;
            initial_ribosomeStates = m.ribosome.states;
            initial_boundMRNAs = m.ribosome.boundMRNAs;
            initial_monomers = m.monomers;
            initial_monomerLengths = m.polypeptide.monomerLengths;

            m.evolveState();

            %no substrate used
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_aminoacylatedTRNAs, m.aminoacylatedTRNAs);

            %all active ribosomes remain active and don't advance
            assertEqual(...
                initial_ribosomeStates(initial_ribosomeStates == m.ribosome.activeValue),...
                m.ribosome.states(initial_ribosomeStates == m.ribosome.activeValue));
            assertEqual(...
                initial_boundMRNAs(initial_ribosomeStates == m.ribosome.activeValue),...
                m.ribosome.boundMRNAs(initial_ribosomeStates == m.ribosome.activeValue));
            assertEqual(...
                initial_boundMRNAs(initial_ribosomeStates == m.ribosome.activeValue),...
                m.polypeptide.boundMRNAs(initial_ribosomeStates == m.ribosome.activeValue));
            
            %no enzyme gain or loss or binding
            initial_enzymes(m.enzymeIndexs_ribosome30SIF3) = 0;
            initial_enzymes(m.enzymeIndexs_ribosome30S) = 0;
            initial_enzymes(m.enzymeIndexs_initiationFactor3) = 0;
            
            assertEqual(initial_enzymes', m.enzymes');
            assertEqual(initial_boundEnzymes', m.boundEnzymes');

            %no additional monomers or nascent monomer elongation
            assertEqual(initial_monomers, m.monomers);
            assertEqual(initial_monomerLengths, m.polypeptide.monomerLengths);
        end

        function testManyTranscripts(this)
            m = this.process;

            %% initial conditions: three mRNAs, components of three ribosomes
            m.stimuli(:) = 0;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_gtp) = 1e6;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.boundEnzymes(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_translationFactors) = 6;
            m.enzymes(m.enzymeIndexs_elongationGFactor) = 12;
            m.enzymes(m.enzymeIndexs_ribosome30S) = 6;
            m.enzymes(m.enzymeIndexs_ribosome50S) = 6;
            m.ribosome.states = zeros(0, 1);
            m.ribosome.boundMRNAs = zeros(0, 1);
            m.ribosome.mRNAPositions = zeros(0, 1);
            m.ribosome.tmRNAPositions = zeros(0, 1);
            m.polypeptide.boundMRNAs = zeros(0, 1);
            m.polypeptide.nascentMonomerLengths = zeros(0, 1);
            m.polypeptide.proteolysisTagLengths = zeros(0, 1);

            m.polypeptide.monomerAASequences = {
                'MMKRTYQPSKLKRAKTH';
                'MGFMARMATAQGRKVLR';
                'MQRRFKNRAQLTVSSER'};
            m.polypeptide.monomerTRNASequences = {
                [10 7 20 17 28 22 21 6 3 20 19 27 13 1 20 4 33]';
                [10 14 12 7 1 15 7 1 4 1 21 18 13 20 29 19 15]';
                [10 21 36 13 12 20 32 13 1 21 34 30 29 24 3 31 13]'};

            m.polypeptide.monomerLengths = [17 17 17]';
            m.monomers = [0 0 0]';

            m.mRNAs = [1 1 1]';
            m.aminoacylatedTRNAs = repmat(100, size(m.rna.matureTRNAIndexs));
            m.aminoacylatedTMRNA = 0;

            % 3 initiated, 3 idle
            m.evolveState();
            assertEqual(3, numel(find(m.ribosome.states==m.ribosome.activeValue)));
            assertTrue(~any(m.polypeptide.nascentMonomerLengths));
            assertEqual(0, sum(m.monomers));

            % 3 elongating, 3 initiated
            m.evolveState();
            assertEqual(6, numel(find(m.ribosome.states==m.ribosome.activeValue)));
            assertEqual([16 16 16 0 0 0]', m.polypeptide.nascentMonomerLengths);
            assertEqual(0, sum(m.monomers));

            % 3 terminated, 3 elongating
            m.evolveState();
            assertEqual(3, numel(find(m.ribosome.states==m.ribosome.activeValue)));
            assertEqual([0 0 0 16 16 16]', m.polypeptide.nascentMonomerLengths);
            assertEqual(3, sum(m.monomers));

            % 3 initiated, 3 terminated
            m.evolveState();
            assertEqual(3, numel(find(m.ribosome.states==m.ribosome.activeValue)));
            assertEqual([0 0 0 0 0 0]', m.polypeptide.nascentMonomerLengths);
            assertEqual(6, sum(m.monomers));

            % 3 elongating, 3 initiated
            m.evolveState();
            assertEqual(6, numel(find(m.ribosome.states==m.ribosome.activeValue)));
            assertEqual([16 16 16 0 0 0]', m.polypeptide.nascentMonomerLengths);
            assertEqual(6, sum(m.monomers));

            assertEqual(6, numel(find(m.ribosome.boundMRNAs)));
            assertEqual(6, numel(find(m.polypeptide.boundMRNAs)));
            assertEqual(100 * numel(m.rna.matureTRNAIndexs) - 150, sum(m.aminoacylatedTRNAs));

            % substrate usage = initiation + elongation + termination
            assertEqual(1e6-(12 + 2 * 150 + 2 * 6), m.substrates(m.substrateIndexs_gtp));
            assertEqual(12 + 2 * 150 + 2 * 6, m.substrates(m.substrateIndexs_gdp));
            assertEqual(12 + 2 * 150 + 2 * 6, m.substrates(m.substrateIndexs_phosphate));
            assertEqual(1e6-(12 + 2*150 + 3 * 6 + 3), m.substrates(m.substrateIndexs_water));
            assertEqual(12 + 3 * 150 + 2 * 6, m.substrates(m.substrateIndexs_hydrogen));

            assertEqual([6;6;6], m.enzymes(m.enzymeIndexs_initiationFactors));
            assertEqual([9;3;3;3], m.enzymes(m.enzymeIndexs_elongationFactors));
            assertEqual([3;3;3;3], m.boundEnzymes(m.enzymeIndexs_elongationFactors));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome30S));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome30SIF3));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome50S));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome70S));
            assertEqual(6, m.boundEnzymes(m.enzymeIndexs_ribosome70S));
        end

        function testBindingProbabilities(this)
            m = this.process;

            %plentiful substrate and enzyme and a small number of ribosomes
            numRibosomes = 50;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.enzymes(m.enzymeIndexs_ribosome30S) = 0;
            m.enzymes(m.enzymeIndexs_ribosome30SIF3) = numRibosomes;
            m.enzymes(m.enzymeIndexs_ribosome50S) = numRibosomes;
            m.enzymes(m.enzymeIndexs_ribosome70S) = 0;

            %no bound ribosomes
            m.boundEnzymes(:) = 0;
            m.ribosome.states = repmat(m.ribosome.notExistValue, numRibosomes, 1);
            m.ribosome.boundMRNAs = zeros(numRibosomes, 1);
            m.polypeptide.boundMRNAs = zeros(numRibosomes, 1);
            m.polypeptide.nascentMonomerLengths = zeros(numRibosomes, 1);
            m.polypeptide.proteolysisTagLengths = zeros(numRibosomes, 1);
            m.ribosome.mRNAPositions = m.polypeptide.nascentMonomerLengths;
            m.ribosome.tmRNAPositions = m.polypeptide.proteolysisTagLengths;

            %binding probabilities
            m.mRNAs(1:numRibosomes) = 1;
            m.mRNAs(numRibosomes+1:end) = 0;

            m.evolveState();
            assertTrue(all(m.ribosome.states == m.ribosome.activeValue), 'not all ribosomes bound');
            assertEqual(1:numRibosomes, unique(m.ribosome.boundMRNAs)', 'not all mRNAs bound');
            assertEqual(1:numRibosomes, unique(m.polypeptide.boundMRNAs)', 'not all mRNAs bound');
        end

        function testRibosomeStalling(this)
            m = this.process;

            %% initial conditions: one ribosome bound to one mRNA
            m.stimuli(:) = 0;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_gtp) = 1e6;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_translationFactors) = 1;
            m.enzymes(m.enzymeIndexs_elongationGFactor) = 2;
            m.enzymes(m.enzymeIndexs_tmRNABindingProtein) = 1;
            m.enzymes(m.enzymeIndexs_peptidylTRNAHydrolase) = 1;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_ribosome70S) = 1;
            m.boundTMRNA = 0;

            m.polypeptide.monomerAASequences = {'MKKKKKVVVVV'};
            m.polypeptide.monomerTRNASequences = {[10 20 20 20 20 20 29 29 29 29 29]'};
            m.polypeptide.monomerLengths = 11;
            m.monomers = 0;

            m.mRNAs = 1;
            m.ribosome.boundMRNAs = 1;
            m.polypeptide.boundMRNAs = 1;
            m.ribosome.states = m.ribosome.activeValue;
            m.ribosome.mRNAPositions = 0;
            m.ribosome.tmRNAPositions = 0;
            m.polypeptide.nascentMonomerLengths = 0;
            m.polypeptide.proteolysisTagLengths = 0;
            m.aminoacylatedTRNAs = repmat(20, size(m.rna.matureTRNAIndexs));
            m.aminoacylatedTRNAs(29) = 0;  % withhold valine
            m.aminoacylatedTMRNA = 1;
            m.tmRNABindingProbability = 1;

            %% elongation
            m.evolveState();
            assertEqual(m.ribosome.activeValue, m.ribosome.states(1));
            assertEqual(6, m.polypeptide.nascentMonomerLengths(1));
            assertEqual(0, m.polypeptide.proteolysisTagLengths(1));
            assertEqual(1, m.polypeptide.boundMRNAs(1));
            assertEqual(1, m.ribosome.boundMRNAs(1));
            assertEqual(0, m.monomers(1));

            %% stalling
            m.evolveState();
            assertEqual(m.ribosome.stalledValue, m.ribosome.states(1));
            assertEqual(6, m.polypeptide.nascentMonomerLengths(1));
            assertEqual(1, m.polypeptide.proteolysisTagLengths(1));
            assertEqual(1, m.polypeptide.boundMRNAs(1));
            assertEqual(1, m.ribosome.boundMRNAs(1));
            assertEqual(0, m.monomers(1));
            assertEqual(0, m.aminoacylatedTMRNA);
            assertEqual(0, m.enzymes(m.enzymeIndexs_tmRNA));
            assertEqual(1, m.boundTMRNA);
            assertEqual(zeros(20,1), m.substrates(m.substrateIndexs_aminoAcids));
            assertEqual(1e6 - 12, m.substrates(m.substrateIndexs_gtp));
            assertEqual(12, m.substrates(m.substrateIndexs_gdp));
            assertEqual(12, m.substrates(m.substrateIndexs_phosphate));
            assertEqual(1e6 - 12 + 6 -6 -1, m.substrates(m.substrateIndexs_water));
            assertEqual(12+6, m.substrates(m.substrateIndexs_hydrogen));

            %% elongation and termination of proteolysis tag
            m.aminoacylatedTRNAs(29) = 20;  % provide valine again
            m.evolveState();
            m.evolveState();
            assertEqual(1, numel(m.polypeptide.abortedSequences));
            assertEqual([1;0;0;0], m.enzymes(m.enzymeIndexs_elongationFactors));
            assertEqual([1;1;1;1], m.boundEnzymes(m.enzymeIndexs_elongationFactors));
            assertEqual(1, m.enzymes(m.enzymeIndexs_ribosome30S));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome30SIF3));
            assertEqual(1, m.enzymes(m.enzymeIndexs_ribosome50S));
            assertEqual(0, m.enzymes(m.enzymeIndexs_ribosome70S));
            assertEqual(0, m.boundEnzymes(m.enzymeIndexs_ribosome70S));
            assertEqual([19 18 18 20 20 20 20 20 20 20 14], m.aminoacylatedTRNAs(10:20)');

            %% freeing the elongation factors
            m.evolveState();
            assertEqual([2;1;1;1], m.enzymes(m.enzymeIndexs_elongationFactors));
            assertEqual([0;0;0;0], m.boundEnzymes(m.enzymeIndexs_elongationFactors));
        end
        
        function testDetailedMassBalance(this)
            m = this.process;
            
            m.mRNAs(:) = 0;
            m.mRNAs(1) = 1;
            
            m.enzymes = m.enzymes + m.boundEnzymes;
            m.boundEnzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_ribosome30S) = 1;
            m.enzymes(m.enzymeIndexs_ribosome30SIF3) = 0;
            m.enzymes(m.enzymeIndexs_ribosome50S) = 1;
            m.enzymes(m.enzymeIndexs_ribosome70S) = 0;
            
            m.ribosome.boundMRNAs = 0;
            m.polypeptide.boundMRNAs = 0;
            m.ribosome.states = m.ribosome.notExistValue;
            m.ribosome.mRNAPositions = 0;
            m.ribosome.tmRNAPositions = 0;
            m.polypeptide.nascentMonomerLengths = 0;
            m.polypeptide.proteolysisTagLengths = 0;
            
            m.substrates(:) = 1e6;
            m.aminoacylatedTRNAs(:) = 1e6;
            
            initSubstrates = m.substrates;
            
            for i = 1:200;
                m.evolveState();
                if m.monomers(1) == 1
                    break;
                end
            end
            
            finalSubstrates = initSubstrates;
            gtpCost = 3 + 2 * m.monomer.lengths(1);
            finalSubstrates(m.substrateIndexs_gtp) = finalSubstrates(m.substrateIndexs_gtp) - gtpCost;
            finalSubstrates(m.substrateIndexs_gdp) = finalSubstrates(m.substrateIndexs_gdp) + gtpCost;
            finalSubstrates(m.substrateIndexs_phosphate) = finalSubstrates(m.substrateIndexs_phosphate) + gtpCost;
            finalSubstrates(m.substrateIndexs_hydrogen) = finalSubstrates(m.substrateIndexs_hydrogen) ...
                + gtpCost ...           %energy
                + m.monomer.lengths(1); %peptidyl-tRNA bond
            finalSubstrates(m.substrateIndexs_water) = finalSubstrates(m.substrateIndexs_water) ...
                - gtpCost ...               %energy
                - m.monomer.lengths(1) ...  %peptidyl-tRNA bond
                + (m.monomer.lengths(1)-1); %peptide bond
            
            assertEqual(1, m.monomers(1));
            assertEqual(0, m.ribosome.boundMRNAs);
            assertEqual(finalSubstrates(setdiff(1:end, m.substrateIndexs_aminoAcids))', ...
                m.substrates(setdiff(1:end, m.substrateIndexs_aminoAcids))');
            assertEqual(finalSubstrates, m.substrates);
        end
        
        function testNoBiasedPausing(this)
            m = this.process;
            rib = m.ribosome;
            pol = m.polypeptide;
            
            %replace tRNA sequences with codons corresponding to limiting
            %tRNA
            pol.monomerTRNASequences = cellfun(@(seq) ones(1, 200), pol.monomerTRNASequences, 'UniformOutput', false);
            
            %all ribosomes
            %- active
            %- translating first mRNA at first codon which we've set to
            %  correspond to the limiting tRNA
            nRibs = 10;
            rib.states = repmat(rib.activeValue, [nRibs 1]);
            rib.boundMRNAs = ones(nRibs, 1);
            rib.mRNAPositions = zeros(nRibs, 1);
            rib.tmRNAPositions = zeros(nRibs, 1);
            
            %synchronize polypepide and ribosome states
            pol.boundMRNAs = rib.boundMRNAs;
            pol.nascentMonomerLengths = rib.mRNAPositions;
            pol.proteolysisTagLengths = rib.tmRNAPositions;
            pol.abortedPolypeptides = zeros(0, 3);
            
            %synchronize ribosome and protein complex states
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_ribosome70S) = nRibs;
            m.enzymes(m.enzymeIndexs_translationFactors) = 100;
            m.boundEnzymes(m.enzymeIndexs_translationFactors) = 100;
            
            %simulate
            iterMax = 100;
            for i = 1:iterMax
                %lot of tRNA, but 1 tRNA limiting
                m.aminoacylatedTRNAs(:) = 100;
                m.aminoacylatedTRNAs(1) = 1;
                m.freeTRNAs(:) = 0;
                
                %simulate
                m.evolveState();
            end
            
            %assert that all ribosomes move
            assertEqual(nRibs, nnz(rib.mRNAPositions));
            assertIn(range(rib.mRNAPositions), [0 20]);
        end
        
        function testGeneEssentiality(this)
            m = this.process;
            
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            numRibosomes = 10;
            m.enzymes(m.enzymeIndexs_ribosome30SIF3) = 0;
            m.enzymes(m.enzymeIndexs_ribosome30S) = numRibosomes;
            m.enzymes(m.enzymeIndexs_ribosome50S) = numRibosomes;
            m.enzymes(m.enzymeIndexs_ribosome70S) = 0;
            m.boundEnzymes(m.enzymeIndexs_ribosome70S) = 0;
            m.aminoacylatedTRNAs = repmat(1e6, size(m.rna.matureTRNAIndexs));
            m.aminoacylatedTMRNA = 1e6;
            m.mRNAs(:) = 1;           
            m.monomers(:) = 0;

            m.ribosome.states = zeros(0, 1);
            m.ribosome.boundMRNAs = zeros(0, 1);
            m.ribosome.mRNAPositions = zeros(0, 1);
            m.ribosome.tmRNAPositions = zeros(0, 1);
            m.polypeptide.boundMRNAs = zeros(0, 1);
            m.polypeptide.nascentMonomerLengths = zeros(0, 1);
            m.polypeptide.proteolysisTagLengths = zeros(0, 1);
            
            this.helpTestGeneEssentiality({
                'MG_173';   %translation initiation factor IF-1
                'MG_142';   %translation initiation factor IF-2
                'MG_196';   %translation initiation factor IF-3
                'MG_089';   %translation elongation factor G
                'MG_026';   %translation elongation factor P
                'MG_451';   %translation elongation factor Tu
                'MG_433';   %translation elongation factor Ts
                'MG_435';   %ribosome recycling factor
                'MG_473';   %ribosomal protein L33 type 2
                'MG_070';   %ribosomal protein S2
                'MG_081';   %ribosomal protein L11
                'MG_082';   %ribosomal protein L1
                'MG_087';   %ribosomal protein S12
                'MG_088';   %ribosomal protein S7
                'MG_090';   %ribosomal protein S6
                'MG_092';   %ribosomal protein S18
                'MG_093';   %ribosomal protein L9
                'MGrrnA5S';   %5S ribosomal rRNA
                'MGrrnA16S';  %16S ribosomal rRNA
                'MGrrnA23S';  %23S ribosomal rRNA
                'MG_150';   %ribosomal protein S10
                'MG_151';   %ribosomal protein L3
                'MG_152';   %ribosomal protein L4/L1 family
                'MG_153';   %ribosomal protein L23
                'MG_154';   %ribosomal protein L2
                'MG_155';   %ribosomal protein S19
                'MG_156';   %ribosomal protein L22
                'MG_157';   %ribosomal protein S3
                'MG_158';   %ribosomal protein L16
                'MG_159';   %ribosomal protein L29
                'MG_160';   %ribosomal protein S17
                'MG_161';   %ribosomal protein L14
                'MG_162';   %ribosomal protein L24
                'MG_163';   %ribosomal protein L5
                'MG_164';   %ribosomal protein S14
                'MG_165';   %ribosomal protein S8
                'MG_166';   %ribosomal protein L6
                'MG_167';   %ribosomal protein L18
                'MG_168';   %ribosomal protein S5
                'MG_169';   %ribosomal protein L15
                'MG_174';   %ribosomal protein L36
                'MG_175';   %ribosomal protein S13
                'MG_176';   %ribosomal protein S11
                'MG_178';   %ribosomal protein L17
                'MG_197';   %ribosomal protein L35
                'MG_198';   %ribosomal protein L20
                'MG_481';   %30S ribosomal protein S21
                'MG_232';   %ribosomal protein L21
                'MG_234';   %ribosomal protein L27
                'MG_257';   %ribosomal protein L31
                'MG_258';   %peptide chain release factor 1
                'MG_311';   %ribosomal protein S4
                'MG_325';   %ribosomal protein L33
                'MG_361';   %ribosomal protein L10
                'MG_362';   %ribosomal protein L7/L12
                'MG_363';   %ribosomal protein L32
                'MG_522';   %ribosomal protein S20
                'MG_417';   %ribosomal protein S9
                'MG_418';   %ribosomal protein L13
                'MG_424';   %ribosomal protein S15
                'MG_426';   %ribosomal protein L28
                'MG_444';   %ribosomal protein L19
                'MG_446';   %ribosomal protein S16
                'MG_466';   %ribosomal protein L34
                },...
                @(m,i) any(m.monomers > i.monomers),...
                struct('lengthSec', 75));
        end
    end

    methods
        function initializeProcessForSuccessfulInitiation(this)
            m = this.process;
            m.substrates(m.substrateIndexs_gtp) = 1e6;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.enzymes(:) = 10;
            m.boundEnzymes(:) = 0;
            m.aminoacylatedTRNAs = repmat(20, size(m.rna.matureTRNAIndexs));
            m.mRNAs(:) = 1;
            m.ribosome.states(:) = m.ribosome.notExistValue;
        end

        function initializeProcessForSuccessfulElongation(this, numAAs)
            m = this.process;
            m.substrates(m.substrateIndexs_gtp) = 2 * numAAs;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.mRNAs(:) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_elongationFactors) = 1;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_ribosome70S) = 1;
            m.ribosome.boundMRNAs = 1;
            m.polypeptide.boundMRNAs = 1;
            m.ribosome.states = m.ribosome.activeValue;
            m.ribosome.mRNAPositions = 0;
            m.ribosome.tmRNAPositions = 0;
            m.polypeptide.nascentMonomerLengths = 0;
            m.polypeptide.proteolysisTagLengths = 0;
            m.aminoacylatedTRNAs(:) = numAAs;
        end
    end
end
