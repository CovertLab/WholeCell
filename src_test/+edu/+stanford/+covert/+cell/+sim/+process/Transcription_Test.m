%Transcription test case
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford
% University
% Last updated: 7/13/2010
classdef Transcription_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = Transcription_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end
    
    %tests
    methods
        function testSimpleFixture(this)
            this.prepareSimpleFixture();
            
            m = this.process;
            
            initial_RNAs = m.RNAs;
            
            for i = 1:1000
                m.evolveState();
                if any(initial_RNAs < m.RNAs)
                    break;
                end
            end
            
            assertTrue(any(initial_RNAs < m.RNAs), 'No RNAs produced during simulation');
            totEnzymes = [3 3 3 3 5 0]';
            
            eIdxs = setdiff(1:numel(m.enzymes), [m.enzymeIndexs_sigmaFactor; m.enzymeIndexs_rnaPolymerase; m.enzymeIndexs_rnaPolymeraseHoloenzyme]);
            assertEqual(totEnzymes(eIdxs), m.enzymes(eIdxs) + m.boundEnzymes(eIdxs));
            
            eIdxs = [m.enzymeIndexs_rnaPolymerase; m.enzymeIndexs_rnaPolymeraseHoloenzyme];
            assertEqual(sum(totEnzymes(eIdxs)), sum(m.enzymes(eIdxs) + m.boundEnzymes(eIdxs)));
            
            eIdxs = [m.enzymeIndexs_sigmaFactor; m.enzymeIndexs_rnaPolymeraseHoloenzyme];
            assertEqual(sum(totEnzymes(eIdxs)), sum(m.enzymes(eIdxs) + m.boundEnzymes(eIdxs)));
        end
        
        %check that transcription units are not so close that RNA
        %polymerase might get stuck by steric hindrance unable to move in
        %opposite directions on adjacent transcription units
        function testConstants(this)
            m = this.process;
            r = m.rnaPolymerases;
            
            assertEqual(13, numel(m.complexIndexs_DnaA_ATP));
            assertEqual(1, numel(m.transcriptionUnitIndexs_DnaAR12345Boxes));
            
            assertEqual(m.enzymeDNAFootprints(m.enzymeIndexs_rnaPolymerase), m.enzymeDNAFootprints(m.enzymeIndexs_rnaPolymeraseHoloenzyme));
            
            [~, d] = eig(m.stateTransitionProbabilities);
            idx = find(abs(diag(d) - 1) < sqrt(eps));
            assertEqual(1, numel(idx));
            assertElementsAlmostEqual(r.stateExpectations, m.stateTransitionProbabilities * r.stateExpectations, ...
                'relative', 1e-1, 1e-1);
        end
        
        function testNoSubstrates(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            
            m.substrates(:) = 0;
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            initial_states = r.states;
            initial_rnaPolymeraseBoundTranscriptionUnits = t.boundTranscriptionUnits;
            initial_RNAs = m.RNAs;
            
            m.evolveState();
            
            %no substrate used
            assertEqual(initial_substrates, m.substrates);
            
            %all actively transcribing polymerases remain actively
            %transcribing, and don't advance
            assertEqual(...
                initial_states(initial_states >= r.activelyTranscribingValue), ...
                r.states(initial_states >= r.activelyTranscribingValue));
            assertEqual(...
                initial_rnaPolymeraseBoundTranscriptionUnits(initial_states >= r.activelyTranscribingValue), ...
                t.boundTranscriptionUnits(initial_states >= r.activelyTranscribingValue));
            
            %elongation, termination, antitermination factors can't bind
            assertEqual(initial_enzymes(m.enzymeIndexs_elongationFactor), m.enzymes(m.enzymeIndexs_elongationFactor));
            assertEqual(initial_enzymes(m.enzymeIndexs_terminationFactor), m.enzymes(m.enzymeIndexs_terminationFactor));
            assertEqual(initial_enzymes(m.enzymeIndexs_antiterminationFactor), m.enzymes(m.enzymeIndexs_antiterminationFactor));
            
            %no net enzyme gain or loss
            eIdxs = setdiff(1:numel(m.enzymes), [m.enzymeIndexs_sigmaFactor; m.enzymeIndexs_rnaPolymerase; m.enzymeIndexs_rnaPolymeraseHoloenzyme]);
            assertEqual(initial_enzymes(eIdxs) + initial_boundEnzymes(eIdxs), m.enzymes(eIdxs) + m.boundEnzymes(eIdxs));
            
            eIdxs = [m.enzymeIndexs_rnaPolymerase; m.enzymeIndexs_rnaPolymeraseHoloenzyme];
            assertEqual(sum(initial_enzymes(eIdxs) + initial_boundEnzymes(eIdxs)), sum(m.enzymes(eIdxs) + m.boundEnzymes(eIdxs)));
            
            eIdxs = [m.enzymeIndexs_sigmaFactor; m.enzymeIndexs_rnaPolymeraseHoloenzyme];
            assertEqual(sum(initial_enzymes(eIdxs) + initial_boundEnzymes(eIdxs)), sum(m.enzymes(eIdxs) + m.boundEnzymes(eIdxs)));
            
            %no additional RNAs
            assertEqual(initial_RNAs, m.RNAs);
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
        end
        
        function testNoWater(this)
            m = this.process;
            r = m.rnaPolymerases;
            
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_water) = 0;
            
            initial_substrates = m.substrates;
            initial_states = r.states;
            initial_RNAs = m.RNAs;
            
            for i = 1:10
                m.evolveState();
            end
            
            %no water used
            assertEqual(initial_substrates(m.substrateIndexs_water), m.substrates(m.substrateIndexs_water));
            
            %all actively transcribing polymerases remain actively
            %transcribing
            assertTrue(all(ismember(...
                find(initial_states >= r.activelyTranscribingValue), ...
                find(r.states >= r.activelyTranscribingValue))));  
            
            %no additional RNAs
            assertEqual(initial_RNAs, m.RNAs);
        end
        
        function testNoEnzymes(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            c.initialize();
            
            m.enzymes(:)=0;
            m.boundEnzymes(:) = 0;
            r.states(:) = r.notExistValue;
            r.positionStrands(:) = 0;
            t.boundTranscriptProgress(:) = t.nullTranscriptValue;
            t.boundTranscriptChromosome(:) = t.nullTranscriptValue;
            t.boundTranscriptionUnits(:) = t.nullTranscriptValue;
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            initial_states = r.states;
            initial_rnaPolymeraseBoundTranscriptionUnits = t.boundTranscriptionUnits;
            initial_RNAs = m.RNAs;
            
            m.evolveState();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(initial_states, r.states);
            assertEqual(initial_rnaPolymeraseBoundTranscriptionUnits, t.boundTranscriptionUnits);
            assertEqual(initial_RNAs, m.RNAs);
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
        end
        
        function testSaturatingSubstrates(this)
            %process
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            c.initialize();
            
            %saturating substrates and enzymes
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            
            %make all RNA polymerases active, and near end of transcript
            numRNAPolymerase = 5;
            m.boundEnzymes(m.enzymeIndexs_rnaPolymerase) = 0;
            m.boundEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = numRNAPolymerase;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            
            t.boundTranscriptionUnits = 2 * m.randStream.randsample(floor(numel(c.transcriptionUnitStrands)/2), numRNAPolymerase, false, ones(floor(numel(c.transcriptionUnitStrands)/2), 1));
            r.states = t.transcriptionUnitLengths(t.boundTranscriptionUnits) - 10;
            r.positionStrands = [...
                t.transcriptionUnitFivePrimeCoordinates(t.boundTranscriptionUnits) + ...
                (1-t.transcriptionUnitDirections(t.boundTranscriptionUnits)) .* ...
                (t.transcriptionUnitLengths(t.boundTranscriptionUnits) - 10) ...
                c.transcriptionUnitStrands(t.boundTranscriptionUnits)];
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptChromosome = ones(numRNAPolymerase, 1);
            
            m.bindProteinToChromosome(r.positionStrands, m.enzymeIndexs_rnaPolymerase);
            
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            
            %store initial state
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            initial_rnaPolymeraseBoundTranscriptionUnits = t.boundTranscriptionUnits;
            initial_RNAs = m.RNAs;
            
            %predict final state
            newRNACounts = histc(initial_rnaPolymeraseBoundTranscriptionUnits, 1:numel(m.RNAs));
            final_RNAs = initial_RNAs;
            final_RNAs = final_RNAs + newRNACounts;
            
            %evolve
            m.evolveState();
            
            %all actively transcribing polymerases terminate
            assertAllEqual(r.freeValue, r.states);
            assertAllEqual(0, t.boundTranscriptionUnits);
            
            %no net enzyme gain or loss
            assertEqual(initial_enzymes + initial_boundEnzymes, m.enzymes + m.boundEnzymes);
            
            %RNAs correctly incremented
            assertEqual(final_RNAs, m.RNAs);
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
        end
        
        function testNoCTP(this)
            %with CTP
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            
            m.evolveState();
            RNAs=m.RNAs;
            substrates = m.substrates;
            
            %no CTP
            this.loadTestFixture();
            
            m.substrates(m.substrateIndexs_ntp(2))=0;
            m.evolveState();
            assertEqual(0, m.substrates(m.substrateIndexs_ntp(2)));
            
            assertTrue(all(RNAs <= m.RNAs));
            assertTrue(all(substrates(m.substrateIndexs_ntp(setdiff(1:end, 2))) >= ...
                m.substrates(m.substrateIndexs_ntp(setdiff(1:end, 2)))));
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
        end
        
        function testNoElongation(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            m.chromosome.initialize();
            
            %set elongation rate to zero, set positions of all bound RNA
            %polymerases to beginning of transciption unit
            m.rnaPolymeraseElongationRate = 0;
            m.enzymes = m.enzymes + m.boundEnzymes;
            m.boundEnzymes(:) = 0;
            
            r.states(:) = r.freeValue;
            r.positionStrands(:) = 0;
            t.boundTranscriptProgress(:) = t.nullTranscriptValue;
            t.boundTranscriptionUnits(:) = t.nullTranscriptValue;
            t.boundTranscriptChromosome(:) = t.nullTranscriptValue;
            
            iTU = find(t.transcriptionUnitDirections == 1, 5, 'first');
            r.states(1:5) = r.activelyTranscribingValue + 1;
            r.positionStrands(1:5, 1) = t.transcriptionUnitFivePrimeCoordinates(iTU) + 1;
            r.positionStrands(1:5, 2) = 1;
            t.boundTranscriptProgress(1:5) = r.states(1:5);
            t.boundTranscriptionUnits(1:5) = iTU;
            t.boundTranscriptChromosome(1:5) = 1;
            
            m.bindProteinToChromosome(r.positionStrands(1:5, :), m.enzymeIndexs_rnaPolymerase);
            
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            
            %store initial state for later comparison
            initial_RNAs = m.RNAs;
            initial_states = r.states;
            initial_rnaPolymeraseBoundTranscriptionUnits = t.boundTranscriptionUnits;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            
            %evolve
            m.evolveState();
            
            %no RNA polymerases advance
            assertEqual(...
                initial_states(initial_states >= r.activelyTranscribingValue), ...
                r.states(initial_states >= r.activelyTranscribingValue));
            assertEqual(...
                initial_rnaPolymeraseBoundTranscriptionUnits(initial_states >= r.activelyTranscribingValue), ...
                t.boundTranscriptionUnits(initial_states >= r.activelyTranscribingValue));
            
            %elongation, termination, antitermination factors can't bind
            assertEqual(initial_enzymes(m.enzymeIndexs_elongationFactor), m.enzymes(m.enzymeIndexs_elongationFactor));
            assertEqual(initial_enzymes(m.enzymeIndexs_terminationFactor), m.enzymes(m.enzymeIndexs_terminationFactor));
            assertEqual(initial_enzymes(m.enzymeIndexs_antiterminationFactor), m.enzymes(m.enzymeIndexs_antiterminationFactor));
            
            %no net enzyme gain or loss
            eIdxs = setdiff(1:numel(m.enzymes), [m.enzymeIndexs_sigmaFactor; m.enzymeIndexs_rnaPolymerase; m.enzymeIndexs_rnaPolymeraseHoloenzyme]);
            assertEqual(initial_enzymes(eIdxs) + initial_boundEnzymes(eIdxs), m.enzymes(eIdxs) + m.boundEnzymes(eIdxs));
            
            eIdxs = [m.enzymeIndexs_rnaPolymerase; m.enzymeIndexs_rnaPolymeraseHoloenzyme];
            assertEqual(sum(initial_enzymes(eIdxs) + initial_boundEnzymes(eIdxs)), sum(m.enzymes(eIdxs) + m.boundEnzymes(eIdxs)));
            
            eIdxs = [m.enzymeIndexs_sigmaFactor; m.enzymeIndexs_rnaPolymeraseHoloenzyme];
            assertEqual(sum(initial_enzymes(eIdxs) + initial_boundEnzymes(eIdxs)), sum(m.enzymes(eIdxs) + m.boundEnzymes(eIdxs)));
            
            %no additional RNAs
            assertEqual(initial_RNAs, m.RNAs);
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
        end
        
        function testFastElongation(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            c.initialize();
            
            %high elongation rate, saturating substrates and enzymes
            m.rnaPolymeraseElongationRate = 1e5;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            
            %start fixed number of RNA polymerases, all in actively bound
            %state at second base of transcription unit
            numRNAPolymerases = 50;
            m.enzymes(m.enzymeIndexs_rnaPolymerase, :) = numRNAPolymerases;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme, :) = 0;
            m.boundEnzymes(m.enzymeIndexs_rnaPolymerase, :) = 0;
            m.boundEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme, :) = 0;
            
            iTU = find(t.transcriptionUnitDirections == 1, numRNAPolymerases * 2, 'first');
            iTU = iTU(1:2:end);
            r.states = repmat(r.activelyTranscribingValue + 1, numRNAPolymerases, 1);
            r.positionStrands = [...
                t.transcriptionUnitFivePrimeCoordinates(iTU) + 1, ...
                ones(numRNAPolymerases, 1)];
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptionUnits = iTU;
            t.boundTranscriptChromosome = ones(numRNAPolymerases, 1);
            
            m.bindProteinToChromosome(r.positionStrands, m.enzymeIndexs_rnaPolymerase);
            
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            
            %compute what the updated counts of RNAs should be after call
            %to evolve state
            newRNACounts = histc(t.boundTranscriptionUnits, 1:numel(m.RNAs));
            final_RNAs = m.RNAs;
            final_RNAs = final_RNAs + newRNACounts;
            
            %evolve
            m.evolveState();
            
            %all actively transcribing RNA polymerases terminate
            assertTrue(all(r.states == r.freeValue));
            assertTrue(all(t.boundTranscriptionUnits == 0));
            
            %RNAs correctly incremented
            assertEqual(final_RNAs, m.RNAs);
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
        end
        
        function testElongationChromosomeIntegration(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            c.initialize();
            
            %saturating substrates and enzymes
            m.rnaPolymeraseElongationRate = 10;
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            
            %start fixed number of RNA polymerases, all in actively bound
            %state at second base of transcription unit
            numRNAPolymerases = 50;
            m.enzymes(m.enzymeIndexs_rnaPolymerase, :) = 0;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme, :) = 0;
            m.boundEnzymes(m.enzymeIndexs_rnaPolymerase, :) = numRNAPolymerases;
            m.boundEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme, :) = 0;
            
            iTU = find(t.transcriptionUnitDirections == 1, numRNAPolymerases, 'first');
            r.states = repmat(r.activelyTranscribingValue + 1, numRNAPolymerases, 1);
            r.positionStrands = [...
                t.transcriptionUnitFivePrimeCoordinates(iTU),...
                ones(numRNAPolymerases,1)];
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptionUnits = iTU;
            t.boundTranscriptChromosome = ones(numRNAPolymerases, 1);
            
            m.bindProteinToChromosome(r.positionStrands, m.enzymeIndexs_rnaPolymerase);
            
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            
            %evolve
            m.evolveState();
            
            rPolyGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymerase);
            initStates = r.states(r.states >= r.activelyTranscribingValue);
            initBoundSites = find(c.complexBoundSites == rPolyGblIdx);
            
            %evolve
            m.evolveState();
            
            %assert proper movement along choromosome
            finalStates = r.states(r.states >= r.activelyTranscribingValue);
            finalBoundSites = find(c.complexBoundSites == rPolyGblIdx);
            deltaBoundSites = abs(finalBoundSites(:,1) - initBoundSites(:,1));
            assertEqual(sum(finalStates-initStates), sum(deltaBoundSites));
        end
        
        function testRNAPolymeraseStateExpectations_free(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            c.initialize();
            
            %put all RNA polymerases in non-specifically bound state
            numRNAPolymerases = 50;
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = numRNAPolymerases;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            m.boundEnzymes(:) = 0;
            
            %put all RNA polymerases in non-specifically bound state
            r.states = repmat(r.nonSpecificallyBoundValue, numRNAPolymerases, 1);
            r.positionStrands = [1000*(1:numRNAPolymerases)' ones(numRNAPolymerases, 1)];
            t.boundTranscriptionUnits = repmat(t.nullTranscriptValue, numRNAPolymerases, 1);
            t.boundTranscriptProgress = repmat(t.nullTranscriptValue, numRNAPolymerases, 1);
            t.boundTranscriptChromosome = repmat(t.nullTranscriptValue, numRNAPolymerases, 1);
            
            m.bindProteinToChromosome(r.positionStrands, m.enzymeIndexs_rnaPolymerase);
            
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            
            %force state transition
            r.stateExpectations = zeros(4, 1);
            r.stateExpectations(r.nonSpecificallyBoundIndex) = 1;
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.nonSpecificallyBoundIndex, :) = 1;
            
            %store initial state
            initial_states = r.states;
            
            %evolve
            m.evolveState();
            
            %all free RNA polymerases advance to non-specifically bound state
            assertTrue(all(r.states(initial_states == r.freeValue) == ...
                r.nonSpecificallyBoundValue));
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
        end
        
        function testRNAPolymeraseStateExpectations_nonSpecificallyBound(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            c.initialize();
            
            %saturating amounts of all enzymes except RNA polymerase
            numRNAPolymerases = 50;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = numRNAPolymerases;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            m.boundEnzymes(:) = 0;
            
            %put all RNA polymerases in non-specifically bound state
            r.states = repmat(r.nonSpecificallyBoundValue, numRNAPolymerases, 1);
            r.positionStrands = [1000*(1:numRNAPolymerases)' ones(numRNAPolymerases, 1)];
            t.boundTranscriptionUnits = repmat(t.nullTranscriptValue, numRNAPolymerases, 1);
            t.boundTranscriptProgress = repmat(t.nullTranscriptValue, numRNAPolymerases, 1);
            t.boundTranscriptChromosome = repmat(t.nullTranscriptValue, numRNAPolymerases, 1);
            
            m.bindProteinToChromosome(r.positionStrands, m.enzymeIndexs_rnaPolymerase);
            
            %set fold changes to 1
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            
            %force state transition
            r.stateExpectations = zeros(4, 1);
            r.stateExpectations(r.specificallyBoundIndex) = 1;
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.specificallyBoundIndex, :) = 1;
            
            %store initial state
            initial_states = r.states;
            
            %evolve
            m.evolveState();
            
            %all non-specifically bound RNA polymerases advance to
            %specifically bound state
            assertTrue(all(r.states(initial_states == r.nonSpecificallyBoundValue) == ...
                r.specificallyBoundValue));
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
            
            %did all RNAPolys bind to the chromosome
            assertEqual(numRNAPolymerases, nnz(c.complexBoundSites));
            assertEqual(0, m.enzymes(m.enzymeIndexs_rnaPolymerase));
            assertEqual(0, m.boundEnzymes(m.enzymeIndexs_rnaPolymerase));
            assertEqual(numRNAPolymerases, m.boundEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme));
            assertEqual(0, m.boundEnzymes(m.enzymeIndexs_sigmaFactor));
        end
        
        function testRNAPolymeraseStateExpectations_specificallyBound(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            c.initialize();           
            
            %put all RNA polymerases in specifically bound state
            numRNAPolymerases = 50;
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = 0;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = numRNAPolymerases;
            m.boundEnzymes(:) = 0;
            
            iTU = find(t.transcriptionUnitDirections == 1, numRNAPolymerases, 'first');
            r.states = repmat(r.specificallyBoundValue, numRNAPolymerases, 1);
            r.positionStrands = [t.transcriptionUnitFivePrimeCoordinates(iTU) ones(numRNAPolymerases, 1)];
            t.boundTranscriptionUnits = iTU;
            t.boundTranscriptProgress = repmat(t.nullTranscriptValue, numRNAPolymerases, 1);
            t.boundTranscriptChromosome = ones(numRNAPolymerases, 1);
            
            m.bindProteinToChromosome(r.positionStrands, m.enzymeIndexs_rnaPolymeraseHoloenzyme);
            
            %set fold changes to 1
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            
            %force state transition
            r.stateExpectations = zeros(4, 1);
            r.stateExpectations(r.specificallyBoundIndex) = 1;
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.specificallyBoundIndex, :) = 1;
            
            %store initial state
            initial_states = r.states;
            
            %evolve
            m.evolveState();
            
            %all specifically bound RNA polymerases remain in specifically
            %bound state
            assertTrue(all(r.states(initial_states == r.specificallyBoundValue) == ...
                r.specificallyBoundValue));
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
        end
        
        function testTranscriptionUnitPolyASequences(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            
            %make all sequences poly A's
            for i = 1:numel(t.transcriptionUnitSequences)
                t.transcriptionUnitSequences{i} = ...
                    repmat('A',1,length(t.transcriptionUnitSequences{i}));
            end
            
            %store initial state for subsequent comparison
            initial_substrates = m.substrates;
            
            %evolve
            m.evolveState();
            
            %check NTPs other than ATP not used
            assertEqual(initial_substrates(m.substrateIndexs_ntp(2:4)), m.substrates(m.substrateIndexs_ntp(2:4)));
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
        end
        
        function testTranscriptionUnitBindingProbabilities(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            
            %number of free RNA polymerases
            numPolymerase = 50;
            
            %saturating substrate and enzyme, except RNA polymerase
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = numPolymerase;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            m.boundEnzymes(:) = 0;
            
            %set Binding Probabilities
            r.transcriptionFactorBindingProbFoldChange(:)=1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            m.transcriptionUnitBindingProbabilities(1:numPolymerase) = m.randStream.rand(numPolymerase, 1);
            m.transcriptionUnitBindingProbabilities(numPolymerase+1:end) = 0;
            m.transcriptionUnitBindingProbabilities = ...
                m.transcriptionUnitBindingProbabilities / ...
                sum(m.transcriptionUnitBindingProbabilities);
            
            %no bound RNA polymeases
            c.initialize();
            m.boundEnzymes(:) = 0;
            r.states = repmat(r.freeValue, numPolymerase, 1);
            r.positionStrands = zeros(numPolymerase, 2);
            t.boundTranscriptionUnits = repmat(t.nullTranscriptValue, numPolymerase, 1);
            t.boundTranscriptProgress = repmat(t.nullTranscriptValue, numPolymerase, 1);
            t.boundTranscriptChromosome = repmat(t.nullTranscriptValue, numPolymerase, 1);
            m.chromosome.initialize();
            
            %set fold changes to 1
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            
            %force RNA polymerase transitions into specifically bound state
            r.stateExpectations = zeros(4, 1);
            r.stateExpectations(r.specificallyBoundIndex) = 1;
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.specificallyBoundIndex, :) = 1;
            
            %evolve
            m.evolveState();
            
            %each transcription unit with non-zero binding probability
            %should be bound by a RNA polymerase, and all transcription
            %units with zero binding probabilities should not be bound by
            %a RNA polymerase
            assertElementsAlmostEqual(length(find(r.states == r.specificallyBoundValue)),...
                numPolymerase, 'absolute', 5);
            assertElementsAlmostEqual(length(unique(t.boundTranscriptionUnits)'),...
                numPolymerase, 'absolute', 5);
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
        end
        
        function testTranscriptionUnitRegulationFoldChange(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            
            %number of free RNA polymerases
            numPolymerase = 50;
            
            %saturating substrate and enzyme other than RNA polymerase
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = numPolymerase;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme)=0;
            m.boundEnzymes(:) = 0;
            
            %set Binding Probabilities
            r.transcriptionFactorBindingProbFoldChange(1:numPolymerase) = m.randStream.rand(numPolymerase, 1);
            r.transcriptionFactorBindingProbFoldChange(numPolymerase+1:end) = 0;
            r.supercoilingBindingProbFoldChange(:) = 1;
            m.transcriptionUnitBindingProbabilities = m.randStream.rand(size(m.transcriptionUnitBindingProbabilities));
            m.transcriptionUnitBindingProbabilities = ...
                m.transcriptionUnitBindingProbabilities / ...
                sum(m.transcriptionUnitBindingProbabilities);
            
            %no bound RNA polymeases
            c.initialize();
            m.boundEnzymes(:)=0;
            r.states = repmat(r.freeValue,numPolymerase,1);
            r.positionStrands = zeros(numPolymerase, 2);
            t.boundTranscriptProgress = repmat(t.nullTranscriptValue, numPolymerase, 1);
            t.boundTranscriptChromosome = repmat(t.nullTranscriptValue, numPolymerase, 1);
            t.boundTranscriptionUnits = repmat(t.nullTranscriptValue, numPolymerase, 1);
            m.chromosome.initialize();
            
            %force RNA polymerase transitions into specifically bound state
            r.stateExpectations = zeros(4, 1);
            r.stateExpectations(r.specificallyBoundIndex) = 1;
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.specificallyBoundIndex, :) = 1;
            
            %evolve
            m.evolveState();
            
            %each transcription unit with non-zero binding probability
            %should be bound by a RNA polymerase, and all transcription
            %units with zero binding probabilities should not be bound by
            %a RNA polymerase
            assertElementsAlmostEqual(length(find(r.states == r.specificallyBoundValue)),...
                numPolymerase, 'absolute', 5);
            assertElementsAlmostEqual(length(unique(t.boundTranscriptionUnits)'),...
                numPolymerase, 'absolute', 5);
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
        end
        
        function testAccessibleTranscriptionUnits(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            
            %number of free RNA polymerases
            numPolymerase = 50;
            
            %saturating substrate and enzyme other than RNA polymerase
            m.substrates(:)=1e6;
            m.enzymes(:)=1e6;
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = numPolymerase;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            m.boundEnzymes(:) = 0;
            
            %set Binding Probabilities
            c.initialize();
            c.polymerizedRegions(1, 1:2) = c.transcriptionUnitStartCoordinates(numPolymerase+1)-1;
            r.transcriptionFactorBindingProbFoldChange = m.randStream.rand(size(r.transcriptionFactorBindingProbFoldChange));
            r.supercoilingBindingProbFoldChange(:) = 1;
            m.transcriptionUnitBindingProbabilities = m.randStream.rand(size(m.transcriptionUnitBindingProbabilities));
            m.transcriptionUnitBindingProbabilities = ...
                m.transcriptionUnitBindingProbabilities / ...
                sum(m.transcriptionUnitBindingProbabilities);
            
            %no bound RNA polymeases
            m.boundEnzymes(:)=0;
            r.states = repmat(r.freeValue, numPolymerase, 1);
            r.positionStrands = zeros(numPolymerase, 2);
            t.boundTranscriptProgress = repmat(t.nullTranscriptValue,numPolymerase,1);
            t.boundTranscriptChromosome = repmat(t.nullTranscriptValue,numPolymerase,1);
            t.boundTranscriptionUnits = repmat(t.nullTranscriptValue,numPolymerase,1);
            
            %set fold changes to 1
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            
            %force RNA polymerase transitions into specifically bound state
            r.stateExpectations = zeros(4, 1);
            r.stateExpectations(r.specificallyBoundIndex) = 1;
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.specificallyBoundIndex, :) = 1;
                        
            %evolve
            m.evolveState();
            
            %each transcription unit with non-zero binding probability
            %should be bound by a RNA polymerase, and all transcription
            %units with zero binding probabilities should not be bound by
            %a RNA polymerase
            assertElementsAlmostEqual(sum(r.states == r.specificallyBoundValue), length(r.states), 'absolute', 5);
            assertEqual(length(unique(t.boundTranscriptionUnits(r.states == r.specificallyBoundValue))'),...
                sum(r.states == r.specificallyBoundValue), 'Not all transcription units with non-zero probability were bound by polymerases');
            
            %transcript and rnaPolymerase states in sync
            assertEqual(t.boundTranscriptProgress, (r.states+abs(r.states))/2);
            assertEqual(length(r.states), length(t.boundTranscriptionUnits));
            assertEqual(length(r.states), length(t.boundTranscriptChromosome));
        end
        
        function testSpecificallyBoundPolysUnbinding(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            c.initialize();
            
            %put all RNA polymerases in specifically bound state
            nPols = sum(r.states ~= r.notExistValue);
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = 0;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = nPols;
            m.enzymes(m.enzymeIndexs_transcriptionFactors(1)) = 0;
            m.boundEnzymes(:) = 0;
            
            iTU = find(t.transcriptionUnitDirections == 1, nPols, 'first');
            iPol = (r.states ~= r.notExistValue);
            r.states(iPol) = r.specificallyBoundValue;
            r.positionStrands(iPol, :) = [...
                t.transcriptionUnitFivePrimeCoordinates(iTU) - (2 * t.transcriptionUnitDirections(iTU) - 1) ...
                ones(size(iTU))];
            t.boundTranscriptionUnits(iPol) = iTU;
            t.boundTranscriptProgress(iPol) = t.nullTranscriptValue;
            t.boundTranscriptChromosome(iPol) = 1;
            
            m.bindProteinToChromosome(r.positionStrands(iPol, :), m.enzymeIndexs_rnaPolymeraseHoloenzyme);
            
            %force state transition
            r.stateExpectations = zeros(4, 1);
            r.stateExpectations(r.specificallyBoundIndex) = 1;
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.specificallyBoundIndex, :) = 1; 
            
            %force binding
            m.evolveState();
            initialPositionStands = r.positionStrands(r.states == r.specificallyBoundValue);
            assertTrue(~isempty(initialPositionStands));
            
            %force state transition
            r.stateExpectations = zeros(4, 1);
            r.stateExpectations(r.nonSpecificallyBoundIndex) = 1;
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.nonSpecificallyBoundIndex, :) = 1; 
            
            %now force unbinding
            m.evolveState();
            assertTrue(~any(r.states == r.specificallyBoundValue));
        end
        
        function testChromosomeBinding(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            
            rPolyGblIdx = [
                m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymerase);
                m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
                ];
            totEnzymes = m.enzymes + m.boundEnzymes;            
            totEnzymes(m.enzymeIndexs_rnaPolymerase) = ...
                totEnzymes(m.enzymeIndexs_rnaPolymerase) + ...
                totEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
            totEnzymes(m.enzymeIndexs_sigmaFactor) = ...
                totEnzymes(m.enzymeIndexs_sigmaFactor) + ...
                totEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
            totEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            
            for seed = round(mod(now, 1) * 1e7)
                %fprintf('seed: %d\n', seed);
                
                %% set seed
                m.seed = seed;
                m.seedRandStream();
                for i = 1:numel(m.states)
                    s = m.states{i};
                    s.seed = seed;
                    s.seedRandStream();
                end
                
                %% initialize to random state with saturating substrates and sufficient enzymes
                %lots of substrates
                m.substrates(:) = 1e6;
                
                %clear bound RNA polymerases
                m.enzymes = totEnzymes;
                m.boundEnzymes(:) = 0;
                c.complexBoundSites(find(c.complexBoundSites == rPolyGblIdx(1))) = 0; %#ok<FNDSB>
                c.complexBoundSites(find(c.complexBoundSites == rPolyGblIdx(2))) = 0; %#ok<FNDSB>
                
                r.states(:) = r.notExistValue;
                r.positionStrands(:) = 0;
                t.boundTranscriptProgress(:) = t.nullTranscriptValue;
                t.boundTranscriptionUnits(:) = t.nullTranscriptValue;
                t.boundTranscriptChromosome(:) = t.nullTranscriptValue;
                
                %initialize RNA polymerases and transcripts to random state
                m.initializeState();
                
                %% assert that chromosome, transcript, and RNA polymerase
                %states remain synchronized
                for i = 1:1000
                    try                        
                        eIdxs = setdiff(1:numel(m.enzymes), [m.enzymeIndexs_sigmaFactor; m.enzymeIndexs_rnaPolymerase; m.enzymeIndexs_rnaPolymeraseHoloenzyme]);
                        assertEqual(totEnzymes(eIdxs), m.enzymes(eIdxs) + m.boundEnzymes(eIdxs));
                        
                        eIdxs = [m.enzymeIndexs_rnaPolymerase; m.enzymeIndexs_rnaPolymeraseHoloenzyme];
                        assertEqual(sum(totEnzymes(eIdxs)), sum(m.enzymes(eIdxs) + m.boundEnzymes(eIdxs)));
                        
                        eIdxs = [m.enzymeIndexs_sigmaFactor; m.enzymeIndexs_rnaPolymeraseHoloenzyme];
                        assertEqual(sum(totEnzymes(eIdxs)), sum(m.enzymes(eIdxs) + m.boundEnzymes(eIdxs)));
                        
                        this.isStateValid();
                    catch exception
                        %detailed status message
                        table = [];
                        table = [table sprintf('\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\n', ...
                            'State', 'State', 'Pos', 'Strnd', 'Prg', 'TU', 'Chr')]; %#ok<*AGROW>
                        table = [table sprintf('\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
                            repmat('=',1,6), repmat('=',1,6), repmat('=',1,6), repmat('=',1,6), repmat('=',1,6), repmat('=',1,6), repmat('=',1,6))];
                        [~, order] = sortrows([min(1, r.states) r.positionStrands], [-1 2 3]);
                        for k = 1:size(r.states, 1)
                            j = order(k);
                            switch r.states(j)
                                case r.notExistValue, stateStr = 'Not Ex';
                                case r.freeValue, stateStr = 'Free';
                                case r.nonSpecificallyBoundValue, stateStr = 'Non-SB';
                                case r.specificallyBoundValue, stateStr = 'SB';
                                otherwise, stateStr = 'Active';
                            end
                            table = [table sprintf('\t%6s\t%6d\t%6d\t%6d\t%6d\t%6d\t%6d\n', ...
                                stateStr, r.states(j), r.positionStrands(j, 1), r.positionStrands(j, 2), ...
                                t.boundTranscriptProgress(j), t.boundTranscriptionUnits(j), t.boundTranscriptChromosome(j))];
                        end
                        throw(MException('Transcription:error', ...
                            'State out of sync\n\tseed: %d\n\titer: %d\n\n%s\n\n%s\n\n\tseed: %d\n\titer: %d', ...
                            seed, i, table, exception.getReport(), seed, i));
                    end
                    
                    m.evolveState();
                end
            end
        end
        
        function testNewPolymerase(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            
            %more RNA polymerase
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = ...
                m.enzymes(m.enzymeIndexs_rnaPolymerase) + 200;
            
            %store current state
            enzymes = m.enzymes;
            boundEnzymes = m.boundEnzymes;
            states = r.states;
            
            %evolve
            m.evolveState();
            
            %test
            assertEqual(enzymes(m.enzymeIndexs_rnaPolymerase) + boundEnzymes(m.enzymeIndexs_rnaPolymerase) + boundEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) - 200, sum(states ~= r.notExistValue));
                        
            eIdxs = setdiff(1:numel(m.enzymes), [m.enzymeIndexs_sigmaFactor; m.enzymeIndexs_rnaPolymerase; m.enzymeIndexs_rnaPolymeraseHoloenzyme]);
            assertEqual(enzymes(eIdxs) + boundEnzymes(eIdxs), m.enzymes(eIdxs) + m.boundEnzymes(eIdxs));
            
            eIdxs = [m.enzymeIndexs_rnaPolymerase; m.enzymeIndexs_rnaPolymeraseHoloenzyme];
            assertEqual(sum(enzymes(eIdxs) + boundEnzymes(eIdxs)), sum(m.enzymes(eIdxs) + m.boundEnzymes(eIdxs)));
            
            eIdxs = [m.enzymeIndexs_sigmaFactor; m.enzymeIndexs_rnaPolymeraseHoloenzyme];
            assertEqual(sum(enzymes(eIdxs) + boundEnzymes(eIdxs)), sum(m.enzymes(eIdxs) + m.boundEnzymes(eIdxs)));
                        
            assertIn(m.enzymes(m.enzymeIndexs_rnaPolymerase), [200 Inf]);
            
            assertAllEqual(r.freeValue, r.states(find(states ~= r.notExistValue, 1, 'last')+(1:200)));
            assertAllEqual(0, r.positionStrands(find(states ~= r.notExistValue, 1, 'last')+(1:200), :));
            assertAllEqual(t.nullTranscriptValue, t.boundTranscriptProgress(find(states ~= r.notExistValue, 1, 'last')+(1:200)));
            assertAllEqual(t.nullTranscriptValue, t.boundTranscriptionUnits(find(states ~= r.notExistValue, 1, 'last')+(1:200)));
            assertAllEqual(t.nullTranscriptValue, t.boundTranscriptChromosome(find(states ~= r.notExistValue, 1, 'last')+(1:200)));
        end
        
        function testCodirectionalPolymerases_OnSameTU(this)
            m = this.process;
            c = m.chromosome;
            r = m.rnaPolymerases;
            t = m.transcripts;
            
            %free all RNA polymerases
            m.enzymes = m.enzymes + m.boundEnzymes;
            m.boundEnzymes(:) = 0;
            r.states(r.states ~= r.notExistValue) = r.freeValue;
            r.positionStrands(:) = 0;
            t.boundTranscriptionUnits(:) = t.nullTranscriptValue;
            t.boundTranscriptProgress(:) = t.nullTranscriptValue;
            t.boundTranscriptChromosome(:) = t.nullTranscriptValue;
            
            c.initialize();
            
            %setup multiple RNA polymerase bound near each other on first
            %transcription unit
            iTU = find(t.transcriptionUnitDirections == 1, 1, 'first');
            
            m.substrates(:) = 1e6;
            r.states(1:4) = 2 + ((1:4)-1) * m.enzymeDNAFootprints(m.enzymeIndexs_rnaPolymerase);
            r.positionStrands(1:4, 1) = t.transcriptionUnitFivePrimeCoordinates(iTU) - 1 +  r.states(1:4);
            r.positionStrands(1:4, 2) = 1;
            t.boundTranscriptionUnits(1:4) = iTU;
            t.boundTranscriptProgress(1:4) = r.states(1:4);
            t.boundTranscriptChromosome(1:4) = 1;
            
            m.bindProteinToChromosome(r.positionStrands(1:4, :), m.enzymeIndexs_rnaPolymerase);
            
            states = r.states;
            positionStrands = r.positionStrands;
            boundTranscriptionUnits = t.boundTranscriptionUnits;
            boundTranscriptProgress = t.boundTranscriptProgress;
            boundTranscriptChromosome = t.boundTranscriptChromosome;
            
            %evolve
            m.evolveState();
            
            %tests
            assertEqual(states(1:3), r.states(1:3));
            assertEqual(positionStrands(1:3, :), r.positionStrands(1:3, :));
            assertEqual(boundTranscriptionUnits(1:3), t.boundTranscriptionUnits(1:3));
            assertEqual(boundTranscriptProgress(1:3), t.boundTranscriptProgress(1:3));
            assertEqual(boundTranscriptChromosome(1:3), t.boundTranscriptChromosome(1:3));
            
            assertEqual(states(4) + m.rnaPolymeraseElongationRate, r.states(4));
            assertEqual(positionStrands(4, 1) + m.rnaPolymeraseElongationRate, r.positionStrands(4, 1));
            assertEqual(positionStrands(4, 2), r.positionStrands(4, 2));
            assertEqual(boundTranscriptionUnits(4), t.boundTranscriptionUnits(4));
            assertEqual(boundTranscriptProgress(4) + m.rnaPolymeraseElongationRate, t.boundTranscriptProgress(4));
            assertEqual(boundTranscriptChromosome(4), t.boundTranscriptChromosome(4));
        end
        
        function testCodirectionalPolymerases_OnDifferentPosStrndTUs(this)
            m = this.process;
            c = m.chromosome;
            r = m.rnaPolymerases;
            t = m.transcripts;
            
            %setup 1 active codirectional RNA polymerases, near end of
            %transcript and near start of next transcription unit, with
            %sufficient NTPs
            m.substrates(:) = 1e6;
            
            nPolymerases = 2;
            m.enzymes = m.enzymes + m.boundEnzymes;
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = nPolymerases;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            m.boundEnzymes(:) = 0;
            m.RNAs(:) = 0;
            
            c.initialize();
            
            tuLens = c.transcriptionUnitLengths;
            tuStarts = c.transcriptionUnitStartCoordinates;
            tuEnds = tuStarts + tuLens;
            tuDirs = t.transcriptionUnitDirections;
            
            iTU = find(...
                tuDirs(1:end-1) == 1 & tuDirs(2:end) == 1 & ...
                tuEnds(1:end-1) > tuStarts(2:end), ...
                1, 'first');
            
            r.states = [tuLens(iTU); r.freeValue];
            r.positionStrands = ...
                [t.transcriptionUnitFivePrimeCoordinates(iTU) + tuLens(iTU) - 1 1; 0 0];
            t.boundTranscriptionUnits = [iTU; t.nullTranscriptValue];
            t.boundTranscriptProgress = [r.states(1); t.nullTranscriptValue];
            t.boundTranscriptChromosome = [1; t.nullTranscriptValue];
            
            m.bindProteinToChromosome(r.positionStrands(1, :), m.enzymeIndexs_rnaPolymerase);
            
            m.transcriptionUnitBindingProbabilities(:) = 0;
            m.transcriptionUnitBindingProbabilities(iTU + 1) = 1;
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            
            %assert 1 RNA produced
            m.evolveState();
            assertEqual([1; 0], m.RNAs(iTU + (0:1)));
            assertTrue(1 <= sum(r.states == r.freeValue));
            assertTrue(0 == sum(r.states == r.specificallyBoundValue));
            assertTrue(0 == sum(r.states >= r.activelyTranscribingValue));
            
            %assert second RNA produced
            for i = 1:ceil(c.transcriptionUnitLengths(iTU+1) / m.rnaPolymeraseElongationRate) + 100
                if m.RNAs(iTU+1)
                    break;
                end
                m.evolveState();
            end
            assertEqual(1, m.RNAs(iTU));
            assertIn(m.RNAs(iTU+1), [1 Inf]);
        end
        
        function testCodirectionalPolymerases_OnDifferentNegStrndTUs(this)
            m = this.process;
            c = m.chromosome;
            r = m.rnaPolymerases;
            t = m.transcripts;
            
            %setup 1 active codirectional RNA polymerases, near end of
            %transcript and near start of next transcription unit, with
            %sufficient NTPs
            m.substrates(:) = 1e6;
            
            nPolymerases = 2;
            m.enzymes = m.enzymes + m.boundEnzymes;
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = nPolymerases;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            m.boundEnzymes(:) = 0;
            m.RNAs(:) = 0;
            
            c.initialize();
            
            tuLens = c.transcriptionUnitLengths;
            tuStarts = c.transcriptionUnitStartCoordinates;
            tuEnds = tuStarts + tuLens;
            tuDirs = t.transcriptionUnitDirections;
            
            iTU = find(...
                tuDirs(1:end-1) == 0 & tuDirs(2:end) == 0 & ...
                tuEnds(1:end-1) > tuStarts(2:end), ...
                1, 'first');
            
            r.states = [tuLens(iTU+1); r.freeValue];
            r.positionStrands = ...
                [t.transcriptionUnitFivePrimeCoordinates(iTU+1) - (tuLens(iTU+1) - 1) 2; 0 0];
            t.boundTranscriptionUnits = [iTU+1; t.nullTranscriptValue];
            t.boundTranscriptProgress = [r.states(1); t.nullTranscriptValue];
            t.boundTranscriptChromosome = [1; t.nullTranscriptValue];
            
            m.bindProteinToChromosome(r.positionStrands(1, :), m.enzymeIndexs_rnaPolymerase);
            
            m.transcriptionUnitBindingProbabilities(:) = 0;
            m.transcriptionUnitBindingProbabilities(iTU) = 1;
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            
            %assert 1 RNA produced
            m.evolveState();
            assertEqual([0; 1], m.RNAs(iTU + (0:1)));
            assertTrue(1 <= sum(r.states == r.freeValue));
            assertTrue(0 == sum(r.states == r.specificallyBoundValue));
            assertTrue(0 == sum(r.states >= r.activelyTranscribingValue));
            
            %assert second RNA produced
            for i = 1:ceil(c.transcriptionUnitLengths(iTU) / m.rnaPolymeraseElongationRate) + 100
                if m.RNAs(iTU)
                    break;
                end
                m.evolveState();
            end
            assertEqual(1, m.RNAs(iTU+1));
            assertIn(m.RNAs(iTU), [1 Inf]);
        end
        
        function testAntiDirectionalPolymerases(this)
            m = this.process;
            c = m.chromosome;
            r = m.rnaPolymerases;
            t = m.transcripts;
            
            %setup two active antidirectional RNA polymerases, both near
            %end of transcript, with sufficient NTPs
            m.substrates(:) = 1e6;
            
            nPolymerases = 2;
            m.enzymes = m.enzymes + m.boundEnzymes;
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = nPolymerases;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            m.boundEnzymes(:) = 0;
            m.RNAs(:) = 0;
            
            c.initialize();
            
            tuLens = c.transcriptionUnitLengths;
            tuStarts = c.transcriptionUnitStartCoordinates;
            tuEnds = tuStarts + tuLens;
            tuDirs = t.transcriptionUnitDirections;
            
            iTU = find(...
                tuDirs(1:end-1) == 1 & tuDirs(2:end) == 0 & ...
                tuEnds(1:end-1) > tuStarts(2:end), ...
                1, 'first');
            
            r.states = tuLens(iTU+(0:1)) - m.enzymeDNAFootprints(m.enzymeIndexs_rnaPolymerase);
            r.positionStrands = ...
                [t.transcriptionUnitFivePrimeCoordinates(iTU+(0:1)) + (2 * tuDirs(iTU+(0:1)) -1) .* (r.states - 1) ...
                2 - tuDirs(iTU+(0:1))];
            t.boundTranscriptionUnits = [iTU; iTU + 1];
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptChromosome = ones(nPolymerases, 1);
            
            m.bindProteinToChromosome(r.positionStrands, m.enzymeIndexs_rnaPolymerase);
            
            %simulate
            for i = 1:5
                m.evolveState();
            end
            
            %assert both RNAs produced
            assertAllEqual(1, m.RNAs(iTU + (0:1)));
        end
        
        function testInitializationConsistency(this)
            %process
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            rPolyGblIdx = [
                m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymerase);
                m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
                ];            
            
            assertEqual(sum(sum([
                r.states == r.nonSpecificallyBoundValue,...
                r.states >  r.activelyTranscribingValue])), nnz(c.complexBoundSites == rPolyGblIdx(1)));
            assertEqual(sum(sum([
                r.states == r.specificallyBoundValue,...
                r.states == r.activelyTranscribingValue])), nnz(c.complexBoundSites == rPolyGblIdx(2)));
            assertEqual(r.states >= r.activelyTranscribingValue, t.boundTranscriptProgress >= r.activelyTranscribingValue)
        end
        
        function testFairBinding(this)
            m = this.process;
            c = m.chromosome;
            t = m.transcripts;
            
            probs = m.computeRNAPolymeraseTUBindingProbabilities();
            probs = probs(:, 1);
            tuDirs = 2 * t.transcriptionUnitDirections - 1;
            tuStrnds = c.transcriptionUnitStrands;
            tu5Coords = t.transcriptionUnitFivePrimeCoordinates;
            
            %randsample
            rnas = zeros(size(m.RNAs));
            for i = 1:1000
                iTU = m.randStream.randsample(numel(tuDirs), 1, false, probs);
                rnas(iTU) = rnas(iTU) + 1;
            end
            
            assertIn(corr(probs, rnas), [0.95 1]);
            assertIn(180 / pi * acos((probs' * rnas) / (sqrt(probs' * probs) * sqrt(rnas' * rnas))), [0 15]);
            
            %bindProteinToChromosome
            rnas = zeros(size(m.RNAs));
            for i = 1:1000
                c.initialize();
                m.enzymes(:) = 1;
                m.boundEnzymes(:) = 0;
                [~, iTU] = m.bindProteinToChromosome( ...
                    [tu5Coords-tuDirs tuStrnds], ...
                    m.enzymeIndexs_rnaPolymerase, 1, probs, ...
                    true, true, 1, false, []);
                rnas(iTU) = rnas(iTU) + 1;
            end
            
            assertIn(corr(probs, rnas), [0.95 1]);
            assertIn(180 / pi * acos((probs' * rnas) / (sqrt(probs' * probs) * sqrt(rnas' * rnas))), [0 20]);
        end
        
        function testExpectations(this)
            import edu.stanford.covert.cell.kb.ssRNA;
            
            %references
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            
            %initialize at beginning of cell cycle
            nPols = ...
                + m.enzymes(m.enzymeIndexs_rnaPolymerase) ...
                + m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) ...
                + m.boundEnzymes(m.enzymeIndexs_rnaPolymerase) ...
                + m.boundEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
            probs = m.computeRNAPolymeraseTUBindingProbabilities();
            probs = probs(:, 1);
            
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_ntp) = 0;
            baseCounts = m.rna.baseCounts(m.rna.nascentIndexs, m.metabolite.nmpIndexs);
            ntpProd = baseCounts' * probs;
            ntpProd = 134.0758 * ntpProd / sum(ntpProd);
            
            m.enzymes = m.enzymes + m.boundEnzymes;
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = ...
                m.enzymes(m.enzymeIndexs_rnaPolymerase) + ...
                m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
            m.enzymes(m.enzymeIndexs_sigmaFactor) = ...
                m.enzymes(m.enzymeIndexs_sigmaFactor) + ...
                m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            m.boundEnzymes(:) = 0;
            m.RNAs(:) = 0;
            totEnzymes = m.enzymes;
            
            c.initialize();
            
            r.states(1:nPols) = r.freeValue;
            r.states(nPols + 1:end) = r.notExistValue;
            r.positionStrands(:) = 0;
            t.boundTranscriptProgress(:) = t.nullTranscriptValue;
            t.boundTranscriptionUnits(:) = t.nullTranscriptValue;
            t.boundTranscriptChromosome(:) = t.nullTranscriptValue;
            m.initializeState();
            assertAllEqual(0, m.RNAs);
            
            iPol = find(r.states >= r.activelyTranscribingValue);
            initTranscriptBaseCounts = zeros(4, 1);
            for i = 1:numel(iPol)
                seq = t.transcriptionUnitSequences{t.boundTranscriptionUnits(iPol(i))}(1:r.states(iPol(i)) - 1);
                initTranscriptBaseCounts = initTranscriptBaseCounts + ssRNA.computeBaseCount(seq, 4, 1:4)';
            end
            
            assertAllEqual(1, sum(m.stateTransitionProbabilities, 1));
            assertEqual(1, sum(r.stateExpectations));
            assertIn(sum(ntpProd) / m.rnaPolymeraseElongationRate / nPols, [0 r.stateExpectations(r.activelyTranscribingIndex)]);
            assertEqual(m.stateTransitionProbabilities(:, r.nonSpecificallyBoundIndex), ...
                m.stateTransitionProbabilities(:, r.freeIndex));
            
            %evolve
            this.isStateValid();
            iterMax = 2000;
            cumNtpProd = zeros(size(ntpProd));
            occ = zeros(4, iterMax);
            for i = 1:iterMax
                %if mod(i, 1000) == 1
                %    fprintf('%d\n',i);
                %end
                tmpNtpProd = m.randStream.stochasticRound(ntpProd);
                cumNtpProd = cumNtpProd + tmpNtpProd;
                m.substrates(m.substrateIndexs_ntp) = ...
                    m.substrates(m.substrateIndexs_ntp) + ...
                    tmpNtpProd;
                
                m.evolveState();
                
                eIdxs = setdiff(1:numel(m.enzymes), [m.enzymeIndexs_sigmaFactor; m.enzymeIndexs_rnaPolymerase; m.enzymeIndexs_rnaPolymeraseHoloenzyme]);
                assertEqual(totEnzymes(eIdxs), m.enzymes(eIdxs) + m.boundEnzymes(eIdxs));
                
                eIdxs = [m.enzymeIndexs_rnaPolymerase; m.enzymeIndexs_rnaPolymeraseHoloenzyme];
                assertEqual(sum(totEnzymes(eIdxs)), sum(m.enzymes(eIdxs) + m.boundEnzymes(eIdxs)));
                
                eIdxs = [m.enzymeIndexs_sigmaFactor; m.enzymeIndexs_rnaPolymeraseHoloenzyme];
                assertEqual(sum(totEnzymes(eIdxs)), sum(m.enzymes(eIdxs) + m.boundEnzymes(eIdxs)));
                
                this.isStateValid();
                
                occ(:, i) = r.stateOccupancies;
            end
            
            %assert binding probabilities and synthesis are 1:1 correlated
            newProbs = m.computeRNAPolymeraseTUBindingProbabilities();
            newProbs = newProbs(:, 1);
            assertEqual(probs(setdiff(1:end, m.transcriptionUnitIndexs_DnaAR12345Boxes)), newProbs(setdiff(1:end, m.transcriptionUnitIndexs_DnaAR12345Boxes)))
            assertIn(corr(probs, m.RNAs), [0.75 1]);
            assertIn(180 / pi * acos((probs' * m.RNAs) / (sqrt(probs' * probs) * sqrt(m.RNAs' * m.RNAs))), [0 40]);
            
            %assert NTPs used
            assertEqual(0, min(m.substrates(m.substrateIndexs_ntp)));
            assertIn((m.substrates(m.substrateIndexs_ntp)' * m.substrateMolecularWeights(m.substrateIndexs_ntp)) / ...
                (m.RNAs' * m.rna.molecularWeights(m.rna.nascentIndexs)), [0 0.10]);
            
            %assert NTPs correctly accounted
            iPol = find(r.states >= r.activelyTranscribingValue);
            transcriptBaseCounts = zeros(4, 1);
            for i = 1:numel(iPol)
                seq = t.transcriptionUnitSequences{t.boundTranscriptionUnits(iPol(i))}(1:r.states(iPol(i)) - 1);
                transcriptBaseCounts = transcriptBaseCounts + ssRNA.computeBaseCount(seq, 4, 1:4)';
            end
            assertEqual(cumNtpProd - m.substrates(m.substrateIndexs_ntp), ...
                m.rna.baseCounts(m.rna.nascentIndexs, m.metabolite.nmpIndexs)' * m.RNAs + ...
                transcriptBaseCounts - initTranscriptBaseCounts);
        end
    end
       
    %test state transitioning
    methods
        function testStateTransitions_Free(this)
            m = this.process;
            r = m.rnaPolymerases;
            
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.freeIndex, :) = 1;
            r.stateExpectations(:) = 0;
            r.stateExpectations(r.freeIndex) = 1;
            
            occ = this.calcStateTransitions(500, true);
            occ = sum(occ, 2) / sum(occ(:));
            assertEqual(r.stateExpectations, occ);
        end
        
        function testStateTransitions_NSB(this)
            m = this.process;
            r = m.rnaPolymerases;
            
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.nonSpecificallyBoundIndex, :) = 1;
            r.stateExpectations(:) = 0;
            r.stateExpectations(r.nonSpecificallyBoundIndex) = 1;
            
            occ = this.calcStateTransitions(500, true);
            occ = sum(occ, 2) / sum(occ(:));
            assertEqual(r.stateExpectations, occ);
        end
        
        function testStateTransitions_SB(this)
            m = this.process;
            r = m.rnaPolymerases;
            
            m.stateTransitionProbabilities = eye(4, 4);
            m.stateTransitionProbabilities(:, r.freeIndex) = m.stateTransitionProbabilities(:, r.nonSpecificallyBoundIndex);
            
            r.stateExpectations(:) = 0;
            r.stateExpectations(r.specificallyBoundIndex) = 1;
            
            occ = this.calcStateTransitions(500, true);
            occ = sum(occ, 2) / sum(occ(:));
            assertElementsAlmostEqual(r.stateExpectations, occ, 'absolute', 0.50);
        end
        
        function testStateTransitions_NSB_Free(this)
            m = this.process;
            r = m.rnaPolymerases;
            
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.activelyTranscribingIndex, r.activelyTranscribingIndex) = 1;
            m.stateTransitionProbabilities(r.specificallyBoundIndex, r.specificallyBoundIndex) = 1;
            m.stateTransitionProbabilities(r.nonSpecificallyBoundIndex, r.nonSpecificallyBoundIndex) = .3;
            m.stateTransitionProbabilities(r.freeIndex, r.nonSpecificallyBoundIndex) = 1 - m.stateTransitionProbabilities(r.nonSpecificallyBoundIndex, r.nonSpecificallyBoundIndex);
            m.stateTransitionProbabilities(r.nonSpecificallyBoundIndex, r.freeIndex) = .3;
            m.stateTransitionProbabilities(r.freeIndex, r.freeIndex) = 1 - m.stateTransitionProbabilities(r.nonSpecificallyBoundIndex, r.freeIndex);
            
            r.stateExpectations(:) = 0;
            r.stateExpectations(r.nonSpecificallyBoundIndex) = 0.3;
            r.stateExpectations(r.freeIndex) = 0.7;
            
            occ = this.calcStateTransitions(500, true);
            occ = sum(occ, 2) / sum(occ(:));
            assertElementsAlmostEqual(r.stateExpectations, occ, 'relative', 5e-2);
        end
        
        function testStateTransitions_Free_NSB(this)
            m = this.process;
            r = m.rnaPolymerases;
            
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.activelyTranscribingIndex, r.activelyTranscribingIndex) = 1;
            m.stateTransitionProbabilities(r.specificallyBoundIndex, r.specificallyBoundIndex) = 1;
            m.stateTransitionProbabilities(r.freeIndex, r.freeIndex) = .3;
            m.stateTransitionProbabilities(r.nonSpecificallyBoundIndex, r.freeIndex) = 1 - m.stateTransitionProbabilities(r.freeIndex, r.freeIndex);
            m.stateTransitionProbabilities(r.freeIndex, r.nonSpecificallyBoundIndex) = .3;
            m.stateTransitionProbabilities(r.nonSpecificallyBoundIndex, r.nonSpecificallyBoundIndex) = 1 - m.stateTransitionProbabilities(r.freeIndex, r.nonSpecificallyBoundIndex);
            
            r.stateExpectations(:) = 0;
            r.stateExpectations(r.nonSpecificallyBoundIndex) = 0.7;
            r.stateExpectations(r.freeIndex) = 0.3;
            
            occ = this.calcStateTransitions(500, true);
            occ = sum(occ, 2) / sum(occ(:));
            assertElementsAlmostEqual(r.stateExpectations, occ, 'relative', 5e-2);
        end
        
        function testStateTransitions_SB_Free(this)
            m = this.process;
            r = m.rnaPolymerases;
            
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.activelyTranscribingIndex, r.activelyTranscribingIndex) = 1;            
            m.stateTransitionProbabilities(r.specificallyBoundIndex, r.specificallyBoundIndex) = .4;
            m.stateTransitionProbabilities(r.freeIndex, r.specificallyBoundIndex) = 1 - m.stateTransitionProbabilities(r.specificallyBoundIndex, r.specificallyBoundIndex);
            m.stateTransitionProbabilities(r.specificallyBoundIndex, r.freeIndex) = (.3-.3*0.4)/0.7;
            m.stateTransitionProbabilities(r.freeIndex, r.freeIndex) = 1 - m.stateTransitionProbabilities(r.specificallyBoundIndex, r.freeIndex);
            m.stateTransitionProbabilities(:, r.nonSpecificallyBoundIndex) = m.stateTransitionProbabilities(:, r.freeIndex);
            
            r.stateExpectations(:) = 0;
            r.stateExpectations(r.specificallyBoundIndex) = 0.3;
            r.stateExpectations(r.freeIndex) = 0.7;
            
            occ = this.calcStateTransitions(500, true);
            occ = sum(occ, 2) / sum(occ(:));
            assertElementsAlmostEqual(r.stateExpectations, occ, 'relative', 5e-2);
        end
        
        function testStateTransitions_Free_SB(this)
            m = this.process;
            r = m.rnaPolymerases;
            
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.activelyTranscribingIndex, r.activelyTranscribingIndex) = 1;
            m.stateTransitionProbabilities(r.freeIndex, r.freeIndex) = .3;
            m.stateTransitionProbabilities(r.specificallyBoundIndex, r.freeIndex) = 1 - m.stateTransitionProbabilities(r.freeIndex, r.freeIndex);
            m.stateTransitionProbabilities(r.freeIndex, r.specificallyBoundIndex) = .3;
            m.stateTransitionProbabilities(r.specificallyBoundIndex, r.specificallyBoundIndex) = 1 - m.stateTransitionProbabilities(r.freeIndex, r.specificallyBoundIndex);
            m.stateTransitionProbabilities(:, r.nonSpecificallyBoundIndex) = m.stateTransitionProbabilities(:, r.freeIndex);
            
            r.stateExpectations(:) = 0;
            r.stateExpectations(r.specificallyBoundIndex) = 0.7;
            r.stateExpectations(r.freeIndex) = 0.3;
            
            occ = this.calcStateTransitions(500, true);
            occ = sum(occ, 2) / sum(occ(:));
            assertElementsAlmostEqual(r.stateExpectations, occ, 'relative', 35e-2);
        end
        
        function testStateTransitions_SB_NSB_Free(this)
            m = this.process;
            r = m.rnaPolymerases;
            
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.activelyTranscribingIndex, r.activelyTranscribingIndex) = 1;
            
            m.stateTransitionProbabilities(r.specificallyBoundIndex, r.specificallyBoundIndex) = .90;
            m.stateTransitionProbabilities(r.nonSpecificallyBoundIndex, r.specificallyBoundIndex) = (.55 - .8*.65)/.35;
            m.stateTransitionProbabilities(r.freeIndex, r.specificallyBoundIndex) = .1 - m.stateTransitionProbabilities(r.nonSpecificallyBoundIndex, r.specificallyBoundIndex);
            
            m.stateTransitionProbabilities(r.specificallyBoundIndex, r.nonSpecificallyBoundIndex) = (.35 - .9*.35) / (.65);
            m.stateTransitionProbabilities(r.nonSpecificallyBoundIndex, r.nonSpecificallyBoundIndex) = .8;
            m.stateTransitionProbabilities(r.freeIndex, r.nonSpecificallyBoundIndex) = .2- m.stateTransitionProbabilities(r.specificallyBoundIndex, r.nonSpecificallyBoundIndex);
            
            m.stateTransitionProbabilities(:, r.freeIndex) = m.stateTransitionProbabilities(:, r.nonSpecificallyBoundIndex);
            
            r.stateExpectations(:) = 0;
            r.stateExpectations(r.specificallyBoundIndex) = 0.35;
            r.stateExpectations(r.nonSpecificallyBoundIndex) = 0.55;
            r.stateExpectations(r.freeIndex) = 0.1;
            
            occ = this.calcStateTransitions(500, true);
            occ = sum(occ, 2) / sum(occ(:));
            assertElementsAlmostEqual(r.stateExpectations, occ, 'relative', 5e-2);
        end
        
        function testStateTransitions_AT_SB_Free(this)
            m = this.process;
            r = m.rnaPolymerases;
            rna = m.rna;
            
            nPols = ...
                + m.enzymes(m.enzymeIndexs_rnaPolymerase) ...
                + m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) ...
                + m.boundEnzymes(m.enzymeIndexs_rnaPolymerase) ...
                + m.boundEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
            probs = m.computeRNAPolymeraseTUBindingProbabilities();
            probs = probs(:, 1);
            
            r.stateExpectations(r.activelyTranscribingIndex) = 0.25;
            r.stateExpectations(r.specificallyBoundIndex) = 0.10;
            r.stateExpectations(r.nonSpecificallyBoundIndex) = 0.0;
            r.stateExpectations(r.freeIndex) = 0.65;
            
            pStProbs = zeros(4, 4);
            pStProbs(r.freeIndex, r.activelyTranscribingIndex) = ...
                134 * 2 / (r.stateExpectations(r.activelyTranscribingIndex) * nPols) / ...
                (probs' * rna.lengths(rna.nascentIndexs));
            pStProbs(r.activelyTranscribingIndex, r.activelyTranscribingIndex) = ...
                1 - pStProbs(r.freeIndex, r.activelyTranscribingIndex);
            pStProbs(r.activelyTranscribingIndex, r.specificallyBoundIndex) = ...
                pStProbs(r.freeIndex, r.activelyTranscribingIndex) * ...
                r.stateExpectations(r.activelyTranscribingIndex) / ...
                r.stateExpectations(r.specificallyBoundIndex);
            pStProbs(r.specificallyBoundIndex, r.specificallyBoundIndex) = ...
                1 - pStProbs(r.activelyTranscribingIndex, r.specificallyBoundIndex);
            pStProbs(r.specificallyBoundIndex, r.freeIndex) = ...
                pStProbs(r.freeIndex, r.activelyTranscribingIndex) * ...
                r.stateExpectations(r.activelyTranscribingIndex) / ...
                r.stateExpectations(r.freeIndex);
            pStProbs(r.freeIndex, r.freeIndex) = ...
                1 - pStProbs(r.specificallyBoundIndex, r.freeIndex);
            pStProbs(:, r.nonSpecificallyBoundIndex) = pStProbs(:, r.freeIndex);
            m.stateTransitionProbabilities = pStProbs;
            
            occ = this.calcStateTransitions(600, true);
            occ = sum(occ, 2) / sum(occ(:));
            assertElementsAlmostEqual(r.stateExpectations, occ, 'relative', 7e-1, 0.1);
            
            assertEqual(0, min(m.substrates(m.substrateIndexs_ntp)));
            assertIn(corr(probs, m.RNAs), [0.60 1]);
            assertIn(180 / pi * acos((probs' * m.RNAs) / (sqrt(probs' * probs) * sqrt(m.RNAs' * m.RNAs))), [0 50]);
            
            %plot(probs, m.RNAs, '.'); 
            %line([0 sum(probs)], [0 sum(m.RNAs)]); 
            %xlim([0 max(probs)]), 
            %ylim([0 max(m.RNAs)])
        end
        
        function testStateTransitions_AT_SB_NSB_Free(this)
            m = this.process;
            r = m.rnaPolymerases;
            rna = m.rna;
            
            m.seed = 3000;
            m.seedRandStream();
            
            nPols = ...
                + m.enzymes(m.enzymeIndexs_rnaPolymerase) ...
                + m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) ...
                + m.boundEnzymes(m.enzymeIndexs_rnaPolymerase) ...
                + m.boundEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
            probs = m.computeRNAPolymeraseTUBindingProbabilities();
            probs = probs(:, 1);
            
            r.stateExpectations(r.activelyTranscribingIndex) = 0.25;
            r.stateExpectations(r.specificallyBoundIndex) = 0.10;
            r.stateExpectations(r.nonSpecificallyBoundIndex) = 0.65;
            r.stateExpectations(r.freeIndex) = 0.0;
            
            pStProbs = zeros(4, 4);
            pStProbs(r.nonSpecificallyBoundIndex, r.activelyTranscribingIndex) = ...
                134 * 2 / (r.stateExpectations(r.activelyTranscribingIndex) * nPols) / ...
                (probs' * rna.lengths(rna.nascentIndexs));
            pStProbs(r.activelyTranscribingIndex, r.activelyTranscribingIndex) = ...
                1 - pStProbs(r.nonSpecificallyBoundIndex, r.activelyTranscribingIndex);
            pStProbs(r.activelyTranscribingIndex, r.specificallyBoundIndex) = ...
                pStProbs(r.nonSpecificallyBoundIndex, r.activelyTranscribingIndex) * ...
                r.stateExpectations(r.activelyTranscribingIndex) / ...
                r.stateExpectations(r.specificallyBoundIndex);
            pStProbs(r.specificallyBoundIndex, r.specificallyBoundIndex) = ...
                1 - pStProbs(r.activelyTranscribingIndex, r.specificallyBoundIndex);
            pStProbs(r.specificallyBoundIndex, r.nonSpecificallyBoundIndex) = ...
                pStProbs(r.nonSpecificallyBoundIndex, r.activelyTranscribingIndex) * ...
                r.stateExpectations(r.activelyTranscribingIndex) / ...
                r.stateExpectations(r.nonSpecificallyBoundIndex);
            pStProbs(r.nonSpecificallyBoundIndex, r.nonSpecificallyBoundIndex) = ...
                1 - pStProbs(r.specificallyBoundIndex, r.nonSpecificallyBoundIndex);
            pStProbs(:, r.freeIndex) = pStProbs(:, r.nonSpecificallyBoundIndex);
            m.stateTransitionProbabilities = pStProbs;
            
            occ = this.calcStateTransitions(500, true);
            occ = sum(occ, 2) / sum(occ(:));
            assertElementsAlmostEqual(r.stateExpectations, occ, 'relative', 7e-1, 0.1);
            
            assertEqual(0, min(m.substrates(m.substrateIndexs_ntp)));
            assertIn(corr(probs, m.RNAs), [0.55 1]);
            assertIn(180 / pi * acos((probs' * m.RNAs) / (sqrt(probs' * probs) * sqrt(m.RNAs' * m.RNAs))), [0 55]);
            
            %plot(probs, m.RNAs, '.'); 
            %line([0 sum(probs)], [0 sum(m.RNAs)]); 
            %xlim([0 max(probs)]), 
            %ylim([0 max(m.RNAs)])
        end
        
        function testStateTransitions(this)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
            m = sim.process('Transcription');
            this.process = m;
            
            fitter = edu.stanford.covert.cell.sim.util.FitConstants(sim);
            paramVec = fitter.constructParameterVectorFromSimulation();
            fitter.applyParameterVectorToSimulation(paramVec);
            
            m = this.process;
            r = m.rnaPolymerases;
            
            occ = this.calcStateTransitions(1000, false);
            assertElementsAlmostEqual(r.stateExpectations, occ(:, end), 'relative', 1e-1, 0.2);
            assertElementsAlmostEqual(r.stateExpectations, sum(occ, 2) / sum(occ(:)), 'relative', 1e-1, 0.1);
        end
        
        function occ = calcStateTransitions(this, iterMax, checkSteadyState)
            %references
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            
            %assert transition probabilties and expectations
            if nargin >= 3 && checkSteadyState
                assertElementsAlmostEqual(m.stateTransitionProbabilities * r.stateExpectations, r.stateExpectations, ...
                    'relative', 1e-8);
            end
            assertEqual(1, sum(r.stateExpectations));
            assertTrue(all(m.stateTransitionProbabilities(:) >= 0));
            assertAllEqual(1, sum(m.stateTransitionProbabilities, 1));
            assertEqual(m.stateTransitionProbabilities(:, r.freeIndex), m.stateTransitionProbabilities(:, r.nonSpecificallyBoundIndex));
            
            %initialize at beginning of cell cycle
            nPols = ...
                + m.enzymes(m.enzymeIndexs_rnaPolymerase) ...
                + m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) ...
                + m.boundEnzymes(m.enzymeIndexs_rnaPolymerase) ...
                + m.boundEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
            probs = m.computeRNAPolymeraseTUBindingProbabilities();
            probs = probs(:, 1);
            
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_ntp) = 0;
            baseCounts = m.rna.baseCounts(m.rna.nascentIndexs, m.metabolite.nmpIndexs);
            ntpProd = baseCounts' * probs;
            ntpProd = 134.0758 * ntpProd / sum(ntpProd);
            
            m.enzymes = max(30, m.enzymes + m.boundEnzymes);
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = ...
                m.enzymes(m.enzymeIndexs_rnaPolymerase) + ...
                m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
            m.enzymes(m.enzymeIndexs_sigmaFactor) = ...
                m.enzymes(m.enzymeIndexs_sigmaFactor) + ...
                m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            m.boundEnzymes(:) = 0;
            m.RNAs(:) = 0;
            
            c.initialize();
            
            r.states(1:nPols) = r.freeValue;
            r.states(nPols + 1:end) = r.notExistValue;
            r.positionStrands(:) = 0;
            t.boundTranscriptProgress(:) = t.nullTranscriptValue;
            t.boundTranscriptionUnits(:) = t.nullTranscriptValue;
            t.boundTranscriptChromosome(:) = t.nullTranscriptValue;
            
            stateTransitionProbabilities = m.stateTransitionProbabilities;
            m.initializeState();
            m.stateTransitionProbabilities = stateTransitionProbabilities;
            
            %evolve
            occ = zeros(4, iterMax);
            for i = 1:iterMax
                tmpNtpProd = m.randStream.stochasticRound(ntpProd);
                m.substrates(m.substrateIndexs_ntp) = ...
                    m.substrates(m.substrateIndexs_ntp) + ...
                    tmpNtpProd;
                
                m.evolveState();
                
                occ(:, i) = r.stateOccupancies;
            end
        end
    end
    
    methods
        function testBindNearDSProteins(this)
            m = this.process;
            c = m.chromosome;
            
            rnaPolGblIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymerase);
            
            %monomer that can't be displaced by RNA polymerase
            idx = setdiff(...
                find(c.monomerDNAFootprintBindingStrandedness == c.dnaStrandedness_dsDNA), ...
                c.reactionBoundMonomer(c.reactionComplexCatalysisMatrix(:, rnaPolGblIdx) ~= 0));
            idx = idx(1);
            c.initialize();
            c.monomerBoundSites([1000 1]) = idx;
            [~,~ , ~, tfs] = c.setSiteProteinBound([1000 2], 1, [], [], rnaPolGblIdx, ...
                [], [], true, false, 1, false, []);
            assertFalse(tfs);
            
            %monomer that can be displaced by RNA polymerase
            idx = intersect(...
                find(c.monomerDNAFootprintBindingStrandedness == c.dnaStrandedness_dsDNA), ...
                c.reactionBoundMonomer(c.reactionComplexCatalysisMatrix(:, rnaPolGblIdx) ~= 0));
            idx = idx(1);
            c.initialize();
            c.monomerBoundSites([1000 1]) = idx;
            [~,~ , ~, tfs] = c.setSiteProteinBound([1000 2], 1, [], [], rnaPolGblIdx, ...
                [], [], true, false, 1, false, []);
            assertTrue(tfs);
            
            %complex that can't be displaced by RNA polymerase
            idx = setdiff(...
                find(c.complexDNAFootprintBindingStrandedness == c.dnaStrandedness_dsDNA), ...
                c.reactionBoundComplex(c.reactionComplexCatalysisMatrix(:, rnaPolGblIdx) ~= 0));
            idx = idx(1);
            c.initialize();
            c.complexBoundSites([1000 1]) = idx;
            [~,~ , ~, tfs] = c.setSiteProteinBound([1000 2], 1, [], [], rnaPolGblIdx, ...
                [], [], true, false, 1, false, []);
            assertFalse(tfs);
            
            %complex that can be displaced by RNA polymerase
            idx = intersect(...
                find(c.complexDNAFootprintBindingStrandedness == c.dnaStrandedness_dsDNA), ...
                c.reactionBoundComplex(c.reactionComplexCatalysisMatrix(:, rnaPolGblIdx) ~= 0));
            idx = idx(1);
            c.initialize();
            c.complexBoundSites([1000 1]) = idx;
            [~,~ , ~, tfs] = c.setSiteProteinBound([1000 2], 1, [], [], rnaPolGblIdx, ...
                [], [], true, false, 1, false, []);
            assertTrue(tfs);
        end
    end
    
    %test displacing
    methods
        function testDisplacingProteinOnOppositeStrand(this)
            m = this.process;
            c = m.chromosome;
            
            %binding position of RNA polymerase
            bindingIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymerase);
            bindingPosStrnds = [131 2];
            
            %already bound protein with large footprint
            boundIdx = 82;
            boundPosStrnds = [579889 1];
            
            %initialize chromosome with 1 bound protein
            c.initialize();
            c.complexBoundSites(boundPosStrnds) = boundIdx;
            
            %assert that bound protein footprint overlaps with binding
            %protein footprint
            assertTrue(boundPosStrnds(1) + c.complexDNAFootprints(boundIdx) - 1 >= bindingPosStrnds(1));
            
            %assert that bound protein can be released by RNA polymerase
            [~, releasableComplexs] = c.getReleasableProteins([], bindingIdx);
            assertTrue(any(boundIdx == releasableComplexs));
            
            %assert that RNA polymerase can't bind and that bound protein
            %not released
            [~, ~, ~, tfs] = c.setSiteProteinBound(bindingPosStrnds, 1, 1, [], bindingIdx, ...
                m.enzymeMonomerGlobalIndexs, m.enzymeComplexGlobalIndexs, ...
                true, false, 1, false, [], false);
            assertTrue(tfs);
            assertEqual(0, c.complexBoundSites(boundPosStrnds));
        end
        
        function testDisplacingProteinOnSameStrand(this)
            m = this.process;
            c = m.chromosome;
            
            %binding position of RNA polymerase
            bindingIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymerase);
            bindingPosStrnds = [131 2];
            
            %already bound protein with large footprint
            boundIdx = 82;
            boundPosStrnds = [579889 2];
            
            %initialize chromosome with 1 bound protein
            c.initialize();
            c.complexBoundSites(boundPosStrnds) = boundIdx;
            
            %assert that bound protein footprint overlaps with binding
            %protein footprint
            assertTrue(boundPosStrnds(1) + c.complexDNAFootprints(boundIdx) - 1 >= bindingPosStrnds(1));
            
            %assert that bound protein can be released by RNA polymerase
            [~, releasableComplexs] = c.getReleasableProteins([], bindingIdx);
            assertTrue(any(boundIdx == releasableComplexs));
            
            %assert that RNA polymerase can't bind and that bound protein
            %not released
            [~, ~, ~, tfs] = c.setSiteProteinBound(bindingPosStrnds, 1, 1, [], bindingIdx, ...
                m.enzymeMonomerGlobalIndexs, m.enzymeComplexGlobalIndexs, ...
                true, false, 1, false, [], false);
            assertTrue(tfs);
            assertEqual(0, c.complexBoundSites(boundPosStrnds));
        end
        
        function testNotDisplacingProteinOnOppositeStrand(this)
            m = this.process;
            c = m.chromosome;
            
            %binding position of RNA polymerase
            bindingIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymerase);
            bindingPosStrnds = [131 1];
            
            %already bound protein with large footprint
            boundIdx = 82;
            boundPosStrnds = [579889 2];
            
            %initialize chromosome with 1 bound protein
            c.initialize();
            c.complexBoundSites(boundPosStrnds) = boundIdx;
            
            %assert that bound protein footprint overlaps with binding
            %protein footprint
            assertTrue(boundPosStrnds(1) + c.complexDNAFootprints(boundIdx) - 1 >= bindingPosStrnds(1));
            
            %assert that bound protein can be released by RNA polymerase
            [~, releasableComplexs] = c.getReleasableProteins([], bindingIdx);
            assertTrue(any(boundIdx == releasableComplexs));
            
            %prevent bound protein from being released by RNA polymerase
            rxnIdx = find(c.reactionBoundComplex == boundIdx & c.reactionComplexCatalysisMatrix(:, bindingIdx));
            c.reactionBoundMonomer(rxnIdx) = [];
            c.reactionBoundComplex(rxnIdx) = [];
            c.reactionMonomerCatalysisMatrix(rxnIdx, :) = [];
            c.reactionComplexCatalysisMatrix(rxnIdx, :) = [];
            c.reactionThresholds = ...
                sum(c.reactionMonomerCatalysisMatrix, 2) + ...
                sum(c.reactionComplexCatalysisMatrix, 2);
            
            [~, releasableComplexs] = c.getReleasableProteins([], bindingIdx);
            assertFalse(any(boundIdx == releasableComplexs));
            
            %assert that RNA polymerase can't bind and that bound protein
            %not released
            [~, ~, ~, tfs] = c.setSiteProteinBound(bindingPosStrnds, 1, 1, [], bindingIdx, ...
                m.enzymeMonomerGlobalIndexs, m.enzymeComplexGlobalIndexs, ...
                true, false, 1, false, [], false);
            assertFalse(tfs);
            assertEqual(boundIdx, c.complexBoundSites(boundPosStrnds));
        end
        
        function testNotDisplacingProteinOnSameStrand(this)
            m = this.process;
            c = m.chromosome;
            
            %binding position of RNA polymerase
            bindingIdx = m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymerase);
            bindingPosStrnds = [131 1];
            
            %already bound protein with large footprint
            boundIdx = 82;
            boundPosStrnds = [579889 1];
            
            %initialize chromosome with 1 bound protein
            c.initialize();
            c.complexBoundSites(boundPosStrnds) = boundIdx;
            
            %assert that bound protein footprint overlaps with binding
            %protein footprint
            assertTrue(boundPosStrnds(1) + c.complexDNAFootprints(boundIdx) - 1 >= bindingPosStrnds(1));
            
            %assert that bound protein can be released by RNA polymerase
            [~, releasableComplexs] = c.getReleasableProteins([], bindingIdx);
            assertTrue(any(boundIdx == releasableComplexs));
            
            %prevent bound protein from being released by RNA polymerase
            rxnIdx = find(c.reactionBoundComplex == boundIdx & c.reactionComplexCatalysisMatrix(:, bindingIdx));
            c.reactionBoundMonomer(rxnIdx) = [];
            c.reactionBoundComplex(rxnIdx) = [];
            c.reactionMonomerCatalysisMatrix(rxnIdx, :) = [];
            c.reactionComplexCatalysisMatrix(rxnIdx, :) = [];
            c.reactionThresholds = ...
                sum(c.reactionMonomerCatalysisMatrix, 2) + ...
                sum(c.reactionComplexCatalysisMatrix, 2);
            
            [~, releasableComplexs] = c.getReleasableProteins([], bindingIdx);
            assertFalse(any(boundIdx == releasableComplexs));
            
            %assert that RNA polymerase can't bind and that bound protein
            %not released
            [~, ~, ~, tfs] = c.setSiteProteinBound(bindingPosStrnds, 1, 1, [], bindingIdx, ...
                m.enzymeMonomerGlobalIndexs, m.enzymeComplexGlobalIndexs, ...
                true, false, 1, false, [], false);
            assertFalse(tfs);
            assertEqual(boundIdx, c.complexBoundSites(boundPosStrnds));
        end
        
        function testBindingOppositeDamage(this)
            m = this.process;
            c = m.chromosome;
            
            %initialize with 1 strand break
            c.initialize();
            c.strandBreaks([414498 1]) = 1;
            
            %assert region opposite strand break is accessible
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions(...
                [], m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymerase));
            assertEqual([414500 1; 1 2], rgnPosStrnds);
            assertEqual([c.sequenceLen-2; c.sequenceLen], rgnLens);
            
            %check that RNA polymerase can bind opposite strand break
            tfs = m.bindProteinToChromosome([414456 2], m.enzymeIndexs_rnaPolymerase, 35, [], [], false, 1, false, [], false);
            assertTrue(tfs);
        end
        
        function testNotBindingOverDamage(this)
            m = this.process;
            c = m.chromosome;
            
            %initialize with 1 strand break
            c.initialize();
            c.strandBreaks([414498 1]) = 1;
            
            %assert region opposite strand break is accessible
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions(...
                [], m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymerase));
            assertEqual([414500 1; 1 2], rgnPosStrnds);
            assertEqual([c.sequenceLen-2; c.sequenceLen], rgnLens);
            
            %check that RNA polymerase can bind opposite strand break
            tfs = m.bindProteinToChromosome([414456 1], m.enzymeIndexs_rnaPolymerase, 35, [], [], false, 1, false, [], false);
            assertFalse(tfs);
        end
    end    
    
    methods
        function testGeneEssentiality(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            
            m.transcriptionUnitBindingProbabilities(:) = 1;
            m.substrates(:) = 1e12;
            m.enzymes(:) = 1e6;
            numRNAPolymerase = 10;
            m.enzymes(m.enzymeIndexs_rnaPolymerase) = numRNAPolymerase;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0;
            m.boundEnzymes(:) = 0;
            
            c.initialize();
            r.states = repmat(r.notExistValue, numRNAPolymerase, 1);
            r.positionStrands = zeros(numRNAPolymerase, 2);
            r.transcriptionFactorBindingProbFoldChange(:) = 1;
            r.supercoilingBindingProbFoldChange(:) = 1;
            t.boundTranscriptionUnits = repmat(t.nullTranscriptValue, numRNAPolymerase, 1);
            t.boundTranscriptProgress = repmat(t.nullTranscriptValue, numRNAPolymerase, 1);
            t.boundTranscriptChromosome = repmat(t.nullTranscriptValue, numRNAPolymerase, 1);
            
            this.helpTestGeneEssentiality({
                'MG_022';
                'MG_141';
                'MG_177';
                'MG_249';
                'MG_282';
                'MG_340';
                'MG_341'},...
                @(m,i) any(m.RNAs > i.RNAs),...
                struct('lengthSec', 100));
        end
    end
    
    %helpers
    methods
        function prepareSimpleFixture(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            
            %metabolites
            m.substrates(m.substrateIndexs_ntp)         = repmat(100,4,1);
            m.substrates(m.substrateIndexs_nmp)         = zeros(4,1);
            m.substrates(m.substrateIndexs_adp)         = 0;
            m.substrates(m.substrateIndexs_diphosphate) = 0;
            m.substrates(m.substrateIndexs_water)       = 1e6;
            m.substrates(m.substrateIndexs_hydrogen)    = 0;
            
            %enzymes
            m.enzymes(m.enzymeIndexs_sigmaFactor)             =3;
            m.enzymes(m.enzymeIndexs_elongationFactor)        =3;
            m.enzymes(m.enzymeIndexs_terminationFactor)       =3;
            m.enzymes(m.enzymeIndexs_antiterminationFactor)   =3;
            m.enzymes(m.enzymeIndexs_rnaPolymerase)           =5;
            m.enzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme) =0;
            
            %bound enzymes
            m.boundEnzymes(:)=0;
            
            %RNA Sequence Properties
            t.transcriptionUnitSequences = {...
                'ACGUUGCAACGUUGCAACGUUGCAACGUUGCAACGUUGCAACGUUGCACA';...
                'UGCAGACGUAACGGUAGCGCUGCCGCGGGCACGUGCACGUUUGACGACGA';...
                'ACUGGCGACGUGCAGGCAGUUGCGCACUAUGCGAGGCACGCGUACGACGG';...
                'ACGUUUGCAACGUAAGCAUUUACGGACGUAACGCGACGUGCGUCGACGCG';...
                'GCGUUCGCAACGGUACGUAGGCAUGGUACGCGCGCCGUCGUGCAUGCACG'};
            t.transcriptionUnitDirections = [1; 1; -1; -1; -1];
            t.transcriptionUnitFivePrimeCoordinates = [500; 6050; 10000; 10500; 37500];
            c.transcriptionUnitStartCoordinates = [500; 6050; 10049; 10549; 37549];
            c.transcriptionUnitStrands = [1; 1; 2; 2; 2];
            t.transcriptionUnitLengths = [50; 50; 50; 50; 50];
            c.transcriptionUnitLengths = t.transcriptionUnitLengths;
            m.transcriptionUnitIndexs_DnaAR12345Boxes = zeros(0, 1);
            %c.polymerizedTranscriptionUnits = [1, 0; 1, 0; 1, 0; 1, 0; 1, 0];
            m.transcriptionUnitBaseCounts = [... %Order: ACGU
                13, 13, 12, 12;...
                10, 14, 18, 8 ;...
                10, 14, 19, 7 ;...
                12, 13, 15, 10;...
                8 , 15, 18, 9 ];
            m.rna.molecularWeights = [15201.7; 15365.8; 15404.8; 15295.7; 15318.7];
            m.rna.nascentIndexs = (1:5)';
            m.transcriptionUnitBindingProbabilities = [0.1; 0.15; 0.2; 0.25; 0.3];
            r.transcriptionFactorBindingProbFoldChange = ones(5, 2);
            r.supercoilingBindingProbFoldChange = ones(5, 2);
            
            %Parameters
            m.rnaPolymeraseElongationRate = 100;
            r.stateExpectations = zeros(4, 1);
            r.stateExpectations(r.activelyTranscribingIndex) = 0.5;
            r.stateExpectations(r.specificallyBoundIndex) = 0.5;
            m.stateTransitionProbabilities = zeros(4, 4);
            m.stateTransitionProbabilities(r.activelyTranscribingIndex, :) = 0.5;
            m.stateTransitionProbabilities(r.specificallyBoundIndex, :) = 0.5;
            
            t.genomeLength = 580076;
            
            %initial conditions
            m.RNAs = zeros(5,1);
            r.states = repmat(r.freeValue,5,1);
            r.positionStrands = zeros(5, 2);
            t.boundTranscriptionUnits = repmat(t.nullTranscriptValue, 5, 1);
            t.boundTranscriptProgress = repmat(t.nullTranscriptValue, 5, 1);
            t.boundTranscriptChromosome = repmat(t.nullTranscriptValue, 5, 1);
            
            c.initialize();
        end
        
        function isStateValid(this)
            m = this.process;
            r = m.rnaPolymerases;
            t = m.transcripts;
            c = m.chromosome;
            rPolyGblIdx = [
                m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymerase);
                m.enzymeGlobalIndexs(m.enzymeIndexs_rnaPolymeraseHoloenzyme);
                ];
            
            %check that state is synchronized
            assertAllEqual(true, m.enzymes >= 0);
            assertAllEqual(true, m.boundEnzymes >= 0);
            assertEqual(nnz(c.complexBoundSites == rPolyGblIdx(1)), m.boundEnzymes(m.enzymeIndexs_rnaPolymerase));
            assertEqual(nnz(c.complexBoundSites == rPolyGblIdx(2)), m.boundEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme));
            assertEqual(sum(r.states == r.specificallyBoundValue | r.states == r.activelyTranscribingValue),...
                m.boundEnzymes(m.enzymeIndexs_rnaPolymeraseHoloenzyme));
            assertEqual(0, m.boundEnzymes(m.enzymeIndexs_sigmaFactor));
            assertEqual(nnz(c.complexBoundSites == rPolyGblIdx(1)), sum(...
                r.states == r.nonSpecificallyBoundValue | ...
                r.states  > r.activelyTranscribingValue));
            assertEqual(nnz(c.complexBoundSites == rPolyGblIdx(2)), sum(...
                r.states == r.specificallyBoundValue | ...
                r.states == r.activelyTranscribingValue));
            assertEqual(m.enzymes(m.enzymeIndexs_rnaPolymerase), ...
                sum(r.states == r.freeValue));
            assertEqual(r.states' >= r.activelyTranscribingValue', t.boundTranscriptProgress' >= r.activelyTranscribingValue');
            
            states = zeros(size(r.states));
            states(r.states >= r.activelyTranscribingValue & isodd(r.positionStrands(:,2))) = ...
                r.positionStrands(r.states >= r.activelyTranscribingValue & isodd(r.positionStrands(:,2)), 1) - ...
                t.transcriptionUnitFivePrimeCoordinates(t.boundTranscriptionUnits(r.states >= r.activelyTranscribingValue & isodd(r.positionStrands(:,2)))) + 1;
            states(r.states >= r.activelyTranscribingValue & iseven(r.positionStrands(:,2))) = ...
                t.transcriptionUnitFivePrimeCoordinates(t.boundTranscriptionUnits(r.states >= r.activelyTranscribingValue & iseven(r.positionStrands(:,2)))) - ...
                r.positionStrands(r.states >= r.activelyTranscribingValue & iseven(r.positionStrands(:,2)), 1) + 1;
            
            assertEqual(states(r.states >= r.activelyTranscribingValue), r.states(r.states >= r.activelyTranscribingValue))
            
            positionStrands = r.positionStrands;
            positionStrands(isodd(positionStrands(:,2)), 1) = ...
                positionStrands(isodd(positionStrands(:,2)), 1) - ...
                m.enzymeDNAFootprints5Prime(m.enzymeIndexs_rnaPolymerase);
            positionStrands(iseven(positionStrands(:,2)), 1) = ...
                positionStrands(iseven(positionStrands(:,2)), 1) - ...
                m.enzymeDNAFootprints3Prime(m.enzymeIndexs_rnaPolymerase);
            assertEqual(find(c.complexBoundSites == rPolyGblIdx(1) | c.complexBoundSites == rPolyGblIdx(2)), sortrows(positionStrands(any(r.positionStrands,2), :), [2 1]))
            
            positionStrands = zeros(size(r.states, 1), 2);
            actPols = r.states >= r.activelyTranscribingValue;
            sbPols = r.states == r.specificallyBoundValue;
            nsbPols = r.states == r.nonSpecificallyBoundValue;
            positionStrands(actPols, :) = [...
                t.transcriptionUnitFivePrimeCoordinates(t.boundTranscriptionUnits(actPols)) + ...
                (2 * t.transcriptionUnitDirections(t.boundTranscriptionUnits(actPols)) -1) .* (r.states(actPols)-1) ...
                2 - t.transcriptionUnitDirections(t.boundTranscriptionUnits(actPols))];
            positionStrands(sbPols, :) = [...
                t.transcriptionUnitFivePrimeCoordinates(t.boundTranscriptionUnits(sbPols)) + ...
                (2 * t.transcriptionUnitDirections(t.boundTranscriptionUnits(sbPols)) - 1) * -1 ...
                2 - t.transcriptionUnitDirections(t.boundTranscriptionUnits(sbPols))];
            positionStrands(nsbPols, :) = r.positionStrands(nsbPols, :);
            
            assertEqual(true(sum(r.states >= r.activelyTranscribingValue), 1), ...
                r.states(r.states >= r.activelyTranscribingValue) <= ...
                t.transcriptionUnitLengths(t.boundTranscriptionUnits(r.states >= r.activelyTranscribingValue)) + 1);
            
            assertTrue(all(ismember(t.boundTranscriptChromosome, [0; 1])));
            
            assertEqual(r.positionStrands, positionStrands)
        end
    end
end
