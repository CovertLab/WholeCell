%RNA production medium test
%
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/10/2010
classdef RNA_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = RNA_Test(name)
            this = this@TestCase(name);
        end
        
        function setUp(this)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'Metabolism'
                'Transcription'
                'RNAProcessing'
                'RNAModification'
                'tRNAAminoacylation'
                'RNADecay'
                });
            
            this.simulation = sim;
        end
        
        function testRnaProduction(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            %% references
            sim = this.simulation;
            %this.seedSimulation(sim, 10);
            %sim.applyOptions('verbosity', 1);
            
            g = sim.gene;
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            rna = sim.state('Rna');
            time = sim.state('Time');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            transcription = sim.process('Transcription');
            
            %% constants
            rnaMWs = rna.molecularWeights(rna.aminoacylatedIndexs);
            
            %% turn off decay
            rna.decayRates(setdiff(1:end, rna.intergenicIndexs)) = 0;
            
            %% keep track of initial state
            initRNAs = ...
                + rna.counts(rna.matureIndexs, sim.compartment.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs, sim.compartment.cytosolIndexs);
            initRNAs(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                initRNAs(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(pc.counts(pc.matureIndexs, :) + pc.counts(pc.boundIndexs, :), 2);
            
            initMonomers = ...
                + pm.counts(pm.matureIndexs, :) ...
                + pm.counts(pm.boundIndexs, :);
            initComplexs = ...
                + pc.counts(pc.matureIndexs, :) ...
                + pc.counts(pc.boundIndexs, :);
            
            %% simulate
            iterMax = 1000;
            
            ntps = zeros(4, iterMax);
            if sim.verbosity > 0
                fprintf('%5s %6s %6s %6s %6s\n', 'Iter ', ' ATP  ', ' CTP  ', ' GTP  ', ' UTP  ');
                fprintf('%5s %6s %6s %6s %6s\n', '=====', '======', '======', '======', '======');
            end
            for i = 1:iterMax
                if sim.verbosity > 0 && mod(i, 100) == 1
                    fprintf('%5d %6d %6d %6d %6d\n', i, ...
                        transcription.substrates(transcription.substrateIndexs_ntp(1)), ...
                        transcription.substrates(transcription.substrateIndexs_ntp(2)), ...
                        transcription.substrates(transcription.substrateIndexs_ntp(3)), ...
                        transcription.substrates(transcription.substrateIndexs_ntp(4)));
                end
                
                %mock protein synthesis
                pm.counts(pm.matureIndexs, :) = ...
                    + pm.counts(pm.matureIndexs, :) ...
                    + pm.randStream.stochasticRound(initMonomers * log(2) / time.cellCycleLength * exp(i * log(2) / time.cellCycleLength));
                pc.counts(pc.matureIndexs, :) = ...
                    + pc.counts(pc.matureIndexs, :) ....
                    + pc.randStream.stochasticRound(initComplexs * log(2) / time.cellCycleLength * exp(i * log(2) / time.cellCycleLength));
                
                %simulate               
                sim.evolveState();
                
                ntps(:, i) = transcription.substrates(transcription.substrateIndexs_ntp);
            end
            
            newRNAs = ...
                - initRNAs ...
                + rna.counts(rna.matureIndexs, sim.compartment.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs, sim.compartment.cytosolIndexs);
            newRNAs(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                newRNAs(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(initComplexs, 2);
            
            %% assertions
            
            %NTPs are limiting
            assertIn(min(min(ntps(:, end-100+1:end))), [0 max(2000, min(ntps(:, 1)))]);
            
            %NTPs are mostly used -- not too many are free
            assertIn(...
                (m.counts(m.ntpIndexs([2 4]), sim.compartment.cytosolIndexs)' * m.molecularWeights(m.ntpIndexs([2 4]))) / ...
                (rna.counts(:, sim.compartment.cytosolIndexs)' * rna.molecularWeights), ...
                [0 0.03]);
            
            %mature/aminoacylated RNA counts increased
            assertAllEqual(true, newRNAs >= 0);
            assertElementsAlmostEqual(newRNAs' * rnaMWs / ConstantUtil.nAvogadro, ...
                (exp(iterMax * log(2) / time.cellCycleLength) - 1) * ...
                mass.initialFractionNTPsInRNAs * mass.cellInitialDryWeight * mass.dryWeightFractionRNA, ...
                'relative', 2.5e-1, 0);
            
            assertIn(sum(sum(rna.counts(rna.nascentIndexs, :))), [0 20]);
            assertIn(sum(sum(rna.counts(rna.processedIndexs, :))), [0 40]);
            assertIn(sum(sum(rna.counts(rna.intergenicIndexs, :))), [0 10]);
            assertAllEqual(0, rna.counts(rna.boundIndexs, :));
            assertAllEqual(0, rna.counts(rna.damagedIndexs, :));
            assertAllEqual(0, rna.counts(rna.misfoldedIndexs, :));
            
            %RNA production matches expectations
            rates = rna.nascentRNAMatureRNAComposition * transcription.computeRNAPolymeraseTUBindingProbabilities();
            rates = rates(:, 1);
            assertIn(corr(rates, newRNAs), [0.50 1]);
            assertIn(180 / pi * acos(rates' * newRNAs / ...
                (sqrt(rates' * rates) * sqrt(newRNAs' * newRNAs))), ...
                [0 60]);
        end
        
        function testRnaProductionAndDecay(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            %% references
            sim = this.simulation;
            %sim.applyOptions('verbosity', 1);
            this.seedSimulation(sim, 1000);
            
            g = sim.gene;
            comp = sim.compartment;
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            rna = sim.state('Rna');
            time = sim.state('Time');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            rnaPol = sim.state('RNAPolymerase');
            transcript = sim.state('Transcript');
            transcription = sim.process('Transcription');
            
            %% keep track of initial state
            initRNAs = ...
                + rna.counts(rna.matureIndexs, comp.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs, comp.cytosolIndexs);
            initRNAs(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                initRNAs(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * ...
                sum(pc.counts(pc.matureIndexs, :) + pc.counts(pc.boundIndexs, :), 2);
            initMonomers = ...
                + pm.counts(pm.matureIndexs, :) ...
                + pm.counts(pm.boundIndexs, :);
            initComplexs = ...
                + pc.counts(pc.matureIndexs, :) ...
                + pc.counts(pc.boundIndexs, :);
            
            %% simulate
            iterMax = 30000;
            
            rnaExp = zeros(size(initRNAs));
            rnaProd = zeros(size(rna.nascentIndexs));
            oldNascentRNAs = transcription.RNAs;
            rnaPolStateOcc = zeros(4, iterMax);
            ntps = zeros(4, iterMax);
            ndps = zeros(4, iterMax);
            nmps = zeros(4, iterMax);
            phosphates = zeros(3, iterMax);
            boundTUs = zeros(numel(rna.nascentIndexs), iterMax);
            rnaDryWt = zeros(1, iterMax);
            transcriptDryWt = zeros(1, iterMax);
            nucDryWt = zeros(1, iterMax);
            if sim.verbosity > 0
                fprintf('%5s %6s %6s %6s %6s\n', 'Iter ', ' ATP  ', ' CTP  ', ' GTP  ', ' UTP  ');
                fprintf('%5s %6s %6s %6s %6s\n', '=====', '======', '======', '======', '======');
            end
            for i = 1:iterMax
                if sim.verbosity > 0 && mod(i, 100) == 1
                    fprintf('%5d %6d %6d %6d %6d\n', i, ...
                        transcription.substrates(transcription.substrateIndexs_ntp(1)), ...
                        transcription.substrates(transcription.substrateIndexs_ntp(2)), ...
                        transcription.substrates(transcription.substrateIndexs_ntp(3)), ...
                        transcription.substrates(transcription.substrateIndexs_ntp(4)));
                end
                
                %mock protein synthesis
                pm.counts(pm.matureIndexs, :) = ...
                    + pm.counts(pm.matureIndexs, :) ...
                    + pm.randStream.stochasticRound(initMonomers * log(2) / time.cellCycleLength * exp(i * log(2) / time.cellCycleLength));
                pc.counts(pc.matureIndexs, :) = ...
                    + pc.counts(pc.matureIndexs, :) + ...
                    + pc.randStream.stochasticRound(initComplexs * log(2) / time.cellCycleLength * exp(i * log(2) / time.cellCycleLength));
                
                %simulate
                sim.evolveState();
                
                %store data
                rnaExp = rnaExp ...
                    + rna.counts(rna.matureIndexs, comp.cytosolIndexs) ...
                    + rna.counts(rna.aminoacylatedIndexs, comp.cytosolIndexs);
                rnaExp(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                    rnaExp(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                    sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(initComplexs, 2);
                rnaProd = rnaProd + max(0, transcription.RNAs - oldNascentRNAs);
                oldNascentRNAs = transcription.RNAs;
                rnaPolStateOcc(:, i) = [
                    rnaPol.nActive
                    rnaPol.nSpecificallyBound
                    rnaPol.nNonSpecificallyBound
                    rnaPol.nFree
                    ];
                ntps(:, i) = transcription.substrates(transcription.substrateIndexs_ntp);
                ndps(:, i) = m.counts(m.ndpIndexs, comp.cytosolIndexs);
                nmps(:, i) = m.counts(m.nmpIndexs, comp.cytosolIndexs);
                phosphates(:, i) = m.counts([m.diphosphateIndexs; m.phosphateIndexs; m.hydrogenIndexs], comp.cytosolIndexs);
                boundTUs(:, i) = histc(transcript.boundTranscriptionUnits(rnaPol.states >= rnaPol.activelyTranscribingValue), (1:numel(rna.nascentIndexs))');
                
                rnaDryWt(:, i) = rna.dryWeight(1);
                transcriptDryWt(:, i) = transcript.dryWeight;
                nucDryWt(:, i) = (...
                    + ntps(:, i)' * m.molecularWeights(m.ntpIndexs) ...
                    + ndps(:, i)' * m.molecularWeights(m.ndpIndexs) ...
                    + nmps(:, i)' * m.molecularWeights(m.nmpIndexs)) / ConstantUtil.nAvogadro;
            end
            
            %% assertions
            
            %NTPs are limiting
            assertIn(min(min(ntps(:, end-100+1:end))), [0 exp(log(2) * iterMax / time.cellCycleLength) * min(ntps(:, 1))]);
            assertIn(min(min(ndps(:, end-100+1:end))), [0 exp(log(2) * iterMax / time.cellCycleLength) * min(ndps(:, 1))]);
            assertIn(min(min(nmps(:, end-100+1:end))), [0 exp(log(2) * iterMax / time.cellCycleLength) * min(nmps(:, 1))]);
            assertIn(min(transcriptDryWt(:, end-100+1:end)), [0 2.5 * exp(log(2) * iterMax / time.cellCycleLength) * transcriptDryWt(1)]);
            
            %NTPs are mostly used -- not too many are free
            assertIn(...
                (m.counts(m.ntpIndexs([2 4]), comp.cytosolIndexs)' * m.molecularWeights(m.ntpIndexs([2 4]))) / ...
                ((rna.dryWeight(1) + transcript.dryWeight) * edu.stanford.covert.util.ConstantUtil.nAvogadro), ...
                [0 0.05]);
            
            rnas = ...
                + rna.counts(rna.matureIndexs, comp.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs, comp.cytosolIndexs);
            rnas(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                rnas(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(initComplexs, 2);
            assertElementsAlmostEqual(rnas' * rna.molecularWeights(rna.aminoacylatedIndexs) / edu.stanford.covert.util.ConstantUtil.nAvogadro, ...
                exp(iterMax * log(2) / time.cellCycleLength) * mass.cellInitialDryWeight * mass.dryWeightFractionRNA, ...
                'relative', 2e-1, 0);
            
            %total NTP incorporation similar to production by metabolism
            assertElementsAlmostEqual(...
                (edu.stanford.covert.util.ComputationUtil.invertCompositionMatrix(rna.nascentRNAMatureRNAComposition) * ...
                (mass.cellInitialDryWeight * mass.dryWeightFractionRNA * edu.stanford.covert.util.ConstantUtil.nAvogadro / ...
                (rna.expression(rna.matureIndexs)' * rna.molecularWeights(rna.matureIndexs)) * rna.expression(rna.matureIndexs) .* ...
                (log(2) / time.cellCycleLength + rna.decayRates(rna.matureIndexs))))' * rna.lengths(rna.nascentIndexs) * ...
                time.cellCycleLength / log(2) * (exp(log(2) / time.cellCycleLength * iterMax)-1), ...
                rnaProd' * rna.lengths(rna.nascentIndexs), ...
                'relative', 0.40);
            
            %RNA expression matches expectations
            assertIn(corr(rna.expression(rna.matureIndexs), rnaExp), [0.9 1]);
            assertIn(180 / pi * acos(rna.expression(rna.matureIndexs)' * rnaExp / ...
                (sqrt(rna.expression(rna.matureIndexs)' * rna.expression(rna.matureIndexs)) * sqrt(rnaExp' * rnaExp))), ...
                [0 20]);
            
            %RNA polymerase occupancy matches expectations
            assertElementsAlmostEqual(rnaPol.stateExpectations, ...
                sum(rnaPolStateOcc, 2) / sum(rnaPolStateOcc(:)), ...
                'relative', 0.75, 0.20); %todo: tighten
            
            %RNAs matured
            assertIn(sum(sum(rna.counts(rna.nascentIndexs, :))), [0 20]);
            assertIn(sum(sum(rna.counts(rna.processedIndexs, :))), [0 30]);
            assertIn(sum(sum(rna.counts(rna.intergenicIndexs, :))), [0 10]);
            assertAllEqual(0, rna.counts(rna.boundIndexs, :));
            assertAllEqual(0, rna.counts(rna.damagedIndexs, :));
            assertAllEqual(0, rna.counts(rna.misfoldedIndexs, :));
        end
        
        function sim = seedSimulation(~, sim, seed)
            if sim.verbosity > 0
                fprintf('Seed: %d\n', seed)
            end
            
            sim.applyOptions('seed', seed);
            sim.seedRandStream();
            
            for i = 1:numel(sim.states)
                o = sim.states{i};
                o.seed = seed;
                o.seedRandStream();
            end
            
            for i = 1:numel(sim.processes)
                o = sim.processes{i};
                o.seed = seed;
                o.seedRandStream();
            end
        end
    end
end
