%Full Length simulation test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/12/2011
classdef Simulation_Large_Test < TestCase
    methods
        function this = Simulation_Large_Test(name)
            this = this@TestCase(name);
        end
    end
    
    methods
        function testRNAAndProteinSynthesisAndDecay(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            %% load fixture
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'Metabolism'
                ...
                'Transcription'
                'RNAProcessing'
                'RNAModification'
                'tRNAAminoacylation'
                'RNADecay'
                ...
                'RibosomeAssembly'
                'Translation'
                'ProteinProcessingI'
                'ProteinTranslocation'
                'ProteinProcessingII'
                'ProteinFolding'
                'ProteinModification'
                'MacromolecularComplexation'
                'TerminalOrganelleAssembly'
                'ProteinActivation'
                'ProteinDecay'});
            %sim.applyOptions('verbosity', 1);
            %this.seedRandStream(sim, round(mod(now, 1) * 1e7));
            
            %% references
            comp = sim.compartment;
            g = sim.gene;
            
            time = sim.state('Time');
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            mr = sim.state('MetabolicReaction');
            pol = sim.state('Polypeptide');
            transcript = sim.state('Transcript');
            rnaPol = sim.state('RNAPolymerase');
            
            dcy = sim.process('ProteinDecay');
            ta = sim.process('tRNAAminoacylation');
            transcription = sim.process('Transcription');
            
            %% constants
            monMWs = pm.molecularWeights(pm.matureIndexs);
            
            metNucComp = edu.stanford.covert.cell.sim.util.FitConstants_Test.getNucleotideCompositionMatrix(m);
            rnaNucBaseComp = (rna.baseCounts * edu.stanford.covert.cell.sim.util.FitConstants_Test.getNucleotideCompositionMatrix(m)')';
            monNucBaseComp = (pm.baseCounts * edu.stanford.covert.cell.sim.util.FitConstants_Test.getNucleotideCompositionMatrix(m)')';
            cpxNucBaseComp = (pc.baseCounts * edu.stanford.covert.cell.sim.util.FitConstants_Test.getNucleotideCompositionMatrix(m)')';
            
            %% keep track of initial state
            initRNAs = ...
                + rna.counts(rna.matureIndexs, sim.compartment.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs, sim.compartment.cytosolIndexs);
            initRNAs(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                initRNAs(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(pc.counts(pc.matureIndexs, :) + pc.counts(pc.boundIndexs, :), 2);
            initMetabolites = m.counts;
            initMonomers = pm.counts;
            initComplexs = pc.counts;
            
            %% simulate
            iterMax = time.replicationInitiationDuration;
            
            rnaExp = zeros(size(initRNAs));
            growth = 0;
            simMrnaExp = zeros(size(rna.matureMRNAIndexs));
            tRNAs = zeros(numel(rna.matureTRNAIndexs), iterMax);
            freedAAs = zeros(numel(dcy.substrateIndexs_aminoAcids), iterMax);
            aaCnts = zeros(numel(ta.substrateIndexs_aminoAcids), iterMax);
            rnaWt = zeros(1, iterMax);
            transcriptWt = zeros(1, iterMax);
            ntpWt = zeros(1, iterMax);
            rnaPolOcc = zeros(4, iterMax);
            nucComp_rna = zeros(4, iterMax);
            nucComp_mon = zeros(4, iterMax);
            nucComp_cpx = zeros(4, iterMax);
            nucComp_met = zeros(4, iterMax);
            nucComp_trs = zeros(4, iterMax);
            nucComp = zeros(4, iterMax);
            if sim.verbosity > 0
                fprintf('%5s %8s %7s %7s %7s %7s %7s %7s %7s %7s %8s %7s\n', 'Iter ', 'Min NTP ', 'Sum AxP', 'Sum CxP', 'Sum GxP', 'Sum UxP', 'Sum NTP', 'Sum NDP', 'Sum NMP', 'Sum NxP', 'Min AA  ', 'Sum AA ');
                fprintf('%5s %8s %7s %7s %7s %7s %7s %7s %7s %7s %8s %7s\n', '=====', '========', '=======', '=======', '=======', '=======', '=======', '=======', '=======', '=======', '========', '=======');
            end
            for iter = 1:iterMax
                if sim.verbosity > 0 && mod(iter, 100) == 1
                    [val1, idx1] = min(m.counts(m.ntpIndexs, comp.cytosolIndexs));
                    [val2, idx2] = min(m.counts(m.aminoAcidIndexs(1:20), comp.cytosolIndexs));
                    fprintf('%5d %3s %4d %7d %7d %7d %7d %7d %7d %7d %7d %3s %4d %7d\n', iter, ...
                        m.wholeCellModelIDs{m.ntpIndexs(idx1)}, val1, ...
                        sum(m.counts([m.ntpIndexs(1); m.ndpIndexs(1); m.nmpIndexs(1)], comp.cytosolIndexs)), ...
                        sum(m.counts([m.ntpIndexs(2); m.ndpIndexs(2); m.nmpIndexs(2)], comp.cytosolIndexs)), ...
                        sum(m.counts([m.ntpIndexs(3); m.ndpIndexs(3); m.nmpIndexs(3)], comp.cytosolIndexs)), ...
                        sum(m.counts([m.ntpIndexs(4); m.ndpIndexs(4); m.nmpIndexs(4)], comp.cytosolIndexs)), ...
                        sum(m.counts(m.ntpIndexs, comp.cytosolIndexs)), ...
                        sum(m.counts(m.ndpIndexs, comp.cytosolIndexs)), ...
                        sum(m.counts(m.nmpIndexs, comp.cytosolIndexs)), ...
                        sum(m.counts([m.ntpIndexs; m.ndpIndexs; m.nmpIndexs], comp.cytosolIndexs)), ...
                        m.wholeCellModelIDs{m.aminoAcidIndexs(idx2)}, val2, ...
                        sum(m.counts(m.aminoAcidIndexs(1:20), comp.cytosolIndexs)));
                end
                
                sim.evolveState();
                
                rnaExp = rnaExp ...
                    + rna.counts(rna.matureIndexs, sim.compartment.cytosolIndexs) ...
                    + rna.counts(rna.aminoacylatedIndexs, sim.compartment.cytosolIndexs);
                rnaExp(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                    rnaExp(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                    sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(initComplexs(pc.matureIndexs, :) + initComplexs(pc.boundIndexs, :), 2);
                
                growth = growth + mr.growth;
                simMrnaExp = simMrnaExp + rna.counts(rna.matureIndexs(rna.matureMRNAIndexs), comp.cytosolIndexs);
                tRNAs(:, iter) = ...
                    + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) ...
                    + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs);
                freedAAs(:, iter) = dcy.substrates(dcy.substrateIndexs_aminoAcids);
                aaCnts(:, iter) = ta.substrates(ta.substrateIndexs_aminoAcids);
                rnaWt(1, iter) = sum(rna.dryWeight);
                transcriptWt(1, iter) = transcript.dryWeight;
                ntpWt(1, iter) = m.counts(m.ntpIndexs, comp.cytosolIndexs)' * m.molecularWeights(m.ntpIndexs) / ConstantUtil.nAvogadro;
                rnaPolOcc(:, iter) = rnaPol.stateOccupancies;
                nucComp_rna(:, iter) = rnaNucBaseComp * sum(rna.counts, 2);
                nucComp_mon(:, iter) = monNucBaseComp * sum(pm.counts, 2);
                nucComp_cpx(:, iter) = cpxNucBaseComp * sum(pc.counts, 2);
                nucComp_met(:, iter) = metNucComp * m.counts(:, comp.cytosolIndexs);
                nucComp_trs(:, iter) = transcript.totalBaseCounts';
                nucComp(:, iter) = ...
                    + nucComp_rna(:, iter) ...
                    + nucComp_mon(:, iter) ...
                    + nucComp_cpx(:, iter) ...
                    + nucComp_met(:, iter) ...
                    + nucComp_trs(:, iter);
                
                this.assertStateIsValid_Protein(sim);
            end
            
            %% assertions - growth
            assertElementsAlmostEqual(log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength), ...
                mr.growth, 'relative', 0.3, 0);
            assertElementsAlmostEqual(exp(log(2) * iter / time.cellCycleLength) - 1, ...
                growth, 'relative', 0.3, 0);
            assertElementsAlmostEqual(mass.cellInitialDryWeight * exp(log(2) * iter / time.cellCycleLength), ...
                sum(mass.cellDry), 'relative', 0.3, 0);
            
            %% assertions -- RNA
            
            %NTPs are limiting
            assertIn(min(transcription.substrates(transcription.substrateIndexs_ntp)), [0 5000]);
            
            %NTPs are mostly used -- not too many are free
            assertIn(...
                (m.counts(m.ntpIndexs([2 4]), sim.compartment.cytosolIndexs)' * m.molecularWeights(m.ntpIndexs([2 4]))) / ...
                (rna.counts(:, sim.compartment.cytosolIndexs)' * rna.molecularWeights), ...
                [0 0.10]);
            
            rnas = ...
                + rna.counts(rna.matureIndexs, sim.compartment.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs, sim.compartment.cytosolIndexs);
            rnas(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                rnas(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(pc.counts(pc.matureIndexs, :) + pc.counts(pc.boundIndexs, :), 2);
            assertElementsAlmostEqual(rnas' * rna.molecularWeights(rna.aminoacylatedIndexs) / edu.stanford.covert.util.ConstantUtil.nAvogadro, ...
                (1+growth) * mass.cellInitialDryWeight * mass.dryWeightFractionRNA, ...
                'relative', 20e-2, 0);
            
            %RNA expression matches expectations
            assertIn(corr(rna.expression(rna.matureIndexs), rnaExp), [0.95 1]);
            assertIn(180 / pi * acos(rna.expression(rna.matureIndexs)' * rnaExp / ...
                (sqrt(rna.expression(rna.matureIndexs)' * rna.expression(rna.matureIndexs)) * sqrt(rnaExp' * rnaExp))), ...
                [0 15]);
            
            %mRNA expression matches expectations
            assertIn(corr(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs)), rnaExp(rna.matureMRNAIndexs)), [0.85 1])
            assertIn(180 / pi * acos(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))' * rnaExp(rna.matureMRNAIndexs) / ...
                (sqrt(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))' * rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))) * ...
                sqrt(rnaExp(rna.matureMRNAIndexs)' * rnaExp(rna.matureMRNAIndexs)))), ...
                [0 25]);
            
            %tRNAs
            assertElementsAlmostEqual(...
                rna.molecularWeights(rna.matureIndexs(rna.matureTRNAIndexs))' * initRNAs(rna.matureTRNAIndexs) * (1 + growth), ...
                sum(rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), :) + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), :), 2)' * ...
                rna.molecularWeights(rna.matureIndexs(rna.matureTRNAIndexs)),...
                'relative', 50e-2); %todo: tighter
            
            %% assertions -- Protein
            expProd = rna.matureRNAGeneComposition(g.mRNAIndexs, rna.matureMRNAIndexs) * simMrnaExp;
            
            newMonomers = pm.counts - initMonomers;
            newComplexs = pc.counts - initComplexs;
            
            totNewMonomers = sum(...
                + newMonomers(pm.nascentIndexs, :) ...
                + newMonomers(pm.processedIIndexs, :) ...
                + newMonomers(pm.processedIIIndexs, :) ...
                + newMonomers(pm.foldedIndexs, :) ...
                + newMonomers(pm.matureIndexs, :) ...
                + newMonomers(pm.boundIndexs, :) ...
                + newMonomers(pm.misfoldedIndexs, :) ...
                + newMonomers(pm.damagedIndexs, :) ...
                + newMonomers(pm.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.nascentIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.matureIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.boundIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.misfoldedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.damagedIndexs, :), ...
                2);
            
            %amino acids are mostly used -- not too many are free
            assertIn(min(min(aaCnts(:, iter-1000+1:iter))), [0 max(5000, min(initMetabolites(m.aminoAcidIndexs(1:20), comp.cytosolIndexs)))]); %todo: tighter
            assertIn(...
                (m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)' * m.molecularWeights(m.aminoAcidIndexs)) / ...
                (sum(pm.dryWeight + pc.dryWeight + pol.dryWeight) * ConstantUtil.nAvogadro), ...
                [0 0.10]);
            
            %mature/bound monomer counts increased
            assertIn(min(totNewMonomers), [-100 Inf]);
            if ~all(expProd)
                assertIn(max(totNewMonomers(expProd == 0)), [-Inf 0]);
            end
            assertElementsAlmostEqual(...
                totNewMonomers' * monMWs / ConstantUtil.nAvogadro, ...
                growth * mass.cellInitialDryWeight * mass.dryWeightFractionProtein, ...
                'relative', 25e-2, 0);
            assertElementsAlmostEqual(...
                (m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)' * m.molecularWeights(m.aminoAcidIndexs) + ...
                totNewMonomers' * monMWs) / ConstantUtil.nAvogadro, ...
                growth * mass.cellInitialDryWeight * mass.dryWeightFractionProtein, ...
                'relative', 25e-2, 0);
            
            %proteins progress through entire synthesis pathway
            assertIn(sum(sum(pm.counts(pm.nascentIndexs, :))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.processedIIndexs, :))), [0 100]);
            assertIn(sum(sum(pm.counts(pm.processedIIIndexs, :))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.signalSequenceIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.foldedIndexs, :))), [0 200]);
            assertIn(sum(sum(pm.counts(pm.misfoldedIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.damagedIndexs, :))), [0 10]);
            
            assertIn(sum(sum(pc.counts(pc.nascentIndexs, :))), [0 10]);
            assertIn(sum(sum(pc.counts(pc.misfoldedIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pc.counts(pc.damagedIndexs, :))), [0 10]);
            
            assertIn(numel(pol.abortedSequences), [0 1000]);
            
            %RNA production / expression matches expectations
            assertIn(corr(expProd, totNewMonomers), [0.90 1]);
            assertIn(180 / pi * acos((expProd' * totNewMonomers) / (sqrt(expProd' * expProd) * sqrt(totNewMonomers' * totNewMonomers))), ...
                [0 20]);
            
            %protein decay
            assertElementsAlmostEqual((...
                + pm.lengths' * max(0, (pm.decayRates .* sum(initMonomers, 2))) ...
                + pm.lengths(pm.matureIndexs)' * sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * ...
                max(0, pc.decayRates(pc.matureIndexs) .* sum(initComplexs(pc.matureIndexs, :) + initComplexs(pc.boundIndexs, :), 2)) ...
                ) * time.cellCycleLength / log(2) * growth, ...
                sum(freedAAs(:)), ...
                'relative', 7.5e-1);
        end
        
        function testReplicationInitiation(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            %% load fixture
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'Metabolism'
                ...
                'ChromosomeCondensation'
                'ChromosomeSegregation'
                'DNADamage'
                'DNARepair'
                'DNASupercoiling'
                'Replication'
                'ReplicationInitiation'
                'TranscriptionalRegulation'
                ...
                'Transcription'
                'RNAProcessing'
                'RNAModification'
                'tRNAAminoacylation'
                'RNADecay'
                ...
                'RibosomeAssembly'
                'Translation'
                'ProteinProcessingI'
                'ProteinTranslocation'
                'ProteinProcessingII'
                'ProteinFolding'
                'ProteinModification'
                'MacromolecularComplexation'
                'TerminalOrganelleAssembly'
                'ProteinActivation'
                'ProteinDecay'});
            %sim.applyOptions('verbosity', 1);
            
            %% references
            comp = sim.compartment;
            g = sim.gene;
            
            time = sim.state('Time');
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            mr = sim.state('MetabolicReaction');
            pol = sim.state('Polypeptide');
            c = sim.state('Chromosome');
            
            dcy = sim.process('ProteinDecay');
            ta = sim.process('tRNAAminoacylation');
            transcription = sim.process('Transcription');
            cc = sim.process('ChromosomeCondensation');
            repInit = sim.process('ReplicationInitiation');
            rep = sim.process('Replication');
            dr = sim.process('DNARepair');
            tr = sim.process('TranscriptionalRegulation');
            
            %% constants
            monMWs = pm.molecularWeights(pm.matureIndexs);
            
            %% initialize
            mr.growth = 0;
            seed = 1;
            while mr.growth < log(2) / time.cellCycleLength || (repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ATP) + repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_1mer_ATP)) < 54
                seed = seed + 1;
                this.seedRandStream(sim, seed);
                sim.initializeState();
            end
            
            %% keep track of initial state
            initRNAs = ...
                + rna.counts(rna.matureIndexs, sim.compartment.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs, sim.compartment.cytosolIndexs);
            initRNAs(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                initRNAs(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(pc.counts(pc.matureIndexs, :) + pc.counts(pc.boundIndexs, :), 2);
            initMetabolites = m.counts;
            initMonomers = pm.counts;
            initComplexs = pc.counts;
            
            [~, vals] = find(c.monomerBoundSites);
            initBoundMonomers = reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1);
            [~, vals] = find(c.complexBoundSites);
            initBoundComplexs = reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1);
            
            %% simulate
            iterMax = ceil(3.0 * time.replicationInitiationDuration);
            
            rnaExp = zeros(size(initRNAs));
            growth = 0;
            simMrnaExp = zeros(size(rna.matureMRNAIndexs));
            tRNAs = zeros(numel(rna.matureTRNAIndexs), iterMax);
            freedAAs = zeros(numel(dcy.substrateIndexs_aminoAcids), iterMax);
            aaCnts = zeros(numel(ta.substrateIndexs_aminoAcids), iterMax);
            if sim.verbosity > 0
                fprintf('%5s %8s %8s %7s %9s\n', 'Iter ', ' Min AA ', ' Sum AA ', ' DnaA  ', '  R1-5   ');
                fprintf('%5s %8s %8s %7s %9s\n', '=====', '========', '========', '=======', '=========');
            end
            for iter = 1:iterMax
                if sim.verbosity > 0 && mod(iter, 100) == 1
                    [val, idx] = min(m.counts(m.aminoAcidIndexs(1:20), comp.cytosolIndexs));
                    dnaAPol = repInit.calculateDnaAR1234Polymerization();
                    fprintf('%5d %3s %4d %8d %3d %3d %1d %1d %1d %1d %1d\n', iter, ...
                        m.wholeCellModelIDs{m.aminoAcidIndexs(idx)}, val, ...
                        sum(m.counts(m.aminoAcidIndexs(1:20), comp.cytosolIndexs)), ...
                        repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ATP), ...
                        repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_1mer_ATP), ...
                        dnaAPol(1, 1), dnaAPol(2, 1), dnaAPol(3, 1), dnaAPol(4, 1), ...
                        any(repInit.calcuateIsDnaAR5Occupied()));
                end
                
                sim.evolveState();
                
                rnaExp = rnaExp ...
                    + rna.counts(rna.matureIndexs, sim.compartment.cytosolIndexs) ...
                    + rna.counts(rna.aminoacylatedIndexs, sim.compartment.cytosolIndexs);
                rnaExp(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                    rnaExp(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                    sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(initComplexs(pc.matureIndexs, :) + initComplexs(pc.boundIndexs, :), 2);
                
                growth = growth + mr.growth;
                simMrnaExp = simMrnaExp + rna.counts(rna.matureIndexs(rna.matureMRNAIndexs), comp.cytosolIndexs);
                tRNAs(:, iter) = ...
                    + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) ...
                    + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs);
                freedAAs(:, iter) = dcy.substrates(dcy.substrateIndexs_aminoAcids);
                aaCnts(:, iter) = ta.substrates(ta.substrateIndexs_aminoAcids);
                
                this.assertStateIsValid_Protein(sim);
                
                %stop once replication initiated
                if rep.isAnyHelicaseBound && ...
                        ~any(any(repInit.calculateDnaAR1234Polymerization())) && ...
                        ~any(repInit.calcuateIsDnaAR5Occupied())
                    break;
                end
            end
            
            %% assertions - growth
            assertElementsAlmostEqual(log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength), ...
                mr.growth, 'relative', 0.75, 0);
            assertElementsAlmostEqual(exp(log(2) * iter / time.cellCycleLength) - 1, ...
                growth, 'relative', 0.75, 0);
            assertElementsAlmostEqual(mass.cellInitialDryWeight * exp(log(2) * iter / time.cellCycleLength), ...
                sum(mass.cellDry), 'relative', 0.75, 0);
            
            %% assertions -- DNA
            %replication initiated
            assertFalse(rep.isDnaAORIComplexAssembled());
            assertTrue(all(rep.helicasePosition));
            assertTrue(all(rep.leadingPolymerasePosition));
            assertAllEqual(0, repInit.calculateDnaAR1234Polymerization());
            assertAllEqual(false, repInit.calcuateIsDnaAR5Occupied());
            assertIn(collapse(c.polymerizedRegions), [2 * c.sequenceLen + 1 4 * c.sequenceLen]);
            
            %time of replication approximately correct
            assertElementsAlmostEqual(iter, time.replicationInitiationDuration, 'relative', 1);
            
            %chromosome supercoiled
            assertIn(collapse(c.supercoiled), [2 Inf]);
            
            %not too much damage, R/M sites fully methylated
            assertIn(nnz(c.damagedSites), [0 10]);
            
            posStrnds = find(c.restrictableMunIRMSites);
            assertIn(sum(c.isRegionDoubleStranded(posStrnds, 1, false)), [0 1]) %#ok<FNDSB>
            posStrnds = find(c.hemiunmethylatedMunIRMSites);
            assertIn(sum(c.isRegionDoubleStranded(posStrnds, 1, false)), [0 1]) %#ok<FNDSB>
            posStrnds = [
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))      ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2))  2 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)];
            [~, ~, posStrnds] = c.isRegionDoubleStranded(posStrnds, 1, false);
            assertAllEqual(m.getIndexs('m6AD'), c.damagedBases(posStrnds));
            
            %appropriate distribution of bound proteins
            %- gyrase, topoisomerase
            %- DnaA
            %- SMCs
            %- RNA polymerase, transcription factors
            %- repair machinery
            [~, vals] = find(c.monomerBoundSites);
            boundMonomers = reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1);
            [~, vals] = find(c.complexBoundSites);
            boundComplexs = reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1);
            
            assertElementsAlmostEqual(initBoundMonomers, boundMonomers, 'relative', 0.5, 5);
            
            assertElementsAlmostEqual(...
                initBoundComplexs([cc.enzymeComplexGlobalIndexs; tr.enzymeComplexGlobalIndexs]), ...
                boundComplexs([cc.enzymeComplexGlobalIndexs; tr.enzymeComplexGlobalIndexs]), ...
                'relative', 0.5, 4);
            assertElementsAlmostEqual(...
                exp(log(2) * iter / time.replicationInitiationDuration) * ...
                initBoundComplexs(setdiff(1:end, [cc.enzymeComplexGlobalIndexs; tr.enzymeComplexGlobalIndexs; rep.enzymeComplexGlobalIndexs; repInit.enzymeComplexGlobalIndexs; pc.rnaPolymeraseIndexs])), ...
                boundComplexs(setdiff(1:end, [cc.enzymeComplexGlobalIndexs; tr.enzymeComplexGlobalIndexs; rep.enzymeComplexGlobalIndexs; repInit.enzymeComplexGlobalIndexs; pc.rnaPolymeraseIndexs])), ...
                'relative', 1, 4);
            assertElementsAlmostEqual(...
                exp(log(2) * iter / time.replicationInitiationDuration) * ...
                sum(initBoundComplexs(pc.rnaPolymeraseIndexs)), ...
                sum(boundComplexs(pc.rnaPolymeraseIndexs)), ...
                'relative', 1, 4);
            
            %DnaA 7mer broken down and regenerated
            assertElementsAlmostEqual(4 * 7 + 1, ...
                + repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ADP) ...
                + repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_1mer_ADP), ...
                'relative', 0.3, 5);
            assertAllEqual(0, repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ATP(2:end)));
            assertAllEqual(0, repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(2:end)));
            assertAllEqual(0, repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_Nmer_ATP(2:end)));
            assertAllEqual(0, repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(2:end)));
            
            %% assertions -- RNA
            
            %NTPs are limiting
            assertIn(min(transcription.substrates(transcription.substrateIndexs_ntp)), [0 1000]);
            
            %NTPs are mostly used -- not too many are free
            assertIn(...
                (m.counts(m.ntpIndexs([2 4]), sim.compartment.cytosolIndexs)' * m.molecularWeights(m.ntpIndexs([2 4]))) / ...
                (rna.counts(:, sim.compartment.cytosolIndexs)' * rna.molecularWeights), ...
                [0 0.10]);
            
            rnas = ...
                + rna.counts(rna.matureIndexs, sim.compartment.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs, sim.compartment.cytosolIndexs);
            rnas(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                rnas(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(pc.counts(pc.matureIndexs, :) + pc.counts(pc.boundIndexs, :), 2);
            assertElementsAlmostEqual(rnas' * rna.molecularWeights(rna.aminoacylatedIndexs) / edu.stanford.covert.util.ConstantUtil.nAvogadro, ...
                (1+growth) * mass.cellInitialDryWeight * mass.dryWeightFractionRNA, ...
                'relative', 20e-2, 0);
            
            %RNA expression matches expectations
            assertIn(corr(rna.expression(rna.matureIndexs), rnaExp), [0.95 1]);
            assertIn(180 / pi * acos(rna.expression(rna.matureIndexs)' * rnaExp / ...
                (sqrt(rna.expression(rna.matureIndexs)' * rna.expression(rna.matureIndexs)) * sqrt(rnaExp' * rnaExp))), ...
                [0 10]);
            
            %mRNA expression matches expectations
            assertIn(corr(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs)), rnaExp(rna.matureMRNAIndexs)), [0.70 1])
            assertIn(180 / pi * acos(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))' * rnaExp(rna.matureMRNAIndexs) / ...
                (sqrt(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))' * rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))) * ...
                sqrt(rnaExp(rna.matureMRNAIndexs)' * rnaExp(rna.matureMRNAIndexs)))), ...
                [0 40]);
            
            %tRNAs
            assertIn(min(min(tRNAs(:, 1:iter))), [3 Inf]);
            assertElementsAlmostEqual(...
                rna.molecularWeights(rna.matureIndexs(rna.matureTRNAIndexs))' * initRNAs(rna.matureTRNAIndexs) * (1 + growth), ...
                sum(rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), :) + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), :), 2)' * ...
                rna.molecularWeights(rna.matureIndexs(rna.matureTRNAIndexs)),...
                'relative', 30e-2);
            
            %% assertions -- Protein
            expProd = rna.matureRNAGeneComposition(g.mRNAIndexs, rna.matureMRNAIndexs) * simMrnaExp;
            
            newMonomers = pm.counts - initMonomers;
            newComplexs = pc.counts - initComplexs;
            
            totNewMonomers = sum(...
                + newMonomers(pm.nascentIndexs, :) ...
                + newMonomers(pm.processedIIndexs, :) ...
                + newMonomers(pm.processedIIIndexs, :) ...
                + newMonomers(pm.foldedIndexs, :) ...
                + newMonomers(pm.matureIndexs, :) ...
                + newMonomers(pm.boundIndexs, :) ...
                + newMonomers(pm.misfoldedIndexs, :) ...
                + newMonomers(pm.damagedIndexs, :) ...
                + newMonomers(pm.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.nascentIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.matureIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.boundIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.misfoldedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.damagedIndexs, :), ...
                2);
            
            %amino acids are mostly used -- not too many are free
            assertIn(min(min(aaCnts(:, iter-1000+1:iter))), [0 max(5000, min(initMetabolites(m.aminoAcidIndexs(1:20), comp.cytosolIndexs)))]);
            assertIn(...
                (m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)' * m.molecularWeights(m.aminoAcidIndexs)) / ...
                (sum(pm.dryWeight + pc.dryWeight + pol.dryWeight) * ConstantUtil.nAvogadro), ...
                [0 0.10]);
            
            %mature/bound monomer counts increased
            assertIn(min(totNewMonomers), [-100 Inf]);
            if ~all(expProd)
                assertIn(max(totNewMonomers(expProd == 0)), [-Inf 0]);
            end
            assertElementsAlmostEqual(...
                totNewMonomers' * monMWs / ConstantUtil.nAvogadro, ...
                growth * mass.cellInitialDryWeight * mass.dryWeightFractionProtein, ...
                'relative', 20e-2, 0);
            assertElementsAlmostEqual(...
                (m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)' * m.molecularWeights(m.aminoAcidIndexs) + ...
                totNewMonomers' * monMWs) / ConstantUtil.nAvogadro, ...
                growth * mass.cellInitialDryWeight * mass.dryWeightFractionProtein, ...
                'relative', 1e-1, 0);
            
            %proteins progress through entire synthesis pathway
            assertIn(sum(sum(pm.counts(pm.nascentIndexs, :))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.processedIIndexs, :))), [0 100]);
            assertIn(sum(sum(pm.counts(pm.processedIIIndexs, :))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.signalSequenceIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.foldedIndexs, :))), [0 200]);
            assertIn(sum(sum(pm.counts(pm.misfoldedIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.damagedIndexs, :))), [0 10]);
            
            assertIn(sum(sum(pc.counts(pc.nascentIndexs, :))), [0 10]);
            assertIn(sum(sum(pc.counts(pc.misfoldedIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pc.counts(pc.damagedIndexs, :))), [0 10]);
            
            assertIn(numel(pol.abortedSequences), [0 200]);
            
            %RNA production / expression matches expectations
            assertIn(corr(expProd, totNewMonomers), [0.90 1]);
            assertIn(180 / pi * acos((expProd' * totNewMonomers) / (sqrt(expProd' * expProd) * sqrt(totNewMonomers' * totNewMonomers))), ...
                [0 25]);
            
            %protein decay
            assertElementsAlmostEqual((...
                + pm.lengths' * max(0, (pm.decayRates .* sum(initMonomers, 2))) ...
                + pm.lengths(pm.matureIndexs)' * sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * ...
                max(0, pc.decayRates(pc.matureIndexs) .* sum(initComplexs(pc.matureIndexs, :) + initComplexs(pc.boundIndexs, :), 2)) ...
                ) * time.cellCycleLength / log(2) * growth, ...
                sum(freedAAs(:)), ...
                'relative', 5e-1);
        end
        
        %test
        %- replication proceeds correctly in presence of SMCs,
        %  DnaA, RNA polymerase, transcription factors, topoisomerases and
        %  gyrases, damage, repair enzymes
        %- metabolism can support replication
        %- R/M sites get fully methylated
        %- equilibrium superhelicity restored
        %- DNA pol stalls on contact with RNA pol
        %- chromosome segregates
        %- FtsZ rings form and contract cell membrane
        %- cell divides
        %- cell doubles in mass and volume, and DNA, RNA, and protein
        %  content
        function testReplicationInitiationReplicationAndCytokinesis(~)
            %% references
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'ChromosomeCondensation'
                'ChromosomeSegregation'
                'DNADamage'
                'DNARepair'
                'DNASupercoiling'
                'Replication'
                'ReplicationInitiation'
                'Transcription'
                'TranscriptionalRegulation'
                ...
                'ChromosomeSegregation';
                'Cytokinesis';
                'FtsZPolymerization';
                });
            %sim.applyOptions('verbosity', 1);
            %this.seedRandStream(sim, round(mod(now, 1) * 1e7));
            
            comp = sim.compartment;
            
            time = sim.state('Time');
            mass = sim.state('Mass');
            met = sim.state('Metabolite');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            c = sim.state('Chromosome');
            rna = sim.state('Rna');
            trnscpt = sim.state('Transcript');
            rnaPol = sim.state('RNAPolymerase');
            geometry = sim.state('Geometry');
            ring = sim.state('FtsZRing');
            
            repInit = sim.process('ReplicationInitiation');
            rep = sim.process('Replication');
            dr = sim.process('DNARepair');
            ftsZPol = sim.process('FtsZPolymerization');
            
            %% constants
            nCompartments = size(pm.counts, 2);
            
            %% store initial state
            initialWeight = mass.cell;
            initialVolume = geometry.volume;
            initRNAs = rna.counts;
            initRNAs(rna.matureIndexs, :) = ...
                + initRNAs(rna.matureIndexs, :) ...
                + initRNAs(rna.boundIndexs, :);
            initRNAs(rna.boundIndexs, :) = 0;
            initMonomers = pm.counts;
            initMonomers(pm.matureIndexs, :) = ...
                + initMonomers(pm.matureIndexs, :) ...
                + initMonomers(pm.boundIndexs, :);
            initMonomers(pm.boundIndexs, :) = 0;
            initComplexs = pc.counts;
            initComplexs(pc.matureIndexs, :) = ...
                + initComplexs(pc.matureIndexs, :) ...
                + initComplexs(pc.boundIndexs, :);
            initComplexs(pc.boundIndexs, :) = 0;
            
            [~, vals] = find(c.monomerBoundSites);
            initBoundMonomers = reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1);
            [~, vals] = find(c.complexBoundSites);
            initBoundComplexs = reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1);
            
            %% simulate
            iterMax = 1.5 * time.cellCycleLength;
            
            pinchedDiameter = geometry.pinchedDiameter;
            nPinchingCycle = 1;
            
            replicationInitiationFinishedIter = NaN;
            replicationFinishedIter = NaN;
            cytokinesisFinishedIter = NaN;
            
            if sim.verbosity > 0
                fprintf('Iter   DnaA1-5  Helicase Pos   Leading Pos   Lagging Pos       FtsZ Ring      Cytokinesis\n');
                fprintf('===== ========= ============= ============= ============= =================== ===========\n');
            end
            
            for iter = 1:iterMax
                if sim.verbosity > 0 && mod(iter, 100) == 1
                    dnaAPol = repInit.calculateDnaAR1234Polymerization();
                    freeFtsZMonomers = ...
                        + ftsZPol.enzymes(ftsZPol.enzymeIndexs_FtsZ) ...
                        + ftsZPol.enzymes(ftsZPol.enzymeIndexs_FtsZ_GDP) ...
                        + ftsZPol.enzymes(ftsZPol.enzymeIndexs_FtsZ_activated)' * (1:9)';
                    totFtsZMonomers = ...
                        + freeFtsZMonomers ...
                        + ftsZPol.boundEnzymes(ftsZPol.enzymeIndexs_FtsZ) ...
                        + ftsZPol.boundEnzymes(ftsZPol.enzymeIndexs_FtsZ_GDP) ...
                        + ftsZPol.boundEnzymes(ftsZPol.enzymeIndexs_FtsZ_activated)' * (1:9)';
                    fprintf('%5d %1d %1d %1d %1d %1d %6d %6d %6d %6d %6d %6d %3d %3d %2d %2d %2d %2d %2d %.1e\n', iter, ...
                        dnaAPol(1, 1), dnaAPol(2, 1), dnaAPol(3, 1), dnaAPol(4, 1), ...
                        any(repInit.calcuateIsDnaAR5Occupied()), ...
                        rep.helicasePosition(1), rep.helicasePosition(2), ...
                        rep.leadingPosition(1), rep.leadingPosition(2), ...
                        rep.laggingPosition(1), rep.laggingPosition(2), ...
                        freeFtsZMonomers, totFtsZMonomers, ...
                        ring.numEdgesOneStraight, ring.numEdgesTwoStraight,...
                        ring.numEdgesTwoBent, ring.numResidualBent, ...
                        nPinchingCycle, geometry.pinchedDiameter);
                end
                
                %mock metabolism
                met.counts = met.counts + ...
                    log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength) * ...
                    (met.biomassProduction - met.byproducts);
                geometry.calculateVolume();
                
                %mock rna synthesis and decay
                newRNAs = rna.randStream.stochasticRound(...
                    + initRNAs * log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength) ...
                    + initRNAs .* min(1, rna.decayRates(:, ones(nCompartments, 1))) * exp(log(2) * iter / time.cellCycleLength) ...
                    );
                decayedRNAs = rna.randStream.stochasticRound(...
                    min(1, rna.decayRates(:, ones(nCompartments, 1))) .* rna.counts);
                notUpdatingRNAs = rna.updateExternalState(-decayedRNAs, true);
                rna.counts = rna.counts ...
                    + newRNAs ...
                    - (decayedRNAs - notUpdatingRNAs);
                met.counts(met.ntpIndexs, comp.cytosolIndexs) = ...
                    + met.counts(met.ntpIndexs, comp.cytosolIndexs) ...
                    - rna.baseCounts(:, met.nmpIndexs)' * sum(newRNAs, 2);
                met.counts(met.nmpIndexs, comp.cytosolIndexs) = ...
                    + met.counts(met.nmpIndexs, comp.cytosolIndexs) ...
                    + rna.baseCounts(:, met.nmpIndexs)' * sum(decayedRNAs - notUpdatingRNAs, 2);
                
                %mock protein synthesis and decay
                newMonomers = pm.randStream.stochasticRound(...
                    + initMonomers * log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength) ...
                    + initMonomers .* min(1, pm.decayRates(:, ones(nCompartments, 1))) * exp(log(2) * iter / time.cellCycleLength) ...
                    );
                newComplexs = pc.randStream.stochasticRound(...
                    + initComplexs * log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength) ...
                    + initComplexs .* min(1, pc.decayRates(:, ones(nCompartments, 1))) * exp(log(2) * iter / time.cellCycleLength) ...
                    );
                decayedMonomers = pm.randStream.stochasticRound(min(1, pm.decayRates(:, ones(nCompartments, 1))) .* pm.counts);
                decayedComplexs = pc.randStream.stochasticRound(min(1, pc.decayRates(:, ones(nCompartments, 1))) .* pc.counts);
                notUpdatingMonomers = pm.updateExternalState(-decayedMonomers, true);
                notUpdatingComplexs = pc.updateExternalState(-decayedComplexs, true);
                pm.counts = pm.counts ...
                    + newMonomers ...
                    - (decayedMonomers - notUpdatingMonomers);
                pc.counts = pc.counts ...
                    + newComplexs ...
                    - (decayedComplexs - notUpdatingComplexs);
                met.counts(met.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    met.counts(met.aminoAcidIndexs, comp.cytosolIndexs) ...
                    - pm.baseCounts(:, met.aminoAcidIndexs)' * sum(newMonomers - (decayedMonomers - notUpdatingMonomers), 2) ...
                    - pc.baseCounts(:, met.aminoAcidIndexs)' * sum(newComplexs - (decayedComplexs - notUpdatingComplexs), 2);
                
                %simulate
                sim.evolveState();
                
                %note when replication initiated
                if isnan(replicationInitiationFinishedIter) && ...
                        rep.isAnyHelicaseBound && ...
                        ~any(any(repInit.calculateDnaAR1234Polymerization())) && ...
                        ~any(repInit.calcuateIsDnaAR5Occupied())
                    replicationInitiationFinishedIter = iter;
                end
                
                %note when replication finished
                if isnan(replicationFinishedIter) && all(rep.strandDuplicated)
                    replicationFinishedIter = iter;
                end
                
                %note when pinched diameter contracts
                if pinchedDiameter ~= geometry.pinchedDiameter
                    nPinchingCycle = nPinchingCycle + 1;
                    pinchedDiameter = geometry.pinchedDiameter;
                end
                
                %stop once divided
                if c.segregated && geometry.pinched && ring.numResidualBent == 0 && geometry.width == -1
                    cytokinesisFinishedIter = iter;
                    break;
                end
            end
            
            %% assert timing
            %time of replication initiation approximately correct
            assertElementsAlmostEqual(time.replicationInitiationDuration, ...
                replicationInitiationFinishedIter, ...
                'relative', 0.8);
            
            %time of replication initiation + replication approximately correct
            assertElementsAlmostEqual(time.replicationInitiationDuration + time.replicationDuration, ...
                replicationFinishedIter, ...
                'relative', 0.8);
            
            %time for cytokinesis approximately correct
            assertElementsAlmostEqual(time.cytokinesisDuration, ...
                cytokinesisFinishedIter - replicationFinishedIter, ...
                'relative', 0.7);
            
            %cell cycle time approximately correct
            assertElementsAlmostEqual(time.cellCycleLength, ...
                cytokinesisFinishedIter, ...
                'relative', 0.7);
            
            %% assert mass, volume, geometry
            %weight, volume doubled
            assertElementsAlmostEqual(sum(2 * initialWeight), sum(mass.cell), 'relative', 0.6, 0); %inaccurate because above "simulation" is not very accurate
            assertElementsAlmostEqual(2 * initialVolume, geometry.volume, 'relative', 0.6, 0); %inaccurate because above "simulation" is not very accurate
            
            %completely pinched
            assertEqual(0, geometry.pinchedDiameter);
            assertTrue(geometry.pinched);
            
            %cell shape undefined, completely pinched
            assertEqual(0, geometry.pinchedDiameter);
            assertEqual(-1, geometry.width);
            assertEqual(-1, geometry.cylindricalLength);
            assertEqual(-1, geometry.surfaceArea);
            assertEqual(-1, geometry.totalLength);
            
            %% assert chromosome duplicated, segregated, supercoiled, undamaged, and appropriately bound by proteins
            %replicated
            assertTrue(all(rep.strandDuplicated));
            
            %chromosome segregated
            assertTrue(c.segregated);
            
            %equilibrium superhelical density restored
            assertEqual([ones(4, 1) (1:4)'], find(c.linkingNumbers));
            assertTrue(all(full(c.supercoiled(1, :))));
            
            %not too much damage
            assertIn(nnz(c.damagedSites), [0 10]);
            
            %R/M sites fully methylated
            assertIn(nnz(c.restrictableMunIRMSites), [0 1]);
            assertIn(nnz(c.hemiunmethylatedMunIRMSites), [0 100]);
            assertEqual([0; met.m6ADIndexs], unique([0; c.damagedBases([
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))      ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2))  2 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))  3 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2))  4 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)])]));
            
            %DNA polymerase stalled, and RNA polymerase aborted resulting
            %in incomplete transcripts
            nActRNAPol = rnaPol.stateExpectations(rnaPol.activelyTranscribingIndex) * (...
                + pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE')), comp.cytosolIndexs) ...
                + pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')), comp.cytosolIndexs) ...
                + pc.counts(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE')), comp.cytosolIndexs) ...
                + pc.counts(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')), comp.cytosolIndexs));
            assertElementsAlmostEqual(nActRNAPol, ...
                numel(trnscpt.abortedSequences), ...
                'relative', 1);
            assertElementsAlmostEqual(iter, time.replicationInitiationDuration + time.replicationDuration, ...
                'relative', 0.8);
            
            %distribution of bound proteins restored
            %- for most proteins this means similar number bound as before
            %- except for SMC which should have the bound density preserved
            [~, vals] = find(c.monomerBoundSites);
            boundMonomers = reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1);
            [~, vals] = find(c.complexBoundSites);
            boundComplexs = reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1);
            
            assertElementsAlmostEqual(...
                exp(log(2) * iter / time.replicationInitiationDuration) * initBoundMonomers, ...
                boundMonomers, ...
                'relative', 0.75, 4);
            assertElementsAlmostEqual(...
                exp(log(2) * iter / time.replicationInitiationDuration) * ...
                initBoundComplexs(setdiff(1:end, [repInit.enzymeComplexGlobalIndexs; pc.rnaPolymeraseIndexs])), ...
                boundComplexs(setdiff(1:end, [repInit.enzymeComplexGlobalIndexs; pc.rnaPolymeraseIndexs])), ...
                'relative', 0.75, 4);
            assertElementsAlmostEqual(...
                exp(log(2) * iter / time.replicationInitiationDuration) * ...
                sum(initBoundComplexs(pc.rnaPolymeraseIndexs)), ...
                sum(boundComplexs(pc.rnaPolymeraseIndexs)), ...
                'relative', 0.75, 4);
            
            %% assert DnaA state
            %DnaA 7mer broken down and regenerated, and not polymerized
            %(too much)
            assertElementsAlmostEqual(zeros(4, 2), repInit.calculateDnaAR1234Polymerization(), 'absolute', 7);            
            
            assertElementsAlmostEqual(0.01 * (4 * 7 + 1), ...
                + repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ADP) ...
                + repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_1mer_ADP), ...
                'relative', 0.2, 3);
            assertAllEqual(0, repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ATP(2:end)));
            assertAllEqual(0, repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(2:end)));
            assertAllEqual(0, repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(2:end)));
        end
        
        %test
        %- replication proceeds correctly in presence of SMCs,
        %  DnaA, RNA polymerase, transcription factors, topoisomerases and
        %  gyrases, damage, repair enzymes
        %- metabolism can support replication
        %- R/M sites get fully methylated
        %- DNA damage doesn't accumulate
        %- equilibrium superhelicity restored
        %- DNA pol stalls on contact with RNA pol
        %- distribution of chromosomally bound proteins is similar at
        %  beginning and end of cell cycle
        %- chromosome segregates
        %- FtsZ rings form and contract cell membrane
        %- cell divides, and geometric properties reflect this
        %- cell doubles in mass and volume, and DNA, RNA, and protein
        %  content
        %- growth rate doubles
        %- timing of replication initiation, replication, and cytokinesis
        %- NTPs, AAs have all been used for transcription, translation
        %- RNA, protein expression matches expectations
        %- no accumulation of immature proteins, RNA
        %- protein decay follows expectations
        function testCellGrowthAndDivision(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            %% references
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'Metabolism'
                ...
                'ChromosomeCondensation'
                'ChromosomeSegregation'
                'DNADamage'
                'DNARepair'
                'DNASupercoiling'
                'Replication'
                'ReplicationInitiation'
                'TranscriptionalRegulation'
                ...
                'Cytokinesis';
                'FtsZPolymerization';
                ...
                'Transcription'
                'RNAProcessing'
                'RNAModification'
                'tRNAAminoacylation'
                'RNADecay'
                ...
                'RibosomeAssembly'
                'Translation'
                'ProteinProcessingI'
                'ProteinTranslocation'
                'ProteinProcessingII'
                'ProteinFolding'
                'ProteinModification'
                'MacromolecularComplexation'
                'TerminalOrganelleAssembly'
                'ProteinActivation'
                'ProteinDecay'
                ...
                'HostInteraction'
                });
            %sim.applyOptions('verbosity', 1);
            %this.seedRandStream(sim, round(mod(now, 1) * 1e7));
            
            comp = sim.compartment;
            g = sim.gene;
            
            time = sim.state('Time');
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            c = sim.state('Chromosome');
            rna = sim.state('Rna');
            trnscpt = sim.state('Transcript');
            rnaPol = sim.state('RNAPolymerase');
            geometry = sim.state('Geometry');
            ring = sim.state('FtsZRing');
            mr = sim.state('MetabolicReaction');
            rib = sim.state('Ribosome');
            pol = sim.state('Polypeptide');
            
            repInit = sim.process('ReplicationInitiation');
            rep = sim.process('Replication');
            dr = sim.process('DNARepair');
            ftsZPol = sim.process('FtsZPolymerization');
            dcy = sim.process('ProteinDecay');
            ta = sim.process('tRNAAminoacylation');
            transcription = sim.process('Transcription');
            
            fthf10Indexs = m.getIndexs({'FMET'; 'FTHF10'; 'FOR'; 'THF'});
            fmetTRNAIndexs = rna.aminoacylatedIndexs(rna.baseCounts(rna.aminoacylatedIndexs, m.aminoAcidIndexs(21)) ~= 0);
            
            %% initialize
            sim.initializeState();
            
            %% store initial state
            initialGrowth = mr.growth;
            initialWeight = mass.cell;
            initialVolume = geometry.volume;
            initMetabolites = m.counts;
            initRNAs = rna.counts;
            initRNAs(rna.matureIndexs, :) = ...
                + initRNAs(rna.matureIndexs, :) ...
                + initRNAs(rna.boundIndexs, :) ...
                + initRNAs(rna.aminoacylatedIndexs, :);
            initRNAs(rna.boundIndexs, :) = 0;
            initRNAs(rna.aminoacylatedIndexs, :) = 0;
            initMonomers = pm.counts;
            initMonomers(pm.matureIndexs, :) = ...
                + initMonomers(pm.matureIndexs, :) ...
                + initMonomers(pm.boundIndexs, :);
            initMonomers(pm.boundIndexs, :) = 0;
            initComplexs = pc.counts;
            initComplexs(pc.matureIndexs, :) = ...
                + initComplexs(pc.matureIndexs, :) ...
                + initComplexs(pc.boundIndexs, :);
            initComplexs(pc.boundIndexs, :) = 0;
            
            [~, vals] = find(c.monomerBoundSites);
            initBoundMonomers = reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1);
            [~, vals] = find(c.complexBoundSites);
            initBoundComplexs = reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1);
            
            %% simulate
            iterMax = 2.0 * time.cellCycleLength;
            
            pinchedDiameter = geometry.pinchedDiameter;
            nPinchingCycle = 1;
            
            replicationInitiationFinishedIter = NaN;
            replicationFinishedIter = NaN;
            cytokinesisFinishedIter = NaN;
            
            growth = 0;
            rnaExp = zeros(size(rna.counts));
            pmExp = zeros(size(pm.counts));
            pcExp = zeros(size(pc.counts));
            freedAAs = zeros(numel(dcy.substrateIndexs_aminoAcids), iterMax);
            ntpCnts = zeros(4, iterMax);
            ndpCnts = zeros(4, iterMax);
            nmpCnts = zeros(4, iterMax);
            rnaCnts = zeros(1, iterMax);
            aaCnts = zeros(numel(ta.substrateIndexs_aminoAcids), iterMax);
            fthf10Cnts = zeros(8, iterMax);
            requirements = zeros(numel(m.counts), numel(sim.processes));
            allocations = zeros(numel(m.counts), numel(sim.processes));
            usages = zeros(numel(m.counts), numel(sim.processes));
            rnaPolOcc = zeros(numel(rnaPol.stateOccupancies), iterMax);
            
            ticker = tic;
            
            if sim.verbosity > 0
                fprintf('Iter  RTime Growth  Mass   RNA   Mon   Cpx   DnaA1-5  Helicase Pos   Leading Pos   Lagging Pos       FtsZ Ring      Cytokinesis\n');
                fprintf('===== ===== ====== ====== ===== ===== ===== ========= ============= ============= ============= =================== ===========\n');
            end
            
            for iter = 1:iterMax
                if sim.verbosity > 0 && mod(iter, 100) == 1
                    dnaAPol = repInit.calculateDnaAR1234Polymerization();
                    freeFtsZMonomers = ...
                        + ftsZPol.enzymes(ftsZPol.enzymeIndexs_FtsZ) ...
                        + ftsZPol.enzymes(ftsZPol.enzymeIndexs_FtsZ_GDP) ...
                        + ftsZPol.enzymes(ftsZPol.enzymeIndexs_FtsZ_activated)' * (1:9)';
                    totFtsZMonomers = ...
                        + freeFtsZMonomers ...
                        + ftsZPol.boundEnzymes(ftsZPol.enzymeIndexs_FtsZ) ...
                        + ftsZPol.boundEnzymes(ftsZPol.enzymeIndexs_FtsZ_GDP) ...
                        + ftsZPol.boundEnzymes(ftsZPol.enzymeIndexs_FtsZ_activated)' * (1:9)';
                    fprintf('%5d %5.1f %.4f %.4f %5d %5d %5d %1d %1d %1d %1d %1d %6d %6d %6d %6d %6d %6d %3d %3d %2d %2d %2d %2d %2d %.1e\n', ...
                        iter, toc(ticker), ...
                        mr.growth / initialGrowth, sum(mass.cellDry) / mass.cellInitialDryWeight, ...
                        sum(rna.counts(:)), sum(pm.counts(:)), sum(pc.counts(:)), ...
                        dnaAPol(1, 1), dnaAPol(2, 1), dnaAPol(3, 1), dnaAPol(4, 1), ...
                        any(repInit.calcuateIsDnaAR5Occupied()), ...
                        rep.helicasePosition(1), rep.helicasePosition(2), ...
                        rep.leadingPosition(1), rep.leadingPosition(2), ...
                        rep.laggingPosition(1), rep.laggingPosition(2), ...
                        freeFtsZMonomers, totFtsZMonomers, ...
                        ring.numEdgesOneStraight, ring.numEdgesTwoStraight,...
                        ring.numEdgesTwoBent, ring.numResidualBent, ...
                        nPinchingCycle, geometry.pinchedDiameter);
                    ticker = tic;
                end
                time.values = iter;
                
                %simulate
                [~, tmp1, tmp2, tmp3] = sim.evolveState();
                
                %note metabolite requirements, allocations, usage
                requirements = requirements + tmp1;
                allocations = allocations + tmp2;
                usages = usages + tmp3;
                
                %error check state
                for i = 1:numel(sim.processes)
                    p = sim.processes{i};
                    
                    %metabolites
                    tfs = p.substrates < 0 | ceil(p.substrates) ~= p.substrates;
                    tfs(p.substrateStimulusLocalIndexs, :) = false;
                    if any(tfs)
                        [idxs, ~] = find(tfs);
                        tmp = cellfun(@(x, y) sprintf('%s\t%f', x, y), p.substrateWholeCellModelIDs(idxs), num2cell(p.substrates(tfs)), 'UniformOutput', false);
                        throw(MException('Simulation_Large_Test:invalidCounts', ...
                            'Invalid substrates in process %s at iter %d\n- %s', ...
                            p.name, iter, strjoin(sprintf('\n- '), tmp{:})));
                    end
                    
                    %enzymes
                    tfs = p.enzymes(:) < 0 | ceil(p.enzymes(:)) ~= p.enzymes(:);
                    if any(tfs)
                        tmp = cellfun(@(x, y) sprintf('%s\t%f', x, y), p.enzymeWholeCellModelIDs(tfs), num2cell(p.enzymes(tfs)), 'UniformOutput', false);
                        throw(MException('Simulation_Large_Test:invalidCounts', ...
                            'Invalid enzymes in process %s at iter %d\n- %s', ...
                            p.name, iter, strjoin(sprintf('\n- '), tmp{:})));
                    end
                    
                    %bound enzymes
                    tfs = p.boundEnzymes(:) < 0 | ceil(p.boundEnzymes(:)) ~= p.boundEnzymes(:);
                    if any(tfs)
                        tmp = cellfun(@(x, y) sprintf('%s\t%f', x, y), p.enzymeWholeCellModelIDs(tfs), num2cell(p.boundEnzymes(tfs)), 'UniformOutput', false);
                        throw(MException('Simulation_Large_Test:invalidCounts', ...
                            'Invalid bound enzymes in process %s at iter %d\n- %s', ...
                            p.name, iter, strjoin(sprintf('\n- '), tmp{:})));
                    end
                end
                
                tfs = m.counts < 0 | ceil(m.counts) ~= m.counts;
                if any(tfs(:))
                    [i, j] = find(tfs);
                    tmp = cellfun(@(x, y) sprintf('%s\t%s\t%f', x, y), m.wholeCellModelIDs(i), comp.wholeCellModelIDs(j), num2cell(m.counts(tfs)), 'UniformOutput', false);
                    throw(MException('Simulation_Large_Test:invalidCounts', ...
                        'Invalid metabolites at iter %d\n- %s', ...
                        iter, strjoin(sprintf('\n- '), tmp{:})));
                end
                
                tfs = rna.counts(:) < 0 | ceil(rna.counts(:)) ~= rna.counts(:);
                if any(tfs)
                    [idxs, ~] = ind2sub(size(rna.counts), find(tfs));
                    tmp = cellfun(@(x, y) sprintf('%s\t%f', x, y), rna.wholeCellModelIDs(idxs), num2cell(rna.counts(tfs)), 'UniformOutput', false);
                    throw(MException('Simulation_Large_Test:invalidCounts', ...
                        'Invalid RNAs at iter %d\n- %s', ...
                        iter, strjoin(sprintf('\n- '), tmp{:})));
                end
                
                tfs = pm.counts(:) < 0 | ceil(pm.counts(:)) ~= pm.counts(:);
                if any(tfs)
                    [idxs, ~] = ind2sub(size(pm.counts), find(tfs));
                    tmp = cellfun(@(x, y) sprintf('%s\t%f', x, y), pm.wholeCellModelIDs(idxs), num2cell(pm.counts(tfs)), 'UniformOutput', false);
                    throw(MException('Simulation_Large_Test:invalidCounts', ...
                        'Invalid monomers at iter %d\n- %s', ...
                        iter, strjoin(sprintf('\n- '), tmp{:})));
                end
                
                tfs = pc.counts(:) < 0 | ceil(pc.counts(:)) ~= pc.counts(:);
                if any(tfs)
                    [idxs, ~] = ind2sub(size(pc.counts), find(tfs));
                    tmp = cellfun(@(x, y) sprintf('%s\t%f', x, y), pc.wholeCellModelIDs(idxs), num2cell(pc.counts(tfs)), 'UniformOutput', false);
                    throw(MException('Simulation_Large_Test:invalidCounts', ...
                        'Invalid complexes at iter %d\n- %s', ...
                        iter, strjoin(sprintf('\n- '), tmp{:})));
                end
                
                %note RNA expression
                growth = growth + mr.growth;
                rnaExp = rnaExp + rna.counts;
                pmExp = pmExp + pm.counts;
                pcExp = pcExp + pc.counts;
                freedAAs(:, iter) = dcy.substrates(dcy.substrateIndexs_aminoAcids);
                rnaCnts(:, iter) = sum(rna.counts(:));
                aaCnts(:, iter) = ta.substrates(ta.substrateIndexs_aminoAcids);
                ntpCnts(:, iter) = m.counts(m.ntpIndexs, comp.cytosolIndexs);
                ndpCnts(:, iter) = m.counts(m.ndpIndexs, comp.cytosolIndexs);
                nmpCnts(:, iter) = m.counts(m.nmpIndexs, comp.cytosolIndexs);
                fthf10Cnts(:, iter) = [
                    m.counts(fthf10Indexs, comp.cytosolIndexs) %[FMET FTHF10 FOR THF]
                    sum(pm.counts(pm.nascentIndexs, comp.cytosolIndexs))
                    rna.counts(fmetTRNAIndexs, comp.cytosolIndexs)
                    sum(rib.mRNAPositions > 0 | rib.tmRNAPositions > 0)
                    size(pol.abortedPolypeptides, 1)
                    ];
                rnaPolOcc(:, iter) = rnaPol.stateOccupancies;
                
                %note when replication initiated
                if isnan(replicationInitiationFinishedIter) && ...
                        rep.isAnyHelicaseBound && ...
                        ~any(any(repInit.calculateDnaAR1234Polymerization())) && ...
                        ~any(repInit.calcuateIsDnaAR5Occupied())
                    replicationInitiationFinishedIter = iter;
                end
                
                %note when replication finished
                if isnan(replicationFinishedIter) && all(rep.strandDuplicated)
                    replicationFinishedIter = iter;
                end
                
                %note when pinched diameter contracts
                if pinchedDiameter ~= geometry.pinchedDiameter
                    nPinchingCycle = nPinchingCycle + 1;
                    pinchedDiameter = geometry.pinchedDiameter;
                end
                
                %stop once divided
                if c.segregated && geometry.pinched && ring.numResidualBent == 0 && geometry.width == -1
                    cytokinesisFinishedIter = iter;
                    break;
                end
            end
            
            %% assert timing
            %time of replication initiation approximately correct
            assertElementsAlmostEqual(time.replicationInitiationDuration, ...
                replicationInitiationFinishedIter, ...
                'relative', 0.8);
            
            %time of replication initiation + replication approximately correct
            assertElementsAlmostEqual(time.replicationInitiationDuration + time.replicationDuration, ...
                replicationFinishedIter, ...
                'relative', 0.8);
            
            %time for cytokinesis approximately correct
            assertElementsAlmostEqual(time.cytokinesisDuration, ...
                cytokinesisFinishedIter - replicationFinishedIter, ...
                'relative', 0.7);
            
            %cell cycle time approximately correct
            assertElementsAlmostEqual(time.cellCycleLength, ...
                cytokinesisFinishedIter, ...
                'relative', 0.7);
            
            %% assertions - growth
            assertElementsAlmostEqual(log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength), ...
                mr.growth, 'relative', 0.75, 0);
            assertElementsAlmostEqual(exp(log(2) * iter / time.cellCycleLength) - 1, ...
                growth, 'relative', 0.75, 0);
            assertElementsAlmostEqual(mass.cellInitialDryWeight * exp(log(2) * iter / time.cellCycleLength), ...
                sum(mass.cellDry), 'relative', 0.75, 0);
            
            %% assert mass, volume, geometry
            %weight, volume doubled
            assertElementsAlmostEqual(sum(2 * initialWeight), sum(mass.cell), 'relative', 0.5, 0);
            assertElementsAlmostEqual(2 * initialVolume, geometry.volume, 'relative', 0.5, 0);
            
            %completely pinched
            assertEqual(0, geometry.pinchedDiameter);
            assertTrue(geometry.pinched);
            
            %cell shape undefined, completely pinched
            assertEqual(0, geometry.pinchedDiameter);
            assertEqual(-1, geometry.width);
            assertEqual(-1, geometry.cylindricalLength);
            assertEqual(-1, geometry.surfaceArea);
            assertEqual(-1, geometry.totalLength);
            
            %% assert chromosome duplicated, segregated, supercoiled, undamaged, and appropriately bound by proteins
            %replicated
            assertTrue(all(rep.strandDuplicated));
            
            %chromosome segregated
            assertTrue(c.segregated);
            
            %equilibrium superhelical density restored
            assertEqual([ones(4, 1) (1:4)'], find(c.linkingNumbers));
            assertTrue(all(full(c.supercoiled(1, :))));
            
            %not too much damage
            assertIn(nnz(c.damagedSites), [0 10]);
            
            %R/M sites fully methylated
            assertIn(nnz(c.restrictableMunIRMSites), [0 1]);
            assertIn(nnz(c.hemiunmethylatedMunIRMSites), [0 1]);
            assertAllEqual(m.m6ADIndexs, c.damagedBases([
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))      ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2))  2 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))  3 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2))  4 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)]));
            
            %DNA polymerase stalled, and RNA polymerase aborted resulting
            %in incomplete transcripts
            nActRNAPol = rnaPol.stateExpectations(rnaPol.activelyTranscribingIndex) * (...
                + pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE')), comp.cytosolIndexs) ...
                + pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')), comp.cytosolIndexs) ...
                + pc.counts(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE')), comp.cytosolIndexs) ...
                + pc.counts(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')), comp.cytosolIndexs));
            assertElementsAlmostEqual(nActRNAPol, ...
                numel(trnscpt.abortedSequences), ...
                'relative', 1);
            assertElementsAlmostEqual(iter, time.replicationInitiationDuration + time.replicationDuration, ...
                'relative', 0.8);
            
            %distribution of bound proteins restored
            %- for most proteins this means similar number bound as before
            %- except for SMC which should have the bound density preserved
            [~, vals] = find(c.monomerBoundSites);
            boundMonomers = reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1);
            [~, vals] = find(c.complexBoundSites);
            boundComplexs = reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1);
            
            assertElementsAlmostEqual(...
                exp(log(2) * iter / time.replicationInitiationDuration) * initBoundMonomers, ...
                boundMonomers, ...
                'relative', 0.75, 4);
            assertElementsAlmostEqual(...
                exp(log(2) * iter / time.replicationInitiationDuration) * ...
                initBoundComplexs(setdiff(1:end, [repInit.enzymeComplexGlobalIndexs; pc.rnaPolymeraseIndexs])), ...
                boundComplexs(setdiff(1:end, [repInit.enzymeComplexGlobalIndexs; pc.rnaPolymeraseIndexs])), ...
                'relative', 0.90, 4);
            assertElementsAlmostEqual(...
                exp(log(2) * iter / time.replicationInitiationDuration) * ...
                sum(initBoundComplexs(pc.rnaPolymeraseIndexs)), ...
                sum(boundComplexs(pc.rnaPolymeraseIndexs)), ...
                'relative', 0.75, 4);
            
            %% assert DnaA state
            %DnaA 7mer broken down and regenerated, and not polymerized
            %(too much)
            assertElementsAlmostEqual(zeros(4, 2), repInit.calculateDnaAR1234Polymerization(), 'absolute', 7);
            
            assertElementsAlmostEqual(.45^((cytokinesisFinishedIter - replicationInitiationFinishedIter)/3600) * (4 * 7 + 1), ...
                + repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ADP) ...
                + repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_1mer_ADP), ...
                'relative', 0.2, 2);
            assertAllEqual(0, repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ATP(2:end)));
            assertAllEqual(0, repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(2:end)));
            assertAllEqual(0, repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(2:end)));
            
            %% assertions -- RNA
            
            %NTPs are limiting
            assertIn(min(transcription.substrates(transcription.substrateIndexs_ntp)), [0 6000]);
            
            %NTPs are mostly used -- not too many are free
            assertIn(...
                (m.counts(m.ntpIndexs([2 4]), sim.compartment.cytosolIndexs)' * m.molecularWeights(m.ntpIndexs([2 4]))) / ...
                (rna.counts(:, sim.compartment.cytosolIndexs)' * rna.molecularWeights), ...
                [0 0.10]);
            
            rnas = ...
                + rna.counts(rna.matureIndexs, sim.compartment.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs, sim.compartment.cytosolIndexs);
            rnas(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                rnas(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(pc.counts(pc.matureIndexs, :) + pc.counts(pc.boundIndexs, :), 2);
            assertElementsAlmostEqual(rnas' * rna.molecularWeights(rna.aminoacylatedIndexs) / edu.stanford.covert.util.ConstantUtil.nAvogadro, ...
                (1 + growth) * mass.cellInitialDryWeight * mass.dryWeightFractionRNA, ...
                'relative', 3e-1, 0);
            
            %RNA expression matches expectations
            matureRNAExpression = ...
                + rnaExp(rna.matureIndexs, comp.cytosolIndexs) ...
                + rnaExp(rna.aminoacylatedIndexs, comp.cytosolIndexs) ...
                + rnaExp(rna.boundIndexs, comp.cytosolIndexs);
            matureRNAExpression(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                + matureRNAExpression(setdiff(1:end, rna.matureMRNAIndexs)) ...
                + sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * sum(...
                + pcExp(pc.nascentIndexs, :) ...
                + pcExp(pc.matureIndexs, :) ...
                + pcExp(pc.boundIndexs, :), 2);
            assertIn(corr(rna.expression(rna.matureIndexs), matureRNAExpression), [0.90 1]);
            assertIn(180 / pi * acos(rna.expression(rna.matureIndexs)' * matureRNAExpression / ...
                (sqrt(rna.expression(rna.matureIndexs)' * rna.expression(rna.matureIndexs)) * ...
                sqrt(matureRNAExpression' * matureRNAExpression))), ...
                [0 20]);
            
            %mRNA expression matches expectations
            mRNAExp = matureRNAExpression(rna.matureMRNAIndexs);
            assertIn(corr(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs)), mRNAExp), [0.90 1]);
            assertIn(180 / pi * acos(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))' * mRNAExp / ...
                (sqrt(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))' * rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))) * ...
                sqrt(mRNAExp' * mRNAExp))), ...
                [0 25]);
            
            %tRNAs
            tRNAExp = matureRNAExpression(rna.matureTRNAIndexs);
            assertIn(min(tRNAExp(:)), [3 Inf]);
            assertElementsAlmostEqual(...
                rna.molecularWeights(rna.matureIndexs(rna.matureTRNAIndexs))' * initRNAs(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) * (1 + growth), ...
                sum(rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), :) + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), :), 2)' * ...
                rna.molecularWeights(rna.matureIndexs(rna.matureTRNAIndexs)),...
                'relative', 60e-2);
            
            %% assertions -- Protein
            expProd = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * matureRNAExpression;
            
            newMonomers = pm.counts - initMonomers;
            newComplexs = pc.counts - initComplexs;
            
            totNewMonomers = sum(...
                + newMonomers(pm.nascentIndexs, :) ...
                + newMonomers(pm.processedIIndexs, :) ...
                + newMonomers(pm.processedIIIndexs, :) ...
                + newMonomers(pm.foldedIndexs, :) ...
                + newMonomers(pm.matureIndexs, :) ...
                + newMonomers(pm.boundIndexs, :) ...
                + newMonomers(pm.misfoldedIndexs, :) ...
                + newMonomers(pm.damagedIndexs, :) ...
                + newMonomers(pm.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.nascentIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.matureIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.boundIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.misfoldedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.damagedIndexs, :), ...
                2);
            totMonExp = sum(...
                + pmExp(pm.nascentIndexs, :) ...
                + pmExp(pm.processedIIndexs, :) ...
                + pmExp(pm.processedIIIndexs, :) ...
                + pmExp(pm.foldedIndexs, :) ...
                + pmExp(pm.matureIndexs, :) ...
                + pmExp(pm.boundIndexs, :) ...
                + pmExp(pm.misfoldedIndexs, :) ...
                + pmExp(pm.damagedIndexs, :) ...
                + pmExp(pm.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * pcExp(pc.nascentIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * pcExp(pc.matureIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * pcExp(pc.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * pcExp(pc.boundIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * pcExp(pc.misfoldedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * pcExp(pc.damagedIndexs, :), ...
                2);
            
            %amino acids are mostly used -- not too many are free
            assertIn(min(min(aaCnts(:, iter-1000+1:iter))), [0 max(4e4, min(initMetabolites(m.aminoAcidIndexs(1:20), comp.cytosolIndexs)))]);
            assertIn(...
                (m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)' * m.molecularWeights(m.aminoAcidIndexs)) / ...
                (sum(pm.dryWeight + pc.dryWeight + pol.dryWeight) * ConstantUtil.nAvogadro), ...
                [0.0 0.10]);
            
            %mature/bound monomer counts increased
            assertIn(min(totNewMonomers), [-100 Inf]);
            assertElementsAlmostEqual(...
                sum(newMonomers' * pm.molecularWeights + newComplexs' * pc.molecularWeights) / ConstantUtil.nAvogadro, ...
                growth * mass.cellInitialDryWeight * mass.dryWeightFractionProtein, ...
                'relative', 2e-1, 0);
            assertElementsAlmostEqual(...
                (m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)' * m.molecularWeights(m.aminoAcidIndexs) + ...
                sum(newMonomers' * pm.molecularWeights + newComplexs' * pc.molecularWeights)) / ConstantUtil.nAvogadro, ...
                growth * mass.cellInitialDryWeight * mass.dryWeightFractionProtein, ...
                'relative', 2e-1, 0);
            
            %proteins progress through entire synthesis pathway
            assertIn(sum(sum(pm.counts(pm.nascentIndexs, :))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.processedIIndexs, :))), [0 150]);
            assertIn(sum(sum(pm.counts(pm.processedIIIndexs, :))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.signalSequenceIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.foldedIndexs, :))), [0 200]);
            assertIn(sum(sum(pm.counts(pm.misfoldedIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.damagedIndexs, :))), [0 10]);
            
            assertIn(sum(sum(pc.counts(pc.nascentIndexs, :))), [0 10]);
            assertIn(sum(sum(pc.counts(pc.misfoldedIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pc.counts(pc.damagedIndexs, :))), [0 10]);
            
            assertIn(numel(pol.abortedSequences), [0 2000]);
            
            %monomer production / mRNA expression matches expectations
            assertIn(corr(expProd, totMonExp), [0.90 1]);
            assertIn(180 / pi * acos((expProd' * totMonExp) / (sqrt(expProd' * expProd) * sqrt(totMonExp' * totMonExp))), ...
                [0 20]);
            
            %protein decay
            assertElementsAlmostEqual((...
                + pm.lengths' * max(0, (pm.decayRates .* sum(initMonomers, 2))) ...
                + pm.lengths(pm.matureIndexs)' * sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * ...
                max(0, pc.decayRates(pc.matureIndexs) .* sum(initComplexs(pc.matureIndexs, :) + initComplexs(pc.boundIndexs, :), 2)) ...
                ) * time.cellCycleLength / log(2) * growth, ...
                sum(freedAAs(:)), ...
                'relative', 6e-1);
        end
    end
    
    methods (Static = true)
        function seedRandStream(sim, seed)
            if sim.verbosity > 0
                fprintf('Seed %d\n', seed);
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
        
        function assertStateIsValid_Protein(sim)
            comp = sim.compartment;
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            rib = sim.state('Ribosome');
            pol = sim.state('Polypeptide');
            
            assertEqual(rib.boundMRNAs, pol.boundMRNAs);
            assertEqual(rib.mRNAPositions, pol.nascentMonomerLengths);
            assertEqual(rib.tmRNAPositions, pol.proteolysisTagLengths);
            
            if any(rib.states == 0)
                assertAllEqual(0, rib.boundMRNAs(rib.states == 0));
                assertAllEqual(0, rib.mRNAPositions(rib.states == 0));
                assertAllEqual(0, rib.tmRNAPositions(rib.states == 0));
            end
            assertTrue(all(rib.mRNAPositions(rib.states ~= 0) <= pol.monomerLengths(rib.boundMRNAs(rib.states ~= 0))));
            assertTrue(all(rib.mRNAPositions >= 0));
            if any(rib.states == rib.activeValue)
                assertAllEqual(0, rib.tmRNAPositions(rib.states == rib.activeValue));
            end
            assertTrue(all(rib.tmRNAPositions >= 0));
            assertTrue(all(rib.tmRNAPositions <= pol.proteolysisTagLength));
            
            assertEqual(pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), comp.cytosolIndexs), sum(rib.states ~= 0));
            assertEqual(rna.counts(rna.boundIndexs(rna.matureTMRNAIndexs), comp.cytosolIndexs), sum(rib.states == rib.stalledValue));
            
            assertTrue(all(rna.counts(:) >= 0));
            assertTrue(all(pm.counts(:) >= 0));
            assertTrue(all(pc.counts(:) >= 0));
        end
    end
end