%Simulation test cases
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef Simulation_Integrated_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = Simulation_Integrated_Test(name)
            this = this@TestCase(name);
        end
        
        function setUp(this)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
            sim.applyOptions('lengthSec', 10);
            this.simulation = sim;
        end
        
        function tearDown(this)
            this.simulation = [];
        end
    end
    
    %allocate memory test
    methods
        function testAllocateMemoryForState(this)
            this.simulation.allocateMemoryForState(1);  %no run-time errors
        end
    end
    
    %initialize state test
    methods
        function testInitializeState(this)
            sim = this.simulation;
            
            %assert process initialization order
            [~, processInitOrderIndexs] = ismember(cellfun(@(x) ['Process_' x], {
                'ProteinDecay'
                'MacromolecularComplexation'
                'RibosomeAssembly'
                'FtsZPolymerization'
                'ProteinFolding'
                'ProteinActivation'
                'Translation'
                'Metabolism'
                'TranscriptionalRegulation'
                'ReplicationInitiation'
                'DNASupercoiling'
                'ChromosomeCondensation'
                'DNARepair'
                'Transcription'
                'ChromosomeSegregation'
                'Cytokinesis'
                'DNADamage'
                'ProteinModification'
                'ProteinProcessingI'
                'ProteinProcessingII'
                'ProteinTranslocation'
                'Replication'
                'RNADecay'
                'RNAModification'
                'RNAProcessing'
                'TerminalOrganelleAssembly'
                'tRNAAminoacylation'
                'HostInteraction'
                }, 'UniformOutput', false), sim.processWholeCellModelIDs);
            assertEqual(processInitOrderIndexs, double(sim.processInitOrderIndexs));
            
            %% initialize state
            sim.allocateMemoryForState(1);
            sim.initializeState();
            
            %% assertions
            this.helpTime_testInitializeState();           %time
            this.helpMass_testInitializeState();           %cell mass
            this.helpGeometry_testInitializeState();       %cell shape
            this.helpChromosome_testInitializeState();     %chromosomes
            this.helpStimuli_testInitializeState();        %stimuli
            this.helpMetabolites_testInitializeState();    %metabolites
            this.helpRNAs_testInitializeState();           %RNA
            this.helpMonomers_testInitializeState();       %protein monomers
            this.helpComplexs_testInitializeState();       %protein complexs
            this.helpFtsZRing_testInitializeState();       %FtsZ ring
            
            this.helpGrowth_testInitializeState();         %growth
            this.helpRNAPolymerases_testInitializeState(); %RNA polymerase
            this.helpRibosomes_testInitializeState();      %ribosome
        end
        
        function helpTime_testInitializeState(this)
            sim = this.simulation;
            assertEqual(0, sim.state('Time').values);
        end
        
        function helpMass_testInitializeState(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            sim = this.simulation;
            g = sim.gene;
            s = sim.state('Mass');
            m = sim.state('Metabolite');
            rna = sim.state('Rna');
            pc = sim.state('ProteinComplex');
            transcript = sim.state('Transcript');
            polypeptide = sim.state('Polypeptide');
            initDry = s.cellInitialDryWeight;
            assertElementsAlmostEqual(initDry / (1-s.fractionWetWeight), sum(s.cell), 'relative', 10e-2, 0);
            
            cIdx = sim.compartment.cytosolIndexs;
            mIdx = sim.compartment.membraneIndexs;
            assertElementsAlmostEqual(initDry * s.dryWeightFractionLipid, ...
                sum(m.counts(m.lipidIndexs, [cIdx mIdx])' * m.molecularWeights(m.lipidIndexs)) / ConstantUtil.nAvogadro, 'relative', ...
                10e-2, 0);
            assertElementsAlmostEqual(initDry * s.dryWeightFractionPolyamine, ...
                sum(m.counts(m.polyamineIndexs, [cIdx mIdx])' * m.molecularWeights(m.polyamineIndexs)) / ConstantUtil.nAvogadro, 'relative', ...
                10e-2, 0);
            assertElementsAlmostEqual(initDry * s.dryWeightFractionCarbohydrate, ...
                sum(m.counts(m.carbohydrateIndexs, [cIdx mIdx])' * m.molecularWeights(m.carbohydrateIndexs)) / ConstantUtil.nAvogadro, 'relative', ...
                10e-2, 0);
            assertElementsAlmostEqual(initDry * s.dryWeightFractionIon, ...
                sum(m.counts(m.ionIndexs, [cIdx mIdx])' * m.molecularWeights(m.ionIndexs)) / ConstantUtil.nAvogadro, 'relative', ...
                1, 0);
            assertElementsAlmostEqual(initDry * s.dryWeightFractionVitamin, ...
                sum(m.counts(m.vitaminIndexs, [cIdx mIdx])' * m.molecularWeights(m.vitaminIndexs)) / ConstantUtil.nAvogadro, 'relative', ...
                10e-2, 0);
            
            assertElementsAlmostEqual(s.dryWeightFractionDNA, ...
                (sum(s.dnaWt) + m.counts(m.dntpIndexs, cIdx)' * m.molecularWeights(m.dntpIndexs) / ConstantUtil.nAvogadro)/ sum(s.cellDry), ...
                'relative', 10e-2, 0);
            assertElementsAlmostEqual(...
                s.dryWeightFractionRNA * s.initialFractionNTPsInRNAs + s.dryWeightFractionProtein * s.initialFractionAAsInMonomers, ...
                (sum(s.rnaWt + s.proteinWt) + transcript.dryWeight + polypeptide.dryWeight) / sum(s.cellDry), 'relative', 5e-2, 0);
            assertElementsAlmostEqual(...
                s.dryWeightFractionRNA * s.initialFractionNTPsInRNAs * sum(s.cellDry), ...
                (sum(s.rnaWt) + transcript.dryWeight + (sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * ...
                sum(pc.counts(pc.matureIndexs, :) + pc.counts(pc.boundIndexs, :), 2))' * ...
                rna.molecularWeights(rna.matureIndexs(setdiff(1:end, rna.matureMRNAIndexs))) / ConstantUtil.nAvogadro), ...
                'relative', 5e-2, 0);
            assertIn(transcript.dryWeight, [0 0.2 * s.dryWeightFractionRNA * sum(s.cellDry)]);
            assertIn(polypeptide.dryWeight, [0 0.01 * s.dryWeightFractionProtein * sum(s.cellDry)]);
        end
        
        function helpGeometry_testInitializeState(this)
            s = this.simulation;
            g = s.state('Geometry');
            
            assertElementsAlmostEqual(1.1910e-017, g.volume, 'relative', 10e-2, 0); %mass = 1.89e-14 * ln(2) g, density = 1100 g/L
            assertElementsAlmostEqual(2.8334e-007, g.width, 'relative', 10e-2, 0);  %mass = 1.89e-14 * ln(2) g, density = 1100 g/L
            assertEqual(g.pinchedDiameter, g.width);
            assertElementsAlmostEqual(0, g.cylindricalLength, 'absolute', 1e-8);
            assertElementsAlmostEqual(4 * pi * (g.width/2)^2, g.surfaceArea, 'relative', 10e-2, 0);
            assertElementsAlmostEqual(g.width, g.totalLength, 'relative', 10e-2, 0);
            assertFalse(g.pinched);
        end
        
        function helpChromosome_testInitializeState(this)
            sim = this.simulation;
            s = sim.state('Chromosome');
            chrLen = size(s.sequence, 1);
            
            %% 1 chromosome polymerized
            assertEqual(2, nnz(s.polymerizedRegions));
            assertAllEqual(size(s.sequence, 1), s.polymerizedRegions([1 1; 1 2]));
            
            %% supercoiled
            assertEqual([true; true], s.supercoiled([1 1; 1 2]));
            assertEqual(2, collapse(s.supercoiled));
            
            %% no proteins bound, except
            %- SMC
            %- transcription factors (spaced approximately 1/7000 bases)
            %- DnaA-ATP
            %- gyrase (topoisomerase IV doesn't bind at sigma = -0.06, and
            %  topoisomerase I isn't processive)
            m = sim.process('ChromosomeCondensation');
            d = sim.process('ReplicationInitiation');
            tr = sim.process('TranscriptionalRegulation');
            c = sim.process('DNASupercoiling');
            t = sim.process('Transcription');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            iMnmr = [
                tr.enzymeMonomerGlobalIndexs];
            iCplx = [
                m.enzymeComplexGlobalIndexs(m.enzymeComplexLocalIndexs == m.enzymeIndexs_SMC_ADP);
                d.enzymeComplexGlobalIndexs(d.enzymeComplexLocalIndexs == d.enzymeIndexs_DnaA_1mer_ATP);
                c.enzymeComplexGlobalIndexs(c.enzymeComplexLocalIndexs == c.enzymeIndexs_gyrase);
                tr.enzymeComplexGlobalIndexs;
                t.enzymeGlobalIndexs(t.enzymeIndexs_rnaPolymerase);
                t.enzymeGlobalIndexs(t.enzymeIndexs_rnaPolymeraseHoloenzyme)];
            [~, bndMnmrs] = find(s.monomerBoundSites);
            [~, bndCplxs] = find(s.complexBoundSites);
            assertTrue(all(ismember(unique(bndMnmrs), iMnmr)));
            assertTrue(all(ismember(unique(bndCplxs), iCplx)));
            
            %SMC
            assertIn(sum(m.enzymes + m.boundEnzymes), [0.75 * chrLen / m.smcSepNt Inf]);
            assertElementsAlmostEqual(chrLen / m.smcSepNt, nnz(s.complexBoundSites == iCplx(1)), 'relative', 0.25, 0);
            
            %DnaA
            %- all bound in DnaA-ATP 1 mer configuration
            %- converged (see replication initiation unit test for more
            %  stringent test)
            assertAllEqual(0, pm.counts(pm.matureIndexs(d.enzymeMonomerGlobalIndexs)));
            assertTrue(all(20 >= pc.counts(pc.matureIndexs(d.enzymeComplexGlobalIndexs))));
            assertAllEqual(0, pm.counts(pm.boundIndexs(d.enzymeMonomerGlobalIndexs)));
            assertAllEqual(0, pc.counts(pc.boundIndexs(setdiff(d.enzymeComplexGlobalIndexs, iCplx(2)))));
            assertIn(pc.counts(pc.boundIndexs(iCplx(2))), [1 Inf]);
            
            monomerBoundSites = s.monomerBoundSites;
            complexBoundSites = s.complexBoundSites;
            d.copyFromState();
            status = d.calculateDnaABoxStatus();
            ssBound8 = sum(status(d.dnaABoxIndexs_8mer, 1) == d.dnaABoxStatus_DnaAATPBound);
            for i = 1:10
                d.evolveState();
            end
            status = d.calculateDnaABoxStatus();
            bound8 = sum(status(d.dnaABoxIndexs_8mer, 1) == d.dnaABoxStatus_DnaAATPBound);
            assertElementsAlmostEqual(ssBound8, bound8, 'relative', 0.3);
            d.copyFromState();
            s.monomerBoundSites = monomerBoundSites;
            s.complexBoundSites = complexBoundSites;
            
            %gyrase
            freGyrase = c.enzymes(c.enzymeIndexs_gyrase);
            bndGyrase = c.boundEnzymes(c.enzymeIndexs_gyrase);
            totGyrase = freGyrase + bndGyrase;
            meanBndGyrase = totGyrase*(1-1/c.gyraseMeanDwellTime/c.stepSizeSec);
            assertIn(bndGyrase, [meanBndGyrase-1*sqrt(meanBndGyrase) totGyrase]);
            
            %% no damage, except methylation at R/M sites
            assertEqual(s.damagedBases, s.damagedSites_nonRedundant);
            
            m = sim.process('DNARepair');
            iMet = sim.state('Metabolite').getIndexs({'m6AD'});
            posStrnds = [
                reshape(m.RM_MunI_RecognitionSites(:, m.RM_MunI_MethylatedPositions), [], 1) ...
                reshape(repmat([1 2], size(m.RM_MunI_RecognitionSites, 1), 1), [], 1)];
            assertEqual(edu.stanford.covert.util.CircularSparseMat(posStrnds, iMet, [chrLen 4], 1), s.damagedBases);
            
            %% chromsomes not segregated
            assertFalse(s.segregated);
            
            %% bound proteins accounted
            [~, vals] = find(s.monomerBoundSites);
            tfs = true(size(pm.boundIndexs));
            tfs(pm.translationFactorIndexs) = false;
            assertEqual(pm.counts(pm.boundIndexs, sim.compartment.cytosolIndexs) .* tfs, reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1));
            
            [~, vals] = find(s.complexBoundSites);
            tfs = true(size(pc.boundIndexs));
            tfs(pc.translationFactorIndexs) = false;
            tfs(pc.ribosome70SIndexs) = false;
            tfs(pc.ftsZGTPIndexs) = false;
            tfs(pc.ftsZGDPIndexs) = false;
            assertEqual(pc.counts(pc.boundIndexs, sim.compartment.cytosolIndexs) .* tfs, reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1));
        end
        
        function helpStimuli_testInitializeState(this)
            import edu.stanford.covert.cell.sim.constant.Condition;
            
            sim = this.simulation;
            stim = sim.state('Stimulus');
            
            %conditions satisfied
            assertEqual(stim.values, Condition.applyConditions(stim.values, stim.setValues, 0));
        end
        
        function helpMetabolites_testInitializeState(this)
            import edu.stanford.covert.cell.sim.constant.Condition;
            import edu.stanford.covert.util.ConstantUtil;
            
            sim = this.simulation;
            met = sim.state('Metabolite');
            mets = met.counts;
            
            %conditions satisfied
            assertEqual(mets, Condition.applyConditions(mets, met.setCounts, 0));
            
            %all non-negative integers
            validateattributes(mets, {'numeric'}, {'nonnegative', 'integer'});
            
            %in correct compartments
            cIdx = sim.compartment.cytosolIndexs;
            mIdx = sim.compartment.membraneIndexs;
            eIdx = sim.compartment.extracellularIndexs;
            m = sim.process('Metabolism');
            hMetIdxs = met.hydrophobicIndexs;
            eMetIdxs = [
                met.setCounts(met.setCounts(:, Condition.compartmentIndexs) == eIdx, Condition.objectIndexs);
                m.substrateMetaboliteGlobalIndexs(ismember(m.substrateMetaboliteLocalIndexs, find(any(m.reactionStoichiometryMatrix(:, m.reactionIndexs_transport, m.compartmentIndexs_extracellular), 2))));
                ];
            assertAllEqual(0, mets(:, setdiff(1:end, [cIdx mIdx eIdx])));
            assertAllEqual(0, mets(hMetIdxs, cIdx));
            assertAllEqual(0, mets(setdiff(1:end, hMetIdxs), mIdx));
            assertAllEqual(0, mets(setdiff(1:end, eMetIdxs), eIdx));
            
            %little free dNTP, NTP, aa in cytosol
            metMWs = met.molecularWeights / ConstantUtil.nAvogadro;
            mass = sim.state('Mass');
            t = sim.state('Time');
            c = sim.state('Chromosome');
            assertElementsAlmostEqual(...
                2 * exp(-log(2) * (t.cellCycleLength - t.cytokinesisDuration) / t.cellCycleLength) - 1, ...
                1 / (2 * c.sequenceLen) * sum(mets(met.dntpIndexs, cIdx)), ...
                'relative', 60e-2, 0);
            assertElementsAlmostEqual(...
                2 * exp(-log(2) * (t.cellCycleLength - t.cytokinesisDuration) / t.cellCycleLength) - 1, ...
                1 / sum(mass.dnaWt) * (metMWs(met.dntpIndexs) - metMWs(met.diphosphateIndexs))' * mets(met.dntpIndexs, cIdx), ...
                'relative', 60e-2, 0);
            assertIn(1/sum(mass.proteinWt) * metMWs(met.aminoAcidIndexs)' * mets(met.aminoAcidIndexs, cIdx), [0 0.12]);
        end
        
        function helpRNAs_testInitializeState(this)
            import edu.stanford.covert.util.ComputationUtil;
            
            sim = this.simulation;
            rna = sim.state('Rna');
            pc = sim.state('ProteinComplex');
            
            %RNAs positive
            assertTrue(all(rna.counts(:) >= 0));
            
            %all RNAs mature (mature: mRNA, rRNA, sRNA; aminoacylated: tRNA, tmRNA)
            assertEqual(nnz(rna.counts), nnz(rna.counts([rna.matureIndexs; rna.aminoacylatedIndexs], :)));
            
            %mature RNAs distributed according expected expression
            rnaExp = ...
                + rna.expression(rna.matureIndexs) ...
                + rna.expression(rna.aminoacylatedIndexs);
            
            complexs = ...
                + pc.counts(pc.matureIndexs, :) ...
                + pc.counts(pc.boundIndexs, :) ...
                + pc.counts(pc.inactivatedIndexs, :);
            rnaSubunits = sum(pc.proteinComplexComposition, 3) * complexs;
            rnaSubunits(sim.gene.mRNAIndexs, :) = 0;
            RNAs = ...
                + rna.counts(rna.matureIndexs, :) ...
                + rna.counts(rna.aminoacylatedIndexs, :) + ...
                + ComputationUtil.invertCompositionMatrix(rna.matureRNAGeneComposition) * rnaSubunits;
            assertElementsAlmostEqual(rnaExp * sum(RNAs(:)), sum(RNAs, 2), 'relative', 0.35, 4);
            
            %in correct compartments
            cIdx = sim.compartment.cytosolIndexs;
            assertEqual(nnz(rna.counts(:, cIdx)), nnz(rna.counts));
            
            %no aborted transcripts
            assertTrue(isempty(sim.state('Transcript').abortedSequences));
        end
        
        function helpMonomers_testInitializeState(this)
            sim = this.simulation;
            t = sim.state('Time');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            %monomers positive
            assertTrue(all(pm.counts(:) >= 0));
            
            %all monomers mature, bound, or inactivated
            assertEqual(nnz(pm.counts), nnz(pm.counts([pm.matureIndexs; pm.boundIndexs; pm.inactivatedIndexs], :)));
            
            %monomers correctly active/inactive
            monomers = pm.counts;
            complexs = pc.counts;
            m = sim.process('ProteinActivation');
            m.copyFromState();
            m.evolveState();
            m.copyToState();
            assertEqual(monomers, pm.counts);
            assertEqual(complexs, pc.counts);
            
            %mature monomers distributed according expected expression
            monExp = (rna.matureRNAGeneComposition(sim.gene.mRNAIndexs, :) * rna.expression(rna.matureIndexs)) ./ ...
                (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
            monExp = monExp / sum(monExp);
            
            complexs = ...
                + pc.counts(pc.matureIndexs, :) ...
                + pc.counts(pc.boundIndexs, :) ...
                + pc.counts(pc.inactivatedIndexs, :);
            monomers = ...
                + pm.counts(pm.matureIndexs, :) ...
                + pm.counts(pm.boundIndexs, :) ...
                + pm.counts(pm.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(sim.gene.mRNAIndexs, :, :), 3) * complexs;
            assertIn(corr(monExp, sum(monomers, 2)), [0.95 1]);
            assertElementsAlmostEqual(monExp * sum(monomers(:)), sum(monomers, 2), 'relative', 0.70, 10);
            
            %in correct compartments
            monomers = ...
                + pm.counts(pm.matureIndexs, :) ...
                + pm.counts(pm.boundIndexs, :) ...
                + pm.counts(pm.inactivatedIndexs, :);
            assertEqual(nnz(monomers), nnz(monomers(sub2ind(size(monomers), (1:size(monomers, 1))', pm.compartments(pm.matureIndexs)))));
            
            %no aborted polypeptides
            m = sim.state('Polypeptide');
            assertTrue(isempty(m.abortedSequences));
        end
        
        function helpComplexs_testInitializeState(this)
            sim = this.simulation;
            
            %few more complexs can form
            ms = {
                sim.process('MacromolecularComplexation');
                sim.process('RibosomeAssembly');
                sim.process('ProteinFolding');
                };
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            RNAs = rna.counts;
            monomers = pm.counts;
            complexs = pc.counts;
            for i = 1:10
                for j = 1:numel(ms)
                    ms{j}.copyFromState();
                    ms{j}.evolveState();
                    ms{j}.copyToState();
                end
            end
            assertElementsAlmostEqual(RNAs, rna.counts, 'relative', 0.05, 2);
            assertElementsAlmostEqual(monomers, pm.counts, 'relative', 0.05, 2);
            assertElementsAlmostEqual(complexs, pc.counts, 'relative', 0.05, 2);
            
            %complexs positive
            assertTrue(all(pc.counts(:) >= 0));
            
            %all complexs mature, bound, or inactivated
            assertEqual(nnz(pc.counts), nnz(pc.counts([pc.matureIndexs; pc.boundIndexs; pc.inactivatedIndexs], :)));
            
            %complexs correctly active/inactive
            complexs = pc.counts;
            m = sim.process('ProteinActivation');
            m.copyFromState();
            m.evolveState();
            m.copyToState();
            assertEqual(complexs, pc.counts);
            
            %in correct compartments
            complexs = ...
                + pc.counts(pc.matureIndexs, :) ...
                + pc.counts(pc.boundIndexs, :) ...
                + pc.counts(pc.inactivatedIndexs, :);
            assertEqual(nnz(complexs), nnz(complexs(sub2ind(size(complexs), (1:size(complexs, 1))', pc.compartments(pc.matureIndexs)))));
        end
        
        function helpFtsZRing_testInitializeState(this)
            sim = this.simulation;
            
            %free polymers
            m = sim.process('FtsZPolymerization');
            initialEnzymes = m.enzymes;
            m.copyFromState();
            for i = 1:5
                m.evolveState();
            end
            m.copyFromState();
            assertElementsAlmostEqual(initialEnzymes, m.enzymes, 'absolute', 3);
            
            %ring
            s = sim.state('FtsZRing');
            assertEqual(0, s.numEdgesOneStraight);
            assertEqual(0, s.numEdgesTwoStraight);
            assertEqual(0, s.numEdgesTwoBent);
            assertEqual(0, s.numResidualBent);
        end
        
        function helpGrowth_testInitializeState(this)
            sim = this.simulation;
            s = sim.state('MetabolicReaction');
            assertElementsAlmostEqual(edu.stanford.covert.cell.sim.process.Metabolism_Test.expected_growthRate, s.growth, 'relative', 0.50, 0);
        end
        
        function helpRNAPolymerases_testInitializeState(this)
            sim = this.simulation;
            m = sim.process('Transcription');
            c = sim.state('Chromosome');
            
            r = m.rnaPolymerases;
            t = m.transcripts;
            
            %RNA polymerase distributed according to expectations
            %- actively transcribing
            %- specifically bound
            %- non-specifically bound
            %- free
            %nPols = m.boundEnzymes(m.enzymeIndexs_rnaPolymerase);
            %assertAllEqual(true, ...
            %    nPols * abs(r.stateExpectations - r.stateOccupancies) <= ...
            %    3 * sqrt(nPols * r.stateExpectations));
            
            %all RNA polymerase bound to accessible transcription units
            atPls = find(r.states >= r.activelyTranscribingValue);
            sbPls = find(r.states == r.specificallyBoundValue);
            atTUs = t.boundTranscriptionUnits(atPls);
            sbTUs = t.boundTranscriptionUnits(sbPls);
            assertAllEqual(true, c.polymerizedTranscriptionUnits([atTUs; sbTUs], 1));
            
            %RNA polymerase distributed according to binding probabilities and
            %lengths
            %probs = ...
            %    t.transcriptionUnitLengths .* ...
            %    m.transcriptionUnitBindingProbabilities .* ...
            %    r.transcriptionFactorBindingProbFoldChange(:,1) .* ...
            %    r.supercoilingBindingProbFoldChange(:,1) .* ...
            %    c.accessibleTranscriptionUnits(:, 1);
            %probs = probs / sum(probs);
            
            %exp = histc(atTUs, 1:numel(t.transcriptionUnitLengths)) / numel(atTUs);
            %assertTrue(norm(numel(atPls) * (probs  - exp(:)), 1) < 12 * sqrt(numel(atPls)));
            
            %exp = histc(sbTUs, 1:numel(t.transcriptionUnitLengths)) / numel(sbTUs);
            %assertTrue(norm(numel(sbPls) * (probs  - exp(:)), 1) < 10 * sqrt(numel(sbPls)));
            
            %RNA polymerase not beyond lengths of TUs
            assertAllEqual(true, r.states(atPls) <= t.transcriptionUnitLengths(atTUs));
            
            %RNA polymerase don't overlap (note: this is more thoroughly
            %ensured by the chromosome class)
            assertEqual(numel(atPls), size(unique([atTUs r.states(atPls)], 'rows'), 1));
            assertEqual(numel(sbTUs), numel(unique(sbTUs)));
        end
        
        function helpRibosomes_testInitializeState(this)
            sim = this.simulation;
            r = sim.state('Ribosome');
            p = sim.state('Polypeptide');
            rna = sim.state('Rna');
            
            %all ribosomes bound, no proteolysis tags being synthesized
            assertEqual(0, r.stateOccupancies(r.stalledIndex));
            assertEqual(0, numel(p.abortedSequences));
            
            %ribosomes bound to expressed mRNAs
            cIdx = sim.compartment.cytosolIndexs;
            bndRibs = r.states == r.activeValue;
            bndRNAs = r.boundMRNAs(bndRibs);
            expRNAs = multiprod(...
                rna.matureRNAGeneComposition(sim.gene.mRNAIndexs, rna.matureMRNAIndexs), ...
                rna.counts(rna.matureIndexs(rna.matureMRNAIndexs), cIdx, :), [1 2], 1);
            assertAllEqual(true, ismember(bndRNAs, find(expRNAs)));
            
            %ribosomes distributed according to mRNA expression and lengths
            probs = p.monomerLengths .* expRNAs;
            probs = probs / sum(probs);
            exp = histc(bndRNAs, 1:numel(p.monomerLengths)) / numel(bndRNAs);
            assertIn(norm(numel(bndRibs) * (probs  - exp), 1), [0 15 * sqrt(numel(bndRibs))]);
            
            %ribosomes not beyond length of monomers
            assertAllEqual(true, p.nascentMonomerLengths(bndRibs) < p.monomerLengths(bndRNAs));
        end
    end
    
    %other tests
    methods
        function testEvolveState(this)
            this.simulation.evolveState();  %no run-time errors
        end
        
        function testCommunicationStateToFromProcesses(this)
            %simulation
            sim = this.simulation;
            met = sim.state('Metabolite');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            stim = sim.state('Stimulus');
            
            %set values of components to unique numbers
            stimuli     = reshape((1:numel(stim.values)),                   size(stim.values));
            metabolites = reshape((1:numel(met.counts)) + stimuli(end),     size(met.counts));
            RNAs        = reshape((1:numel(rna.counts)) + metabolites(end), size(rna.counts));
            monomers    = reshape((1:numel(pm.counts))  + RNAs(end),        size(pm.counts));
            complexs    = reshape((1:numel(pc.counts))  + monomers(end),    size(pc.counts));
            
            %set values of simulation components
            stim.values  = stimuli;
            met.counts   = metabolites;
            rna.counts   = RNAs;
            pm.counts    = monomers;
            pc.counts    = complexs;
            
            %communicate values of components to/from processes
            for i = 1:length(sim.processes)
                process = sim.processes{i};
                process.copyFromState();
                process.copyToState();
            end
            
            %assert that values of components haven't changed as a result
            %of communication to and then from processes
            assertEqual(stimuli,     stim.values,  'State <==> Process communication error');
            assertEqual(metabolites, met.counts,   'State <==> Process communication error');
            assertEqual(RNAs,        rna.counts,   'State <==> Process communication error');
            assertEqual(monomers,    pm.counts,    'State <==> Process communication error');
            assertEqual(complexs,    pc.counts,    'State <==> Process communication error');
        end
        
        %test communication simulation ==> processes
        function testStateToProcessCommunication(this)
            sim = this.simulation;
            met = sim.state('Metabolite');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            stim = sim.state('Stimulus');
            
            %set values of components to unique numbers
            stim.values(:)  = (1:numel(stim.values));
            met.counts(:)   = (1:numel(met.counts)) + stim.values(end);
            rna.counts(:)   = (1:numel(rna.counts)) + met.counts(end);
            pm.counts(:)    = (1:numel(pm.counts))  + rna.counts(end);
            pc.counts(:)    = (1:numel(pc.counts))  + pm.counts(end);
            
            for i = 1:length(sim.processes)
                process = sim.processes{i};
                
                %clear process state
                process.allocateMemoryForState(1);
                
                %communicate simulation => process
                process.copyFromState();
                
                %stimuli
                for j = 1:length(process.stimuliWholeCellModelIDs)
                    properCommunication = false;
                    
                    %stimuli
                    idxs = find(strcmp(process.stimuliWholeCellModelIDs{j}, stim.wholeCellModelIDs));
                    for k = 1:numel(idxs)
                        if all(ismember(process.stimuli(j,:), stim.values(idxs(k),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %Metabolites
                    idxs = find(strcmp(process.stimuliWholeCellModelIDs{j}, met.wholeCellModelIDs));
                    for k = 1:numel(idxs)
                        if all(ismember(process.stimuli(j,:), met.counts(idxs(k),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %RNAs
                    idxs = find(strcmp(process.stimuliWholeCellModelIDs{j}, rna.wholeCellModelIDs(rna.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.stimuli(j,:), rna.counts(rna.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %monomers
                    idxs = find(strcmp(process.stimuliWholeCellModelIDs{j}, pm.wholeCellModelIDs(pm.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.stimuli(j,:), pm.counts(pm.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %complexs
                    idxs = find(strcmp(process.stimuliWholeCellModelIDs{j}, pc.wholeCellModelIDs(pc.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.stimuli(j,:), pc.counts(pc.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    assertTrue(false, sprintf('%s %s stimuli values not correctly communicated from simulation.', ...
                        process.wholeCellModelID, process.stimuliWholeCellModelIDs{j}));
                end
                
                %substrates
                for j=1:length(process.substrateWholeCellModelIDs)
                    properCommunication = false;
                    
                    %stimuli
                    idxs = find(strcmp(process.substrateWholeCellModelIDs{j}, stim.wholeCellModelIDs));
                    for k = 1:numel(idxs)
                        if all(ismember(process.substrates(j,:), stim.values(idxs(k),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %Metabolites
                    idxs = find(strcmp(process.substrateWholeCellModelIDs{j}, met.wholeCellModelIDs));
                    for k = 1:numel(idxs)
                        if all(ismember(process.substrates(j,:), met.counts(idxs(k),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %RNAs
                    idxs = find(strcmp(process.substrateWholeCellModelIDs{j}, rna.wholeCellModelIDs(rna.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.substrates(j,:), rna.counts(rna.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %monomers
                    idxs = find(strcmp(process.substrateWholeCellModelIDs{j}, pm.wholeCellModelIDs(pm.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.substrates(j,:), pm.counts(pm.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %complexs
                    idxs = find(strcmp(process.substrateWholeCellModelIDs{j}, pc.wholeCellModelIDs(pc.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.substrates(j,:), pc.counts(pc.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    assertTrue(false, sprintf('%s %s substrate values not correctly communicated from simulation.',...
                        process.wholeCellModelID, process.substrateWholeCellModelIDs{j}));
                end
                
                %enzymes
                for j=1:length(process.enzymeWholeCellModelIDs)
                    properCommunication = false;
                    
                    %stimuli
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, stim.wholeCellModelIDs));
                    for k = 1:numel(idxs)
                        if all(ismember(process.enzymes(j,:), stim.values(idxs(k),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %Metabolites
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, met.wholeCellModelIDs));
                    for k = 1:numel(idxs)
                        if all(ismember(process.enzymes(j,:), met.counts(idxs(k),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %RNAs
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, rna.wholeCellModelIDs(rna.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.enzymes(j,:), rna.counts(rna.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %monomers
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, pm.wholeCellModelIDs(pm.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.enzymes(j,:), pm.counts(pm.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %complexs
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, pc.wholeCellModelIDs(pc.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.enzymes(j,:), pc.counts(pc.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    assertTrue(false, sprintf('%s %s enzyme values not correctly communicated from simulation.',...
                        process.wholeCellModelID, process.enzymeWholeCellModelIDs{j}));
                end
                
                %bound enzymes
                for j=1:length(process.enzymeWholeCellModelIDs)
                    properCommunication = false;
                    
                    %RNAs
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, rna.wholeCellModelIDs(rna.boundIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.boundEnzymes(j,:), rna.counts(rna.boundIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %monomers
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, pm.wholeCellModelIDs(pm.boundIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.boundEnzymes(j,:), pm.counts(pm.boundIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %complexs
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, pc.wholeCellModelIDs(pc.boundIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.boundEnzymes(j,:), pc.counts(pc.boundIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    assertTrue(false, sprintf('%s %s bound enzyme values not correctly communicated from simulation.', ...
                        process.wholeCellModelID, process.enzymeWholeCellModelIDs{j}));
                end
            end
        end
        
        function testProcessToStateCommunication(this)
            sim = this.simulation;
            met = sim.state('Metabolite');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            stim = sim.state('Stimulus');
            
            %loop over processes
            for i = 1:length(sim.processes)
                sim.allocateMemoryForState(1);
                
                %set state to unique values
                stim.values(:)  = (1:numel(stim.values));
                met.counts(:)   = (1:numel(met.counts))   + stim.values(end);
                rna.counts(:)   = (1:numel(rna.counts))   + met.counts(end);
                pm.counts(:)    = (1:numel(pm.counts))    + rna.counts(end);
                pc.counts(:)    = (1:numel(pc.counts))    + pm.counts(end);
                
                %process
                process = sim.processes{i};
                
                %communicate state => process
                process.copyFromState();
                
                %set process state to unique numbers
                assertTrue(all(process.stimuli(:)));
                assertTrue(all(process.substrates(:)));
                assertTrue(all(process.enzymes(:)));
                assertTrue(all(process.boundEnzymes(:)));
                
                %communicate process => state
                process.copyToState();
                
                %substrates
                for j=1:length(process.substrateWholeCellModelIDs)
                    properCommunication = false;
                    
                    %stimuli
                    if any(strcmp(process.substrateWholeCellModelIDs{j}, stim.wholeCellModelIDs))
                        break;
                    end
                    
                    %Metabolites
                    idxs = find(strcmp(process.substrateWholeCellModelIDs{j}, sim.state('Metabolite').wholeCellModelIDs));
                    for k = 1:numel(idxs)
                        if all(ismember(process.substrates(j,:), met.counts(idxs(k),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %RNAs
                    idxs = find(strcmp(process.substrateWholeCellModelIDs{j}, rna.wholeCellModelIDs(rna.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.substrates(j,:), rna.counts(rna.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %monomers
                    idxs = find(strcmp(process.substrateWholeCellModelIDs{j}, pm.wholeCellModelIDs(pm.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.substrates(j,:), pm.counts(pm.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %complexs
                    idxs = find(strcmp(process.substrateWholeCellModelIDs{j}, pc.wholeCellModelIDs(pc.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.substrates(j,:), pc.counts(pc.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    assertTrue(false, sprintf('%s %s substrate values not correctly communicated to simulation.', ...
                        process.wholeCellModelID, process.substrateWholeCellModelIDs{j}));
                end
                
                %enzymes
                for j=1:length(process.enzymeWholeCellModelIDs)
                    properCommunication = false;
                    
                    %RNAs
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, rna.wholeCellModelIDs(rna.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.enzymes(j,:), rna.counts(rna.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %monomers
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, pm.wholeCellModelIDs(pm.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.enzymes(j,:), pm.counts(pm.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %complexs
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, pc.wholeCellModelIDs(pc.matureIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.enzymes(j,:), pc.counts(pc.matureIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    assertTrue(false, sprintf('%s %s enzyme values not correctly communicated to simulation.', ...
                        process.wholeCellModelID, process.enzymeWholeCellModelIDs{j}));
                end
                
                %bound enzymes
                for j=1:length(process.enzymeWholeCellModelIDs)
                    properCommunication = false;
                    
                    %RNAs
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, rna.wholeCellModelIDs(rna.boundIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.boundEnzymes(j,:), rna.counts(rna.boundIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %monomers
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, pm.wholeCellModelIDs(pm.boundIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.boundEnzymes(j,:), pm.counts(pm.boundIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    %complexs
                    idxs = find(strcmp(process.enzymeWholeCellModelIDs{j}, pc.wholeCellModelIDs(pc.boundIndexs)));
                    for k = 1:numel(idxs)
                        if all(ismember(process.boundEnzymes(j,:), pc.counts(pc.boundIndexs(idxs(k)),:)))
                            properCommunication = true;
                            break;
                        end
                    end
                    if properCommunication; break; end;
                    
                    assertTrue(false, sprintf('%s %s bound enzyme values not correctly communicated to simulation.', ...
                        process.wholeCellModelID, process.enzymeWholeCellModelIDs{j}));
                end
            end
        end
    end
    
    %test individual processes
    methods
        function testMacromolecularComplexation(this)
            %simulation
            sim = this.simulation;
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            %subunits available for complexation
            rna.counts(rna.matureIndexs, :) = 10;
            pm.counts(pm.matureIndexs, :) = 10;
            
            %hold onto initial state
            initial_complexs = pc.counts;
            
            %run sim
            sim.evolveState();
            
            %check complexes being produced
            assertTrue(sum(initial_complexs(:)) < sum(pc.counts(:)), ...
                'Macromolecular complexation process didn''t produce any complexes');
        end
        
        function testMetabolism(this)
            sim = this.simulation;
            m = sim.process('Metabolism');
            s = sim.state('MetabolicReaction');
            
            initial_substrates = m.substrates;
            
            sim.evolveState();
            
            %check growth rate
            assertElementsAlmostEqual(edu.stanford.covert.cell.sim.process.Metabolism_Test.expected_growthRate, ...
                s.growth, 'relative', 35e-2, 0);
            
            %check substrates being produced
            assertFalse(all(initial_substrates(:) == m.substrates(:)), ...
                'Metabolism process didn''t produce any substrates')
        end
        
        function testTranscription(this)
            sim = this.simulation;
            t = sim.process('Transcription');
            m = sim.process('Metabolism');
            met = sim.state('Metabolite');
            
            r = t.rnaPolymerases;
            transcripts = t.transcripts;
            
            %lots of metabolites
            met.counts = 10 * met.counts;
            
            %more metabolic enzymes
            m.enzymes = 10 * m.enzymes;
            m.copyToState();
            
            %more RNA polymerase
            numRNAPolymerase = 10;
            t.enzymes = t.enzymes + t.boundEnzymes;
            t.enzymes(t.enzymeIndexs_rnaPolymerase, :) = numRNAPolymerase;
            t.enzymes(t.enzymeIndexs_rnaPolymeraseHoloenzyme, :) = 0;
            t.boundEnzymes(:) = 0;
            r.states = repmat(r.freeValue, numRNAPolymerase, 1);
            r.positionStrands(:) = r.notExistValue;
            transcripts.boundTranscriptionUnits = zeros(numRNAPolymerase, 1);
            transcripts.boundTranscriptProgress = zeros(numRNAPolymerase, 1);
            transcripts.boundTranscriptChromosome = zeros(numRNAPolymerase, 1);
            transcripts.abortedTranscripts = zeros(0, 2);
            t.copyToState();
            
            %hold on to initial numbers of RNAs
            initial_RNAs = t.RNAs;
            
            %evolve state
            for i = 1:100
                sim.evolveState();
                
                if any(t.RNAs > initial_RNAs)
                    break;
                end
            end
            
            %check RNAs being produced
            assertTrue(any(t.RNAs > initial_RNAs), 'RNAs not produced');
        end
    end
    
    %test groups of processes
    methods
        function testMetabolismTranscriptionCommunication(this)
            sim = this.simulation;
            metabolism = sim.process('Metabolism');
            transcription = sim.process('Transcription');
            
            [~, idxs_transcription, idxs_metabolism] = intersect(...
                transcription.substrateWholeCellModelIDs, ...
                metabolism.substrateWholeCellModelIDs);
            
            initial_substrates = metabolism.substrates;
            
            sim.evolveState();
            metabolism.copyFromState();
            transcription.copyFromState();
            
            assertFalse(all(initial_substrates(:) == metabolism.substrates(:)),...
                'Substrates didn''t change during simulation');
            assertEqual(...
                metabolism.substrates(idxs_metabolism, metabolism.compartmentIndexs_cytosol), ...
                transcription.substrates(idxs_transcription, :),...
                'Error in communication of metabolites between metabolism and transcription');
        end
        
        function testMassBalance(this)
            import edu.stanford.covert.cell.sim.constant.Condition;
            
            sim = this.simulation;
            
            sim.state('Metabolite').setCounts = zeros(0, size(sim.state('Metabolite').setCounts, 2));
            
            mass = sum(sim.state('Mass').total);
            for i = 1:100
                sim.evolveState();
            end
            assertElementsAlmostEqual(mass, sum(sim.state('Mass').total), 'relative', 1e-6, 0);
        end
        
        function testMediaSupportsGrowth(~)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'Metabolism'});
            
            time = sim.state('Time');
            m = sim.state('Metabolite');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            mr = sim.state('MetabolicReaction');
            met = sim.process('Metabolism');
            
            initialGrowth = mr.growth;
            initialMetabolites = m.counts;
            initialMonomers = pm.counts;
            initialComplexs = pc.counts;
            
            sim.evolveState();
            tmpIdxs = met.substrateMetaboliteGlobalCompartmentIndexs(:, met.compartmentIndexs_extracellular);
            neededMetabolites = 1/initialGrowth * (initialMetabolites(tmpIdxs) - m.counts(tmpIdxs));
            
            assertTrue(all(initialMetabolites(tmpIdxs) >= neededMetabolites), ...
                sprintf(sprintf('Insufficient metabolites:\n- %s', strjoin('\n- ', ...
                m.wholeCellModelIDs{met.substrateMetaboliteGlobalIndexs(initialMetabolites(tmpIdxs) < neededMetabolites, 1)}))));
            
            for i = 1:time.cellCycleLength
                pm.counts = ceil(initialMonomers * exp(i * log(2) / time.cellCycleLength));
                pc.counts = ceil(initialComplexs * exp(i * log(2) / time.cellCycleLength));
                sim.evolveState();
            end
            
            assertElementsAlmostEqual(2 * initialGrowth, mr.growth, 'relative', 1e-2, 0);
        end
        
        function testRNADecayIntegration(this)
            sim = this.simulation;
            c = sim.compartment;
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            rib = sim.state('Ribosome');
            pol = sim.state('Polypeptide');
            prtDcy = sim.process('ProteinDecay');
            
            %% decay actively translating mRNA
            %setup by with active ribosome, decayRates = 0, decay enzymes = 0
            pm.counts(prtDcy.enzymeMonomerGlobalCompartmentIndexs) = 0;
            pc.counts(prtDcy.enzymeComplexGlobalCompartmentIndexs) = 0;
            
            r.counts(r.matureIndexs, c.cytosolIndexs) = 0;
            r.counts(r.matureIndexs(r.matureMRNAIndexs(1:2)), c.cytosolIndexs) = 1;
            r.counts(r.boundIndexs(r.matureTMRNAIndexs), c.cytosolIndexs) = 0;
            r.counts(r.matureIndexs(r.matureTMRNAIndexs), c.cytosolIndexs) = 1;
            
            pc.counts(pc.matureIndexs(pc.ribosome30SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.ribosome50SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.ribosome70SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs) = 2;
            pm.counts(pm.matureIndexs(pm.translationFactorIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.translationFactorIndexs), c.cytosolIndexs) = 0;
            pm.counts(pm.boundIndexs(pm.translationFactorIndexs), c.cytosolIndexs) = 2;
            pc.counts(pc.boundIndexs(pc.translationFactorIndexs), c.cytosolIndexs) = 2;
            
            rib.states = [rib.activeValue; rib.activeValue];
            rib.boundMRNAs = [1; 2];
            rib.mRNAPositions = [10; 0];
            rib.tmRNAPositions = [0; 0];
            pol.boundMRNAs = rib.boundMRNAs;
            pol.nascentMonomerLengths = rib.mRNAPositions;
            pol.proteolysisTagLengths = rib.tmRNAPositions;
            pol.abortedPolypeptides = zeros(0, 3);
            
            m.counts(:) = 0;
            m.counts(m.waterIndexs, c.cytosolIndexs) = 1e7;
            
            r.decayRates(:) = 0;
            r.decayRates(r.matureIndexs(r.matureMRNAIndexs(1))) = 10;
            r.decayRates(r.matureIndexs(r.matureMRNAIndexs(2))) = 10;
            
            %evolve and assert that mRNA degraded, tmRNA + ribosome released
            sim.evolveState();
            
            assertEqual(0, r.counts(r.matureIndexs(r.matureMRNAIndexs(1)), c.cytosolIndexs));
            assertEqual(0, r.counts(r.boundIndexs(r.matureTMRNAIndexs), c.cytosolIndexs));
            assertEqual(1, r.counts(r.matureIndexs(r.matureTMRNAIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.matureIndexs(pc.ribosome70SIndexs), c.cytosolIndexs));
            assertEqual(2, pc.counts(pc.matureIndexs(pc.ribosome30SIndexs), c.cytosolIndexs)+pc.counts(pc.matureIndexs(pc.getIndexs('RIBOSOME_30S_IF3')), c.cytosolIndexs));
            assertEqual(2, pc.counts(pc.matureIndexs(pc.ribosome50SIndexs), c.cytosolIndexs));
            assertEqual(2, pm.counts(pm.matureIndexs(pm.translationFactorIndexs), c.cytosolIndexs) + pm.counts(pm.boundIndexs(pm.translationFactorIndexs), c.cytosolIndexs));
            assertAllEqual(2, pc.counts(pc.matureIndexs(pc.translationFactorIndexs), c.cytosolIndexs) + pc.counts(pc.boundIndexs(pc.translationFactorIndexs), c.cytosolIndexs));
            assertEqual(1, numel(pol.abortedSequences));
            assertEqual(10, numel(pol.abortedSequences{1}));
            
            %add protein decay enzymes, evolve, and assert that polypeptide is degraded
            pm.counts(prtDcy.enzymeMonomerGlobalCompartmentIndexs) = 10;
            pc.counts(prtDcy.enzymeComplexGlobalCompartmentIndexs) = 10;
            
            sim.evolveState();
            assertEqual(0, numel(pol.abortedSequences));
            
            %% decay stalled mRNA
            %setup by with 1 stalled ribosome, decayRates = 0, decay enzymes = 0
            pm.counts(prtDcy.enzymeMonomerGlobalCompartmentIndexs) = 0;
            pc.counts(prtDcy.enzymeComplexGlobalCompartmentIndexs) = 0;
            
            r.counts(r.matureIndexs, c.cytosolIndexs) = 0;
            r.counts(r.matureIndexs(r.matureMRNAIndexs(1)), c.cytosolIndexs) = 1;
            r.counts(r.boundIndexs(r.matureTMRNAIndexs), c.cytosolIndexs) = 1;
            
            pc.counts(pc.matureIndexs(pc.getIndexs('RIBOSOME_30S_IF3')), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.ribosome30SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.ribosome50SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.ribosome70SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs) = 1;
            pm.counts(pm.matureIndexs(pm.translationFactorIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.translationFactorIndexs), c.cytosolIndexs) = 0;
            pm.counts(pm.boundIndexs(pm.translationFactorIndexs), c.cytosolIndexs) = 1;
            pc.counts(pc.boundIndexs(pc.translationFactorIndexs), c.cytosolIndexs) = 1;
            
            rib.states = rib.stalledValue;
            rib.boundMRNAs = 1;
            rib.mRNAPositions = 10;
            rib.tmRNAPositions = 5;
            pol.boundMRNAs = rib.boundMRNAs;
            pol.nascentMonomerLengths = rib.mRNAPositions;
            pol.proteolysisTagLengths = rib.tmRNAPositions;
            pol.abortedPolypeptides = zeros(0, 3);
            
            m.counts(:) = 0;
            m.counts(m.waterIndexs, c.cytosolIndexs) = 1e7;
            
            r.decayRates(:) = 0;
            r.decayRates(r.matureIndexs(r.matureMRNAIndexs(1))) = 10;
            
            %evolve and assert that mRNA degraded, tmRNA + ribosome released
            sim.evolveState();
            
            assertAllEqual(0, r.counts(r.matureIndexs(r.matureMRNAIndexs(1:2)), c.cytosolIndexs));
            assertEqual(0, r.counts(r.boundIndexs(r.matureTMRNAIndexs), c.cytosolIndexs));
            assertEqual(1, r.counts(r.matureIndexs(r.matureTMRNAIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.matureIndexs(pc.ribosome70SIndexs), c.cytosolIndexs));
            assertEqual(1, pc.counts(pc.matureIndexs(pc.ribosome30SIndexs), c.cytosolIndexs)+pc.counts(pc.matureIndexs(pc.getIndexs('RIBOSOME_30S_IF3')), c.cytosolIndexs));
            assertEqual(1, pc.counts(pc.matureIndexs(pc.ribosome50SIndexs), c.cytosolIndexs));
            assertEqual(1, pm.counts(pm.matureIndexs(pm.translationFactorIndexs), c.cytosolIndexs) + pm.counts(pm.boundIndexs(pm.translationFactorIndexs), c.cytosolIndexs));
            assertAllEqual(1, pc.counts(pc.matureIndexs(pc.translationFactorIndexs), c.cytosolIndexs) + pc.counts(pc.boundIndexs(pc.translationFactorIndexs), c.cytosolIndexs));
            assertEqual(1, numel(pol.abortedSequences));
            assertEqual(15, numel(pol.abortedSequences{1}));
            
            %add protein decay enzymes, evolve, and assert that polypeptide is degraded
            pm.counts(prtDcy.enzymeMonomerGlobalCompartmentIndexs) = 10;
            pc.counts(prtDcy.enzymeComplexGlobalCompartmentIndexs) = 10;
            
            sim.evolveState();
            assertEqual(0, numel(pol.abortedSequences));
            
            %% decay tmRNA
            %setup by with 1 stalled ribosome, decayRates = 0, decay enzymes = 0
            pm.counts(prtDcy.enzymeMonomerGlobalCompartmentIndexs) = 0;
            pc.counts(prtDcy.enzymeComplexGlobalCompartmentIndexs) = 0;
            
            r.counts(r.matureIndexs(r.matureMRNAIndexs(1)), c.cytosolIndexs) = 1;
            r.counts(r.boundIndexs(r.matureTMRNAIndexs), c.cytosolIndexs) = 1;
            r.counts(r.matureIndexs(r.matureTMRNAIndexs), c.cytosolIndexs) = 0;
            
            pc.counts(pc.matureIndexs(pc.getIndexs('RIBOSOME_30S_IF3')), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.ribosome30SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.ribosome50SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.ribosome70SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs) = 1;
            pm.counts(pm.matureIndexs(pm.translationFactorIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.translationFactorIndexs), c.cytosolIndexs) = 0;
            pm.counts(pm.boundIndexs(pm.translationFactorIndexs), c.cytosolIndexs) = 1;
            pc.counts(pc.boundIndexs(pc.translationFactorIndexs), c.cytosolIndexs) = 1;
            
            rib.states = rib.stalledValue;
            rib.boundMRNAs = 1;
            rib.mRNAPositions = 10;
            rib.tmRNAPositions = 5;
            pol.boundMRNAs = rib.boundMRNAs;
            pol.nascentMonomerLengths = rib.mRNAPositions;
            pol.proteolysisTagLengths = rib.tmRNAPositions;
            pol.abortedPolypeptides = zeros(0, 3);
            
            m.counts(:) = 0;
            m.counts(m.waterIndexs, c.cytosolIndexs) = 1e7;
            
            r.decayRates(:) = 0;
            r.decayRates(r.boundIndexs(r.matureTMRNAIndexs)) = 10;
            
            %evolve and assert that tmRNA degraded, mRNA + ribosome released
            sim.evolveState();
            
            assertEqual(1, r.counts(r.matureIndexs(r.matureMRNAIndexs(1)), c.cytosolIndexs));
            assertEqual(0, r.counts(r.boundIndexs(r.matureTMRNAIndexs), c.cytosolIndexs));
            assertEqual(0, r.counts(r.matureIndexs(r.matureTMRNAIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.matureIndexs(pc.ribosome70SIndexs), c.cytosolIndexs));
            assertEqual(1, pc.counts(pc.matureIndexs(pc.ribosome30SIndexs), c.cytosolIndexs)+pc.counts(pc.matureIndexs(pc.getIndexs('RIBOSOME_30S_IF3')), c.cytosolIndexs));
            assertEqual(1, pc.counts(pc.matureIndexs(pc.ribosome50SIndexs), c.cytosolIndexs));
            assertEqual(1, pm.counts(pm.matureIndexs(pm.translationFactorIndexs), c.cytosolIndexs) + pm.counts(pm.boundIndexs(pm.translationFactorIndexs), c.cytosolIndexs));
            assertAllEqual(1, pc.counts(pc.matureIndexs(pc.translationFactorIndexs), c.cytosolIndexs) + pc.counts(pc.boundIndexs(pc.translationFactorIndexs), c.cytosolIndexs));
            assertEqual(1, numel(pol.abortedSequences));
            assertEqual(15, numel(pol.abortedSequences{1}));
            
            %add protein decay enzymes, evolve, and assert that polypeptide is degraded
            pm.counts(prtDcy.enzymeMonomerGlobalCompartmentIndexs) = 10;
            pc.counts(prtDcy.enzymeComplexGlobalCompartmentIndexs) = 10;
            
            sim.evolveState();
            assertEqual(0, numel(pol.abortedSequences));
        end
        
        function testProteinDecayIntegration(this)
            sim = this.simulation;
            c = sim.compartment;
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            rib = sim.state('Ribosome');
            pol = sim.state('Polypeptide');
            chr = sim.state('Chromosome');
            rnaPol = sim.state('RNAPolymerase');
            trnscpt = sim.state('Transcript');
            ring = sim.state('FtsZRing');
            rnaDcy = sim.process('RNADecay');
            prtDcy = sim.process('ProteinDecay');
            ftsZPol = sim.process('FtsZPolymerization');
            transcription = sim.process('Transcription');
            
            %% decay ribosome
            %setup by with 1 stalled ribosome, decayRates = 0, decay enzymes = 0
            pc.counts(pc.matureIndexs(prtDcy.enzymeGlobalIndexs(prtDcy.enzymeIndexs_ftsHProtease)), :) = 0;
            
            r.counts(r.matureIndexs(r.matureMRNAIndexs(1:2)), c.cytosolIndexs) = 1;
            r.counts(r.boundIndexs(r.matureTMRNAIndexs), c.cytosolIndexs) = 1;
            r.counts(r.matureIndexs(r.matureTMRNAIndexs), c.cytosolIndexs) = 0;
            
            pc.counts(pc.matureIndexs(pc.ribosome30SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.ribosome50SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.ribosome70SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs) = 2;
            pm.counts(pm.matureIndexs(pm.translationFactorIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.translationFactorIndexs), c.cytosolIndexs) = 0;
            pm.counts(pm.boundIndexs(pm.translationFactorIndexs), c.cytosolIndexs) = 2;
            pc.counts(pc.boundIndexs(pc.translationFactorIndexs), c.cytosolIndexs) = 2;
            
            rib.states = [rib.stalledValue; rib.activeValue];
            rib.boundMRNAs = [1; 2];
            rib.mRNAPositions = [100; 0];
            rib.tmRNAPositions = [5; 0];
            pol.boundMRNAs = rib.boundMRNAs;
            pol.nascentMonomerLengths = rib.mRNAPositions;
            pol.proteolysisTagLengths = rib.tmRNAPositions;
            pol.abortedPolypeptides = zeros(0, 3);
            
            m.counts(:) = 0;
            m.counts(m.waterIndexs, c.cytosolIndexs) = 1e7;
            
            r.decayRates(:) = 0;
            pc.decayRates(:) = 0;
            pc.decayRates(pc.boundIndexs(pc.ribosome70SIndexs)) = 10;
            
            %evolve and assert that tmRNA degraded, mRNA + ribosome released
            sim.evolveState();
            
            assertAllEqual(1, r.counts(r.matureIndexs(r.matureMRNAIndexs(1:2)), c.cytosolIndexs));
            assertEqual(0, r.counts(r.boundIndexs(r.matureTMRNAIndexs), c.cytosolIndexs));
            assertEqual(1, r.counts(r.matureIndexs(r.matureTMRNAIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.matureIndexs(pc.ribosome70SIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.matureIndexs(pc.ribosome30SIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.matureIndexs(pc.ribosome50SIndexs), c.cytosolIndexs));
            assertEqual(2, pm.counts(pm.matureIndexs(pm.translationFactorIndexs), c.cytosolIndexs) + pm.counts(pm.boundIndexs(pm.translationFactorIndexs), c.cytosolIndexs));
            assertAllEqual(2, pc.counts(pc.matureIndexs(pc.translationFactorIndexs), c.cytosolIndexs) + pc.counts(pc.boundIndexs(pc.translationFactorIndexs), c.cytosolIndexs));
            assertEqual(1, numel(pol.abortedSequences));
            assertEqual(105, numel(pol.abortedSequences{1}));
            
            %add protein decay enzymes, evolve, and assert that polypeptide is degraded
            m.counts(m.atpIndexs, c.cytosolIndexs) = 1e6;
            m.counts(m.waterIndexs, c.cytosolIndexs) = 1e6;
            pc.counts(pc.matureIndexs(prtDcy.enzymeGlobalIndexs(prtDcy.enzymeIndexs_ftsHProtease)), :) = 1000;
            
            sim.evolveState();
            assertEqual(0, numel(pol.abortedSequences));
            
            %% decay translation elongation factor
            %setup by with active ribosome, decayRates = 0, decay enzymes = 0
            pc.counts(pc.matureIndexs(prtDcy.enzymeGlobalIndexs(prtDcy.enzymeIndexs_ftsHProtease)), :) = 0;
            
            r.counts(r.matureIndexs(r.matureMRNAIndexs(1)), c.cytosolIndexs) = 1;
            r.counts(r.boundIndexs(r.matureTMRNAIndexs), c.cytosolIndexs) = 0;
            r.counts(r.matureIndexs(r.matureTMRNAIndexs), c.cytosolIndexs) = 0;
            
            pc.counts(pc.matureIndexs(pc.ribosome30SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.ribosome50SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.ribosome70SIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs) = 1;
            pm.counts(pm.matureIndexs(pm.translationFactorIndexs), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.translationFactorIndexs), c.cytosolIndexs) = 0;
            pm.counts(pm.boundIndexs(pm.translationFactorIndexs), c.cytosolIndexs) = 1;
            pc.counts(pc.boundIndexs(pc.translationFactorIndexs), c.cytosolIndexs) = 1;
            
            rib.states = rib.activeValue;
            rib.boundMRNAs = 1;
            rib.mRNAPositions = 50;
            rib.tmRNAPositions = 0;
            pol.boundMRNAs = rib.boundMRNAs;
            pol.nascentMonomerLengths = rib.mRNAPositions;
            pol.proteolysisTagLengths = rib.tmRNAPositions;
            pol.abortedPolypeptides = zeros(0, 3);
            
            m.counts(:) = 0;
            m.counts(m.ntpIndexs([1 3]), c.cytosolIndexs) = 1e6;
            m.counts(m.waterIndexs, c.cytosolIndexs) = 1e6;
            
            r.decayRates(:) = 0;
            pc.decayRates(:) = 0;
            pm.decayRates(:) = 0;
            pm.decayRates(pm.boundIndexs(pm.translationFactorIndexs)) = 10;
            
            %evolve and assert that tmRNA degraded, mRNA + ribosome released
            sim.evolveState();
            
            assertEqual(1, r.counts(r.matureIndexs(r.matureMRNAIndexs(1)), c.cytosolIndexs));
            assertEqual(0, r.counts(r.boundIndexs(r.matureTMRNAIndexs), c.cytosolIndexs));
            assertEqual(0, r.counts(r.matureIndexs(r.matureTMRNAIndexs), c.cytosolIndexs));
            assertEqual(1, pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.matureIndexs(pc.ribosome70SIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.matureIndexs(pc.ribosome30SIndexs), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.matureIndexs(pc.ribosome50SIndexs), c.cytosolIndexs));
            assertEqual(0, pm.counts(pm.matureIndexs(pm.translationFactorIndexs), c.cytosolIndexs) + pm.counts(pm.boundIndexs(pm.translationFactorIndexs), c.cytosolIndexs));
            assertAllEqual(1, pc.counts(pc.matureIndexs(pc.translationFactorIndexs), c.cytosolIndexs) + pc.counts(pc.boundIndexs(pc.translationFactorIndexs), c.cytosolIndexs));
            assertEqual(0, numel(pol.abortedSequences));
            assertEqual(rib.activeValue, rib.states);
            
            %% decay FtsZ
            pm.counts(ftsZPol.enzymeMonomerGlobalCompartmentIndexs) = 0;
            pc.counts(ftsZPol.enzymeComplexGlobalCompartmentIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.getIndexs('MG_224_9MER_GTP')), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.getIndexs('MG_224_9MER_GDP')), c.cytosolIndexs) = 2;
            pc.counts(pc.boundIndexs(pc.getIndexs('MG_224_9MER_GTP')), c.cytosolIndexs) = 3;
            pc.counts(pc.boundIndexs(pc.getIndexs('MG_224_9MER_GDP')), c.cytosolIndexs) = 4;
            ring.numEdgesOneStraight = 1;
            ring.numEdgesTwoStraight = 1;
            ring.numEdgesTwoBent = 2;
            ring.numResidualBent = 0;
            
            m.counts(:) = 0;
            m.counts(m.ntpIndexs([1 3]), c.cytosolIndexs) = 1e6;
            m.counts(m.waterIndexs, c.cytosolIndexs) = 1e6;
            
            pc.decayRates(:) = 0;
            pc.decayRates(pc.boundIndexs(pc.getIndexs('MG_224_9MER_GTP'))) = 10;
            pc.decayRates(pc.boundIndexs(pc.getIndexs('MG_224_9MER_GDP'))) = 10;
            sim.evolveState();
            
            assertEqual(0, pc.counts(pc.matureIndexs(pc.getIndexs('MG_224_9MER_GTP')), c.cytosolIndexs));
            assertEqual(2, pc.counts(pc.matureIndexs(pc.getIndexs('MG_224_9MER_GDP')), c.cytosolIndexs));
            assertEqual(3, pc.counts(pc.boundIndexs(pc.getIndexs('MG_224_9MER_GTP')), c.cytosolIndexs));
            assertEqual(4, pc.counts(pc.boundIndexs(pc.getIndexs('MG_224_9MER_GDP')), c.cytosolIndexs));
            assertEqual(1, ring.numEdgesOneStraight);
            assertEqual(1, ring.numEdgesTwoStraight);
            assertEqual(2, ring.numEdgesTwoBent);
            assertEqual(0, ring.numResidualBent);
            
            %% decay RNA polymerase
            chr.initialize();
            sim.process('DNARepair').initializeState();
            
            pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE')), c.cytosolIndexs) = 0;
            pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')), c.cytosolIndexs) = 0;
            pc.counts(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE')), c.cytosolIndexs) = 1;
            pc.counts(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')), c.cytosolIndexs) = 0;
            
            pm.counts(rnaDcy.enzymeMonomerGlobalCompartmentIndexs(rnaDcy.enzymeIndexs_ribonucleaseR)) = 0;
            
            m.counts(:) = 0;
            m.counts(m.ntpIndexs([1 3]), c.cytosolIndexs) = 1e6;
            m.counts(m.waterIndexs, c.cytosolIndexs) = 1e6;
            
            iTU = find(trnscpt.transcriptionUnitDirections == 1, 1, 'first');
            trnscpt.boundTranscriptionUnits = iTU;
            rnaPol.states = 10;
            rnaPol.positionStrands = [...
                trnscpt.transcriptionUnitFivePrimeCoordinates(iTU) + 10 - 1 ...
                chr.transcriptionUnitStrands(iTU)];
            trnscpt.boundTranscriptProgress = rnaPol.states;
            trnscpt.boundTranscriptChromosome = 1;
            transcription.bindProteinToChromosome(rnaPol.positionStrands, transcription.enzymeIndexs_rnaPolymerase);
            transcription.rnaPolymeraseElongationRate = 0;
            trnscpt.abortedTranscripts = zeros(0, 2);
            
            pc.decayRates(:) = 0;
            pc.decayRates(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE'))) = 10;
            
            %evolve and assert that RNA polymerase degraded, transcript
            %released
            sim.evolveState();
            
            assertEqual(0, pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE')), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE')), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')), c.cytosolIndexs));
            assertEqual(0, nnz(chr.complexBoundSites == pc.getIndexs('RNA_POLYMERASE')));
            assertEqual(rnaPol.notExistValue, rnaPol.states);
            assertEqual([0 0], rnaPol.positionStrands);
            assertEqual(0, trnscpt.boundTranscriptProgress);
            assertEqual(0, trnscpt.boundTranscriptionUnits);
            assertEqual(0, trnscpt.boundTranscriptChromosome);
            assertEqual(1, numel(trnscpt.abortedSequences));
            assertIn(numel(trnscpt.abortedSequences{1}), [9 trnscpt.transcriptionUnitLengths(iTU)+1]);
            
            %add RNA decay enzymes and check that aborted transcript is
            %degraded
            pm.counts(rnaDcy.enzymeMonomerGlobalCompartmentIndexs(rnaDcy.enzymeIndexs_ribonucleaseR)) = 10;
            
            sim.evolveState();
            
            assertEqual(0, numel(trnscpt.abortedSequences));
            
            %% decay other chromosomally bound protein
            monIdx = find(pm.counts(pm.boundIndexs, c.cytosolIndexs) == 0, 10, 'first');
            monIdx = monIdx(end);
            cpxIdx = find(pc.counts(pc.boundIndexs, c.cytosolIndexs) == 0, 10, 'first');
            cpxIdx = cpxIdx(end);
            
            pm.counts(pm.boundIndexs(monIdx), c.cytosolIndexs) = 1;
            pc.counts(pc.boundIndexs(cpxIdx), c.cytosolIndexs) = 1;
            
            chr.monomerBoundSites(1000, 1) = monIdx;
            chr.complexBoundSites(2000, 1) = cpxIdx;
            
            pm.decayRates(:) = 0;
            pc.decayRates(:) = 0;
            pm.decayRates(pm.boundIndexs(monIdx)) = 1;
            pc.decayRates(pc.boundIndexs(cpxIdx)) = 1;
            
            m.counts(:) = 0;
            m.counts(m.ntpIndexs([1 3]), c.cytosolIndexs) = 1e6;
            m.counts(m.waterIndexs, c.cytosolIndexs) = 1e6;
            
            sim.evolveState();
            
            assertEqual(0, pm.counts(pm.boundIndexs(monIdx), c.cytosolIndexs));
            assertEqual(0, pc.counts(pc.boundIndexs(cpxIdx), c.cytosolIndexs));
            assertEqual(0, nnz(chr.monomerBoundSites == monIdx));
            assertEqual(0, nnz(chr.complexBoundSites == cpxIdx));
        end
    end
    
    methods
        %Test that process uses rand streams correctly
        %- only changes state of process (and possibly state) rand
        %  stream(s)
        %- doesn't change state of global rand stream
        function testRandStreamUsage(this)
            lengthSec = 100;
            seed = 10;
            
            %% hold state of default stream
            defaultStream = RandStream.getDefaultStream;
            defaultStream.reset(1);
            
            defaultStreamInit = struct();
            defaultStreamInit.Type = defaultStream.Type;
            defaultStreamInit.Seed = defaultStream.Seed;
            defaultStreamInit.NumStreams = defaultStream.NumStreams;
            defaultStreamInit.StreamIndex = defaultStream.StreamIndex;
            defaultStreamInit.State = defaultStream.State;
            defaultStreamInit.Substream = defaultStream.Substream;
            defaultStreamInit.RandnAlg = defaultStream.RandnAlg;
            defaultStreamInit.Antithetic = defaultStream.Antithetic;
            defaultStreamInit.FullPrecision = defaultStream.FullPrecision;
            
            %% simulate
            %run-1
            this.setUp();
            sim = this.simulation;
            sim.applyOptions('lengthSec', lengthSec, 'seed', seed);
            streams1 = cell(1 + numel(sim.states) + numel(sim.processes), 1);
            streams1{1} = sim.getForTest('randStream');
            for i = 1:numel(sim.states)
                streams1{i+1} = sim.states{i}.randStream;
            end
            for i = 1:numel(sim.processes)
                streams1{i+1+numel(sim.states)} = sim.processes{i}.randStream;
            end
            sim.run();
            timeCourses1 = sim.getTimeCourses();
            this.tearDown();
            checkDefualtStreamUnchanged(defaultStream, defaultStreamInit);
            
            %run-2
            this.setUp();
            sim = this.simulation;
            sim.applyOptions('lengthSec', lengthSec, 'seed', seed);
            streams2 = cell(1 + numel(sim.states) + numel(sim.processes), 1);
            streams2{1} = sim.getForTest('randStream');
            for i = 1:numel(sim.states)
                streams2{i+1} = sim.states{i}.randStream;
            end
            for i = 1:numel(sim.processes)
                streams2{i+1+numel(sim.states)} = sim.processes{i}.randStream;
            end
            sim.run();
            timeCourses2 = sim.getTimeCourses();
            this.tearDown();
            checkDefualtStreamUnchanged(defaultStream, defaultStreamInit);
            
            %% assert state of simulation, state, and process rand streams updated the same in both runs
            for i = 1:numel(streams1)
                assertTrue(streams1{i} == streams2{i});
            end
            assertEqual(timeCourses1, timeCourses2);

            %% assert state of default stream unchanged
            checkDefualtStreamUnchanged(defaultStream, defaultStreamInit);
            function checkDefualtStreamUnchanged(defaultStream, defaultStreamInit)
                isequal(defaultStream, RandStream.getDefaultStream);
                assertEqual(defaultStreamInit.Type, defaultStream.Type);
                assertEqual(defaultStreamInit.Seed, defaultStream.Seed);
                assertEqual(defaultStreamInit.NumStreams, defaultStream.NumStreams);
                assertEqual(defaultStreamInit.StreamIndex, defaultStream.StreamIndex);
                assertEqual(defaultStreamInit.State, defaultStream.State);
                assertEqual(defaultStreamInit.Substream, defaultStream.Substream);
                assertEqual(defaultStreamInit.RandnAlg, defaultStream.RandnAlg);
                assertEqual(defaultStreamInit.Antithetic, defaultStream.Antithetic);
                assertEqual(defaultStreamInit.FullPrecision, defaultStream.FullPrecision);
            end
        end
    end
end
