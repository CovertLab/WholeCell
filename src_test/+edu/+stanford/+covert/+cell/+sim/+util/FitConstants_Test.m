%FitConstants test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updaetd 2/3/2011
classdef FitConstants_Test < TestCase
    properties
        fitter
        simulation
        warningStatus
    end
    
    methods
        function this = FitConstants_Test(name)
            this = this@TestCase(name);
        end
        
        function setUp(this)
            %turn warning off
            this.warningStatus = warning('query', 'WholeCell:warning');
            warning('off', 'WholeCell:warning');
            
            %create simulation, fitter
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
            sim.applyOptions('lengthSec', 10);
            this.simulation = sim;
            this.fitter = edu.stanford.covert.cell.sim.util.FitConstants(sim, struct('verbosity', 0));
        end
        
        function tearDown(this)
            %reset warning state
            warning(this.warningStatus.state, 'WholeCell:warning');
        end
    end
    
    methods
        function testConstants(this)
            fitter = this.fitter;
            sim = this.simulation;
            g = sim.gene;
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            %base counts calculated correctly
            assertEqual(r.lengths(r.matureIndexs), sum(fitter.getRnaNMPCounts(), 2));
            assertEqual(pm.lengths(pm.matureIndexs), sum(fitter.getMonomerAACounts(), 2));            
            
            %gyrase genes in transcription unit, and same stoichiometry in
            %complex
            assertEqual(1, nnz(any(r.matureRNAGeneComposition(ismember(g.wholeCellModelIDs, {'MG_003', 'MG_004'}), :),1)));
            pcComp = sum(pc.proteinComplexComposition(ismember(g.wholeCellModelIDs, {'MG_003', 'MG_004'}), :, :),3);
            assertEqual(1, nnz(any(pcComp, 1)));
            assertEqual(2, numel(unique(pcComp)));
            
            %topoisomerase I is monomer
            assertAllEqual(0, pc.proteinComplexComposition(ismember(g.wholeCellModelIDs, 'MG_122'), :, :));
        end
        
        function disabled_testFitAnalytically(this)
            %references
            fitter = this.fitter;
            
            %experimentally observed values
            paramVec0 = fitter.initializeFittedConstants();
            
            %fit
            fitter.method = 'analytic';
            fitter.run();
            paramVecAnal = fitter.constructParameterVector(paramVec0);
            this.areJointParameterConstraintsSatisfied();
            
            %check fit is better than heuristic method
            fitter.method = 'heuristic';
            fitter.run();
            paramVecHeur = fitter.constructParameterVector(paramVec0);
            this.areJointParameterConstraintsSatisfied();
            
            assertIn(sum((paramVecAnal{1} - paramVec0{1}).^2), [0 sum((paramVecHeur{1} - paramVec0{1}).^2) * (1+5e-3)]);
        end
        
        function testFitHeuristically(this)
            %fit
            this.fitter.method = 'heuristic';
            this.fitter.run();
            
            %check joint parameter constraints are satisfied
            this.areJointParameterConstraintsSatisfied();
        end
        
        function areJointParameterConstraintsSatisfied(this)
            %classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            % references
            fitter = this.fitter;
            sim = this.simulation;
            g = sim.gene;
            t = sim.state('Time');
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            ring = sim.state('FtsZRing');
            geom = sim.state('Geometry');
            mass = sim.state('Mass');
            sc = sim.process('DNASupercoiling');
            
            %indices
            iCyt = sim.compartment.cytosolIndexs;
            
            %constants
            rnaMWs = r.molecularWeights(r.matureIndexs) / ConstantUtil.nAvogadro;
            paramVec = fitter.constructParameterVectorFromSimulation();
            [rnaExp, ~, ~, rnaWtFracs, rnaDecayRates, bmComp, bmProd, byProd, unaccEComp] = ...
                fitter.extractParameterVector(paramVec);
            monExp = (r.matureRNAGeneComposition(g.mRNAIndexs, :) * rnaExp) ./ ...
                (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
            monExp = monExp / sum(monExp);
            nmpComp = fitter.getRnaNMPCounts()' * rnaExp;
            nmpComp = nmpComp / sum(nmpComp);
            aaComp = fitter.getMonomerAACounts()' * monExp;
            aaComp = aaComp / sum(aaComp);
            paramVec = fitter.constructParameterVector(rnaExp, nmpComp, aaComp, ...
                rnaWtFracs, rnaDecayRates, bmComp, bmProd, byProd, unaccEComp);
            
            %RNA expression matches RNA production + half lives
            rnaProd = sim.process('Transcription').computeRNAPolymeraseTUBindingProbabilities(false);
            rnaProd = rnaProd(:, 1) .* fitter.calcEffectiveMeanTranscriptionUnitCopyNumbers();
            tmpRnaExp = (r.nascentRNAMatureRNAComposition * rnaProd) ./ ...
                (log(2) / t.cellCycleLength + r.decayRates(r.matureIndexs));
            tmpRnaExp = tmpRnaExp / sum(tmpRnaExp);
            assertElementsAlmostEqual(rnaExp, tmpRnaExp, 'relative', 1e-8, 0);
            
            %RNA production matches NTP production
            invMat = ComputationUtil.invertCompositionMatrix(r.nascentRNAMatureRNAComposition);
            [~, ~, ~, ~, ~, totRnas, ~, ~, rnaDecays] = fitter.calcMacromolecularCounts(paramVec);
            ntpProd = r.baseCounts(r.nascentIndexs, m.nmpIndexs)' * invMat * (totRnas + rnaDecays);
            assertElementsAlmostEqual(...
                bmProd(m.ntpIndexs([2 4])), ...
                ntpProd([2 4]), ...
                'relative', 15e-2, 0);
            
            %assert all RNA expression is in mature state
            assertElementsAlmostEqual(1, sum(rnaExp), 'relative', 1e-12, 0);
            
            %production rates of genes in transcription units is equal
            assertElementsAlmostEqual(rnaExp, (r.nascentRNAMatureRNAComposition * ...
                ComputationUtil.invertCompositionMatrix(r.nascentRNAMatureRNAComposition) * ...
                (rnaExp .* (log(2) / t.cellCycleLength + rnaDecayRates))) ./ (log(2) / t.cellCycleLength + rnaDecayRates), ...
                'relative', 1e-12, 0);
            assertElementsAlmostEqual(totRnas + rnaDecays, ...
                r.nascentRNAMatureRNAComposition * (invMat * (totRnas + rnaDecays)), ...
                'relative', 1e-8);
            
            %RNA expression matches weight fractions (mrstRNA)
            assertElementsAlmostEqual(1, sum(r.expectedWeightFractions), 'relative', 1e-8, 0);
            assertElementsAlmostEqual(1, sum([...
                rnaMWs(r.matureMRNAIndexs)' * rnaExp(r.matureMRNAIndexs)
                rnaMWs(r.matureRibosomalRRNAIndexs) .* rnaExp(r.matureRibosomalRRNAIndexs)
                rnaMWs([r.matureSRNAIndexs; r.matureTRNAIndexs])' * rnaExp([r.matureSRNAIndexs; r.matureTRNAIndexs])] / ...
                (rnaMWs' * rnaExp)), ...
                'relative', 1e-8, 0);
            
            %NMP composition matches expression, sequences
            tmp = sum(bmComp(m.nmpIndexs, :), 2);
            assertElementsAlmostEqual(tmp / sum(tmp), nmpComp, 'relative', 2e-2, 0);
            
            tmp = bmComp(m.nmpIndexs, iCyt);
            tmp(m.nmpIndexs == m.getIndexs('AMP')) = tmp(m.nmpIndexs == m.getIndexs('AMP')) + bmComp(m.getIndexs('m62AMP'), iCyt);
            tmp(m.nmpIndexs == m.getIndexs('CMP')) = tmp(m.nmpIndexs == m.getIndexs('CMP')) + bmComp(m.getIndexs('k2CMP'), iCyt);
            tmp(m.nmpIndexs == m.getIndexs('GMP')) = tmp(m.nmpIndexs == m.getIndexs('GMP')) + bmComp(m.getIndexs('GmMP'), iCyt);
            tmp(m.nmpIndexs == m.getIndexs('GMP')) = tmp(m.nmpIndexs == m.getIndexs('GMP')) + bmComp(m.getIndexs('m1GMP'), iCyt);
            tmp(m.nmpIndexs == m.getIndexs('GMP')) = tmp(m.nmpIndexs == m.getIndexs('GMP')) + bmComp(m.getIndexs('m2GMP'), iCyt);
            tmp(m.nmpIndexs == m.getIndexs('GMP')) = tmp(m.nmpIndexs == m.getIndexs('GMP')) + bmComp(m.getIndexs('m7GMP'), iCyt);
            tmp(m.nmpIndexs == m.getIndexs('UMP')) = tmp(m.nmpIndexs == m.getIndexs('UMP')) + bmComp(m.getIndexs('cmnm5s2UMP'), iCyt);
            tmp(m.nmpIndexs == m.getIndexs('UMP')) = tmp(m.nmpIndexs == m.getIndexs('UMP')) + bmComp(m.getIndexs('PSIURIMP'), iCyt);
            tmp(m.nmpIndexs == m.getIndexs('UMP')) = tmp(m.nmpIndexs == m.getIndexs('UMP')) + bmComp(m.getIndexs('s2UMP'), iCyt);
            tmp(m.nmpIndexs == m.getIndexs('UMP')) = tmp(m.nmpIndexs == m.getIndexs('UMP')) + bmComp(m.getIndexs('s4UMP'), iCyt);
            tmp(m.nmpIndexs == m.getIndexs('UMP')) = tmp(m.nmpIndexs == m.getIndexs('UMP')) + bmComp(m.getIndexs('UmMP'), iCyt);
            assertElementsAlmostEqual(tmp / sum(tmp), nmpComp, 'relative', 20e-2, 0);
            
            assertElementsAlmostEqual(...
                this.getNucleotideCompositionMatrix(m) * bmComp, ...
                this.getNucleotideCompositionMatrix(m) * (bmProd - byProd), ...
                'relative', 1e-3);
            
            %amino acid composition matches expression, sequences
            tmp = sum(bmComp(m.aminoAcidIndexs, :), 2);
            assertElementsAlmostEqual(tmp / sum(tmp), aaComp, 'relative', 1e-2, 0);
            
            tmp = bmComp(m.aminoAcidIndexs, iCyt);
            tmp(m.aminoAcidIndexs == m.getIndexs('CYS')) = tmp(m.aminoAcidIndexs == m.getIndexs('CYS')) + bmComp(m.getIndexs('diacylglycerolCys'), iCyt);
            tmp(m.aminoAcidIndexs == m.getIndexs('LYS')) = tmp(m.aminoAcidIndexs == m.getIndexs('LYS')) + bmComp(m.getIndexs('LIPOYLLYS'), iCyt);
            tmp(m.aminoAcidIndexs == m.getIndexs('SER')) = tmp(m.aminoAcidIndexs == m.getIndexs('SER')) + bmComp(m.getIndexs('pSER'), iCyt);
            tmp(m.aminoAcidIndexs == m.getIndexs('THR')) = tmp(m.aminoAcidIndexs == m.getIndexs('THR')) + bmComp(m.getIndexs('pTHR'), iCyt);
            tmp(m.aminoAcidIndexs == m.getIndexs('TYR')) = tmp(m.aminoAcidIndexs == m.getIndexs('TYR')) + bmComp(m.getIndexs('pTYR'), iCyt);
            assertElementsAlmostEqual(tmp / sum(tmp), aaComp, 'relative', 5e-2, 0);
            
            %RNA, protein expression matches cell cycle length (> calculated min expression)
            nonLargeCpxMonIdxs = find(~any(any(pc.proteinComplexComposition(g.mRNAIndexs, pc.getIndexs({'MG_271_272_273_274_192MER'; 'MG_271_272_273_274_192MER'; 'MG_398_399_400_401_402_403_404_405_22MER'}), :), 3), 2));
            sim.process('Metabolism').tolerance = 1e-1;
            [~, ~, ~, ~, ~, rnaCnt, monCnt] = fitter.calcMacromolecularCounts(paramVec);
            [~, ~, ~, ~, minRNACnt, minMonCnt, maxRNACnt, maxMonCnt] = fitter.calcResourceRequirements(paramVec);
            assertIn(min(rnaCnt ./ minRNACnt), [1 - 2e-3 Inf]);
            assertIn(max(rnaCnt ./ maxRNACnt), [0 1]);
            assertIn(min(monCnt ./ minMonCnt), [1 - 2e-2 Inf]); 
            assertIn(max(monCnt(nonLargeCpxMonIdxs) ./ maxMonCnt(nonLargeCpxMonIdxs)), [0 1 + 2e-2]);
            assertIn(min(monCnt(minMonCnt > 0)), [(1-fitter.tolerance) * pm.minimumAverageExpression Inf]);
            assertIn(min(monCnt(setdiff(1:end, [fitter.monomerIndex_DnaA; fitter.monomerIndex_FtsZ])) ./ ...
                max(sum(pc.proteinComplexComposition(g.mRNAIndexs(...
                setdiff(1:end, [fitter.monomerIndex_DnaA; fitter.monomerIndex_FtsZ])), :, :), ...
                3), [], 2)), [(1-fitter.tolerance) * pc.minimumAverageExpression Inf]);
            assertIn(min(rnaCnt(r.matureSRNAIndexs(any(any(pc.proteinComplexComposition(g.sRNAIndexs, :, :), 3), 2)))), [12 Inf]);
            
            %FtsZ, DnaA expression are preserved so that replication initiation and cytokinesis
            %durations don't need to be refit/recalculated
            assertElementsAlmostEqual(fitter.dnaACnt, ...
                monCnt(fitter.monomerIndex_DnaA), ...
                'absolute', 1);
            
            numEdges = ring.calcNumEdges(...
                geom.calculateWidth(mass.cellInitialDryWeight / (1 - mass.fractionWetWeight) / geom.density), ...
                ring.filamentLengthInNm);
            assertElementsAlmostEqual(fitter.ftsZCnt, ...
                monCnt(fitter.monomerIndex_FtsZ), ...
                'absolute', 5);
            assertElementsAlmostEqual(ceil(4.0 / exp(log(2) * (1-t.cytokinesisDuration / t.cellCycleLength)) * numEdges * ring.numFtsZSubunitsPerFilament), ...
                monCnt(fitter.monomerIndex_FtsZ), ...
                'absolute', 5);
            
            %topoisomerase I / gyrase activity balanced
            assertElementsAlmostEqual(...
                monCnt(fitter.monomerIndex_topoisomeraseI) * -sc.topoIDeltaLK * sc.topoIActivityRate, ...
                min(monCnt(fitter.monomerIndex_gyrase) / fitter.gryStoichiometry) * sc.gyraseDeltaLK * sc.gyraseActivityRate, ...
                'relative', 1e-2);
            
            %test objective function
            [~, ~, ~, ~, ~, paramVec{1}(fitter.rnaExpIdxs)] = fitter.calcMacromolecularCounts(paramVec);
            H = eye(numel(paramVec{1}));
            [val, gradient, hessian] = fitter.objFunc(paramVec{1}, paramVec{1}, H);
            assertEqual(0, val + paramVec{1}' * paramVec{1});
            assertEqual(size(paramVec{1}), size(gradient));
            assertEqual([numel(paramVec{1}) numel(paramVec{1})], size(hessian));
            
            %all linear constraints satisfied
            [Aineq, bineq, Aeq, beq] = fitter.linConstraintFunc(paramVec);
            assertEqual(true, all(Aineq * paramVec{1} <= bineq + 5e-2));
            assertElementsAlmostEqual(beq, Aeq * paramVec{1}, 'relative', 1e-10);
            
            %all non-linear constraints satisfied
            [c, ceq, gradC, gradCeq] = fitter.nlinConstraintFunc(paramVec{1});
            assertEqual([], c);
            assertEqual([], gradC);
            assertElementsAlmostEqual(zeros(size(ceq))', ceq', 'absolute', 1e-2);
            assertEqual([numel(paramVec{1}) numel(ceq)], size(gradCeq));
            
            %check sizes of hessian and hessian * Y
            assertEqual(size(H), size(fitter.hessFunc(paramVec{1}, struct('ineqnonlin', [], 'eqnonlin', ones(size(ceq))), H)));
            assertEqual([size(H, 1) 2], size(fitter.hessMultFunc(H, ones(size(H, 2), 2))));
        end
        
        function testCalcResourceRequirements(this)
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %references
            sim = this.simulation;
            g = sim.gene;
            mass = sim.state('Mass');
            t = sim.state('Time');
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            %indices
            iCyt = sim.compartment.cytosolIndexs;
            iMem = sim.compartment.membraneIndexs;
            
            %constants
            metMWs = m.molecularWeights;
            
            %calculate biomass composition, maintenance energy, and byproducts
            this.fitter.method = 'heuristic';
            this.fitter.run();
            
            bmComp = m.biomassComposition;
            bmProd = m.biomassProduction;
            byProd = m.byproducts;
            
            %initialize
            c.initialize();
            sim.process('DNARepair').initializeState();
            
            %% assertions
            %all positive reals, correct size
            validateattributes(bmComp, {'numeric'}, {'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(bmComp(~any(pc.baseCounts < 0, 1), :), {'numeric'}, {'nonnegative', 'real'});
            validateattributes(bmProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(byProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            
            %metabolites all in cytosol and membrane
            assertAllEqual(0, bmComp(:, setdiff(1:end, [iCyt; iMem])));
            assertAllEqual(0, bmProd(:, setdiff(1:end, [iCyt; iMem])));
            assertAllEqual(0, byProd(:, setdiff(1:end, [iCyt; iMem])));
            
            %hydrophobic metabolites not in cytosol
            assertAllEqual(0, bmComp(m.hydrophobicIndexs, iCyt));
            assertAllEqual(0, bmProd(m.hydrophobicIndexs, iCyt));
            assertAllEqual(0, byProd(m.hydrophobicIndexs, iCyt));
            
            %hydrophilic metabolites not in membrane
            assertAllEqual(0, bmComp(setdiff(1:end, m.hydrophobicIndexs), iMem));
            assertAllEqual(0, bmProd(setdiff(1:end, m.hydrophobicIndexs), iMem));
            assertAllEqual(0, byProd(setdiff(1:end, m.hydrophobicIndexs), iMem));
            
            %correct biomass composition weight fractions
            baseCount = getBaseCounts(c.sequence);
            modBaseCount = getBaseCounts(c.sequence);
            modBaseCount(1) = modBaseCount(1) - 2 * size(sim.process('DNARepair').RM_MunI_RecognitionSites, 1);
            modBaseCount(5) = 2 * size(sim.process('DNARepair').RM_MunI_RecognitionSites, 1);
            ntpCorrFactor = 2 * exp(-log(2) * (t.cellCycleLength - t.cytokinesisDuration) / t.cellCycleLength);
            chrWt = modBaseCount' * (m.molecularWeights([m.dnmpIndexs; m.getIndexs('m6dAMP')]) - (ConstantUtil.elements.O + ConstantUtil.elements.H));
            corrFactor = ((ntpCorrFactor - 1) * baseCount' * m.molecularWeights(m.dntpIndexs) + chrWt) /  chrWt;
            
            assertElementsAlmostEqual(c.dryWeight / mass.cellInitialDryWeight * corrFactor, ...
                mass.dryWeightFractionDNA, 'relative', 20e-2, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionDNA, ...
                bmComp([m.dnmpIndexs; m.getIndexs({'m6dAMP'}); m.dntpIndexs], iCyt)' / ConstantUtil.nAvogadro * ...
                [metMWs([m.dnmpIndexs; m.getIndexs({'m6dAMP'})]) - (ConstantUtil.elements.O + ConstantUtil.elements.H); metMWs(m.dntpIndexs)] / ...
                mass.cellInitialDryWeight, 'relative', 10e-2, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionRNA, ...
                bmComp(m.nmpIndexs, iCyt)' / ConstantUtil.nAvogadro * ...
                (metMWs(m.nmpIndexs)-(ConstantUtil.elements.O + ConstantUtil.elements.H)) / ...
                mass.cellInitialDryWeight, 'relative', 20e-2, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionProtein, ...
                bmComp(m.aminoAcidIndexs, iCyt)' / ConstantUtil.nAvogadro * ...
                (metMWs(m.aminoAcidIndexs)-(ConstantUtil.elements.O +  2 * ConstantUtil.elements.H)) / ...
                mass.cellInitialDryWeight, 'relative', 10e-2, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionLipid, ...
                bmComp(m.lipidIndexs, iMem)' / ConstantUtil.nAvogadro * ...
                metMWs(m.lipidIndexs) / mass.cellInitialDryWeight, 'relative', 10e-2, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionPolyamine, ...
                bmComp(m.polyamineIndexs, iCyt)' / ConstantUtil.nAvogadro * ...
                metMWs(m.polyamineIndexs) / mass.cellInitialDryWeight, 'relative', 10e-2, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionCarbohydrate, ...
                bmComp(m.carbohydrateIndexs, iCyt)' / ConstantUtil.nAvogadro * ...
                metMWs(m.carbohydrateIndexs) / mass.cellInitialDryWeight, 'relative', 10e-2, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionVitamin, ...
                bmComp(m.vitaminIndexs, iCyt)' / ConstantUtil.nAvogadro * ...
                metMWs(m.vitaminIndexs) / mass.cellInitialDryWeight, 'relative', 10e-2, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionIon, ...
                bmComp(m.ionIndexs, iCyt)' / ConstantUtil.nAvogadro * ...
                metMWs(m.ionIndexs) / mass.cellInitialDryWeight, 'relative', 50e-2, 0);
            
            %correct biomass normalization
            assertElementsAlmostEqual(1, ...
                + mass.dryWeightFractionRNA ...
                + mass.dryWeightFractionDNA ...
                + mass.dryWeightFractionProtein ...
                + mass.dryWeightFractionLipid ...
                + mass.dryWeightFractionPolyamine ...
                + mass.dryWeightFractionCarbohydrate ...
                + mass.dryWeightFractionVitamin ...
                + mass.dryWeightFractionIon ...
                + mass.dryWeightFractionNucleotide, ...
                'relative', 1e-8);
            
            effMetMWs = m.molecularWeights;
            effMetMWs(m.dnmpIndexs) = effMetMWs(m.dnmpIndexs) - (ConstantUtil.elements.O + ConstantUtil.elements.H);
            effMetMWs(m.nmpIndexs) = effMetMWs(m.nmpIndexs) - (ConstantUtil.elements.O + ConstantUtil.elements.H);
            effMetMWs(m.aminoAcidIndexs) = effMetMWs(m.aminoAcidIndexs) - (ConstantUtil.elements.O + 2 * ConstantUtil.elements.H);
            assertElementsAlmostEqual(mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), ...
                sum(bmComp, 2)' * effMetMWs / ConstantUtil.nAvogadro, ...
                'relative', 1e-2, 0);
            
            assertElementsAlmostEqual(sum(bmProd - byProd, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 5e-2, 0);
            
            %correct biomass production, byproducts -- DNA
            assertElementsAlmostEqual(bmComp(m.dnmpIndexs, iCyt) + [bmComp(m.getIndexs({'m6dAMP'}), iCyt); 0; 0; 0] + bmComp(m.dntpIndexs, iCyt), ...
                bmProd(m.dntpIndexs, iCyt), 'relative', 1e-2, 0);
            
            %correct biomass production, byproducts -- RNA
            invMat = edu.stanford.covert.util.ComputationUtil.invertCompositionMatrix(r.nascentRNAMatureRNAComposition);
            rnas0 = mass.dryWeightFractionRNA * mass.cellInitialDryWeight * r.expression(r.matureIndexs) / ...
                (r.expression(r.matureIndexs)' * r.molecularWeights(r.matureIndexs) / ConstantUtil.nAvogadro);
            rnaProd = rnas0 + rnas0 .* r.decayRates(r.matureIndexs) * t.cellCycleLength / log(2);
            decayedNMPs = r.baseCounts(r.nascentIndexs, m.nmpIndexs)' * invMat * ...
                (rnas0 .* r.decayRates(r.matureIndexs)) * t.cellCycleLength / log(2);
            prodNMPs = r.baseCounts(r.nascentIndexs, m.nmpIndexs)' * invMat * ...
                (rnas0 + (rnas0 .* r.decayRates(r.matureIndexs)) * ...
                t.cellCycleLength / log(2));
            
            assertElementsAlmostEqual(...
                r.baseCounts(r.nascentIndexs, m.nmpIndexs([2 4]))' * invMat * rnaProd,...
                bmProd(m.ntpIndexs([2 4]), iCyt), 'relative', 3e-1, 0);
            assertElementsAlmostEqual(bmComp(m.nmpIndexs([2 4]), iCyt) + decayedNMPs([2 4]), bmProd(m.ntpIndexs([2 4]), iCyt), 'relative', 3e-1, 0);
            assertElementsAlmostEqual(prodNMPs([2 4]), bmProd(m.ntpIndexs([2 4]), iCyt), 'relative', 3e-1, 0);
            
            assertElementsAlmostEqual((r.baseCounts(r.matureIndexs, m.nmpIndexs([2 4]))' * rnas0)', ...
                (bmProd(m.ntpIndexs([2 4]), iCyt) - byProd(m.nmpIndexs([2 4]), iCyt))', 'relative', 6e-1, 0);
            
            %correct biomass production, byproducts -- Protein
            monExp = (r.matureRNAGeneComposition(g.mRNAIndexs, :) * r.expression(r.matureIndexs)) ./ ...
                (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
            monExp = monExp / sum(monExp);
            mons0 = mass.dryWeightFractionProtein * mass.cellInitialDryWeight * monExp / ...
                (monExp' * pm.molecularWeights(pm.matureIndexs) / ConstantUtil.nAvogadro);
            monProd = mons0 + mons0 .* pm.decayRates(pm.matureIndexs) * t.cellCycleLength / log(2);
            
            subunits = zeros(size(g.wholeCellModelIDs));
            subunits(g.mRNAIndexs) = mons0;
            subunits(g.rRNAIndexs) = rnas0(r.matureRRNAIndexs);
            subunits(g.sRNAIndexs) = rnas0(r.matureSRNAIndexs);
            subunits(g.tRNAIndexs) = rnas0(r.matureTRNAIndexs);
            cpxs0 = min(repmat(subunits, 1, numel(pc.matureIndexs)) ./ sum(pc.proteinComplexComposition, 3), [], 1)';
            cpxs0(max((sum(pc.proteinComplexComposition, 3)' * sum(pc.proteinComplexComposition, 3)) .* (1-eye(numel(cpxs0))), [], 2) >= sum(sum(pc.proteinComplexComposition, 3),1)') = 0;
            cpxs0(pc.formationProcesses(pc.matureIndexs) == sim.processIndex('Metabolism')) = 0;
            
            subunits = zeros(size(g.wholeCellModelIDs));
            subunits(g.mRNAIndexs) = monProd;
            subunits(g.rRNAIndexs) = rnaProd(r.matureRRNAIndexs);
            subunits(g.sRNAIndexs) = rnaProd(r.matureSRNAIndexs);
            subunits(g.tRNAIndexs) = rnaProd(r.matureTRNAIndexs);
            cpxProd = min(repmat(subunits, 1, numel(pc.matureIndexs)) ./ sum(pc.proteinComplexComposition, 3), [], 1)';
            cpxProd(max((sum(pc.proteinComplexComposition, 3)' * sum(pc.proteinComplexComposition, 3)) .* (1-eye(numel(cpxProd))), [], 2) >= sum(sum(pc.proteinComplexComposition, 3),1)') = 0;
            cpxProd(pc.formationProcesses(pc.matureIndexs) == sim.processIndex('Metabolism')) = 0;
            
            aaProd = bmProd(m.aminoAcidIndexs, iCyt);
            aaProd(m.getIndexs({'GLU'}) == m.aminoAcidIndexs) = ...
                aaProd(m.getIndexs({'GLU'}) == m.aminoAcidIndexs) - ...
                aaProd(m.getIndexs({'GLN'}) == m.aminoAcidIndexs);
            aaProd(m.getIndexs({'CYS'}) == m.aminoAcidIndexs) = ...
                aaProd(m.getIndexs({'CYS'}) == m.aminoAcidIndexs) + ...
                sim.process('RNAModification').reactionStoichiometryMatrix(sim.process('RNAModification').substrateIndexs_cys, :) * ...
                sim.process('RNAModification').reactionModificationMatrix * rnaProd;
            assertElementsAlmostEqual(...
                (pm.baseCounts(pm.matureIndexs, m.aminoAcidIndexs)' * monProd)', ...
                aaProd', 'relative', 20e-2, 0);
            
            aaProd = bmProd(m.aminoAcidIndexs, iCyt);
            aaProd(m.getIndexs({'CYS'}) == m.aminoAcidIndexs) = ...
                aaProd(m.getIndexs({'CYS'}) == m.aminoAcidIndexs) + ...
                sim.process('RNAModification').reactionStoichiometryMatrix(sim.process('RNAModification').substrateIndexs_cys, :) * ...
                sim.process('RNAModification').reactionModificationMatrix * rnaProd;
            assertElementsAlmostEqual((pm.baseCounts(pm.matureIndexs, m.aminoAcidIndexs(1:20))' * mons0)', ...
                (aaProd(1:20) - byProd(m.aminoAcidIndexs(1:20), iCyt))', 'relative', 1.5e-1, 0);
            
            %correct biomass production, byproducts -- Lipids, polyamines, carbohydrates, vitamins, and ions
            lipidProd = bmComp(m.lipidIndexs, :);
            lipidProd(m.lipidIndexs == m.getIndexs({'PG160'}), sim.compartment.membraneIndexs) = ...
                lipidProd(m.lipidIndexs == m.getIndexs({'PG160'}), sim.compartment.membraneIndexs) + ...
                sum(monProd(sim.process('ProteinProcessingII').lipoproteinMonomerIndexs));
            assertElementsAlmostEqual(lipidProd, bmProd(m.lipidIndexs, :), 'relative', 2e-1, 0);
            assertAllEqual(0, byProd(m.lipidIndexs, :));
            
            assertElementsAlmostEqual(bmComp(m.polyamineIndexs, :), bmProd(m.polyamineIndexs, :), 'relative', 1e-3, 0);
            assertAllEqual(0, byProd(m.polyamineIndexs, :));
            
            assertElementsAlmostEqual(bmComp(m.carbohydrateIndexs, :), bmProd(m.carbohydrateIndexs, :), 'relative', 1e-3, 0);
            assertElementsAlmostEqual(zeros(numel(m.carbohydrateIndexs), sim.compartment.count), byProd(m.carbohydrateIndexs, :), 'absolute', 2, 0);
            
            vitBmProd = bmComp(m.vitaminIndexs, :);
            vitByProd = zeros(size(vitBmProd));
            vitBmProd(m.vitaminIndexs == m.getIndexs({'FTHF10'}), sim.compartment.cytosolIndexs) = ...
                vitBmProd(m.vitaminIndexs == m.getIndexs({'FTHF10'}), sim.compartment.cytosolIndexs) + ...
                sum(monProd);
            vitByProd(m.vitaminIndexs == m.getIndexs({'THF'}), sim.compartment.cytosolIndexs) = ...
                vitByProd(m.vitaminIndexs == m.getIndexs({'THF'}), sim.compartment.cytosolIndexs) + ...
                sum(monProd);
            vitBmProd(m.vitaminIndexs == m.getIndexs({'FTHF5'}), sim.compartment.cytosolIndexs) = ...
                vitBmProd(m.vitaminIndexs == m.getIndexs({'FTHF5'}), sim.compartment.cytosolIndexs) - ...
                sim.process('RNAModification').reactionStoichiometryMatrix(strcmp(sim.process('RNAModification').substrateWholeCellModelIDs, 'FTHF5'), :) * ...
                sim.process('RNAModification').reactionModificationMatrix * rnaProd;
            vitByProd(m.vitaminIndexs == m.getIndexs({'THF'}), sim.compartment.cytosolIndexs) = ...
                vitByProd(m.vitaminIndexs == m.getIndexs({'THF'}), sim.compartment.cytosolIndexs) + ...
                sim.process('RNAModification').reactionStoichiometryMatrix(strcmp(sim.process('RNAModification').substrateWholeCellModelIDs, 'THF'), :) * ...
                sim.process('RNAModification').reactionModificationMatrix * rnaProd;
            vitBmProd(m.vitaminIndexs == m.getIndexs({'AMET'}), sim.compartment.cytosolIndexs) = ...
                vitBmProd(m.vitaminIndexs == m.getIndexs({'AMET'}), sim.compartment.cytosolIndexs) - ...
                sim.process('RNAModification').reactionStoichiometryMatrix(strcmp(sim.process('RNAModification').substrateWholeCellModelIDs, 'AMET'), :) * ...
                sim.process('RNAModification').reactionModificationMatrix * rnaProd;                
            vitBmProd(m.vitaminIndexs == m.getIndexs({'AMET'}), sim.compartment.cytosolIndexs) = ...
                vitBmProd(m.vitaminIndexs == m.getIndexs({'AMET'}), sim.compartment.cytosolIndexs) + ...
                 2 * size(sim.process('DNARepair').RM_MunI_RecognitionSites, 1);
            assertElementsAlmostEqual(vitBmProd, bmProd(m.vitaminIndexs, :), 'relative', 5e-1, 0);
            assertElementsAlmostEqual(vitByProd, byProd(m.vitaminIndexs, :), 'relative', 2e-1, 0);
            
            ionBmProd = zeros(numel(m.ionIndexs), sim.compartment.count);
            ionByProd = zeros(numel(m.ionIndexs), sim.compartment.count);
            ionBmProd(:, sim.compartment.cytosolIndexs) = ...
                m.experimentalBiomassCompositionMolFractions(m.ionIndexs) * ...
                mass.dryWeightFractionIon / (m.experimentalBiomassCompositionMolFractions(m.ionIndexs)' * m.molecularWeights(m.ionIndexs)) * ...
                mass.cellInitialDryWeight * ConstantUtil.nAvogadro;
            [tfs, idxs] = ismember(m.wholeCellModelIDs(m.ionIndexs), sim.process('ProteinFolding').substrateWholeCellModelIDs);
            ionBmProd(tfs, sim.compartment.cytosolIndexs) = ...
                ionBmProd(tfs, sim.compartment.cytosolIndexs) + ...
                sim.process('ProteinFolding').proteinProstheticGroupMatrix(:, idxs(tfs))' * [monProd; cpxProd];
            ionByProd(tfs, sim.compartment.cytosolIndexs) = ...
                ionByProd(tfs, sim.compartment.cytosolIndexs) + ...
                sim.process('ProteinFolding').proteinProstheticGroupMatrix(:, idxs(tfs))' * ...
                [mons0 .* pm.decayRates(pm.matureIndexs) * t.cellCycleLength / log(2);
                cpxs0 .* pc.decayRates(pc.matureIndexs) * t.cellCycleLength / log(2)];
            
            tfs = ~ismember(m.wholeCellModelIDs(m.ionIndexs), {'H'; 'PI'; 'PPI'});
            assertElementsAlmostEqual(ionBmProd(tfs, :), bmProd(m.ionIndexs(tfs), :), 'relative', 10e-2, 0);
            assertElementsAlmostEqual(ionByProd(tfs, :), byProd(m.ionIndexs(tfs), :), 'relative', 20e-2, 0);
            
            %correct biomass production, byproducts -- energy
            assertIn(sim.process('Metabolism').unaccountedEnergyConsumption, ...
                [0 sim.process('Metabolism').growthAssociatedMaintenance * ConstantUtil.nAvogadro / 1000 * mass.cellInitialDryWeight])
            
            %correct biomass production and byproducts normalization
            assertElementsAlmostEqual(mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), ...
                sum((bmProd - byProd), 2)' * metMWs / ConstantUtil.nAvogadro, ...
                'relative', 5e-2, 0);
            
            assertElementsAlmostEqual(mass.dryWeightFractionDNA * mass.cellInitialDryWeight, ...
                sum(bmProd(m.dntpIndexs, :) - byProd(m.dntpIndexs, :), 2)' * ...
                (metMWs(m.dntpIndexs) - metMWs(m.diphosphateIndexs)) / ConstantUtil.nAvogadro, ...
                'relative', 10e-2, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionProtein * mass.cellInitialDryWeight, ...
                sum(bmProd(m.aminoAcidIndexs, :) - byProd(m.aminoAcidIndexs, :), 2)' * metMWs(m.aminoAcidIndexs) / ConstantUtil.nAvogadro, ...
                'relative', 2e-1, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionLipid * mass.cellInitialDryWeight, ...
                sum(bmProd(m.lipidIndexs, :) - byProd(m.lipidIndexs, :), 2)' * metMWs(m.lipidIndexs) / ConstantUtil.nAvogadro, ...
                'relative', 10e-2, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionPolyamine * mass.cellInitialDryWeight, ...
                sum(bmProd(m.polyamineIndexs, :) - byProd(m.polyamineIndexs, :), 2)' * metMWs(m.polyamineIndexs) / ConstantUtil.nAvogadro, ...
                'relative', 10e-2, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionCarbohydrate * mass.cellInitialDryWeight, ...
                sum(bmProd(m.carbohydrateIndexs, :) - byProd(m.carbohydrateIndexs, :), 2)' * metMWs(m.carbohydrateIndexs) / ConstantUtil.nAvogadro, ...
                'relative', 10e-2, 0);
            assertElementsAlmostEqual(mass.dryWeightFractionVitamin * mass.cellInitialDryWeight, ...
                sum(bmProd(m.vitaminIndexs, :) - byProd(m.vitaminIndexs, :), 2)' * metMWs(m.vitaminIndexs) / ConstantUtil.nAvogadro, ...
                'relative', 5e-1, 0);
            
            %biomass production, byproducts match biomass composition
            assertElementsAlmostEqual(bmComp(m.dnmpIndexs, :) + [[bmComp(m.getIndexs({'m6dAMP'})); 0; 0; 0] zeros(4, sim.compartment.count-1)] + bmComp(m.dntpIndexs, :),  ...
                bmProd(m.dntpIndexs, :) - byProd(m.dntpIndexs, :), ...
                'relative', 1e-2, 0);
            
            assertElementsAlmostEqual(bmComp(m.nmpIndexs, :), ...
                bmProd(m.ntpIndexs, :) - byProd(m.nmpIndexs, :)- byProd(m.ndpIndexs, :), ...
                'relative', 7e-1, 0);
            
            assertElementsAlmostEqual(bmComp(m.aminoAcidIndexs(1:20), 1)', ...
                bmProd(m.aminoAcidIndexs(1:20), 1)' - byProd(m.aminoAcidIndexs(1:20), 1)', ...
                'relative', 5e-1, 0);
            
            lipidProd = bmProd(m.lipidIndexs, :);
            lipidProd(m.lipidIndexs == m.getIndexs({'PG160'}), sim.compartment.membraneIndexs) = ...
                lipidProd(m.lipidIndexs == m.getIndexs({'PG160'}), sim.compartment.membraneIndexs) - ...
                sum(monProd(sim.process('ProteinProcessingII').lipoproteinMonomerIndexs));
            assertElementsAlmostEqual(bmComp(m.lipidIndexs, :), ...
                lipidProd - byProd(m.lipidIndexs, :), ...
                'relative', 1e-2, 1e2);
            
            assertElementsAlmostEqual(bmComp(m.polyamineIndexs, :), ...
                bmProd(m.polyamineIndexs, :) - byProd(m.polyamineIndexs, :), ...
                'relative', 1e-3, 0);
            
            p = sim.process('DNARepair');
            constants = edu.stanford.covert.util.StructUtil.catstruct(...
                sim.getFixedConstants(), sim.getFittedConstants());
            constants.processes.DNADamage.substrateGlobalIndexs = sim.process('DNADamage').substrateGlobalIndexs;
            constants.processes.DNADamage.substrateStimulusGlobalIndexs = sim.process('DNADamage').substrateStimulusGlobalIndexs;
            constants.processes.DNADamage.substrateStimulusCompartmentIndexs = sim.process('DNADamage').substrateStimulusCompartmentIndexs;
            constants.processes.DNADamage.substrateStimulusGlobalCompartmentIndexs = sim.process('DNADamage').substrateStimulusGlobalCompartmentIndexs;
            constants.processes.DNADamage.substrateStimulusLocalIndexs = sim.process('DNADamage').substrateStimulusLocalIndexs;
            
            dnaRepairRates = p.calcExpectedReactionRates(constants);
            carbByProd = byProd(m.carbohydrateIndexs, :);
            assertElementsAlmostEqual(bmComp(m.carbohydrateIndexs, :), ...
                bmProd(m.carbohydrateIndexs, :) - carbByProd, ...
                'relative', 1e-1, 2);
            
            vitBmProd = bmProd(m.vitaminIndexs, :);
            vitByProd = byProd(m.vitaminIndexs, :);
            vitBmProd(m.vitaminIndexs == m.getIndexs({'FTHF10'}), sim.compartment.cytosolIndexs) = ...
                vitBmProd(m.vitaminIndexs == m.getIndexs({'FTHF10'}), sim.compartment.cytosolIndexs) - ...
                sum(monProd);
            vitByProd(m.vitaminIndexs == m.getIndexs({'THF'}), sim.compartment.cytosolIndexs) = ...
                vitByProd(m.vitaminIndexs == m.getIndexs({'THF'}), sim.compartment.cytosolIndexs) - ...
                sum(monProd);
            vitBmProd(m.vitaminIndexs == m.getIndexs({'FTHF5'}), sim.compartment.cytosolIndexs) = ...
                vitBmProd(m.vitaminIndexs == m.getIndexs({'FTHF5'}), sim.compartment.cytosolIndexs) + ...
                sim.process('RNAModification').reactionStoichiometryMatrix(strcmp(sim.process('RNAModification').substrateWholeCellModelIDs, 'FTHF5'), :) * ...
                sim.process('RNAModification').reactionModificationMatrix * rnaProd;
            vitByProd(m.vitaminIndexs == m.getIndexs({'THF'}), sim.compartment.cytosolIndexs) = ...
                vitByProd(m.vitaminIndexs == m.getIndexs({'THF'}), sim.compartment.cytosolIndexs) - ...
                sim.process('RNAModification').reactionStoichiometryMatrix(strcmp(sim.process('RNAModification').substrateWholeCellModelIDs, 'THF'), :) * ...
                sim.process('RNAModification').reactionModificationMatrix * rnaProd;
            vitBmProd(m.vitaminIndexs == m.getIndexs({'NAD'}), sim.compartment.cytosolIndexs) = ...
                vitBmProd(m.vitaminIndexs == m.getIndexs({'NAD'}), sim.compartment.cytosolIndexs) ...
                - p.reactionSmallMoleculeStoichiometryMatrix(strcmp(p.substrateWholeCellModelIDs, 'NAD'), :) * dnaRepairRates ...
                - (2 + sum(cellfun(@numel, sim.process('Replication').primaseBindingLocations)));
            vitBmProd(m.vitaminIndexs == m.getIndexs({'AMET'}), sim.compartment.cytosolIndexs) = ...
                vitBmProd(m.vitaminIndexs == m.getIndexs({'AMET'}), sim.compartment.cytosolIndexs) + ...
                sim.process('RNAModification').reactionStoichiometryMatrix(strcmp(sim.process('RNAModification').substrateWholeCellModelIDs, 'AMET'), :) * ...
                sim.process('RNAModification').reactionModificationMatrix * rnaProd;
            vitBmProd(m.vitaminIndexs == m.getIndexs({'AMET'}), sim.compartment.cytosolIndexs) = ...
                vitBmProd(m.vitaminIndexs == m.getIndexs({'AMET'}), sim.compartment.cytosolIndexs) - ...
                2 * size(sim.process('DNARepair').RM_MunI_RecognitionSites, 1);
            assertElementsAlmostEqual(bmComp(m.vitaminIndexs, 1)', ...
                vitBmProd(:,1)' - vitByProd(:,1)', ...
                'relative', 10e-1, 1e5);
            
            assertElementsAlmostEqual(... %ntps accounted for
                this.getNucleotideCompositionMatrix(m) * bmComp, ...
                this.getNucleotideCompositionMatrix(m) * (bmProd - byProd), ...
                'relative', 1e-3);
        end
        
        function testCalcResourceRequirements_DNAWeightFraction(this)
            %import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %references
            sim = this.simulation;
            fitter = this.fitter; %#ok<*PROP>
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            
            %indices
            cIdx = sim.compartment.cytosolIndexs;
            
            %constants
            metMWs = m.molecularWeights;
            
            %set weight fractions
            mass.dryWeightFractionRNA          = 0;
            mass.dryWeightFractionDNA          = 1;
            mass.dryWeightFractionProtein      = 0;
            mass.dryWeightFractionLipid        = 0;
            mass.dryWeightFractionPolyamine    = 0;
            mass.dryWeightFractionCarbohydrate = 0;
            mass.dryWeightFractionVitamin      = 0;
            mass.dryWeightFractionIon          = 0;
            mass.dryWeightFractionNucleotide   = 0;
            
            %calculate biomass composition, maintenance energy, and byproducts
            paramVec = fitter.initializeFittedConstants();
            sim.state('MetabolicReaction').growth0 = 0;
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            [bmComp, bmProd, byProd] = fitter.calcResourceRequirements(paramVec);
            
            %assertions
            validateattributes(bmComp, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(bmProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(byProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            
            tmp_bmProd = bmProd;
            tmp_byProd = byProd;
            
            p = sim.process('DNADamage');
            rates = p.calcExpectedReactionRates();
            tmp_bmProd(p.substrateMetaboliteGlobalIndexs, cIdx) = ...
                tmp_bmProd(p.substrateMetaboliteGlobalIndexs, cIdx) - ...
                max(0, -p.reactionSmallMoleculeStoichiometryMatrix(p.substrateMetaboliteLocalIndexs, :) * rates);
            tmp_byProd(p.substrateMetaboliteGlobalIndexs, cIdx) = ...
                tmp_byProd(p.substrateMetaboliteGlobalIndexs, cIdx) - ...
                max(0, p.reactionSmallMoleculeStoichiometryMatrix(p.substrateMetaboliteLocalIndexs, :) * rates);
            
            p = sim.process('DNARepair');
            constants = edu.stanford.covert.util.StructUtil.catstruct(...
                sim.getFixedConstants(), sim.getFittedConstants());
            constants.processes.DNADamage.substrateGlobalIndexs = sim.process('DNADamage').substrateGlobalIndexs;
            constants.processes.DNADamage.substrateStimulusGlobalIndexs = sim.process('DNADamage').substrateStimulusGlobalIndexs;
            constants.processes.DNADamage.substrateStimulusCompartmentIndexs = sim.process('DNADamage').substrateStimulusCompartmentIndexs;
            constants.processes.DNADamage.substrateStimulusGlobalCompartmentIndexs = sim.process('DNADamage').substrateStimulusGlobalCompartmentIndexs;
            constants.processes.DNADamage.substrateStimulusLocalIndexs = sim.process('DNADamage').substrateStimulusLocalIndexs;
            rates = p.calcExpectedReactionRates(constants);
            tmp_bmProd(p.substrateMetaboliteGlobalIndexs, cIdx) = ...
                tmp_bmProd(p.substrateMetaboliteGlobalIndexs, cIdx) - ...
                max(0, -p.reactionSmallMoleculeStoichiometryMatrix(p.substrateMetaboliteLocalIndexs, :) * rates);
            tmp_byProd(p.substrateMetaboliteGlobalIndexs, cIdx) = ...
                tmp_byProd(p.substrateMetaboliteGlobalIndexs, cIdx) - ...
                max(0, p.reactionSmallMoleculeStoichiometryMatrix(p.substrateMetaboliteLocalIndexs, :) * rates);
            
            tmp_bmProd(m.getIndexs('ATP'), :) = tmp_bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('ATP'), :);
            tmp_bmProd(m.getIndexs('GTP'), :) = tmp_bmProd(m.getIndexs('GTP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_bmProd(m.getIndexs('H2O'), :) = tmp_bmProd(m.getIndexs('H2O'), :) - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('ADP'), :) = tmp_byProd(m.getIndexs('ADP'), :) - bmProd(m.getIndexs('ATP'), :);
            tmp_byProd(m.getIndexs('GDP'), :) = tmp_byProd(m.getIndexs('GDP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('PI'), :)  = tmp_byProd(m.getIndexs('PI'), :)  - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('H'), :)   = tmp_byProd(m.getIndexs('H'), :)   - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            
            nLigations = ...
                sum(cellfun(@(x) ~isempty(x), sim.process('Replication').primaseBindingLocations)) + ...
                sum(cellfun(@numel, sim.process('Replication').primaseBindingLocations));
            tmp_bmProd(m.getIndexs({'NAD'}), cIdx) = tmp_bmProd(m.getIndexs({'NAD'}), cIdx) - nLigations;
            tmp_byProd(m.getIndexs({'NMN'}), cIdx) = tmp_byProd(m.getIndexs({'NMN'}), cIdx) - nLigations;
            tmp_byProd(m.getIndexs({'AMP'}), cIdx) = tmp_byProd(m.getIndexs({'AMP'}), cIdx) - nLigations;
            tmp_byProd(m.getIndexs({'H'  }), cIdx) = tmp_byProd(m.getIndexs({'H'  }), cIdx) - nLigations;
            
            tmp_bmProd(m.dnmpIndexs, cIdx) = tmp_bmProd(m.dnmpIndexs, cIdx) + getBaseCounts(c.sequence);
            tmp_bmProd(m.dntpIndexs, cIdx) = tmp_bmProd(m.dntpIndexs, cIdx) - getBaseCounts(c.sequence);
            tmp_byProd(m.getIndexs('PPI'), cIdx) = tmp_byProd(m.getIndexs('PPI'), cIdx) - 2 * c.sequenceLen;
            
            tmp_bmProd(m.getIndexs('m6dAMP'), cIdx) = tmp_bmProd(m.getIndexs('m6dAMP'), cIdx) + 2*size(sim.process('DNARepair').RM_MunI_RecognitionSites, 1);
            tmp_bmProd(m.getIndexs('DAMP'),   cIdx) = tmp_bmProd(m.getIndexs('DAMP'),   cIdx) - 2*size(sim.process('DNARepair').RM_MunI_RecognitionSites, 1);
            
            tmp_bmProd(m.hydrogenIndexs, cIdx) = tmp_bmProd(m.hydrogenIndexs, cIdx) - tmp_byProd(m.hydrogenIndexs, cIdx);
            tmp_byProd(m.hydrogenIndexs, cIdx) = 0;
            
            atp = tmp_bmProd(m.getIndexs('ATP'), :);
            tmp_bmProd(m.getIndexs('ATP'), :) = tmp_bmProd(m.getIndexs('ATP'), :) - atp;
            tmp_bmProd(m.getIndexs('H2O'), :) = tmp_bmProd(m.getIndexs('H2O'), :) - atp;
            tmp_byProd(m.getIndexs('ADP'), :) = tmp_byProd(m.getIndexs('ADP'), :) - atp;
            tmp_byProd(m.getIndexs('PI'),  :) = tmp_byProd(m.getIndexs('PI'),  :) - atp;
            tmp_byProd(m.getIndexs('H'),   :) = tmp_byProd(m.getIndexs('H'),   :) - atp;
            
            h = tmp_bmProd(m.getIndexs('H'), :);
            tmp_bmProd(m.getIndexs('H'), :) = tmp_bmProd(m.getIndexs('H'), :) - h;
            tmp_byProd(m.getIndexs('H'), :) = tmp_byProd(m.getIndexs('H'), :) - h;
            
            metIdxs = [m.ntpIndexs; m.ndpIndexs; m.nmpIndexs; m.getIndexs({'H2O'; 'H'; 'PI'; 'PPI'; 'm6dAMP'; 'AHCYS'; 'NMN'; 'FTHF10'})];
            
            assertTrue(all(all(bmComp(setdiff(1:end, [m.dntpIndexs; m.dnmpIndexs; metIdxs]), :) < 100)));
            assertElementsAlmostEqual(...
                zeros(size(tmp_bmProd(setdiff(1:end, [m.dntpIndexs; m.dnmpIndexs; metIdxs]), :))), ...
                tmp_bmProd(setdiff(1:end, [m.dntpIndexs; m.dnmpIndexs; metIdxs]), :), ...
                'absolute', 100);
            assertElementsAlmostEqual(...
                zeros(size(tmp_byProd(setdiff(1:end, metIdxs), :))), ...
                tmp_byProd(setdiff(1:end, metIdxs), :), ...
                'absolute', 100);
            assertAllEqual(0, tmp_byProd(m.hydrogenIndexs, setdiff(1:end, cIdx)));
            assertElementsAlmostEqual(0, tmp_byProd(m.hydrogenIndexs, cIdx), 'absolute', 2e4);
            
            tmpMetMWs = metMWs;
            tmpMetMWs([m.dnmpIndexs; m.getIndexs({'H2O'; 'ADP'; 'GDP'; 'PI'; 'H'; 'm6dAMP'})]) = ...
                tmpMetMWs([m.dnmpIndexs; m.getIndexs({'H2O'; 'ADP'; 'GDP'; 'PI'; 'H'; 'm6dAMP'})]) - ...
                ConstantUtil.elements.O - ConstantUtil.elements.H;
            assertElementsAlmostEqual(sum(bmComp, 2)' * tmpMetMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight, 'relative', 2e-1, 0);
            assertElementsAlmostEqual(sum(bmProd - byProd, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-2, 0);
            
            assertElementsAlmostEqual(... %ntps accounted for
                (bmComp(m.ntpIndexs) + bmComp(m.ndpIndexs) + bmComp(m.nmpIndexs)), ...
                (bmProd(m.ntpIndexs) + bmProd(m.ndpIndexs) + bmProd(m.nmpIndexs)) - ...
                (byProd(m.ntpIndexs) + byProd(m.ndpIndexs) + byProd(m.nmpIndexs)) + ...
                bmProd(m.getIndexs('NAD')) * [1; 0; 0; 0], ...
                'relative', 1e-3);
        end
        
        function testCalcResourceRequirements_RNAWeightFraction(this)
            %import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %references
            sim = this.simulation;
            fitter = this.fitter; %#ok<*PROP>
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            
            %constants
            metMWs = m.molecularWeights;
            rnaLens = r.lengths(r.matureIndexs);
            
            %set weight fractions
            mass.dryWeightFractionRNA          = 1;
            mass.dryWeightFractionDNA          = 0;
            mass.dryWeightFractionProtein      = 0;
            mass.dryWeightFractionLipid        = 0;
            mass.dryWeightFractionPolyamine    = 0;
            mass.dryWeightFractionCarbohydrate = 0;
            mass.dryWeightFractionVitamin      = 0;
            mass.dryWeightFractionIon          = 0;
            mass.dryWeightFractionNucleotide   = 0;
            
            sim.process('DNADamage').reactionBounds(:) = 0;
            sim.process('Replication').primaseBindingLocations = {[] []};
            sim.process('DNARepair').RM_MunI_RecognitionSites = sim.process('DNARepair').RM_MunI_RecognitionSites([], :);
            
            %calculate biomass composition, maintenance energy, and byproducts
            paramVec = fitter.initializeFittedConstants();
            sim.state('MetabolicReaction').growth0 = 0;
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            [bmComp, bmProd, byProd] = fitter.calcResourceRequirements(paramVec);
            
            %assertions
            validateattributes(bmComp, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(bmProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(byProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            
            ts = sim.process('Transcription');
            rmod = sim.process('RNAModification');
            ta = sim.process('tRNAAminoacylation');
            
            baseCounts = zeros(numel(r.nascentIndexs), numel(m.wholeCellModelIDs));
            baseCounts(:, m.nmpIndexs) = baseCounts(:, m.nmpIndexs) + ts.transcriptionUnitBaseCounts(:, ts.substrateIndexs_nmp); %Transcription
            modMat = r.nascentRNAMatureRNAComposition' * rmod.reactionModificationMatrix' * rmod.reactionCatalysisMatrix; %RNA modification
            baseCounts(:, m.getIndexs('s2UMP')) = baseCounts(:, m.getIndexs('s2UMP')) - modMat(:, rmod.enzymeIndexs('MG_008_379_TETRAMER'));
            baseCounts(:, m.getIndexs('cmnm5s2UMP')) = baseCounts(:, m.getIndexs('cmnm5s2UMP')) + modMat(:, rmod.enzymeIndexs('MG_008_379_TETRAMER'));
            baseCounts(:, m.getIndexs('CMP')) = baseCounts(:, m.getIndexs('CMP')) - modMat(:, rmod.enzymeIndexs('MG_084_TETRAMER'));
            baseCounts(:, m.getIndexs('k2CMP')) = baseCounts(:, m.getIndexs('k2CMP')) + modMat(:, rmod.enzymeIndexs('MG_084_TETRAMER'));
            baseCounts(:, m.getIndexs('UMP')) = baseCounts(:, m.getIndexs('UMP')) - modMat(:, rmod.enzymeIndexs('MG_182_DIMER'));
            baseCounts(:, m.getIndexs('PSIURIMP')) = baseCounts(:, m.getIndexs('PSIURIMP')) + modMat(:, rmod.enzymeIndexs('MG_182_DIMER'));
            baseCounts(:, m.getIndexs('UMP')) = baseCounts(:, m.getIndexs('UMP')) - modMat(:, rmod.enzymeIndexs('MG_209_MONOMER'));
            baseCounts(:, m.getIndexs('PSIURIMP')) = baseCounts(:, m.getIndexs('PSIURIMP')) + modMat(:, rmod.enzymeIndexs('MG_209_MONOMER'));
            baseCounts(:, m.getIndexs('GMP')) = baseCounts(:, m.getIndexs('GMP')) - modMat(:, rmod.enzymeIndexs('MG_252_DIMER'));
            baseCounts(:, m.getIndexs('GmMP')) = baseCounts(:, m.getIndexs('GmMP')) + modMat(:, rmod.enzymeIndexs('MG_252_DIMER'));
            baseCounts(:, m.getIndexs('UMP')) = baseCounts(:, m.getIndexs('UMP')) - modMat(:, rmod.enzymeIndexs('MG_295_MONOMER'));
            baseCounts(:, m.getIndexs('s2UMP')) = baseCounts(:, m.getIndexs('s2UMP')) + modMat(:, rmod.enzymeIndexs('MG_295_MONOMER'));
            baseCounts(:, m.getIndexs('UMP')) = baseCounts(:, m.getIndexs('UMP')) - modMat(:, rmod.enzymeIndexs('MG_346_DIMER'));
            baseCounts(:, m.getIndexs('UmMP')) = baseCounts(:, m.getIndexs('UmMP')) + modMat(:, rmod.enzymeIndexs('MG_346_DIMER'));
            baseCounts(:, m.getIndexs('GMP')) = baseCounts(:, m.getIndexs('GMP')) - modMat(:, rmod.enzymeIndexs('MG_347_DIMER'));
            baseCounts(:, m.getIndexs('m7GMP')) = baseCounts(:, m.getIndexs('m7GMP')) + modMat(:, rmod.enzymeIndexs('MG_347_DIMER'));
            baseCounts(:, m.getIndexs('UMP')) = baseCounts(:, m.getIndexs('UMP')) - modMat(:, rmod.enzymeIndexs('MG_370_MONOMER'));
            baseCounts(:, m.getIndexs('PSIURIMP')) = baseCounts(:, m.getIndexs('PSIURIMP')) + modMat(:, rmod.enzymeIndexs('MG_370_MONOMER'));
            baseCounts(:, m.getIndexs('UMP')) = baseCounts(:, m.getIndexs('UMP')) - modMat(:, rmod.enzymeIndexs('MG_372_DIMER'));
            baseCounts(:, m.getIndexs('s4UMP')) = baseCounts(:, m.getIndexs('s4UMP')) + modMat(:, rmod.enzymeIndexs('MG_372_DIMER'));
            baseCounts(:, m.getIndexs('GMP')) = baseCounts(:, m.getIndexs('GMP')) - modMat(:, rmod.enzymeIndexs('MG_380_MONOMER'));
            baseCounts(:, m.getIndexs('m2GMP')) = baseCounts(:, m.getIndexs('m2GMP')) + modMat(:, rmod.enzymeIndexs('MG_380_MONOMER'));
            baseCounts(:, m.getIndexs('GMP')) = baseCounts(:, m.getIndexs('GMP')) - modMat(:, rmod.enzymeIndexs('MG_445_DIMER'));
            baseCounts(:, m.getIndexs('m1GMP')) = baseCounts(:, m.getIndexs('m1GMP')) + modMat(:, rmod.enzymeIndexs('MG_445_DIMER'));
            baseCounts(:, m.getIndexs('AMP')) = baseCounts(:, m.getIndexs('AMP')) - modMat(:, rmod.enzymeIndexs('MG_463_MONOMER')) / 2;
            baseCounts(:, m.getIndexs('m62AMP')) = baseCounts(:, m.getIndexs('m62AMP')) + modMat(:, rmod.enzymeIndexs('MG_463_MONOMER')) / 2;
            baseCounts = baseCounts - r.intergenicRNAMatrix' * r.baseCounts(r.intergenicIndexs, :); %RNA processing
            assertEqual(r.nascentRNAMatureRNAComposition' * r.baseCounts(r.matureIndexs, :), baseCounts);
            baseCounts(:, m.aminoAcidIndexs) = baseCounts(:, m.aminoAcidIndexs) - ... %tRNA aminoacylation
                r.nascentRNAMatureRNAComposition(ta.aminoacylatedRNAGlobalIndexs, :)' * ...
                ta.reactionModificationMatrix' * ta.reactionStoichiometryMatrix([ta.substrateIndexs_aminoAcids; ta.substrateIndexs_fmethionine], :)';
            idx = strcmp(r.wholeCellModelIDs(r.nascentIndexs), 'TU_139');
            baseCounts(idx, m.getIndexs('MET'))  = baseCounts(idx, m.getIndexs('MET'))  - 1;
            baseCounts(idx, m.getIndexs('FMET')) = baseCounts(idx, m.getIndexs('FMET')) + 1;
            assertEqual(r.nascentRNAMatureRNAComposition' * r.baseCounts(r.aminoacylatedIndexs, :), baseCounts);
            
            tmpMetMWs = metMWs;
            tmpMetMWs([m.nmpIndexs; m.getIndexs({'s2UMP'; 'cmnm5s2UMP'; 'k2CMP'; 'PSIURIMP'; 'GmMP'; 's2UMP'; 'UmMP'; 'm7GMP'; 's4UMP'; 'm2GMP'; 'm1GMP'; 'm62AMP'})]) = ...
                tmpMetMWs([m.nmpIndexs; m.getIndexs({'s2UMP'; 'cmnm5s2UMP'; 'k2CMP'; 'PSIURIMP'; 'GmMP'; 's2UMP'; 'UmMP'; 'm7GMP'; 's4UMP'; 'm2GMP'; 'm1GMP'; 'm62AMP'})]) - ...
                (ConstantUtil.elements.O + ConstantUtil.elements.H) * (1 - 1/(r.expression(r.matureIndexs)' * rnaLens));
            assertElementsAlmostEqual(sum(bmComp, 2)' * tmpMetMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-3, 0);
            assertElementsAlmostEqual(sum(bmProd - byProd, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 2e-2, 0);
            
            assertElementsAlmostEqual(... %ntps accounted for
                this.getNucleotideCompositionMatrix(m) * bmComp, ...
                this.getNucleotideCompositionMatrix(m) * (bmProd - byProd), ...
                'relative', 1e-3);
        end
        
        function testCalcResourceRequirements_ProteinWeightFraction(this)
            %import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %references
            sim = this.simulation;
            fitter = this.fitter; %#ok<*PROP>
            g = sim.gene;
            t = sim.state('Time');
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            %constants
            metMWs = m.molecularWeights;
            rnaExp = r.expression(r.matureIndexs);
            monExp = (r.matureRNAGeneComposition(g.mRNAIndexs, :) * rnaExp) ./ ...
                (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
            monExp = monExp / sum(monExp);
            monLens = pm.lengths(pm.matureIndexs);
            
            %set weight fractions
            mass.dryWeightFractionRNA          = 0;
            mass.dryWeightFractionDNA          = 0;
            mass.dryWeightFractionProtein      = 1;
            mass.dryWeightFractionLipid        = 0;
            mass.dryWeightFractionPolyamine    = 0;
            mass.dryWeightFractionCarbohydrate = 0;
            mass.dryWeightFractionVitamin      = 0;
            mass.dryWeightFractionIon          = 0;
            mass.dryWeightFractionNucleotide   = 0;
            
            sim.process('DNADamage').reactionBounds(:) = 0;
            sim.process('Replication').primaseBindingLocations = {[] []};
            sim.process('DNARepair').RM_MunI_RecognitionSites = sim.process('DNARepair').RM_MunI_RecognitionSites([], :);
            
            %calculate biomass composition, maintenance energy, and byproducts
            paramVec = fitter.initializeFittedConstants();
            sim.state('MetabolicReaction').growth0 = 1 / sim.state('Mass').timeAveragedCellWeight / sim.state('Time').cellCycleLength;
            sim.state('ProteinMonomer').minimumAverageExpression = 0;
            [bmComp, bmProd, byProd] = fitter.calcResourceRequirements(paramVec);
            
            %assertions
            validateattributes(bmComp, {'numeric'}, {'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(bmComp(~any(pc.baseCounts < 0, 1), :), {'numeric'}, {'nonnegative', 'real'});
            validateattributes(bmProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(byProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            
            ta = sim.process('tRNAAminoacylation');
            pd = sim.process('ProteinDecay');
            pmod = sim.process('ProteinModification');
            pf = sim.process('ProteinFolding');
            pp1 = sim.process('ProteinProcessingI');
            pp2 = sim.process('ProteinProcessingII');
            
            baseCounts = zeros(numel(pm.matureIndexs), numel(m.wholeCellModelIDs));
            [~, idxs] = ismember(ta.aminoacylatedRNAtRNAIndexs, ta.aminoacylatedRNAGlobalIndexs); %tRNA aminoacylation
            baseCounts(:, m.aminoAcidIndexs) = baseCounts(:, m.aminoAcidIndexs) - (...
                ta.reactionStoichiometryMatrix([ta.substrateIndexs_aminoAcids; ta.substrateIndexs_fmethionine], :) * ...
                ta.reactionModificationMatrix(:, idxs) * ta.monomerTRNACounts'...
                )';
            baseCounts(:, m.getIndexs({'MET'})) = baseCounts(:, m.getIndexs({'MET'})) - pp1.nascentMonomerNTerminalMethionineCleavages; %protein processing I
            baseCounts(pp2.lipoproteinMonomerIndexs, pp2.substrateGlobalIndexs(pp2.substrateIndexs_diacylglycerolCys)) = ... %Protein processing II
                baseCounts(pp2.lipoproteinMonomerIndexs, pp2.substrateGlobalIndexs(pp2.substrateIndexs_diacylglycerolCys)) + 1;
            baseCounts(pp2.lipoproteinMonomerIndexs, m.getIndexs({'CYS'})) = ...
                baseCounts(pp2.lipoproteinMonomerIndexs, m.getIndexs({'CYS'})) - 1;
            baseCounts(:, m.aminoAcidIndexs) = ... %Protein processing II / decay
                baseCounts(:, m.aminoAcidIndexs) - ...
                pd.monomerDecayReactions(pd.substrateIndexs_aminoAcids, pd.monomer.signalSequenceIndexs)';
            baseCounts(:, pf.substrateMetaboliteGlobalIndexs) = ... %Protein folding
                baseCounts(:, pf.substrateMetaboliteGlobalIndexs) + ...
                pf.proteinProstheticGroupMatrix(1:numel(pm.matureIndexs), pf.substrateMetaboliteLocalIndexs);
            baseCounts(:, m.getIndexs('GLU')) = ... %Protein modification
                baseCounts(:, m.getIndexs('GLU')) - ...
                (pmod.reactionStoichiometryMatrix(pmod.substrateIndexs('GLU'), :) * pmod.reactionModificationMatrix)';
            baseCounts(:, m.getIndexs({'SER';'THR';'TYR';'LYS'})) = ...
                baseCounts(:, m.getIndexs({'SER';'THR';'TYR';'LYS'})) ...
                - pm.baseCounts(pm.matureIndexs, m.getIndexs({'pSER';'pTHR';'pTYR';'LIPOYLLYS'}));
            baseCounts(:, m.getIndexs({'pSER';'pTHR';'pTYR';'LIPOYLLYS'})) = ...
                baseCounts(:, m.getIndexs({'pSER';'pTHR';'pTYR';'LIPOYLLYS'})) + ...
                pm.baseCounts(pm.matureIndexs, m.getIndexs({'pSER';'pTHR';'pTYR';'LIPOYLLYS'}));
            assertEqual(pm.baseCounts(pm.matureIndexs, :), baseCounts);
            
            tmpMetMWs = metMWs;
            tmpMetMWs(m.aminoAcidIndexs) = ...
                tmpMetMWs(m.aminoAcidIndexs) - ...
                (ConstantUtil.elements.O + 2 * ConstantUtil.elements.H) * (1 - 1/(monExp' * monLens));
            assertElementsAlmostEqual(sum(bmComp, 2)' * tmpMetMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 5e-2, 0);
            assertElementsAlmostEqual(sum(bmProd - byProd, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-2, 0);
            
            assertElementsAlmostEqual(... %ntps accounted for
                this.getNucleotideCompositionMatrix(m) * bmComp, ...
                this.getNucleotideCompositionMatrix(m) * (bmProd - byProd), ...
                'relative', 1e-3);
        end
        
        function testCalcResourceRequirements_LipidWeightFraction(this)
            %import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %references
            sim = this.simulation;
            fitter = this.fitter; %#ok<*PROP>
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            
            %constants
            metMWs = m.molecularWeights;
            
            %set weight fractions
            mass.dryWeightFractionRNA          = 0;
            mass.dryWeightFractionDNA          = 0;
            mass.dryWeightFractionProtein      = 0;
            mass.dryWeightFractionLipid        = 1;
            mass.dryWeightFractionPolyamine    = 0;
            mass.dryWeightFractionCarbohydrate = 0;
            mass.dryWeightFractionVitamin      = 0;
            mass.dryWeightFractionIon          = 0;
            mass.dryWeightFractionNucleotide   = 0;
            
            sim.process('DNADamage').reactionBounds(:) = 0;
            sim.process('Replication').primaseBindingLocations = {[] []};
            sim.process('DNARepair').RM_MunI_RecognitionSites = sim.process('DNARepair').RM_MunI_RecognitionSites([], :);
            
            %calculate biomass composition, maintenance energy, and byproducts
            paramVec = fitter.initializeFittedConstants();
            sim.state('MetabolicReaction').growth0 = 0;
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            [bmComp, bmProd, byProd] = fitter.calcResourceRequirements(paramVec);
            
            %assertions
            validateattributes(bmComp, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(bmProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(byProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            
            tmp_bmProd = bmProd;
            tmp_byProd = byProd;
            tmp_bmProd(m.getIndexs('ATP'), :) = tmp_bmProd(m.getIndexs('ATP'), :) - byProd(m.getIndexs('ADP'), :);
            tmp_bmProd(m.getIndexs('GTP'), :) = tmp_bmProd(m.getIndexs('GTP'), :) - byProd(m.getIndexs('GDP'), :);
            tmp_bmProd(m.getIndexs('H2O'), :) = tmp_bmProd(m.getIndexs('H2O'), :) - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('ADP'), :) = tmp_byProd(m.getIndexs('ADP'), :) - bmProd(m.getIndexs('ATP'), :);
            tmp_byProd(m.getIndexs('GDP'), :) = tmp_byProd(m.getIndexs('GDP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('PI'), :)  = tmp_byProd(m.getIndexs('PI'), :)  - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('H'), :)   = tmp_byProd(m.getIndexs('H'), :)   - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            
            metIdxs = [m.ntpIndexs; m.ndpIndexs; m.nmpIndexs; m.getIndexs({'H2O'; 'H'; 'PI'; 'PPI'; 'FTHF10'})];
            
            assertAllEqual(0, bmComp(setdiff(1:end, [m.lipidIndexs; metIdxs]), :));
            assertElementsAlmostEqual(zeros(size(tmp_bmProd(setdiff(1:end, [m.lipidIndexs; metIdxs]), :))), ...
                tmp_bmProd(setdiff(1:end, [m.lipidIndexs; metIdxs]), :), 'absolute', 1e-4);
            assertElementsAlmostEqual(zeros(size(tmp_byProd(setdiff(1:end, [m.lipidIndexs; metIdxs]), :))), ...
                tmp_byProd(setdiff(1:end, [m.lipidIndexs; metIdxs]), :), 'absolute', 1e-4);
            
            assertElementsAlmostEqual(sum(bmComp, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-8, 0);
            assertElementsAlmostEqual(sum(bmProd - byProd, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-8, 0);
            
            assertElementsAlmostEqual(bmComp(setdiff(1:end, metIdxs), :), ...
                + tmp_bmProd(setdiff(1:end, metIdxs), :) ...
                - tmp_byProd(setdiff(1:end, metIdxs), :), ...
                'relative', 1e-3);
            
            assertElementsAlmostEqual(... %ntps accounted for
                this.getNucleotideCompositionMatrix(m) * bmComp, ...
                this.getNucleotideCompositionMatrix(m) * (bmProd - byProd), ...
                'relative', 1e-3);
        end
        
        function testCalcResourceRequirements_PolyamineWeightFraction(this)
            %import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %references
            sim = this.simulation;
            fitter = this.fitter; %#ok<*PROP>
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            
            %constants
            metMWs = m.molecularWeights;
            
            %set weight fractions
            mass.dryWeightFractionRNA          = 0;
            mass.dryWeightFractionDNA          = 0;
            mass.dryWeightFractionProtein      = 0;
            mass.dryWeightFractionLipid        = 0;
            mass.dryWeightFractionPolyamine    = 1;
            mass.dryWeightFractionCarbohydrate = 0;
            mass.dryWeightFractionVitamin      = 0;
            mass.dryWeightFractionIon          = 0;
            mass.dryWeightFractionNucleotide   = 0;
            
            sim.process('Replication').primaseBindingLocations = {[] []};
            sim.process('DNADamage').reactionBounds(:) = 0;
            sim.process('DNARepair').RM_MunI_RecognitionSites = sim.process('DNARepair').RM_MunI_RecognitionSites([], :);
            
            %calculate biomass composition, maintenance energy, and byproducts
            paramVec = fitter.initializeFittedConstants();
            sim.state('MetabolicReaction').growth0 = 0;
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            [bmComp, bmProd, byProd] = fitter.calcResourceRequirements(paramVec);
            
            %assertions
            validateattributes(bmComp, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(bmProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(byProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            
            tmp_bmProd = bmProd;
            tmp_byProd = byProd;
            tmp_bmProd(m.getIndexs('ATP'), :) = tmp_bmProd(m.getIndexs('ATP'), :) - byProd(m.getIndexs('ADP'), :);
            tmp_bmProd(m.getIndexs('GTP'), :) = tmp_bmProd(m.getIndexs('GTP'), :) - byProd(m.getIndexs('GDP'), :);
            tmp_bmProd(m.getIndexs('H2O'), :) = tmp_bmProd(m.getIndexs('H2O'), :) - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('ADP'), :) = tmp_byProd(m.getIndexs('ADP'), :) - bmProd(m.getIndexs('ATP'), :);
            tmp_byProd(m.getIndexs('GDP'), :) = tmp_byProd(m.getIndexs('GDP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('PI'), :)  = tmp_byProd(m.getIndexs('PI'), :)  - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('H'), :)   = tmp_byProd(m.getIndexs('H'), :)   - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            
            metIdxs = [m.ntpIndexs; m.ndpIndexs; m.nmpIndexs; m.getIndexs({'H2O'; 'H'; 'PI'; 'PPI'; 'FTHF10'})];
            
            assertAllEqual(0, bmComp(setdiff(1:end, [m.polyamineIndexs; metIdxs]), :));
            assertElementsAlmostEqual(zeros(size(tmp_bmProd(setdiff(1:end, [m.polyamineIndexs; metIdxs]), :))), ...
                tmp_bmProd(setdiff(1:end, [m.polyamineIndexs; metIdxs]), :), 'absolute', 1e-4);
            assertElementsAlmostEqual(zeros(size(tmp_byProd(setdiff(1:end, [m.polyamineIndexs; metIdxs]), :))), ...
                tmp_byProd(setdiff(1:end, [m.polyamineIndexs; metIdxs]), :), 'absolute', 1e-4);
            
            assertElementsAlmostEqual(sum(bmComp, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-8, 0);
            assertElementsAlmostEqual(sum(bmProd - byProd, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-8, 0);
            
            assertElementsAlmostEqual(bmComp(setdiff(1:end, metIdxs), :), ...
                tmp_bmProd(setdiff(1:end, metIdxs), :) - ...
                tmp_byProd(setdiff(1:end, metIdxs), :), ...
                'relative', 1e-3);
            
            assertElementsAlmostEqual(... %ntps accounted for
                this.getNucleotideCompositionMatrix(m) * bmComp, ...
                this.getNucleotideCompositionMatrix(m) * (bmProd - byProd), ...
                'relative', 1e-3);
        end
        
        function testCalcResourceRequirements_CarbohydrateWeightFraction(this)
            %import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %references
            sim = this.simulation;
            fitter = this.fitter; %#ok<*PROP>
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
                        
            %constants
            metMWs = m.molecularWeights;
            
            %set weight fractions
            mass.dryWeightFractionRNA          = 0;
            mass.dryWeightFractionDNA          = 0;
            mass.dryWeightFractionProtein      = 0;
            mass.dryWeightFractionLipid        = 0;
            mass.dryWeightFractionPolyamine    = 0;
            mass.dryWeightFractionCarbohydrate = 1;
            mass.dryWeightFractionVitamin      = 0;
            mass.dryWeightFractionIon          = 0;
            mass.dryWeightFractionNucleotide   = 0;
            
            sim.process('Replication').primaseBindingLocations = {[] []};
            sim.process('DNADamage').reactionBounds(:) = 0;
            sim.process('DNARepair').RM_MunI_RecognitionSites = sim.process('DNARepair').RM_MunI_RecognitionSites([], :);
            
            %calculate biomass composition, maintenance energy, and byproducts
            paramVec = fitter.initializeFittedConstants();
            sim.state('MetabolicReaction').growth0 = 0;
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            [bmComp, bmProd, byProd] = fitter.calcResourceRequirements(paramVec);
            
            %assertions
            validateattributes(bmComp, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(bmProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(byProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            
            tmp_bmProd = bmProd;
            tmp_byProd = byProd;
            tmp_bmProd(m.getIndexs('ATP'), :) = tmp_bmProd(m.getIndexs('ATP'), :) - byProd(m.getIndexs('ADP'), :);
            tmp_bmProd(m.getIndexs('GTP'), :) = tmp_bmProd(m.getIndexs('GTP'), :) - byProd(m.getIndexs('GDP'), :);
            tmp_bmProd(m.getIndexs('H2O'), :) = tmp_bmProd(m.getIndexs('H2O'), :) - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('ADP'), :) = tmp_byProd(m.getIndexs('ADP'), :) - bmProd(m.getIndexs('ATP'), :);
            tmp_byProd(m.getIndexs('GDP'), :) = tmp_byProd(m.getIndexs('GDP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('PI'), :)  = tmp_byProd(m.getIndexs('PI'), :)  - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('H'), :)   = tmp_byProd(m.getIndexs('H'), :)   - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            
            metIdxs = [m.ntpIndexs; m.ndpIndexs; m.nmpIndexs; m.getIndexs({'H2O'; 'H'; 'PI'; 'PPI'; 'FTHF10'})];
            
            assertAllEqual(0, bmComp(setdiff(1:end, [m.carbohydrateIndexs; metIdxs]), :));
            assertElementsAlmostEqual(zeros(size(tmp_bmProd(setdiff(1:end, [m.carbohydrateIndexs; metIdxs]), :))), ...
                tmp_bmProd(setdiff(1:end, [m.carbohydrateIndexs; metIdxs]), :), 'absolute', 1e-4);
            assertElementsAlmostEqual(zeros(size(tmp_byProd(setdiff(1:end, [m.carbohydrateIndexs; metIdxs]), :))), ...
                tmp_byProd(setdiff(1:end, [m.carbohydrateIndexs; metIdxs]), :), 'absolute', 1e-4);
            
            assertElementsAlmostEqual(sum(bmComp, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-8, 0);
            assertElementsAlmostEqual(sum(bmProd - byProd, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-8, 0);
            
            assertElementsAlmostEqual(bmComp(setdiff(1:end, metIdxs), :), ...
                tmp_bmProd(setdiff(1:end, metIdxs), :) - ...
                tmp_byProd(setdiff(1:end, metIdxs), :), ...
                'relative', 1e-3);
            
            assertElementsAlmostEqual(... %ntps accounted for
                this.getNucleotideCompositionMatrix(m) * bmComp, ...
                this.getNucleotideCompositionMatrix(m) * (bmProd - byProd), ...
                'relative', 1e-3);
        end
        
        function testCalcResourceRequirements_VitaminWeightFraction(this)
            %import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %references
            sim = this.simulation;
            fitter = this.fitter; %#ok<*PROP>
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            
            %constants
            metMWs = m.molecularWeights;
            
            %set weight fractions
            mass.dryWeightFractionRNA          = 0;
            mass.dryWeightFractionDNA          = 0;
            mass.dryWeightFractionProtein      = 0;
            mass.dryWeightFractionLipid        = 0;
            mass.dryWeightFractionPolyamine    = 0;
            mass.dryWeightFractionCarbohydrate = 0;
            mass.dryWeightFractionVitamin      = 1;
            mass.dryWeightFractionIon          = 0;
            mass.dryWeightFractionNucleotide   = 0;
            
            sim.process('DNADamage').reactionBounds(:) = 0;
            sim.process('Replication').primaseBindingLocations = {[] []};
            sim.process('DNARepair').RM_MunI_RecognitionSites = sim.process('DNARepair').RM_MunI_RecognitionSites([], :);
            
            %calculate biomass composition, maintenance energy, and byproducts
            paramVec = fitter.initializeFittedConstants();
            sim.state('MetabolicReaction').growth0 = 0;
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            [bmComp, bmProd, byProd] = fitter.calcResourceRequirements(paramVec);
            
            %assertions
            validateattributes(bmComp, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(bmProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(byProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            
            tmp_bmProd = bmProd;
            tmp_byProd = byProd;
            tmp_bmProd(m.getIndexs('ATP'), :) = tmp_bmProd(m.getIndexs('ATP'), :) - byProd(m.getIndexs('ADP'), :);
            tmp_bmProd(m.getIndexs('GTP'), :) = tmp_bmProd(m.getIndexs('GTP'), :) - byProd(m.getIndexs('GDP'), :);
            tmp_bmProd(m.getIndexs('H2O'), :) = tmp_bmProd(m.getIndexs('H2O'), :) - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('ADP'), :) = tmp_byProd(m.getIndexs('ADP'), :) - bmProd(m.getIndexs('ATP'), :);
            tmp_byProd(m.getIndexs('GDP'), :) = tmp_byProd(m.getIndexs('GDP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('PI'), :)  = tmp_byProd(m.getIndexs('PI'), :)  - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('H'), :)   = tmp_byProd(m.getIndexs('H'), :)   - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            
            metIdxs = [m.ntpIndexs; m.ndpIndexs; m.nmpIndexs; m.getIndexs({'H2O'; 'H'; 'PI'; 'PPI'; 'FTHF10'})];
            
            assertAllEqual(0, bmComp(setdiff(1:end, [m.vitaminIndexs; metIdxs]), :));
            assertElementsAlmostEqual(zeros(size(tmp_bmProd(setdiff(1:end, [m.vitaminIndexs; metIdxs]), :))), ...
                tmp_bmProd(setdiff(1:end, [m.vitaminIndexs; metIdxs]), :), 'absolute', 1e-4);
            assertElementsAlmostEqual(zeros(size(tmp_byProd(setdiff(1:end, [m.vitaminIndexs; metIdxs]), :))), ...
                tmp_byProd(setdiff(1:end, [m.vitaminIndexs; metIdxs]), :), 'absolute', 1e-4);
            
            assertElementsAlmostEqual(sum(bmComp, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-8, 0);
            assertElementsAlmostEqual(sum(bmProd - byProd, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-8, 0);
            
            assertElementsAlmostEqual(bmComp(setdiff(1:end, metIdxs), :), ...
                tmp_bmProd(setdiff(1:end, metIdxs), :) - ...
                tmp_byProd(setdiff(1:end, metIdxs), :), ...
                'relative', 1e-3);
            
            assertElementsAlmostEqual(... %ntps accounted for
                this.getNucleotideCompositionMatrix(m) * bmComp, ...
                this.getNucleotideCompositionMatrix(m) * (bmProd - byProd), ...
                'relative', 1e-3);
        end
        
        function testCalcResourceRequirements_IonWeightFraction(this)
            %import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %references
            sim = this.simulation;
            fitter = this.fitter; %#ok<*PROP>
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            
            %constants
            metMWs = m.molecularWeights;
            
            %set weight fractions
            mass.dryWeightFractionRNA          = 0;
            mass.dryWeightFractionDNA          = 0;
            mass.dryWeightFractionProtein      = 0;
            mass.dryWeightFractionLipid        = 0;
            mass.dryWeightFractionPolyamine    = 0;
            mass.dryWeightFractionCarbohydrate = 0;
            mass.dryWeightFractionVitamin      = 0;
            mass.dryWeightFractionIon          = 1;
            mass.dryWeightFractionNucleotide   = 0;
            
            sim.process('DNADamage').reactionBounds(:) = 0;
            sim.process('Replication').primaseBindingLocations = {[] []};
            sim.process('DNARepair').RM_MunI_RecognitionSites = sim.process('DNARepair').RM_MunI_RecognitionSites([], :);
            
            %calculate biomass composition, maintenance energy, and byproducts
            paramVec = fitter.initializeFittedConstants();
            sim.state('MetabolicReaction').growth0 = 0;
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            [bmComp, bmProd, byProd] = fitter.calcResourceRequirements(paramVec);
            
            %assertions
            validateattributes(bmComp, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(bmProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            validateattributes(byProd, {'numeric'}, {'nonnegative', 'real', 'size', [numel(m.wholeCellModelIDs) sim.compartment.count]});
            
            tmp_bmProd = bmProd;
            tmp_byProd = byProd;
            tmp_bmProd(m.getIndexs('ATP'), :) = tmp_bmProd(m.getIndexs('ATP'), :) - byProd(m.getIndexs('ADP'), :);
            tmp_bmProd(m.getIndexs('GTP'), :) = tmp_bmProd(m.getIndexs('GTP'), :) - byProd(m.getIndexs('GDP'), :);
            tmp_bmProd(m.getIndexs('H2O'), :) = tmp_bmProd(m.getIndexs('H2O'), :) - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('ADP'), :) = tmp_byProd(m.getIndexs('ADP'), :) - bmProd(m.getIndexs('ATP'), :);
            tmp_byProd(m.getIndexs('GDP'), :) = tmp_byProd(m.getIndexs('GDP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('PI'), :)  = tmp_byProd(m.getIndexs('PI'), :)  - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            tmp_byProd(m.getIndexs('H'), :)   = tmp_byProd(m.getIndexs('H'), :)   - bmProd(m.getIndexs('ATP'), :) - bmProd(m.getIndexs('GTP'), :);
            
            metIdxs = [m.ntpIndexs; m.ndpIndexs; m.nmpIndexs; m.getIndexs({'H2O'; 'H'; 'PI'; 'PPI'; 'FTHF10'})];
            
            assertAllEqual(0, bmComp(setdiff(1:end, [m.ionIndexs; metIdxs]), :));
            assertElementsAlmostEqual(zeros(size(tmp_bmProd(setdiff(1:end, [m.ionIndexs; metIdxs]), :))), ...
                tmp_bmProd(setdiff(1:end, [m.ionIndexs; metIdxs]), :), 'absolute', 1e-3);
            assertElementsAlmostEqual(zeros(size(tmp_byProd(setdiff(1:end, [m.ionIndexs; metIdxs]), :))), ...
                tmp_byProd(setdiff(1:end, [m.ionIndexs; metIdxs]), :), 'absolute', 1e-10);
            
            assertElementsAlmostEqual(sum(bmComp, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-8, 0);
            assertElementsAlmostEqual(sum(bmProd - byProd, 2)' * metMWs / ConstantUtil.nAvogadro, ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight), 'relative', 1e-8, 0);
            
            assertElementsAlmostEqual(bmComp(setdiff(1:end, metIdxs), :), ...
                tmp_bmProd(setdiff(1:end, metIdxs), :) - ...
                tmp_byProd(setdiff(1:end, metIdxs), :), ...
                'relative', 1e-3);
            
            assertElementsAlmostEqual(... %ntps accounted for
                this.getNucleotideCompositionMatrix(m) * bmComp, ...
                this.getNucleotideCompositionMatrix(m) * (bmProd - byProd), ...
                'relative', 1e-3);
        end
    end
    
    %helper methods
    methods (Static = true)
        function A = getNucleotideCompositionMatrix(m)
            A = zeros(4, numel(m.wholeCellModelIDs));
            A(sub2ind(size(A), 1:4, m.ntpIndexs')) = 1;
            A(sub2ind(size(A), 1:4, m.ndpIndexs')) = 1;
            A(sub2ind(size(A), 1:4, m.nmpIndexs')) = 1;
            A(4, m.getIndexs('cmnm5s2UMP')) = 1;
            A(3, m.getIndexs('GmMP')) = 1;
            A(2, m.getIndexs('k2CMP')) = 1;
            A(3, m.getIndexs('m1GMP')) = 1;
            A(3, m.getIndexs('m2GMP')) = 1;
            A(1, m.getIndexs('m62AMP')) = 1;
            A(3, m.getIndexs('m7GMP')) = 1;
            A(4, m.getIndexs('PSIURIMP')) = 1;
            A(4, m.getIndexs('s2UMP')) = 1;
            A(4, m.getIndexs('s4UMP')) = 1;
            A(4, m.getIndexs('UmMP')) = 1;
            A(1, m.getIndexs('NAD')) = 1;
            A(1, m.getIndexs('LIPOYLAMP')) = 1;
        end
    end
end