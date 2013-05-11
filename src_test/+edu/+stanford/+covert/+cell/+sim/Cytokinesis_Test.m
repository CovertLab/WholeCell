%Cytokinesis medium test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/29/2010
classdef Cytokinesis_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = Cytokinesis_Test(name)
            this = this@TestCase(name);
        end
        
        %used to estimate the duration of cytokinesis
        function testCytokinesisDuration(this)
            %% simulate
            durations = this.estimateMeanCytokinesisDurationDistribution(5, 0);
            
            %% assert
            %total time for cytokinesis near cytokinesis duration
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {});
            time = sim.state('Time');
            assertElementsAlmostEqual(time.cytokinesisDuration, mean(durations), 'relative', 0.25);
        end
        
        function testSufficientFtsZ(~)
            import edu.stanford.covert.cell.sim.ProteinGrowth_Test;
            import edu.stanford.covert.util.ConstantUtil;
            
            %% simulate
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'FtsZPolymerization'});
            g = sim.gene;
            time = sim.state('Time');
            mass = sim.state('Mass');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            ring = sim.state('FtsZRing');
            ftsZPol = sim.process('FtsZPolymerization');
            
            ftsZMonIdxs = ftsZPol.enzymeGlobalIndexs(ftsZPol.enzymeIndexs_FtsZ);
            ftsZCpxIdxs = ftsZPol.enzymeGlobalIndexs([ftsZPol.enzymeIndexs_FtsZ_GDP; ftsZPol.enzymeIndexs_FtsZ_activated]);
            maxTime = time.cellCycleLength - time.cytokinesisDuration;
            minFtsZ = 2.25 * ring.numEdges * ring.numFtsZSubunitsPerFilament;
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            monExp = mRNAExp ./ (log(2) / time.cellCycleLength + pm.decayRates(pm.matureIndexs));
            monExp = monExp / sum(monExp);
            monCnts = mass.cellInitialDryWeight * mass.dryWeightFractionProtein * ConstantUtil.nAvogadro / ...
                (monExp' * pm.molecularWeights(pm.matureIndexs)) * monExp;
            %meanFtsZCopyNumber = 3.45 * ring.numEdges * ring.numFtsZSubunitsPerFilament;
            meanFtsZCopyNumber = monCnts(ftsZPol.enzymeMonomerGlobalIndexs);
            
            tfs = rna.matureRNAGeneComposition(g.mRNAIndexs(ftsZPol.enzymeMonomerGlobalIndexs), :) ~= 0;
            rnaExp = rna.expression;
            rnaExp(rna.matureIndexs(tfs)) = rnaExp(rna.matureIndexs(tfs)) * ...
                meanFtsZCopyNumber / monCnts(ftsZPol.enzymeMonomerGlobalIndexs);
            rnaExp = rnaExp / sum(rnaExp);
            
            clear sim g time mass rna pm ring ftsZPol;
            
            [~, ~, ~, matureMonCounts, matureCpxCounts] = ProteinGrowth_Test.calcExponentialGrowthDistribution(...
                100, maxTime, 100, 100, 1, 1, [], [], [], rnaExp);
            
            %% assert
            ftsZMonomers = ...
                + matureMonCounts(ftsZMonIdxs, :)' ...
                + matureCpxCounts(ftsZCpxIdxs, :)' * [1 1:9]';
            assertIn(quantile(ftsZMonomers, 0.10), [minFtsZ Inf]);
        end
    end
    
    methods (Static = true)
        function durations = estimateMeanCytokinesisDurationDistribution(nIter, verbosity)
            import edu.stanford.covert.util.ConstantUtil;
            
            %% warnings
            warningState = warning('query', 'WholeCell:warning');
            warning('off', 'WholeCell:warning');
            
            %% setup
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'FtsZPolymerization'});
            g = sim.gene;
            time = sim.state('Time');
            mass = sim.state('Mass');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            ring = sim.state('FtsZRing');
            ftsZPol = sim.process('FtsZPolymerization');
            
            cellCycleLength = time.cellCycleLength;
            expectedDuration = time.cytokinesisDuration;
            maxTime = 5 * expectedDuration;
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            monExp = mRNAExp ./ (log(2) / time.cellCycleLength + pm.decayRates(pm.matureIndexs));
            monExp = monExp / sum(monExp);
            monCnts = mass.initialFractionAAsInMonomers * mass.cellInitialDryWeight * ...
                mass.dryWeightFractionProtein * ConstantUtil.nAvogadro / ...
                (monExp' * pm.molecularWeights(pm.matureIndexs)) * monExp;
            
            meanFtsZCopyNumber = monCnts(ftsZPol.enzymeGlobalIndexs(ftsZPol.enzymeIndexs_FtsZ));
            assertElementsAlmostEqual(...
                4.00 * ring.numEdges * ring.numFtsZSubunitsPerFilament / exp(log(2) * (1 - expectedDuration / cellCycleLength)), ...
                meanFtsZCopyNumber, ...
                'relative', 5e-2);
            
            %clear simulation
            clear sim g time mass rna pm ring ftsZPol;
            
            %% simulate
            durations = zeros(nIter, 1);
            for i = 1:nIter
                if verbosity > 0
                    fprintf('Iter %d\n', i);
                end
                durations(i) = edu.stanford.covert.cell.sim.Cytokinesis_Test.sampleCytokinesisDuration(meanFtsZCopyNumber, expectedDuration, maxTime, i, verbosity);
            end
            
            %% plot
            %hist(durations/3600);
            %xlabel('Time (h)');
            %ylabel('Frequency');
            
            %% warning
            warning(warningState.state, 'WholeCell:warning');
        end
        
        function cytokinesisDuration = sampleCytokinesisDuration(meanFtsZCopyNumber, expectedDuration, maxIter, seed, verbosity)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'ChromosomeSegregation';
                'Cytokinesis';
                'FtsZPolymerization';
                });
            sim.applyOptions('verbosity', verbosity);
            edu.stanford.covert.cell.sim.Cytokinesis_Test.seedRandStream(sim, seed);
            
            time        = sim.state('Time');
            mass        = sim.state('Mass');
            geometry    = sim.state('Geometry');
            met         = sim.state('Metabolite');
            chr         = sim.state('Chromosome');
            rna         = sim.state('Rna');
            pm          = sim.state('ProteinMonomer');
            pc          = sim.state('ProteinComplex');
            ring        = sim.state('FtsZRing');
            segregation = sim.process('ChromosomeSegregation');
            cytokinesis = sim.process('Cytokinesis');
            ftsZPol     = sim.process('FtsZPolymerization');
            
            %% store initial state
            initialWeight = mass.cell;
            initialVolume = geometry.volume;
            
            %% initialize simulation to end of replication
            %cell mass
            f = exp(log(2) * (1 - expectedDuration / time.cellCycleLength));
            met.counts = f * met.counts;
            rna.counts = f * rna.counts;
            pm.counts  = f * pm.counts;
            pc.counts  = f * pc.counts;
            
            %chromosome
            chr.initialize();
            chr.polymerizedRegions(1, :) = chr.polymerizedRegions(1, 1);
            chr.linkingNumbers(1, :)     = chr.linkingNumbers(1, 1);
            
            chr.monomerBoundSites(:, :) = 0;
            chr.complexBoundSites(:, :) = 0;
            
            chr.gapSites(:, :)               = 0;
            chr.abasicSites(:, :)            = 0;
            chr.damagedSugarPhosphates(:, :) = 0;
            chr.damagedBases(:, :)           = 0;
            chr.intrastrandCrossLinks(:, :)  = 0;
            chr.strandBreaks(:, :)           = 0;
            chr.hollidayJunctions(:, :)      = 0;
            
            chr.segregated = false;
            
            %unbind proteins in processes
            for i = 1:numel(sim.processes)
                m = sim.processes{i};
                
                m.copyFromState();
                
                m.enzymes(:) = m.enzymes(:) + m.boundEnzymes(:);
                m.boundEnzymes(:) = 0;
                
                m.copyToState();
            end
            
            %sufficient energy
            met.counts(met.getIndexs('GTP'), sim.compartment.cytosolIndexs) = 1e6;
            met.counts(met.getIndexs('H2O'), sim.compartment.cytosolIndexs) = 1e6;
            
            %sufficient FtsZ
            ftsZPol.copyFromState();
            ftsZPol.enzymes(:) = 0;
            ftsZPol.enzymes(ftsZPol.enzymeIndexs_FtsZ) = ceil(f * meanFtsZCopyNumber);
            ftsZPol.boundEnzymes(:) = 0;
            ftsZPol.copyToState();
            
            cytokinesis.copyFromState();
            
            %% assertions
            %weight, volume double
            f = exp(log(2) * (1 - expectedDuration / time.cellCycleLength));
            assertElementsAlmostEqual(f * initialWeight, mass.cell);
            assertElementsAlmostEqual(f * initialVolume, geometry.volume);
            
            %chromosome not segregated
            assertEqual(false, chr.segregated);
            
            %cell shape defined, not pinched
            assertEqual(geometry.width, geometry.pinchedDiameter);
            assertTrue(geometry.width > 0);
            assertElementsAlmostEqual(0, geometry.cylindricalLength, 'absolute', 1e-7, 0);
            assertTrue(geometry.surfaceArea > 0);
            assertTrue(geometry.totalLength > 0);
            
            %sufficient FtsZ for cytokinesis
            ftsZPolEnzymes = ftsZPol.enzymes + ftsZPol.boundEnzymes;
            ftsZMonomers = ...
                + ftsZPolEnzymes(ftsZPol.enzymeIndexs_FtsZ) ...
                + ftsZPolEnzymes(ftsZPol.enzymeIndexs_FtsZ_GDP) ...
                + ftsZPolEnzymes(ftsZPol.enzymeIndexs_FtsZ_activated)' * (1:9)';
            assertIn(ftsZMonomers, [2.25 * ring.numEdges * ring.numFtsZSubunitsPerFilament Inf]);
            
            %sufficient energy for chromosome segregation and FtsZ polymerization
            [~, ftsZGTP] = cytokinesis.calcRequiredPinchingCycles(geometry.width, ring.numFtsZSubunitsPerFilament, ring.filamentLengthInNm);
            assertIn(met.counts(met.getIndexs('GTP'), sim.compartment.cytosolIndexs), [segregation.gtpCost + ftsZGTP Inf]);
            
            %% run
            pinchedDiameter = geometry.pinchedDiameter;
            nCycle = 1;
            cytokinesisDuration = NaN;
            
            if sim.verbosity > 0
                fprintf('\t%4s %4s %4s %4s %7s %4s %4s %4s %4s %4s %11s\n', ...
                    'Time', 'Segr', ...
                    'FtsZ', '9mer', 'Bnd9mer', ...
                    '1Srt', '2Srt', '2Bnt', 'RBnt', ...
                    'Cycl', 'Diameter   ');
                fprintf('\t%4s %4s %4s %4s %7s %4s %4s %4s %4s %4s %11s\n', ...
                    '====', '====', ...
                    '====', '====', '=======', ...
                    '====', '====', '====', '====', ...
                    '====', '===========');
                edu.stanford.covert.cell.sim.Cytokinesis_Test.printState(sim, 0, 1);
            end
            for iter = 1:maxIter
                sim.evolveState();
                
                if pinchedDiameter ~= geometry.pinchedDiameter
                    nCycle = nCycle + 1;
                end
                
                pinchedDiameter = geometry.pinchedDiameter;
                
                if sim.verbosity > 0 && mod(iter, 100) == 1
                    edu.stanford.covert.cell.sim.Cytokinesis_Test.printState(sim, iter, nCycle);
                end
                
                if chr.segregated && geometry.pinched && ring.numResidualBent == 0 && geometry.width == -1
                    cytokinesisDuration = iter;
                    if sim.verbosity > 0
                        fprintf('\tCell pinched in %d s\n', iter);
                    end
                    break;
                end
            end
            
            %% assert
            %chromosome segregated
            assertTrue(chr.segregated);
            
            %completely pinched
            assertEqual(0, geometry.pinchedDiameter);
            assertTrue(geometry.pinched);
            
            %cell shape undefined, completely pinched
            assertEqual(0, geometry.pinchedDiameter);
            assertEqual(-1, geometry.width);
            assertEqual(-1, geometry.cylindricalLength);
            assertEqual(-1, geometry.surfaceArea);
            assertEqual(-1, geometry.totalLength);
        end
    end
    
    %helper methods
    methods (Static = true)
        function printState(sim, i, nCycle)
            chr         = sim.state('Chromosome');
            geometry    = sim.state('Geometry');
            ring        = sim.state('FtsZRing');
            cytokinesis = sim.process('Cytokinesis');
            ftsZPol     = sim.process('FtsZPolymerization');
            
            ftsZPolEnzymes = ftsZPol.enzymes + ftsZPol.boundEnzymes;
            ftsZMonomers = ...
                + ftsZPolEnzymes(ftsZPol.enzymeIndexs_FtsZ) ...
                + ftsZPolEnzymes(ftsZPol.enzymeIndexs_FtsZ_GDP) ...
                + ftsZPolEnzymes(ftsZPol.enzymeIndexs_FtsZ_activated(1:end-1))' * (1:8)';
            fprintf('\t%4d %4d %4d %4d %7d %4d %4d %4d %4d %4d %.4e\n', ...
                i, chr.segregated, ...
                ftsZMonomers, ...
                cytokinesis.enzymes(cytokinesis.enzymeIndexs_ftsZ_GTP_polymer),...
                cytokinesis.boundEnzymes(cytokinesis.enzymeIndexs_ftsZ_GTP_polymer), ...
                ring.numEdgesOneStraight, ring.numEdgesTwoStraight,...
                ring.numEdgesTwoBent, ring.numResidualBent, ...
                nCycle, geometry.pinchedDiameter);
        end
        
        function seedRandStream(sim, seed)
            if sim.verbosity > 0
                fprintf('\tSeed: %d\n', seed);
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
