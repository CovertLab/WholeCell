%DNA medium test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 4/19/2011
classdef DNADamageRepair_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = DNADamageRepair_Test(name)
            this = this@TestCase(name);
        end
        
        function setUp(this)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'DNADamage'
                'DNARepair'
                });
            sim.applyOptions('verbosity', 0, 'seed', 1);
            
            this.simulation = sim;
        end
        
        function testExpectedDamageRates_sampleAccessibleSites(this)
            %references
            sim = this.simulation;
            c = sim.state('Chromosome');
            
            %test sampleAccessibleSites
            iterMax = 1000;
            cnts = zeros(iterMax, 1);
            prob = 1e-6;
            for i = 1:iterMax
                cnts(i) = size(c.sampleAccessibleSites(prob, Inf, 'T'), 1);
            end
            assertEqual([0; 1], unique(cnts));
            assertElementsAlmostEqual(...
                sum(cnts), ...
                iterMax * prob * 2 * c.sequenceLen * (1-c.sequenceGCContent)/2, ...
                'relative', 0.10);
        end
        
        function testExpectedDamageRates_setSiteDamaged(this)
            %references
            sim = this.simulation;
            c = sim.state('Chromosome');
            
            %test setSiteDamaged
            iterMax = 1000;
            cnts = zeros(iterMax, 1);
            prob = 1e-6;
            for i = 1:iterMax
                c.initialize();
                cnts(i) = size(c.setSiteDamaged('damagedBases', 1, prob, Inf, 'T', 'sequence'), 1);
            end
            assertEqual([0; 1], unique(cnts));
            assertElementsAlmostEqual(...
                sum(cnts), ...
                iterMax * prob * 2 * c.sequenceLen * (1-c.sequenceGCContent)/2, ...
                'relative', 0.10);
        end
        
        function testExpectedDamageRates_evolveState_simpleReactionSet(this)
            %references
            sim = this.simulation;
            c = sim.state('Chromosome');
            d = sim.process('DNADamage');
            
            %test evolve state with simplified reaction set
            iterMax = 1000;
            
            d.reactionWholeCellModelIDs = d.reactionWholeCellModelIDs(1:4);
            d.reactionStoichiometryMatrix = d.reactionStoichiometryMatrix(:, 1:4);
            d.reactionDamageTypes = {'damagedBases'; 'damagedBases'; 'intrastrandCrossLinks'; 'intrastrandCrossLinks'};
            d.reactionBounds = [zeros(4, 1) repmat(1e-6, 4, 1)];
            d.reactionRadiation = zeros(4, 1);
            d.reactionDNAProduct = (1:4)';
            d.reactionVulnerableMotifs = {'A'; 'C'; 'G'; 'T'};
            d.reactionVulnerableMotifTypes = repmat({'sequence'}, 4, 1);
            
            c.initialize();
            expRates = d.calcExpectedReactionRates() * iterMax;
            
            d.substrates(any(d.reactionSmallMoleculeStoichiometryMatrix,2)) = 1e6;
            for i = 1:iterMax
                d.evolveState();
            end
            [~, vals] = find(c.damagedBases);
            assertElementsAlmostEqual(histc(vals, d.reactionDNAProduct), expRates .* [1; 1; 0; 0], ...
                'relative', 0.10);
            [~, vals] = find(c.intrastrandCrossLinks);
            assertElementsAlmostEqual(histc(vals, d.reactionDNAProduct), expRates .* [0; 0; 1; 1], ...
                'relative', 0.10);
        end
        
        function testExpectedDamageRates_evolveState_lessSimpleReactionSet(this)
            %references
            sim = this.simulation;
            c = sim.state('Chromosome');
            m = sim.state('Metabolite');
            d = sim.process('DNADamage');
            
            %simulate damage
            iterMax = 1000;
            
            d.reactionBounds(d.reactionBounds(:, 2) < 1000, 2) = 1e-6;
            d.substrates(d.reactionRadiation(d.reactionRadiation ~= 0)) = 1;
            d.reactionBounds(~cellfun(@ischar, d.reactionVulnerableMotifs), 2) = 0;
            
            expRates = d.calcExpectedReactionRates() * iterMax;
            
            assertFalse(any(d.reactionDNAProduct == m.m6ADIndexs & strcmp(d.reactionDamageTypes, 'damagedBases')));
            
            d.substrates(any(d.reactionSmallMoleculeStoichiometryMatrix,2)) = 1e6;
            for i = 1:iterMax
                d.evolveState();
            end
            
            %assert damage rates ~= expected
            this.assertDamageRateIsExpected(c, d, expRates, iterMax, 0.4, 0.5);
        end
        
        function testExpectedDamageRates_evolveState_fullReactionSet(this)
            %references
            sim = this.simulation;
            c = sim.state('Chromosome');
            m = sim.state('Metabolite');
            s = sim.state('Stimulus');
            d = sim.process('DNADamage');
            r = sim.process('DNARepair');
            
            %simulate damage
            iterMax = 1000;
            
            d.reactionBounds(d.reactionBounds(:, 2) < 1000, 2) = 1e-6;
            d.substrates(:) = 1e12;
            d.substrates(d.reactionRadiation(d.reactionRadiation ~= 0)) = 1;
            s.setValues(s.getIndexs(d.substrateWholeCellModelIDs(d.reactionRadiation(d.reactionRadiation ~= 0))), 3) = 1;
            
            expRates = d.calcExpectedReactionRates() * iterMax;
            
            assertFalse(any(d.reactionDNAProduct == m.m6ADIndexs & ...
                strcmp(d.reactionDamageTypes, 'damagedBases')));
            
            d.reactionSmallMoleculeStoichiometryMatrix(:) = 0;
            
            r.reactionBounds(:, 2) = 0;
            r.enzymeBounds(:, 2) = 0;
            
            sim.applyOptions('lengthSec', iterMax);
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            sim.run();
            
            %assert damage rates ~= expected
            this.assertDamageRateIsExpected(c, d, expRates, iterMax, 0.1, 10);
        end
        
        function testAllDamageCanBeRepaired(this)
            sim = this.simulation;
            c = sim.state('Chromosome');
            d = sim.process('DNADamage');
            r = sim.process('DNARepair');
            
            repairRates = 10000 * r.enzymeBounds;
            
            for i = 1:numel(d.reactionWholeCellModelIDs)
                c.initialize();
                d.initializeState();
                r.initializeState();
                
                sim.state('Metabolite').counts(:, sim.compartment.cytosolIndexs) = 1e6;
                sim.state('Stimulus').values(:) = 0;
                sim.state('Stimulus').setValues(:, 3) = 0;
                
                d.reactionBounds(:, 2) = 0;
                
                switch d.reactionVulnerableMotifTypes{i}
                    case 'sequence'
                        d.reactionBounds(i, 2) = 3e-5 * 4 ^ (numel(d.reactionVulnerableMotifs{i})-1);
                    case {'damagedSugarPhosphate', 'damagedBases', 'intrastrandCrossLinks'}
                        c.(d.reactionVulnerableMotifTypes{i})((1:5)*10000, 1) = ...
                            d.reactionVulnerableMotifs{i};
                        c.(d.reactionVulnerableMotifTypes{i})((6:10)*10000, 2) = ...
                            d.reactionVulnerableMotifs{i};
                        d.reactionBounds(i, 2) = 1;
                    case {'gapSites', 'abasicSites', 'strandBreaks', 'hollidayJunctions'}
                        c.(d.reactionVulnerableMotifTypes{i})((1:5) * 10000, 1) = 1;
                        c.(d.reactionVulnerableMotifTypes{i})((6:10) * 10000, 2) = 1;
                        d.reactionBounds(i, 2) = 1;
                    otherwise, throw(MException('DNADamageRepair_Test:error', ...
                            'Unimplemented vulnerable motif %s', d.reactionVulnerableMotifTypes{i}));
                end
                
                if d.reactionRadiation(i)
                    sim.state('Stimulus').setValues(sim.state('Stimulus').setValues(:,1) == ...
                        d.substrateGlobalIndexs(d.reactionRadiation(i)), 3) = 1;
                end
                
                r.enzymeBounds(:) = 0;
                
                %stochastically damage
                for j = 1:2
                    sim.evolveState();
                end
                
                switch d.reactionDamageTypes{i}
                    case {'damagedSugarPhosphate', 'damagedBases', 'intrastrandCrossLinks'}
                        reactionDNAProduct = d.reactionDNAProduct(i);
                    case {'gapSites', 'abasicSites', 'strandBreaks', 'hollidayJunctions'}
                        reactionDNAProduct = 1;
                    otherwise, throw(MException('DNADamageRepair_Test:error', ...
                            'Unimplemented damage type %s', d.reactionDamageTypes{i}));
                end
                assertIn(nnz(c.(d.reactionDamageTypes{i}) == reactionDNAProduct), [2 50], ...
                    sprintf('%d %s %s\n', i, d.reactionWholeCellModelIDs{i}, d.reactionDamageTypes{i}));
                
                %repair
                d.reactionBounds(:, 2) = 0;
                r.enzymeBounds = repairRates;
                
                posStrnds = find(any(c.(d.reactionDamageTypes{i}) == reactionDNAProduct, 2));
                nDamagesMax = sum(diff([posStrnds(:, 1); posStrnds(1,1) + c.sequenceLen]) <= max(r.enzymeDNAFootprints));
                
                initDamagedSites = nnz(c.damagedSites);
                
                for j = 1:350
                    sim.evolveState();
                    
                    damagedSites = c.damagedSites;
                    damagedSites(find(c.damagedBases)) = 0; %#ok<FNDSB>
                    damagedSites(c.shiftPositionsStrandsBase3Prime(find(c.damagedBases))) = 0; %#ok<FNDSB>
                    damagedSites(c.shiftPositionsStrandsBase5Prime(find(c.damagedBases))) = 0; %#ok<FNDSB>
                    if nnz(damagedSites) <= nDamagesMax
                        break;
                    end
                end
                
                assertIn(nnz(damagedSites), [0 nDamagesMax], ...
                    sprintf('%d %s %d\n', i, d.reactionWholeCellModelIDs{i}, initDamagedSites));
            end
        end
        
        function testDNARepairOutweighsDamage(this)
            %% references
            sim = this.simulation;
            c = sim.state('Chromosome');
            d = sim.process('DNADamage');
            r = sim.process('DNARepair');
            
            c.initialize();
            r.initializeState();
            
            %% indices of proteins which should never bind
            monGblIdxs = [];
            cpxGblIdxs = [];
            for i = 1:numel(sim.processes)
                p = sim.processes{i};
                monGblIdxs = [monGblIdxs; ...
                    p.enzymeGlobalIndexs(p.enzymeMonomerLocalIndexs)]; %#ok<*AGROW>
                cpxGblIdxs = [cpxGblIdxs; ...
                    p.enzymeGlobalIndexs(p.enzymeComplexLocalIndexs)];
            end
            monGblIdxs = sort(monGblIdxs);
            cpxGblIdxs = sort(cpxGblIdxs);
            
            monGblIdxs = setdiff(monGblIdxs, r.enzymeGlobalIndexs(r.enzymeIndexs_recombinationStrandExchange));
            cpxGblIdxs = setdiff(cpxGblIdxs, r.enzymeGlobalIndexs(r.enzymeIndexs_DisA));
            
            %% indices of allowed damage
            dmgSPIdxs = [];
            dmgBaseIdxs = sort([
                d.reactionDNAProduct(ismember(d.reactionDamageTypes, {'damagedBases'}))
                r.substrateGlobalIndexs(r.substrateIndexs_m6AD)
                ]);
            dmgCLIdxs = sort(...
                d.reactionDNAProduct(ismember(d.reactionDamageTypes, {'intrastrandCrossLinks'})));
            
            %% simulate
            iterMax = 250;
            
            d.reactionBounds = 1e4 * d.reactionBounds;
            r.enzymeBounds = 1e4 * r.enzymeBounds;
            
            sim.state('Metabolite').counts(:, sim.compartment.cytosolIndexs) = 1e6;
            
            nDamages = zeros(iterMax, 1);
            for i = 1:iterMax
                this.assertSecondChromosomeIgnored(c);
                this.assert1ChromsomePolymerized(c);
                this.assertCorrectProteinsBound(c, monGblIdxs, cpxGblIdxs);
                this.assertCorrectDamages(c, dmgSPIdxs, dmgBaseIdxs, dmgCLIdxs);
                
                sim.evolveState();
                
                nDamages(i) = nnz(c.damagedSites);
            end
            
            this.assertDamageDoesntAccumulate(c, r);
            
            d.reactionBounds(:, 2) = 0;
            
            for i = 1:20
                sim.evolveState();
            end
            
            assertEqual(0, nnz(c.damagedSites));
        end
    end
    
    %helper methods
    methods
        %check that second chromosome is never modified
        function assertSecondChromosomeIgnored(~, c)
            posStrnds = find(c.polymerizedRegions);     assertFalse(any(posStrnds(:, 2) > 2));
            posStrnds = find(c.linkingNumbers);         assertFalse(any(posStrnds(:, 2) > 2));
            
            posStrnds = find(c.gapSites);               assertFalse(any(posStrnds(:, 2) > 2));
            posStrnds = find(c.abasicSites);            assertFalse(any(posStrnds(:, 2) > 2));
            posStrnds = find(c.damagedSugarPhosphates); assertFalse(any(posStrnds(:, 2) > 2));
            posStrnds = find(c.damagedBases);           assertFalse(any(posStrnds(:, 2) > 2));
            posStrnds = find(c.strandBreaks);           assertFalse(any(posStrnds(:, 2) > 2));
            posStrnds = find(c.intrastrandCrossLinks);  assertFalse(any(posStrnds(:, 2) > 2));
            posStrnds = find(c.hollidayJunctions);      assertFalse(any(posStrnds(:, 2) > 2));
            
            posStrnds = find(c.monomerBoundSites);      assertFalse(any(posStrnds(:, 2) > 2));
            posStrnds = find(c.complexBoundSites);      assertFalse(any(posStrnds(:, 2) > 2));
            
            assertFalse(c.segregated);
        end
        
        %check that only first chromosome is polymerized
        function assert1ChromsomePolymerized(~, c)
            import edu.stanford.covert.util.CircularSparseMat;
            
            lk0 = c.sequenceLen / c.relaxedBasesPerTurn * (1 + c.equilibriumSuperhelicalDensity);
            assertEqual(CircularSparseMat([1 1; 1 2], [c.sequenceLen; c.sequenceLen], [c.sequenceLen c.nCompartments], 1), c.polymerizedRegions);
            assertEqual(CircularSparseMat([1 1; 1 2], [lk0; lk0], [c.sequenceLen c.nCompartments], 1), c.linkingNumbers);
        end
        
        %check that only allowed proteins ever bind chromosome
        function assertCorrectProteinsBound(~, c, monGblIdxs, cpxGblIdxs)
            [~, monVals] = find(c.monomerBoundSites);
            [~, cpxVals] = find(c.complexBoundSites);
            
            assertFalse(any(ismembc(monVals, monGblIdxs)));
            assertFalse(any(ismembc(cpxVals, cpxGblIdxs)));
        end
        
        %check damage is of allowed types
        function assertCorrectDamages(~, c, dmgSPIdxs, dmgBaseIdxs, dmgCLIdxs)
            [~, vals] = find(c.gapSites);               assertTrue(all(ismembc(vals, [false; true])));
            [~, vals] = find(c.abasicSites);            assertTrue(all(ismembc(vals, [false; true])));
            [~, vals] = find(c.damagedSugarPhosphates); assertTrue(all(ismembc(vals, dmgSPIdxs)));
            [~, vals] = find(c.damagedBases);           assertTrue(all(ismembc(vals, dmgBaseIdxs)));
            [~, vals] = find(c.strandBreaks);           assertTrue(all(ismembc(vals, [false; true])));
            [~, vals] = find(c.intrastrandCrossLinks);  assertTrue(all(ismembc(vals, dmgCLIdxs)));
            [~, vals] = find(c.hollidayJunctions);      assertTrue(all(ismembc(vals, [false; true])));
        end
        
        %check damage doesn't accumulate -- damage that occurs is repaired quickly
        function assertDamageDoesntAccumulate(~, c, r)
            assertIn(nnz(c.gapSites), [0 10 + r.NER_UvrABC_IncisionMargin3 + r.NER_UvrABC_IncisionMargin5]);
            assertIn(nnz(c.abasicSites), [0 10 + r.NER_UvrABC_IncisionMargin3 + r.NER_UvrABC_IncisionMargin5]);
            assertIn(nnz(c.damagedSugarPhosphates), [0 10]);
            assertIn(nnz(c.damagedBases), [0 10 + 2 * size(r.RM_MunI_RecognitionSites, 1)]);
            assertIn(nnz(c.strandBreaks), [0 10]);
            assertIn(nnz(c.intrastrandCrossLinks), [0 20]);
            assertIn(nnz(c.hollidayJunctions), [0 20]);
        end
        
        %check that damage rate follow expectations
        function assertDamageRateIsExpected(~, c, d, expRates, iterMax, relTol, absTol)
            [types, ~, idxs] = unique(d.reactionDamageTypes);
            reactionDNAProduct = d.reactionDNAProduct;
            reactionDNAProduct(ismember(d.reactionDamageTypes, ...
                {'strandBreaks','abasicSites','gapSites','hollidayJunctions'})) = 1;
            dmgTypesProds = unique([idxs reactionDNAProduct], 'rows');
            calcRates = zeros(size(dmgTypesProds, 1), 1);
            obsRates = zeros(size(dmgTypesProds, 1), 1);
            for i = 1:size(dmgTypesProds, 1)
                rxnTfs = ...
                    ismember(d.reactionDamageTypes, types(dmgTypesProds(i,1))) & ...
                    reactionDNAProduct == dmgTypesProds(i, 2);
                calcRates(i) = sum(expRates(rxnTfs));
                obsRates(i) = nnz(c.(types{dmgTypesProds(i, 1)}) == dmgTypesProds(i, 2));
            end
            
            reactionVulnerableMotifs = d.reactionVulnerableMotifs;
            reactionVulnerableMotifs(cellfun(@ischar, reactionVulnerableMotifs)) = {0};
            reactionVulnerableMotifs(cellfun(@islogical, reactionVulnerableMotifs)) = {1};
            reactionVulnerableMotifs = cell2mat(reactionVulnerableMotifs);
            rates = d.reactionBounds(:, 2);
            rates(d.reactionRadiation~=0) = ...
                rates(d.reactionRadiation~=0) .* ...
                d.substrates(d.reactionRadiation(d.reactionRadiation~=0));
            for i = 1:size(dmgTypesProds, 1)
                j = find(...
                    rates > 0 & ...
                    ismember(reactionVulnerableMotifs, dmgTypesProds(i, 2)) & ...
                    ismember(d.reactionVulnerableMotifTypes, types(dmgTypesProds(i, 1))));
                
                [~, k] = ismember([idxs(j) reactionDNAProduct(j)], dmgTypesProds, 'rows');
                
                if isempty(j)
                    continue;
                end
                
                calcRates(k) = calcRates(k) + min(1, iterMax * rates(j)) * calcRates(i);
                if all(dmgTypesProds(k,1) == dmgTypesProds(i,1))
                    calcRates(i) = max(0, 1 - iterMax * sum(rates(j))) * calcRates(i);
                end
            end
            
            %metIds = d.metabolite.wholeCellModelIDs(dmgTypesProds(:,2));
            %tfs =ismember(types(dmgTypesProds(:,1)), {'strandBreaks','abasicSites','gapSites','hollidayJunctions'});
            %metIds(tfs) = cell(sum(tfs), 1);
            %[types(dmgTypesProds(:,1)) metIds num2cell(calcRates) num2cell(obsRates)]
            
            assertElementsAlmostEqual(calcRates', obsRates', 'relative', relTol, absTol);
        end
    end
end
