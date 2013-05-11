function this = initializeState(this)
% Initialize state of organism to expectation values.
% - Time
% - Numbers of metabolites, RNAs, protein monomers, protein complexes
% - Resource allocation
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/17/2009

%Todo: better match initial conditions with final conditions of simulations
%- 46 immature (folded) molecules MG-272 monomer
%- 1 immature (folded) molecule MG-166 monomer
%- 2.2 nascent RNA
%- 2.2 processed RNA
%- Fewer RNA pol complexes (more uncomplexed monomers)

%% import classes
import edu.stanford.covert.util.ComputationUtil;
import edu.stanford.covert.util.ConstantUtil;
import edu.stanford.covert.cell.sim.constant.Condition;
import edu.stanford.covert.cell.sim.Simulation;

%% data
%states
g = this.gene;
t = this.state('Time');
s = this.state('Stimulus');
m = this.state('Metabolite');
c = this.state('Chromosome');
r = this.state('Rna');
pm = this.state('ProteinMonomer');
pc = this.state('ProteinComplex');
mass = this.state('Mass');
mr = this.state('MetabolicReaction');
geom = this.state('Geometry');
ftsZRing = this.state('FtsZRing');
transcript = this.state('Transcript');
polypeptide = this.state('Polypeptide');
dnaRepair = this.process('DNARepair');

%indices
cIdx = this.compartment.cytosolIndexs;

%create local variables for convenience
rnaMWs = 1 / ConstantUtil.nAvogadro * r.molecularWeights;
monMWs = 1 / ConstantUtil.nAvogadro * pm.molecularWeights;

chrComp = [getBaseCounts(c.sequence); 0];
if ~isempty(dnaRepair)
    m6dAMPComp = 2 * size(dnaRepair.RM_MunI_RecognitionSites, 1);
    chrComp(5) = m6dAMPComp;
    chrComp(1) = chrComp(1) - m6dAMPComp;
end
rnaComp = r.baseCounts;
monComp = pm.baseCounts;

rnaExp = r.expression;

%% seed random number generator
this.seedRandStream();

%% clear simulated properties
this.allocateMemoryForState(1);

%% time
t.initialize();

%% total cell mass
totCellWt = this.randStream.random('norm', 1, mass.cellInitialMassVariation);

bComp = totCellWt * m.biomassComposition;
rnaWt = totCellWt * mass.cellInitialDryWeight * mass.dryWeightFractionRNA;
prtWt = totCellWt * mass.cellInitialDryWeight * mass.dryWeightFractionProtein;

%% metabolites
m.counts = max(0, this.randStream.stochasticRound(bComp));

%% DNA
%Initialize the free dNTP pool. Chromosome is initialized below in
%initialize states section.
c.initialize();
m.counts([m.dnmpIndexs; m.getIndexs('m6dAMP')], cIdx) = ...
    m.counts([m.dnmpIndexs; m.getIndexs('m6dAMP')], cIdx) - chrComp;

%% RNAs
%initialize RNAs based on expected expression
switch this.macromoleculeStateInitialization
    case 'expected'
        r.counts(:, cIdx) = this.randStream.stochasticRound(...
            mass.initialFractionNTPsInRNAs * rnaWt * rnaExp / (rnaMWs' * rnaExp));
    case 'multinomial'
        prob = rnaExp;
        idxs = this.randStream.randsample(numel(prob), ceil(2 * rnaWt / (rnaMWs' * prob)), true, prob);
        idxs = idxs(1:find(cumsum(rnaMWs(idxs)) > mass.initialFractionNTPsInRNAs * rnaWt, 1, 'first') - 1);
        r.counts(:, cIdx) = this.randStream.stochasticRound(histc(idxs, 1:numel(prob)));
    otherwise, throw(MException('Simulation:error', 'Unsupported macromoleculeStateInitialization %s', this.macromoleculeStateInitialization));
end
m.counts = m.counts - rnaComp' * r.counts;

%% protein monomers
%1. initialize protein monomers based on expected expression
%2. sort to appropriate compartment
switch this.macromoleculeStateInitialization
    case 'expected'
        monExp = zeros(size(pm.wholeCellModelIDs));
        monExp(pm.matureIndexs) = (r.matureRNAGeneComposition(g.mRNAIndexs, :) * r.expression(r.matureIndexs)) ./ ...
            (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
        monExp = monExp / sum(monExp);
        
        pm.counts(sub2ind(size(pm.counts), (1:size(pm.counts,1))', pm.compartments)) = ...
            this.randStream.stochasticRound(mass.initialFractionAAsInMonomers * prtWt * monExp / (monMWs' * monExp));
    case 'multinomial'
        n = pm.macromoleculeStateInitializationVariation;
        
        if rnaWt > 0
            mRNAExp = rnaExp(r.matureIndexs(r.matureMRNAIndexs));
            mRNAExp = mRNAExp / sum(mRNAExp);
            mRNAProd = mRNAExp .* r.decayRates(r.matureIndexs(r.matureMRNAIndexs));
            mRNAProd = mRNAProd / sum(mRNAProd);
            mRNACnt = mRNAExp * rnaWt * r.weightFractionMRNA * ConstantUtil.nAvogadro / ...
                (mRNAExp' * r.molecularWeights(r.matureIndexs(r.matureMRNAIndexs)));
            totMRNAProd = (sum(mRNACnt) * log(2) / t.cellCycleLength + r.decayRates(r.matureIndexs(r.matureMRNAIndexs))' * mRNACnt) * t.cellCycleLength / log(2);
            
            sample_mRNAProd = this.randStream.mnrnd(round(n * totMRNAProd), mRNAProd)';
            sample_mRNAExp = zeros(size(sample_mRNAProd));
            for i = 1:numel(sample_mRNAExp)
                if sample_mRNAProd(i) > 0
                    sample_mRNAExp(i) = sample_mRNAProd(i) * ...
                        mean(this.randStream.random('exponential', ...
                        1 / r.decayRates(r.matureIndexs(r.matureMRNAIndexs(i))), ...
                        [sample_mRNAProd(i) 1]));
                end
            end
        else
            sample_mRNAExp = r.expression(r.matureIndexs(r.matureMRNAIndexs));
            sample_mRNAExp = sample_mRNAExp / sum(sample_mRNAExp);
        end
        
        if prtWt > 0
            monProd = r.matureRNAGeneComposition(g.mRNAIndexs, r.matureMRNAIndexs) * sample_mRNAExp;
            monProd = monProd / sum(monProd);
            monExp = monProd ./ (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
            monExp = monExp / sum(monExp);
            monCnt = prtWt * ConstantUtil.nAvogadro / (monExp' * pm.molecularWeights(pm.matureIndexs)) * monExp;
            totMonProd = monCnt' * (1 + t.cellCycleLength / log(2) * pm.decayRates(pm.matureIndexs));
            
            sample_monExp = this.randStream.mnrnd(round(n * totMonProd), monExp)';
            sample_monExp = this.randStream.stochasticRound(n * sample_monExp * prtWt * mass.initialFractionAAsInMonomers * ...
                ConstantUtil.nAvogadro / (sample_monExp' * pm.molecularWeights(pm.matureIndexs)));
            sample_monExp2 = zeros(size(sample_monExp));
            tmp = cumsum(sample_monExp);
            order = this.randStream.randperm(tmp(end));
            for i = 1:tmp(end)/2
                idx = find(order(i) <= tmp, 1, 'first');
                sample_monExp2(idx) = sample_monExp2(idx) + 1;
            end
            pm.counts(sub2ind(size(pm.counts), pm.matureIndexs, pm.compartments(pm.matureIndexs))) = ...
                this.randStream.stochasticRound(2/n * sample_monExp2);
        end
    otherwise, throw(MException('Simulation:error', 'Unsupported macromoleculeStateInitialization %s', this.macromoleculeStateInitialization));
end
m.counts(:, cIdx) = m.counts(:, cIdx) - monComp' * sum(pm.counts, 2);

if strcmp(this.macromoleculeStateInitialization, 'expected')
    pcComp = sum(pc.proteinComplexComposition, 3);
    nonFormCpxIdxs = find(pc.formationProcesses(pc.matureIndexs) ~= this.processIndex('MacromolecularComplexation'));
    
    subunits = zeros(size(g.wholeCellModelIDs));
    subunits(g.mRNAIndexs) = sum(pm.counts(pm.matureIndexs, :), 2);
    subunits(g.rRNAIndexs) = sum(r.counts(r.matureRRNAIndexs, :), 2);
    subunits(g.sRNAIndexs) = sum(r.counts(r.matureSRNAIndexs, :), 2);
    subunits(g.tRNAIndexs) = sum(r.counts(r.matureTRNAIndexs, :), 2);
    cpxs = floor(min(repmat(subunits, 1, numel(pc.matureIndexs)) ./ pcComp, [], 1))';
    cpxs(nonFormCpxIdxs) = 0;
    
    [subs2Nets, cpxs2Nets, nets] = edu.stanford.covert.util.findNonInteractingRowsAndColumns(pcComp);
    for i = 1:numel(nets)
        if size(nets{i}, 2) <= 1
            continue;
        end
        
        tmpSubIdxs = find(subs2Nets == i);
        tmpCpxIdxs = find(cpxs2Nets == i);
        tmpPcComp = pcComp(tmpSubIdxs, tmpCpxIdxs);
        tmpSubs = subunits(tmpSubIdxs);
        tmpCpxs = zeros(size(tmpCpxIdxs));
        while true
            if all(sum(tmpCpxIdxs, 1) == sum(tmpCpxIdxs(:, 1)))
                tmpRates = prod(repmat(tmpSubs / mean(tmpSubs), 1, numel(tmpCpxIdxs)) .^ tmpPcComp, 1)';
            else
                tmpRates = prod(repmat(tmpSubs / sum(totMons), 1, numel(tmpCpxIdxs)) .^ tmpPcComp, 1)';
            end
            tmpRates(ismember(tmpCpxIdxs, nonFormCpxIdxs)) = 0;
            if ~any(tmpRates) || max(min(repmat(tmpSubs, 1, numel(tmpCpxIdxs)) ./ tmpPcComp, [], 1)) < 1e-3
                break;
            end
            tmpCpxs = tmpCpxs + tmpRates * min(tmpSubs ./ (tmpPcComp * tmpRates));
            tmpSubs = subunits(tmpSubIdxs) - tmpPcComp * tmpCpxs;
        end
        cpxs(tmpCpxIdxs) = floor(tmpCpxs);
    end
    
    pc.counts(sub2ind(...
        size(pc.counts), ...
        pc.matureIndexs, ...
        pc.compartments(pc.matureIndexs))) = ...
        cpxs;
    
    r.counts(r.matureIndexs(setdiff(1:end, r.matureMRNAIndexs)), :) = ...
        r.counts(r.matureIndexs(setdiff(1:end, r.matureMRNAIndexs)), :) - ...
        sum(permute(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), [1 3 2]) .* ...
        repmat(permute(cpxs, [2 3 1]), [numel(g.rRNAIndexs)+numel(g.sRNAIndexs)+numel(g.tRNAIndexs) size(pc.proteinComplexComposition, 3) 1]), 3);
    pm.counts(pm.matureIndexs, :) = ...
        pm.counts(pm.matureIndexs, :) - ...
        sum(permute(pc.proteinComplexComposition(g.mRNAIndexs, :, :), [1 3 2]) .* ...
        repmat(permute(cpxs, [2 3 1]), [numel(g.mRNAIndexs) size(pc.proteinComplexComposition, 3) 1]), 3);
end

%% nucleotides
V = totCellWt * mass.cellInitialDryWeight / (1 - mass.fractionWetWeight) / geom.density; %volume (L)
m.counts(m.ntpIndexs([1 3]), cIdx)  = round(1e-3 * ConstantUtil.nAvogadro * V * m.meanNTPConcentration);
m.counts(m.ntpIndexs([2 4]), cIdx)  = 0;
m.counts(m.ndpIndexs([1 3]), cIdx)  = round(1e-3 * ConstantUtil.nAvogadro * V * m.meanNDPConcentration);
m.counts(m.ndpIndexs([2 4]), cIdx)  = 0;
m.counts(m.nmpIndexs, cIdx)         = round(1e-3 * ConstantUtil.nAvogadro * V * m.meanNMPConcentration);
m.counts(m.phosphateIndexs, cIdx)   = sum(m.counts(m.ndpIndexs, cIdx));
m.counts(m.diphosphateIndexs, cIdx) = sum(m.counts(m.nmpIndexs, cIdx));
m.counts(m.hydrogenIndexs, cIdx)    = sum(m.counts(m.ndpIndexs, cIdx)) + sum(m.counts(m.nmpIndexs, cIdx));

%% protein complexes: initialized by processes (see initialize processes section below)
%1. Macromolecular complexation
%2. Ribosome assembly
%3. Protein folding
%4. Protein activation

%% initialize states
mass.initialize();
geom.initialize();
ftsZRing.initialize();
transcript.initialize();
polypeptide.initialize();

%% apply media and stimulation conditions
m.counts = Condition.applyConditions(m.counts, m.setCounts, t.values); %media
s.values = Condition.applyConditions(s.values, s.setValues, t.values); %stimulation

%% Assert metabolites non-negative
negIdxs = find(any(m.counts < min(0, -bComp) * 15e-2, 2));
if ~isempty(negIdxs)
    warnStr = cellfun(@(x, y) sprintf('- %s: %.3e', x, y), m.wholeCellModelIDs(negIdxs), num2cell(m.counts(negIdxs)), 'UniformOutput', false);
    warning('WholeCell:warning:initializeNegativeCounts', ...
        ['Metabolites should be non-negative:\n' repmat('%s\n', 1, numel(negIdxs))], warnStr{:});
end
m.counts = max(0, m.counts);

%% initialize processes
%- Initialize "own" state while preserving the cell mass
%- Processes "own" state includes, for example,
%  - RNA polymerases in transcription
%  - Ribosomes in translation
%- In particular, when processes incorporate free metabolites into the cell's mass
%  they must alter the counts of simulation's free metabolites
%- Processes should assume that energy is not limited, and should not account for
%  energy usage. Specifically, processes should not alter the simulation's
%  metabolite counts to account for energy usage.
processes = this.processesInInitOrder;
for i = 1:length(processes)
    process = processes{i};
    process.simulationStateSideEffects = [];
    process.copyFromState();
    process.initializeState();
    process.copyToState();
    if ~isempty(process.simulationStateSideEffects)
        process.simulationStateSideEffects.updateSimulationState(this);
    end
end
this.state('Polypeptide').abortedPolypeptides = this.state('Polypeptide').abortedPolypeptides([], :);
this.state('Transcript').abortedTranscripts = this.state('Transcript').abortedTranscripts([], :);

%% Assert metabolites non-negative
negIdxs = find(any(m.counts < min(0, -bComp) * 15e-2, 2));
if ~isempty(negIdxs)
    warnStr = cellfun(@(x, y) sprintf('- %s: %.3e', x, y), m.wholeCellModelIDs(negIdxs), num2cell(m.counts(negIdxs)), 'UniformOutput', false);
    warning('WholeCell:warning', ...
        ['Metabolites should be non-negative:\n' repmat('%s\n', 1, numel(negIdxs))], warnStr{:});
end
m.counts = max(0, m.counts);

%% update cell volume
mass.calcMass();
geom.calculateVolume();

%% apply media and stimulation conditions
m.counts = Condition.applyConditions(m.counts, m.setCounts, t.values); %media
s.values = Condition.applyConditions(s.values, s.setValues, t.values); %stimulation

%% Synchronize process states with simulation
for i = 1:length(this.processes)
    this.processes{i}.copyFromState();
end

%% If cell is dead, rerun initialize state
if abs(mr.growth - mr.meanInitialGrowthRate) / mr.meanInitialGrowthRate > mr.initialGrowthFilterWidth && ~isempty(this.seed)
    this.seed = this.randStream.randi([0 2^32-1], 1);
    this.initializeState();
end
