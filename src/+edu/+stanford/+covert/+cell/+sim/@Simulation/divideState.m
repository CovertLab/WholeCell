function [daughters, state_daughters] = divideState(this)
% Divide state of organism into 2 daughter cells.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/22/2012

import edu.stanford.covert.cell.sim.constant.Condition;
import edu.stanford.covert.util.CircularSparseMat;
import edu.stanford.covert.util.SparseMat;

c = this.compartment;
m = this.state('Metabolite');
r = this.state('Rna');
pm = this.state('ProteinMonomer');
pc = this.state('ProteinComplex');
ring = this.state('FtsZRing');
rnaPol = this.state('RNAPolymerase');
transcript = this.state('Transcript');
chr = this.state('Chromosome');
rib = this.state('Ribosome');
pol = this.state('Polypeptide');
trl = this.process('Translation');
geom = this.state('Geometry');

%% check if cell has divided
if chr.ploidy ~= 2 || ~chr.segregated || geom.pinchedDiameter ~= 0
    throw(MException('Simulation:divideState:notReplicated', ...
        'Cannot properly divide state because cell has not completely replicated.'));
end

%% initiate daughters
daughters = cell(2, 1);
for i = 1:2
    tmp = load('data/Simulation_fitted.mat');
    tmp.constructRandStream();
    
    daughters{i} = tmp.simulation;
    
    daughters{i}.allocateMemoryForState(1);
    
    daughters{i}.applyOptions(this.getOptions());
    daughters{i}.applyParameters(this.getParameters());
    
    daughters{i}.applyPerturbations();
end

%% time
for i = 1:2
    daughters{i}.state('Time').values = 0;
end

%% partition molecules
m_cnts = m.counts;
r_cnts = r.counts;
pm_cnts = pm.counts;
pc_cnts = pc.counts;

%chromosomes
d1_chr = daughters{1}.state('Chromosome');
d2_chr = daughters{2}.state('Chromosome');

props = setdiff(chr.stateNames, 'segregated');
for i = 1:numel(props)
    [posStrnd, val] = find(chr.(props{i}));
    tfs = posStrnd(:, 2) <= 2;
    
    posStrnd1 = posStrnd( tfs, :);
    posStrnd2 = posStrnd(~tfs, :);
    posStrnd2(:, 2) = posStrnd2(:, 2) - 2;
    val1 = val( tfs, :);
    val2 = val(~tfs, :);
    d1_chr.(props{i}) = CircularSparseMat(posStrnd1, val1, [chr.sequenceLen, chr.nCompartments], 1);
    d2_chr.(props{i}) = CircularSparseMat(posStrnd2, val2, [chr.sequenceLen, chr.nCompartments], 1);
end

d1_chr.segregated = false;
d2_chr.segregated = false;

%nascent transcripts, RNA polymerases
d1_rnaPol = daughters{1}.state('RNAPolymerase');
d2_rnaPol = daughters{2}.state('RNAPolymerase');
d1_transcript = daughters{1}.state('Transcript');
d2_transcript = daughters{2}.state('Transcript');

nFree = sum(rnaPol.states == rnaPol.freeValue);
nFree1 = rnaPol.randStream.binornd(nFree, 0.5);
nFree2 = nFree - nFree1;

tfs1 = rnaPol.states ~= rnaPol.notExistValue & rnaPol.positionStrands(:, 2) <= 2 & rnaPol.positionStrands(:, 2) >= 1;
tfs2 = rnaPol.states ~= rnaPol.notExistValue & rnaPol.positionStrands(:, 2)  > 2;
nPad1 = sum(tfs1) + nFree1;
nPad2 = sum(tfs2) + nFree2;
d1_rnaPol.states = [rnaPol.states(tfs1, :); repmat(rnaPol.freeValue, nFree1, 1); zeros(nPad1, 1)];
d2_rnaPol.states = [rnaPol.states(tfs2, :); repmat(rnaPol.freeValue, nFree2, 1); zeros(nPad2, 1)];
d1_rnaPol.positionStrands = [rnaPol.positionStrands(tfs1, 1) rnaPol.positionStrands(tfs1, 2);   zeros(nFree1 + nPad1, 2)];
d2_rnaPol.positionStrands = [rnaPol.positionStrands(tfs2, 1) rnaPol.positionStrands(tfs2, 2)-2; zeros(nFree2 + nPad2, 2)];
d1_transcript.boundTranscriptionUnits = [transcript.boundTranscriptionUnits(tfs1, :); zeros(nFree1 + nPad1, 1)];
d2_transcript.boundTranscriptionUnits = [transcript.boundTranscriptionUnits(tfs2, :); zeros(nFree2 + nPad2, 1)];
d1_transcript.boundTranscriptProgress = [transcript.boundTranscriptProgress(tfs1, :); zeros(nFree1 + nPad1, 1)];
d2_transcript.boundTranscriptProgress = [transcript.boundTranscriptProgress(tfs2, :); zeros(nFree2 + nPad2, 1)];
d1_transcript.boundTranscriptChromosome = [transcript.boundTranscriptChromosome(tfs1) == 1; zeros(nFree1 + nPad1, 1)];
d2_transcript.boundTranscriptChromosome = [transcript.boundTranscriptChromosome(tfs2) == 2; zeros(nFree2 + nPad2, 1)];

d1_rnaPol.transcriptionFactorBindingProbFoldChange = [rnaPol.transcriptionFactorBindingProbFoldChange(:, 1) zeros(size(rnaPol.transcriptionFactorBindingProbFoldChange, 1), 1)];
d2_rnaPol.transcriptionFactorBindingProbFoldChange = [zeros(size(rnaPol.transcriptionFactorBindingProbFoldChange, 1), 1) rnaPol.transcriptionFactorBindingProbFoldChange(:, 1)];

d1_rnaPol.supercoilingBindingProbFoldChange = [rnaPol.supercoilingBindingProbFoldChange(:, 1) zeros(size(rnaPol.supercoilingBindingProbFoldChange, 1), 1)];
d2_rnaPol.supercoilingBindingProbFoldChange = [zeros(size(rnaPol.supercoilingBindingProbFoldChange, 1), 1) rnaPol.supercoilingBindingProbFoldChange(:, 1)];

n = size(transcript.abortedTranscripts, 1);
k = transcript.randStream.binornd(n, 0.5);
tfs = false(n, 1);
tfs(pol.randStream.randomlySelectNRows((1:n)', k)) = true;
daughters{1}.state('Polypeptide').abortedPolypeptides = transcript.abortedTranscripts( tfs, :);
daughters{2}.state('Polypeptide').abortedPolypeptides = transcript.abortedTranscripts(~tfs, :);

%ribosomes, nascent polypeptides
d1_rib = daughters{1}.state('Ribosome');
d2_rib = daughters{2}.state('Ribosome');

n = size(rib.states, 1);
k = rib.randStream.binornd(n, 0.5);
tfs = false(n, 1);
tfs(rib.randStream.randomlySelectNRows((1:n)', k)) = true;

nPad = sum(tfs);
d1_rib.states = [rib.states( tfs, :); zeros(nPad, 1)];
d2_rib.states = [rib.states(~tfs, :); zeros(nPad, 1)];
d1_rib.boundMRNAs = [rib.boundMRNAs( tfs, :); zeros(nPad, 1)];
d2_rib.boundMRNAs = [rib.boundMRNAs(~tfs, :); zeros(nPad, 1)];
d1_rib.mRNAPositions = [rib.mRNAPositions( tfs, :); zeros(nPad, 1)];
d2_rib.mRNAPositions = [rib.mRNAPositions(~tfs, :); zeros(nPad, 1)];
d1_rib.tmRNAPositions = [rib.tmRNAPositions( tfs, :); zeros(nPad, 1)];
d2_rib.tmRNAPositions = [rib.tmRNAPositions(~tfs, :); zeros(nPad, 1)];

for i = 1:numel(daughters)
    d_pol = daughters{i}.state('Polypeptide');
    d_rib = daughters{i}.state('Ribosome');
    
    d_pol.boundMRNAs = d_rib.boundMRNAs;
    d_pol.nascentMonomerLengths = d_rib.mRNAPositions;
    d_pol.proteolysisTagLengths = d_rib.tmRNAPositions;
end

n = size(pol.abortedPolypeptides, 1);
k = pol.randStream.binornd(n, 0.5);
tfs = false(n, 1);
tfs(pol.randStream.randomlySelectNRows((1:n)', k)) = true;
daughters{1}.state('Polypeptide').abortedPolypeptides = pol.abortedPolypeptides( tfs, :);
daughters{2}.state('Polypeptide').abortedPolypeptides = pol.abortedPolypeptides(~tfs, :);

%FtsZ ring
pc_cnts(pc.boundIndexs(pc.getIndexs('MG_224_9MER_GTP')), c.cytosolIndexs) = ...
    + pc_cnts(pc.boundIndexs(pc.getIndexs('MG_224_9MER_GTP')), c.cytosolIndexs) ...
    - ring.numEdgesOneStraight ...
    - 2 * ring.numEdgesTwoStraight;
pc_cnts(pc.matureIndexs(pc.getIndexs('MG_224_9MER_GTP')), c.cytosolIndexs) = ...
    + pc_cnts(pc.matureIndexs(pc.getIndexs('MG_224_9MER_GTP')), c.cytosolIndexs) ...
    + ring.numEdgesOneStraight ...
    + 2 * ring.numEdgesTwoStraight;

pc_cnts(pc.boundIndexs(pc.getIndexs('MG_224_9MER_GDP')), c.cytosolIndexs) = ...
    + pc_cnts(pc.boundIndexs(pc.getIndexs('MG_224_9MER_GDP')), c.cytosolIndexs) ...
    - ring.numResidualBent ...
    - 2 * ring.numEdgesTwoBent;
pc_cnts(pc.matureIndexs(pc.getIndexs('MG_224_9MER_GDP')), c.cytosolIndexs) = ...
    + pc_cnts(pc.matureIndexs(pc.getIndexs('MG_224_9MER_GDP')), c.cytosolIndexs) ...
    + ring.numResidualBent ...
    + 2 * ring.numEdgesTwoBent;

for i = 1:numel(daughters)    
    d_ring = daughters{i}.state('FtsZRing');
    d_ring.numEdgesOneStraight = 0;
    d_ring.numEdgesTwoStraight = 0;
    d_ring.numEdgesTwoBent = 0;
    d_ring.numResidualBent = 0;
end

%metabolites
daughters{1}.state('Metabolite').counts = ...
    + m.randStream.binornd(m_cnts .* (m_cnts < 1e4), 0.5) ...
    + min(m_cnts, m.randStream.poissrnd(m_cnts .* (m_cnts >= 1e4) * 0.5));
daughters{2}.state('Metabolite').counts = m_cnts - daughters{1}.state('Metabolite').counts;

%RNA
d1_r = daughters{1}.state('Rna');
d2_r = daughters{2}.state('Rna');
d1_r.counts =  r.randStream.binornd(r_cnts, 0.5);
d2_r.counts =  r_cnts - d1_r.counts;

d1_r.counts(r.boundIndexs(trl.tmRNAGlobalIndex), c.cytosolIndexs) = ...
    d1_rib.nStalled;
d2_r.counts(r.boundIndexs(trl.tmRNAGlobalIndex), c.cytosolIndexs) = ...
    + r_cnts(r.boundIndexs(trl.tmRNAGlobalIndex), c.cytosolIndexs) ...
    - d1_r.counts(r.boundIndexs(trl.tmRNAGlobalIndex), c.cytosolIndexs);

%protein monomers
d1_pm = daughters{1}.state('ProteinMonomer');
d2_pm = daughters{2}.state('ProteinMonomer');
d1_pm.counts =  pm.randStream.binornd(pm_cnts, 0.5);
d2_pm.counts =  pm_cnts - d1_pm.counts;

[~, vals1] = find(d1_chr.monomerBoundSites);
[~, vals2] = find(d2_chr.monomerBoundSites);
d1_boundMon = histc(vals1, (1:numel(pm.boundIndexs))');
d2_boundMon = histc(vals2, (1:numel(pm.boundIndexs))');
idxs = setdiff((1:numel(pm.boundIndexs))', [
    pm.getIndexs(trl.enzymeWholeCellModelIDs(trl.enzymeIndexs_elongationFactors))
    ]); %#ok<NBRAK>
d1_pm.counts(pm.boundIndexs(idxs), c.cytosolIndexs) = d1_boundMon(idxs);
d2_pm.counts(pm.boundIndexs(idxs), c.cytosolIndexs) = d2_boundMon(idxs);

%protein complexes
d1_pc = daughters{1}.state('ProteinComplex');
d2_pc = daughters{2}.state('ProteinComplex');
d1_pc.counts = pc.randStream.binornd(pc_cnts, 0.5);
d2_pc.counts = pc_cnts - d1_pc.counts;

[~, vals1] = find(d1_chr.complexBoundSites);
[~, vals2] = find(d2_chr.complexBoundSites);
d1_boundCpx = histc(vals1, (1:numel(pc.boundIndexs))');
d2_boundCpx = histc(vals2, (1:numel(pc.boundIndexs))');
idxs = setdiff((1:numel(pc.boundIndexs))', [
    pc.getIndexs(trl.enzymeWholeCellModelIDs(trl.enzymeIndexs_elongationFactors))
    pc.getIndexs({'MG_224_9MER_GTP'; 'MG_224_9MER_GDP'})
    pc.ribosome70SIndexs
    ]);
d1_pc.counts(pc.boundIndexs(idxs), c.cytosolIndexs) = d1_boundCpx(idxs);
d2_pc.counts(pc.boundIndexs(idxs), c.cytosolIndexs) = d2_boundCpx(idxs);

d1_pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs(1)), c.cytosolIndexs) = ...
    min(d1_rnaPol.nFree, pc_cnts(pc.matureIndexs(pc.rnaPolymeraseIndexs(1)), c.cytosolIndexs));
d1_pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs(2)), c.cytosolIndexs) = ...
    + d1_rnaPol.nFree ...
    - d1_pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs(1)), c.cytosolIndexs);
d2_pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), c.cytosolIndexs) =  ...
    + pc_cnts(pc.matureIndexs(pc.rnaPolymeraseIndexs), c.cytosolIndexs) ...
    - d1_pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), c.cytosolIndexs);

d1_pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs) = d1_rib.nActive + d1_rib.nStalled;
d2_pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), c.cytosolIndexs) = d2_rib.nActive + d2_rib.nStalled;

%% Split process requirements, allocations, usages
mm = this.state('Metabolite');
dm1 = daughters{1}.state('Metabolite');
dm2 = daughters{2}.state('Metabolite');

mm.processRequirements = full(mm.processRequirements);
mm.processAllocations = full(mm.processAllocations);
mm.processUsages = full(mm.processUsages);

dm1.processRequirements = mm.randStream.stochasticRound(mm.processRequirements / 2);
dm1.processAllocations = mm.randStream.stochasticRound(mm.processAllocations / 2);
dm1.processUsages = mm.randStream.stochasticRound(mm.processUsages / 2);

dm2.processRequirements = mm.processRequirements - dm1.processRequirements;
dm2.processAllocations = mm.processAllocations - dm1.processAllocations;
dm2.processUsages = mm.processUsages - dm1.processUsages;

for i = 1:numel(daughters)
    dm = daughters{i}.state('Metabolite');
    dm.processRequirements = SparseMat(dm.processRequirements);
    dm.processAllocations = SparseMat(dm.processAllocations);
    dm.processUsages = SparseMat(dm.processUsages);
end

%% Set external conditions
for i = 1:numel(daughters)
    t = daughters{i}.state('Time');
    
    %media
    m = daughters{i}.state('Metabolite');
    m.counts = Condition.applyConditions(m.counts, m.setCounts, t.values);
    
    %stimuli
    s = daughters{i}.state('Stimulus');
    s.values = Condition.applyConditions(s.values, s.setValues, t.values);
end

%% recalculate mass, geometry, growth, metabolic fluxes
for i = 1:numel(daughters)
    daughters{i}.state('Mass').initialize();
    daughters{i}.state('Geometry').initialize();
    
    daughters{i}.process('HostInteraction').copyFromState();
    daughters{i}.process('HostInteraction').initializeState();
    
    daughters{i}.process('Metabolism').copyFromState();
    daughters{i}.process('Metabolism').initializeState();
end

%% Extract state
state_daughters = struct();
for i = 1:numel(daughters)
    for j = 1:numel(daughters{i}.states)
        s = daughters{i}.states{j};
        sid = s.wholeCellModelID(7:end);
        state_daughters(i).(sid) = struct();
        
        for k = 1:numel(s.stateNames)
            state_daughters(i).(sid).(s.stateNames{k}) = ...
                s.(s.stateNames{k});
        end
        
        for k = 1:numel(s.dependentStateNames)
            state_daughters(i).(sid).(s.dependentStateNames{k}) = ...
                s.(s.dependentStateNames{k});
        end
    end
end
