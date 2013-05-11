function [this, requirements, allocations, usages] = evolveState(this)
% Evolves state of organism:
% 1. Increments time
% 2. Evaluates stimuli
% 3. Evaluates processes
% 4. Applies media conditions
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanfod.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/15/2010

import edu.stanford.covert.cell.sim.constant.Condition;

time = this.state_time;
mets = this.state_metabolite;
stim = this.state_stimulus;

time.values = time.values + this.stepSizeSec;

%% evaluate/apply stimuli
stim.values = Condition.applyConditions(stim.values, stim.setValues, time.values);

%% estimate metabolic requirements of processes
processes = this.processes;
nProcesses = length(processes);
requirements = zeros([numel(mets.counts) nProcesses]);
for i = 1:nProcesses
    mod = processes{i};
    mod.copyFromState();
    r = mod.calcResourceRequirements_Current();
    requirements(mod.substrateMetaboliteGlobalCompartmentIndexs, i) = ...
        reshape(r(mod.substrateMetaboliteLocalIndexs, :), [], 1);
end
requirements = max(0, requirements);
tmp = mets.counts(:) ./ max(1, sum(requirements, 2));
allocations = max(0, fix(requirements .* tmp(:, ones(nProcesses, 1))));
if nargout >= 4
    usages = zeros(size(allocations));
end

%% run cell processes (a.k.a. processes):
% - update them with latest state
% - allocate metabolites fairly
% - run the cell processes
% - update simulation with latest state

%determine order of evaluation with constraint that tRNA aminoacylation
%always occurs in same order with respect to translation
while true
    processEvalOrderIndexs = this.randStream.randperm(nProcesses);
    idx1 = find(processEvalOrderIndexs == this.processIndex_tRNAAminoacylation, 1);
    idx2 = find(processEvalOrderIndexs == this.processIndex_translation, 1);
    if isempty(idx1) || isempty(idx2) || idx1 < idx2
        break;
    end
end

%simulate processes
for i = 1:nProcesses
    mod = processes{processEvalOrderIndexs(i)};
    
    allocation = reshape(allocations(mod.substrateMetaboliteGlobalCompartmentIndexs, processEvalOrderIndexs(i)), ...
        size(mod.substrateMetaboliteGlobalCompartmentIndexs));
    counts = mets.counts(mod.substrateMetaboliteGlobalCompartmentIndexs);
    
    mod.simulationStateSideEffects = [];
    mod.copyFromState();
    mod.substrates(mod.substrateMetaboliteLocalIndexs, :) = allocation;
    mod.evolveState();
    mod.copyToState();
    mets.counts(mod.substrateMetaboliteGlobalCompartmentIndexs) = counts + ...
        mod.substrates(mod.substrateMetaboliteLocalIndexs, :) - allocation;
    if nargout >= 4
        usages(mod.substrateMetaboliteGlobalCompartmentIndexs, processEvalOrderIndexs(i)) = ...
            reshape(allocation - mod.substrates(mod.substrateMetaboliteLocalIndexs, :), [], 1);
    end
    if ~isempty(mod.simulationStateSideEffects)
        mod.simulationStateSideEffects.updateSimulationState(this);
    end
end

%% apply media conditions
mets.counts = Condition.applyConditions(mets.counts, mets.setCounts, time.values);
