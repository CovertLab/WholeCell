%Runs the simulation, and optionally logs results. varargin optionally
%contains two input arguments:
%- adf
%- single instance of SimulationLogger or cell array of a SimulationLogger instances
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 3/24/2011
function [this, loggers] = run(this, varargin)

%process options
loggers = {};
ic = struct();
for i = 1:numel(varargin)
    if isa(varargin{i}, 'edu.stanford.covert.cell.sim.util.Logger') || ...
            (iscell(varargin{i}) && isa(varargin{i}{1}, 'edu.stanford.covert.cell.sim.util.Logger'))
        loggers = varargin{i};
        if ~iscell(loggers)
            loggers = {loggers};
        end
    else
        ic = varargin{i};
    end
end

%references
g = this.state('Geometry');
met = this.state('Metabolite');

%allocate memory
this.allocateMemoryForState(1);

%initialize state
this.initializeState();
if isstruct(ic) && isfield(ic, 'states') && numel(fieldnames(ic.states)) > 0
    %override default initial conditions
    fields = fieldnames(ic.states);
    for i = 1:numel(fields)
        s = this.state(regexprep(fields{i}, '^State_', ''));
        
        subfields = fieldnames(ic.states.(fields{i}));
        for j = 1:numel(subfields)
            tmp = ic.states.(fields{i}).(subfields{j});
            if isnumeric(tmp)
                s.(subfields{j})(~isnan(tmp)) = tmp(~isnan(tmp));
            else
                s.(subfields{j}) = tmp;
            end
        end
    end
    
    %synchronize calculated state
    this.state('Time').values = 0;
    this.state('Mass').initialize();
    this.state('Geometry').initialize();
   
    %synchronize processes
    for i = 1:numel(this.processes)
        proc = this.processes{i};
        proc.copyFromState();
    end
end

%apply perturbations
this.applyPerturbations();

%evolve state
for j = 1:numel(loggers)
    loggers{j}.initialize(this);
end

try    
    for i = 1:this.getNumSteps
        [~, requirements, allocations, usages] = this.evolveState();
        met.processRequirements = edu.stanford.covert.util.SparseMat(requirements);
        met.processAllocations = edu.stanford.covert.util.SparseMat(allocations);
        met.processUsages = edu.stanford.covert.util.SparseMat(usages);
        
        for j = 1:numel(loggers)
            loggers{j}.append(this);
        end
        if ~isempty(g) && g.pinched
            break;
        end
    end
catch exception
    for j = 1:numel(loggers)
        loggers{j}.finalize(this);
    end
    exception.rethrow();
end

for j = 1:numel(loggers)
    loggers{j}.finalize(this);
end