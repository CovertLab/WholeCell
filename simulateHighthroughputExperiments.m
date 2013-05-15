%simulateHighthroughputExperiments
% Runs whole-cell simulation using specified parameter values and returns
% simulated high-throughput experimental data (see
% HighthroughputExperimentsLogger for more information).
%
% Inputs (apply using key, value pairs)
% - seed [integer]: random number generator seed
% - parameterVals [struct]: struct containing desired values of simulation
%   parameters. Initialize using the getAllParameters method of the
%   simulation class to get default parameter values. Overwrite struct
%   elements to set parameter values.
% - parameterValsPath [.mat/.xml file path]: .mat file path for stored
%   struct of parameter values or .xml file describing parameter values.
%   Use http://wholecell.stanford.edu/simulation/runSimulations.php to
%   generate XML file.
% - initialConditions [.
% - simPath [.mat file path]: Desired file path for simulated in silico
%   experimental data (see HighthroughputExperimentsLogger for information
%   about simulated data)
% - verbosity [integer]: Desired verbosity level. Zero supresses output.
%   Higher value prints more output.
%
% Output
% - If simPath is set, saves simulated in silico experimental data (see
%   HighthroughputExperimentsLogger for information about simulated data)
%   to path specified by simPath
%
% Examples:
%   %Run simulation using struct of parameter values
%   sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
%   parameterVals = sim.getAllParameters(); %get default parameters
%   parameterVals.lengthSec = 1000;         %override defaults
%   simulateHighthroughputExperiments(...
%       'seed', 1, ...
%       'parameterVals', parameterVals, ...
%       'simPath', 'sim1.mat' ...
%       );
%
%   %Run simulation using mat file of parameter values
%   sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
%   parameterVals = sim.getAllParameters(); %get default parameters
%   parameterVals.lengthSec = 1000;         %override defaults
%   parameterValsPath = 'sim1-parameters.mat';
%   save(parameterValsPath, '-struct', 'parameterVals');
%   simulateHighthroughputExperiments(...
%       'seed', 1, ...
%       'parameterValsPath', parameterValsPath, ...
%       'simPath', 'sim1.mat' ...
%       );
%
% See also:
% - edu.stanford.covert.cell.sim.Simulation
% - edu.stanford.covert.cell.sim.util.HighthroughputExperimentsLogger
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
function simulateHighthroughputExperiments(varargin)

%% import classes
import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
import edu.stanford.covert.cell.sim.util.ConditionSet;
import edu.stanford.covert.cell.sim.util.HighthroughputExperimentsLogger;
import edu.stanford.covert.cell.sim.util.SummaryLogger;
import edu.stanford.covert.util.StructUtil;

%% process arguments
ip = inputParser();

ip.addParamValue('seed', []);
ip.addParamValue('geneticKnockouts', [], @(x) ischar(x) || iscell(x));
ip.addParamValue('parameterVals', [], @(x) isstruct(x));
ip.addParamValue('parameterValsPath', '', @(x) exist(x, 'file'));
ip.addParamValue('simPath', '', @(x) ischar(x));
ip.addParamValue('initialConditions', struct(), @(x) isstruct(x));
ip.addParamValue('initialConditionsPath', '', @(x) exist(x, 'file'));
ip.addParamValue('verbosity', 1);

ip.parse(varargin{:});

seed                  = ip.Results.seed;
geneticKnockouts      = ip.Results.geneticKnockouts;
parameterVals         = ip.Results.parameterVals;
parameterValsPath     = ip.Results.parameterValsPath;
initialConditions     = ip.Results.initialConditions;
initialConditionsPath = ip.Results.initialConditionsPath;
simPath               = ip.Results.simPath;
verbosity             = ip.Results.verbosity;

if ischar(seed)
    seed = str2double(seed);
    validateattributes(seed, {'numeric'}, {'integer'});
end

if ischar(geneticKnockouts)
    geneticKnockouts = {geneticKnockouts};
end

if numel(fields(initialConditions)) == 0 && ~isempty(initialConditionsPath)
    initialConditions = load(initialConditionsPath);
end

if ischar(verbosity)
    verbosity = str2double(verbosity);
end
validateattributes(verbosity, {'numeric'}, {'integer'});

%% load simulation
sim = CachedSimulationObjectUtil.load();

%% set parameter values
%parameters
if ~isempty(parameterVals)
    sim.applyAllParameters(parameterVals);
elseif ~isempty(parameterValsPath)
    [~, ~, ext] = fileparts(parameterValsPath);
    switch ext
        case '.mat'
            parameterVals = load(parameterValsPath);
        case '.xml'
            tmp = ConditionSet.parseConditionSet(sim, parameterValsPath);
            parameterVals = struct();
            parameterVals = StructUtil.catstruct(parameterVals, tmp.options);
            parameterVals = StructUtil.catstruct(parameterVals, tmp.perturbations);
            parameterVals = StructUtil.catstruct(parameterVals, tmp.parameters);
            parameterVals = StructUtil.catstruct(parameterVals, tmp.fittedConstants);
            parameterVals = StructUtil.catstruct(parameterVals, tmp.fixedConstants);
        otherwise
            throw(MException('simulateHighthroughputExperiments:invalidFileFormat', 'Invalid input file format "%s"', ext));
    end
    sim.applyAllParameters(parameterVals);
end

%verbosity
sim.applyOptions('verbosity', verbosity);

%knockout perturbations
if iscell(geneticKnockouts)
    sim.applyOptions('geneticKnockouts', geneticKnockouts);
end

%seed random number generator
if ~isempty(seed)
    sim.applyOptions('seed', seed);
end

%% setup loggers
loggers = {SummaryLogger(1, verbosity)};
if ~isempty(simPath)
    loggers = [loggers; {HighthroughputExperimentsLogger(simPath)}];
end

%% run
sim.run(initialConditions, loggers);