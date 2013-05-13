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
% - parameterValsPath [.mat file path]: .mat file path for stored struct of
%   parameter values
% - outPath [.mat file path]: Desired file path for simulated in silico
%   experimental data (see HighthroughputExperimentsLogger for information
%   about simulated data)
%
% Output
% - If outPath is set, saves simulated in silico experimental data (see
%   HighthroughputExperimentsLogger for information about simulated data)
%   to path specified by outPath
%
% Examples:
%   %Run simulation using struct of parameter values
%   sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
%   parameterVals = sim.getAllParameters(); %get default parameters
%   parameterVals.lengthSec = 1000;         %override defaults
%   simulateHighthroughputExperiments(...
%       'seed', 1, ...
%       'parameterVals', parameterVals, ...
%       'outPath', 'sim1.mat' ...
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
%       'outPath', 'sim1.mat' ...
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

%import classes
import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
import edu.stanford.covert.cell.sim.util.HighthroughputExperimentsLogger;
import edu.stanford.covert.cell.sim.util.SummaryLogger;

%process arguments
ip = inputParser();

ip.addParamValue('seed', 0);
ip.addParamValue('parameterVals', [], @(x) isstruct(x));
ip.addParamValue('parameterValsPath', '', @(x) exist(x, 'file'));
ip.addParamValue('outPath', '', @(x) ischar(x));

ip.parse(varargin{:});

seed              = ip.Results.seed;
parameterVals     = ip.Results.parameterVals;
parameterValsPath = ip.Results.parameterValsPath;
outPath           = ip.Results.outPath;

if ischar(seed)
    seed = str2double(seed);
end
validateattributes(seed, {'numeric'}, {'integer'});

%load simulation
sim = CachedSimulationObjectUtil.load();

%seed random number generator
sim.applyOptions('seed', seed);

%set parameter values
if ~isempty(parameterVals)
    sim.applyAllParameters(parameterVals);
elseif ~isempty(parameterValsPath)
    parameterVals = load(parameterValsPath);
    sim.applyAllParameters(parameterVals);
end

%setup loggers
loggers = {SummaryLogger(1, 1)};
if ~isempty(outPath)
    loggers = [loggers; {HighthroughputExperimentsLogger(outPath)}];
end

%run
sim.run(loggers);