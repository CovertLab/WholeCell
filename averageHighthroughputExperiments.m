%averageHighthroughputExperiments
% Averages multiple in silico high-throughput experiments.
%
% Inputs (apply using key, value pairs)
% - inPathPattern [.mat file path pattern]: file path to .mat files
%   created by the HighthroughputExperimentsLogger.
%
% Output
% - avgVals [struct]: Struct containing average value of in silico
%   experimental data
% - labels [struct]: Struct containing row and column labels for avgVals
% - If outPath is set, saves average simulated in silico experimental data
%   to .mat file
%
% Example:
%   sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
%   parameterVals = sim.getAllParameters(); %get default parameters
%   parameterVals.lengthSec = 1000;         %override defaults
%   for i = 1:10
%     simulateHighthroughputExperiments(...
%           'seed', seed, ...
%           'parameterVals', parameterVals, ...
%           'outPath', sprintf('output/dream-sim-%d.mat', i) ...
%           );
%   end
%   averageHighthroughputExperiments(...
%       'inPathPattern', 'output/dream-sim-*.mat', ...
%       'outPath', 'output/dream-sim-avg.mat' ...
%       );
%
% See also:
% - edu.stanford.covert.cell.sim.Simulation
% - edu.stanford.covert.cell.sim.util.HighthroughputExperimentsLogger
% - simulateHighthroughputExperiments
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
function [avgVals, labels] = averageHighthroughputExperiments(varargin)

%import classes
import edu.stanford.covert.cell.sim.util.HighthroughputExperimentsLogger;

%process arguments
ip = inputParser();

ip.addParamValue('inPathPattern', '', @(x) ischar(x));
ip.addParamValue('outPath', '', @(x) ischar(x));
ip.addParamValue('verbosity', 1);

ip.parse(varargin{:});

inPathPattern = ip.Results.inPathPattern;
outPath       = ip.Results.outPath;
verbosity     = ip.Results.verbosity;

if ischar(verbosity)
    verbosity = str2double(verbosity);
end
validateattributes(verbosity, {'numeric'}, {'integer'});

%average experiments
[avgVals, labels] = HighthroughputExperimentsLogger.average(inPathPattern, verbosity);
if ~isempty(outPath)
    tmp = avgVals;
    tmp.labels = labels;
    save(outPath, '-struct', 'tmp');
    clear tmp;
end