%averageHighthroughputExperimentsAndCalcErrors
%
% Inputs (apply using key, value pairs)
% - parameterVals [struct]: struct containing desired values of simulation
%   parameters. Initialize using the getAllParameters method of the
%   simulation class to get default parameter values. Overwrite struct
%   elements to set parameter values.
% - parameterValsPath [.mat/.xml file path]: .mat file path for stored
%   struct of parameter values or .xml file describing parameter values.
%   Use http://wholecell.stanford.edu/simulation/runSimulations.php to
%   generate XML file.
% - simPathPattern [.mat file path pattern]: file path to .mat files
%   created by the HighthroughputExperimentsLogger.
% - avgValsPath [.mat file path]: Location where average values should be
%   saved in .mat format
% - refParameterVals: same as parameterVals, but for reference parameter
%   values
% - refParameterValsPath: same as parameterValsPath, but for reference parameter
%   values
% - refAvgVals: same as avgVals, but for reference parameter
%   values
% - refAvgValsPath: same as avgValsPath:, but for reference parameter
%   values
% - verbosity [integer]: Desired verbosity level. Zero supresses output.
%   Higher value prints more output.
%
% Output
% - dists [struct]: struct with two fields (parameter, prediction)
%   containing the cacluated parameter and prediction error from comparison
%   to the reference 
% - avgVals [struct]: Struct containing average value of in silico
%   experimental data
% - avgValLabels [struct]: Struct containing row and column labels for avgVals
% - If avgValsPath is set, saves average simulated in silico experimental data
%   to .mat file
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
function [dists, avgVals, avgValLabels] = averageHighthroughputExperimentsAndCalcErrors(varargin)

%% parse inputs
ip = inputParser();

ip.addParamValue('parameterVals', [], @(x) isstruct(x));
ip.addParamValue('parameterValsPath', '', @(x) exist(x, 'file'));
ip.addParamValue('simPathPattern', '', @(x) ischar(x));
ip.addParamValue('avgValsPath', '', @(x) ischar(x));
ip.addParamValue('distsPath', '', @(x) ischar(x));
ip.addParamValue('refParameterVals', [], @(x) isstruct(x));
ip.addParamValue('refParameterValsPath', '', @(x) exist(x, 'file'));
ip.addParamValue('refAvgVals', [], @(x) isstruct(x));
ip.addParamValue('refAvgValsPath', '', @(x) exist(x, 'file'));
ip.addParamValue('verbosity', 1);

ip.parse(varargin{:});

parameterVals        = ip.Results.parameterVals;
parameterValsPath    = ip.Results.parameterValsPath;
simPathPattern       = ip.Results.simPathPattern;
avgValsPath          = ip.Results.avgValsPath;
distsPath            = ip.Results.distsPath;
refParameterVals     = ip.Results.refParameterVals;
refParameterValsPath = ip.Results.refParameterValsPath;
refAvgVals           = ip.Results.refAvgVals;
refAvgValsPath       = ip.Results.refAvgValsPath;
verbosity            = ip.Results.verbosity;

%load parameter values
refParameterVals = loadParameterVals(refParameterVals, refParameterValsPath);
parameterVals = loadParameterVals(parameterVals, parameterValsPath);

%load average values
refAvgVals = loadAvgVals(refAvgVals, refAvgValsPath);

%parse verbosity
if ischar(verbosity)
    verbosity = str2double(verbosity);
end
validateattributes(verbosity, {'numeric'}, {'integer'});


%% average simulation
[avgVals, avgValLabels] = averageHighthroughputExperiments(...
    'simPathPattern', simPathPattern, ...
    'avgValsPath', avgValsPath, ...
    'verbosity', verbosity ...
    );

%% calculate parameter and prediction errors
dists = calcParametersAndPredictionErrors(...    
    'parameterVals', parameterVals, ...
    'avgVals', avgVals, ...
    'refParameterVals', refParameterVals, ...
    'refAvgVals', refAvgVals, ...
    'distsPath', distsPath ...
    );

%% send data to Synapse (TODO)
%- requires user id
%- simulation id (timestamp?)

function parameterVals = loadParameterVals(parameterVals, parameterValsPath)
if ~isempty(parameterVals)
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
end

function avgVals = loadAvgVals(avgVals, avgValsPath)
if ~isempty(avgVals)
elseif ~isempty(avgValsPath)
    avgVals = load(avgValsPath);
end
