%calcParametersAndPredictionErrors
% 
% Input
% - parameterVals [struct]: struct containing desired values of simulation
%   parameters. Initialize using the getAllParameters method of the
%   simulation class to get default parameter values. Overwrite struct
%   elements to set parameter values.
% - parameterValsPath [.mat/.xml file path]: .mat file path for stored
%   struct of parameter values or .xml file describing parameter values.
%   Use http://wholecell.stanford.edu/simulation/runSimulations.php to
%   generate XML file.
% - avgVals [struct]: Struct containing average value of in silico
%   experimental data
% - avgValsPath [.mat file path]: .mat file path to struct containing
%   average value of in silico experimental data
% - refParameterVals: same as parameterVals, but for reference parameter
%   values
% - refParameterValsPath: same as parameterValsPath, but for reference parameter
%   values
% - refAvgVals: same as avgVals, but for reference parameter
%   values
% - refAvgValsPath: same as avgValsPath:, but for reference parameter
%   values
%
% Output
% - dists [struct]: struct with two fields (parameter, prediction)
%   containing the cacluated parameter and prediction error from comparison
%   to the reference 
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
function dists = calcParametersAndPredictionErrors(varargin)

%% process arguments
ip = inputParser();

ip.addParamValue('parameterVals', [], @(x) isstruct(x));
ip.addParamValue('parameterValsPath', '', @(x) exist(x, 'file'));
ip.addParamValue('avgVals', [], @(x) isstruct(x));
ip.addParamValue('avgValsPath', '', @(x) ischar(x));
ip.addParamValue('refParameterVals', [], @(x) isstruct(x));
ip.addParamValue('refParameterValsPath', '', @(x) exist(x, 'file'));
ip.addParamValue('refAvgVals', [], @(x) isstruct(x));
ip.addParamValue('refAvgValsPath', '', @(x) exist(x, 'file'));

ip.parse(varargin{:});

parameterVals        = ip.Results.parameterVals;
parameterValsPath    = ip.Results.parameterValsPath;
avgVals              = ip.Results.avgVals;
avgValsPath          = ip.Results.avgValsPath;
refParameterVals     = ip.Results.refParameterVals;
refParameterValsPath = ip.Results.refParameterValsPath;
refAvgVals           = ip.Results.refAvgVals;
refAvgValsPath       = ip.Results.refAvgValsPath;

%load parameter values
refParameterVals = loadParameterVals(refParameterVals, refParameterValsPath);
parameterVals = loadParameterVals(parameterVals, parameterValsPath);

%load average values
avgVals = loadAvgVals(avgVals, avgValsPath);
refAvgVals = loadAvgVals(refAvgVals, refAvgValsPath);

%% distances
dists = struct('parameter', [], 'prediction', []);

%parameters
dists.parameter = calcAvgSumSquaredLogRatio(...
    getParameterVector(parameterVals), ...
    getParameterVector(refParameterVals));

%prediction
dists.prediction = calcAvgNormSquaredDiff(...
    getPredictionVector(avgVals), ...
    getPredictionVector(refAvgVals));

%Get struct of parameter values
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

%Get struct of averages of high-throughput in silico experiments
function avgVals = loadAvgVals(avgVals, avgValsPath)
if ~isempty(avgVals)
elseif ~isempty(avgValsPath)
    avgVals = load(avgValsPath);
end

%Create feature vector from full set of parameters
function paramVec = getParameterVector(paramStruct)
sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
rna = sim.state('Rna');
met = sim.process('Metabolism');
trn = sim.process('Transcription');

halfLives = rna.halfLives;
enzBounds = met.enzymeBounds;
tuBindProbs = trn.transcriptionUnitBindingProbabilities;

if isfield(paramStruct, 'states') && isfield(paramStruct.states, 'Rna') && isfield(paramStruct.states.Rna, 'halfLives')
    halfLives = paramStruct.states.Rna.halfLives;
end
if isfield(paramStruct, 'processes') && isfield(paramStruct.processes, 'Metabolism') && isfield(paramStruct.processes.Metabolism, 'enzymeBounds')
    enzBounds = paramStruct.processes.Metabolism.enzymeBounds;
end
if isfield(paramStruct, 'processes') && isfield(paramStruct.processes, 'Transcription') && isfield(paramStruct.processes.Transcription, 'transcriptionUnitBindingProbabilities')
    tuBindProbs = paramStruct.processes.Transcription.transcriptionUnitBindingProbabilities;
end

paramVec = [
    halfLives(rna.matureIndexs)
    enzBounds(:, 1)
    enzBounds(:, 2)
    tuBindProbs
    ];

%Create feature vector from predicted high-throughput data
function predVec = getPredictionVector(predStruct)
sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
time = sim.state('Time');
cellCycleLength = time.cellCycleLength;

predVec = struct('mean', [], 'std', []);
predVec.mean = [
    mean(nanmean(predStruct.growth, 2), 1)
    mean(nanmean(predStruct.mass, 2), 1)
    mean(nanmean(predStruct.volume, 2), 1)
    
    mean(min(cellCycleLength, predStruct.repInitTime))
    mean(min(cellCycleLength, predStruct.repTermTime))
    mean(min(cellCycleLength, predStruct.cellCycleLen))
    
    predStruct.metConcs.mean
    predStruct.dnaSeq.mean
    predStruct.rnaSeq.mean
    predStruct.chipSeq.mean(:)
    predStruct.rnaArray.mean
    predStruct.protArray.mean
    predStruct.rxnFluxes.mean
    ];
predVec.std = [
    sqrt(mean(nanvar(predStruct.growth, 1, 2)))
    sqrt(mean(nanvar(predStruct.mass, 1, 2)))
    sqrt(mean(nanvar(predStruct.volume, 1, 2)))
    
    std(min(cellCycleLength, predStruct.repInitTime))
    std(min(cellCycleLength, predStruct.repTermTime))
    std(min(cellCycleLength, predStruct.cellCycleLen))
    
    predStruct.metConcs.std
    predStruct.dnaSeq.std
    predStruct.rnaSeq.std
    predStruct.chipSeq.std(:)
    predStruct.rnaArray.std
    predStruct.protArray.std
    predStruct.rxnFluxes.std
    ];

%Calculate <(log(test/ref))^2>
function dist = calcAvgSumSquaredLogRatio(test, ref)
test = max(-realmax, min(realmax, test));
ref  = max(-realmax, min(realmax, ref));

dist = mean(log(test ./ ref) .^ 2);

dist = full(dist);

%Calculate <((mean(ref) - mean(test))^2) / var(ref)>
function dist = calcAvgNormSquaredDiff(test, ref)
test.mean = max(-realmax, min(realmax, test.mean));
ref.mean  = max(-realmax, min(realmax, ref.mean));
test.std  = max(-realmax, min(realmax, test.std));
ref.std   = max(-realmax, min(realmax, ref.std));

dist = mean(((ref.mean - test.mean) .^ 2) ./ (ref.std .^ 2));

dist = full(dist);