%DreamScoring
% Functions for calculating DREAM distances and scoring
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
classdef DreamScoring
    methods (Static = true)
        %calcScores
        % Calculates score and rank of each submission by the procedure established
        % in the DREAM 7 parameter estimation challenge (see
        % http://www.the-dream-project.org/sites/the-dream-project.org/files/documents/DREAM6/DREAM7/Descriptions/dream7c1_scoring.pdf)
        %
        % Input:
        % - dists [struct array]: struct array of distances calculated by
        %   calcParametersAndPredictionError
        %
        % Output
        % - Scores [double array]: score of each submission
        % - Ranks [int array]: rank of each submission
        function [submissionList, distances, pValues, scores, ranks] = calcScores(varargin)
            import edu.stanford.covert.cell.sim.util.DreamScoring;
            
            %% parse inputs
            ip = inputParser();
            
            ip.addParamValue('refParameterValsPath', '', @(x) exist(x, 'file'));
            ip.addParamValue('refAvgValsPath', '', @(x) exist(x, 'file'));
            ip.addParamValue('submissionsPathBase', '', @(x) exist(x, 'dir'));
            ip.addParamValue('nullValsDistribSize', 10, @(x) isnumeric(x) && x == ceil(x));
            ip.addParamValue('nullDistanceDistribSize', 1e4, @(x) isnumeric(x) && x == ceil(x));
            
            ip.parse(varargin{:});
            
            refParameterValsPath    = ip.Results.refParameterValsPath;
            refAvgValsPath          = ip.Results.refAvgValsPath;
            submissionsPathBase     = ip.Results.submissionsPathBase;
            nullValsDistribSize     = ip.Results.nullValsDistribSize;
            nullDistanceDistribSize = ip.Results.nullDistanceDistribSize;
            
            %% list of submissions
            submissionDistances = dir(fullfile(submissionsPathBase, '*.distances.*'));
            submissionParameters = dir(fullfile(submissionsPathBase, '*.parameters.*'));
            submissionPredictions = dir(fullfile(submissionsPathBase, '*.predictions.*'));
            
            nSubmissions = numel(submissionDistances);
            nullValsDistribSize = min(nullValsDistribSize, nSubmissions);
            
            submissionList = repmat(struct('user', [], 'parameterSetTimestamp', []), nSubmissions, 1);
            for i = 1:nSubmissions
                tmp = strsplit('_', submissionDistances(i).name(1:end-14));
                submissionList(i).user = tmp{1};
                submissionList(i).parameterSetTimestamp = tmp{2};
            end
            
            %% build table of distances
            pValues = repmat(struct('parameter', [], 'prediction', []), nSubmissions, 1);
            distances = repmat(struct('parameter', [], 'prediction', []), nSubmissions, 1);
            for i = 1:nSubmissions
                distances(i) = load(fullfile(submissionsPathBase, submissionDistances(i).name));
            end
            
            %% parameters
            refParameterVals = DreamScoring.getParameterVector(DreamScoring.loadParameterVals([], refParameterValsPath));
            nParameters = numel(refParameterVals);
            
            [~, bestDistanceIdxs] = sort([distances.parameter]);
            
            nullParameterVals = zeros(nParameters, nullValsDistribSize);
            for i = 1:nullValsDistribSize
                nullParameterVals(:, i) = DreamScoring.getParameterVector(DreamScoring.loadParameterVals([], fullfile(submissionsPathBase, submissionParameters(bestDistanceIdxs(i)).name)));
            end
            
            nullDistances = zeros(1, nullDistanceDistribSize);
            for i = 1:nullDistanceDistribSize
                nullDistances(i) = DreamScoring.calcAvgSumSquaredLogRatio(...
                    refParameterVals, ...
                    nullParameterVals(sub2ind(size(nullParameterVals), (1:nParameters)', randi(nullValsDistribSize, [nParameters 1]))));
            end
            
            tmp = DreamScoring.calcPVal([distances.parameter]', nullDistances');
            for i = 1:nSubmissions
                pValues(i).parameter = tmp(i);
            end
            
            %% predictions
            refAvgVals = DreamScoring.getPredictionVector(DreamScoring.loadAvgVals([], refAvgValsPath));
            nAvgVals = numel(refAvgVals.mean);
            
            [~, bestDistanceIdxs] = sort([distances.prediction]);
            
            nullAvgVals = struct('mean', zeros(nAvgVals, nullValsDistribSize), 'std', zeros(nAvgVals, nullValsDistribSize));
            for i = 1:nullValsDistribSize
                tmp = DreamScoring.getPredictionVector(DreamScoring.loadAvgVals([], fullfile(submissionsPathBase, submissionPredictions(bestDistanceIdxs(i)).name)));
                nullAvgVals.mean(:, i) = tmp.mean;
                nullAvgVals.std(:, i) = tmp.std;
            end
            
            nullDistances = zeros(1, nullDistanceDistribSize);
            for i = 1:nullDistanceDistribSize
                tmp = sub2ind(size(nullAvgVals.mean), (1:nAvgVals)', randi(nullValsDistribSize, [nAvgVals 1]));
                
                nullDistances(i) = DreamScoring.calcAvgNormSquaredDiff(...
                    refAvgVals, ...
                    struct('mean', nullAvgVals.mean(tmp), 'std', nullAvgVals.std(tmp)) ...
                    );
            end
            
            tmp = DreamScoring.calcPVal([distances.prediction]', nullDistances');
            for i = 1:nSubmissions
                pValues(i).prediction = tmp(i);
            end
            
            %% Calculate scores
            scores = -log([pValues.parameter]' .* [pValues.prediction]');
            
            %% calculate ranks
            [~, idxs] = sort(scores);
            ranks = zeros(size(scores));
            ranks(idxs) = numel(scores):-1:1;
        end
        
        %calcErrors
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
        function dists = calcErrors(varargin)
            import edu.stanford.covert.cell.sim.util.DreamScoring;
            
            %% process arguments
            ip = inputParser();
            
            ip.addParamValue('parameterVals', [], @(x) isstruct(x));
            ip.addParamValue('parameterValsPath', '', @(x) exist(x, 'file'));
            ip.addParamValue('avgVals', [], @(x) isstruct(x));
            ip.addParamValue('avgValsPath', '', @(x) ischar(x));
            ip.addParamValue('distsPath', '', @(x) ischar(x));
            ip.addParamValue('refParameterVals', [], @(x) isstruct(x));
            ip.addParamValue('refParameterValsPath', '', @(x) exist(x, 'file'));
            ip.addParamValue('refAvgVals', [], @(x) isstruct(x));
            ip.addParamValue('refAvgValsPath', '', @(x) exist(x, 'file'));
            
            ip.parse(varargin{:});
            
            parameterVals        = ip.Results.parameterVals;
            parameterValsPath    = ip.Results.parameterValsPath;
            avgVals              = ip.Results.avgVals;
            avgValsPath          = ip.Results.avgValsPath;
            distsPath            = ip.Results.distsPath;
            refParameterVals     = ip.Results.refParameterVals;
            refParameterValsPath = ip.Results.refParameterValsPath;
            refAvgVals           = ip.Results.refAvgVals;
            refAvgValsPath       = ip.Results.refAvgValsPath;
            
            %load parameter values
            refParameterVals = DreamScoring.loadParameterVals(refParameterVals, refParameterValsPath);
            parameterVals = DreamScoring.loadParameterVals(parameterVals, parameterValsPath);
            
            %load average values
            avgVals = DreamScoring.loadAvgVals(avgVals, avgValsPath);
            refAvgVals = DreamScoring.loadAvgVals(refAvgVals, refAvgValsPath);
            
            %% distances
            dists = struct('parameter', [], 'prediction', []);
            
            %parameters
            dists.parameter = DreamScoring.calcAvgSumSquaredLogRatio(...
                DreamScoring.getParameterVector(parameterVals), ...
                DreamScoring.getParameterVector(refParameterVals));
            
            %prediction
            dists.prediction = DreamScoring.calcAvgNormSquaredDiff(...
                DreamScoring.getPredictionVector(avgVals), ...
                DreamScoring.getPredictionVector(refAvgVals));
            
            if ~isempty(distsPath)
                save(distsPath, '-struct', 'dists');
            end
        end
    end
    
    %helper methods
    methods (Static = true, Access = protected)        
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
        end
        
        
        %Get struct of averages of high-throughput in silico experiments
        function avgVals = loadAvgVals(avgVals, avgValsPath)
            if ~isempty(avgVals)
            elseif ~isempty(avgValsPath)
                avgVals = load(avgValsPath);
            end
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
        end
        
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
        end
        
        %Calculate <(log(test/ref))^2>
        function dist = calcAvgSumSquaredLogRatio(test, ref)
            test = max(-realmax, min(realmax, test));
            ref  = max(-realmax, min(realmax, ref));
            
            dist = mean(log10(test ./ ref) .^ 2);
            
            dist = full(dist);
        end
        
        %Calculate <((mean(ref) - mean(test))^2) / var(ref)>
        function dist = calcAvgNormSquaredDiff(test, ref)
            test.mean(isinf(test.mean)) = sign(test.mean(isinf(test.mean))) * Inf;
            ref.mean(isinf(ref.mean)) = sign(ref.mean(isinf(ref.mean))) * Inf;
            test.std(isinf(test.std)) = sign(test.std(isinf(test.std))) * Inf;
            ref.std(isinf(ref.std)) = sign(ref.std(isinf(ref.std))) * Inf;            
            
            dist = mean(((ref.mean - test.mean) .^ 2) ./ (ref.std .^ 2));
            
            dist = full(dist);
        end        
        
        %calculate p-value
        function p = calcPVal(tests, distrib)
            [sortedTestAndDistrib, idx] = sortrows([
                tests zeros(size(tests))
                distrib ones(size(distrib))
                ], [-1 -2]);
            percentiles = zeros(size(sortedTestAndDistrib, 1), 1);
            percentiles(idx) = cumsum(sortedTestAndDistrib(:, 2)) / numel(distrib);
            
            p = percentiles(1:numel(tests));
            p(p == 0) = 1 / numel(distrib);
        end
    end
end