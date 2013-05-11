% SimulationEnsemble
%   Retrieves data from an ensemble (a batch) of simulations.
%   Data is stored and annotated
%   Statistics can be computed on the data
%
% Author: Derek Macklin, macklin@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Creation Date: 7/12/2011
% Last Updated: 8/1/2011

classdef SimulationEnsemble < handle
    properties
        stateData
    end
    
    properties (SetAccess = protected)
        directory
        timeSpan
        selectedSimulations
    end
    
    properties (Constant = true)
        summaryLogToEnsemble = {
            'runTime'           'runTime'               1
            'mass'              'mass'                  1
            'growth_rate'       'growth'                1
            'ploidy'            'ploidy'                1
            'atp'               'metabolites'           1
            'adp'               'metabolites'           2
            'amp'               'metabolites'           3
            'ntps'              'metabolites'           4
            'amino_acids'       'metabolites'           7
            'lipids'            'metabolites'          10
            'polysaccharides'   'metabolites'          11
            'dntps'             'metabolites'          12
            'antibiotics'       'metabolites'          13
            'rnas'              'rnas'                  1
            'mrnas'             'rnas'                  2
            'rrnas'             'rnas'                  3
            'srnas'             'rnas'                  4
            'trnas'             'rnas'                  5
            'immatureRnas'      'rnas'                  6
            'damagedRnas'       'rnas'                  7
            'monomers'          'proteins'              1
            'matureMonomers'    'proteins'              2
            'immatureMonomers'  'proteins'              3
            'complexs'          'proteins'              4
            'matureComplexs'    'proteins'              5
            'immatureComplexs'  'proteins'              6
            'damagedProteins'   'proteins'              7
            'proteins'          'proteins'              [1 4]
            'rnaPolymerases'    'rnaPolymerases'        1
            'ribosomes'         'ribosomes'             1
            'dnaA'              'dnaA'                  1:2
            'dnaAATP_free'      'dnaA'                  1
            'dnaA_total'        'dnaA'                  2
            'dnaA_box1'         'dnaAFunctionalBoxes'   1
            'dnaA_box2'         'dnaAFunctionalBoxes'   2
            'dnaA_box3'         'dnaAFunctionalBoxes'   3
            'dnaA_box4'         'dnaAFunctionalBoxes'   4
            'dnaA_box5'         'dnaAFunctionalBoxes'   5
            'dnaA_boxes'        'dnaAFunctionalBoxes'   1:5
            'helicase1'         'replisome'             1
            'helicase2'         'replisome'             2
            'leadingPol1'       'replisome'             3
            'leadingPol2'       'replisome'             4
            'laggingPol1'       'replisome'             5
            'laggingPol2'       'replisome'             6
            'ftsZ'              'ftsZ'                  2
            'ftsZRing1st'       'ftsZRing'              1
            'ftsZRing2st'       'ftsZRing'              2
            'ftsZRing2bt'       'ftsZRing'              3
            'ftsZRingRbt'       'ftsZRing'              4
            'pinchedDiameter'   'pinchedDiameter'       1
            'supercoiled'       'supercoiled'           1
            'superhelicity'     'supercoiled'           2
            'atpUsage'          'metaboliteUsages'      1
            'gtpUsage'          'metaboliteUsages'      2
            'atpProduction'     'metaboliteProductions' 1
            'gtpProduction'     'metaboliteProductions' 2
            };
    end
        
    methods (Static = true)
        function incompSims = remainingSimulations(time, endTimes)
            % Data for a Kaplan-Meier curve
            incompSims = zeros(1, length(time));
            endTimes = sort(endTimes(~isnan(endTimes)));
            nSims = length(endTimes);
            
            r1 = 1;
            r2 = find(time >= endTimes(1), 1);
            incompSims(r1:r2) = nSims;
            nSims = nSims - 1;
            
            for i = 1:length(endTimes) - 1;
                r1 = find(time >= endTimes(i), 1);
                r2 = find(time >= endTimes(i+1), 1);
                incompSims(r1:r2) = nSims;
                nSims = nSims - 1;
            end
            
            incompSims(find(time >= endTimes(end), 1):end) = nSims;
        end
        
        function sData = createStateDataStruct(numTimePoints, nSims, selectedStates)
            
            numProperties = size(selectedStates, 1);
            
            sData = struct;
            sData.time = nan(1, numTimePoints);
            sData.downsampleStepSec = NaN;
            sData.values = nan(numProperties, 1, numTimePoints, nSims);
            sData.simulationNumber = nan(nSims, 1);
            sData.simulationStatus = nan(nSims, 1);
            sData.simulationEndTimes = nan(nSims, 1);
            sData.properties = selectedStates(:, 1);
            if size(selectedStates, 2) == 1
                sData.properties = selectedStates(:, 1);
            else
                sData.descriptions = selectedStates(:, 2);
            end
        end
    end
    
    methods
        function this = SimulationEnsemble(simID, selectedStates, timeSpan, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            if nargin <= 1
                throw(MException('SimulationEnsemble:error', 'Insufficient arguments'));
            end
            
            simDir = SimulationDiskUtil.getSimulation([simID filesep '1']);
            this.directory = SimulationDiskUtil.getSimulationBatchDir(simDir);
            
            this.timeSpan = [0 Inf];
            this.selectedSimulations = SimulationDiskUtil.getCompleteSimulations(this.directory);
            
            if nargin >= 3 && ~isempty(timeSpan)
                this.timeSpan = timeSpan;
            end
            
            if nargin >= 4
                this.selectedSimulations = selectedSimulations;
            end
            
            this.stateData = SimulationEnsemble.createStateDataStruct(...
                this.maximumSimulationTimePoints(), numel(this.selectedSimulations), selectedStates);
            
            this.initializeStateData();
        end
        
        function initializeStateData(this)
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            %identify properties which have (and haven't) been saved in summary
            [tfs, ensemblePropIdxs] = ismember(this.stateData.properties, this.summaryLogToEnsemble(:, 1));
            
            if ~isempty(this.selectedSimulations)
                summaryFileName = sprintf('%s/%d/summary.mat', this.directory, this.selectedSimulations(1));
                if exist(summaryFileName, 'file')
                    summary = load(summaryFileName);
                    tfs(tfs) = ismember(this.summaryLogToEnsemble(ensemblePropIdxs(tfs), 2), fieldnames(summary));
                    clear summary;
                end
            end
            
            if ~all(tfs)
                warning('WholeCell:warning:undefinedSummaryProperty', 'Undefined properties:\n- %s', strjoin(sprintf('\n- '), this.stateData.properties{~tfs}));
            end
            
            %get data
            for i = 1:length(this.selectedSimulations)
                this.stateData.simulationNumber(i) = this.selectedSimulations(i);
                summaryFileName = sprintf('%s/%d/summary.mat', this.directory, this.selectedSimulations(i));
                
                if ~exist(sprintf('%s/%d/state-0.mat', this.directory, this.selectedSimulations(i)), 'file')
                    this.stateData.simulationStatus(i) = SummaryLogger.SIM_STATUS_DIDNT_START;
                elseif ~exist(sprintf('%s/%d/err.mat', this.directory, this.selectedSimulations(i)), 'file') && ~exist(sprintf('%s/%d/summary.mat', this.directory, this.selectedSimulations(i)), 'file')
                    this.stateData.simulationStatus(i) = SummaryLogger.SIM_STATUS_STILL_RUNNING;
                elseif ~exist(summaryFileName, 'file')
                    this.stateData.simulationStatus(i) = SummaryLogger.SIM_STATUS_TERMINATED_WITH_ERROR;
                else                    
                    summary = load(summaryFileName);
                    
                    if exist(sprintf('%s/%d/err.mat', this.directory, this.selectedSimulations(i)), 'file')
                        this.stateData.simulationStatus(i) = SummaryLogger.SIM_STATUS_TERMINATED_WITH_ERROR;
                    elseif summary.pinchedDiameter(end) == 0
                        this.stateData.simulationStatus(i) = SummaryLogger.SIM_STATUS_COMPLETED_WITH_DIVISION;
                    else
                        this.stateData.simulationStatus(i) = SummaryLogger.SIM_STATUS_COMPLETED_WITHOUT_DIVISION;
                    end
                    
                    indStart = find(summary.time >= this.timeSpan(1), 1, 'first');
                    indStop = find(summary.time <= this.timeSpan(2), 1, 'last');
                    
                    if numel(summary.time(indStart:indStop)) > sum(~isnan(this.stateData.time))
                        this.stateData.time(1:(indStop - indStart + 1)) = summary.time(indStart:indStop);
                    end
                    this.stateData.downsampleStepSec = diff(this.stateData.time(1:2));
                    
                    this.stateData.simulationEndTimes(i) = summary.time(end);
                    
                    for j = 1:length(this.stateData.properties)
                        if ~tfs(j)
                            continue;
                        end
                        
                        results = double(summary.(this.summaryLogToEnsemble{ensemblePropIdxs(j), 2})(this.summaryLogToEnsemble{ensemblePropIdxs(j), 3}, indStart:indStop));
                        if size(results, 1) > 1
                            results = sum(results, 1);
                        end
                        
                        if isequal(this.summaryLogToEnsemble{ensemblePropIdxs(j), 2}, 'replisome')
                            if isodd(this.summaryLogToEnsemble{ensemblePropIdxs(j), 3})
                                results(results == 0) = NaN;
                                a = find(results < 580076 / 3);     % Ideally use 2, but want room for noise (like footprint size)--this is overkill
                                results(a) = results(a) + 580076;
                            else
                                results(results == 0) = NaN;
                                a = find(results > 580076 * 2 / 3); % Ideally use 1/2, but want room for noise (like footprint size)---this is overkill
                                results(a) = results(a) - 580076;
                            end
                        end
                        this.stateData.values(j, 1, 1:length(results), i) = results;
                    end
                    
                    clear summary;
                end
            end
        end
        
        function value = maximumSimulationTimePoints(this)
            value = 0;
            for i = 1:length(this.selectedSimulations)
                try %#ok<TRYNC>
                    summary = load([this.directory filesep num2str(this.selectedSimulations(i)) filesep 'summary.mat'], 'time');
                    value = max(value, sum(summary.time >= this.timeSpan(1) & summary.time <= this.timeSpan(2)));
                end
            end
        end
        
        % 'props' is a cell array of properties
        function value = getPropertyIndices(this, props)
            [~, value] = ismember(props, this.stateData.properties);
            value(value == 0) = NaN;
        end
        
        function value = getCompletedWithReplication(this)
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            value = sum(this.stateData.simulationStatus == SummaryLogger.SIM_STATUS_COMPLETED_WITH_DIVISION);
        end
        
        function value = getCompletedWithoutError(this)
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            value = sum(...
                this.stateData.simulationStatus == SummaryLogger.SIM_STATUS_COMPLETED_WITH_DIVISION | ...
                this.stateData.simulationStatus == SummaryLogger.SIM_STATUS_COMPLETED_WITHOUT_DIVISION);
        end
        
        function value = getEarliestTermination(this)
            value = min(this.stateData.simulationEndTimes);
        end
        
        function value = getTimeStamp(this)
            value = edu.stanford.covert.cell.sim.util.SimulationDiskUtil.getSimulationTimeStamp([this.directory filesep '1']);
        end
    end
    
    methods (Static = true)
        function [states, metadata, options, parameters, fittedConstants] = load(simGroupDir, stateNames, ...
                initTime, finTime, downsampleStepSec, ...
                downsampleType, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            if ~((isunix && simGroupDir(1) == '/') || (ispc && simGroupDir(2) == ':'))
                simGroupDir = [SimulationDiskUtil.getBaseDir() filesep simGroupDir];
            end
            
            simDir = [simGroupDir filesep '1'];
            metadata = DiskLogger.loadMetadata(simDir); metadata.downsampleStepSec = downsampleStepSec;
            options = DiskLogger.loadOptions(simDir);
            parameters = DiskLogger.loadParameters(simDir);
            fittedConstants = DiskLogger.loadFittedConstants(simDir);
            
            if ~exist('selectedSimulations', 'var') || isempty(selectedSimulations)
                selectedSimulations = SimulationDiskUtil.getCompleteSimulations(simGroupDir);
            end
            nSims = numel(selectedSimulations);
            
            lengthSecs = zeros(nSims, 1);
            states = struct();
            for j = 1:size(stateNames, 1)
                for i = 1:nSims
                    simDir = [simGroupDir filesep num2str(selectedSimulations(i))];
                    if j == 1
                        tmp = DiskLogger.loadMetadata(simDir, 'lengthSec');
                        lengthSecs(i) = tmp.lengthSec;
                    end
                    metadata.lengthSec = lengthSecs(i);
                    tmpStates = DiskLogger.loadTimecourses(simDir, stateNames, initTime, finTime, downsampleStepSec, ...
                        downsampleType, options, metadata);
                    
                    if i == 1
                        states.(stateNames{j, 1}).(stateNames{j, 2}) = tmpStates.(stateNames{j, 1}).(stateNames{j, 2});
                        continue;
                    end
                    
                    tmp1 = states.(stateNames{j, 1}).(stateNames{j, 2});
                    tmp2 = tmpStates.(stateNames{j, 1}).(stateNames{j, 2});
                    siz1 = size(tmp1); siz1 = [siz1 ones(1, 4-numel(siz1))]; %#ok<AGROW>
                    siz2 = size(tmp2); siz2 = [siz2 ones(1, 4-numel(siz2))]; %#ok<AGROW>
                    if isa(states.(stateNames{j, 1}).(stateNames{j, 2}), 'edu.stanford.covert.util.SparseMat')
                        if ~isequal(siz1(1:3), siz2(1:3))
                            states.(stateNames{j, 1}).(stateNames{j, 2}) = cat(1, states.(stateNames{j, 1}).(stateNames{j, 2}), DiskLogger.allocateData(class(tmp1), [max(0, siz2(1) - siz1(1)) siz1(2:end)]));
                            states.(stateNames{j, 1}).(stateNames{j, 2}) = cat(2, states.(stateNames{j, 1}).(stateNames{j, 2}), DiskLogger.allocateData(class(tmp1), [max(siz2(1), siz1(1)) max(0, siz2(2) - siz1(2)) siz1(3:end)]));
                            states.(stateNames{j, 1}).(stateNames{j, 2}) = cat(3, states.(stateNames{j, 1}).(stateNames{j, 2}), DiskLogger.allocateData(class(tmp1), [max(siz2(1:2), siz1(1:2)) max(0, siz2(3) - siz1(3)) siz1(4:end)]));
                            
                            tmpStates.(stateNames{j, 1}).(stateNames{j, 2}) = cat(1, tmpStates.(stateNames{j, 1}).(stateNames{j, 2}), DiskLogger.allocateData(class(tmp1), [max(0, siz1(1) - siz2(1)) siz2(2:end)]));
                            tmpStates.(stateNames{j, 1}).(stateNames{j, 2}) = cat(2, tmpStates.(stateNames{j, 1}).(stateNames{j, 2}), DiskLogger.allocateData(class(tmp1), [max(siz1(1), siz2(1)) max(0, siz1(2) - siz2(2)) siz2(3:end)]));
                            tmpStates.(stateNames{j, 1}).(stateNames{j, 2}) = cat(3, tmpStates.(stateNames{j, 1}).(stateNames{j, 2}), DiskLogger.allocateData(class(tmp1), [max(siz1(1:2), siz2(1:2)) max(0, siz1(3) - siz2(3)) siz2(4:end)]));
                        end
                        
                        states.(stateNames{j, 1}).(stateNames{j, 2}) = cat(4, ...
                            states.(stateNames{j, 1}).(stateNames{j, 2}), tmpStates.(stateNames{j, 1}).(stateNames{j, 2}));
                    else
                        states.(stateNames{j, 1}).(stateNames{j, 2}) = DiskLogger.allocateData(class(tmp1), [max(siz1(1:3), siz2(1:3)) siz1(4)+1], nnz(tmp1) / numel(tmp1));
                        states.(stateNames{j, 1}).(stateNames{j, 2})(1:siz1(1), 1:siz1(2), 1:siz1(3), 1:siz1(4)) = tmp1;
                        states.(stateNames{j, 1}).(stateNames{j, 2})(1:siz2(1), 1:siz2(2), 1:siz2(3), end) = tmp2;
                    end
                    
                    clear tmpStates tmp1 tmp2 siz1 siz2;
                end
            end
        end
    end
end
