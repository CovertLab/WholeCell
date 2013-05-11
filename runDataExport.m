% runDataExport.m
%   Exports simulation data in JSON format
%
% Author: Ruby Lee
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 8/22/2012
function runDataExport(simBatch, simIdx)
setWarnings();
setPath();
setPreferences();
warning('off', 'MATLAB:MKDIR:DirectoryExists');

%import classes
import edu.stanford.covert.cell.sim.util.DiskLogger;
import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
import edu.stanford.covert.util.SparseMat;

%options
verbosity = 1;
exportProps = {
    'Chromosome', 'abasicSites'
    'Chromosome', 'complexBoundSites'
    'Chromosome', 'damagedBases'
    'Chromosome', 'damagedSugarPhosphates'
    'Chromosome', 'doubleStrandedRegions'
    'Chromosome', 'gapSites'
    'Chromosome', 'hollidayJunctions'
    'Chromosome', 'intrastrandCrossLinks'
    'Chromosome', 'linkingNumbers'
    'Chromosome', 'monomerBoundSites'
    'Chromosome', 'ploidy'
    'Chromosome', 'polymerizedRegions'
    'Chromosome', 'segregated'
    'Chromosome', 'singleStrandedRegions'
    'Chromosome', 'strandBreaks'
    'Chromosome', 'superhelicalDensity'
    'FtsZRing', 'numEdges'
    'FtsZRing', 'numEdgesOneStraight'
    'FtsZRing', 'numEdgesTwoBent'
    'FtsZRing', 'numEdgesTwoStraight'
    'FtsZRing', 'numResidualBent'
    'Geometry', 'chamberVolume'
    'Geometry', 'cylindricalLength'
    'Geometry', 'pinchedDiameter'
    'Geometry', 'surfaceArea'
    'Geometry', 'totalLength'
    'Geometry', 'volume'
    'Geometry', 'width'
    'Host', 'isBacteriumAdherent'
    'Host', 'isInflammatoryResponseActivated'
    'Host', 'isNFkBActivated'
    'Host', 'isTLRActivated'
    'Mass', 'cell'
    'Mass', 'cellDry'
    'Mass', 'dnaWt'
    'Mass', 'media'
    'Mass', 'metaboliteWt'
    'Mass', 'proteinWt'
    'Mass', 'rnaWt'
    'Mass', 'total'
    'Mass', 'waterWt'
    'MetabolicReaction', 'doublingTime'
    'MetabolicReaction', 'fluxs'
    'MetabolicReaction', 'growth'
    'Metabolite', 'counts'
    'ProteinComplex', 'counts'
    'ProteinMonomer', 'counts'
    'Ribosome', 'boundMRNAs'
    'Ribosome', 'mRNAPositions'
    'Ribosome', 'nActive'
    'Ribosome', 'nNotExist'
    'Ribosome', 'nStalled'
    'Ribosome', 'stateOccupancies'
    'Ribosome', 'states'
    'Ribosome', 'tmRNAPositions'
    'Rna', 'counts'
    'RNAPolymerase', 'nActive'
    'RNAPolymerase', 'nFree'
    'RNAPolymerase', 'nNonSpecificallyBound'
    'RNAPolymerase', 'nSpecificallyBound'
    'RNAPolymerase', 'positionStrands'
    'RNAPolymerase', 'stateOccupancies'
    'RNAPolymerase', 'states'
    'RNAPolymerase', 'supercoilingBindingProbFoldChange'
    'RNAPolymerase', 'transcriptionFactorBindingProbFoldChange'
    'Stimulus', 'values'
    'Transcript', 'abortedTranscripts'
    'Transcript', 'boundTranscriptChromosome'
    'Transcript', 'boundTranscriptionUnits'
    'Transcript', 'boundTranscriptProgress'
    'Transcript', 'rnaBoundRNAPolymerases'
    };

%sanitize arguments
if ischar(simIdx)
    simIdx = str2double(simIdx);
end

%load simulation object which contains meta data
[~, ~, sim] = SimulationDiskUtil.getSimulation([simBatch filesep num2str(simIdx)]);

%export data to JSON
errorLog = [];
outDir = [SimulationDiskUtil.getBaseDir() filesep simBatch filesep num2str(simIdx) filesep 'json'];

%create output directory
mkdir(outDir);

%summary data
if verbosity > 0, fprintf('Exporting summary ...'), end
summary = load(sprintf('%s%s%s%s%d%ssummary.mat', SimulationDiskUtil.getBaseDir(), filesep, simBatch, filesep, simIdx, filesep));
mkdir([outDir filesep 'Summary']);
fields = fieldnames(summary);
for i = 1:numel(fields)
    data = summary.(fields{i});
    for j = 1:size(data, 1)
        tmp = struct('label', [], 'data', []);
        tmp.label = [fields{i} '_' num2str(j)];
        tmp.data = compressData(data(j, :), 1, 1000);
        
        writeJSON(tmp, [outDir filesep 'Summary' filesep fields{i} '_' num2str(j) '.json']);
    end
    
    clear data;
end
clear summary;
if verbosity > 0, fprintf(' done\n'), end

%complete data
for i = 1:numel(sim.states)
    stateID = sim.states{i}.wholeCellModelID(7:end);
    propIDs = [sim.states{i}.stateNames; sim.states{i}.dependentStateNames];
    if verbosity > 0, fprintf('Exporting %s ...\n', stateID), end
    mkdir([outDir filesep stateID]);
    for j = 1:numel(propIDs)
        if ~any(strcmp(exportProps(:, 1), stateID) & strcmp(exportProps(:, 2), propIDs{j}))
            continue;
        end
        
        if verbosity > 0, fprintf('\t%s ...', propIDs{j}), end
        %try
        %retrieve data
        result = SimulationEnsemble.load(simBatch, {stateID propIDs{j}}, ...
            [], [], 1, 'extract', simIdx);
        result = result.(stateID).(propIDs{j});
        
        %output data in JSON format
        if isa(sim.states{i}, 'edu.stanford.covert.cell.sim.state.Chromosome') && isa(result, 'SparseMat')
            %chromosome properties
            [subs, vals] = find(result);
            tfs1 = subs(:, 3) < size(result, 3);
            tfs2 = subs(:, 3) > 1;
            
            result1 = SparseMat(subs(tfs1, :), vals(tfs1, :), size(result) - [0 0 1]);
            result2 = SparseMat([subs(tfs2, 1:2) subs(tfs2, 3)-1], vals(tfs2, :), size(result) - [0 0 1]);
            
            times = downsampleData([1; 1 + find(full(permute(any(any(result1 ~= result2, 1), 2), [3 1 2]))); size(result, 3)], 1000);
            if numel(times) > 1000
                times = downsampleData((1:size(result, 3))', 1000);
            end
            
            tmp = struct('time', [], 'strand', [], 'pos', [], 'val', []);
            tmp = tmp(ones(numel(times), 1), 1);
            for k = 1:numel(times)
                [tmpsubs, tmpvals] = find(result(:, :, times(k)));
                
                tmp(k).time = times(k);
                tmp(k).pos = tmpsubs(:, 1);
                tmp(k).strand = tmpsubs(:, 2);
                tmp(k).val = tmpvals;
            end
            
            writeJSON(tmp, [outDir filesep stateID filesep propIDs{j} '.json']);
        else
            %generic output
            
            %sum over compartments
            if ...
                    (isa(sim.states{i}, 'edu.stanford.covert.cell.sim.MoleculeCountState') && strcmp(propIDs{j}, 'counts')) || ...
                    (isa(sim.states{i}, 'edu.stanford.covert.cell.sim.state.Stimulus') && strcmp(propIDs{j}, 'values'))
                result = sum(result, 2);
            end
            
            %output
            if size(result, 2) == 1
                %save rows separately
                for k = 1:size(result, 1)
                    tmp = struct('label', [], 'data', []);
                    tmp.label = [stateID '_' propIDs{j} '_' num2str(k)];
                    tmp.data = compressData(result(k, :, :), 0, 1000);
                    
                    writeJSON(tmp, [outDir filesep stateID filesep propIDs{j} '_' num2str(k) '.json']);
                end
            else
                %save rows and columns separately
                for k = 1:size(result, 1)
                    for l = 1:size(result, 2)
                        tmp = struct('label', [], 'data', []);
                        tmp.label = [stateID '_' propIDs{j} '_' num2str(k) '_' num2str(l)];
                        tmp.data = compressData(result(k, l, :), 0, 1000);
                        
                        writeJSON(tmp, [outDir filesep stateID filesep propIDs{j} '_' num2str(k) '_' num2str(l) '.json']);
                    end
                end
            end
        end
        
        clear result tmp;
        %catch exception
        %    errorLog = [errorLog sprintf('Error in exporting %s in simulation %s/%s\n%s', ...
        %        propIDs{j}, simBatch, simIdx, exception.getReport())]; %#ok<AGROW>
        %end
        if verbosity > 0, fprintf(' done\n'), end
    end
    if verbosity > 0, fprintf('done\n'), end
end

%print error log
if ~isempty(errorLog)
    throw(MException('runAnalysis:error', errorLog));
end

function result = compressData(data, timeOffset, maxSamples)
data = full(double(data));
data = data(:);

%get time points where value changes
idxs = unique([1; find(diff(data(1:end-1))); 1 + find(diff(data(1:end-1))); numel(data)]);

if numel(idxs) > maxSamples
    %downsample
    idxs = downsampleData((1:numel(data))', maxSamples);
end

%return as time, value pairs
result = [idxs-timeOffset  data(idxs)];

function idxs = downsampleData(idxs, maxSamples)
if numel(idxs) > maxSamples
    tmp = idxs(1:floor(numel(idxs)/maxSamples):end);
    if tmp(end) ~= idxs(end)
        tmp = [tmp; idxs(end)];
    end
    idxs = tmp;
end

function writeJSON(data, fileName)
fid = fopen(fileName, 'w+');
writeJSON_data(fid, data);
fclose(fid);

function writeJSON_data(fid, data)
if ischar(data)
    writeJSON_char(fid, data)
elseif isnumeric(data)
    writeJSON_numeric(fid, data)
elseif isstruct(data)
    writeJSON_struct(fid, data)
else
    throw(MException('writeJSON_data:error', 'Unsupported data type %s', class(data)))
end

function writeJSON_struct(fid, data)
if numel(data) == 1
    fprintf(fid, '{');
    fields = fieldnames(data);
    for i = 1:numel(fields)
        fprintf(fid, '"%s": ', fields{i});
        writeJSON_data(fid, data.(fields{i}));
        if i < numel(fields)
            fprintf(fid, ', ');
        end
    end
    fprintf(fid, '}');
else
    fprintf(fid, '[');
    for i = 1:numel(data)
        writeJSON_struct(fid, data(i));
        if i < numel(data)
            fprintf(fid, ',');
        end
    end
    fprintf(fid, ']');
end

function writeJSON_char(fid, data)
fprintf(fid, '"%s"', data);

function writeJSON_numeric(fid, data)
fprintf(fid, '[');

if ndims(data) > 2
    throw(MException('writeJSON:error', 'Data with more than 2 dimensions is not supported'))
end

if isempty(data)
elseif size(data, 2) == 1
    if all(abs(data - round(data)) < 1e-10)
        fprintf(fid, '%d', round(data(1)));
        if numel(data) > 1
            fprintf(fid, ',%d', round(data(2:end)));
        end
    else
        fprintf(fid, '%e', data(1));
        if numel(data) > 1
            fprintf(fid, ',%e', data(2:end));
        end
    end
elseif size(data, 2) == 2 && all(abs(data(:, 1) - round(data(:, 1))) < 1e-10)
    for j = 1:size(data, 1)
        fprintf(fid, '[');
        fprintf(fid, '%d', round(data(j, 1)));
        fprintf(fid, ',%e', data(j, 2:end));
        fprintf(fid, ']');
        if j < size(data, 1)
            fprintf(fid, ',');
        end
    end
else
    for j = 1:size(data, 1)
        fprintf(fid, '[');
        fprintf(fid, '%e', data(j, 1));
        fprintf(fid, ',%e', data(j, 2:end));
        fprintf(fid, ']');
        if j < size(data, 1)
            fprintf(fid, ',');
        end
    end
end

fprintf(fid, ']');