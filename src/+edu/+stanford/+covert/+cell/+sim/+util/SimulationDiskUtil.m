%SimulationDiskUtil
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Derek Macklin, macklin@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 7/6/2011
classdef SimulationDiskUtil
    methods (Static = true)
        function [simDir, revision, sim] = getSimulation(simDir)
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            if ~((isunix && simDir(1) == '/') || (ispc && simDir(2) == ':'))
                simDir = [edu.stanford.covert.cell.sim.util.SimulationDiskUtil.getBaseDir() filesep simDir];
            end
            
            if nargout >= 2
                metadata = DiskLogger.loadMetadata(simDir);
                revision = metadata.revision;
            end
            
            if nargout >= 3
                sim = CachedSimulationObjectUtil.load(revision);
            end
        end
        
        function [simGroupDir, revision, sim] = getLatestSimulationGroup()
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            baseDir = edu.stanford.covert.cell.sim.util.SimulationDiskUtil.getBaseDir(); %#ok<*PROP>
            
            files = dir([baseDir filesep '201*']);
            files = files([files.isdir]);
            fileIdx = [];
            for i = numel(files):-1:1
                if ~SimulationDiskUtil.getNumSimulations(files(i).name)
                    continue;
                end
                
                tmp = true;
                for j = 1:SimulationDiskUtil.getNumSimulations(files(i).name)
                    if ~exist([baseDir filesep files(i).name filesep num2str(j) filesep 'metadata.mat'], 'file')
                        tmp = false;
                        break;
                    end
                end
                
                if ~tmp
                    continue;
                end
                
                for j = 1:SimulationDiskUtil.getNumSimulations(files(i).name)
                    try
                        load([baseDir filesep files(i).name filesep num2str(j) filesep 'metadata.mat'])
                    catch %#ok<CTCH>
                        tmp = false;
                        break;
                    end
                end
                
                if ~tmp
                    continue;
                end
                
                fileIdx = i;
                break;
            end
            if isempty(fileIdx)
                simGroupDir = [];
                revision = [];
                sim = [];
                return;
            end
            simTimeStamp = files(fileIdx).name; %simulation group
            simGroupDir = sprintf('%s%s%s', baseDir, filesep, simTimeStamp);
            
            if nargout >= 2
                simIdx = 1; %simulation within group
                simDir = sprintf('%s%s%d', simGroupDir, filesep, simIdx);
                metadata = DiskLogger.loadMetadata(simDir);
                revision = metadata.revision;
            end
            
            if nargout >= 3
                sim = CachedSimulationObjectUtil.load(revision);
            end
        end
        
        function [simDir, revision, sim] = getLatestSimulation()
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            simGroupDir = SimulationDiskUtil.getLatestSimulationGroup();
            
            simIdx = 1; %simulation within group
            simDir = sprintf('%s%s%d', simGroupDir, filesep, simIdx);
            
            if nargout >= 2
                metadata = DiskLogger.loadMetadata(simDir);
                revision = metadata.revision;
            end
            
            if nargout >= 3
                sim = CachedSimulationObjectUtil.load(revision);
            end
        end
        
        function simGroup = getLatestWildTypeSimulationGroup()
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            baseDir = SimulationDiskUtil.getBaseDir();
            files = dir([baseDir filesep '201*']);
            files = files([files.isdir]);
            [~, order] = sort({files.name}); 
            files = files(order(end:-1:1));
            for i = 1:numel(files)
                for j = 1:SimulationDiskUtil.getNumSimulations(files(i).name)
                    simDir = [baseDir filesep files(i).name filesep num2str(j)];
                    simGroup = files(i).name;
                    
                    if ...
                            ~exist([simDir filesep 'metadata.mat'], 'file') || ...
                            ~exist([simDir filesep 'options.mat'], 'file') || ...
                            ~exist([simDir filesep 'parameters.mat'], 'file') || ...
                            exist([simDir filesep 'skipAnalysis'], 'file')
                        continue;
                    end
                    
                    options = load([simDir filesep 'options.mat']);
                    if isempty(options.geneticKnockouts)
                        return;
                    end
                end
            end
        end
        
        function simGroup = getLatestSingleGeneDeletionSimulationGroup()
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            baseDir = SimulationDiskUtil.getBaseDir();
            files = dir([baseDir filesep '201*']);
            files = files([files.isdir]);
            [~, order] = sort({files.name}); 
            files = files(order(end:-1:1));
            for i = 1:numel(files)
                for j = 1:SimulationDiskUtil.getNumSimulations(files(i).name)
                    simDir = [baseDir filesep files(i).name filesep num2str(j)];
                    simGroup = files(i).name;
                    
                    if ...
                            ~exist([simDir filesep 'metadata.mat'], 'file') || ...
                            ~exist([simDir filesep 'options.mat'], 'file') || ...
                            ~exist([simDir filesep 'parameters.mat'], 'file') || ...
                            exist([simDir filesep 'skipAnalysis'], 'file')
                        continue;
                    end
                    
                    options = load([simDir filesep 'options.mat']);
                    if isscalar(options.geneticKnockouts)
                        return;
                    end
                end
            end
        end
        
        function [metaData, options, parameters] = getSimulations(simBatchDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            baseDir = SimulationDiskUtil.getBaseDir(); %#ok<*PROP>
            
            if nargin >= 1
                simDir = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
                simBatchDir = simDir(1:find(simDir == filesep, 1, 'last')-1);
                files = struct('name', simBatchDir, 'isdir', ~~exist(simBatchDir, 'dir'));
            else
                files = dir([baseDir filesep '201*']);
            end
            files = files([files.isdir]);
            metaData = repmat(struct(...
                'firstName', [], ...
                'lastName', [], ...
                'email', [], ...
                'affiliation', [], ...
                'userName', [], ...
                'hostName', [], ...
                'ipAddress', [], ...
                'revision', [], ...
                'differencesFromRevision', [], ...
                'shortDescription', [], ...
                'longDescription', [], ...
                'knowledgeBaseWID', [], ...
                'startTime', [], ...
                'endTime', [], ...
                'lengthSec', [], ...
                'outputDirectory', [], ...
                'segmentSizeStep', [], ...
                'downsampleStepSec', [], ...
                'baseDirectory', [], ...
                'directory', [], ...
                'simGroup', [], ...
                'simIdx', [], ...
                'name', [], ...
                'stateNames', [], ...
                'dependentStateNames' ,[]), 0, 1);
            options = repmat(struct, 0, 1);
            parameters = repmat(struct, 0, 1);
            for i = 1:numel(files)
                if nargin < 2
                    selectedSimulations = 1:SimulationDiskUtil.getNumSimulations(files(i).name);
                end
                for j = 1:numel(selectedSimulations)
                    if (isunix && files(i).name(1) == '/') || (ispc && length(files(i).name) >= 2 && files(i).name(2) == ':')
                        simDir = [files(i).name filesep num2str(selectedSimulations(j))];
                        simDir = strrep(simDir, '/', filesep);
                        simDir = strrep(simDir, '\', filesep);
                        idxs = find(files(i).name == '/' | files(i).name == '\');
                        simGroup = files(i).name(idxs(end)+1:end);
                    else
                        simDir = [baseDir filesep files(i).name filesep num2str(selectedSimulations(j))];
                        simGroup = files(i).name;
                    end
                    
                    if exist([simDir filesep 'skipAnalysis'], 'file')
                        continue;
                    end
                    
                    if nargout >= 1 && exist([simDir filesep 'metadata.mat'], 'file')
                        try
                            tmpMetadata = load([simDir filesep 'metadata.mat']);
                        catch %#ok<CTCH>
                            continue;
                        end
                        tmpMetadata.baseDirectory = baseDir;
                        tmpMetadata.directory = simDir;
                        tmpMetadata.simGroup = simGroup;
                        tmpMetadata.simIdx = selectedSimulations(j);
                        tmpMetadata.name = [files(i).name ' - ' num2str(selectedSimulations(j)) ' - ' tmpMetadata.shortDescription];
                        
                        tmp = setdiff(fieldnames(metaData), fieldnames(tmpMetadata));
                        for k = 1:numel(tmp)
                            tmpMetadata.(tmp{k}) = [];
                        end
                        
                        metaData = [metaData; tmpMetadata]; %#ok<AGROW>
                    end
                    
                    if nargout >= 2 && exist([simDir filesep 'options.mat'], 'file')
                        tmpOptions = load([simDir filesep 'options.mat']);
                        
                        if isempty(options)
                            options = tmpOptions;
                        else
                            tmp = setdiff(fieldnames(options), fieldnames(tmpOptions));
                            for k = 1:numel(tmp)
                                tmpOptions.(tmp{k}) = [];
                            end
                            
                            tmp = setdiff(fieldnames(tmpOptions), fieldnames(options));
                            for k = 1:numel(tmp)
                                for l = 1:numel(options)
                                    options(l).(tmp{k}) = [];
                                end
                            end
                            
                            options = [options; tmpOptions]; %#ok<AGROW>
                        end
                    end
                    
                    if nargout >= 3 && exist([simDir filesep 'parameters.mat'], 'file')
                        tmpParameters = load([simDir filesep 'parameters.mat']);
                        
                        if isempty(parameters)
                            parameters = tmpParameters;
                        else
                            tmp = setdiff(fieldnames(parameters), fieldnames(tmpParameters));
                            for k = 1:numel(tmp)
                                tmpParameters.(tmp{k}) = [];
                            end
                            
                            tmp = setdiff(fieldnames(tmpParameters), fieldnames(parameters));
                            for k = 1:numel(tmp)
                                for l = 1:numel(parameters)
                                    parameters(l).(tmp{k}) = [];
                                end
                            end
                            
                            parameters = [parameters; tmpParameters]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
        
        function value = getBaseDir()
            config = getConfig();
            value = config.outputPath;
        end
        
        function value = getSimulationIndex(simDirAbs)
            % Get rid of any trailing slashes
            while simDirAbs(end) == '/'
                simDirAbs(end)=[];
            end
            
            p = find(simDirAbs == filesep, 1, 'last');
            value = str2double(simDirAbs(p+1:end));
        end
        
        function value = getSimulationTimeStamp(simDirAbs)
            % Get rid of any trailing slashes
            while simDirAbs(end) == filesep
                simDirAbs(end) = [];
            end
                
            p = find(simDirAbs == filesep, 2, 'last');
            value = simDirAbs(p(1) + 1:p(2) - 1);
			
			value = strrep(value, '_', ':');
			value(11) = ' ';
        end
        
        function value = getSimulationBatchDir(simDirAbs)
            % Get rid of any trailing slashes
            while simDirAbs(end) == '/'
                simDirAbs(end) = [];
            end
            
            p = find(simDirAbs == filesep, 2, 'last');
            value = simDirAbs(1:p(2) - 1);
            while value(end) == filesep
                value(end) = [];
            end
        end
        
        function value = getNumSimulations(simBatchDir)
			if (isunix && simBatchDir(1) == '/') || (ispc && simBatchDir(2) == ':')
				baseDir = simBatchDir;
			else
				baseDir = [edu.stanford.covert.cell.sim.util.SimulationDiskUtil.getBaseDir filesep simBatchDir];
			end
            files = dir(baseDir);
            files = files([files.isdir] & ~isnan(str2double({files.name})));
            value = numel(files);
        end
        
        function value = getCompleteSimulations(simBatchDir)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if (isunix && simBatchDir(1) == '/') || (ispc && simBatchDir(2) == ':')
                baseDir = simBatchDir;
            else
                baseDir = [edu.stanford.covert.cell.sim.util.SimulationDiskUtil.getBaseDir filesep simBatchDir];
            end
            
            nSims = SimulationDiskUtil.getNumSimulations(simBatchDir);
            tfs = false(nSims, 1);
            for i = 1:nSims
                tfs(i) = exist([baseDir filesep num2str(i) filesep 'summary.mat'], 'file');
            end
            value = find(tfs);
        end

        function warnings = getSimulationWarnings(simBatchDir, simIdx)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            fid = fopen([SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(simIdx)]) filesep 'out.log'], 'r');
            if fid == -1
                throw(MException('SimulationDiskUtil:error', 'Unable to open file'))
            end
            
            warningTemplate = struct();
            warningTemplate.message = [];
            warningTemplate.files = {};
            warningTemplate.lineNumbers = [];
            warningTemplate.times = [];
            warnings = repmat(warningTemplate, 0, 1);
            isWarn = false;
            isWarnTrace = false;
            time = 0;
            while ~feof(fid)
                line = fgetl(fid);
                
                if isWarnTrace && (length(line) < 5 || ~isequal(line(1:5), '  In '))
                    isWarn = false;
                    isWarnTrace = false;
                    
                    tmpWarning.message = tmpWarning.message(12:end-4);
                    
                    repeatedWarning = 0;
                    for j = 1:numel(warnings)
                        if isequal(warnings(j).message, tmpWarning.message) && ...
                                isequal(warnings(j).files, tmpWarning.files) && ...
                                isequal(warnings(j).lineNumbers, tmpWarning.lineNumbers)
                            repeatedWarning = j;
                            break;
                        end
                    end
                    
                    if repeatedWarning
                        warnings(repeatedWarning).times = ...
                            [warnings(repeatedWarning).times; tmpWarning.times];
                    else
                        warnings = [warnings; tmpWarning]; %#ok<AGROW>
                    end
                end
                if numel(line) >= 11 && isequal(line(3:11), 'Warning: ')
                    isWarn = true;
                    isWarnTrace = false;
                    tmpWarning = warningTemplate;
                end
                if isWarn && line(1) == '>'
                    isWarnTrace = true;
                end
                
                if ~isWarn && length(line) >= 1 && ~all(line == ' ') && any(line == ' ')
                    idx1 = find(line ~= ' ', 1, 'first');
                    idx2 = find(line(idx1:end) == ' ', 1, 'first') - idx1;
                    if ~isempty(idx1) && ~isempty(idx2) && idx2 >= idx1
                        time = max(time, str2double(line(idx1:idx2)));
                    end
                elseif isWarn && ~isWarnTrace
                    tmpWarning.times = time;
                    tmpWarning.message = [tmpWarning.message line ' '];
                elseif isWarnTrace
                    tmp = line(6:end);
                    tmpWarning.files{end+1, 1} = tmp(1:strfind(tmp, ' at ')-1);
                    tmpWarning.lineNumbers(end+1, 1) = str2double(tmp(strfind(tmp, ' at ')+4:end));
                end
            end
            
            assert(~isWarn);
            assert(~isWarnTrace);
            
            fclose(fid);
        end
    end
end
