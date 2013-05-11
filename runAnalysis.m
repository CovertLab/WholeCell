% Analyzes whole cell simulations.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 6/30/2011
function runAnalysis(varargin)
setWarnings();
setPath();
setPreferences();

import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;

outputDirectory = varargin{1};
if numel(varargin) < 2
    iJob = 1;
    nJobs = 1;
else
	iJob = varargin{2};
    nJobs = varargin{3};	
	if ischar(iJob)
		iJob = str2double(iJob);
		nJobs = str2double(nJobs);
	end
end

simNum = 1;
simDir = [outputDirectory filesep num2str(simNum)];
[~, ~, sim] = SimulationDiskUtil.getSimulation(simDir);

%cleanup any xls files to avoid slyk warning messages on windows
if ispc
    delete([outputDirectory filesep '*.xls']);
end

%choose analyses to run
[~, options] = SimulationDiskUtil.getSimulations(outputDirectory);
allSingleGeneDeletions = true;
for i = 1:numel(options)
    if ...
            isempty(isscalar(options(i).geneticKnockouts)) || ...
            ~isscalar(options(i).geneticKnockouts) || ...
            ~isempty(options(i).media) || ...
            ~isempty(options(i).stimulus)
        allSingleGeneDeletions = false;
    end
end

if allSingleGeneDeletions
    analysisFunc = {
        %single gene deletions
        'analysis'  'SingleGeneDeletions'           'run'                   {outputDirectory, false, iJob, nJobs}
        };
else
    analysisFunc = {
        %static
        'analysis'  'BiomassCompositionProduction'  'run'                   {sim, [outputDirectory filesep 'BiomassCompositionProduction']}
        'analysis'  'CellState'                     'run'                   {sim, [outputDirectory filesep 'CellState']}
        'analysis'  'Constants'                     'run'                   {sim, [outputDirectory filesep 'Constants']}
        'analysis'  'DNADamage'                     'run'                   {sim, [outputDirectory filesep 'DNADamage']}
        'analysis'  'FBA'                           'run'                   {sim, [outputDirectory filesep 'FBA']}
        'analysis'  'SimulationStructure'           'run'                   {sim, [outputDirectory filesep 'SimulationStructure']}
        
        %single cell
        'analysis'  'TranslationAnalysis'           'run'                   {simDir, [outputDirectory filesep 'TranslationAnalysis']}
        'analysis'  'ChromosomeSpaceTimePlot'       'run'                   {simDir, [], [outputDirectory filesep 'ChromosomeSpaceTimePlot']}
        'analysis'  'ChromosomePositionHistogram'   'run'                   {simDir, [outputDirectory filesep 'ChromosomePositionHistogram']}
        
        %population
        'util'      'SummaryLogger'                 'summarizeSimulations'  {outputDirectory}
        'analysis'  'CellOverview'                  'run'                   {outputDirectory, [outputDirectory filesep 'summary']}
        'analysis'  'ProcessMetaboliteUsage'        'run'                   {outputDirectory, [], [outputDirectory filesep 'processMetaboliteUsage']}
        };
    analysisFunc = analysisFunc(iJob:nJobs:end, :);
    
    md1 = meta.class.fromName('edu.stanford.covert.cell.sim.analysis.SingleCell');
    md2 = meta.class.fromName('edu.stanford.covert.cell.sim.analysis.Population');
    mdNames1 = cellfun(@(md) md.Name, md1.Methods, 'UniformOutput', false);
    mdNames2 = cellfun(@(md) md.Name, md2.Methods, 'UniformOutput', false);
    analysisFunc = [
        analysisFunc; {
        'analysis'  'SingleCell'  'run'  {outputDirectory, simNum, [outputDirectory filesep 'singleCell'], mdNames1(iJob:nJobs:end)}
        'analysis'  'Population'  'run'  {outputDirectory, [outputDirectory filesep 'population'], [], mdNames2(iJob:nJobs:end)}
        }];
end

%run analyses and log any errors
errorLog = [];
for i = 1:size(analysisFunc, 1)
    try
        edu.stanford.covert.cell.sim.(analysisFunc{i, 1}).(analysisFunc{i, 2}).(analysisFunc{i, 3})(analysisFunc{i, 4}{:});
    catch exception
        errorLog = [errorLog sprintf('Error in %s.%s analysis\n%s', analysisFunc{i, 2}, analysisFunc{i, 3}, exception.getReport())]; %#ok<AGROW>
    end
end

%print error log
if ~isempty(errorLog)
    throw(MException('runAnalysis:error', errorLog));
end
