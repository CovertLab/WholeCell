% Generates an animation
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/6/2011
function runAnimation(varargin)
setWarnings();
setPath();
setPreferences();

%import classes
import edu.stanford.covert.cell.sim.analysis.OverviewAnimation;

%options
if nargin >= 1
    simBatchId = varargin{1};
end
if nargin >= 2
	movieSubFolder = varargin{2};
else 
	movieSubFolder = 'movie';
end
if nargin >= 3
    renderFrames = varargin{3};
    if ischar(renderFrames)
        renderFrames = str2double(renderFrames);
    end
end
if nargin >= 4
    jobIdx = varargin{4};
    if ischar(jobIdx)
        jobIdx = str2double(jobIdx);
    end
end
if nargin >= 5
    nJobs = varargin{5};
    if ischar(nJobs)
        nJobs = str2double(nJobs);
    end
end

%make animation
outDir = [edu.stanford.covert.cell.sim.util.SimulationDiskUtil.getBaseDir() filesep simBatchId filesep movieSubFolder];
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
animator = OverviewAnimation(simBatchId, 1:10, [outDir filesep 'movie.avi'], 50);
animator.tmpDirectory = outDir;
if renderFrames
    animator.run([jobIdx nJobs], true, true, false, false);
else
    animator.run([], false, false, true, false);
	
	animator.codec.ffmpeg = 'wmv2';
	animator.movieFileName = [outDir filesep 'movie.wmv'];
	animator.renderMovie();
end
