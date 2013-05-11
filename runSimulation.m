% Executes whole cell simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/9/2011
function sim = runSimulation(varargin)
%parse options
if nargin >= 1 && isstruct(varargin{1})
    options = varargin{1};
else
    validateattributes(numel(varargin), {'numeric'}, {'even'});
    
    options = struct(varargin{:});
end

if isfield(options, 'logToDisk') && ischar(options.logToDisk)
    options.logToDisk = strcmp(options.logToDisk, 'true');
end

if isfield(options, 'logToDb') && ischar(options.logToDb)
    options.logToDb = strcmp(options.logToDb, 'true');
end

if isfield(options, 'runner')
    runnerClassName = options.runner;
else
    runnerClassName = 'SimulationRunner';
end

%run simulation
runner = edu.stanford.covert.cell.sim.runners.(runnerClassName)(options);
sim = runner.run();