% Load whole cell simulation from disk/database.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/12/2011
function [simulation, states, options, parameters, fittedConstants, metadata] = loadSimulation(directoryOrWID, timeInit, timeFin, downsampleStepSec, downsampleType)
%initialize
setWarnings();
setPath();
setPreferences();

%import
import edu.stanford.covert.cell.sim.Simulation;
import edu.stanford.covert.cell.sim.util.DatabaseLogger;
import edu.stanford.covert.cell.sim.util.DiskLogger;
import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
import edu.stanford.covert.db.MySQLDatabase;

%load metadata
if ischar(directoryOrWID)
    directoryOrWID = SimulationDiskUtil.getSimulation(directoryOrWID);
    metadata = DiskLogger.loadMetadata(directoryOrWID);
else
    if ~exist('database', 'var'); database = MySQLDatabase(config); end
    metadata = DatabaseLogger.loadMetadata(database, directoryOrWID);
end

%construct simulation object
simulation = CachedSimulationObjectUtil.load(metadata.revision);

%load options, parameters, fitted constants based on data stored in disk/database
if ischar(directoryOrWID)
    if nargin < 2
        timeInit = [];
    end
    if nargin < 3
        timeFin = [];
    end
    if nargin < 4
        downsampleStepSec = [];
    end
    if nargin < 5
        downsampleType = [];
    end
    [states, metadata, options, parameters, fittedConstants, randStreamStates] = ...
        DiskLogger.load(directoryOrWID, '-independent', timeInit, timeFin, downsampleStepSec, downsampleType);
else
    if ~exist('database', 'var')
        database = MySQLDatabase(config);
    end
    if nargin > 1
        warning('WholeCell:warning', 'Additional options ignored');
    end
    [states, metadata, options, parameters, fittedConstants, randStreamStates] = ...
        DatabaseLogger.load(simulation, database, directoryOrWID);
end

%apply options, parameters, constants, time courses to simulation
simulation.applyOptions(options);
simulation.applyParameters(parameters);
simulation.applyFittedConstants(fittedConstants);
simulation.applyRandStreamStates(randStreamStates);
simulation.allocateMemoryForState(numel(states.Time.values));

for i = 1:numel(simulation.states)
    state = simulation.states{i};
    stateID = state.wholeCellModelID(7:end);
    for j = 1:numel(state.stateNames)
        name = state.stateNames{j};
        try
            state.(name) = states.(stateID).(name);
        catch %#ok<CTCH>
            warning('WholeCell:warning', 'Data not provided for %s %s state', state.name, name);
        end
    end
end

for i = 1:numel(simulation.processes)
    process = simulation.processes{i};
    process.copyFromState();
end

simulation.applyPerturbationsToConstants();

%cleanup 
if exist('database', 'var')
    database.close();
end