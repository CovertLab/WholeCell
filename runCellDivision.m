% Divide cell into 2 daughter cells.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/24/2012
function [sim, daughters] = runCellDivision(simBatch, simIdx, simIdxChild1, simIdxChild2)
import edu.stanford.covert.cell.sim.util.DiskLogger;
import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;

baseOutDir = [SimulationDiskUtil.getBaseDir() filesep simBatch];
simDir = [baseOutDir filesep num2str(simIdx)];

md = DiskLogger.loadMetadata(simDir);
sim = loadSimulation(simDir, md.lengthSec, md.lengthSec, 1, 'extract');

[daughters, daughterInitialConditions] = sim.divideState();
daughterInitialCondition1 = daughterInitialConditions(1);
daughterInitialCondition2 = daughterInitialConditions(2);

save([baseOutDir filesep num2str(simIdxChild1) filesep 'initialConditions.mat'], '-struct', 'daughterInitialCondition1');
save([baseOutDir filesep num2str(simIdxChild2) filesep 'initialConditions.mat'], '-struct', 'daughterInitialCondition2');