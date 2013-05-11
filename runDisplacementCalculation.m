% Calculates DNA-bound protein displacements
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 4/6/2012
function runDisplacementCalculation(simBatch, simIdx, outDir)
setWarnings;
setPath;

%import classes
import edu.stanford.covert.cell.sim.analysis.Figures34;
import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;

%parse inputs
if ischar(simIdx)
	simIdx = str2num(simIdx);
end

%run calculation
simBatchDir = [SimulationDiskUtil.getBaseDir() filesep simBatch filesep];
Figures34.calculateDNABoundProteinDisplacements(outDir, simBatchDir, simIdx, 1, 1);