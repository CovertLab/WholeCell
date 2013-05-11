% Reindexes whole cell simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 6/30/2011
function runReindexing(outputDirectory)
setWarnings();
setPath();
setPreferences();

import edu.stanford.covert.cell.sim.util.DiskLogger;

DiskLogger.reindexTimeCourses(outputDirectory);