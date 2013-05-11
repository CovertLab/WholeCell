% Runs the large tests for the whole cell simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/13/2010
function runLargeTests()
%% initialize
warning('off','WholeCell:warning');
setWarnings();
setPath();
setPreferences();

%% run tests
import edu.stanford.covert.test.XMLTestRunDisplay;
import edu.stanford.covert.test.runtests;
monitor = XMLTestRunDisplay('Large Tests','Whole cell simulation large tests', 'output/runLargeTests/tests.xml');
runtests(monitor, {
    'edu.stanford.covert.cell.sim.Simulation_Large_Test'    
    });
