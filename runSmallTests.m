% Runs the small tests for the whole cell simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/13/2010
function runSmallTests()
%% initialize
warning('off', 'WholeCell:warning');
setWarnings();
setPath();
setPreferences();

%% run tests
import edu.stanford.covert.test.XMLTestRunDisplay;
import edu.stanford.covert.test.runtests;
monitor = XMLTestRunDisplay('Small Tests', 'Whole cell simulation small tests', 'output/runSmallTests/tests.xml');
runtests(monitor, {
    'edu.stanford.covert.cell.sim.constant'
    'edu.stanford.covert.cell.sim.process'
    'edu.stanford.covert.cell.sim.state'
    'edu.stanford.covert.cell.sim.util'
    'edu.stanford.covert.cell.kb'
    'edu.stanford.covert.db'
    'edu.stanford.covert.io'
    'edu.stanford.covert.test'
    'edu.stanford.covert.util'
    });
