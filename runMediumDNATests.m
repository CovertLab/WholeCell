% Runs the medium DNA tests for the whole cell simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/13/2010
function runMediumDNATests()
%% initialize
warning('off','WholeCell:warning');
setWarnings();
setPath();
setPreferences();

%% run tests
import edu.stanford.covert.test.XMLTestRunDisplay;
import edu.stanford.covert.test.runtests;
monitor = XMLTestRunDisplay('Medium DNA Tests','Whole cell simulation medium DNA tests', 'output/runMediumDNATests/tests.xml');
runtests(monitor, {
    'edu.stanford.covert.cell.sim.DNA_Test'
    });

