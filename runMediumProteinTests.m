% Runs the medium protein tests for the whole cell simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/13/2010
function runMediumProteinTests()
%% initialize
warning('off','WholeCell:warning');
setWarnings();
setPath();
setPreferences();

%% run tests
import edu.stanford.covert.test.XMLTestRunDisplay;
import edu.stanford.covert.test.runtests;
monitor = XMLTestRunDisplay('Medium Protein Tests','Whole cell simulation medium protein tests', 'output/runMediumProteinTests/tests.xml');
runtests(monitor, {
    'edu.stanford.covert.cell.sim.Protein_Test'
    });

