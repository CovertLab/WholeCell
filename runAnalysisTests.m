% Runs the analysis tests for the whole cell simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 3/23/2011
function runAnalysisTests()
%% initialize
warning('off','WholeCell:warning');
setWarnings();
setPath();
setPreferences();

%% run tests
import edu.stanford.covert.test.XMLTestRunDisplay;
import edu.stanford.covert.test.runtests;
monitor = XMLTestRunDisplay('Analysis Tests','Whole cell simulation analysis tests', 'output/runAnalysisTests/tests.xml');
runtests(monitor, {
    'edu.stanford.covert.cell.sim.analysis'
    });