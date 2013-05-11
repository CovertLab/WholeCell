% Runs the medium tests for the whole cell simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/13/2010
function runMediumTests()
%% initialize
warning('off','WholeCell:warning');
setWarnings();
setPath();
setPreferences();

%% run tests
import edu.stanford.covert.test.XMLTestRunDisplay;
import edu.stanford.covert.test.runtests;
monitor = XMLTestRunDisplay('Medium Tests','Whole cell simulation medium tests', 'output/runMediumTests/tests.xml');
runtests(monitor, {
    'edu.stanford.covert.cell.sim.CellCycle_Test'
    'edu.stanford.covert.cell.sim.Cytokinesis_Test'
    'edu.stanford.covert.cell.sim.DNADamageRepair_Test'
    'edu.stanford.covert.cell.sim.ProteinGrowth_Test'
    'edu.stanford.covert.cell.sim.RNA_Test'
    'edu.stanford.covert.cell.sim.Simulation_Integrated_Test'
    'edu.stanford.covert.cell.sim.Simulation_KnowledgeBase_Test'
    'edu.stanford.covert.cell.sim.Simulation_Test'
    'edu.stanford.covert.cell.sim.SimulationStateSideEffect_Test'    
    });

