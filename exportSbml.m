%Set warnings, load code
setWarnings();
setPath();

%load model
load('data\Simulation_fitted.mat');
sim = simulation;

%export model to SBML
edu.stanford.covert.cell.sim.util.ExportSbml.run(sim, 'output/sbml');