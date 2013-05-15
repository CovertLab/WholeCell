% Generates test fixtures for whole cell simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/31/2011
function generateTestFixtures(useCachedKb)
import edu.stanford.covert.cell.sim.runners.SimulationRunner;

if nargin < 1
	useCachedKb = true;
end

%create runner which helps do the tasks below
runner = SimulationRunner(...
    'useCachedKb', useCachedKb, ...
    'useCachedSim', false, ...
    'cacheKb', true, ...
    'cacheSimulation', true);

%construction KB and simulation
sim = runner.constructKbAndSimulation();

% verify initial growth rate distribution
if ~runner.isSimulationFitted(sim)
    throw(MException('Simulation:notFitted', 'Simulation not properly fitted'));
end

%save simulation parameters in JSON format
runner.cacheSimulationParametersJson(sim);

%generate simulation test fixtures
runner.cacheSimulationTestFixtures(sim);