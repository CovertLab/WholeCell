%Summary logger test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updaetd 2/3/2011
classdef SummaryLogger_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = SummaryLogger_Test(name)
            this = this@TestCase(name);
        end
        
        function setUp(this)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load();
            sim.applyOptions('lengthSec', 10);
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            this.simulation = sim;
            
            if ~exist('output/runSmallTests/SummaryLogger', 'dir')
                mkdir('output/runSmallTests/SummaryLogger');
            end
        end
    end
    
    methods
        function testLogging(this)
            sim = this.simulation;
            stepSizeSec = 2;
            logger = edu.stanford.covert.cell.sim.util.SummaryLogger(stepSizeSec, 2, 'output/runSmallTests/SummaryLogger');
            
            sim.run(logger);
            close(logger.figureHandle);
            
            nTimepoints = sim.lengthSec / stepSizeSec + 1;
            assertEqual(19, numel(fieldnames(logger.log)));
            assertEqual([1 nTimepoints], size(logger.log.runTime));
            assertEqual([1 nTimepoints], size(logger.log.time));
            assertEqual([1 nTimepoints], size(logger.log.mass));
            assertEqual([1 nTimepoints], size(logger.log.growth));
            assertEqual([13 nTimepoints], size(logger.log.metabolites));
            assertEqual([7 nTimepoints], size(logger.log.rnas));
            assertEqual([7 nTimepoints], size(logger.log.proteins));
        end
        
        function testEarlyTermination(this)
            sim = this.simulation;
            stepSizeSec = 2;
            logger = edu.stanford.covert.cell.sim.util.SummaryLogger(stepSizeSec, 2, 'output/runSmallTests/SummaryLogger');
            logger.initialize(sim);
            
            for i = 1:sim.getNumSteps / 2
                sim.evolveState();
                logger.append(sim);
            end
            
            logger.finalize(sim);
            
            fields = fieldnames(logger.log);
            for i = 1:numel(fields)
                assertEqual(ceil(sim.getNumSteps / 2 / stepSizeSec + 1), size(logger.log.(fields{i}), 2));
            end
            
            close(logger.figureHandle);
        end
    end
end