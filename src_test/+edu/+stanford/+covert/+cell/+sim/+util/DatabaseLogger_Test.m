%Database logger test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updaetd 2/3/2011
classdef DatabaseLogger_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = DatabaseLogger_Test(name)
            this = this@TestCase(name);
        end
        
        function setUp(this)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
            sim.applyOptions('lengthSec', 10);
            this.simulation = sim;
        end
    end
    
    methods
        function testLogging(this)
            import edu.stanford.covert.cell.sim.util.DatabaseLogger;
            
            dbConnectionParameters = getConfig();
            database = edu.stanford.covert.db.MySQLDatabase(dbConnectionParameters);
            knowledgeBaseWID = edu.stanford.covert.cell.kb.KnowledgeBaseUtil.selectLatestKnowledgeBase(database);
            
            sim = this.simulation;
            
            metadata = edu.stanford.covert.cell.sim.util.ConditionSet.parseConditionSet(sim, 'src_test/+edu/+stanford/+covert/+cell/+sim/+util/fixtures/runSimulations.xml').metadata;
            metadata.knowledgeBaseWID = knowledgeBaseWID;
            metadata.revision = edu.stanford.covert.util.revision;
                        
            downsampleStepSec = 2;
            logger = DatabaseLogger(database, downsampleStepSec, metadata);
            logger.addMetadata(metadata);
            sim.run(logger);
                        
            metadatafields = {'shortDescription'; 'longDescription';
                'email'; 'lastName'; 'firstName'; 'affiliation'; 'knowledgeBaseWID';
                'startTime'; 'endTime'; 'lengthSec'; 'revision'; 'differencesFromRevision';
                'userName'; 'hostName'; 'ipAddress';
                'outputDirectory'; 'segmentSizeStep'; 'downsampleStepSec'};
            assertEqual(sort(metadatafields), sort(fieldnames(logger.metadata)));
            
            %load simulation data back from database
            [states, metadata, options, parameters, fittedConstants, randStreamStates] = ...
                DatabaseLogger.load(sim, database, logger.simulationWID);
            
            %assert states
            nTimepoints = sim.lengthSec / downsampleStepSec;
            assertEqual(sort(cellfun(@(x) x(7:end), sim.stateWholeCellModelIDs, 'UniformOutput', false)), sort(fieldnames(states)));                    
            assertEqual([1 1 nTimepoints], size(states.Time.values));
            assertEqual(permute(downsampleStepSec:downsampleStepSec:sim.lengthSec, [1 3 2]), states.Time.values);
            assertEqual([numel(sim.getRandStreamStates().simulation) numel(states.Time.values)], size(randStreamStates.simulation));
            
            %assert metadata
            assertEqual(logger.metadata, metadata);
            
            %assert options
            assertEqual(sim.getOptions, options);
                       
            %assert parameters
            assertEqual(sort(fieldnames(sim.getParameters)), sort(fieldnames(parameters)));
            assertEqual(sort(fieldnames(sim.getParameters.states)), sort(fieldnames(parameters.states)));
            assertEqual(sort(fieldnames(sim.getParameters.processes)), sort(fieldnames(parameters.processes)));
            for i=1:numel(sim.states)
                id = sim.states{1}.wholeCellModelID(7:end);
                assertEqual(sort(fieldnames(sim.getParameters.states.(id))), sort(fieldnames(parameters.states.(id))));
                
                fields = fieldnames(parameters.states.(id));
                for j = 1:numel(fields)
                    switch class(parameters.states.(id).(fields{j}))
                        case {'char', 'cell', 'struct'}
                            assertEqual(sim.getParameters.states.(id).(fields{j}), parameters.states.(id).(fields{j}));
                        otherwise
                            assertElementsAlmostEqual(sim.getParameters.states.(id).(fields{j}), parameters.states.(id).(fields{j}), 'relative', 1e-4);
                    end
                end
            end
            for i=1:numel(sim.states)
                id = sim.states{1}.wholeCellModelID(7:end);
                assertEqual(sort(fieldnames(sim.getParameters.states.(id))), sort(fieldnames(parameters.states.(id))));
                
                fields = fieldnames(parameters.states.(id));
                for j = 1:numel(fields)
                    switch class(parameters.states.(id).(fields{j}))
                        case {'char', 'cell', 'struct'}
                            assertEqual(sim.getParameters.states.(id).(fields{j}), parameters.states.(id).(fields{j}));
                        otherwise
                            assertElementsAlmostEqual(sim.getParameters.states.(id).(fields{j}), parameters.states.(id).(fields{j}), 'relative', 1e-4);
                    end
                end
            end
            
            %assert fitted constants
            assertEqual(sort(fieldnames(sim.getFittedConstants)), sort(fieldnames(fittedConstants)));
            assertEqual(sort(fieldnames(sim.getFittedConstants.states)), sort(fieldnames(fittedConstants.states)));
            assertEqual(sort(fieldnames(sim.getFittedConstants.processes)), sort(fieldnames(fittedConstants.processes)));
            for i=1:numel(sim.states)
                id = sim.states{1}.wholeCellModelID(7:end);
                assertEqual(sort(fieldnames(sim.getFittedConstants.states.(id))), sort(fieldnames(fittedConstants.states.(id))));
                
                fields = fieldnames(fittedConstants.states.(id));
                for j = 1:numel(fields)
                    switch class(fittedConstants.states.(id).(fields{j}))
                        case {'char', 'cell', 'struct'}
                            assertEqual(sim.getFittedConstants.states.(id).(fields{j}), fittedConstants.states.(id).(fields{j}));
                        otherwise
                            assertElementsAlmostEqual(sim.getFittedConstants.states.(id).(fields{j}), fittedConstants.states.(id).(fields{j}), 'relative', 1e-4);
                    end
                end
            end
            for i=1:numel(sim.states)
                id = sim.states{1}.wholeCellModelID(7:end);
                assertEqual(sort(fieldnames(sim.getFittedConstants.states.(id))), sort(fieldnames(fittedConstants.states.(id))));
                
                fields = fieldnames(fittedConstants.states.(id));
                for j = 1:numel(fields)
                    switch class(fittedConstants.states.(id).(fields{j}))
                        case {'char', 'cell', 'struct'}
                            assertEqual(sim.getFittedConstants.states.(id).(fields{j}), fittedConstants.states.(id).(fields{j}));
                        otherwise
                            assertElementsAlmostEqual(sim.getFittedConstants.states.(id).(fields{j}), fittedConstants.states.(id).(fields{j}), 'relative', 1e-4);
                    end
                end
            end
            
            %test loadSimulation
            [sim2, states, ~, ~, ~, metadata] = loadSimulation(logger.simulationWID);
            metadata.downsampleStepSec = logger.metadata.downsampleStepSec;
            
            assertEqual(options, sim2.getOptions());
            assertStructElementsAlmostEqual(parameters, sim2.getParameters(), 'relative', 1e-4, 0);
            assertStructElementsAlmostEqual(fittedConstants, sim2.getFittedConstants(), 'relative', 1e-8, 0);
            assertEqual(metadata, logger.metadata);
            
            states2 = struct;
            for i = 1:numel(sim2.states)
                state = sim2.states{i};
                stateID = state.wholeCellModelID(7:end);
                states2.(stateID) = struct;
                for j = 1:numel(state.stateNames)
                    name = state.stateNames{j};
                    states2.(stateID).(name) = state.(name);
                end
            end
            assertEqual(states, states2);
            
            %cleanup -- remove simulation from database
            database.prepareStatement('CALL delete_simulation("{Si}")', logger.simulationWID);
            database.query();
            database.close();
        end
        
        function testEarlyTermination(this)
            import edu.stanford.covert.cell.sim.util.DatabaseLogger;
            
            %setup
            dbConnectionParameters = getConfig();
            database = edu.stanford.covert.db.MySQLDatabase(dbConnectionParameters);
            knowledgeBaseWID = edu.stanford.covert.cell.kb.KnowledgeBaseUtil.selectLatestKnowledgeBase(database);
            
            sim = this.simulation;
            
            metadata = edu.stanford.covert.cell.sim.util.ConditionSet.parseConditionSet(sim, 'src_test/+edu/+stanford/+covert/+cell/+sim/+util/fixtures/runSimulations.xml').metadata;
            metadata.knowledgeBaseWID = knowledgeBaseWID;
                        
            downsampleStepSec = 2;
            logger = DatabaseLogger(database, downsampleStepSec, metadata);
            logger.addMetadata(metadata);
            
            %run simulation
            logger.initialize(sim);
            
            for i = 1:sim.getNumSteps / 2
                sim.evolveState();
                logger.append(sim);
            end
            
            logger.finalize(sim);
            
            %test
            [states, metadata] = logger.load(sim, database, logger.simulationWID);            
            assertEqual(sim.getNumSteps / 2, metadata.lengthSec);
            assertEqual(permute([2 4 5], [1 3 2]), states.Time.values);
            
            %cleanup -- remove simulation from database
            database.prepareStatement('CALL delete_simulation("{Si}")', logger.simulationWID);
            database.query();
            database.close();
        end
    end
end