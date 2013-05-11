%Disk logger test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updaetd 2/3/2011
classdef DiskLogger_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = DiskLogger_Test(name)
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
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            sim = this.simulation;
            pm = sim.state('ProteinMonomer');
            
            metadata = edu.stanford.covert.cell.sim.util.ConditionSet.parseConditionSet(...
                sim, 'src_test/+edu/+stanford/+covert/+cell/+sim/+util/fixtures/runSimulations.xml').metadata;
            metadata.knowledgeBaseWID = 145;
            metadata.revision = edu.stanford.covert.util.revision;
            
            outputDirectory = 'output/runSmallTests/DiskLogger';
            delete([outputDirectory '/*']);
            segmentSizeSec = 2;
            logger = DiskLogger(outputDirectory, segmentSizeSec, metadata);
            logger.addMetadata(metadata);
            sim.run(logger);
            
            DiskLogger.reindexTimeCourses(outputDirectory);
            
            metadatafields = {'shortDescription'; 'longDescription';
                'email'; 'lastName'; 'firstName'; 'affiliation'; 'knowledgeBaseWID';
                'startTime'; 'endTime'; 'lengthSec'; 'revision'; 'differencesFromRevision';
                'userName'; 'hostName'; 'ipAddress';
                'outputDirectory'; 'segmentSizeStep'; 'downsampleStepSec';
                'stateNames'; 'dependentStateNames'};
            assertEqual(sort(metadatafields), sort(fieldnames(logger.metadata)));
            
            %% assert created files
            files = dir(outputDirectory);
            fileNames = {'.'; '..'; 'fittedConstants.mat'; 'metadata.mat'; 'options.mat'; 'parameters.mat';
                'state-0.mat'; 'state-1.mat'; 'state-2.mat'; 'state-3.mat'; 'state-4.mat'; 'state-5.mat'; 
                'randStreamStates.mat'};
            assertTrue(all(ismember(fileNames, {files.name}')));
            
            %% assert loading back from disk
            downsampleStepSec = 2;
            nTimepoints = sim.lengthSec / downsampleStepSec;
            
            %extraction downsampling
            [states, metadata, options, parameters, fittedConstants, randStreamStates] = ...
                DiskLogger.load(outputDirectory, '-all', [], [], downsampleStepSec, 'extract');
            metadata.downsampleStepSec = logger.metadata.downsampleStepSec;
            metadata.segmentSizeStep = logger.metadata.segmentSizeStep;
            
            assertEqual(sort(cellfun(@(x) x(7:end), sim.stateWholeCellModelIDs, 'UniformOutput', false)), sort(fieldnames(states)));
            assertEqual([1 1 nTimepoints], size(states.Time.values));
            assertEqual(permute([1:downsampleStepSec:sim.lengthSec-3 sim.lengthSec], [1 3 2]), states.Time.values);
            assertEqual([numel(sim.getRandStreamStates().simulation) sim.lengthSec+1], size(randStreamStates.simulation));
            
            assertEqual(logger.metadata, metadata);
            assertEqual(sim.getOptions, options);
            assertEqual(sim.getParameters, parameters);
            assertEqual(sim.getFittedConstants, fittedConstants);
            
            pmCounts = states.ProteinMonomer.counts;
            
            states = DiskLogger.load(outputDirectory, {'ProteinMonomer' 'counts' '-sum' ':'}, [], [], downsampleStepSec, 'extract');
            assertEqual(sum(pmCounts, 1), states.ProteinMonomer.counts);
            
            states = DiskLogger.load(outputDirectory, {'ProteinMonomer' 'counts' '-sum' '-sum'}, [], [], downsampleStepSec, 'extract');
            assertEqual(sum(sum(pmCounts, 1), 2), states.ProteinMonomer.counts);
            
            states = DiskLogger.load(outputDirectory, {'ProteinMonomer' 'counts' ':' [3 4]}, [], [], downsampleStepSec, 'extract');
            assertEqual(sum(pmCounts(:, [3 4], :), 2), states.ProteinMonomer.counts);
            
            states = DiskLogger.load(outputDirectory, {'ProteinMonomer' 'counts' [pm.matureIndexs pm.boundIndexs] [1 4]}, [], [], downsampleStepSec, 'extract');
            assertEqual(sum(pmCounts(pm.matureIndexs, [1 4], :) + pmCounts(pm.boundIndexs, [1 4], :), 2), states.ProteinMonomer.counts);
            
            %mean downsampling
            [states, metadata, options, parameters, fittedConstants] = ...
                DiskLogger.load(outputDirectory, '-all', [], [], downsampleStepSec, 'mean');
            metadata.downsampleStepSec = logger.metadata.downsampleStepSec;
            
            assertEqual(sort(cellfun(@(x) x(7:end), sim.stateWholeCellModelIDs, 'UniformOutput', false)), sort(fieldnames(states)));
            assertEqual([1 1 nTimepoints], size(states.Time.values));
            assertEqual(permute(0.5+(1:downsampleStepSec:sim.lengthSec), [1 3 2]), states.Time.values);
            
            assertEqual(logger.metadata, metadata);
            assertEqual(sim.getOptions, options);
            assertEqual(sim.getParameters, parameters);
            assertEqual(sim.getFittedConstants, fittedConstants);
            
            states = DiskLogger.load(outputDirectory, '-independent', [], [], downsampleStepSec, 'extract'); %#ok<NASGU>
            states = DiskLogger.load(outputDirectory, '-dependent', [], [], downsampleStepSec, 'extract'); %#ok<NASGU>
            
            %test loadSimulation
            [sim2, states, ~, ~, ~, metadata] = loadSimulation([pwd filesep outputDirectory]);
            metadata.downsampleStepSec = logger.metadata.downsampleStepSec;
            
            assertEqual(options, sim2.getOptions());
            assertEqual(parameters, sim2.getParameters());
            assertEqual(fittedConstants, sim2.getFittedConstants());
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
        end
        
        function testEarlyTermination(this)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            %setup
            sim = this.simulation;
            
            metadata = edu.stanford.covert.cell.sim.util.ConditionSet.parseConditionSet(sim, 'src_test/+edu/+stanford/+covert/+cell/+sim/+util/fixtures/runSimulations.xml').metadata;
            metadata.knowledgeBaseWID = 145;
            
            outputDirectory = 'output/runSmallTests/DiskLogger';
            delete([outputDirectory '/*']);
            segmentSizeSec = 2;
            logger = DiskLogger(outputDirectory, segmentSizeSec, metadata);
            logger.addMetadata(metadata);
            
            %run simulation
            logger.initialize(sim);
            
            for i = 1:sim.getNumSteps / 2
                sim.evolveState();
                logger.append(sim);
            end
            
            logger.finalize(sim);
            
            %test
            files = dir(outputDirectory);
            fileNames = {'.'; '..'; 'fittedConstants.mat'; 'metadata.mat'; 'options.mat'; 'parameters.mat';
                'state-0.mat'; 'state-1.mat'; 'state-2.mat'; 'state-3.mat'};
            assertTrue(all(ismember(fileNames, {files.name}')));
            
            [states, metadata] = logger.load(outputDirectory, {'Time' 'values'}, [], [], 2, 'extract');
            assertEqual(sim.getNumSteps / 2, metadata.lengthSec);
            assertEqual(permute([1 3 5], [1 3 2]), states.Time.values);
            
            [states, metadata] = logger.load(outputDirectory, {'Time' 'values'}, [], [], 2, 'mean');
            assertEqual(sim.getNumSteps / 2, metadata.lengthSec);
            assertEqual(permute([1.5 3.5 5], [1 3 2]), states.Time.values);
        end
    end
end