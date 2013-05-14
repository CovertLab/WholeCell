%DreamCompetitionTest
% Test new classes and functions created for DREAM parameter estimation
% competition.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
classdef DreamCompetitionTest < TestCase
    %constructor
    methods
        function this = DreamCompetitionTest(name)
            this = this@TestCase(name);
        end
    end
    
    %test getting, setting parameters
    methods
        function test_getApplyAllParameters(~)
            %import
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            
            %load simulation object
            sim = CachedSimulationObjectUtil.load();
            
            %set, get parameter values
            sim.applyAllParameters(struct('lengthSec', 1));
            assertEqual(1, sim.getAllParameters().lengthSec);
            
            sim.applyAllParameters(struct('processes', struct('Metabolism', struct('cellCycleLength', 1000))));
            assertEqual(1000, sim.getAllParameters().processes.Metabolism.cellCycleLength);
        end
        
        function test_getMetabolicReactionKinetics(~)
            %import
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            
            %load simulation object
            sim = CachedSimulationObjectUtil.load();
            met = sim.process('Metabolism');
            
            %get kinetics
            kinetics = sim.getMetabolicReactionKinetics();
            assertEqual(met.reactionWholeCellModelIDs, fields(kinetics));
            assertEqual(met.enzymeBounds(met.reactionIndexs_fba, :), ...
                met.fbaEnzymeBounds(met.fbaReactionIndexs_metabolicConversion, :));
        end
        
        function test_setMetabolicReactionKinetics(~)
            %import
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            
            %load simulation object
            sim = CachedSimulationObjectUtil.load();
            met = sim.process('Metabolism');
            
            %set kinetics
            rxnId = met.reactionWholeCellModelIDs{1};
            kinetics.(rxnId).for = 1;
            sim.setMetabolicReactionKinetics(kinetics);
            assertEqual(kinetics.(rxnId).for, sim.getMetabolicReactionKinetics().(rxnId).for);
            assertEqual(met.enzymeBounds(met.reactionIndexs_fba, :), ...
                met.fbaEnzymeBounds(met.fbaReactionIndexs_metabolicConversion, :));
        end
    end
    
    %test in silico experiment logger
    methods
        %Run simulation using struct of parameter values
        function test_simulateHighthroughputExperiments1(~)
            sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
            
            parameterVals = sim.getAllParameters(); %get default parameters
            parameterVals.lengthSec = 10;          %override defaults
            
            simulateHighthroughputExperiments(...
                'seed', 1, ...
                'parameterVals', parameterVals, ...
                'outPath', 'output/dream-sim-1.mat' ...
                );
            
            experimentalData = load('output/dream-sim-1.mat');
            assertEqual(0:10, experimentalData.time);
        end
        
        %Run simulation using .mat file of parameter values
        function test_simulateHighthroughputExperiments2(~)
            sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
            
            parameterVals = sim.getAllParameters(); %get default parameters
            parameterVals.lengthSec = 5;          %override defaults
            
            parameterValsPath = 'output/dream-sim-parameters-2.mat';
            save(parameterValsPath, '-struct', 'parameterVals');
            
            simulateHighthroughputExperiments(...
                'seed', 1, ...
                'parameterValsPath', parameterValsPath, ...
                'outPath', 'output/dream-sim-2.mat' ...
                );
            
            experimentalData = load('output/dream-sim-2.mat');
            assertEqual(0:5, experimentalData.time);
        end
        
        %Run simulation using XML file of parameter values
        function test_simulateHighthroughputExperiments3(~)
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.ConditionSet;
            import edu.stanford.covert.util.StructUtil;
            
            sim = CachedSimulationObjectUtil.load();
            sim.applyAllParameters('lengthSec', 5);
            
            parameterValsPath = 'output/dream-sim-parameters-3.xml';
            
            metadata = struct();
            metadata.firstName = 'Jonathan';
            metadata.lastName = 'Karr';
            metadata.email = 'jkarr@stanford.edu';
            metadata.affiliation = 'Stanford University';
            metadata.userName = 'jkarr';
            metadata.hostName = 'covertlab-jkarr.stanford.edu';
            metadata.ipAddress = '171.65.92.146';
            metadata.revision = 1;
            metadata.differencesFromRevision = '';
            
            condition.shortDescription = 'Test simulation';
            condition.longDescription = 'Test simulation';
            condition.replicates = 1;
            condition.options = sim.getOptions();
            condition.parameters = sim.getParameters();
            
            ConditionSet.generateConditionSet(sim, metadata, condition, parameterValsPath);
            
            %verify condition set
            data = ConditionSet.parseConditionSet(sim, parameterValsPath);
            tmp = StructUtil.catstruct(metadata, struct(...
                'shortDescription', condition.shortDescription, ...
                'longDescription', condition.longDescription ...
                ));
            assertEqual(tmp, data.metadata);
            
            assertEqual(condition.options, data.options);
            
            for i = 1:numel(sim.states)
                s = sim.states{i};
                id = s.wholeCellModelID(7:end);
                parameters = s.getParameters();
                if numel(fields(parameters)) > 0
                    assertEqual(condition.parameters.states.(id), data.parameters.states.(id));
                else
                    assertEqual(struct(), condition.parameters.states.(id));
                    assertFalse(isfield(data.parameters.states, id));
                end
            end
            
            for i = 1:numel(sim.processes)
                p = sim.processes{i};
                id = p.wholeCellModelID(9:end);
                parameters = p.getParameters();
                if numel(fields(parameters)) > 0
                    assertEqual(condition.parameters.processes.(id), data.parameters.processes.(id));
                else
                    assertEqual(struct(), condition.parameters.processes.(id));
                    assertFalse(isfield(data.parameters.processes, id));
                end
            end
            
            %simulate using XML
            simulateHighthroughputExperiments(...
                'seed', 1, ...
                'parameterValsPath', parameterValsPath, ...
                'outPath', 'output/dream-sim-3.mat' ...
                );
            
            experimentalData = load('output/dream-sim-3.mat');
            assertEqual(0:5, experimentalData.time);
        end
        
        %test knocking out
        function test_simulateHighthroughputExperiments4(~)
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            
            sim = CachedSimulationObjectUtil.load();
            sim.applyAllParameters(...
                'geneticKnockouts', {'MG_001'}, ...
                'lengthSec', 5 ...
                );
            
            simulateHighthroughputExperiments(...
                'parameterVals', sim.getAllParameters(), ...
                'outPath', 'output/dream-sim-4.mat' ...
                );
            
            experimentalData = load('output/dream-sim-4.mat');
            assertEqual(0:5, experimentalData.time);
        end
        
        %test setting initial conditions
        function test_simulateHighthroughputExperiments5(~)
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.util.StructUtil;
            
            sim = CachedSimulationObjectUtil.load();
            sim.applyAllParameters(...
                'lengthSec', 5 ...
                );
            
            initialConditions = sim.getTimeCourses();
            simulateHighthroughputExperiments(...
                'parameterVals', sim.getAllParameters(), ...
                'initialConditions', initialConditions, ...
                'outPath', 'output/dream-sim-5.mat' ...
                );
            
            experimentalData = load('output/dream-sim-5.mat');
            assertEqual(0:5, experimentalData.time);
        end
        
        %test setting initial conditions
        function test_simulateHighthroughputExperiments6(~)
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.util.StructUtil;
            
            sim = CachedSimulationObjectUtil.load();
            sim.applyAllParameters(...
                'lengthSec', 10 ...
                );
            
            initialConditionsPath = 'output/dream-sim-initial-conditions-6.mat';
            initialConditions = sim.getTimeCourses(); %#ok<NASGU>
            save(initialConditionsPath, '-struct', 'initialConditions');
            
            simulateHighthroughputExperiments(...
                'parameterVals', sim.getAllParameters(), ...
                'initialConditionsPath', initialConditionsPath, ...
                'outPath', 'output/dream-sim-6.mat' ...
                );
            
            experimentalData = load('output/dream-sim-6.mat');
            assertEqual(0:10, experimentalData.time);
        end
        
        function test_averageHighthroughputExperiments(~)
            %simulate
            sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
            parameterVals = sim.getAllParameters(); %get default parameters
            for i = 1:2
                parameterVals.lengthSec = i * 10; %override defaults
                simulateHighthroughputExperiments(...
                    'seed', i, ...
                    'parameterVals', parameterVals, ...
                    'outPath', sprintf('output/dream-sim-batch-%d.mat', i) ...
                    );
            end
            
            %average
            averageHighthroughputExperiments(...
                'inPathPattern', 'output/dream-sim-batch-*.mat', ...
                'outPath', 'output/dream-sim-avg.mat' ...
                );
            
            %assert
            assertEqual(2, exist('output/dream-sim-avg.mat', 'file'));
        end
    end
end