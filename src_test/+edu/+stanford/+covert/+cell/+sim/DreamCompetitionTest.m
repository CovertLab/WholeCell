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
            
            sim.applyAllParameters(struct('processes', struct('Metabolism', struct('unaccountedEnergyConsumption', 1000))));
            assertEqual(1000, sim.getAllParameters().processes.Metabolism.unaccountedEnergyConsumption);
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
                'simPath', 'output/dream-sim-1.mat', ...
                'verbosity', 0 ...
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
                'simPath', 'output/dream-sim-2.mat' ...
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
            condition.fittedConstants = sim.getFittedConstants();
            %condition.fixedConstants = sim.getFixedConstants();
            
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
                
                fittedConstants = s.getFittedConstants();
                if numel(fields(fittedConstants)) > 0
                    assertEqual(condition.fittedConstants.states.(id), data.fittedConstants.states.(id));
                else
                    assertEqual(struct(), condition.fittedConstants.states.(id));
                    assertFalse(isfield(data.fittedConstants.states, id));
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
                
                fittedConstants = p.getFittedConstants();
                if numel(fields(fittedConstants)) > 0
                    assertEqual(condition.fittedConstants.processes.(id), data.fittedConstants.processes.(id));
                else
                    assertEqual(struct(), condition.fittedConstants.processes.(id));
                    assertFalse(isfield(data.fittedConstants.processes, id));
                end
            end
            
            %simulate using XML
            simulateHighthroughputExperiments(...
                'seed', 1, ...
                'parameterValsPath', parameterValsPath, ...
                'simPath', 'output/dream-sim-3.mat' ...
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
                'simPath', 'output/dream-sim-4.mat' ...
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
                'simPath', 'output/dream-sim-5.mat' ...
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
                'simPath', 'output/dream-sim-6.mat' ...
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
                    'simPath', sprintf('output/dream-sim-batch-%d.mat', i) ...
                    );
            end
            
            %average
            averageHighthroughputExperiments(...
                'simPathPattern', 'output/dream-sim-batch-*.mat', ...
                'avgValsPath', 'output/dream-sim-avg.mat' ...
                );
            
            %assert
            assertEqual(2, exist('output/dream-sim-avg.mat', 'file'));
        end
        
        function test_averageHighthroughputExperimentsAndCalcErrors(~)
            %% simulate
            sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
            mass = sim.state('Mass');
            
            parameterVals = sim.getAllParameters(); %get default parameters
            nSimulations = 2;
            for i = 1:nSimulations
                parameterVals.lengthSec = i * 10; %override defaults
                simulateHighthroughputExperiments(...
                    'seed', i, ...
                    'parameterVals', parameterVals, ...
                    'simPath', sprintf('output/dream-sim-batch-%d.mat', i) ...
                    );
                
                % validate sim
                simVals = load(sprintf('output/dream-sim-batch-%d.mat', i));
                
                assertEqual(0:i*10, simVals.time);
                assertEqual([1 i*10+1], size(simVals.growth))
                assertEqual([1 i*10+1], size(simVals.mass))
                assertEqual([1 i*10+1], size(simVals.volume))
                assertTrue(all(isreal(simVals.growth) & isfinite(simVals.growth)))
                assertElementsAlmostEqual(mass.cellInitialDryWeight / (1 - mass.fractionWetWeight) * ones(1, i*10+1), simVals.mass);
                assertTrue(all(isreal(simVals.volume) & isfinite(simVals.volume)))
                
                assertTrue(isscalar(simVals.repInitTime))
                assertTrue(isscalar(simVals.repTermTime))
                assertTrue(isscalar(simVals.cellCycleLen))
                
                assertTrue(all(isreal(simVals.metConcs.mean) & isfinite(simVals.metConcs.mean) & simVals.metConcs.mean >= 0))
                assertTrue(all(isreal(simVals.metConcs.std) & isfinite(simVals.metConcs.std) & simVals.metConcs.std >= 0))
                
                assertTrue(all(isreal(simVals.dnaSeq.mean) & isfinite(simVals.dnaSeq.mean) & simVals.dnaSeq.mean >= 0))
                assertTrue(all(isreal(simVals.dnaSeq.std) & isfinite(simVals.dnaSeq.std) & simVals.dnaSeq.std >= 0))
                assertEqual(2, min(simVals.dnaSeq.mean))
                assertEqual(0, min(simVals.dnaSeq.std))
                
                assertTrue(all(isreal(simVals.rnaSeq.mean) & isfinite(simVals.rnaSeq.mean) & simVals.rnaSeq.mean >= 0))
                assertTrue(all(isreal(simVals.rnaSeq.std) & isfinite(simVals.rnaSeq.std) & simVals.rnaSeq.std >= 0))
                
                assertTrue(all(all(isreal(simVals.chipSeq.mean) & isfinite(simVals.chipSeq.mean) & simVals.chipSeq.mean >= 0)))
                assertTrue(all(all(isreal(simVals.chipSeq.std) & isfinite(simVals.chipSeq.std) & simVals.chipSeq.std >= 0)))
                
                assertTrue(all(isreal(simVals.rnaArray.mean) & isfinite(simVals.rnaArray.mean) & simVals.rnaArray.mean >= 0))
                assertTrue(all(isreal(simVals.rnaArray.std) & isfinite(simVals.rnaArray.std) & simVals.rnaArray.std >= 0))
                
                assertTrue(all(isreal(simVals.protArray.mean) & isfinite(simVals.protArray.mean) & simVals.protArray.mean >= 0))
                assertTrue(all(isreal(simVals.protArray.std) & isfinite(simVals.protArray.std) & simVals.protArray.std >= 0))
                
                assertTrue(all(isreal(simVals.rxnFluxes.mean) & isfinite(simVals.rxnFluxes.mean)))
                assertTrue(all(isreal(simVals.rxnFluxes.std) & isfinite(simVals.rxnFluxes.std) & simVals.rxnFluxes.std >= 0))
            end
            
            %% average
            refParameterVals = parameterVals;
            refAvgVals = averageHighthroughputExperiments('simPathPattern', 'output/dream-sim-batch-*.mat');
            
            assertEqual([nSimulations 21], size(refAvgVals.growth))
            assertEqual([nSimulations 21], size(refAvgVals.mass))
            assertEqual([nSimulations 21], size(refAvgVals.volume))
            assertEqual([nSimulations 1], size(refAvgVals.repInitTime))
            assertEqual([nSimulations 1], size(refAvgVals.repTermTime))
            assertEqual([nSimulations 1], size(refAvgVals.cellCycleLen))
            assertAllEqual(false, isnan(refAvgVals.growth(nSimulations, :)))
            assertAllEqual(false, isnan(refAvgVals.mass(nSimulations, :)))
            assertAllEqual(false, isnan(refAvgVals.volume(nSimulations, :)))
            assertAllEqual(NaN, refAvgVals.repTermTime)
            assertAllEqual(NaN, refAvgVals.cellCycleLen)
            
            propNames = {'metConcs'; 'dnaSeq'; 'rnaSeq'; 'chipSeq'; 'rnaArray'; 'protArray'; 'rxnFluxes'};
            for i = 1:numel(propNames)
                assertTrue(all(all(isreal(refAvgVals.(propNames{i}).mean))));
                assertTrue(all(all(isfinite(refAvgVals.(propNames{i}).mean))));
                if ~strcmp(propNames{i}, 'rxnFluxes')
                    assertTrue(all(all(refAvgVals.(propNames{i}).mean >= 0)));
                end
                
                assertTrue(all(all(isreal(refAvgVals.(propNames{i}).std))));
                assertTrue(all(all(isfinite(refAvgVals.(propNames{i}).std))));
                assertTrue(all(all(refAvgVals.(propNames{i}).std >= 0)));
            end
            
            %% calculate distances from "reference"
            [dists, avgVals] = averageHighthroughputExperimentsAndCalcErrors(...
                'parameterVals', parameterVals, ...
                'simPathPattern', 'output/dream-sim-batch-*.mat', ...
                'avgValsPath', 'output/dream-sim-avg.mat', ...
                'refParameterVals', refParameterVals, ...
                'refAvgVals', refAvgVals, ...
                'verbosity', 0 ...
                );
            
            %assertions            
            assertEqual(refAvgVals, avgVals);
            assertEqual(0, dists.parameter);
            assertEqual(0, dists.prediction);
        end
        
        function test_calcParameterAndPredictionScoring(~)
            refParameterValsPath = 'output/1_1.parameters.mat';
            refAvgValsPath = 'output/1_1.predictions.mat';
            
            sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
            parameterVals = sim.getAllParameters(); %get default parameters
            parameterVals.lengthSec = 10;
            for j = 1:3
                %save parameters
                parameterValsPath = sprintf('output/%d_1.parameters.mat', j);
                save(parameterValsPath, '-struct', 'parameterVals');
                
                for i = 1:2
                    %fprintf('Running simulation %d for batch %d\n', i, j);
                    
                    %simulate
                    simulateHighthroughputExperiments(...
                        'seed', j * 100 + i, ...
                        'parameterValsPath', parameterValsPath, ...
                        'simPath', sprintf('output/%d_1.predictions_%d.mat', j, i), ...
                        'verbosity', 0 ...
                        );
                end
                
                % calculate "reference"
                if j == 1
                    averageHighthroughputExperiments(...
                        'simPathPattern', sprintf('output/%d_1.predictions_*.mat', j), ...
                        'avgValsPath', sprintf('output/%d_1.predictions.mat', j), ...
                        'verbosity', 0 ...
                        );
                end
                
                %average
                averageHighthroughputExperimentsAndCalcErrors(...
                    'parameterValsPath', parameterValsPath, ...
                    'simPathPattern', sprintf('output/%d_1.predictions_*.mat', j), ...
                    'avgValsPath', sprintf('output/%d_1.predictions.mat', j), ...
                    'distsPath', sprintf('output/%d_1.distances.mat', j), ...
                    'refParameterValsPath', refParameterValsPath, ...
                    'refAvgValsPath', refAvgValsPath, ...
                    'verbosity', 0 ...
                    );
            end
            
            %calculate scores
            [submissionList, distances, pValues, scores, ranks] = calcParameterAndPredictionScoring(...
                'refParameterValsPath', refParameterValsPath, ...
                'refAvgValsPath', refAvgValsPath, ...
                'submissionsPathBase', 'output', ...
                'nullDistanceDistribSize', 20 ...
                );
            
            assertEqual([
                struct('user', '1', 'parameterSetTimestamp', '1')
                struct('user', '2', 'parameterSetTimestamp', '1')
                struct('user', '3', 'parameterSetTimestamp', '1')
                ], submissionList);
            
            assertEqual([3 1], size(distances))
            assertEqual([3 1], size(pValues))
            assertEqual([3 1], size(scores))
            assertEqual([3 1], size(ranks))
            
            assertAllEqual(0, [distances.parameter])
            assertAllEqual(1/20, [pValues.parameter])
            
            assertEqual(0, distances(1).prediction)
            assertEqual(1/20, pValues(1).prediction)
            assertAllEqual(true, [distances(2:end).prediction] > 0)
            assertAllEqual(true, [pValues(2:end).prediction] > 0)
            assertAllEqual(true, scores(2:end) < scores(1))
            
            assertEqual(1, ranks(1))
        end
    end
end