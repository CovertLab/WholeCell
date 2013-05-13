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
        
        %Run simulation using mat file of parameter values
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
        
        function test_averageHighthroughputExperiments(~)
            %simulate
            sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
            parameterVals = sim.getAllParameters(); %get default parameters
            for i = 1:2
                parameterVals.lengthSec = i * 10; %override defaults
                simulateHighthroughputExperiments(...
                    'seed', i, ...
                    'parameterVals', parameterVals, ...
                    'outPath', sprintf('output/dream-sim-%d.mat', i) ...
                    );
            end
            
            %average
            averageHighthroughputExperiments(...
                'inPathPattern', 'output/dream-sim-*.mat', ...
                'outPath', 'output/dream-sim-avg.mat' ...
                );
            
            %assert
            assertEqual(2, exist('output/dream-sim-avg.mat', 'file'));
        end
    end
end