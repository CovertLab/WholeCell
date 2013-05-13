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
    
    %tests
    methods
        %Run simulation using struct of parameter values
        function test_runHighthroughputExperimentsSimulation1(~)
            sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
            
            parameterVals = sim.getAllParameters(); %get default parameters
            parameterVals.lengthSec = 10;          %override defaults
            
            runHighthroughputExperimentsSimulation(...
                'seed', 1, ...
                'parameterVals', parameterVals, ...
                'outPath', 'output/dream-sim-1.mat' ...
                );
            
            experimentalData = load('output/dream-sim-1.mat');
            assertEqual(0:10, experimentalData.time);
        end
        
        %Run simulation using mat file of parameter values
        function test_runHighthroughputExperimentsSimulation2(~)
            sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
            
            parameterVals = sim.getAllParameters(); %get default parameters
            parameterVals.lengthSec = 5;          %override defaults
            
            parameterValsPath = 'output/dream-sim-2-parameters.mat';
            save(parameterValsPath, '-struct', 'parameterVals');
            
            runHighthroughputExperimentsSimulation(...
                'seed', 1, ...
                'parameterValsPath', parameterValsPath, ...
                'outPath', 'output/dream-sim-2.mat' ...
                );
            
            experimentalData = load('output/dream-sim-2.mat');
            assertEqual(0:5, experimentalData.time);
        end
    end
end