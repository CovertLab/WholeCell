%SimulationStateSideEffect_Test
% SimulationStateSideEffect, SimulationStateSideEffectItem test class
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/17/2010
classdef SimulationStateSideEffect_Test < TestCase
    properties
        simulation
    end
    
    methods 
        function this = SimulationStateSideEffect_Test(name)
            this = this@TestCase(name);
        end
    end
    
    methods
        function setUp(this)
            this.simulation = edu.stanford.covert.cell.sim.SimulationFixture.load();
        end
    end
    
    %tests
    methods
        function testSimulationStateSideEffects(this)
            import edu.stanford.covert.cell.sim.SimulationStateSideEffect;
            import edu.stanford.covert.cell.sim.SimulationStateSideEffectItem;

            sim = this.simulation;
 
            sim.state('Metabolite').molecularWeights = [1;2;3];
            sim.state('Metabolite').counts = repmat(10, 3, 2);
            sim.state('ProteinMonomer').counts = repmat(10, 3, 2);
            sim.state('ProteinComplex').counts = repmat(10, 3, 2);
            sim.state('ProteinMonomer').molecularWeights = [1;2;3];
            sim.state('ProteinComplex').molecularWeights = [1;2;3];            
            
            sideEffects = [
                SimulationStateSideEffect([
                    SimulationStateSideEffectItem('Metabolite', 'counts', [], 2, 2, -1);
                    SimulationStateSideEffectItem('Metabolite', 'counts', [], 1, 1, 2)]);
                SimulationStateSideEffect([
                    SimulationStateSideEffectItem('ProteinMonomer', 'counts', [], 1, 2, 3);
                    SimulationStateSideEffectItem('ProteinComplex', 'counts', [], 3, 1, -1);]);
                SimulationStateSideEffect([
                    SimulationStateSideEffectItem('ProteinMonomer', 'counts', [], 1, 2, 1);
                    SimulationStateSideEffectItem('ProteinMonomer', 'counts', [], 2, 1, 1);
                    SimulationStateSideEffectItem('ProteinComplex', 'counts', [], 3, 2, -1)])];
                
            sideEffects.updateSimulationState(sim);
            
            assertEqual(10 + [2 0; 0 -1; 0 0], sim.state('Metabolite').counts);
            assertEqual(10 + [0 4; 1 0; 0 0], sim.state('ProteinMonomer').counts);
            assertEqual(10 + [0 0; 0 0; -1 -1], sim.state('ProteinComplex').counts);                   
        end
        
        function testDisp(~)
            import edu.stanford.covert.cell.sim.SimulationStateSideEffect;
            import edu.stanford.covert.cell.sim.SimulationStateSideEffectItem;            
            
            sideEffects = [
                SimulationStateSideEffect([
                    SimulationStateSideEffectItem('Metabolite', 'counts', [], 2, 2, -1);
                    SimulationStateSideEffectItem('Metabolite', 'counts', [], 1, 1, 2)]);
                SimulationStateSideEffect([
                    SimulationStateSideEffectItem('ProteinMonomer', 'counts', [], 1, 2, 3);
                    SimulationStateSideEffectItem('ProteinComplex', 'counts', [], 3, 1, -1);]);
                SimulationStateSideEffect([
                    SimulationStateSideEffectItem('ProteinMonomer', 'counts', [], 1, 2, 1);
                    SimulationStateSideEffectItem('ProteinMonomer', 'counts', [], 2, 1, 1);
                    SimulationStateSideEffectItem('ProteinComplex', 'counts', [], 3, 2, -1)])];
                
            disp(sideEffects);
        end
    end
end