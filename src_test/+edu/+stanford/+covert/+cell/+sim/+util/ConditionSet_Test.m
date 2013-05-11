%ConditionSet test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updaetd 2/3/2011
classdef ConditionSet_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = ConditionSet_Test(name)
            this = this@TestCase(name);
        end
        
        function setUp(this)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
            sim.applyOptions('lengthSec', 10);
            this.simulation = sim;
        end
    end
    
    methods
        function testGenerateConditionSet(this)
            import edu.stanford.covert.cell.sim.util.ConditionSet;
            
            sim = this.simulation;
            fileName = 'output/runSmallTests/conditions.xml';
            
            expectedConditionSet = this.getTestConditionSet();
            expectedCondition = rmfield(expectedConditionSet, 'metadata');
            expectedCondition.shortDescription = expectedConditionSet.metadata.shortDescription;
            expectedCondition.longDescription = expectedConditionSet.metadata.longDescription;
            expectedCondition.replicates = 1;
            ConditionSet.generateConditionSet(sim, ...
                expectedConditionSet.metadata, ...
                expectedCondition, fileName);
            
            assertEqual(expectedConditionSet, ConditionSet.parseConditionSet(sim, fileName));
        end
        
        function testGenerateSingleGeneDeletionConditionSet(this)
            import edu.stanford.covert.cell.sim.util.ConditionSet;
            
            sim = this.simulation;
            replicates = 10;
            fileName = 'output/runSmallTests/conditions.xml';
            tmp = this.getTestConditionSet();
            metadata = tmp.metadata;
            geneWholeCellModelID = 'MG_001';
            geneName = sim.gene.names{strcmp(sim.gene.wholeCellModelIDs, geneWholeCellModelID)};
            metadata.shortDescription = sprintf('Single-gene (%s; %s) deletion simulation set', geneWholeCellModelID, geneName);
            metadata.longDescription = sprintf('Single-gene (%s; %s) deletion simulation set with %d replicates', geneWholeCellModelID, geneName, replicates);
            
            %single gene
            ConditionSet.generateSingleGeneDeletionConditionSet(sim, ...
                metadata, geneWholeCellModelID, replicates, fileName);
            expectedConditionSet = struct;
            expectedConditionSet.metadata = metadata;
            expectedConditionSet.options = struct('states', struct(), 'processes', struct());
            expectedConditionSet.parameters = struct('states', struct(), 'processes', struct());
            expectedConditionSet.perturbations = struct('geneticKnockouts', [], 'stimulus', [], 'media', []);
            expectedConditionSet.perturbations.geneticKnockouts = {geneWholeCellModelID};
            expectedConditionSet.perturbations.stimulus = zeros(0, 6);
            expectedConditionSet.perturbations.media = zeros(0, 6);
            assertEqual(expectedConditionSet, ConditionSet.parseConditionSet(sim, fileName));
            
            %all genes
            ConditionSet.generateSingleGeneDeletionConditionSet(sim, ...
                metadata, '-all', replicates, fileName);
        end
        
        function testParseConditionSet(this)
            expectedConditionSet = this.getTestConditionSet();
            fileName = 'src_test/+edu/+stanford/+covert/+cell/+sim/+util/fixtures/runSimulations.xml';
            conditionSet = edu.stanford.covert.cell.sim.util.ConditionSet.parseConditionSet(this.simulation, fileName);
            
            assertEqual(expectedConditionSet.metadata, conditionSet.metadata);
            
            assertEqual(expectedConditionSet.options.states, conditionSet.options.states);
            assertEqual(expectedConditionSet.options.processes, conditionSet.options.processes);
            assertEqual(expectedConditionSet.options, conditionSet.options);
            
            assertEqual(expectedConditionSet.parameters.states, conditionSet.parameters.states);
            assertEqual(expectedConditionSet.parameters.processes, conditionSet.parameters.processes);
            assertEqual(expectedConditionSet.parameters, conditionSet.parameters);
            
            assertEqual(expectedConditionSet.perturbations.geneticKnockouts, conditionSet.perturbations.geneticKnockouts);
            assertEqual(expectedConditionSet.perturbations.stimulus, conditionSet.perturbations.stimulus);
            assertEqual(expectedConditionSet.perturbations.media, conditionSet.perturbations.media);
            assertEqual(expectedConditionSet.perturbations, conditionSet.perturbations);
            
            assertEqual(expectedConditionSet, conditionSet);
        end
        
        function value = getTestConditionSet(this)
            sim = this.simulation;
            
            value = struct(...
                'metadata', struct('firstName','Jonathan', ...
                'lastName', 'Karr', 'email', 'jkarr@stanford.edu', ...
                'affiliation','Covert Lab, Department of Bioengineering, Stanford University', ...
                'userName','jkarr','hostName','covertlab-server1.stanford.edu','ipAddress','171.65.102.24','revision',1151,'differencesFromRevision','',...
                'shortDescription', 'Test condition 1', 'longDescription', 'Test condition 1'), ...
                'options', struct('states',struct('Metabolite',struct('option',0)),'processes',struct), ...
                'parameters', struct('states',struct,'processes',struct('Metabolism',struct('parameter',0))), ...
                'perturbations', struct(...
                    'geneticKnockouts', [], ...
                    'stimulus', [], ...
                    'media', []));
            value.perturbations.geneticKnockouts = {'MG_001'; 'MG_003'};
            value.perturbations.stimulus = [
                sim.state('Stimulus').getIndexs('alpha_radiation') sim.compartment.extracellularIndexs 0 0 Inf 0
                sim.state('Stimulus').getIndexs('gamma_radiation') sim.compartment.extracellularIndexs 0 0 Inf 0
                ];
            value.perturbations.media = [
                sim.state('Metabolite').getIndexs('G6P') sim.compartment.extracellularIndexs 0 0 Inf 0
                sim.state('Metabolite').getIndexs('G1P') sim.compartment.extracellularIndexs 0 0 Inf 0
                ];
            value.perturbations.stimulus(:, end) = sub2ind(...
                [numel(sim.state('Stimulus').wholeCellModelIDs) sim.compartment.count], ...
                value.perturbations.stimulus(:, 1), ...
                value.perturbations.stimulus(:, 2));
            value.perturbations.media(:, end) = sub2ind(...
                [numel(sim.state('Metabolite').wholeCellModelIDs) sim.compartment.count], ...
                value.perturbations.media(:, 1), ...
                value.perturbations.media(:, 2));
        end
    end
end