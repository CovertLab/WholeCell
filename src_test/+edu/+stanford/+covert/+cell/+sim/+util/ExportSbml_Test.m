classdef ExportSbml_Test < TestCase
    methods
        function this = ExportSbml_Test(name)
            this = this@TestCase(name);
        end
    end
    
    methods
        function testEncodeRNADecay(~)
            import edu.stanford.covert.cell.sim.util.ExportSbml;
            
            %load model
            load('data\Simulation_fitted.mat');
            sim = simulation;
            
            %export model to SBML
            [mdl, species, enzymes, reactions] = ExportSbml.constructRNADecaySimBiologyModel(sim);
            
            outDir = 'output/sbml';
            ExportSbml.exportModel(mdl, species, enzymes, reactions, fullfile(outDir, 'RNADecay.xml'), fullfile(outDir, 'RNADecay.xls'));            
        end
        
        function testExport(~)
            import edu.stanford.covert.cell.sim.util.ExportSbml;
            
            %load model
            load('data\Simulation_fitted.mat');
            sim = simulation;
            
            %export model to SBML
            ExportSbml.run(sim, 'output/sbml');
        end
    end
end