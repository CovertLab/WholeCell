%Rna_Test
% RNA test class.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 1/4/2011
classdef Rna_Test < edu.stanford.covert.cell.sim.CellStateTestCase
    methods
        function this = Rna_Test(name)
            this = this@edu.stanford.covert.cell.sim.CellStateTestCase(name);
        end
    end
    
    methods
        function testMolecularWeights(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            s = this.state;
            m = edu.stanford.covert.cell.sim.CellStateFixture.load(edu.stanford.covert.cell.sim.state.Metabolite('', ''));
            molecularWeights = s.baseCounts * m.molecularWeights - ...
                (ConstantUtil.elements.O + ConstantUtil.elements.H) * max(0, s.lengths-1);
            molecularWeights(s.aminoacylatedIndexs) = ...
                molecularWeights(s.aminoacylatedIndexs) - ...
                (ConstantUtil.elements.O + ConstantUtil.elements.H) * ...
                (molecularWeights(s.aminoacylatedIndexs) > molecularWeights(s.matureIndexs));
            assertElementsAlmostEqual(s.molecularWeights, ...
                molecularWeights, ...
                'relative', 1e-8, 0);
        end
    end
end