%Host interaction process test case
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/27/2011
classdef HostInteraction_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = HostInteraction_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end
    
    methods
        function testProteinComparatments(this)
            m = this.process;
            c = m.compartment;
            
            assertTrue(all(ismember(m.enzymeMonomerCompartmentIndexs(ismember(m.enzymeMonomerLocalIndexs, m.enzymeIndexs_terminalOrganelle)), ...
                [c.terminalOrganelleCytosolIndexs; c.terminalOrganelleMembraneIndexs])));
            assertTrue(all(ismember(m.enzymeComplexCompartmentIndexs(ismember(m.enzymeComplexLocalIndexs, m.enzymeIndexs_terminalOrganelle)), ...
                [c.terminalOrganelleCytosolIndexs; c.terminalOrganelleMembraneIndexs])));
            
            assertFalse(any(ismember(m.enzymeMonomerCompartmentIndexs(~ismember(m.enzymeMonomerLocalIndexs, m.enzymeIndexs_terminalOrganelle)), ...
                [c.terminalOrganelleCytosolIndexs; c.terminalOrganelleMembraneIndexs])));
            assertFalse(any(ismember(m.enzymeComplexCompartmentIndexs(~ismember(m.enzymeComplexLocalIndexs, m.enzymeIndexs_terminalOrganelle)), ...
                [c.terminalOrganelleCytosolIndexs; c.terminalOrganelleMembraneIndexs])));
        end
        
        function testAdherent(this)
            m = this.process();
            h = m.host;
            
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_terminalOrganelle) = 1;
            
            m.evolveState();
            
            assertEqual(true, h.isBacteriumAdherent);
            assertAllEqual(false, h.isTLRActivated);
            assertEqual(false, h.isNFkBActivated);
            assertEqual(false, h.isInflammatoryResponseActivated);
        end
        
        function testNotAdherent(this)
            m = this.process();
            h = m.host;
            
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_terminalOrganelle(1:end-1)) = 1;
            
            m.evolveState();
            
            assertEqual(false, h.isBacteriumAdherent);
            assertAllEqual(false, h.isTLRActivated);
            assertEqual(false, h.isNFkBActivated);
            assertEqual(false, h.isInflammatoryResponseActivated);
        end
        
        function testTLR26Activation(this)
            m = this.process();
            h = m.host;
            
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_terminalOrganelle) = 1;
            m.enzymes(m.enzymeIndexs_tlr26Ligand) = 1;
            
            m.evolveState();
            
            assertEqual(true, h.isBacteriumAdherent);
            assertEqual(false, h.isTLRActivated(h.tlrIndexs_1));
            assertEqual(true, h.isTLRActivated(h.tlrIndexs_2));
            assertEqual(true, h.isTLRActivated(h.tlrIndexs_6));
            assertEqual(true, h.isNFkBActivated);
            assertEqual(true, h.isInflammatoryResponseActivated);
        end
        
        function testNoTLR26Activation(this)
            m = this.process();
            h = m.host;
            
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_terminalOrganelle(1:end-1)) = 1;
            m.enzymes(m.enzymeIndexs_tlr26Ligand) = 1;
            
            m.evolveState();
            
            assertEqual(false, h.isBacteriumAdherent);
            assertEqual(false, h.isTLRActivated(h.tlrIndexs_1));
            assertEqual(false, h.isTLRActivated(h.tlrIndexs_2));
            assertEqual(false, h.isTLRActivated(h.tlrIndexs_6));
            assertEqual(false, h.isNFkBActivated);
            assertEqual(false, h.isInflammatoryResponseActivated);
        end
        
        function testGeneEssentiality(this)
            this.helpTestGeneEssentiality({}, @(~,~) true);
        end
    end
end
