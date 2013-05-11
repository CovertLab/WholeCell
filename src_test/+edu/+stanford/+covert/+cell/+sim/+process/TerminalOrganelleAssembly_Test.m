%Terminal organelle assembly test cases
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef TerminalOrganelleAssembly_Test < edu.stanford.covert.cell.sim.ReactionProcessTestCase
    methods
        function this = TerminalOrganelleAssembly_Test(name)
            this = this@edu.stanford.covert.cell.sim.ReactionProcessTestCase(name);
        end

        function testSimpleFixture(this)
            this.loadSimpleTestFixture();
            m = this.process;

            m.evolveState();
            assertAllEqual(0, m.substrates(:, m.compartmentIndexs_unincorporated));
            assertAllEqual(20, m.substrates(:, m.compartmentIndexs_incorporated));
        end

        function testNoIncorporatedProteins(this)
            m = this.process;
            m.substrates(:, m.compartmentIndexs_incorporated) = 0;
            m.substrates(:, m.compartmentIndexs_unincorporated) = 10;

            m.evolveState();
            assertAllEqual(0, m.substrates(:, m.compartmentIndexs_unincorporated));
            assertAllEqual(10, m.substrates(:, m.compartmentIndexs_incorporated));
        end

        function testNoUnincorporatedProteins(this)
            m = this.process;
            m.substrates(:, m.compartmentIndexs_incorporated) = 10;
            m.substrates(:, m.compartmentIndexs_unincorporated) = 0;

            m.evolveState();
            assertAllEqual(0, m.substrates(:, m.compartmentIndexs_unincorporated));
            assertAllEqual(10, m.substrates(:, m.compartmentIndexs_incorporated));
        end

        function testNoHMW12(this)
            m = this.process;

            m.substrates(:,m.compartmentIndexs_unincorporated) = 1e6;
            m.substrates(:,m.compartmentIndexs_incorporated) = 0;
            m.substrates(m.substrateIndexs_HMW1,:) = 0;
            m.substrates(m.substrateIndexs_HMW2,:) = 0;
            initialState = struct('substrates', m.substrates);
            m.evolveState();
            assertFalse(someUseOfAllSubstrates(m, initialState));

            m.substrates(:) = 100;
            m.substrates(m.substrateIndexs_HMW1, m.compartmentIndexs_incorporated) = 0;
            m.substrates(m.substrateIndexs_HMW2, m.compartmentIndexs_incorporated) = 0;
            initialState = struct('substrates', m.substrates);
            m.evolveState();
            assertTrue(someUseOfAllSubstrates(m, initialState));

            m.substrates(:) = 100;
            m.substrates(m.substrateIndexs_HMW1, m.compartmentIndexs_unincorporated) = 0;
            m.substrates(m.substrateIndexs_HMW2, m.compartmentIndexs_unincorporated) = 0;
            initialState = struct('substrates', m.substrates);
            m.evolveState();
            assertTrue(someUseOfAllSubstrates(m, initialState));
        end

        function testNoHMW1(this)
            m = this.process;

            m.substrates(:) = 100;
            m.substrates(m.substrateIndexs_HMW1,:)=0;
            initialState = struct('substrates', m.substrates);
            m.evolveState();
            assertFalse(someUseOfAllSubstrates(m, initialState));

            m.substrates(:) = 100;
            m.substrates(m.substrateIndexs_HMW1, m.compartmentIndexs_incorporated) = 0;
            initialState = struct('substrates', m.substrates);
            m.evolveState();
            assertTrue(someUseOfAllSubstrates(m, initialState));

            m.substrates(:) = 100;
            m.substrates(m.substrateIndexs_HMW1, m.compartmentIndexs_unincorporated) = 0;
            initialState = struct('substrates', m.substrates);
            m.evolveState();
            assertTrue(someUseOfAllSubstrates(m, initialState));
        end

        function testNoHMW2(this)
            m = this.process;

            m.substrates(:) = 100;
            m.substrates(m.substrateIndexs_HMW2,:) = 0;
            initialState = struct('substrates', m.substrates);
            m.evolveState();
            assertFalse(someUseOfAllSubstrates(m, initialState));

            m.substrates(:) = 100;
            m.substrates(m.substrateIndexs_HMW2, m.compartmentIndexs_incorporated) = 0;
            initialState = struct('substrates', m.substrates);
            m.evolveState();
            assertTrue(someUseOfAllSubstrates(m, initialState));

            m.substrates(:) = 100;
            m.substrates(m.substrateIndexs_HMW2, m.compartmentIndexs_unincorporated) = 0;
            initialState = struct('substrates', m.substrates);
            m.evolveState();
            assertTrue(someUseOfAllSubstrates(m, initialState));
        end

        function testGeneEssentiality(this)
            m = this.process;
            m.substrates(:) = 100;
            this.helpTestGeneEssentiality({
                'MG_218';     %hmw2, high molecular weight cytadherence accessory protein 2
                'MG_312';     %hmw1, high molecular weight cytadherence accessory protein 1
                'MG_317';     %hmw3, high molecular weight cytadherence accessory protein 3
                'MG_318'},... %p32, P32 adhesin
                @someUseOfAllSubstrates);
        end
    end
    
    methods (Access = private)
        function loadSimpleTestFixture(this)
            m = this.process;

            enzymeSubstrateIndexs = [
                m.enzymeIndexs_HMW1;
                m.enzymeIndexs_HMW2;
                m.enzymeIndexs_HMW3;
                m.enzymeIndexs_P32];

            m.substrateWholeCellModelIDs = sort(m.substrateWholeCellModelIDs__);
            m.enzymeWholeCellModelIDs = m.substrateWholeCellModelIDs(enzymeSubstrateIndexs);

            m.reactionWholeCellModelIDs = {
                'MG_191_TerminalOrganelleAssembly_01';
                'MG_192_TerminalOrganelleAssembly_01';
                'MG_217_TerminalOrganelleAssembly_01';
                'MG_218_TerminalOrganelleAssembly_01';
                'MG_218_TerminalOrganelleAssembly_02';
                'MG_312_TerminalOrganelleAssembly_01';
                'MG_312_TerminalOrganelleAssembly_02';
                'MG_317_TerminalOrganelleAssembly_01';
                'MG_318_TerminalOrganelleAssembly_01';
                'MG_386_TerminalOrganelleAssembly_01'};
            m.reactionNames = m.reactionWholeCellModelIDs;

            nSubstrates = length(m.substrateWholeCellModelIDs);
            nReactions  = length(m.reactionWholeCellModelIDs);

            m.localizationSubstrates = zeros(nReactions, nSubstrates);
            m.localizationSubstrates(1, substrateIdx('MG_191_MONOMER')) = 1;
            m.localizationSubstrates(2, substrateIdx('MG_192_MONOMER')) = 1;
            m.localizationSubstrates(3, substrateIdx('MG_217_MONOMER')) = 1;
            m.localizationSubstrates(4, substrateIdx('MG_218_MONOMER')) = 1;
            m.localizationSubstrates(5, substrateIdx('MG_218_MONOMER')) = 1;
            m.localizationSubstrates(6, substrateIdx('MG_312_MONOMER')) = 1;
            m.localizationSubstrates(7, substrateIdx('MG_312_MONOMER')) = 1;
            m.localizationSubstrates(8, substrateIdx('MG_317_MONOMER')) = 1;
            m.localizationSubstrates(9, substrateIdx('MG_318_MONOMER')) = 1;
            m.localizationSubstrates(10, substrateIdx('MG_386_MONOMER')) = 1;

            i = m.compartmentIndexs_incorporated;
            u = m.compartmentIndexs_unincorporated;
            
            m.localizationReactions = zeros(nReactions, nSubstrates, 2);
            m.localizationReactions(1, substrateIdx('MG_312_MONOMER'), i) = 1;
            m.localizationReactions(3, substrateIdx('MG_318_MONOMER'), i) = 1;
            m.localizationReactions(4, substrateIdx('MG_312_MONOMER'), u) = 1;
            m.localizationReactions(5, substrateIdx('MG_312_MONOMER'), i) = 1;
            m.localizationReactions(6, substrateIdx('MG_218_MONOMER'), u) = 1;
            m.localizationReactions(7, substrateIdx('MG_218_MONOMER'), i) = 1;
            m.localizationReactions(8, substrateIdx('MG_312_MONOMER'), i) = 1;
            m.localizationReactions(9, substrateIdx('MG_317_MONOMER'), i) = 1;
            m.localizationReactions(10, substrateIdx('MG_312_MONOMER'), i) = 1;

            m.localizationThreshold = sum(sum(m.localizationReactions,3),2);

            m.substrates = repmat(10, length(m.substrateWholeCellModelIDs), 2);
            m.enzymes    = m.substrates(enzymeSubstrateIndexs, :);
            
            function i = substrateIdx(id)
                [~,i] = ismember(id, m.substrateWholeCellModelIDs);
            end
        end
    end
end

function result = someUseOfAllSubstrates(m, i)
    result = all(...
        i.substrates(:,m.compartmentIndexs_unincorporated)==0 | ...
        i.substrates(:,m.compartmentIndexs_unincorporated) > ...
        m.substrates(:,m.compartmentIndexs_unincorporated));
end
