%Protein folding process test case
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef ProteinFolding_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = ProteinFolding_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end

    %fixtures
    methods
        function loadSimpleTestFixture(this)
            %process
            m=this.process;

            %whole cell model IDs
            m.substrateWholeCellModelIDs       = {'ADP'; 'ATP'; 'FE2'; 'H'; 'H2O'; 'K'; 'MG'; 'MN'; 'NA'; 'PI'; 'ZN'};
            m.enzymeWholeCellModelIDs          = {'MG_019_DIMER'; 'MG_201_DIMER'; 'MG_238_MONOMER'; 'MG_305_MONOMER'; 'MG_392_393_21MER'};
            m.unfoldedMonomerWholeCellModelIDs = {'MG_001_MONOMER'; 'MG_003_MONOMER'; 'MG_004_MONOMER'; 'MG_006_MONOMER'; 'MG_394_MONOMER'};
            m.unfoldedComplexWholeCellModelIDs = {'DNA_GYRASE'; 'MG_394_TETRAMER';};
            m.foldedMonomerWholeCellModelIDs   = m.unfoldedMonomerWholeCellModelIDs;
            m.foldedComplexWholeCellModelIDs   = m.unfoldedComplexWholeCellModelIDs;

            %names
            m.substrateNames = m.substrateWholeCellModelIDs;
            m.enzymeNames    = m.enzymeWholeCellModelIDs;

            %indices
            m.substrateIndexs_adp        = m.substrateIndexs({'ADP'});
            m.substrateIndexs_atp        = m.substrateIndexs({'ATP'});
            m.substrateIndexs_fe2        = m.substrateIndexs({'FE2'});
            m.substrateIndexs_hydrogen   = m.substrateIndexs({'H'});
            m.substrateIndexs_water      = m.substrateIndexs({'H2O'});
            m.substrateIndexs_k          = m.substrateIndexs({'K'});
            m.substrateIndexs_mg         = m.substrateIndexs({'MG'});
            m.substrateIndexs_mn         = m.substrateIndexs({'MN'});
            m.substrateIndexs_sodium     = m.substrateIndexs({'NA'});
            m.substrateIndexs_phosphate  = m.substrateIndexs({'PI'});
            m.substrateIndexs_zinc       = m.substrateIndexs({'ZN'});

            m.enzymeIndexs_triggerFactor = m.enzymeIndexs({'MG_238_MONOMER'});
            m.enzymeIndexs_dnaK          = m.enzymeIndexs({'MG_305_MONOMER'});
            m.enzymeIndexs_dnaJ          = m.enzymeIndexs({'MG_019_DIMER'});
            m.enzymeIndexs_grpE          = m.enzymeIndexs({'MG_201_DIMER'});
            m.enzymeIndexs_groELES       = m.enzymeIndexs({'MG_392_393_21MER'});

            %compartments
            m.monomerCompartments = ones(size(m.unfoldedMonomerWholeCellModelIDs));
            m.complexCompartments = ones(size(m.unfoldedComplexWholeCellModelIDs));

            %prosthetic groups required for folding
            monomerProstheticGroupMatrix = zeros(length(m.unfoldedMonomerWholeCellModelIDs), length(m.substrateWholeCellModelIDs));
            complexProstheticGroupMatrix = zeros(length(m.unfoldedComplexWholeCellModelIDs), length(m.substrateWholeCellModelIDs));

            monomerProstheticGroupMatrix(strcmp('MG_006_MONOMER',  m.unfoldedMonomerWholeCellModelIDs), m.substrateIndexs_mg)=1;
            monomerProstheticGroupMatrix(strcmp('MG_394_MONOMER',  m.unfoldedMonomerWholeCellModelIDs), m.substrateIndexs_zinc)=1;

            m.proteinProstheticGroupMatrix = [monomerProstheticGroupMatrix; complexProstheticGroupMatrix];

            %chaperones required for folding
            monomerChaperoneMatrix = zeros(length(m.unfoldedMonomerWholeCellModelIDs), length(m.enzymeWholeCellModelIDs));
            complexChaperoneMatrix = zeros(length(m.unfoldedComplexWholeCellModelIDs), length(m.enzymeWholeCellModelIDs));

            monomerChaperoneMatrix(strcmp('MG_001_MONOMER',  m.unfoldedMonomerWholeCellModelIDs), [])=1;
            monomerChaperoneMatrix(strcmp('MG_003_MONOMER',  m.unfoldedMonomerWholeCellModelIDs), [])=1;
            monomerChaperoneMatrix(strcmp('MG_004_MONOMER',  m.unfoldedMonomerWholeCellModelIDs), [m.enzymeIndexs_dnaK; m.enzymeIndexs_groELES])=1;
            monomerChaperoneMatrix(strcmp('MG_006_MONOMER',  m.unfoldedMonomerWholeCellModelIDs), [])=1;
            monomerChaperoneMatrix(strcmp('MG_394_MONOMER',  m.unfoldedMonomerWholeCellModelIDs), m.enzymeIndexs_groELES)=1;
            complexChaperoneMatrix(strcmp('MG_394_TETRAMER', m.unfoldedComplexWholeCellModelIDs), [])=1;
            complexChaperoneMatrix(strcmp('DNA_GYRASE',      m.unfoldedComplexWholeCellModelIDs), [])=1;

            m.proteinChaperoneMatrix = [monomerChaperoneMatrix; complexChaperoneMatrix];
            m.proteinChaperoneMatrix(1:size(monomerChaperoneMatrix,1), m.enzymeIndexs_triggerFactor)=1;
            m.proteinChaperoneMatrix(:,m.enzymeIndexs_dnaJ) = m.proteinChaperoneMatrix(:,m.enzymeIndexs_dnaK);
            m.proteinChaperoneMatrix(:,m.enzymeIndexs_grpE) = m.proteinChaperoneMatrix(:,m.enzymeIndexs_dnaK);
            
            m.initializeSpeciesNetwork();

            %molecular weights
            m.substrateMolecularWeights       = [424.1769  503.1489         0    1.0079   18.0152   39.0983         0   14.0067         0   95.9793   14.0067]';
            m.enzymeMolecularWeights          = 1e5 * [0.8756    0.5005    0.5094    0.6507    9.0263]';
            m.unfoldedMonomerMolecularWeights = 1e5 * (1:length(m.unfoldedMonomerWholeCellModelIDs))';
            m.unfoldedComplexMolecularWeights = 1e6 * (1:length(m.unfoldedComplexWholeCellModelIDs))';
            m.foldedMonomerMolecularWeights   = m.unfoldedMonomerMolecularWeights + monomerProstheticGroupMatrix*m.substrateMolecularWeights;
            m.foldedComplexMolecularWeights   = m.unfoldedComplexMolecularWeights + complexProstheticGroupMatrix*m.substrateMolecularWeights;

            %initial state
            m.substrates       = repmat(1e6,length(m.substrateWholeCellModelIDs),  1);
            m.enzymes          = ones(length(m.enzymeWholeCellModelIDs),           1);
            m.boundEnzymes     = zeros(length(m.enzymeWholeCellModelIDs),          1);
            m.unfoldedMonomers = ones(length(m.unfoldedMonomerWholeCellModelIDs),  1);
            m.unfoldedComplexs = ones(length(m.unfoldedComplexWholeCellModelIDs),  1);
            m.foldedMonomers   = zeros(length(m.foldedMonomerWholeCellModelIDs),   1);
            m.foldedComplexs   = zeros(length(m.foldedComplexWholeCellModelIDs),   1);
        end
    end

    %tests
    methods
        function testConstants(this)
             m = this.process;
             
             %all monomers fold
             assertTrue(isempty(m.monomerIndexs_notFolding));
        end
        
        function testOneMonomerWithoutProstheticGroup(this)
            m = this.process;
            [~,i] = ismember('MG_001_MONOMER', m.unfoldedMonomerWholeCellModelIDs);
            assertFalse(isempty(i));
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_triggerFactor) = 1;
            m.unfoldedMonomers(:) = 0;
            m.unfoldedMonomers(i) = 1;
            m.unfoldedComplexs(:) = 0;
            m.foldedMonomers(:) = 0;
            m.foldedComplexs(:) = 0;

            m.evolveState();
            assertEqual(zeros(size(m.unfoldedMonomers)), m.unfoldedMonomers);
            assertEqual(1, m.foldedMonomers(i));
            assertEqual(1, sum(m.foldedMonomers));
            assertEqual(zeros(size(m.unfoldedComplexs)), m.unfoldedComplexs);
            assertEqual(zeros(size(m.foldedComplexs)), m.foldedComplexs);
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(1, m.enzymes(m.enzymeIndexs_triggerFactor));
        end

        function testNoMonomersFoldWithoutTriggerFactor(this)
            m = this.process;
            m.substrates(:) = 1e6;
            m.enzymes(:) = 0;
            m.unfoldedMonomers(:) = 1;
            m.unfoldedComplexs(:) = 0;
            m.foldedMonomers(:) = 0;
            m.foldedComplexs(:) = 0;

            m.evolveState();
            assertEqual(ones(size(m.unfoldedMonomers)), m.unfoldedMonomers);
            assertEqual(zeros(size(m.foldedMonomers)), m.foldedMonomers);
            assertEqual(zeros(size(m.unfoldedComplexs)), m.unfoldedComplexs);
            assertEqual(zeros(size(m.foldedComplexs)), m.foldedComplexs);
            assertEqual(1e6 * ones(size(m.substrates)), m.substrates);
        end

        function testOneMonomerWithProstheticGroup(this)
            m = this.process;
            [~,i] = ismember('MG_006_MONOMER', m.unfoldedMonomerWholeCellModelIDs);
            assertFalse(isempty(i));
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_mg) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_triggerFactor) = 1;
            m.unfoldedMonomers(:) = 0;
            m.unfoldedMonomers(i) = 1;
            m.unfoldedComplexs(:) = 0;
            m.foldedMonomers(:) = 0;
            m.foldedComplexs(:) = 0;

            m.evolveState();
            assertEqual(zeros(size(m.unfoldedMonomers)), m.unfoldedMonomers);
            assertEqual(1, m.foldedMonomers(i));
            assertEqual(1, sum(m.foldedMonomers));
            assertEqual(zeros(size(m.unfoldedComplexs)), m.unfoldedComplexs);
            assertEqual(zeros(size(m.foldedComplexs)), m.foldedComplexs);
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(1, m.enzymes(m.enzymeIndexs_triggerFactor));
        end

        function testOneMonomerWithProstheticGroup_noIon(this)
            m = this.process;
            [~,i] = ismember('MG_006_MONOMER', m.unfoldedMonomerWholeCellModelIDs);
            assertFalse(isempty(i));
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_triggerFactor) = 1;
            m.unfoldedMonomers(:) = 0;
            m.unfoldedMonomers(i) = 1;
            m.unfoldedComplexs(:) = 0;
            m.foldedMonomers(:) = 0;
            m.foldedComplexs(:) = 0;

            m.evolveState();
            assertEqual(1, m.unfoldedMonomers(i));
            assertEqual(zeros(size(m.foldedMonomers)), m.foldedMonomers);
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(1, m.enzymes(m.enzymeIndexs_triggerFactor));
        end

        function testSelfFoldingComplex(this)
            m = this.process;
            [~,i] = ismember('DNA_GYRASE', m.unfoldedComplexWholeCellModelIDs);
            assertFalse(isempty(i));
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.unfoldedMonomers(:) = 0;
            m.unfoldedComplexs(:) = 0;
            m.unfoldedComplexs(i) = 1;
            m.foldedMonomers(:) = 0;
            m.foldedComplexs(:) = 0;

            m.evolveState();
            assertEqual(zeros(size(m.unfoldedMonomers)), m.unfoldedMonomers);
            assertEqual(zeros(size(m.foldedMonomers)), m.foldedMonomers);
            assertEqual(zeros(size(m.unfoldedComplexs)), m.unfoldedComplexs);
            assertEqual(1, m.foldedComplexs(i));
            assertEqual(1, sum(m.foldedComplexs));
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(zeros(size(m.enzymes)), m.enzymes);
        end

        function testOneComplexRequiringMultipleChaperones(this)
            m = this.process;
            [~,i] = ismember('MG_394_TETRAMER', m.unfoldedComplexWholeCellModelIDs);
            assertFalse(isempty(i));
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_triggerFactor) = 1;
            m.enzymes(m.enzymeIndexs_dnaK) = 1;
            m.enzymes(m.enzymeIndexs_dnaJ) = 1;
            m.enzymes(m.enzymeIndexs_grpE) = 1;
            m.enzymes(m.enzymeIndexs_groELES) = 1;
            m.unfoldedMonomers(:) = 0;
            m.unfoldedComplexs(:) = 0;
            m.unfoldedComplexs(i) = 1;
            m.foldedMonomers(:) = 0;
            m.foldedComplexs(:) = 0;

            m.evolveState();
            assertEqual(zeros(size(m.unfoldedMonomers)), m.unfoldedMonomers);
            assertEqual(zeros(size(m.foldedMonomers)), m.foldedMonomers);
            assertEqual(zeros(size(m.unfoldedComplexs)), m.unfoldedComplexs);
            assertEqual(1, m.foldedComplexs(i));
            assertEqual(1, sum(m.foldedComplexs));
            assertEqual(zeros(size(m.substrates)), m.substrates);
        end

        function testComplexFoldingWithMissingChaperone(this)
            m = this.process;
            [~,i] = ismember('MG_394_TETRAMER', m.unfoldedComplexWholeCellModelIDs);
            assertFalse(isempty(i));
            chaperones = [
                m.enzymeIndexs_triggerFactor
                m.enzymeIndexs_dnaK
                m.enzymeIndexs_dnaJ
                m.enzymeIndexs_grpE
                m.enzymeIndexs_groELES];
            for j = 1:length(chaperones)
                m.substrates(:) = 0;
                m.enzymes(:) = 0;
                m.enzymes(chaperones) = 1;
                m.enzymes(chaperones(j)) = 0;
                m.unfoldedMonomers(:) = 0;
                m.unfoldedComplexs(:) = 0;
                m.unfoldedComplexs(i) = 1;
                m.foldedMonomers(:) = 0;
                m.foldedComplexs(:) = 0;

                m.evolveState();
                assertEqual(1, m.foldedComplexs(i),...
                    ['expected folding, j=' num2str(j)]);
            end
        end

        function testFoldEveryKindOfProtein(this)
            n = 1;  %how many of each protein
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(any(m.proteinProstheticGroupMatrix)) = 1e6;
            m.enzymes(:)          = 1;
            m.unfoldedMonomers(:) = n;
            m.unfoldedComplexs(:) = n;
            m.foldedMonomers(:)   = 0;
            m.foldedComplexs(:)   = 0;

            m.evolveState();
            assertEqual(zeros(size(m.unfoldedMonomers)), m.unfoldedMonomers);
            assertEqual(n * ones(size(m.foldedMonomers)), m.foldedMonomers);
            assertEqual(zeros(size(m.unfoldedComplexs)), m.unfoldedComplexs);
            assertEqual(n * ones(size(m.foldedComplexs)), m.foldedComplexs);
        end

        function testGeneEssentiality(this)
            % Provide everything necessary to fold one of each kind of protein.
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(any(m.proteinProstheticGroupMatrix)) = 1e6;
            m.enzymes(:)          = 1;
            m.unfoldedMonomers(:) = 1;
            m.unfoldedComplexs(:) = 1;
            m.foldedMonomers(:)   = 0;
            m.foldedComplexs(:)   = 0;

            this.helpTestGeneEssentiality({
                'MG_019';     %chaperone protein DnaJ
                'MG_201';     %co-chaperone GrpE
                'MG_238';     %trigger factor
                'MG_305';     %chaperone
                'MG_392';     %chaperonin
                'MG_393'},... %chaperonin, 10 kDa
                @(m,i) all(i.foldedMonomers < m.foldedMonomers) && ...
                       all(i.foldedComplexs < m.foldedComplexs));
        end
    end
end
