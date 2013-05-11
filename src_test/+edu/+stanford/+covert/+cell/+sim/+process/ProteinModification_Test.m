%Protein modification process test case
%
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/29/2010
classdef ProteinModification_Test < edu.stanford.covert.cell.sim.ReactionProcessTestCase

    %monomer indices
    properties
        MG_001_MONOMER
        MG_070_MONOMER
        MG_166_MONOMER
        MG_217_MONOMER
        MG_272_MONOMER
    end

    methods
        function this = ProteinModification_Test(name)
            this = this@edu.stanford.covert.cell.sim.ReactionProcessTestCase(name);
        end
    end

    methods
        function loadTestFixture(this)
            this.loadTestFixture@edu.stanford.covert.cell.sim.ReactionProcessTestCase();

            m = this.process;
            
            this.MG_001_MONOMER = find(strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_001_MONOMER'),1,'first');
            this.MG_070_MONOMER = find(strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_070_MONOMER'),1,'first');
            this.MG_166_MONOMER = find(strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_166_MONOMER'),1,'first');
            this.MG_217_MONOMER = find(strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_217_MONOMER'),1,'first');
            this.MG_272_MONOMER = find(strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_272_MONOMER'),1,'first');
        end
    end

    methods
        function loadSimpleTestFixture(this)
            %process
            m = this.process;

            %Whole Cell model IDs
            m.substrateWholeCellModelIDs = {'ADP' 'AMP' 'ATP' 'GLU' 'H' 'LIPOYLAMP' 'LIPOYLLYS' 'LYS' 'PI' 'pSER' 'pTHR' 'pTYR' 'SER' 'THR' 'TYR'}';
            m.enzymeWholeCellModelIDs = {'MG_012_MONOMER' 'MG_270_MONOMER' 'MG_109_DIMER'}';
            m.reactionWholeCellModelIDs = {...
                'MG_070_MONOMER_MODIFICATION_001';...
                'MG_070_MONOMER_MODIFICATION_002';...
                'MG_070_MONOMER_MODIFICATION_003';...
                'MG_070_MONOMER_MODIFICATION_004';...
                'MG_070_MONOMER_MODIFICATION_005';...
                'MG_166_MONOMER_MODIFICATION_001';...
                'MG_217_MONOMER_MODIFICATION_001';...
                'MG_272_MONOMER_MODIFICATION_001'};
            m.unmodifiedMonomerWholeCellModelIDs = {'MG_001_MONOMER'; 'MG_070_MONOMER'; 'MG_166_MONOMER'; 'MG_217_MONOMER'; 'MG_272_MONOMER'};
            m.modifiedMonomerWholeCellModelIDs = m.unmodifiedMonomerWholeCellModelIDs;

            %names
            m.substrateNames = m.substrateWholeCellModelIDs;
            m.enzymeNames    = m.enzymeWholeCellModelIDs;
            m.reactionNames  = m.reactionWholeCellModelIDs;

            %types
            m.reactionTypes = mat2cell(repmat('adduction',8,1), repmat(1,8,1), length('adduction'));
            m.reactionTypes{6}='ligation';

            %indices
            m.substrateIndexs_aminoAcids         = m.substrateIndexs({'GLU';'LYS';'SER';'THR';'TYR';});
            m.substrateIndexs_modifiedAminoAcids = m.substrateIndexs({'LIPOYLLYS';'pSER';'pTHR';'pTYR'});
            m.substrateIndexs_atp                = m.substrateIndexs({'ATP'});
            m.substrateIndexs_adp                = m.substrateIndexs({'ADP'});
            m.substrateIndexs_amp                = m.substrateIndexs({'AMP'});
            m.substrateIndexs_hydrogen           = m.substrateIndexs({'H'});
            m.substrateIndexs_phosphate          = m.substrateIndexs({'PI'});
            m.substrateIndexs_glutamate          = m.substrateIndexs({'GLU'});
            m.substrateIndexs_lipoylAmp          = m.substrateIndexs({'LIPOYLAMP'});
            m.substrateIndexs_lipoylLys          = m.substrateIndexs({'LIPOYLLYS'});

            m.enzymeIndexs_serineThreonineKinase = m.enzymeIndexs({'MG_109_DIMER'});
            m.enzymeIndexs_lipoylTransferase     = m.enzymeIndexs({'MG_270_MONOMER'});
            m.enzymeIndexs_glutamateLigase       = m.enzymeIndexs({'MG_012_MONOMER'});

            m.reactionIndexs_adduction           = find(strcmp(m.reactionTypes,'adduction'));
            m.reactionIndexs_ligation            = find(strcmp(m.reactionTypes,'ligation'));

            m.monomerCompartments                = ones(size(m.unmodifiedMonomerWholeCellModelIDs));

            %other physical properties
            m.matureMonomerLengths = ones(size(m.unmodifiedMonomerWholeCellModelIDs));

            %proteins modified by each reaction
            m.reactionModificationMatrix = zeros(length(m.reactionWholeCellModelIDs), length(m.unmodifiedMonomerWholeCellModelIDs));
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG_070_MONOMER_MODIFICATION_001'), strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_070_MONOMER'))=1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG_070_MONOMER_MODIFICATION_002'), strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_070_MONOMER'))=1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG_070_MONOMER_MODIFICATION_003'), strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_070_MONOMER'))=1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG_070_MONOMER_MODIFICATION_004'), strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_070_MONOMER'))=1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG_070_MONOMER_MODIFICATION_005'), strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_070_MONOMER'))=1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG_166_MONOMER_MODIFICATION_001'), strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_166_MONOMER'))=1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG_217_MONOMER_MODIFICATION_001'), strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_217_MONOMER'))=1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG_272_MONOMER_MODIFICATION_001'), strcmp(m.unmodifiedMonomerWholeCellModelIDs,'MG_272_MONOMER'))=1;

            %enzyme which catalyze each reaction
            m.reactionCatalysisMatrix = zeros(length(m.reactionWholeCellModelIDs), length(m.enzymeWholeCellModelIDs));
            m.reactionCatalysisMatrix(:,m.enzymeIndexs_serineThreonineKinase)=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG_166_MONOMER_MODIFICATION_001'),:)=0;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG_272_MONOMER_MODIFICATION_001'),:)=0;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG_166_MONOMER_MODIFICATION_001'),m.enzymeIndexs_glutamateLigase)=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG_272_MONOMER_MODIFICATION_001'),m.enzymeIndexs_lipoylTransferase)=1;

            %reactions
            m.reactionStoichiometryMatrix = zeros(length(m.substrateWholeCellModelIDs), length(m.reactionWholeCellModelIDs));
            m.reactionStoichiometryMatrix([m.substrateIndexs_atp;m.substrateIndexs_adp;m.substrateIndexs_hydrogen], m.reactionCatalysisMatrix(:,m.enzymeIndexs_serineThreonineKinase)>0)=repmat([-1 1 1]', 1, sum(m.reactionCatalysisMatrix(:,m.enzymeIndexs_serineThreonineKinase)));
            m.reactionStoichiometryMatrix([m.substrateIndexs_atp;m.substrateIndexs_glutamate;m.substrateIndexs_adp;m.substrateIndexs_hydrogen;m.substrateIndexs_phosphate], m.reactionCatalysisMatrix(:,m.enzymeIndexs_glutamateLigase)>0)=repmat([-4 -4 4 4 4]', 1, sum(m.reactionCatalysisMatrix(:,m.enzymeIndexs_glutamateLigase)));
            m.reactionStoichiometryMatrix([m.substrateIndexs_lipoylAmp;m.substrateIndexs_amp;m.substrateIndexs_hydrogen], m.reactionCatalysisMatrix(:,m.enzymeIndexs_lipoylTransferase)>0)=repmat([-1 1 2]', 1, sum(m.reactionCatalysisMatrix(:,m.enzymeIndexs_lipoylTransferase)));

            %reaction kinetics
            m.enzymeBounds = repmat([-Inf Inf], length(m.reactionWholeCellModelIDs), 2);
            m.enzymeBounds(m.reactionCatalysisMatrix(:,m.enzymeIndexs_serineThreonineKinase)>0, 2) = 0.0427;
            m.enzymeBounds(m.reactionCatalysisMatrix(:,m.enzymeIndexs_glutamateLigase)>0,       2) = 3.9321e-004;
            m.enzymeBounds(m.reactionCatalysisMatrix(:,m.enzymeIndexs_lipoylTransferase)>0,     2) = 31.6347;
            
            m.initializeSpeciesNetwork();

            %molecular weights
            m.substrateMolecularWeights          = [424.1769 345.2049 503.1489 146.1210  1.0079 534.5226 334.4968 147.1949  95.9793 183.0564 197.0829 259.1522 105.0923 119.1188 181.1881]';
            m.enzymeMolecularWeights             = 1e4 * [3.2768    3.8947    8.9271]';
            m.unmodifiedMonomerMolecularWeights  = 1e5 * (1:size(m.unmodifiedMonomerWholeCellModelIDs,1))';
            m.modifiedMonomerMolecularWeights    = m.unmodifiedMonomerMolecularWeights - ...
                m.reactionModificationMatrix'*m.reactionStoichiometryMatrix'*m.substrateMolecularWeights;

            %initial state
            m.substrates         = zeros(length(m.substrateWholeCellModelIDs),         1);
            m.enzymes            = zeros(length(m.enzymeWholeCellModelIDs),            1);
            m.boundEnzymes       = zeros(length(m.enzymeWholeCellModelIDs),            1);
            m.unmodifiedMonomers = zeros(length(m.unmodifiedMonomerWholeCellModelIDs), 1);
            m.modifiedMonomers   = zeros(length(m.modifiedMonomerWholeCellModelIDs),   1);
        end
    end

    %tests
    methods
        function testConstants(this)
            m = this.process;
            assertFalse(any(ismember(m.substrateWholeCellModelIDs, 'H2O')));
        end
        
        function testMonomerWithNoModifications(this)
            m = this.process;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(this.MG_001_MONOMER) = 1;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(0, m.unmodifiedMonomers(this.MG_001_MONOMER));
            assertEqual(1, m.modifiedMonomers(this.MG_001_MONOMER));
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(zeros(size(m.enzymes)), m.enzymes);
        end

        function testOnePhosphorylation(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_serineThreonineKinase) = 1e6;
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(this.MG_217_MONOMER) = 1;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(0,   m.unmodifiedMonomers(this.MG_217_MONOMER));
            assertEqual(1,   m.modifiedMonomers(this.MG_217_MONOMER));
            assertEqual(0,   m.substrates(m.substrateIndexs_atp));
            assertEqual(1,   m.substrates(m.substrateIndexs_adp));
            assertEqual(1,   m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_serineThreonineKinase));
        end

        function testOnePhosphorylation_noEnyzme(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.enzymes(:) = 0;
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(this.MG_217_MONOMER) = 1;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(1, m.unmodifiedMonomers(this.MG_217_MONOMER));
            assertEqual(0, m.modifiedMonomers(this.MG_217_MONOMER));
            assertEqual(1, m.substrates(m.substrateIndexs_atp));
            assertEqual(0, m.substrates(m.substrateIndexs_adp));
            assertEqual(0, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(0, m.enzymes(m.enzymeIndexs_serineThreonineKinase));
        end

        function testOnePhosphorylation_noATP(this)
            m = this.process;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_serineThreonineKinase) = 1;
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(this.MG_217_MONOMER) = 1;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(1, m.unmodifiedMonomers(this.MG_217_MONOMER));
            assertEqual(0, m.modifiedMonomers(this.MG_217_MONOMER));
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(1, m.enzymes(m.enzymeIndexs_serineThreonineKinase));
        end

        function testMonomerRequiringMultiplePhosphorylations(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 5;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_serineThreonineKinase) = 1e6;
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(this.MG_070_MONOMER) = 1;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(0,   m.unmodifiedMonomers(this.MG_070_MONOMER));
            assertEqual(1,   m.modifiedMonomers(this.MG_070_MONOMER));
            assertEqual(0,   m.substrates(m.substrateIndexs_atp));
            assertEqual(5,   m.substrates(m.substrateIndexs_adp));
            assertEqual(5,   m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_serineThreonineKinase));
        end

        function testMonomerRequiringMultiplePhosphorylations_oneTooFewATP(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 4;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_serineThreonineKinase) = 1;
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(this.MG_070_MONOMER) = 1;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(1, m.unmodifiedMonomers(this.MG_070_MONOMER));
            assertEqual(0, m.modifiedMonomers(this.MG_070_MONOMER));
            assertEqual(4, m.substrates(m.substrateIndexs_atp));
            assertEqual(0, m.substrates(m.substrateIndexs_adp));
            assertEqual(0, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1, m.enzymes(m.enzymeIndexs_serineThreonineKinase));
        end

        function testGlutamateAdduction(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 4;
            m.substrates(m.substrateIndexs_glutamate) = 4;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_glutamateLigase) = 4e4;
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(this.MG_166_MONOMER) = 1;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(0, m.unmodifiedMonomers(this.MG_166_MONOMER));
            assertEqual(1, m.modifiedMonomers(this.MG_166_MONOMER));
            assertEqual(0, m.substrates(m.substrateIndexs_atp));
            assertEqual(0, m.substrates(m.substrateIndexs_glutamate));
            assertEqual(4, m.substrates(m.substrateIndexs_adp));
            assertEqual(4, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(4, m.substrates(m.substrateIndexs_phosphate));
            assertEqual(4e4, m.enzymes(m.enzymeIndexs_glutamateLigase));
        end

        function testGlutamateAdduction_oneTooFewSubstrate(this)
            m = this.process;
            substrates = [
                m.substrateIndexs_atp
                m.substrateIndexs_glutamate];
            for i = 1:numel(substrates)
                m.substrates(:) = 0;
                m.substrates(m.substrateIndexs_atp) = 4;
                m.substrates(substrates) = 4;
                m.substrates(substrates(i)) = m.substrates(substrates(i)) - 1;
                m.enzymes(:) = 0;
                m.enzymes(m.enzymeIndexs_glutamateLigase) = 4e4;
                m.unmodifiedMonomers(:) = 0;
                m.unmodifiedMonomers(this.MG_166_MONOMER) = 1;
                m.modifiedMonomers(:) = 0;

                m.evolveState();
                assertEqual(1, m.unmodifiedMonomers(this.MG_166_MONOMER));
                assertEqual(0, m.modifiedMonomers(this.MG_166_MONOMER));
                assertEqual(...
                    4 * ones(numel(substrates)-1), ...
                    m.substrates(substrates([1:i-1 i+1:end])));
                assertEqual(3, m.substrates(substrates(i)));
                assertEqual(0, m.substrates(m.substrateIndexs_adp));
                assertEqual(0, m.substrates(m.substrateIndexs_hydrogen));
                assertEqual(0, m.substrates(m.substrateIndexs_phosphate));
                assertEqual(4e4, m.enzymes(m.enzymeIndexs_glutamateLigase));
            end
        end

        function testGlutamateAdduction_insufficientEnzyme(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 4;
            m.substrates(m.substrateIndexs_glutamate) = 4;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_glutamateLigase) = 100;  %(probability is low)
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(this.MG_166_MONOMER) = 1;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(1, m.unmodifiedMonomers(this.MG_166_MONOMER));
            assertEqual(0, m.modifiedMonomers(this.MG_166_MONOMER));
            assertEqual(4, m.substrates(m.substrateIndexs_atp));
            assertEqual(4, m.substrates(m.substrateIndexs_glutamate));
            assertEqual(0, m.substrates(m.substrateIndexs_adp));
            assertEqual(0, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(0, m.substrates(m.substrateIndexs_phosphate));
            assertEqual(100, m.enzymes(m.enzymeIndexs_glutamateLigase));
        end

        function testLipoateLigation(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_lipoylAmp) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lipoylTransferase) = 1;
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(this.MG_272_MONOMER) = 1;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(0, m.unmodifiedMonomers(this.MG_272_MONOMER));
            assertEqual(1, m.modifiedMonomers(this.MG_272_MONOMER));
            assertEqual(0, m.substrates(m.substrateIndexs_lipoylAmp));
            assertEqual(2, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1, m.substrates(m.substrateIndexs_amp));
            assertEqual(0, m.substrates(m.substrateIndexs_lipoylLys));
            assertEqual(1, m.enzymes(m.enzymeIndexs_lipoylTransferase));
        end

        function testLipoateLigation_noLipoylAMP(this)
            m = this.process;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_lipoylTransferase) = 1;
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(this.MG_272_MONOMER) = 1;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(1, m.unmodifiedMonomers(this.MG_272_MONOMER));
            assertEqual(0, m.modifiedMonomers(this.MG_272_MONOMER));
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(1, m.enzymes(m.enzymeIndexs_lipoylTransferase));
        end

        function testLipoateLigation_insufficientEnzyme(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_lipoylAmp) = 1;
            m.enzymes(:) = 0;
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(this.MG_272_MONOMER) = 1;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(1, m.unmodifiedMonomers(this.MG_272_MONOMER));
            assertEqual(0, m.modifiedMonomers(this.MG_272_MONOMER));
            assertEqual(1, m.substrates(m.substrateIndexs_lipoylAmp));
            assertEqual(zeros(size(m.enzymes)), m.enzymes);
        end

        % 61 total phosphorylations, 1 glutamate adduction, 1 lipoate ligation
        % Reactants:
        %    ATP: 61 + 4
        %    GLU: 4
        %    LipoylAMP: 1
        % Products:
        %    ADP: 61 + 4
        %    H: 61 + 4 + 2
        %    PI: 4
        %    AMP: 1
        function testOneOfEveryMonomer(this)
            m = this.process;

            m.substrates = sum(max(0,-m.reactionStoichiometryMatrix),2);
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_serineThreonineKinase) = 1e6;
            m.enzymes(m.enzymeIndexs_glutamateLigase) = 4e4;
            m.enzymes(m.enzymeIndexs_lipoylTransferase) = 1;
            m.unmodifiedMonomers(:) = 1;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(zeros(size(m.unmodifiedMonomers)), m.unmodifiedMonomers);
            assertEqual(ones(size(m.modifiedMonomers)), m.modifiedMonomers);
            assertEqual(sum(max(0,m.reactionStoichiometryMatrix),2),m.substrates);
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_serineThreonineKinase));
            assertEqual(4e4, m.enzymes(m.enzymeIndexs_glutamateLigase));
            assertEqual(1, m.enzymes(m.enzymeIndexs_lipoylTransferase));
        end

        % Start with a lot of each of two kinds of monomers requiring
        % phosphorylation and insufficient ATP to modify all of both kinds, then
        % verify that we modify them in quantities proportional to how many
        % phosphorylations each requires.
        function testFairAllocationOfResources(this)
            numATP = 800;
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = numATP;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_serineThreonineKinase) = 1e10;
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(this.MG_217_MONOMER) = 1000;
            m.unmodifiedMonomers(this.MG_070_MONOMER) = 1000;
            m.modifiedMonomers(:) = 0;

            m.evolveState();
            assertEqual(0, m.substrates(m.substrateIndexs_atp));
            assertEqual(...
                numATP,...
                m.modifiedMonomers(this.MG_217_MONOMER) + ...
                m.modifiedMonomers(this.MG_070_MONOMER) * 5);
            assertTrue(...
                abs(m.modifiedMonomers(this.MG_217_MONOMER) - ...
                    m.modifiedMonomers(this.MG_070_MONOMER) * 5) < ...
                    0.20 * numATP, ...
                sprintf('significant monomer imbalance: %d, %d', ...
                    m.modifiedMonomers(this.MG_217_MONOMER), ...
                    m.modifiedMonomers(this.MG_070_MONOMER) * 5));
        end

        function testGeneEssentiality(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 65;
            m.substrates(m.substrateIndexs_glutamate) = 4;
            m.substrates(m.substrateIndexs_lipoylAmp) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_serineThreonineKinase) = 1e6;
            m.enzymes(m.enzymeIndexs_glutamateLigase) = 4e4;
            m.enzymes(m.enzymeIndexs_lipoylTransferase) = 1;
            
            indices = [
                this.MG_001_MONOMER
                this.MG_070_MONOMER
                this.MG_166_MONOMER
                this.MG_217_MONOMER
                this.MG_272_MONOMER];
            m.unmodifiedMonomers(:) = 0;
            m.unmodifiedMonomers(indices) = 1;
            m.modifiedMonomers(:) = 0;

            this.helpTestGeneEssentiality({
                'MG_012';     %alpha-L-glutamate ligase
                'MG_109';     %serine/threonine protein kinase, putative
                'MG_270'},... %lipoyltransferase/lipoate-protein ligase, putative
                @(m,i) all(m.modifiedMonomers(indices) - i.modifiedMonomers(indices))); 
        end
    end
end
