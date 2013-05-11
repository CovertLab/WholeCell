% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/30/2010
classdef ProteinTranslocation_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = ProteinTranslocation_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end

    %fixtures
    methods
        function loadSimpleFixture(this)
            m = this.process;

            %IDs
            m.compartmentWholeCellModelIDs = {'c';'e';'m';'tc';'tm'};
            m.substrateWholeCellModelIDs = {'ATP';'GTP';'ADP';'GDP';'H';'H2O';'PI'};
            m.enzymeWholeCellModelIDs = {
                'MG_0001_048';                   %signal recognition particle
                'MG_297_MONOMER';                %signal recognition particle receptor
                'MG_072_DIMER';                  %preprotein translocase, SecA subunit
                'MG_055_170_277_464_476_20MER'}; %preprotein translocase
            m.monomerWholeCellModelIDs = {'1';'2';'3';'4'};

            %names
            m.substrateNames = m.substrateWholeCellModelIDs;
            m.enzymeNames    = m.enzymeWholeCellModelIDs;

            %indices
            m.compartment.cytosolIndexs       = find(strcmp(m.compartmentWholeCellModelIDs, 'c'), 1, 'first');
            m.compartment.membraneIndexs      = find(strcmp(m.compartmentWholeCellModelIDs, 'm'), 1, 'first');
            m.compartment.extracellularIndexs = find(strcmp(m.compartmentWholeCellModelIDs, 'e'), 1, 'first');

            m.substrateIndexs_atp       = 1;
            m.substrateIndexs_gtp       = 2;
            m.substrateIndexs_adp       = 3;
            m.substrateIndexs_gdp       = 4;
            m.substrateIndexs_hydrogen  = 5;
            m.substrateIndexs_water     = 6;
            m.substrateIndexs_phosphate = 7;

            m.enzymeIndexs_signalRecognitionParticle         = 1;
            m.enzymeIndexs_signalRecognitionParticleReceptor = 2;
            m.enzymeIndexs_translocaseATPase                 = 3;
            m.enzymeIndexs_translocasePore                   = 4;

            %rates
            m.translocaseSpecificRate = 2.710e12;
            m.preproteinTranslocase_aaTranslocatedPerATP    = 35;
            m.SRP_GTPUsedPerMonomer = 2;

            %weights
            m.substrateMolecularWeights = [503.1489 519.148 424.1769 440.1760 1.0079 18.0152 95.9793]';
            m.enzymeMolecularWeights    = 1e6 * (1:length(m.enzymeWholeCellModelIDs))';
            m.monomerMolecularWeights   = 1e5 * (1:length(m.monomerWholeCellModelIDs))';

            %other physical properties
            m.monomerCompartments  = [m.compartment.membraneIndexs; m.compartment.membraneIndexs; m.compartment.extracellularIndexs; m.compartment.cytosolIndexs];
            m.monomerLengths       = [100 200 150 125]';
            m.matureMonomerLengths = m.monomerLengths;
            m.monomerSRPPathways   = [true;false;false;false];

            m.monomerIndexs_translocating    = find(m.monomerCompartments~=m.compartment.cytosolIndexs);
            m.monomerIndexs_nontranslocating = find(m.monomerCompartments==m.compartment.cytosolIndexs);

            %initial state
            m.substrates   = zeros(length(m.substrateWholeCellModelIDs), 1);
            m.enzymes      = zeros(length(m.enzymeWholeCellModelIDs),    1);
            m.boundEnzymes = zeros(length(m.enzymeWholeCellModelIDs),    1);
            m.monomers     = zeros(length(m.monomerWholeCellModelIDs),   length(m.compartmentWholeCellModelIDs));
        end
    end

    %tests
    methods
        function testOneProteinThatDoesNotRequireTranslocation(this)
            m = this.process;

            m.monomerCompartments(:) = m.compartment.cytosolIndexs;
            m.monomerSRPPathways(:)=false;
            m.monomerIndexs_nontranslocating = 1:length(m.monomerCompartments);
            m.monomerIndexs_translocating = [];

            m.monomers(:) = 0;
            m.monomers(:, m.compartment.cytosolIndexs) = 1;

            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;

            initial_substrates = m.substrates;
            initial_monomers = m.monomers;
            m.evolveState();

            assertEqual(initial_monomers, m.monomers);
            assertEqual(initial_substrates, m.substrates);
        end

        function testOneTranslocationToMembrane_IntegralMembrane(this)
            m = this.process;

            idx=1;
            m.monomerCompartments(:) = m.compartment.membraneIndexs;
            m.monomerSRPPathways(:)=false;
            m.monomerSRPPathways(1)=true;
            m.monomerIndexs_translocating = (1:length(m.monomerCompartments))';
            m.monomerIndexs_nontranslocating = [];

            m.monomers(:) = 0;
            m.monomers(idx, m.compartment.cytosolIndexs) = 1;

            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp)   = ceil(m.monomerLengths(idx)/m.preproteinTranslocase_aaTranslocatedPerATP);
            m.substrates(m.substrateIndexs_gtp)   = m.SRP_GTPUsedPerMonomer;
            m.substrates(m.substrateIndexs_water) = sum(m.substrates([m.substrateIndexs_atp m.substrateIndexs_gtp]));
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalRecognitionParticle)         = 1;
            m.enzymes(m.enzymeIndexs_signalRecognitionParticleReceptor) = 1;
            m.enzymes(m.enzymeIndexs_translocaseATPase) = ceil(m.monomerLengths(idx)/m.preproteinTranslocase_aaTranslocatedPerATP)/m.translocaseSpecificRate+sqrt(eps);
            m.enzymes(m.enzymeIndexs_translocasePore) = m.enzymes(m.enzymeIndexs_translocaseATPase);

            initial_substrates = m.substrates;
            m.evolveState();

            monomers=zeros(size(m.monomers(idx, :)));
            monomers(m.compartment.membraneIndexs)=1;
            assertEqual(monomers, m.monomers(idx, :));
            assertEqual(0, m.substrates(m.substrateIndexs_atp));
            assertEqual(0, m.substrates(m.substrateIndexs_gtp));
            assertEqual(0, m.substrates(m.substrateIndexs_water));
            assertEqual(initial_substrates(m.substrateIndexs_atp), m.substrates(m.substrateIndexs_adp));
            assertEqual(initial_substrates(m.substrateIndexs_gtp), m.substrates(m.substrateIndexs_gdp));
            assertEqual(initial_substrates(m.substrateIndexs_water), m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(initial_substrates(m.substrateIndexs_water), m.substrates(m.substrateIndexs_phosphate));
        end

        function testOneTranslocationToMembrane_Lipoprotein(this)
            m = this.process;

            idx=1;
            m.monomerCompartments(:) = m.compartment.cytosolIndexs;
            m.monomerCompartments(1) = m.compartment.membraneIndexs;
            m.monomerSRPPathways(:)=false;
            m.monomerSRPPathways(1)=false;
            m.monomerIndexs_translocating = 1;
            m.monomerIndexs_nontranslocating = (2:length(m.monomerCompartments))';

            m.monomers(:) = 0;
            m.monomers(idx, m.compartment.cytosolIndexs) = 1;

            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp)   = ceil(m.monomerLengths(idx)/m.preproteinTranslocase_aaTranslocatedPerATP);
            m.substrates(m.substrateIndexs_gtp)   = 0;
            m.substrates(m.substrateIndexs_water) = sum(m.substrates([m.substrateIndexs_atp m.substrateIndexs_gtp]));
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalRecognitionParticle)         = 0;
            m.enzymes(m.enzymeIndexs_signalRecognitionParticleReceptor) = 0;
            m.enzymes(m.enzymeIndexs_translocaseATPase) = ceil(m.monomerLengths(idx)/m.preproteinTranslocase_aaTranslocatedPerATP)/m.translocaseSpecificRate+sqrt(eps);
            m.enzymes(m.enzymeIndexs_translocasePore) = m.enzymes(m.enzymeIndexs_translocaseATPase);

            initial_substrates = m.substrates;
            m.evolveState();

            monomers=zeros(size(m.monomers(idx, :)));
            monomers(m.compartment.membraneIndexs)=1;
            assertEqual(monomers, m.monomers(idx, :));
            assertEqual(0, m.substrates(m.substrateIndexs_atp));
            assertEqual(0, m.substrates(m.substrateIndexs_gtp));
            assertEqual(0, m.substrates(m.substrateIndexs_water));
            assertEqual(initial_substrates(m.substrateIndexs_atp), m.substrates(m.substrateIndexs_adp));
            assertEqual(initial_substrates(m.substrateIndexs_gtp), m.substrates(m.substrateIndexs_gdp));
            assertEqual(initial_substrates(m.substrateIndexs_water), m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(initial_substrates(m.substrateIndexs_water), m.substrates(m.substrateIndexs_phosphate));
        end

        function testOneTranslocationToExterior(this)
            m = this.process;

            idx=1;
            m.monomerCompartments(:) = m.compartment.cytosolIndexs;
            m.monomerCompartments(1) = m.compartment.extracellularIndexs;
            m.monomerSRPPathways(:)=false;
            m.monomerSRPPathways(1)=false;
            m.monomerIndexs_translocating = 1;
            m.monomerIndexs_nontranslocating = (2:length(m.monomerCompartments))';

            m.monomers(:) = 0;
            m.monomers(idx, m.compartment.cytosolIndexs) = 1;

            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp)   = ceil(m.monomerLengths(idx)/m.preproteinTranslocase_aaTranslocatedPerATP);
            m.substrates(m.substrateIndexs_gtp)   = 0;
            m.substrates(m.substrateIndexs_water) = sum(m.substrates([m.substrateIndexs_atp m.substrateIndexs_gtp]));
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalRecognitionParticle)         = 0;
            m.enzymes(m.enzymeIndexs_signalRecognitionParticleReceptor) = 0;
            m.enzymes(m.enzymeIndexs_translocaseATPase) = ceil(m.monomerLengths(idx)/m.preproteinTranslocase_aaTranslocatedPerATP)/m.translocaseSpecificRate+sqrt(eps);
            m.enzymes(m.enzymeIndexs_translocasePore) = m.enzymes(m.enzymeIndexs_translocaseATPase);

            initial_substrates = m.substrates;
            m.evolveState();

            monomers=zeros(size(m.monomers(idx, :)));
            monomers(m.compartment.extracellularIndexs)=1;
            assertEqual(monomers, m.monomers(idx, :));
            assertEqual(0, m.substrates(m.substrateIndexs_atp));
            assertEqual(0, m.substrates(m.substrateIndexs_gtp));
            assertEqual(0, m.substrates(m.substrateIndexs_water));
            assertEqual(initial_substrates(m.substrateIndexs_atp), m.substrates(m.substrateIndexs_adp));
            assertEqual(initial_substrates(m.substrateIndexs_gtp), m.substrates(m.substrateIndexs_gdp));
            assertEqual(initial_substrates(m.substrateIndexs_water), m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(initial_substrates(m.substrateIndexs_water), m.substrates(m.substrateIndexs_phosphate));
        end

        %lots of two kinds of monomers and not enough energy to translocate
        %all of both
        function testFairness(this)
            %process
            m = this.process;

            %setup network
            m.monomerCompartments(:) = m.compartment.cytosolIndexs;
            m.monomerCompartments(1:2) = m.compartment.membraneIndexs;
            m.monomerSRPPathways(:)=false;
            m.monomerSRPPathways(1:2)=true;
            m.monomerIndexs_translocating = (1:2)';
            m.monomerIndexs_nontranslocating = (3:length(m.monomerCompartments))';
            m.monomerLengths(:) = 500;

            %lots of all monomers
            nMonomers = 10000;
            m.monomers(:) = 0;
            m.monomers(:, m.compartment.cytosolIndexs) = nMonomers;

            %sufficient substrates and enzymes to translocate all proteins
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp)   = sum(ceil(m.monomerLengths(m.monomerIndexs_translocating)/m.preproteinTranslocase_aaTranslocatedPerATP));
            m.substrates(m.substrateIndexs_gtp)   = sum(m.monomerSRPPathways) * m.SRP_GTPUsedPerMonomer;
            m.substrates(m.substrateIndexs_water) = sum(m.substrates([m.substrateIndexs_atp m.substrateIndexs_gtp]));

            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalRecognitionParticle)         = length(m.monomerIndexs_translocating);
            m.enzymes(m.enzymeIndexs_signalRecognitionParticleReceptor) = m.enzymes(m.enzymeIndexs_signalRecognitionParticle);
            m.enzymes(m.enzymeIndexs_translocaseATPase) = sum(ceil(m.monomerLengths(m.monomerIndexs_translocating)/m.preproteinTranslocase_aaTranslocatedPerATP)/m.translocaseSpecificRate)+sqrt(eps);
            m.enzymes(m.enzymeIndexs_translocasePore) = m.enzymes(m.enzymeIndexs_translocaseATPase);

            %cut back substrates and enzymes so insufficient among to
            %translocate everything
            m.substrates = nMonomers * 0.6 * m.substrates;
            m.enzymes    = nMonomers * 0.6 * m.enzymes;

            %evolve state
            m.evolveState();

            %assert that energy, enzymes allocated fairly
            monomers = m.monomers(1:2, m.compartment.membraneIndexs);
            assertTrue( abs(diff(monomers)) < max(monomers)/100);
        end

        %one of every kind of protein can be processed in one time step
        %- membrane, lipoprotein, extracellular
        function testAllProteinsCanBeProcessed(this)
            m = this.process;

            m.monomerCompartments(1)=m.compartment.membraneIndexs;
            m.monomerCompartments(2)=m.compartment.membraneIndexs;
            m.monomerCompartments(3)=m.compartment.extracellularIndexs;
            m.monomerCompartments(4:end) = m.compartment.cytosolIndexs;
            m.monomerSRPPathways(1)=true;
            m.monomerSRPPathways(2)=false;
            m.monomerSRPPathways(3)=false;
            m.monomerSRPPathways(4:end)=false;
            m.monomerIndexs_translocating = (1:3)';
            m.monomerIndexs_nontranslocating = (4:length(m.monomerCompartments))';

            m.monomers(:) = 0;
            m.monomers(:, m.compartment.cytosolIndexs) = 1;

            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp)   = sum(ceil(m.monomerLengths(m.monomerIndexs_translocating)/m.preproteinTranslocase_aaTranslocatedPerATP));
            m.substrates(m.substrateIndexs_gtp)   = sum(m.monomerSRPPathways) * m.SRP_GTPUsedPerMonomer;
            m.substrates(m.substrateIndexs_water) = sum(m.substrates([m.substrateIndexs_atp m.substrateIndexs_gtp]));
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalRecognitionParticle)         = length(m.monomerIndexs_translocating);
            m.enzymes(m.enzymeIndexs_signalRecognitionParticleReceptor) = m.enzymes(m.enzymeIndexs_signalRecognitionParticle);
            m.enzymes(m.enzymeIndexs_translocaseATPase) = sum(ceil(m.monomerLengths(m.monomerIndexs_translocating)/m.preproteinTranslocase_aaTranslocatedPerATP)/m.translocaseSpecificRate)+sqrt(eps);
            m.enzymes(m.enzymeIndexs_translocasePore) = m.enzymes(m.enzymeIndexs_translocaseATPase);

            initial_substrates = m.substrates;
            m.evolveState();

            monomers=zeros(size(m.monomers));
            monomers(sub2ind(size(monomers), (1:length(m.monomerCompartments))', m.monomerCompartments))=1;
            assertEqual(monomers, m.monomers);
            assertEqual(0, m.substrates(m.substrateIndexs_atp));
            assertEqual(0, m.substrates(m.substrateIndexs_gtp));
            assertEqual(0, m.substrates(m.substrateIndexs_water));
            assertEqual(initial_substrates(m.substrateIndexs_atp), m.substrates(m.substrateIndexs_adp));
            assertEqual(initial_substrates(m.substrateIndexs_gtp), m.substrates(m.substrateIndexs_gdp));
            assertEqual(initial_substrates(m.substrateIndexs_water), m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(initial_substrates(m.substrateIndexs_water), m.substrates(m.substrateIndexs_phosphate));
        end

        function testGeneEssentiality(this)
            m = this.process;

            m.monomers(:)=0;
            m.monomers(:,m.compartment.cytosolIndexs) = 1;

            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            m.substrates(m.substrateIndexs_gtp) = 1e6;
            m.substrates(m.substrateIndexs_water) = 1e6;

            m.substrates(:)=1e10;

            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_signalRecognitionParticle) = 1e6;
            m.enzymes(m.enzymeIndexs_signalRecognitionParticleReceptor) = 1e6;
            m.enzymes(m.enzymeIndexs_translocaseATPase) = 1e6;
            m.enzymes(m.enzymeIndexs_translocasePore) = 1e6;

            m.enzymes(:) = 1e10;

            this.helpTestGeneEssentiality({
                'MG_048';     %signal recognition particle protein
                'MG_055';     %preprotein translocase, SecE subunit
                'MG_072';     %preprotein translocase, SecA subunit
                'MG_476';     %preprotein translocase, SecG subunit
                'MG_170';     %preprotein translocase, SecY subunit
                'MG_0001';    %scRNA, signal recognition particle 4.5S RNA, MCS1
                'MG_277';     %Protein-export membrane protein secDF
                'MG_297';     %signal recognition particle receptor
                'MG_464'},... %inner-membrane protein insertion factor
                @(m,~) ~any(m.monomers(m.monomerIndexs_translocating, m.compartment.cytosolIndexs)));
        end
    end
end
