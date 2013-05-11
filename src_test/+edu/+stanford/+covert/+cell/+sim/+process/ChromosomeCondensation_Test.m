%ChromosomeCondensation process test case
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/9/2010
classdef ChromosomeCondensation_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = ChromosomeCondensation_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end
    
    %fixtures
    methods
        function loadSimpleTestFixture(this)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            import edu.stanford.covert.cell.sim.state.Chromosome;      
            
            %% process
            m = this.process;
            
            %% parameters
            m.enzymeDNAFootprints(m.enzymeIndexs_SMC_ADP) = 630;
            m.smcSepNt = 7130;     %1 SMC complex per x bp; Jensen, Shapiro 2003
            m.smcSepProbCenter = 2800;
            
            %% substrates
            m.substrateWholeCellModelIDs = {
                'ATP'; 'ADP'; 'PI'; 'H2O'; 'H'};
            m.substrateNames = m.substrateWholeCellModelIDs;
            
            m.substrateIndexs_atp       = 1;
            m.substrateIndexs_adp       = 2;
            m.substrateIndexs_phosphate = 3;
            m.substrateIndexs_water     = 4;
            m.substrateIndexs_hydrogen  = 5;
            
            m.substrateMetaboliteLocalIndexs = (1:numel(m.substrateWholeCellModelIDs))';
            m.substrateMetaboliteGlobalIndexs = m.substrateMetaboliteLocalIndexs;
            
            m.substrateMolecularWeights = [503.1489; 424.1769; 95.9793; 18.0152; 1.0079];
            
            %% enzymes
            m.enzymeWholeCellModelIDs = {
                'MG_213_214_298_6MER';      %Chromosome Segregation Protein SMC with SCP Proteins
                'MG_213_214_298_6MER_ADP';  %Chromosome Segregation Protein SMC with SCP Proteins-ADP
                };
            m.enzymeNames = m.enzymeWholeCellModelIDs;
            
            m.enzymeIndexs_SMC     = 1;
            m.enzymeIndexs_SMC_ADP = 2;
            
            m.enzymeMonomerLocalIndexs = zeros(0, 1);
            m.enzymeComplexLocalIndexs = (1:numel(m.enzymeWholeCellModelIDs))';
            m.enzymeMonomerGlobalIndexs = m.enzymeMonomerLocalIndexs;
            m.enzymeComplexGlobalIndexs = m.enzymeComplexLocalIndexs;
            
            m.enzymeMolecularWeights = [378486.032; 378910.208900];
            
            %% chromosome state            
            c = Chromosome([], []);
            
            c.relaxedBasesPerTurn = 10.5;
            c.equilibriumSuperhelicalDensity = -0.06;
            c.supercoiledSuperhelicalDensityTolerance = 0.1;
            
            bases = 'ACGT';
            seq = bases(randi(4, 580076, 1));
            c.sequence = ChromosomeSequence(seq);
            c.sequenceLen = numel(seq);
                                    
            m.states = {c};
            m.chromosome = c;
                        
            c.metabolite.molecularWeights = [503.1489; 424.1769; 95.9793; 18.0152; 1.0079; 329.2055; 305.1808; 345.2049; 320.1921; 212.0942; 149.1530];
            c.metabolite.dnmpIndexs = (6:9)';
            c.metabolite.waterIndexs = 4;
            c.metabolite.dr5pIndexs = 10;
            c.metabolite.m6ADIndexs = 11;
            
            %protein release reactions
            c.monomerDNAFootprints = ones(size(m.enzymeMonomerLocalIndexs));
            c.complexDNAFootprints = ones(size(m.enzymeComplexLocalIndexs));
            c.complexDNAFootprints(m.enzymeToComplex(m.enzymeIndexs_SMC_ADP)) = 630;
            
            c.monomerDNAFootprintBindingStrandedness = repmat(c.dnaStrandedness_dsDNA, size(m.enzymeMonomerLocalIndexs));
            c.complexDNAFootprintBindingStrandedness = repmat(c.dnaStrandedness_dsDNA, size(m.enzymeComplexLocalIndexs));
            c.monomerDNAFootprintRegionStrandedness = repmat(c.dnaStrandedness_dsDNA, size(m.enzymeMonomerLocalIndexs));
            c.complexDNAFootprintRegionStrandedness = repmat(c.dnaStrandedness_dsDNA, size(m.enzymeComplexLocalIndexs));
            
            nReactions = 0;
            c.reactionBoundMonomer = zeros(nReactions, 1);
            c.reactionBoundComplex = zeros(nReactions, 1);
            c.reactionMonomerCatalysisMatrix = zeros(nReactions, numel(m.enzymeMonomerGlobalIndexs));
            c.reactionComplexCatalysisMatrix = zeros(nReactions, numel(m.enzymeComplexGlobalIndexs));
            c.reactionThresholds = ...
                sum(c.reactionMonomerCatalysisMatrix, 2) + ...
                sum(c.reactionComplexCatalysisMatrix, 2);
            
            %% initial state
            m.substrates    = zeros(length(m.substrateWholeCellModelIDs), 1);
            m.enzymes       = zeros(length(m.enzymeWholeCellModelIDs),    1);
            m.boundEnzymes  = zeros(length(m.enzymeWholeCellModelIDs),    1);
            
            m.enzymes(m.enzymeIndexs_SMC) = 1e6;
            m.substrates(m.substrateIndexs_atp) = m.enzymes(m.enzymeIndexs_SMC);
            m.substrates(m.substrateIndexs_water) = m.enzymes(m.enzymeIndexs_SMC);
                        
            c.initialize();
        end
    end

    %tests
    methods
        function testConstants(this)
            m = this.process;
            
            %check all enzymes are complexes
            assertEqual((1:numel(m.enzymeComplexLocalIndexs))', m.enzymeComplexLocalIndexs);
            
            %check SMC complex DNA footprint is 630 nt
            assertEqual(630, m.enzymeDNAFootprints(m.enzymeIndexs_SMC_ADP));
            
            %smcSepProbCenter = 2800 nt
            assertEqual(2800, m.smcSepProbCenter);
            
            %smcDNAFootprint, smcSepProbCenter are even
            assertTrue(iseven(m.enzymeDNAFootprints(m.enzymeIndexs_SMC_ADP)));
            assertTrue(iseven(m.smcSepProbCenter));            
            assertIn(m.smcSepNt + m.smcSepProbCenter, [m.enzymeDNAFootprints(m.enzymeIndexs_SMC_ADP) Inf]);
        end
                
        function testBindSeveralSMCs(this)
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 10;
            m.substrates(m.substrateIndexs_water) = 10;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_SMC) = 10;
            m.boundEnzymes(:) = 0;
            
            initial_polymerizedRegions = c.polymerizedRegions;
            
            m.evolveState();
            
            assertEqual(0, m.substrates(m.substrateIndexs_atp));
            assertEqual(0, m.substrates(m.substrateIndexs_water));
            assertEqual(0, m.substrates(m.substrateIndexs_adp));
            assertEqual(10, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(10, m.substrates(m.substrateIndexs_phosphate));            
            
            assertEqual(0, nnz(m.enzymes));
            assertEqual(0, m.boundEnzymes(m.enzymeIndexs_SMC));
            assertEqual(10, m.boundEnzymes(m.enzymeIndexs_SMC_ADP));
            assertEqual(0, nnz(c.monomerBoundSites));
            assertEqual(10, nnz(c.complexBoundSites));
            assertEqual(0, nnz(c.getDamagedSites(true, true, true, true, false, false, true)));
            assertEqual(initial_polymerizedRegions, c.polymerizedRegions);            
        end
        
        function testInitializeStateConverges(this)
            m = this.process;
            c = m.chromosome;
            
            c.initialize();         
            
            m.substrates(:) = 0;            
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_SMC) = 1e6;            
            m.substrates(m.substrateIndexs_atp) = m.enzymes(m.enzymeIndexs_SMC);
            m.substrates(m.substrateIndexs_water) = m.enzymes(m.enzymeIndexs_SMC);
            
            m.initializeState();
            assertIn(m.boundEnzymes(m.enzymeIndexs_SMC_ADP), [0.9 1.1] * size(c.sequence, 1) / m.smcSepNt);
            
            boundSMCs = m.boundEnzymes(m.enzymeIndexs_SMC_ADP);
            for i = 1:1000
                m.evolveState();
            end
            assertElementsAlmostEqual(boundSMCs, m.boundEnzymes(m.enzymeIndexs_SMC_ADP), 'relative', 0.05);
        end
        
        function testBindSMCs_NoRoomForMoreSMCs(this)
            m = this.process;
            c = m.chromosome;            
            
            smcADPIdx = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_SMC_ADP);            
            
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.enzymes(m.enzymeIndexs_SMC_ADP) = 0;
            m.boundEnzymes(:) = 0;
            
            c.complexBoundSites(1:m.smcSepNt:end, 1) = smcADPIdx;
            m.boundEnzymes(m.enzymeIndexs_SMC_ADP) = nnz(c.complexBoundSites);
            
            m.evolveState();
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_SMC));
                
            for i = 1:1000
                m.evolveState();                
            end
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_SMC));
        end

        function testBindSMCs_RoomForOneMoreSMC(this)
            m = this.process;
            c = m.chromosome;
            c.initialize();      

            smcADPIdx = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_SMC_ADP);

            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.enzymes(m.enzymeIndexs_SMC_ADP) = 0;
            m.boundEnzymes(:) = 0;

            c.complexBoundSites(1:m.smcSepNt:end-2*m.smcSepNt, 1) = smcADPIdx;
            m.boundEnzymes(m.enzymeIndexs_SMC_ADP) = nnz(c.complexBoundSites);

            m.evolveState();
            assertEqual(1e6-1, m.enzymes(m.enzymeIndexs_SMC));

            for i = 1:10
                m.evolveState();
            end 
            assertEqual(1e6-1, m.enzymes(m.enzymeIndexs_SMC));
        end
        
        function testBindSMCs_NoATP(this)            
            m = this.process;
            c = m.chromosome;
            c.initialize();
            
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_atp) = 0;
            m.enzymes(:) = 1e6;
            m.enzymes(m.enzymeIndexs_SMC_ADP) = 0;
            m.boundEnzymes(:) = 0;
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            initial_polymerizedRegions = c.polymerizedRegions;
            
            m.evolveState();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(0, nnz(c.monomerBoundSites));
            assertEqual(0, nnz(c.complexBoundSites));            
            assertEqual(initial_polymerizedRegions, c.polymerizedRegions);
        end
        
        function testBindSMCs_NoWater(this)            
            m = this.process;
            c = m.chromosome;
            c.initialize();
            
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_water) = 0;
            m.enzymes(:) = 1e6;
            m.enzymes(m.enzymeIndexs_SMC_ADP) = 0;
            m.boundEnzymes(:) = 0;
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            initial_polymerizedRegions = c.polymerizedRegions;
            
            m.evolveState();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(0, nnz(c.monomerBoundSites));
            assertEqual(0, nnz(c.complexBoundSites));            
            assertEqual(initial_polymerizedRegions, c.polymerizedRegions);
        end
      
        function testBindSMCs_NoEnzyme(this)
            m = this.process;
            c = m.chromosome;
            c.initialize();
            
            m.substrates(:) = 1e6;
            m.enzymes(:) = 1e6;
            m.enzymes(m.enzymeIndexs_SMC) = 0;
            m.enzymes(m.enzymeIndexs_SMC_ADP) = 0;
            m.boundEnzymes(:) = 0;
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            initial_polymerizedRegions = c.polymerizedRegions;
            
            m.evolveState();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(0, nnz(c.monomerBoundSites));
            assertEqual(0, nnz(c.complexBoundSites));
            assertEqual(initial_polymerizedRegions, c.polymerizedRegions);
        end
        
        function testFreeSMCADPDissociation(this)
            m = this.process;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_SMC_ADP) = 10;
            m.boundEnzymes(:) = 0;
            
            m.evolveState();
            
            assertEqual(10, m.substrates(m.substrateIndexs_adp));
            assertEqual(1, nnz(m.substrates));
            assertEqual(10, m.enzymes(m.enzymeIndexs_SMC));
            assertEqual(1, nnz(m.enzymes));
            assertAllEqual(0, m.boundEnzymes);
        end
        
        % Gene essentiality
        function testGeneEssentiality(this)
            m = this.process;
            c = this.process.chromosome;
            
            c.initialize();            
            m.enzymes(:) = 1e6;
            m.boundEnzymes(:) = 0;
            m.substrates(:) = 1e6;

            this.helpTestGeneEssentiality({
                'MG_213';      %segregation and condensation protein A
                'MG_214';      %segregation and condensation protein B
                'MG_298'},...  %chromosome segregation protein SMC
                @this.isProperlyFunctioning,...
                struct('lengthSec', 20));
        end
    end

    %helper methods
    methods (Access = private)
        function result = isProperlyFunctioning(~, m, ~)
            nExpectedBound = size(m.chromosome.sequence, 1) / m.smcSepNt;
            result = ...
                0.9 * nExpectedBound < m.boundEnzymes(m.enzymeIndexs_SMC_ADP) & ...
                1.1 * nExpectedBound > m.boundEnzymes(m.enzymeIndexs_SMC_ADP);
        end
    end
end
