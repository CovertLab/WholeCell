%Replication Initiation process test case
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/12/2010
classdef ReplicationInitiation_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = ReplicationInitiation_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end
    end
    
    methods
        function loadTestFixture(this)
            this.loadTestFixture@edu.stanford.covert.cell.sim.ProcessTestCase();
            
            m = this.process;
            c = m.chromosome;
            c.initialize();
        end
        
        function loadSimpleTestFixture(this)
            %import classes
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            import edu.stanford.covert.cell.sim.state.Chromosome;
            
            %% process
            m = this.process;
            
            %% parameters
            m.enzymeDNAFootprints(m.enzymeIndexs_DnaA_1mer_ATP)        = 11; %nt
            
            m.siteCooperativity    = 70;       %
            m.stateCooperativity   = 20;       %
            
            m.kb1ATP               = 25;       %nM/h
            m.kb2ATP               = 0.6;      %nM/h
            m.kd1ATP               = 20;       %1/h
            m.kb1ADP               = 2.5;      %nM/h
            m.kb2ADP               = 0.61;     %nM/h
            m.kd1ADP               = 20;       %1/h
            m.k_Regen              = 2.315e-6; %1/s
            m.K_Regen_P4           = 0.018;    %g/L
            m.k_inact              = 4.24e14;  %1/s
            
            %% metabolites
            m.substrateWholeCellModelIDs = {'ATP'; 'ADP'; 'PI'; 'H2O'; 'H'};
            m.substrateNames = m.substrateWholeCellModelIDs;
            
            m.substrateIndexs_atp        = 1;
            m.substrateIndexs_adp        = 2;
            m.substrateIndexs_phosphate  = 3;
            m.substrateIndexs_water      = 4;
            m.substrateIndexs_hydrogen   = 5;
            
            m.substrateMetaboliteLocalIndexs = (1:numel(m.substrateWholeCellModelIDs))';
            m.substrateMetaboliteGlobalIndexs = m.substrateMetaboliteLocalIndexs;
            
            m.substrateMolecularWeights = [503.1489; 424.1769; 95.9793; 18.0152; 1.0079];
            
            %% enzymes
            m.enzymeWholeCellModelIDs = {
                'MG_469_1MER_ADP'
                'MG_469_1MER_ATP'
                'MG_469_2MER_1ATP_ADP'
                'MG_469_2MER_ATP'
                'MG_469_3MER_2ATP_ADP'
                'MG_469_3MER_ATP'
                'MG_469_4MER_3ATP_ADP'
                'MG_469_4MER_ATP'
                'MG_469_5MER_4ATP_ADP'
                'MG_469_5MER_ATP'
                'MG_469_6MER_5ATP_ADP'
                'MG_469_6MER_ATP'
                'MG_469_7MER_6ATP_ADP'
                'MG_469_7MER_ATP'
                'MG_469_MONOMER'
                };
            m.enzymeNames = m.enzymeWholeCellModelIDs;
            m.enzymeIndexs_DnaA             = 15;
            m.enzymeIndexs_DnaA_1mer_ADP    = 1;
            m.enzymeIndexs_DnaA_1mer_ATP    = 2;
            m.enzymeIndexs_DnaA_Nmer_ATP    = (2:2:14)';
            m.enzymeIndexs_DnaA_Nmer_ADP    = (1:2:13)';
            m.enzymeIndexs_DnaA_polymer_ATP = (4:2:14)';
            m.enzymeIndexs_DnaA_polymer_ADP = (3:2:13)';
            
            m.enzymeMonomerLocalIndexs = m.enzymeIndexs_DnaA;
            m.enzymeMonomerGlobalIndexs = 1;
            m.enzymeComplexLocalIndexs = (1:numel(m.enzymeWholeCellModelIDs)-1)';
            m.enzymeComplexGlobalIndexs = m.enzymeComplexLocalIndexs;
            
            m.enzymeMolecularWeights = zeros(size(m.enzymeWholeCellModelIDs));
            m.enzymeMolecularWeights(m.enzymeIndexs_DnaA) = 50766.528;
            m.enzymeMolecularWeights(m.enzymeIndexs_DnaA_Nmer_ATP) = (1:7)' * (m.enzymeMolecularWeights(m.enzymeIndexs_DnaA) + m.substrateMolecularWeights(m.substrateIndexs_atp));
            m.enzymeMolecularWeights(m.enzymeIndexs_DnaA_Nmer_ADP) = (1:7)' * m.enzymeMolecularWeights(m.enzymeIndexs_DnaA) + (0:6)' * m.substrateMolecularWeights(m.substrateIndexs_atp) + m.substrateMolecularWeights(m.substrateIndexs_adp);
            
            %% chromosome
            c = Chromosome([], []);
            m.states = {c};
            m.chromosome = c;
            
            %sequence
            bases = 'ACGT';
            seq = bases(randi(4, 580076, 1));
            c.sequence = ChromosomeSequence(seq);
            c.sequenceLen = length(seq);
            c.sequenceGCContent = getGCContent(c.sequence);
            
            %metabolites
            c.metabolite.molecularWeights = [503.1489; 424.1769; 95.9793; 18.0152; 1.0079; 329.2055; 305.1808; 345.2049; 320.1921; 212.0942; 149.1530];
            c.metabolite.dnmpIndexs = (6:9)';
            c.metabolite.waterIndexs = 4;
            c.metabolite.dr5pIndexs = 10;
            c.metabolite.m6ADIndexs = 11;
            
            %protein footprints
            c.monomerDNAFootprints = ones(size(m.enzymeMonomerLocalIndexs));
            c.complexDNAFootprints = [
                repmat(m.enzymeDNAFootprints(m.enzymeIndexs_DnaA_1mer_ATP), size(m.enzymeComplexLocalIndexs));
                20];
            
            c.monomerDNAFootprintBindingStrandedness = repmat(c.dnaStrandedness_dsDNA, size(m.enzymeMonomerLocalIndexs));
            c.complexDNAFootprintBindingStrandedness = [
                repmat(c.dnaStrandedness_dsDNA, size(m.enzymeComplexLocalIndexs));
                c.dnaStrandedness_ssDNA];
            c.monomerDNAFootprintRegionStrandedness = repmat(c.dnaStrandedness_dsDNA, size(m.enzymeMonomerLocalIndexs));
            c.complexDNAFootprintRegionStrandedness = [
                repmat(c.dnaStrandedness_dsDNA, size(m.enzymeComplexLocalIndexs));
                c.dnaStrandedness_xsDNA];
            
            %protein release reactions -- helicase can release all kinds of DnaA
            %complexes
            nReactions = numel(m.enzymeComplexLocalIndexs);
            c.reactionBoundMonomer = zeros(nReactions, 1);
            c.reactionBoundComplex = (1:numel(m.enzymeComplexLocalIndexs))';
            c.reactionMonomerCatalysisMatrix = zeros(nReactions, numel(m.enzymeMonomerGlobalIndexs));
            c.reactionComplexCatalysisMatrix = zeros(nReactions, numel(m.enzymeComplexGlobalIndexs)+1);
            c.reactionComplexCatalysisMatrix(:, end) = 1;
            c.reactionThresholds = ...
                sum(c.reactionMonomerCatalysisMatrix, 2) + ...
                sum(c.reactionComplexCatalysisMatrix, 2);
            
            %% DnaA boxes
            nBoxes = 2283;
            
            m.dnaABoxIndexs_9mer = [1:148 nBoxes - 3]';
            m.dnaABoxIndexs_7mer = nBoxes - 4;
            m.dnaABoxIndexs_8mer = [149:nBoxes-5 nBoxes-2 nBoxes-1 nBoxes]';
            m.dnaABoxIndexs_R12345 = (nBoxes:-1:nBoxes-4)';
            
            m.dnaABoxStartPositions = ceil((1+(1:nBoxes)') * size(c.sequence, 1) / (nBoxes+100));
            
            m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345) = [
                579011 + ceil(8/2 - m.enzymeDNAFootprints(m.enzymeIndexs_DnaA_1mer_ATP)/2);
                578924 + ceil(8/2 - m.enzymeDNAFootprints(m.enzymeIndexs_DnaA_1mer_ATP)/2);
                578896 + ceil(8/2 - m.enzymeDNAFootprints(m.enzymeIndexs_DnaA_1mer_ATP)/2);
                578792 + ceil(9/2 - m.enzymeDNAFootprints(m.enzymeIndexs_DnaA_1mer_ATP)/2);  
                578597 + ceil(7/2 - m.enzymeDNAFootprints(m.enzymeIndexs_DnaA_1mer_ATP)/2)]; 
            
            %% initial state
            m.geometry.volume               = 2e-21;
            
            m.metabolite.counts(:, sim.membraneCompartmentIndexs) = 0;
            m.metabolite.counts(1, sim.membraneCompartmentIndexs) = 3.15e-15 / m.metabolite.molecularWeights(1) * ...
                edu.stanford.covert.util.ConstantUtil.nAvogadro;
            
            c.initialize();
            c.segregated = false;
            
            m.substrates   = zeros(size(m.substrateWholeCellModelIDs));
            m.enzymes      = zeros(size(m.enzymeWholeCellModelIDs));
            m.boundEnzymes = zeros(size(m.enzymeWholeCellModelIDs));
        end
    end
    
    % tests
    methods
        function testConstants(this)
            m = this.process;
            c = m.chromosome;
            
            %assert that indices of complexes into enzymeComplexLocalIndexs are
            %the same as their local indexs so that for example
            %  this.enzymeComplexGlobalIndexs(this.enzymeComplexLocalIndexs == this.enzymeIndexs_DnaA_1mer_ATP)
            %is equivalent to
            %  this.enzymeComplexGlobalIndexs(this.enzymeIndexs_DnaA_1mer_ATP)
            %and thus we can use the 2nd shorter syntax
            assertEqual((1:numel(m.enzymeComplexLocalIndexs))', m.enzymeComplexLocalIndexs);
            
            %assert that all DnaA-AxP polymers have same DNA footprint size
            assertEqual(1, numel(unique(c.complexDNAFootprints(m.enzymeComplexGlobalIndexs([m.enzymeIndexs_DnaA_Nmer_ATP; m.enzymeIndexs_DnaA_Nmer_ADP])))));
            
            %assert DnaA polymer molecular weights computed correctly
            enzymeMolecularWeights = zeros(size(m.enzymeWholeCellModelIDs));
            enzymeMolecularWeights(m.enzymeIndexs_DnaA) = m.enzymeMolecularWeights(m.enzymeIndexs_DnaA);
            enzymeMolecularWeights(m.enzymeIndexs_DnaA_Nmer_ATP) = (1:7)' * (m.enzymeMolecularWeights(m.enzymeIndexs_DnaA) + m.substrateMolecularWeights(m.substrateIndexs_atp));
            enzymeMolecularWeights(m.enzymeIndexs_DnaA_Nmer_ADP) = (1:7)' * m.enzymeMolecularWeights(m.enzymeIndexs_DnaA) + (0:6)' * m.substrateMolecularWeights(m.substrateIndexs_atp) + m.substrateMolecularWeights(m.substrateIndexs_adp);
            assertElementsAlmostEqual(enzymeMolecularWeights', m.enzymeMolecularWeights');
            
            %assert that DnaA boxes don't overlap
            overlapIdxs = zeros(0, 2);
            for i = 1:numel(m.dnaABoxStartPositions)-1
                tmp = find(m.dnaABoxStartPositions(i) == m.dnaABoxStartPositions(i+1:end));
                overlapIdxs = [overlapIdxs; repmat(i, size(tmp)) i+tmp]; %#ok<AGROW>
            end
            assertTrue(isempty(overlapIdxs), sprintf('DnaA boxes cannot overlap. Overlapping DnaA boxes:\n\tBox 1\tBox 2\n\t=====\t=====\n%s', sprintf('\t%d\t%d\n', overlapIdxs)));
            
            %assert positions of R1-5 functional DnaA boxes
            assertEqual([
                579011 + ceil(8/2 - m.enzymeDNAFootprints(m.enzymeIndexs_DnaA_1mer_ATP)/2);
                578924 + ceil(8/2 - m.enzymeDNAFootprints(m.enzymeIndexs_DnaA_1mer_ATP)/2);
                578896 + ceil(8/2 - m.enzymeDNAFootprints(m.enzymeIndexs_DnaA_1mer_ATP)/2);
                578792 + ceil(9/2 - m.enzymeDNAFootprints(m.enzymeIndexs_DnaA_1mer_ATP)/2);
                578597 + ceil(7/2 - m.enzymeDNAFootprints(m.enzymeIndexs_DnaA_1mer_ATP)/2)], ...
                m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345));
            
            %asert kd1ATP = kd1ADP
            assertEqual(m.kd1ATP, m.kd1ADP);
            
            %assert enzymeComplexGlobalIndexs sorted allowing ismembc to be used
            %in place of ismember
            assertTrue(issorted(m.enzymeComplexGlobalIndexs));
            assertTrue(issorted(m.enzymeGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP)));
            assertTrue(issorted(m.enzymeGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ADP)));
        end
        
        function testInitializedToSteadyState(this)
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            
            totalDnaA = ...
                + m.enzymes(m.enzymeIndexs_DnaA) ...
                + (1:7) * m.enzymes(m.enzymeIndexs_DnaA_Nmer_ATP) ...
                + (1:7) * m.enzymes(m.enzymeIndexs_DnaA_Nmer_ADP) ...
                + (1:7) * m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ATP) ...
                + (1:7) * m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ADP);
            
            %% initialize state
            m.initializeState();
            
            %% assertions
            %assert conservation of DnaA monomers
            assertElementsAlmostEqual(totalDnaA, ...
                + m.enzymes(m.enzymeIndexs_DnaA) ...
                + (1:7) * m.enzymes(m.enzymeIndexs_DnaA_Nmer_ATP) ...
                + (1:7) * m.enzymes(m.enzymeIndexs_DnaA_Nmer_ADP) ...
                + (1:7) * m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ATP) ...
                + (1:7) * m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ADP));
            
            %assert no free DnaA, DnaA bound
            assertAllEqual(0, m.enzymes(setdiff(1:end, m.enzymeIndexs_DnaA_1mer_ATP)));
            assertIn(m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP), [0 1]);
            assertEqual(1, nnz(m.boundEnzymes));
            assertEqual(totalDnaA, m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP) + m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            
            status = m.calculateDnaABoxStatus();
            ssBound8 = sum(status(m.dnaABoxIndexs_8mer, 1) == m.dnaABoxStatus_DnaAATPBound);
            
            %% evolve state
            for i = 1:100
                m.evolveState();
                
                assertAllEqual(0, m.enzymes(m.enzymeIndexs_DnaA));
                assertAllEqual(0, m.enzymes(m.enzymeIndexs_DnaA_Nmer_ATP(2:end)));
                assertAllEqual(true, 1 >= m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ATP(2:end)));
                assertAllEqual(0, m.enzymes(m.enzymeIndexs_DnaA_Nmer_ADP));
                assertAllEqual(0, m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ADP));
            end
            
            %% assertions
            %assert conservation of DnaA monomers
            assertElementsAlmostEqual(totalDnaA, ...
                + m.enzymes(m.enzymeIndexs_DnaA) ...
                + (1:7) * m.enzymes(m.enzymeIndexs_DnaA_Nmer_ATP) ...
                + (1:7) * m.enzymes(m.enzymeIndexs_DnaA_Nmer_ADP) ...
                + (1:7) * m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ATP) ...
                + (1:7) * m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ADP));
            
            %assert no free DnaA, DnaA bound
            status = m.calculateDnaABoxStatus();
            bound8 = sum(status(m.dnaABoxIndexs_8mer, 1) == m.dnaABoxStatus_DnaAATPBound);
            
            assertElementsAlmostEqual(0, m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP), 'absolute', 3);
            assertElementsAlmostEqual(bound8, ssBound8, 'relative', 0.3);
        end
                
        function testDnaAAtpBindingTo9mer(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA) = 1;
            m.boundEnzymes(:) = 0;
            
            m.kb2ATP = 0;  %disable DnaA-ATP 8mer binding
            m.kd1ATP = 0;  %disable DnaA-ATP unbinding
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.k_inact = 0; %disable DnaA-ATP inactivation
            
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertEqual(1, sum(dnaABoxStatus(m.dnaABoxIndexs_9mer, 1) == m.dnaABoxStatus_DnaAATPBound));
            assertEqual(size(dnaABoxStatus, 1) -1, sum(dnaABoxStatus(:, 1) == m.dnaABoxStatus_NotBound));
            assertEqual(size(dnaABoxStatus, 1), sum(dnaABoxStatus(:, 2) == m.dnaABoxStatus_NotExist));
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(zeros(size(m.enzymes)), m.enzymes);
            assertEqual(1, sum(abs(m.boundEnzymes)));
            assertEqual(1, m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ATP));
        end
        
        function testDnaAAtpBindingTo9mer_LimitedByDnaA(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA) = 0;
            m.boundEnzymes(:) = 0;
            
            m.kb2ATP = 0;  %disable DnaA-ATP 8mer binding
            m.kd1ATP = 0;  %disable DnaA-ATP unbinding
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.k_inact = 0; %disable DnaA-ATP inactivation
            
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertEqual(0, sum(dnaABoxStatus(m.dnaABoxIndexs_9mer, 1) == m.dnaABoxStatus_DnaAATPBound));
            assertEqual(size(dnaABoxStatus, 1), sum(dnaABoxStatus(:, 1) == m.dnaABoxStatus_NotBound));
            assertEqual(size(dnaABoxStatus, 1), sum(dnaABoxStatus(:, 2) == m.dnaABoxStatus_NotExist));
            assertEqual(1, m.substrates(m.substrateIndexs_atp));
            assertEqual(1, sum(abs(m.substrates)));
            assertEqual(zeros(size(m.enzymes)), m.enzymes);
            assertEqual(zeros(size(m.boundEnzymes)), m.boundEnzymes);
        end
        
        function testDnaAAtpBindingTo9mer_LimitedByATP(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA) = 1;
            m.boundEnzymes(:) = 0;
            
            m.kb2ATP = 0;  %disable DnaA-ATP 8mer binding
            m.kd1ATP = 0;  %disable DnaA-ATP unbinding
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.k_inact = 0; %disable DnaA-ATP inactivation
            
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertEqual(0, sum(dnaABoxStatus(m.dnaABoxIndexs_9mer, 1) == m.dnaABoxStatus_DnaAATPBound));
            assertEqual(size(dnaABoxStatus, 1), sum(dnaABoxStatus(:, 1) == m.dnaABoxStatus_NotBound));
            assertEqual(size(dnaABoxStatus, 1), sum(dnaABoxStatus(:, 2) == m.dnaABoxStatus_NotExist));
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(1, m.enzymes(m.enzymeIndexs_DnaA));
            assertEqual(1, sum(abs(m.enzymes)));
            assertEqual(zeros(size(m.boundEnzymes)), m.boundEnzymes);
        end
        
        function testDnaAAtpBindingTo8mer(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA) = 1;
            m.boundEnzymes(:) = 0;
            
            m.kb1ATP = 0;  %disable DnaA-ATP 9mer binding
            m.kd1ATP = 0;  %disable DnaA-ATP unbinding
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.k_inact = 0; %disable DnaA-ATP inactivation
            
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertEqual(1, sum(dnaABoxStatus(m.dnaABoxIndexs_8mer, 1) == m.dnaABoxStatus_DnaAATPBound));
            assertEqual(size(dnaABoxStatus, 1) -1, sum(dnaABoxStatus(:, 1) == m.dnaABoxStatus_NotBound));
            assertEqual(size(dnaABoxStatus, 1), sum(dnaABoxStatus(:, 2) == m.dnaABoxStatus_NotExist));
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(zeros(size(m.enzymes)), m.enzymes);
            assertEqual(1, m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            assertEqual(1, sum(abs(m.boundEnzymes)));
        end
        
        function testDnaAAtpBindingInOriCRegion(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions, 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ATP);
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), 1) = 0;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA) = 1;
            m.boundEnzymes(:) = 0;
            
            m.kb1ATP = 1e6;
            m.kb2ATP = 0;  %disable DnaA-ATP 8mer binding
            m.kd1ATP = 0;  %disable DnaA-ATP unbinding
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.k_inact = 0; %disable DnaA-ATP inactivation
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertEqual(1, sum(sum(dnaABoxStatus(m.dnaABoxIndexs_9mer, :) == m.dnaABoxStatus_NotBound)));
            assertEqual(3, sum(sum(dnaABoxStatus(m.dnaABoxIndexs_8mer, :) == m.dnaABoxStatus_NotBound)));
            assertEqual(1, sum(sum(dnaABoxStatus(m.dnaABoxIndexs_7mer, :) == m.dnaABoxStatus_NotBound)));
            
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertEqual(1, sum(dnaABoxStatus(m.dnaABoxIndexs_R12345, 1) == m.dnaABoxStatus_DnaAATPBound));
            assertEqual(4, sum(dnaABoxStatus(:, 1) == m.dnaABoxStatus_NotBound));
            assertEqual(size(dnaABoxStatus, 1), sum(dnaABoxStatus(:, 2) == m.dnaABoxStatus_NotExist));
            
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(zeros(size(m.enzymes)), m.enzymes);
            assertEqual(1, m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            assertEqual(1, sum(abs(m.boundEnzymes)));
        end
        
        function testDnaAAtpUnbinding(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions, 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ATP);
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ATP) = size(m.dnaABoxStartPositions, 1);
            
            m.kb1ATP = 0;  %disable DnaA-ATP 9mer binding
            m.kb2ATP = 0;  %disable DnaA-ATP 8mer binding
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.kd1ATP = 1e3; %boost DnaA-ATP unbinding
            m.kd1ADP = 1e3; %boost DnaA-ATP unbinding
            
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertFalse(any(m.enzymes(setdiff(1:end, m.enzymeIndexs_DnaA_1mer_ATP))));
            assertFalse(any(m.boundEnzymes(setdiff(1:end, m.enzymeIndexs_DnaA_1mer_ATP))));
            assertTrue(m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP) > 0);
            assertEqual(sum(dnaABoxStatus(:, 1) == m.dnaABoxStatus_NotBound), m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            assertEqual(sum(dnaABoxStatus(:, 1) == m.dnaABoxStatus_DnaAATPBound), m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ATP))
            assertElementsAlmostEqual(...
                sum(dnaABoxStatus(m.dnaABoxIndexs_9mer, 1) == m.dnaABoxStatus_NotBound) / sum(dnaABoxStatus(m.dnaABoxIndexs_8mer, 1) == m.dnaABoxStatus_NotBound), ...
                numel(m.dnaABoxIndexs_9mer) / numel(m.dnaABoxIndexs_8mer), ...
                'relative', 0.2);
        end
        
        function testDnaAAtpPolymerUnbinding(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions, 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ATP);
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R1234(1:3)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(2));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R1234(4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ATP) = size(m.dnaABoxStartPositions, 1) - 3;
            m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ATP(2)) = 3;
            
            m.kb1ATP = 0;  %disable DnaA-ATP 9mer binding
            m.kb2ATP = 0;  %disable DnaA-ATP 8mer binding
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.kd1ATP = 1e6; %boost DnaA-ATP unbinding
            m.kd1ADP = 1e6; %boost DnaA-ATP unbinding
            
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertAllEqual(0, m.substrates);
            assertFalse(any(m.enzymes(setdiff(1:end, m.enzymeIndexs_DnaA_1mer_ATP))));
            assertFalse(any(m.boundEnzymes(setdiff(1:end, m.enzymeIndexs_DnaA_1mer_ATP))));
            assertEqual(numel(m.dnaABoxStartPositions)-1, m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            assertEqual(4, m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            assertAllEqual(m.dnaABoxStatus_NotBound, dnaABoxStatus(setdiff(1:end, m.dnaABoxIndexs_R1234), 1));
            assertAllEqual(m.dnaABoxStatus_DnaAATPBound, dnaABoxStatus(m.dnaABoxIndexs_R1234, 1));
        end
        
        function testDnaAAtpComplexNotUnbinding(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions, 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ATP);
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R1234), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(2));
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ATP) = size(m.dnaABoxStartPositions, 1) - 4;
            m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ATP(2)) = 4;
            
            m.kb1ATP = 0;  %disable DnaA-ATP 9mer binding
            m.kb2ATP = 0;  %disable DnaA-ATP 8mer binding
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.kd1ATP = 1e6; %boost DnaA-ATP unbinding
            m.kd1ADP = 1e6; %boost DnaA-ATP unbinding
            
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertAllEqual(0, m.substrates);
            assertFalse(any(m.enzymes(setdiff(1:end, m.enzymeIndexs_DnaA_1mer_ATP))));
            assertFalse(any(m.boundEnzymes(setdiff(1:end, m.enzymeIndexs_DnaA_Nmer_ATP(2)))));
            assertEqual(numel(m.dnaABoxStartPositions)-4, m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            assertEqual(4, m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ATP(2)));
            assertAllEqual(m.dnaABoxStatus_NotBound, dnaABoxStatus(setdiff(1:end, m.dnaABoxIndexs_R1234), 1));
            assertAllEqual(m.dnaABoxStatus_DnaAATPBound, dnaABoxStatus(m.dnaABoxIndexs_R1234, 1));
        end
        
        function testDnaAAtpCompletedComplexNotUnbinding(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions, 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ATP);
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R1234), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ATP) = size(m.dnaABoxStartPositions, 1) - 4;
            m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ATP(7)) = 4;
            
            m.kb1ATP = 0;  %disable DnaA-ATP 9mer binding
            m.kb2ATP = 0;  %disable DnaA-ATP 8mer binding
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.kd1ATP = 1e6; %boost DnaA-ATP unbinding
            m.kd1ADP = 1e6; %boost DnaA-ATP unbinding
            
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertAllEqual(0, m.substrates);
            assertFalse(any(m.enzymes(setdiff(1:end, m.enzymeIndexs_DnaA_1mer_ATP))));
            assertFalse(any(m.boundEnzymes(setdiff(1:end, m.enzymeIndexs_DnaA_Nmer_ATP([1 7])))));
            assertEqual(numel(m.dnaABoxStartPositions)-5, m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            assertEqual(1, m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ATP(1)));
            assertEqual(4, m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ATP(7)));
            assertAllEqual(m.dnaABoxStatus_NotBound, dnaABoxStatus(setdiff(1:end, m.dnaABoxIndexs_R12345), 1));
            assertAllEqual(m.dnaABoxStatus_DnaAATPBound, dnaABoxStatus(m.dnaABoxIndexs_R12345, 1));
        end
        
        function testDnaAAtpUnbindingInOriCRegion(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345([1 4])), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ATP);
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ATP) = 2;
            
            m.kb1ATP = 0;  %disable DnaA-ATP 9mer binding
            m.kb2ATP = 0;  %disable DnaA-ATP 8mer binding
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.kd1ATP = 1e4; %boost DnaA-ATP unbinding
            m.kd1ADP = 1e4; %boost DnaA-ATP unbinding
            
            m.evolveState();
            
            assertEqual([0; 0], c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345([1 4])) ones(2, 1)]));
            
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertFalse(any(m.enzymes(setdiff(1:end, m.enzymeIndexs_DnaA_1mer_ATP))));
            assertEqual(2, m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            assertEqual(zeros(size(m.boundEnzymes)), m.boundEnzymes);
        end
        
        function testDnaAAdpBindingTo9mer(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA_1mer_ADP) = 1;
            m.boundEnzymes(:) = 0;
            
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.kd1ATP = 0;
            m.kd1ADP = 0;
            m.k_Regen = 0;
            m.K_Regen_P4 = 0;
            
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertEqual(1, sum(dnaABoxStatus(m.dnaABoxIndexs_9mer, 1) == m.dnaABoxStatus_DnaAADPBound));
            assertEqual(size(dnaABoxStatus, 1) -1, sum(dnaABoxStatus(:, 1) == m.dnaABoxStatus_NotBound));
            assertEqual(size(dnaABoxStatus, 1), sum(dnaABoxStatus(:, 2) == m.dnaABoxStatus_NotExist));
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(zeros(size(m.enzymes)), m.enzymes);
            assertEqual(1, m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ADP));
            assertEqual(1, sum(abs(m.boundEnzymes)));
        end
        
        function testDnaAAdpBindingTo8mer(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA_1mer_ADP) = 1;
            m.boundEnzymes(:) = 0;
            
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kd1ATP = 0;
            m.kd1ADP = 0;
            m.k_Regen = 0;
            m.K_Regen_P4 = 0;
            
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertEqual(1, sum(dnaABoxStatus(m.dnaABoxIndexs_8mer, 1) == m.dnaABoxStatus_DnaAADPBound));
            assertEqual(size(dnaABoxStatus, 1) -1, sum(dnaABoxStatus(:, 1) == m.dnaABoxStatus_NotBound));
            assertEqual(size(dnaABoxStatus, 1), sum(dnaABoxStatus(:, 2) == m.dnaABoxStatus_NotExist));
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(zeros(size(m.enzymes)), m.enzymes);
            assertEqual(1, m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ADP));
            assertEqual(1, sum(abs(m.boundEnzymes)));
        end
        
        function testDnaAAdpBindingInOriCRegion(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions, 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ADP);
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), 1) = 0;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA_1mer_ADP) = 1;
            m.boundEnzymes(:) = 0;
            
            m.kb1ADP = 0;   %disable DnaA-ADP 9mer binding
            m.kb2ADP = 1e8; %increase binding rate
            m.kd1ATP = 0;
            m.kd1ADP = 0;
            m.k_Regen = 0;
            m.K_Regen_P4 = 0;
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertEqual(1, sum(sum(dnaABoxStatus(m.dnaABoxIndexs_9mer, :) == m.dnaABoxStatus_NotBound)));
            assertEqual(3, sum(sum(dnaABoxStatus(m.dnaABoxIndexs_8mer, :) == m.dnaABoxStatus_NotBound)));
            assertEqual(1, sum(sum(dnaABoxStatus(m.dnaABoxIndexs_7mer, :) == m.dnaABoxStatus_NotBound)));
            
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertEqual(1, sum(dnaABoxStatus(m.dnaABoxIndexs_R12345, 1) == m.dnaABoxStatus_DnaAADPBound));
            assertEqual(4, sum(dnaABoxStatus(:, 1) == m.dnaABoxStatus_NotBound));
            assertEqual(size(dnaABoxStatus, 1), sum(dnaABoxStatus(:, 2) == m.dnaABoxStatus_NotExist));
            
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(zeros(size(m.enzymes)), m.enzymes);
            assertEqual(1, m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ADP));
            assertEqual(1, sum(abs(m.boundEnzymes)));
        end
        
        function testDnaAAdpUnbinding(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions, 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ADP);
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ADP) = size(m.dnaABoxStartPositions, 1);
            
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kd1ATP = 1e3;
            m.kd1ADP = 1e3;
            m.k_Regen = 0;
            m.K_Regen_P4 = 0;
                        
            m.evolveState();
            
            dnaABoxStatus = m.calculateDnaABoxStatus();
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertTrue(m.enzymes(m.enzymeIndexs_DnaA_1mer_ADP) > 0);
            assertFalse(any(m.enzymes(setdiff(1:end, m.enzymeIndexs_DnaA_1mer_ADP))));
            assertEqual(size(m.dnaABoxStartPositions, 1), m.enzymes(m.enzymeIndexs_DnaA_1mer_ADP) + m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ADP));
            assertFalse(any(m.boundEnzymes(setdiff(1:end, m.enzymeIndexs_DnaA_1mer_ADP))));
            assertEqual(sum(dnaABoxStatus(:, 1) == m.dnaABoxStatus_NotBound), m.enzymes(m.enzymeIndexs_DnaA_1mer_ADP));
            assertEqual(sum(dnaABoxStatus(:, 1) == m.dnaABoxStatus_DnaAADPBound), m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ADP))
            assertElementsAlmostEqual(...
                sum(dnaABoxStatus(m.dnaABoxIndexs_9mer, 1) == m.dnaABoxStatus_NotBound) / sum(dnaABoxStatus(m.dnaABoxIndexs_8mer, 1) == m.dnaABoxStatus_NotBound), ...
                numel(m.dnaABoxIndexs_9mer) / numel(m.dnaABoxIndexs_8mer), ...
                'relative', 0.2);
        end
        
        function testDnaARegeneration(this)
            m = this.process;
            c = m.chromosome;
            
            m.kb1ADP = 0;   %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;   %disable DnaA-ADP 8mer binding
            m.kb1ATP = 0;   %disable DnaA-ATP 9mer binding
            m.kb2ATP = 0;   %disable DnaA-ATP 8mer binding
            m.k_Regen = 1e6;
            m.K_Regen_P4 = 0;
            
            c.complexBoundSites(:, :) = 0;
            
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA_1mer_ADP) = 1e6;
            m.boundEnzymes(:) = 0;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e6;
            
            m.evolveState();
            
            assertEqual(1e6, m.substrates(m.substrateIndexs_adp));
            assertFalse(any(m.substrates(setdiff(1:end, m.substrateIndexs_adp))));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            assertFalse(any(m.enzymes(setdiff(1:end, m.enzymeIndexs_DnaA_1mer_ATP))));
            assertEqual(zeros(size(m.boundEnzymes)), m.boundEnzymes);
        end
        
        function testDnaAComplexGrowth(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions, 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ATP);
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(3));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2:4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(2));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(5)), 1) = 0;
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1e4;
            m.enzymes(m.enzymeIndexs_DnaA) = 1e4;
            
            m.kd1ATP = 0;   %disable DnaA-ATP unbinding
            m.kb1ADP = 0;
            m.kb2ADP = 0;
            m.k_inact = 0;
            
            m.evolveState();
            
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertAllEqual(3, polATPs(:,1));
            assertAllEqual(0, polATPs(:,2));
            assertAllEqual(0, polADPs);
            assertEqual([3 0], m.calculateDnaAR1234ComplexSize(polATPs, polADPs));
            assertEqual([false false], m.calculateIsDnaAORIComplexAssembled());
        end
        
        function testNoDnaAComplexGrowthWithDnaAAdpBoundToRSite(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions, 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ATP);
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345([1 2 4])), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ATP);
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(3)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ADP);
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(5)), 1) = 0;
            
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP) = 1e4;
            
            m.kd1ATP = 0;   %disable DnaA-ATP unbinding
            m.kd1ADP = 0;   %disable DnaA-ADP unbinding
            m.k_inact = 0;
            
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([1; 1; 0; 1], polATPs(:,1));
            assertAllEqual(0, polATPs(:,2));
            assertEqual([0; 0; 1; 0], polADPs(:,1));
            assertAllEqual(0, polADPs(:,2));
            assertEqual([0 0], m.calculateDnaAR1234ComplexSize(polATPs, polADPs));
            
            m.evolveState();
            
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([1; 1; 0; 1], polATPs(:,1));
            assertAllEqual(0, polATPs(:,2));
            assertEqual([0; 0; 1; 0], polADPs(:,1));
            assertAllEqual(0, polADPs(:,2));
            assertEqual([0 0], m.calculateDnaAR1234ComplexSize(polATPs, polADPs));
            assertEqual([false false], m.calculateIsDnaAORIComplexAssembled());
        end
        
        % Verifies that a DnaA-ATP will bind R5 when the complex is 7 deep,
        % and that replication initiation will then be triggered.
        function testDnaAComplexCompletion(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1:4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA) = 1;
            m.boundEnzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ATP) = 28;
                       
            m.kb1ATP = 0;   %disable DnaA-ATP 9mer binding
            m.kb2ATP = 0;   %disable DnaA-ATP 8mer binding
            
            m.evolveState();
            
            assertEqual(0, m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            assertEqual(29, m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            assertEqual([true false], m.calcuateIsDnaAR5Occupied, 'site R5 not occupied');
            
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertAllEqual(7, polATPs(:,1));
            assertAllEqual(0, polATPs(:,2));
            assertAllEqual(0, polADPs);
            assertEqual([7 0], m.calculateDnaAR1234ComplexSize(polATPs, polADPs));
            assertEqual([true false], m.calculateIsDnaAORIComplexAssembled, 'replication not initiated');
        end
        
        function testDnaAComplexDissociation(this)
            m = this.process;
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 1e6;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA_Nmer_ATP(7)) = 4;
            m.enzymes(m.enzymeIndexs_DnaA_Nmer_ATP(1)) = 1;
            m.boundEnzymes(:) = 0;
            
            m.kb1ATP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ATP = 0;  %disable DnaA-ADP 8mer binding
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.k_Regen = 0; %disable DnaA-ATP regeneration
            m.k_inact = 1e21;
            
            m.evolveState();
            
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertAllEqual(0, polATPs);
            assertAllEqual(0, polADPs);
            assertEqual([0 0], m.calculateDnaAR1234ComplexSize(polATPs, polADPs));
            assertEqual([false false], m.calcuateIsDnaAR5Occupied(), 'R5 still occupied');
            
            assertEqual(28, m.substrates(m.substrateIndexs_phosphate));
            assertEqual(28, m.substrates(m.substrateIndexs_hydrogen));
            assertEqual(1e6-28, m.substrates(m.substrateIndexs_water));
            assertEqual(0, m.substrates(m.substrateIndexs_atp));
            assertEqual(0, m.substrates(m.substrateIndexs_adp));
            
            assertEqual(28, m.enzymes(m.enzymeIndexs_DnaA_1mer_ADP));
            assertEqual(1, m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP));
            assertEqual(2, nnz(m.enzymes));
            assertAllEqual(0, m.boundEnzymes);
        end
        
        function testDnaAComplexDissociation_NoWater(this)
            m = this.process;
            
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA_Nmer_ATP(7)) = 4;
            m.enzymes(m.enzymeIndexs_DnaA_Nmer_ATP(1)) = 1;
            m.boundEnzymes(:) = 0;
            
            m.kb1ATP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ATP = 0;  %disable DnaA-ADP 8mer binding
            m.kb1ADP = 0;  %disable DnaA-ADP 9mer binding
            m.kb2ADP = 0;  %disable DnaA-ADP 8mer binding
            m.k_Regen = 0; %disable DnaA-ATP regeneration
            m.k_inact = 1e21;
            
            m.evolveState();
            
            assertEqual(4, m.enzymes(m.enzymeIndexs_DnaA_Nmer_ATP(7)));
        end
        
        % Verifies that nothing changes when no enzymes are present.
        function testNoEnzymes(this)
            m = this.process;
            c = m.chromosome;
            
            m.substrates(:) = 1e6;
            m.enzymes(:) = 0;
            m.boundEnzymes(:) = 0;
            
            initial_substrates = m.substrates;
            initial_enzymes = m.enzymes;
            initial_boundEnzymes = m.boundEnzymes;
            initial_monomerBoundSites = c.monomerBoundSites;
            initial_complexBoundSites = c.complexBoundSites;
            
            m.evolveState();
            
            assertEqual(initial_substrates, m.substrates);
            assertEqual(initial_enzymes, m.enzymes);
            assertEqual(initial_boundEnzymes, m.boundEnzymes);
            assertEqual(initial_monomerBoundSites, c.monomerBoundSites);
            assertEqual(initial_complexBoundSites, c.complexBoundSites);
        end

        function testComplexFormation(this)
            m = this.process;
            c = m.chromosome;
            c.monomerBoundSites(:, :) = 0;
            c.complexBoundSites(:, :) = 0;
            
            DnaA = ...
                + (1:7) * m.enzymes(m.enzymeIndexs_DnaA_Nmer_ATP) ...
                + (1:7) * m.enzymes(m.enzymeIndexs_DnaA_Nmer_ADP) ...
                + (1:7) * m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ATP) ...
                + (1:7) * m.boundEnzymes(m.enzymeIndexs_DnaA_Nmer_ADP) ...
                + m.enzymes(m.enzymeIndexs_DnaA) ...
                + m.boundEnzymes(m.enzymeIndexs_DnaA);
            
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA) = DnaA;
            m.boundEnzymes(:) = 0; 
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 40000;
            m.substrates(m.substrateIndexs_water) = 40000;
            
            m.stateCooperativity = 1e4;
            m.siteCooperativity = 1e4;
            
            m.initializeState();
            
            assertAllEqual(true, m.calculateDnaAR1234Polymerization <= 1);
            for i = 1:5000
                m.evolveState();
                
                if any(m.calculateIsDnaAORIComplexAssembled)
                    break;
                end
            end
            
            assertEqual([true false], m.calculateIsDnaAORIComplexAssembled);
        end
        
        function testDnaAADPRegenerationFollowReplicationInitiation(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = this.process;
            c = m.chromosome;
            
            %setup just after replication initiates
            nDnaAADP = 4 * 7;
            m.substrates(:) = 1e6;
            m.enzymes = 2 * m.enzymes + m.boundEnzymes;
            m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP) = ...
                + m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP) ...
                - nDnaAADP;
            m.enzymes(m.enzymeIndexs_DnaA_Nmer_ATP(7)) = ...
                + m.enzymes(m.enzymeIndexs_DnaA_Nmer_ATP(7)) ...
                + 4;
            assertIn(m.enzymes(m.enzymeIndexs_DnaA_1mer_ATP), [0 Inf]);
            
            %simulate during replication phase
            iterMax = 5 * 3600;
            for i = 1:iterMax
                l = ceil(c.sequenceLen * i / iterMax / 2);
                m.chromosome.polymerizedRegions = CircularSparseMat(...
                    [1 1; 1 2; 1 3; 1 4; c.sequenceLen-l+1 3; c.sequenceLen-l+1 4], ...
                    [repmat(c.sequenceLen, 2, 1); repmat(l, 4, 1)], ...
                    [c.sequenceLen c.nCompartments], 1);
                m.chromosome.linkingNumbers = CircularSparseMat(...
                    [1 1; 1 2; 1 3; 1 4; c.sequenceLen-l+1 3; c.sequenceLen-l+1 4], ...
                    (1 + c.equilibriumSuperhelicalDensity) * [repmat(c.sequenceLen, 2, 1); repmat(l, 4, 1)] / c.relaxedBasesPerTurn, ...
                    [c.sequenceLen c.nCompartments], 1);
                m.evolveState();
            end
            m.chromosome.polymerizedRegions = CircularSparseMat(...
                [1 1; 1 2; 1 3; 1 4], ...
                repmat(c.sequenceLen, 4, 1), ...
                [c.sequenceLen c.nCompartments], 1);
            m.chromosome.linkingNumbers = CircularSparseMat(...
                [1 1; 1 2; 1 3; 1 4], ...
                (1 + c.equilibriumSuperhelicalDensity) * repmat(c.sequenceLen, 4, 1) / c.relaxedBasesPerTurn, ...
                [c.sequenceLen c.nCompartments], 1);
            
            assertElementsAlmostEqual(0.45^5 * nDnaAADP, ...
                + m.enzymes(m.enzymeIndexs_DnaA_1mer_ADP) ...
                + m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ADP), ...
                'relative', 0.2, 3);
            
            %simulate during cytokinesis phase
            for i = 1:3600
                m.evolveState();
            end
            assertElementsAlmostEqual(.45^6 * nDnaAADP, ...
                + m.enzymes(m.enzymeIndexs_DnaA_1mer_ADP) ...
                + m.boundEnzymes(m.enzymeIndexs_DnaA_1mer_ADP), ...
                'relative', 0.2, 2);
        end
        
        function testInactivateFreeDnaAATP(this)
            m = this.process;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA_polymer_ATP) = 10;
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = (2:7) * m.enzymes(m.enzymeIndexs_DnaA_polymer_ATP);
            
            m.inactivateFreeDnaAATP();
            
            assertAllEqual(0, m.enzymes(m.enzymeIndexs_DnaA_polymer_ATP));
            assertEqual(0, m.substrates(m.substrateIndexs_water))
        end
        
        function testInactivateFreeDnaAATP_LimitedWater(this)
            m = this.process;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA_polymer_ATP) = 10;
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_water) = 20;
            
            m.inactivateFreeDnaAATP();
            
            assertIn((2:7) * (10 - m.enzymes(m.enzymeIndexs_DnaA_polymer_ATP)), [14 20]);
            assertIn(m.substrates(m.substrateIndexs_water), [0 6])
        end
        
        function testDuration(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            %default parameter values
            m = this.process;
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load();
            g = sim.gene;
            time = sim.state('Time');
            mass = sim.state('Mass');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            
            siteCooperativity  = m.siteCooperativity;
            stateCooperativity = m.stateCooperativity;
            cellCycleLength = time.cellCycleLength;
            expectedDuration = time.replicationInitiationDuration;
            nDnaABoxes = numel(m.dnaABoxStartPositions);
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            monExp = mRNAExp ./ (log(2) / time.cellCycleLength + pm.decayRates(pm.matureIndexs));
            monExp = monExp / sum(monExp);
            meanDnaACopyNumber = mass.cellInitialDryWeight * mass.dryWeightFractionProtein * ConstantUtil.nAvogadro / ...
                (monExp' * pm.molecularWeights(pm.matureIndexs)) * monExp(m.enzymeMonomerGlobalIndexs);
            
            clear sim g time mass rna pm;
            
            %estimate mean replication initiation duration based on default
            %parameter values
            duration = edu.stanford.covert.cell.sim.process.ReplicationInitiation_Test.sampleDuration(...
                meanDnaACopyNumber, siteCooperativity, stateCooperativity, nDnaABoxes, ...
                1.5 * cellCycleLength, 100, 1, 1, [], 0);
            
            %assert that estimated mean replication initation duration is equal
            %to the target replication initiation duration
            assertElementsAlmostEqual(expectedDuration, duration, 'relative', 0.75);
        end
        
        function testGeneEssentiality(this)
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
            c.complexBoundSites(:, :) = 0;
            
            %remove all binding sites other than R1-5 sites
            m.dnaABoxStartPositions = m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345);
            m.dnaABoxIndexs_R12345 = (1:5)';
            m.dnaABoxIndexs_R1234 = m.dnaABoxIndexs_R12345(1:4,1);
            m.dnaABoxIndexs_R5    = m.dnaABoxIndexs_R12345(5,1);
            m.dnaABoxIndexs_9mer = 4;
            m.dnaABoxIndexs_8mer = (1:3)';
            m.dnaABoxIndexs_7mer = 5;
            
            %bump up the quantity of free enzyme and ATP
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_DnaA) = 1e6;
            m.substrates(m.substrateIndexs_atp) = 1e6;
                        
            this.helpTestGeneEssentiality({
                'MG_469'},... %chromosomal replication initiator protein dnaA
                @(m,i) all([true false] == m.calculateIsDnaAORIComplexAssembled()),...
                struct('lengthSec', 8));
        end
    end
    
    %test helper helper functions
    methods
        %tests 
        %- calculateDnaAR1234Polymerization
        %- calculateDnaAR1234ComplexSize
        %- calculateDnaAR1234ATPPolymerizationRates
        %- calculateDnaAR1234ATPPolymerizationCooperativity        
        function testCalculateDnaAR1234ATPPolymerization(this)
            m = this.process;
            c = m.chromosome;
            sc = 50;
            sC = 2;
            m.siteCooperativity = sc;
            m.stateCooperativity = sC;
            
            nAvo = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            r9 = m.kb1ATP * 1e9 / 3600 / nAvo / m.geometry.volume * m.stepSizeSec;
            r8 = m.kb2ATP * 1e9 / 3600 / nAvo / m.geometry.volume * m.stepSizeSec;
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([0 0 0 0; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([0 0], pol);
            assertEqual([1 1 1 1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 0 0 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            %1 DnaA-ATP 1 mer
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([1 0 0 0; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([0 0], pol);
            assertEqual([1 1 1 sc; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 0 0 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([0 1 0 0; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([0 0], pol);
            assertEqual([1 1 1 sc; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 0 0 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(3)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([0 0 1 0; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([0 0], pol);
            assertEqual([1 1 1 sc; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 0 0 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([0 0 0 1; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([0 0], pol);
            assertEqual([sc 1 1 1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 0 0 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            %multiple Dna-ATP 1mers
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345([1 4])), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([1 0 0 1; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([0 0], pol);
            assertEqual([1 sc sc 1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 0 0 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2:4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([0 1 1 1; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([0 0], pol);
            assertEqual([sc 1 1 1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 0 0 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ADP(1));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2:4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([0 1 1 1; 0 0 0 0]', polATPs);
            assertEqual([1 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([0 0], pol);
            assertEqual([1 1 1 1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 0 0 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            %complex size 1
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1:4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([1 1 1 1; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([1 0], pol);
            assertEqual([1 1 1 sC^1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([r8 r8 r8 sC^1*r9; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(2));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2:4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([2 1 1 1; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([1 0], pol);
            assertEqual([1 1 1 sc+sC*1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 r8 r8 (sc + sC*1)*r9; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(2));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1:3)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([1 1 1 2; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([1 0], pol);
            assertEqual([sc 1 1 1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([sc*r8 r8 r8 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345([1 3])), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345([2 4])), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(2));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([1 2 1 2; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([1 0], pol);
            assertEqual([sc 1 1 1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([sc*r8 0 r8 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2:4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(2));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([1 2 2 2; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([1 0], pol);
            assertEqual([sc 1 1 1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([sc*r8 0 0 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(3)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345([1 2 4])), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(2));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([2 2 1 2; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([1 0], pol);
            assertEqual([1 1 sc 1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 0 sc*r8 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345([2 3])), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345([1 4])), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(2));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([2 1 1 2; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([1 0], pol);
            assertEqual([1 sc sc 1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 sc*r8 sc*r8 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            %complex size 5
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2:4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(5));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(6));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([6 5 5 5; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([5 0], pol);
            assertEqual([1 1 1 sc+sC*5; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 r8 r8 (sc+sC*5)*r9; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(3)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ADP(5));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345([2 4])), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(5));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(6));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([6 5 0 5; 0 0 0 0]', polATPs);
            assertEqual([0 0 5 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([4 0], pol);
            assertEqual([1 1 1 1; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 0 0 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
            
            %complex size 7
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345), [1 3]) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1:4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual([7 7 7 7; 0 0 0 0]', polATPs);
            assertEqual([0 0 0 0; 0 0 0 0]', polADPs);
            polATPs = polATPs(1:4, :);
            polADPs = polADPs(1:4, :);
            pol = m.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            assertEqual([7 0], pol);
            assertEqual([1 1 1 sC*7; 1 1 1 1]', m.calculateDnaAR1234ATPPolymerizationCooperativity(polATPs, polADPs, pol));
            assertEqual([0 0 0 0; 0 0 0 0]', m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
        end
        
        function testIsDnaAORIComplexAssembled(this)
            m = this.process;
            c = m.chromosome;
            
            %ex 1
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(6));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2)) 1]) = m.enzymeComplexGlobalIndexs(end) + 1;
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(3)) 1]) = 0;
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(4)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(5)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            
            assertEqual([false false], m.calculateIsDnaAORIComplexAssembled());
            
            %ex 2
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(3)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(4)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(5)) 1]) = 0;
            
            assertEqual([false false], m.calculateIsDnaAORIComplexAssembled());
            
            %ex 3
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(3)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(4)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(5)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            
            assertEqual([true false], m.calculateIsDnaAORIComplexAssembled());
        end
        
        function testCalculateDnaAR12345Polymerization(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(6));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ADP);
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(3)) 1]) = 0;
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(4)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(5)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            
            assertEqual([[6; 0; 0; 7] zeros(4, 1)], m.calculateDnaAR1234Polymerization());
            assertEqual([true false], m.calcuateIsDnaAR5Occupied());
        end
        
        function testIsDnaAR5Occupied(this)
            m = this.process;
            c = m.chromosome;
            
            %ex 1
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(6));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2)) 1]) = m.enzymeComplexGlobalIndexs(end) + 1;
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(3)) 1]) = 0;
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(4)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(5)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            
            assertEqual([true false], m.calcuateIsDnaAR5Occupied());
            
            %ex 2
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(6));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(2));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(3)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(3));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(4)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(5)) 1]) = 0;
            
            assertEqual([false false], m.calcuateIsDnaAR5Occupied());
        end
        
        function testGetDnaABoxStatus(this)
            m = this.process;
            c = m.chromosome;
            
            c.complexBoundSites([m.dnaABoxStartPositions(1) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ADP(1));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(6));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(2)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_1mer_ADP);
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(3)) 1]) = 0;
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(4)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(7));
            c.complexBoundSites([m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(5)) 1]) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(1));
            
            dnaABoxStatus = repmat([m.dnaABoxStatus_NotBound m.dnaABoxStatus_NotExist], size(m.dnaABoxStartPositions));
            dnaABoxStatus(m.dnaABoxIndexs_R12345([1 4 5]), 1) = m.dnaABoxStatus_DnaAATPBound;
            dnaABoxStatus(m.dnaABoxIndexs_R12345(2), 1) = m.dnaABoxStatus_DnaAADPBound;
            dnaABoxStatus(1, 1) = m.dnaABoxStatus_DnaAADPBound;
            
            assertEqual(dnaABoxStatus, m.calculateDnaABoxStatus());
        end
        
        function testcalculateDnaAR1234ATPPolymerizationRates_nothingBound(this)
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
            c.complexBoundSites(:, :) = 0;
            
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertEqual(zeros(4,2), m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
        end
        
        function testcalculateDnaAR1234ATPPolymerizationRates_complexEffects(this)
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
            c.polymerizedRegions(1, :) = size(c.sequence, 1);
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1:4)), 1:2:end) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(5));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(3)), 3) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(6));
            
            Navo = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            polRates = zeros(4, 2);
            polRates(4,   :) = m.kb1ATP * 1e9 / 3600 / Navo / m.geometry.volume * m.stepSizeSec;
            polRates(1:3, :) = m.kb2ATP * 1e9 / 3600 / Navo / m.geometry.volume * m.stepSizeSec;
            polRates(3, 2) = 0;
            polRates(4, 1) = polRates(4, 1) * m.stateCooperativity*5;
            polRates(4, 2) = polRates(4, 2) * (m.stateCooperativity*5 + m.siteCooperativity);
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertElementsAlmostEqual(polRates, m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
        end
        
        function testcalculateDnaAR1234ATPPolymerizationRates_R4Bound(this)
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
            c.polymerizedRegions(1, :) = size(c.sequence, 1);
            c.complexBoundSites(:, :) = 0;
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1:4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(3));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(1:4)), 3) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(2));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345(4)), 1) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(4));
            c.complexBoundSites(m.dnaABoxStartPositions(m.dnaABoxIndexs_R12345([1 2 4])), 3) = m.enzymeComplexGlobalIndexs(m.enzymeIndexs_DnaA_Nmer_ATP(3));
            
            Navo = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            polRates = zeros(4, 2);
            polRates(4,   :) = m.kb1ATP * 1e9 / 3600 / Navo / m.geometry.volume * m.stepSizeSec;
            polRates(1:3, :) = m.kb2ATP * 1e9 / 3600 / Navo / m.geometry.volume * m.stepSizeSec;
            polRates(4, 1) = 0;
            polRates([1 2 4], 2) = 0;
            polRates(1, 1) = polRates(1, 1) * m.siteCooperativity;
            polRates(3, 2) = polRates(3, 2) * m.siteCooperativity;
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            assertElementsAlmostEqual(polRates, m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs));
        end
        
        function testCalculateDnaAAxPBindingRates_nothingBound(this)
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
            c.complexBoundSites(:, :) = 0;
            
            %ATP
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            [positionsStrands, rates] = m.calculateDnaAATPBindingRates(polATPs, polADPs);
            
            Navo = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            assertEqual(m.dnaABoxStartPositions, positionsStrands(:,1));
            assertAllEqual(1, positionsStrands(:,2));
            assertAllEqual(m.kb1ATP * 1e9 / 3600 / Navo / m.geometry.volume * m.stepSizeSec, rates(m.dnaABoxIndexs_9mer, 1));
            assertAllEqual(m.kb2ATP * 1e9 / 3600 / Navo / m.geometry.volume * m.stepSizeSec, rates(m.dnaABoxIndexs_8mer, 1));
            assertAllEqual(0, rates(m.dnaABoxIndexs_7mer, 1));
            
            %ADP
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            [positionsStrands, rates] = m.calculateDnaAADPBindingRates(polATPs, polADPs);
            
            Navo = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            assertEqual(m.dnaABoxStartPositions, positionsStrands(:,1));
            assertAllEqual(1, positionsStrands(:,2));
            assertAllEqual(m.kb1ADP * 1e9 / 3600 / Navo / m.geometry.volume * m.stepSizeSec, rates(m.dnaABoxIndexs_9mer, 1));
            assertAllEqual(m.kb2ADP * 1e9 / 3600 / Navo / m.geometry.volume * m.stepSizeSec, rates(m.dnaABoxIndexs_8mer, 1));
            assertAllEqual(0, rates(m.dnaABoxIndexs_7mer, 1));
        end
        
        function testCalculateDnaAAxPBindingRates_chromosome2(this)
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
            c.polymerizedRegions(1, :) = size(c.sequence, 1);
            c.complexBoundSites(:, :) = 0;
            
            %% ATP
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            [positionsStrands, rates] = m.calculateDnaAATPBindingRates(polATPs, polADPs);
            
            %binding just as likely on chromosome 2
            assertEqual(positionsStrands(1:end/2,1), positionsStrands(end/2+1:end,1));
            assertEqual(rates(1:end/2), rates(end/2+1:end));
            
            %% ADP
            [polATPs, polADPs] = m.calculateDnaAR1234Polymerization();
            [positionsStrands, rates] = m.calculateDnaAADPBindingRates(polATPs, polADPs);
            
            %binding just as likely on chromosome 2
            assertEqual(positionsStrands(1:end/2,1), positionsStrands(end/2+1:end,1));
            assertEqual(rates(1:end/2), rates(end/2+1:end));
        end
        
        function testCalculateDnaAAxPBPolymerizationRates_chromosome2(this)
            m = this.process;
            c = m.chromosome;
            
            c.initialize();
            c.polymerizedRegions(1, :) = size(c.sequence, 1);
            c.complexBoundSites(:, :) = 0;
            
            %% ATP - 1mers
            polATPs = ones(4, 2);
            polADPs = zeros(4, 2);
            
            %binding just as likely on chromosome 2
            [positionsStrands, rates] = m.calculateDnaAATPBindingRates(polATPs, polADPs);
            positionsStrands = cat(3, positionsStrands(1:end/2, :), positionsStrands(end/2+1:end, :));
            rates = cat(3, rates(1:end/2, :), rates(end/2+1:end, :));
            
            assertEqual(positionsStrands(:, 1, 1), positionsStrands(:, 1, 2));
            assertEqual(rates(m.dnaABoxIndexs_R12345, 1), rates(m.dnaABoxIndexs_R12345, 2));
            assertEqual(rates(:, 1), rates(:, 2));
            
            %polymerization just as likely on chromosome 2
            rates = m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs);
            assertEqual(rates(:, 1), rates(:, 2));
            
            %% ATP - 2mers
            polATPs = 2*ones(4, 2);
            polADPs = zeros(4, 2);
            
            %binding just as likely on chromosome 2
            [positionsStrands, rates] = m.calculateDnaAATPBindingRates(polATPs, polADPs);
            positionsStrands = cat(3, positionsStrands(1:end/2, :), positionsStrands(end/2+1:end, :));
            rates = cat(3, rates(1:end/2, :), rates(end/2+1:end, :));
            
            assertEqual(positionsStrands(:, 1, 1), positionsStrands(:, 1, 2));
            assertEqual(rates(m.dnaABoxIndexs_R12345, 1), rates(m.dnaABoxIndexs_R12345, 2));
            assertEqual(rates(:, 1), rates(:, 2));
            
            %polymerization just as likely on chromosome 2
            rates = m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs);
            assertEqual(rates(:, 1), rates(:, 2));
            
            %% ATP - 2/3mers
            polATPs = [2 2; 3 3; 2 2; 3 3];
            polADPs = zeros(4, 2);
            
            %binding just as likely on chromosome 2
            [positionsStrands, rates] = m.calculateDnaAATPBindingRates(polATPs, polADPs);
            positionsStrands = cat(3, positionsStrands(1:end/2, :), positionsStrands(end/2+1:end, :));
            rates = cat(3, rates(1:end/2, :), rates(end/2+1:end, :));
            
            assertEqual(positionsStrands(:, 1, 1), positionsStrands(:, 1, 2));
            assertEqual(rates(m.dnaABoxIndexs_R12345, 1), rates(m.dnaABoxIndexs_R12345, 2));
            assertEqual(rates(:, 1), rates(:, 2));
            
            %polymerization just as likely on chromosome 2
            rates = m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs);
            assertEqual(rates(:, 1), rates(:, 2));
            
            %% ATP - 6mers
            polATPs = 6*ones(4, 2);
            polADPs = zeros(4, 2);
            
            %binding just as likely on chromosome 2
            [positionsStrands, rates] = m.calculateDnaAATPBindingRates(polATPs, polADPs);
            positionsStrands = cat(3, positionsStrands(1:end/2, :), positionsStrands(end/2+1:end, :));
            rates = cat(3, rates(1:end/2, :), rates(end/2+1:end, :));
            
            assertEqual(positionsStrands(:, 1, 1), positionsStrands(:, 1, 2));
            assertEqual(rates(m.dnaABoxIndexs_R12345, 1), rates(m.dnaABoxIndexs_R12345, 2));
            assertEqual(rates(:, 1), rates(:, 2));
            
            %polymerization just as likely on chromosome 2
            rates = m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs);
            assertEqual(rates(:, 1), rates(:, 2));
            
            %% ATP - 7mers
            polATPs = 7*ones(4, 2);
            polADPs = zeros(4, 2);
            
            %binding just as likely on chromosome 2
            [positionsStrands, rates] = m.calculateDnaAATPBindingRates(polATPs, polADPs);
            positionsStrands = cat(3, positionsStrands(1:end/2, :), positionsStrands(end/2+1:end, :));
            rates = cat(3, rates(1:end/2, :), rates(end/2+1:end, :));
            
            assertEqual(positionsStrands(:, 1, 1), positionsStrands(:, 1, 2));
            assertEqual(rates(m.dnaABoxIndexs_R12345, 1), rates(m.dnaABoxIndexs_R12345, 2));
            assertEqual(rates(:, 1), rates(:, 2));
            
            %polymerization just as likely on chromosome 2
            rates = m.calculateDnaAR1234ATPPolymerizationRates(polATPs, polADPs);
            assertEqual(rates(:, 1), rates(:, 2));
            
            %% ADP
            polATPs = 2*ones(4, 2);
            polADPs = zeros(4, 2);
            
            %binding just as likely on chromosome 2
            [positionsStrands, rates] = m.calculateDnaAADPBindingRates(polATPs, polADPs);
            positionsStrands = cat(3, positionsStrands(1:end/2, :), positionsStrands(end/2+1:end, :));
            rates = cat(3, rates(1:end/2, :), rates(end/2+1:end, :));
            
            assertEqual(positionsStrands(:, 1, 1), positionsStrands(:, 1, 2));
            assertEqual(rates(m.dnaABoxIndexs_R12345, 1), rates(m.dnaABoxIndexs_R12345, 2));
            assertEqual(rates(:, 1), rates(:, 2));
            
            %polymerization just as likely on chromosome 2
            rates = m.calculateDnaAR1234ADPPolymerizationRates(polATPs, polADPs);
            assertEqual(rates(:, 1), rates(:, 2));
            
            %% ATP - 2/3mers
            polATPs = [2 2; 3 3; 2 2; 3 3];
            polADPs = zeros(4, 2);
            
            %binding just as likely on chromosome 2
            [positionsStrands, rates] = m.calculateDnaAADPBindingRates(polATPs, polADPs);
            positionsStrands = cat(3, positionsStrands(1:end/2, :), positionsStrands(end/2+1:end, :));
            rates = cat(3, rates(1:end/2, :), rates(end/2+1:end, :));
            
            assertEqual(positionsStrands(:, 1, 1), positionsStrands(:, 1, 2));
            assertEqual(rates(m.dnaABoxIndexs_R12345, 1), rates(m.dnaABoxIndexs_R12345, 2));
            assertEqual(rates(:, 1), rates(:, 2));
            
            %polymerization just as likely on chromosome 2
            rates = m.calculateDnaAR1234ADPPolymerizationRates(polATPs, polADPs);
            assertEqual(rates(:, 1), rates(:, 2));
            
            %% ATP - 7mers
            polATPs = 7*ones(4, 2);
            polADPs = zeros(4, 2);
            
            %binding just as likely on chromosome 2
            [positionsStrands, rates] = m.calculateDnaAADPBindingRates(polATPs, polADPs);
            positionsStrands = cat(3, positionsStrands(1:end/2, :), positionsStrands(end/2+1:end, :));
            rates = cat(3, rates(1:end/2, :), rates(end/2+1:end, :));
            
            assertEqual(positionsStrands(:, 1, 1), positionsStrands(:, 1, 2));
            assertEqual(rates(m.dnaABoxIndexs_R12345, 1), rates(m.dnaABoxIndexs_R12345, 2));
            assertEqual(rates(:, 1), rates(:, 2));
            
            %polymerization just as likely on chromosome 2
            rates = m.calculateDnaAR1234ADPPolymerizationRates(polATPs, polADPs);
            assertEqual(rates(:, 1), rates(:, 2));
        end
    end
    
    %methods for fitting replication initiation duration
    methods (Static)
        function [meanDnaACopyNumber, siteCooperativity, stateCooperativity, nDnaABoxes] = ...
                fitDuration(meanDnaACopyNumber, siteCooperativitys, stateCooperativity, nDnaABoxes)
            import edu.stanford.covert.cell.sim.process.ReplicationInitiation_Test;
            import edu.stanford.covert.util.ComputationUtil;
            
            nTrials = 32;
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {'ReplicationInitiation'});
            expectedDuration = sim.state('Time').replicationInitiationDuration;
            cellCycleLength = sim.state('Time').cellCycleLength;
            maxTime = 3 * expectedDuration;
            stepSizeSec = 100;
            useCommonInitialState = false;
            verbosity = 0;
            if nargin < 1
                meanDnaACopyNumber = 54;
            end
            if nargin < 2
                siteCooperativitys = sim.process('ReplicationInitiation').siteCooperativity * [0.9 1.1];
            end
            if nargin < 3
                stateCooperativity = sim.process('ReplicationInitiation').stateCooperativity;
            end
            if nargin < 4
                nDnaABoxes = numel(sim.process('ReplicationInitiation').dnaABoxStartPositions);
            end
            paralellize = true;
            %growthFunc = @(iter) log(2) / cellCycleLength * exp(log(2) * iter / cellCycleLength);
            growthFunc = [];
            clear sim;
            
            %clear saved files
            %delete('src_test/+edu/+stanford/+covert/+cell/+sim/+process/ReplicationInitiationDuration*.mat')            
            
            %fit
            if paralellize
                matlabpool('open');
            end
            try
                warningState = warning('query', 'WholeCell:warning');
                warning('off', 'WholeCell:warning');
                [siteCooperativity, ~, exitFlag, output] = ComputationUtil.fzero(@(siteCooperativity) ...
                    ReplicationInitiation_Test.fitDuration_diff(...
                    expectedDuration, meanDnaACopyNumber, siteCooperativity, stateCooperativity, nDnaABoxes, ...
                    nTrials, maxTime, stepSizeSec, growthFunc, useCommonInitialState, verbosity), ...
                    siteCooperativitys, optimset('TolX', 1e-2));
            catch exception
                if paralellize
                    matlabpool('close');
                end
                warning(warningState.state, 'WholeCell:warning');
                rethrow(exception);
            end
            if paralellize
                matlabpool('close');
            end
            warning(warningState.state, 'WholeCell:warning');
            
            if exitFlag ~= 1
                throw(MException('ReplicationInitiation_Test:error', output.message));
            end
        end
        
        function diff = fitDuration_diff(expectedDuration, meanDnaACopyNumber, siteCooperativity, stateCooperativity, nDnaABoxes, ...
                nTrials, maxTime, stepSizeSec, growthFunc, useCommonInitialState, verbosity)
            fprintf('Mean DnaA: %d, Site Cooperativity: %f, State Cooperativity: %f, No. DnaA boxes %d\n', ...
                meanDnaACopyNumber, siteCooperativity, stateCooperativity, nDnaABoxes);
            
            [durations, dnaACopyNumbers] = edu.stanford.covert.cell.sim.process.ReplicationInitiation_Test.meanDuration(...
                meanDnaACopyNumber, siteCooperativity, stateCooperativity, nDnaABoxes, ...
                nTrials, maxTime, stepSizeSec, growthFunc, useCommonInitialState, verbosity);
            
            diff = expectedDuration - mean(min(maxTime, durations));
            
            fprintf('Mean DnaA: %d, Site Cooperativity: %f, State Cooperativity: %f, No. DnaA boxes %d, Mean duration: %f\n', ...
                meanDnaACopyNumber, siteCooperativity, stateCooperativity, nDnaABoxes, mean(durations)/3600);
            
            save(sprintf('output/runSmallTests/ReplicationInitiationDuration-%d-%f-%f-%d.mat', ...
                meanDnaACopyNumber, siteCooperativity, stateCooperativity, nDnaABoxes), ...
                'durations', 'dnaACopyNumbers', 'siteCooperativity', 'stateCooperativity', 'maxTime', 'stepSizeSec');
        end
        
        function [durations, dnaACopyNumbers] = meanDuration(meanDnaACopyNumber, siteCooperativity, stateCooperativity, nDnaABoxes, ...
                nTrials, iterMax, stepSizeSec, growthFunc, useCommonInitialState, verbosity)
            durations = zeros(nTrials, 1);
            dnaACopyNumbers = zeros(nTrials, 1);
            if useCommonInitialState
                initSeeds = ones(nTrials, 1);
            else
                initSeeds = 10 * (1:nTrials)';
            end
            seeds = (1:nTrials)';
            parfor k = 1:nTrials
                [durations(k), dnaACopyNumbers(k)] = edu.stanford.covert.cell.sim.process.ReplicationInitiation_Test.sampleDuration(...
                    meanDnaACopyNumber, siteCooperativity, stateCooperativity, nDnaABoxes, ...
                    iterMax, stepSizeSec, initSeeds(k), seeds(k), growthFunc, verbosity);
            end
        end
        
        function [replicationInitiationTime, dnaACopyNumber] = sampleDuration(meanDnaACopyNumber, siteCooperativity, stateCooperativity, nDnaABoxes, ...
                iterMax, stepSizeSec, initSeed, seed, growthFunc, verbosity)
            import edu.stanford.covert.util.ConstantUtil;
            
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'MacromolecularComplexation'
                'Metabolism'
                'ReplicationInitiation'
                });
            sim.applyOptions('verbosity', verbosity);
            
            g = sim.gene;
            comp = sim.compartment;
            time = sim.state('Time');
            mass = sim.state('Mass');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            mr = sim.state('MetabolicReaction');
            met = sim.process('Metabolism');
            repInit = sim.process('ReplicationInitiation');
            
            sim.state('MetabolicReaction').initialGrowthFilterWidth = Inf;
            mass.initialFractionNTPsInRNAs = 1;
            mass.initialFractionAAsInMonomers = 1;
            
            pc.formationProcesses(~ismember(pc.formationProcesses, [9 10 19])) = 0;
            pc.formationProcesses(pc.formationProcesses == 9) = 1;  %Macromolecular complexatopm
            pc.formationProcesses(pc.formationProcesses == 19) = 3; %Replication initiation
            pc.formationProcesses(pc.formationProcesses == 10) = 2; %Metabolism
            
            %% adjust constants
            repInit.siteCooperativity = siteCooperativity;
            repInit.stateCooperativity = stateCooperativity;
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            monExp = mRNAExp ./ (log(2) / time.cellCycleLength + pm.decayRates(pm.matureIndexs));
            monExp = monExp / sum(monExp);
            monCnts = mass.cellInitialDryWeight * mass.dryWeightFractionProtein * ConstantUtil.nAvogadro / ...
                (monExp' * pm.molecularWeights(pm.matureIndexs)) * monExp;
            
            tfs = rna.matureRNAGeneComposition(g.mRNAIndexs(repInit.enzymeMonomerGlobalIndexs), :) ~= 0;
            rna.expression(rna.matureIndexs(tfs)) = rna.expression(rna.matureIndexs(tfs)) * ...
                meanDnaACopyNumber / monCnts(repInit.enzymeMonomerGlobalIndexs);
            rna.expression(rna.matureIndexs) = rna.expression(rna.matureIndexs) / ...
                sum(rna.expression(rna.matureIndexs));
            
            if nDnaABoxes ~= round(nDnaABoxes)
                nDnaABoxes = round(nDnaABoxes);
                warning('WholeCell:warning', 'Non-integer number of DnaA boxes');
            end
            repInit.sampleDnaABoxes(nDnaABoxes);
            
            %% constants needed to mock transcription, translation
            rnaProd = edu.stanford.covert.util.ComputationUtil.invertCompositionMatrix(rna.nascentRNAMatureRNAComposition) * ...
                (rna.expression(rna.matureIndexs) .* ...
                (log(2) / time.cellCycleLength + rna.decayRates(rna.matureIndexs)));
            rnaProd = rna.nascentRNAMatureRNAComposition * rnaProd / sum(rnaProd) * ...
                (rna.expression(rna.matureIndexs) * mass.cellInitialDryWeight * mass.dryWeightFractionRNA * ConstantUtil.nAvogadro / ...
                (rna.molecularWeights(rna.matureIndexs)' * rna.expression(rna.matureIndexs)))' * (...
                + log(2) / time.cellCycleLength ...
                + rna.decayRates(rna.matureIndexs));
            rnaProd = 0.3929 * rnaProd / sum(rnaProd);
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            mRNAExp = mRNAExp / sum(mRNAExp);
            nMonProd = mass.cellInitialDryWeight * mass.dryWeightFractionProtein * ConstantUtil.nAvogadro / (pm.molecularWeights(pm.matureIndexs)' * mRNAExp) * mRNAExp' * (...
                + log(2) / time.cellCycleLength ...
                + pm.decayRates(pm.matureIndexs));
                       
            matureRNADecayRates = min(1, rna.decayRates(rna.matureIndexs));
            matureRNAMRNAComposition = rna.matureRNAGeneComposition(g.mRNAIndexs, :);
            
            matureMonDecayRates = pm.decayRates(pm.matureIndexs);
            matureCpxDecayRates = pc.decayRates(pc.matureIndexs);
            
            nComplexs = numel(pc.matureIndexs);
            pcMonComp = sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3);
            notFormComplexTfs = pc.formationProcesses(pc.matureIndexs) ~= 1 | ...
                any(any(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3), 1)';
            
            metGasIdxs = met.substrateIndexs({'O2'; 'CO2'});
            
            metSubstrateMonomerLocalIndexs = met.substrateMonomerLocalIndexs(:, 1);
            metSubstrateComplexLocalIndexs = met.substrateComplexLocalIndexs(:, 1);
            metSubstrateMonomerGlobalIndexs = met.substrateMonomerGlobalIndexs(:, 1);
            metSubstrateComplexGlobalIndexs = met.substrateComplexGlobalIndexs(:, 1);
            
            %% initialize
            if ~isempty(growthFunc)
                edu.stanford.covert.cell.sim.process.ReplicationInitiation_Test.seedRandStream(sim, initSeed);
                sim.initializeState();
                pc.counts(sub2ind(size(pc.counts), pc.matureIndexs, pc.compartments(pc.matureIndexs))) = ...
                    + pc.counts(sub2ind(size(pc.counts), pc.matureIndexs, pc.compartments(pc.matureIndexs))) ...
                    + pc.counts(sub2ind(size(pc.counts), pc.nascentIndexs, pc.compartments(pc.nascentIndexs)));
                pc.counts(pc.nascentIndexs, :) = 0;
            else
                i = 0;
                while true
                    if initSeed
                        edu.stanford.covert.cell.sim.process.ReplicationInitiation_Test.seedRandStream(sim, initSeed + i);
                    end
                    
                    sim.initializeState();
                    pc.counts(sub2ind(size(pc.counts), pc.matureIndexs, pc.compartments(pc.matureIndexs))) = ...
                        + pc.counts(sub2ind(size(pc.counts), pc.matureIndexs, pc.compartments(pc.matureIndexs))) ...
                        + pc.counts(sub2ind(size(pc.counts), pc.nascentIndexs, pc.compartments(pc.nascentIndexs)));
                    pc.counts(pc.nascentIndexs, :) = 0;
                    sim.evolveState();
                    
                    if mr.growth > 0.2 * log(2) / time.cellCycleLength
                        break;
                    end
                    
                    i = i + 1;
                end
            end
            repInit.substrates(:) = 0;
            repInit.substrates(repInit.substrateIndexs_atp) = 1e6;
            
            %% seed
            edu.stanford.covert.cell.sim.process.ReplicationInitiation_Test.seedRandStream(sim, seed);
            
            %% keep track of initial state
            matureRNACounts = rna.counts(rna.matureIndexs, comp.cytosolIndexs);
            matureMonCounts = sum(pm.counts(pm.matureIndexs, :), 2);
            boundMonCounts = sum(pm.counts(pm.boundIndexs, :), 2);
            matureCpxCounts = sum(pc.counts(pc.matureIndexs, :), 2);
            boundCpxCounts = sum(pc.counts(pc.boundIndexs, :), 2);
            tfs = any(any(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3), 1);
            matureRNACounts(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                matureRNACounts(setdiff(1:end, rna.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), tfs, :), 3) * (matureCpxCounts(tfs) + boundCpxCounts(tfs));
            matureMonCounts = ...
                matureMonCounts + ...
                sum(pc.proteinComplexComposition(g.mRNAIndexs, tfs, :), 3) * (matureCpxCounts(tfs) + boundCpxCounts(tfs));
            matureCpxCounts(tfs) = 0;
            boundCpxCounts(tfs) = 0;
            
            initGas = met.substrates(metGasIdxs, met.compartmentIndexs_extracellular);
            
            %% simulate
            replicationInitiationTime = NaN;
            
            assertElementsAlmostEqual(...
                mass.cellInitialDryWeight * (...
                + mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA ...
                + mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein ...
                ) * edu.stanford.covert.util.ConstantUtil.nAvogadro, ...
                + matureRNACounts' * rna.molecularWeights(rna.matureIndexs) ...
                + (matureMonCounts+boundMonCounts)' * pm.molecularWeights(pm.matureIndexs) ...
                + (matureCpxCounts+boundCpxCounts)' * pc.molecularWeights(pc.matureIndexs), ...
                'relative', 5e-2)           
            
            if sim.verbosity > 0
                fprintf('\t%4s %5s %9s %7s %13s\n', 'Seed', 'Iter ', ' Growth  ', ' DnaA  ', '    R1-5     ');
                fprintf('\t%4s %5s %9s %7s %13s\n', '====', '=====', '=========', '=======', '=============');
            end
            for iter = 1:iterMax
                if sim.verbosity > 0 && mod(iter, 1000) == 1
                    [polATPs, polADPs] = repInit.calculateDnaAR1234Polymerization();
                    s1 = repInit.calculateDnaAR1234ComplexSize(polATPs, polADPs);
                    b1 = repInit.calculateDnaABoxStatus();
                    o1 = repInit.calcuateIsDnaAR5Occupied();
                    s1 = s1(1);
                    b1 = b1(repInit.dnaABoxIndexs_R12345, 1);
                    o1 = o1(1);
                    
                    fprintf('\t%4d %5d %1.2e %3d %3d %d %d %d %d %d %d %d\n',...
                        seed, iter, mr.growth, ...
                        repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ATP), ...
                        repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_1mer_ATP), ...
                        s1, b1(1), b1(2), b1(3), b1(4), b1(5), o1);
                end
                
                %% time
                time.values = iter;
                
                if mod(iter, stepSizeSec) == 1
                    %% mock transcription, RNA maturation, and decay
                    newRNAs = rna.randStream.stochasticRound(rnaProd * mr.growth * time.cellCycleLength / log(2) * stepSizeSec);
                    decayedRNAs = rna.randStream.stochasticRound(matureRNACounts .* matureRNADecayRates * stepSizeSec);
                    matureRNACounts = ...
                        + matureRNACounts ...
                        + newRNAs ...
                        - decayedRNAs;
                    
                    %% mock translation, protein maturation, and decay
                    mRNAExp = matureRNAMRNAComposition * matureRNACounts;
                    newMonomers = pm.randStream.stochasticRound(nMonProd * mr.growth * time.cellCycleLength / log(2) * max(0, mRNAExp / sum(mRNAExp)) * stepSizeSec);
                    decayedMatureMonomers = pm.randStream.stochasticRound(matureMonCounts .* matureMonDecayRates * stepSizeSec);
                    decayedMatureComplexs = pc.randStream.stochasticRound(matureCpxCounts .* matureCpxDecayRates * stepSizeSec);
                    decayedBoundMonomers = pm.randStream.stochasticRound(boundMonCounts .* matureMonDecayRates * stepSizeSec);
                    decayedBoundComplexs = pc.randStream.stochasticRound(boundCpxCounts .* matureCpxDecayRates * stepSizeSec);
                    
                    tmp = zeros(size(pc.counts));
                    tmp(pc.boundIndexs, comp.cytosolIndexs) = -decayedBoundComplexs;
                    notUpdatingComplexs = pc.updateExternalState(tmp, true);
                    decayedBoundComplexs = decayedBoundComplexs - notUpdatingComplexs(pc.boundIndexs, comp.cytosolIndexs);
                    
                    matureMonCounts = ...
                        + matureMonCounts ...
                        + newMonomers ...
                        - decayedMatureMonomers;
                    matureCpxCounts = ...
                        + matureCpxCounts ...
                        - decayedMatureComplexs;
                    boundMonCounts = ...
                        + boundMonCounts ...
                        - decayedBoundMonomers;
                    boundCpxCounts = ...
                        + boundCpxCounts ...
                        - decayedBoundComplexs;
                    
                    %% mock protein complexation
                    newComplexs = max(0, floor(min(matureMonCounts(:, ones(nComplexs, 1)) ./ pcMonComp, [], 1))');
                    newComplexs(notFormComplexTfs) = 0;
                    matureMonCounts = ...
                        + matureMonCounts ...
                        - pcMonComp * newComplexs;
                    matureCpxCounts = ...
                        + matureCpxCounts ...
                        + newComplexs;
                    
                    %% mock growth
                    if ~isempty(growthFunc)
                        mr.growth = growthFunc(iter);
                    else
                        met.substrates(metGasIdxs, met.compartmentIndexs_extracellular) = initGas;
                        met.substrates(metSubstrateMonomerLocalIndexs, :) = matureMonCounts(metSubstrateMonomerGlobalIndexs, ones(3, 1));
                        met.substrates(metSubstrateComplexLocalIndexs, :) = matureCpxCounts(metSubstrateComplexGlobalIndexs, ones(3, 1));
                        
                        met.enzymes(met.enzymeMonomerLocalIndexs) = matureMonCounts(met.enzymeMonomerGlobalIndexs);
                        met.enzymes(met.enzymeComplexLocalIndexs) = matureCpxCounts(met.enzymeComplexGlobalIndexs);
                        
                        met.evolveState();
                    end
                end
                
                %% simulate DnaA-ATP dynamics
                repInit.enzymes(repInit.enzymeMonomerLocalIndexs) = matureMonCounts(repInit.enzymeMonomerGlobalIndexs);
                repInit.enzymes(repInit.enzymeComplexLocalIndexs) = matureCpxCounts(repInit.enzymeComplexGlobalIndexs);
                repInit.boundEnzymes(repInit.enzymeMonomerLocalIndexs) = boundMonCounts(repInit.enzymeMonomerGlobalIndexs);
                repInit.boundEnzymes(repInit.enzymeComplexLocalIndexs) = boundCpxCounts(repInit.enzymeComplexGlobalIndexs);
                
                repInit.evolveState();
                
                matureMonCounts(repInit.enzymeMonomerGlobalIndexs) = repInit.enzymes(repInit.enzymeMonomerLocalIndexs);
                matureCpxCounts(repInit.enzymeComplexGlobalIndexs) = repInit.enzymes(repInit.enzymeComplexLocalIndexs);
                boundMonCounts(repInit.enzymeMonomerGlobalIndexs) = repInit.boundEnzymes(repInit.enzymeMonomerLocalIndexs);
                boundCpxCounts(repInit.enzymeComplexGlobalIndexs) = repInit.boundEnzymes(repInit.enzymeComplexLocalIndexs);
                
                if any(repInit.calculateIsDnaAORIComplexAssembled)
                    replicationInitiationTime = iter;
                    if sim.verbosity > 0
                        fprintf('\tReplication Initiated at %d s\n', iter);
                    end
                    break;
                end
            end
            
            dnaACopyNumber = ...
                + repInit.enzymes(repInit.enzymeIndexs_DnaA) ...
                + (1:7) * repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ATP) ...
                + (1:7) * repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP) ....
                + repInit.boundEnzymes(repInit.enzymeIndexs_DnaA) ...
                + (1:7) * repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_Nmer_ATP) ...
                + (1:7) * repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP);
            
            if sim.verbosity > 0 && snan(replicationInitiationTime)
                fprintf('\tReplication did not initiate within %d s\n', iterMax);
            end
        end
        
        function seedRandStream(sim, seed)
            if sim.verbosity > 0
                fprintf('\tSeed %d\n', seed);
            end
            
            sim.applyOptions('seed', seed);
            sim.seedRandStream();
            for i = 1:numel(sim.states)
                o = sim.states{i};
                o.seed = seed;
                o.seedRandStream();
            end
            
            for i = 1:numel(sim.processes)
                o = sim.processes{i};
                o.seed = seed;
                o.seedRandStream();
            end
        end
    end
end
