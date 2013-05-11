%RNA modification process test case
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef RNAModification_Test < edu.stanford.covert.cell.sim.ReactionProcessTestCase
    %constructor
    methods
        function this = RNAModification_Test(name)
            this = this@edu.stanford.covert.cell.sim.ReactionProcessTestCase(name);
        end
    end

    %RNA indices
    properties
        TU_001
        MG471
        MG472
        MG475
    end

    %fixtures
    methods
        function loadTestFixture(this)
            this.loadTestFixture@edu.stanford.covert.cell.sim.ReactionProcessTestCase();

            m = this.process;
            this.TU_001 = find(strcmp(m.unmodifiedRNAWholeCellModelIDs, 'TU_001'), 1, 'first');
            this.MG471  = find(strcmp(m.unmodifiedRNAWholeCellModelIDs, 'MG471'),  1, 'first');
            this.MG472  = find(strcmp(m.unmodifiedRNAWholeCellModelIDs, 'MG472'),  1, 'first');
            this.MG475  = find(strcmp(m.unmodifiedRNAWholeCellModelIDs, 'MG475'),  1, 'first');
        end

        function loadSimpleTestFixture(this)
            %process
            m = this.process;

            %whole cell model IDs
            m.reactionWholeCellModelIDs = {...
                'MG471_MODIFICATION_001'; 'MG471_MODIFICATION_002'; 'MG471_MODIFICATION_003'; 'MG472_MODIFICATION_001'; 'MG472_MODIFICATION_002'; ...
                'MG475_MODIFICATION_001'};
            m.substrateWholeCellModelIDs = {...
                'AHCYS'; 'AMET'; 'AMP'; 'ATP'; 'cmnm5s2UMP'; 'CMP'; 'CYS'; 'FAD'; 'FTHF5'; 'GLY'; 'GmMP'; 'GMP'; 'H'; ...
                'H2O'; 'k2CMP'; 'LYS'; 'm1GMP'; 'm2GMP'; 'm62AMP'; 'm6AMP'; 'm7GMP'; 'PPI'; 'PSIURIMP'; 's2UMP'; 's4UMP'; ...
                'SER'; 'THF'; 'UmMP'; 'UMP'};
            m.enzymeWholeCellModelIDs = {...
                'MG_209_MONOMER'; 'MG_295_MONOMER'; 'MG_370_MONOMER'; 'MG_380_MONOMER'; 'MG_463_MONOMER'; 'MG_008_379_TETRAMER'; ...
                'MG_084_TETRAMER'; 'MG_182_DIMER'; 'MG_252_DIMER'; 'MG_346_DIMER'; 'MG_347_DIMER'; 'MG_372_DIMER'; 'MG_445_DIMER'};
            m.unmodifiedRNAWholeCellModelIDs = {...
                'TU_001'; 'MG471'; 'MG472'; 'MG475'};
            m.modifiedRNAWholeCellModelIDs = m.unmodifiedRNAWholeCellModelIDs;

            %names
            m.reactionNames  = m.reactionWholeCellModelIDs;
            m.substrateNames = m.substrateWholeCellModelIDs;
            m.enzymeNames    = m.enzymeWholeCellModelIDs;

            %indices
            m.substrateIndexs_amet        = m.substrateIndexs({'AMET'});
            m.substrateIndexs_fthf5       = m.substrateIndexs({'FTHF5'});
            m.substrateIndexs_cys         = m.substrateIndexs({'CYS'});
            m.substrateIndexs_gly         = m.substrateIndexs({'GLY'});
            m.substrateIndexs_lys         = m.substrateIndexs({'LYS'});
            m.substrateIndexs_atp         = m.substrateIndexs({'ATP'});
            m.substrateIndexs_water       = m.substrateIndexs({'H2O'});
            m.substrateIndexs_nmp         = m.substrateIndexs({'AMP';'CMP';'GMP';'UMP'});
            m.substrateIndexs_modifiedNMP = m.substrateIndexs({'PSIURIMP';'m62AMP';'m2GMP';'m7GMP';'GmMP';'UmMP';'m1GMP';'k2CMP';'s4UMP';'cmnm5s2UMP'});

            m.enzymeIndexs_rRNA16SDimethyladenosineTransferase = m.enzymeIndexs({'MG_463_MONOMER'});      %dimethyladenosine transferase
            m.enzymeIndexs_rRNA16SMethyltransferaseGidB        = m.enzymeIndexs({'MG_380_MONOMER'});      %methyltransferase GidB
            m.enzymeIndexs_rRNA23SMethyltransferaseI           = m.enzymeIndexs({'MG_252_DIMER'});        %23S rRNA methyltransferase; G2251
            m.enzymeIndexs_rRNA23SMethyltransferaseII          = m.enzymeIndexs({'MG_346_DIMER'});        %23S rRNA methyltransferase; U2552
            m.enzymeIndexs_rRNA23SPseudouridineSynthaseI       = m.enzymeIndexs({'MG_209_MONOMER'});      %23S rRNA pseudouridine synthase; U955, U2504, U2580
            m.enzymeIndexs_rRNA23SPseudouridineSynthaseII      = m.enzymeIndexs({'MG_370_MONOMER'});      %23S rRNA pseudouridine synthase; U1911, U1915, U1917
            m.enzymeIndexs_tRNAGuanineN1Methyltransferase      = m.enzymeIndexs({'MG_445_DIMER'});        %tRNA (guanine-N1)-methyltransferase
            m.enzymeIndexs_tRNAGuanineN7Methyltransferase      = m.enzymeIndexs({'MG_347_DIMER'});        %tRNA (guanine-N(7)-)-methyltransferase
            m.enzymeIndexs_tRNALysidineSynthetase              = m.enzymeIndexs({'MG_084_TETRAMER'});     %tRNA(Ile)-lysidine synthetase
            m.enzymeIndexs_tRNAPseudouridineSynthase           = m.enzymeIndexs({'MG_182_DIMER'});        %tRNA pseudouridine synthase A
            m.enzymeIndexs_tRNAUracil2Sulfurtransferase        = m.enzymeIndexs({'MG_295_MONOMER'});      %tRNA U34 sulfurtransferase
            m.enzymeIndexs_tRNAUracil4Sulfurtransferase        = m.enzymeIndexs({'MG_372_DIMER'});        %thiamine biosynthesis/tRNA modification protein ThiI
            m.enzymeIndexs_tRNAUracil5Carboxymethylaminomethyl = m.enzymeIndexs({'MG_008_379_TETRAMER'}); %tRNA uridine 5-carboxymethylaminomethyl modification enzyme

            %other physical properties
            m.matureRNALengths = repmat(500, size(m.unmodifiedRNAWholeCellModelIDs));

            %RNA modified by each reaction
            m.reactionModificationMatrix = zeros(length(m.reactionWholeCellModelIDs), length(m.unmodifiedRNAWholeCellModelIDs));
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_001'), strcmp(m.unmodifiedRNAWholeCellModelIDs,'MG471'))=1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_002'), strcmp(m.unmodifiedRNAWholeCellModelIDs,'MG471'))=1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_003'), strcmp(m.unmodifiedRNAWholeCellModelIDs,'MG471'))=1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG472_MODIFICATION_001'), strcmp(m.unmodifiedRNAWholeCellModelIDs,'MG472'))=1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG472_MODIFICATION_002'), strcmp(m.unmodifiedRNAWholeCellModelIDs,'MG472'))=1;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG475_MODIFICATION_001'), strcmp(m.unmodifiedRNAWholeCellModelIDs,'MG475'))=1;

            %reaction catalysis
            m.reactionCatalysisMatrix = zeros(length(m.reactionWholeCellModelIDs), length(m.enzymeWholeCellModelIDs));
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_001'), strcmp(m.enzymeWholeCellModelIDs,'MG_347_DIMER'))=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_002'), strcmp(m.enzymeWholeCellModelIDs,'MG_182_DIMER'))=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_003'), strcmp(m.enzymeWholeCellModelIDs,'MG_372_DIMER'))=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG472_MODIFICATION_001'), strcmp(m.enzymeWholeCellModelIDs,'MG_347_DIMER'))=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG472_MODIFICATION_002'), strcmp(m.enzymeWholeCellModelIDs,'MG_372_DIMER'))=1;
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG475_MODIFICATION_001'), strcmp(m.enzymeWholeCellModelIDs,'MG_372_DIMER'))=1;

            %reaction stoichiometry
            m.reactionStoichiometryMatrix = zeros(length(m.substrateWholeCellModelIDs), length(m.reactionWholeCellModelIDs));

            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'AMET'), m.reactionCatalysisMatrix(:,strcmp(m.enzymeWholeCellModelIDs,'MG_347_DIMER'))>0)=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'AHCYS'), m.reactionCatalysisMatrix(:,strcmp(m.enzymeWholeCellModelIDs,'MG_347_DIMER'))>0)=1;

            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ATP'), m.reactionCatalysisMatrix(:,strcmp(m.enzymeWholeCellModelIDs,'MG_372_DIMER'))>0)=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'CYS'), m.reactionCatalysisMatrix(:,strcmp(m.enzymeWholeCellModelIDs,'MG_372_DIMER'))>0)=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'H2O'), m.reactionCatalysisMatrix(:,strcmp(m.enzymeWholeCellModelIDs,'MG_372_DIMER'))>0)=-1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'AMP'), m.reactionCatalysisMatrix(:,strcmp(m.enzymeWholeCellModelIDs,'MG_372_DIMER'))>0)=1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'H'),   m.reactionCatalysisMatrix(:,strcmp(m.enzymeWholeCellModelIDs,'MG_372_DIMER'))>0)=1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'PPI'), m.reactionCatalysisMatrix(:,strcmp(m.enzymeWholeCellModelIDs,'MG_372_DIMER'))>0)=1;
            m.reactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'SER'), m.reactionCatalysisMatrix(:,strcmp(m.enzymeWholeCellModelIDs,'MG_372_DIMER'))>0)=1;            
            
            %reaction kinetics
            m.enzymeBounds = repmat([-Inf Inf], length(m.reactionWholeCellModelIDs), 1);
            m.enzymeBounds(m.reactionCatalysisMatrix(:,strcmp(m.enzymeWholeCellModelIDs,'MG_182_DIMER'))>0,2)=0.1800;
            m.enzymeBounds(m.reactionCatalysisMatrix(:,strcmp(m.enzymeWholeCellModelIDs,'MG_347_DIMER'))>0,2)=3.2661;
            m.enzymeBounds(m.reactionCatalysisMatrix(:,strcmp(m.enzymeWholeCellModelIDs,'MG_372_DIMER'))>0,2)=0.0066;
            
            m.initializeSpeciesNetwork();

            %molecular weights
            m.substrateMolecularWeights = [...
                384.4100, 399.4444, 345.2049, 503.1489, 425.3077, 321.1802, 121.1579, 783.5327, 473.4384,  75.0664, 375.2308, 361.2043,   1.0079, ...
                18.0152, 449.3520, 147.1949, 375.2308, 375.2308, 373.2579, 359.2314, 376.2387, 174.9513, 322.1650, 338.2306, 338.2306, 105.0923, ...
                443.4125, 336.1915, 322.1650]';
            m.enzymeMolecularWeights = 1e5 * [...
                0.3501 0.4190 0.3798 0.2218 0.2997 2.3975 1.3860 0.5523 0.5544 0.3864 0.4961 0.8837 0.5296]';
            m.unmodifiedRNAMolecularWeights = 1e6 * (1:length(m.unmodifiedRNAWholeCellModelIDs))';
            m.modifiedRNAMolecularWeights    = m.unmodifiedRNAMolecularWeights - ...
                m.reactionModificationMatrix'*m.reactionStoichiometryMatrix'*m.substrateMolecularWeights;

            %initial state
            m.substrates     = zeros(length(m.substrateWholeCellModelIDs),    1);
            m.enzymes        = zeros(length(m.enzymeWholeCellModelIDs),       1);
            m.boundEnzymes   = zeros(length(m.enzymeWholeCellModelIDs),       1);
            m.unmodifiedRNAs = zeros(length(m.unmodifiedRNAWholeCellModelIDs),1);
            m.modifiedRNAs   = zeros(length(m.modifiedRNAWholeCellModelIDs),  1);
        end
    end

    %tests
    methods
        function testRNAWithNoModifications(this)
            m = this.process;
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.unmodifiedRNAs(:) = 0;
            m.unmodifiedRNAs(this.TU_001) = 1;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertEqual(0, m.unmodifiedRNAs(this.TU_001));
            assertEqual(1, m.modifiedRNAs(this.TU_001));
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(zeros(size(m.enzymes)), m.enzymes);
        end

        function testOneTRNAU8SulfurTransferase(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.substrates(m.substrateIndexs_cys) = 1;
            m.substrates(m.substrateIndexs_water) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_tRNAUracil4Sulfurtransferase) = 1e6;
            m.unmodifiedRNAs(:) = 0;
            m.unmodifiedRNAs(this.MG475) = 1;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertEqual(0,   m.unmodifiedRNAs(this.MG475));
            assertEqual(1,   m.modifiedRNAs(this.MG475));
            assertEqual(0,   m.substrates(m.substrateIndexs_atp));
            assertEqual(0,   m.substrates(m.substrateIndexs_cys));
            assertEqual(0,   m.substrates(m.substrateIndexs_water));
            assertEqual(1,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'AMP')));
            assertEqual(1,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'H')));
            assertEqual(1,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'PPI')));
            assertEqual(1,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'SER')));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_tRNAUracil4Sulfurtransferase));
        end

        function testOneTRNAU8SulfurTransferase_noEnyzme(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.substrates(m.substrateIndexs_cys) = 1;
            m.substrates(m.substrateIndexs_water) = 1;
            m.enzymes(:) = 0;
            m.unmodifiedRNAs(:) = 0;
            m.unmodifiedRNAs(this.MG475) = 1;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertEqual(1, m.unmodifiedRNAs(this.MG475));
            assertEqual(0, m.modifiedRNAs(this.MG475));
            assertEqual(1, m.substrates(m.substrateIndexs_atp));
            assertEqual(1, m.substrates(m.substrateIndexs_cys));
            assertEqual(1, m.substrates(m.substrateIndexs_water));
            assertEqual(0, m.substrates(strcmp(m.substrateWholeCellModelIDs,'AMP')));
            assertEqual(0, m.substrates(strcmp(m.substrateWholeCellModelIDs,'H')));
            assertEqual(0, m.substrates(strcmp(m.substrateWholeCellModelIDs,'PPI')));
            assertEqual(0, m.substrates(strcmp(m.substrateWholeCellModelIDs,'SER')));
            assertEqual(0, m.enzymes(m.enzymeIndexs_tRNAUracil4Sulfurtransferase));
        end

        function testOneTRNAU8SulfurTransferase_noATP(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 0;
            m.substrates(m.substrateIndexs_cys) = 1;
            m.substrates(m.substrateIndexs_water) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_tRNAUracil4Sulfurtransferase) = 1e6;
            m.unmodifiedRNAs(:) = 0;
            m.unmodifiedRNAs(this.MG475) = 1;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertEqual(1,   m.unmodifiedRNAs(this.MG475));
            assertEqual(0,   m.modifiedRNAs(this.MG475));
            assertEqual(0,   m.substrates(m.substrateIndexs_atp));
            assertEqual(1,   m.substrates(m.substrateIndexs_cys));
            assertEqual(1,   m.substrates(m.substrateIndexs_water));
            assertEqual(0,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'AMP')));
            assertEqual(0,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'H')));
            assertEqual(0,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'PPI')));
            assertEqual(0,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'SER')));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_tRNAUracil4Sulfurtransferase));
        end

        function testRNARequiringMultipleModifications(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 1;
            m.substrates(m.substrateIndexs_cys) = 1;
            m.substrates(m.substrateIndexs_water) = 1;
            m.substrates(m.substrateIndexs_amet) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_tRNAPseudouridineSynthase)      = 1e6;
            m.enzymes(m.enzymeIndexs_tRNAGuanineN7Methyltransferase) = 1e6;
            m.enzymes(m.enzymeIndexs_tRNAUracil4Sulfurtransferase)   = 1e6;
            m.unmodifiedRNAs(:) = 0;
            m.unmodifiedRNAs(this.MG471) = 1;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertEqual(0,   m.unmodifiedRNAs(this.MG471));
            assertEqual(1,   m.modifiedRNAs(this.MG471));
            assertEqual(0,   m.substrates(m.substrateIndexs_atp));
            assertEqual(0,   m.substrates(m.substrateIndexs_cys));
            assertEqual(0,   m.substrates(m.substrateIndexs_water));
            assertEqual(0,   m.substrates(m.substrateIndexs_amet));
            assertEqual(1,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'AMP')));
            assertEqual(1,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'H')));
            assertEqual(1,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'PPI')));
            assertEqual(1,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'SER')));
            assertEqual(1,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'AHCYS')));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_tRNAPseudouridineSynthase));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_tRNAGuanineN7Methyltransferase));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_tRNAUracil4Sulfurtransferase));
        end

        function testRNARequiringMultipleModifications_oneTooFewATP(this)
            m = this.process;
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = 0;
            m.substrates(m.substrateIndexs_cys) = 1;
            m.substrates(m.substrateIndexs_water) = 1;
            m.substrates(m.substrateIndexs_amet) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_tRNAPseudouridineSynthase)      = 1e6;
            m.enzymes(m.enzymeIndexs_tRNAGuanineN7Methyltransferase) = 1e6;
            m.enzymes(m.enzymeIndexs_tRNAUracil4Sulfurtransferase)   = 1e6;
            m.unmodifiedRNAs(:) = 0;
            m.unmodifiedRNAs(this.MG471) = 1;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertEqual(1,   m.unmodifiedRNAs(this.MG471));
            assertEqual(0,   m.modifiedRNAs(this.MG471));
            assertEqual(0,   m.substrates(m.substrateIndexs_atp));
            assertEqual(1,   m.substrates(m.substrateIndexs_cys));
            assertEqual(1,   m.substrates(m.substrateIndexs_water));
            assertEqual(1,   m.substrates(m.substrateIndexs_amet));
            assertEqual(0,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'AMP')));
            assertEqual(0,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'H')));
            assertEqual(0,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'PPI')));
            assertEqual(0,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'SER')));
            assertEqual(0,   m.substrates(strcmp(m.substrateWholeCellModelIDs,'AHCYS')));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_tRNAPseudouridineSynthase));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_tRNAGuanineN7Methyltransferase));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_tRNAUracil4Sulfurtransferase));
        end

        function testTRNAGuanineN7MethylTransferase(this)
            m = this.process;

            m.reactionModificationMatrix(m.reactionModificationMatrix(:,this.MG471)==1,this.MG472)=1;
            m.reactionModificationMatrix(m.reactionModificationMatrix(:,this.MG471)==1,this.MG471)=0;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_001'),:)=0;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_001'),this.MG471)=1;
            
            m.initializeSpeciesNetwork();

            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_amet) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_tRNAGuanineN7Methyltransferase) = 1e6;
            m.unmodifiedRNAs(:) = 0;
            m.unmodifiedRNAs(this.MG471) = 1;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertEqual(0, m.unmodifiedRNAs(this.MG471));
            assertEqual(1, m.modifiedRNAs(this.MG471));
            assertEqual(0, m.substrates(m.substrateIndexs_amet));
            assertEqual(1, m.substrates(strcmp(m.substrateWholeCellModelIDs,'AHCYS')));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_tRNAGuanineN7Methyltransferase));
        end

        function testTRNAGuanineN7MethylTransferase_oneTooFewSubstrate(this)
            m = this.process;

            m.reactionModificationMatrix(m.reactionModificationMatrix(:,this.MG471)==1,this.MG472)=1;
            m.reactionModificationMatrix(m.reactionModificationMatrix(:,this.MG471)==1,this.MG471)=0;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_001'),:)=0;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_001'),this.MG471)=1;
            
            m.initializeSpeciesNetwork();

            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_amet) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_tRNAGuanineN7Methyltransferase) = 1e6;
            m.unmodifiedRNAs(:) = 0;
            m.unmodifiedRNAs(this.MG471) = 1;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertEqual(1, m.unmodifiedRNAs(this.MG471));
            assertEqual(0, m.modifiedRNAs(this.MG471));
            assertEqual(0, m.substrates(m.substrateIndexs_amet));
            assertEqual(0, m.substrates(strcmp(m.substrateWholeCellModelIDs,'AHCYS')));
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_tRNAGuanineN7Methyltransferase));
        end

        function testTRNAGuanineN7MethylTransferase_insufficientEnzyme(this)
            m = this.process;

            m.reactionModificationMatrix(m.reactionModificationMatrix(:,this.MG471)==1,this.MG472)=1;
            m.reactionModificationMatrix(m.reactionModificationMatrix(:,this.MG471)==1,this.MG471)=0;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_001'),:)=0;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_001'),this.MG471)=1;
            
            m.initializeSpeciesNetwork();

            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_amet) = 1;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_tRNAGuanineN7Methyltransferase) = 0;
            m.unmodifiedRNAs(:) = 0;
            m.unmodifiedRNAs(this.MG471) = 1;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertEqual(1, m.unmodifiedRNAs(this.MG471));
            assertEqual(0, m.modifiedRNAs(this.MG471));
            assertEqual(1, m.substrates(m.substrateIndexs_amet));
            assertEqual(0, m.substrates(strcmp(m.substrateWholeCellModelIDs,'AHCYS')));
            assertEqual(0, m.enzymes(m.enzymeIndexs_tRNAGuanineN7Methyltransferase));
        end

        function testTRNAPseudouridineSynthaseA(this)
            m = this.process;

            m.reactionModificationMatrix(m.reactionModificationMatrix(:,this.MG471)==1,this.MG472)=1;
            m.reactionModificationMatrix(m.reactionModificationMatrix(:,this.MG471)==1,this.MG471)=0;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_002'),:)=0;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_002'),this.MG471)=1;
            
            m.initializeSpeciesNetwork();

            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_tRNAPseudouridineSynthase) = 1e6;
            m.unmodifiedRNAs(:) = 0;
            m.unmodifiedRNAs(this.MG471) = 1;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertEqual(0, m.unmodifiedRNAs(this.MG471));
            assertEqual(1, m.modifiedRNAs(this.MG471));
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(1e6, m.enzymes(m.enzymeIndexs_tRNAPseudouridineSynthase));
        end

        function testTRNAPseudouridineSynthaseA_insufficientEnzyme(this)
            m = this.process;

            m.reactionModificationMatrix(m.reactionModificationMatrix(:,this.MG471)==1,this.MG472)=1;
            m.reactionModificationMatrix(m.reactionModificationMatrix(:,this.MG471)==1,this.MG471)=0;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_002'),:)=0;
            m.reactionModificationMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_002'),this.MG471)=1;
            
            m.initializeSpeciesNetwork();

            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_tRNAPseudouridineSynthase) = 0;
            m.unmodifiedRNAs(:) = 0;
            m.unmodifiedRNAs(this.MG471) = 1;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertEqual(1, m.unmodifiedRNAs(this.MG471));
            assertEqual(0, m.modifiedRNAs(this.MG471));
            assertEqual(zeros(size(m.substrates)), m.substrates);
            assertEqual(0, m.enzymes(m.enzymeIndexs_tRNAPseudouridineSynthase));
        end

        function testOneOfEveryRNA(this)
            m = this.process;

            m.substrates = sum(max(0,-m.reactionStoichiometryMatrix),2);
            m.enzymes(:) = 1e6;
            m.unmodifiedRNAs(:) = 1;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertEqual(zeros(size(m.unmodifiedRNAs)), m.unmodifiedRNAs);
            assertEqual(ones(size(m.modifiedRNAs)), m.modifiedRNAs);
            assertEqual(sum(max(0,m.reactionStoichiometryMatrix),2),m.substrates);
            assertEqual(repmat(1e6,size(m.enzymes)), m.enzymes);
        end

        % Start with a lot of each of two kinds of RNAs requiring
        % same modifications, and insufficient resources to modify all RNAs
        % verify that we modify them fairly.
        function testFairAllocationOfResources(this)
            m = this.process;

            n=2000;

            m.reactionStoichiometryMatrix(:,strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_001'))=m.reactionStoichiometryMatrix(:,strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_003'));
            m.reactionStoichiometryMatrix(:,strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_002'))=m.reactionStoichiometryMatrix(:,strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_003'));
            m.reactionStoichiometryMatrix(:,strcmp(m.reactionWholeCellModelIDs,'MG472_MODIFICATION_001'))=m.reactionStoichiometryMatrix(:,strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_003'));
            m.reactionStoichiometryMatrix(:,strcmp(m.reactionWholeCellModelIDs,'MG472_MODIFICATION_002'))=m.reactionStoichiometryMatrix(:,strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_003'));

            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_001'),:)=m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_003'),:);
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_002'),:)=m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_003'),:);
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG472_MODIFICATION_001'),:)=m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_003'),:);
            m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG472_MODIFICATION_002'),:)=m.reactionCatalysisMatrix(strcmp(m.reactionWholeCellModelIDs,'MG471_MODIFICATION_003'),:);
            
            m.initializeSpeciesNetwork();

            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_atp) = n;
            m.substrates(m.substrateIndexs_cys) = n;
            m.substrates(m.substrateIndexs_water) = n;
            m.enzymes(:) = 0;
            m.enzymes(m.enzymeIndexs_tRNAUracil4Sulfurtransferase) = 1e12;
            m.unmodifiedRNAs(:) = 0;
            m.unmodifiedRNAs(this.MG471) = n;
            m.unmodifiedRNAs(this.MG472) = n;
            m.modifiedRNAs(:) = 0;

            m.evolveState();
            assertTrue(1 >= m.substrates(m.substrateIndexs_atp));
            assertTrue(1 >= m.substrates(m.substrateIndexs_cys));
            assertTrue(1 >= m.substrates(m.substrateIndexs_water));
            assertTrue(n-1 <= m.substrates(strcmp(m.substrateWholeCellModelIDs,'AMP')));
            assertTrue(n-1 <= m.substrates(strcmp(m.substrateWholeCellModelIDs,'H')));
            assertTrue(n-1 <= m.substrates(strcmp(m.substrateWholeCellModelIDs,'PPI')));
            assertTrue(n-1 <= m.substrates(strcmp(m.substrateWholeCellModelIDs,'SER')));
            assertEqual(1e12, m.enzymes(m.enzymeIndexs_tRNAUracil4Sulfurtransferase));

            modifiedRNAs = [m.modifiedRNAs(this.MG471)*3/2; m.modifiedRNAs(this.MG472)];
            assertTrue(abs(diff(modifiedRNAs)) < 0.15 * max(modifiedRNAs), ...
                sprintf('significant monomer imbalance: %d, %d', ...
                    modifiedRNAs(1), modifiedRNAs(2)));
        end

        function testGeneEssentiality(this)
            m = this.process;
            m.substrates = sum(max(0,-m.reactionStoichiometryMatrix),2);
            m.enzymes(:) = 1e12;
            m.unmodifiedRNAs(:) = 1;
            m.modifiedRNAs(:) = 0;

            this.helpTestGeneEssentiality({
                'MG_008'; %tRNA uridine 5-carboxymethylaminomethyl modification enzyme
                'MG_084'; %tRNA(Ile)-lysidine synthetase
                'MG_182'; %tRNA pseudouridine synthase A
                'MG_209'; %23S rRNA pseudouridine synthase; U955, U2504, U2580
                'MG_252'; %23S rRNA methyltransferase; G2251
                'MG_295'; %tRNA U34 sulfurtransferase
                'MG_346'; %23S rRNA methyltransferase; U2552
                'MG_347'; %tRNA (guanine-N(7)-)-methyltransferase
                'MG_370'; %23S rRNA pseudouridine synthase; U1911, U1915, U1917
                'MG_372'; %thiamine biosynthesis/tRNA modification protein ThiI
                'MG_379'; %tRNA uridine 5-carboxymethylaminomethyl modification enzyme
                'MG_380'; %methyltransferase GidB
                'MG_445'; %tRNA (guanine-N1)-methyltransferase
                'MG_463'; %dimethyladenosine transferase
                },...
                @(m,i) all((i.unmodifiedRNAs > m.unmodifiedRNAs) | i.unmodifiedRNAs==0));
        end
    end
end
