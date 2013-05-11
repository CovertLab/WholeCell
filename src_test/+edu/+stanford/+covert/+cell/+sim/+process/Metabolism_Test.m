%Metabolism test cases
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef Metabolism_Test < edu.stanford.covert.cell.sim.ReactionProcessTestCase
    %constants
    properties (Constant = true)
        expected_growthRate = 1 / (9.0 * 3600) * log(2);
    end

    %constructor
    methods
        function this = Metabolism_Test(name)
            this = this@edu.stanford.covert.cell.sim.ReactionProcessTestCase(name);
        end
    end
    
    %simple test fixture
    methods
        function testSimpleFixture(this)
            %load dummy network
            this.loadSimpleTestFixture();

            %process
            m = this.process;
            s = m.metabolicReaction;
            m.substrates(m.substrateIndexs_internalExchangedMetabolites) = 10;
            m.substrates(:, m.compartmentIndexs_extracellular) = 10;

            %test model
            m.evolveState();

            %assert growth
            assertElementsAlmostEqual(10, s.growth, 'relative', 1e-6);
        end

        function loadSimpleTestFixture(this)
            %% process
            m = this.process;

            %% numbers of components
            numSubstrates   = 25;
            numReactions    = 67;
            numCompartments = 3;

            %% IDs
            m.stimuliWholeCellModelIDs   = {};
            m.substrateWholeCellModelIDs = {...
                'glucose';'glucose-6-phosphate';'fructose-6-phosphate';...
                'fructose-1,6-diphosphate';'dihydroxyacetone phosphate';...
                'glyceraldehyde-3-phosphate';'glycerate-1,3-diphosphate';...
                '3-phosphoglycerate';'2-phosphoglycerate';'phosphoenolpyruvate';...
                'pyruvate';'6-phophogluconolactone';'6-phospho-gluconate';...
                'ribulose-5-phosphate';'ribose-5-phosphate';'xylulose-5-phosphate';...
                'ATP';'ADP';'PI';'NAD';'NADH';'H2O';'CO2';'H'};
            m.enzymeWholeCellModelIDs    = {};
            m.reactionWholeCellModelIDs = {...
                'PTS,glk';'pgi';'pfk';'fbp';'fda';'tpi';'gap';'pgk';'gpm';...
                'eno';'pyk';'zwf';'pgl';'gnd';'rpe';'rpi';'tkt,tal';...
                'TX_glucose';'TX_glucose-6-phosphate';'TX_fructose-6-phosphate';...
                'TX_fructose-1,6-diphosphate';'TX_dihydroxyacetone phosphate';...
                'TX_glyceraldehyde-3-phosphate';'TX_glycerate-1,3-diphosphate';...
                'TX_3-phosphoglycerate';'TX_2-phosphoglycerate';'TX_phosphoenolpyruvate';...
                'TX_pyruvate';'TX_6-phophogluconolactone';'TX_6-phospho-gluconate';...
                'TX_ribulose-5-phosphate';'TX_ribose-5-phosphate';'TX_xylulose-5-phosphate';...
                'TX_ATP';'TX_ADP';'TX_Pi';'TX_NAD';'TX_NADH';'TX_H2O';'TX_CO2';'TX_H'};

            %% names
            m.stimuliNames               = {};
            m.substrateNames             = m.substrateWholeCellModelIDs;
            m.enzymeNames                = {};
            m.reactionNames              = m.reactionWholeCellModelIDs;

            %% indices
            m.compartmentIndexs_cytosol          = 1;
            m.compartmentIndexs_extracellular    = 2;
            m.compartmentIndexs_membrane         = 3;

            m.substrateIndexs_amp                = [];
            m.substrateIndexs_adp                = find(strcmp(m.substrateWholeCellModelIDs, 'ADP'));
            m.substrateIndexs_atp                = find(strcmp(m.substrateWholeCellModelIDs, 'ATP'));
            m.substrateIndexs_diphosphate        = [];
            m.substrateIndexs_phosphate          = find(strcmp(m.substrateWholeCellModelIDs, 'PI'));
            m.substrateIndexs_water              = find(strcmp(m.substrateWholeCellModelIDs, 'H2O'));
            m.substrateIndexs_hydrogen           = find(strcmp(m.substrateWholeCellModelIDs, 'H'));
            [~, m.substrateIndexs_atpHydrolysis] = ismember({'ATP'; 'H2O'; 'ADP'; 'PI'; 'H'}, m.substrateWholeCellModelIDs);

            m.substrateIndexs_fba                = (1:length(m.substrateWholeCellModelIDs))';
            m.fbaSubstrateIndexs_substrates      = sub2ind([25 3], ...
                repmat((1:length(m.substrateWholeCellModelIDs))', 1, 3), ...
                repmat((1:3), 24, 1));
            m.substrateIndexs_limitableProteins = zeros(0, 1);
            m.substrateIndexs_externalExchangedMetabolites = m.substrateIndexs_fba;
            m.substrateIndexs_internalExchangedMetabolites = zeros(0, 1);
            m.substrateIndexs_internalExchangedLimitedMetabolites = zeros(0, 1);

            m.fbaSubstrateIndexs_biomass = length(m.substrateWholeCellModelIDs)+1;
            
            m.proteinLimitableProteinComposition    = zeros(0, 0);
            
            m.reactionIndexs_chemical               = (1:17)';
            m.reactionIndexs_transport              = (18:41)';
            m.reactionIndexs_fba                    = (1:41)';
            m.fbaReactionIndexs_metaboliteExternalExchange = (42:65)';
            m.fbaReactionIndexs_metaboliteInternalExchange = zeros(0, 1);
            m.fbaReactionIndexs_metaboliteInternalLimitedExchange = zeros(0, 1);
            m.fbaReactionIndexs_metaboliteInternalUnlimitedExchange = zeros(0, 1);
            m.fbaReactionIndexs_metabolicConversion = [m.reactionIndexs_chemical; m.reactionIndexs_transport];
            m.fbaReactionIndexs_biomassProduction   = max(m.fbaReactionIndexs_metaboliteExternalExchange) + 1;
            m.fbaReactionIndexs_biomassExchange     = max(m.fbaReactionIndexs_metaboliteExternalExchange) + 2;

            m.stimulusStimulusLocalIndexs               = [];
            m.stimulusStimulusGlobalIndexs              = [];
            m.stimulusStimulusCompartmentIndexs         = [];
            m.stimulusMetaboliteLocalIndexs             = [];
            m.stimulusMetaboliteGlobalIndexs            = [];
            m.stimulusMetaboliteCompartmentIndexs       = [];
            m.stimulusRNALocalIndexs                    = [];
            m.stimulusRNAGlobalIndexs                   = [];
            m.stimulusRNACompartmentIndexs              = [];
            m.stimulusMonomerLocalIndexs                = [];
            m.stimulusMonomerGlobalIndexs               = [];
            m.stimulusMonomerCompartmentIndexs          = [];
            m.stimulusComplexLocalIndexs                = [];
            m.stimulusComplexGlobalIndexs               = [];
            m.stimulusComplexCompartmentIndexs          = [];

            m.substrateStimulusLocalIndexs              = [];
            m.substrateStimulusGlobalIndexs             = [];
            m.substrateStimulusCompartmentIndexs        = [];
            m.substrateMetaboliteLocalIndexs            = (1:length(m.substrateWholeCellModelIDs))';
            m.substrateMetaboliteGlobalIndexs           = [];
            m.substrateMetaboliteCompartmentIndexs      = repmat(1:numCompartments, length(m.substrateWholeCellModelIDs),1);
            m.substrateRNALocalIndexs                   = [];
            m.substrateRNAGlobalIndexs                  = [];
            m.substrateRNACompartmentIndexs             = [];
            m.substrateMonomerLocalIndexs               = [];
            m.substrateMonomerGlobalIndexs              = [];
            m.substrateMonomerCompartmentIndexs         = [];
            m.substrateComplexLocalIndexs               = [];
            m.substrateComplexGlobalIndexs              = [];
            m.substrateComplexCompartmentIndexs         = [];         
            m.substrateMolecularWeights                 = [];

            m.enzymeStimulusLocalIndexs                 = [];
            m.enzymeStimulusGlobalIndexs                = [];
            m.enzymeStimulusCompartmentIndexs           = [];
            m.enzymeMetaboliteLocalIndexs               = [];
            m.enzymeMetaboliteGlobalIndexs              = [];
            m.enzymeMetaboliteCompartmentIndexs         = [];
            m.enzymeRNALocalIndexs                      = [];
            m.enzymeRNAGlobalIndexs                     = [];
            m.enzymeRNACompartmentIndexs                = [];
            m.enzymeMonomerLocalIndexs                  = [];
            m.enzymeMonomerGlobalIndexs                 = [];
            m.enzymeMonomerCompartmentIndexs            = [];
            m.enzymeComplexLocalIndexs                  = [];
            m.enzymeComplexGlobalIndexs                 = [];
            m.enzymeComplexCompartmentIndexs            = [];
            m.enzymeMolecularWeights                    = [];

            %% FBA stoichiometry matrix
            fbaReactionStoichiometryMatrix = zeros(numSubstrates, numReactions, numCompartments);

            %chemical reactions
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'glucose'),                    strcmp(m.reactionWholeCellModelIDs, 'PTS,glk'), m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'glucose-6-phosphate'),        strcmp(m.reactionWholeCellModelIDs, 'PTS,glk'), m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'glucose-6-phosphate'),        strcmp(m.reactionWholeCellModelIDs, 'pgi'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'glucose-6-phosphate'),        strcmp(m.reactionWholeCellModelIDs, 'zwf'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'fructose-6-phosphate'),       strcmp(m.reactionWholeCellModelIDs, 'pgi'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'fructose-6-phosphate'),       strcmp(m.reactionWholeCellModelIDs, 'pfk'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'fructose-6-phosphate'),       strcmp(m.reactionWholeCellModelIDs, 'fbp'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'fructose-6-phosphate'),       strcmp(m.reactionWholeCellModelIDs, 'tkt,tal'), m.compartmentIndexs_cytosol) =  2;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'fructose-1,6-diphosphate'),   strcmp(m.reactionWholeCellModelIDs, 'pfk'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'fructose-1,6-diphosphate'),   strcmp(m.reactionWholeCellModelIDs, 'fbp'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'fructose-1,6-diphosphate'),   strcmp(m.reactionWholeCellModelIDs, 'fda'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'dihydroxyacetone phosphate'), strcmp(m.reactionWholeCellModelIDs, 'fda'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'dihydroxyacetone phosphate'), strcmp(m.reactionWholeCellModelIDs, 'tpi'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'glyceraldehyde-3-phosphate'), strcmp(m.reactionWholeCellModelIDs, 'tpi'),     m.compartmentIndexs_cytosol) =  2;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'glyceraldehyde-3-phosphate'), strcmp(m.reactionWholeCellModelIDs, 'gap'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'glyceraldehyde-3-phosphate'), strcmp(m.reactionWholeCellModelIDs, 'tkt,tal'), m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'glycerate-1,3-diphosphate'),  strcmp(m.reactionWholeCellModelIDs, 'gap'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'glycerate-1,3-diphosphate'),  strcmp(m.reactionWholeCellModelIDs, 'pgk'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'3-phosphoglycerate'),         strcmp(m.reactionWholeCellModelIDs, 'pgk'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'3-phosphoglycerate'),         strcmp(m.reactionWholeCellModelIDs, 'gpm'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'2-phosphoglycerate'),         strcmp(m.reactionWholeCellModelIDs, 'gpm'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'2-phosphoglycerate'),         strcmp(m.reactionWholeCellModelIDs, 'eno'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'phosphoenolpyruvate'),        strcmp(m.reactionWholeCellModelIDs, 'PTS,glk'), m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'phosphoenolpyruvate'),        strcmp(m.reactionWholeCellModelIDs, 'eno'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'phosphoenolpyruvate'),        strcmp(m.reactionWholeCellModelIDs, 'pyk'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'pyruvate'),                   strcmp(m.reactionWholeCellModelIDs, 'PTS,glk'), m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'pyruvate'),                   strcmp(m.reactionWholeCellModelIDs, 'pyk'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'6-phophogluconolactone'),     strcmp(m.reactionWholeCellModelIDs, 'zwf'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'6-phophogluconolactone'),     strcmp(m.reactionWholeCellModelIDs, 'pgl'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'6-phospho-gluconate'),        strcmp(m.reactionWholeCellModelIDs, 'pgl'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'6-phospho-gluconate'),        strcmp(m.reactionWholeCellModelIDs, 'gnd'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ribulose-5-phosphate'),       strcmp(m.reactionWholeCellModelIDs, 'gnd'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ribulose-5-phosphate'),       strcmp(m.reactionWholeCellModelIDs, 'rpe'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ribulose-5-phosphate'),       strcmp(m.reactionWholeCellModelIDs, 'rpi'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ribose-5-phosphate'),         strcmp(m.reactionWholeCellModelIDs, 'rpi'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ribose-5-phosphate'),         strcmp(m.reactionWholeCellModelIDs, 'tkt,tal'), m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'xylulose-5-phosphate'),       strcmp(m.reactionWholeCellModelIDs, 'rpe'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'xylulose-5-phosphate'),       strcmp(m.reactionWholeCellModelIDs, 'tkt,tal'), m.compartmentIndexs_cytosol) = -2;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ATP'),                        strcmp(m.reactionWholeCellModelIDs, 'PTS,glk'), m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ATP'),                        strcmp(m.reactionWholeCellModelIDs, 'pfk'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ATP'),                        strcmp(m.reactionWholeCellModelIDs, 'pgk'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ATP'),                        strcmp(m.reactionWholeCellModelIDs, 'pyk'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ADP'),                        strcmp(m.reactionWholeCellModelIDs, 'PTS,glk'), m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ADP'),                        strcmp(m.reactionWholeCellModelIDs, 'pfk'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ADP'),                        strcmp(m.reactionWholeCellModelIDs, 'pgk'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'ADP'),                        strcmp(m.reactionWholeCellModelIDs, 'pyk'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'PI'),                         strcmp(m.reactionWholeCellModelIDs, 'fbp'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'PI'),                         strcmp(m.reactionWholeCellModelIDs, 'gap'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'NAD'),                        strcmp(m.reactionWholeCellModelIDs, 'gap'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'NADH'),                       strcmp(m.reactionWholeCellModelIDs, 'gap'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'H2O'),                        strcmp(m.reactionWholeCellModelIDs, 'eno'),     m.compartmentIndexs_cytosol) =  1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'H2O'),                        strcmp(m.reactionWholeCellModelIDs, 'pgl'),     m.compartmentIndexs_cytosol) = -1;
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs,'CO2'),                        strcmp(m.reactionWholeCellModelIDs, 'gnd'),     m.compartmentIndexs_cytosol) =  1;

            %transport reactions
            fbaReactionStoichiometryMatrix(1:24, m.reactionIndexs_transport, m.compartmentIndexs_extracellular) = ...
                -eye(length(m.substrateWholeCellModelIDs));
            fbaReactionStoichiometryMatrix(1:24, m.reactionIndexs_transport, m.compartmentIndexs_cytosol) = ...
                eye(length(m.substrateWholeCellModelIDs));

            %exchange reactions
            fbaReactionStoichiometryMatrix(1:24, m.fbaReactionIndexs_metaboliteExternalExchange, m.compartmentIndexs_extracellular) = ...
                eye(length(m.substrateWholeCellModelIDs));

            %biomass production
            fbaReactionStoichiometryMatrix(strcmp(m.substrateWholeCellModelIDs, 'pyruvate'), m.fbaReactionIndexs_biomassProduction, m.compartmentIndexs_cytosol) = ...
                -1;
            fbaReactionStoichiometryMatrix(m.fbaSubstrateIndexs_biomass, m.fbaReactionIndexs_biomassProduction, m.compartmentIndexs_cytosol) = ...
                1;

            %biomass efflux
            fbaReactionStoichiometryMatrix(m.fbaSubstrateIndexs_biomass, m.fbaReactionIndexs_biomassExchange, m.compartmentIndexs_cytosol) = ...
                -1;

            %reshape
            m.fbaReactionStoichiometryMatrix = reshape(...
                permute(fbaReactionStoichiometryMatrix, [2 1 3]),...
                size(fbaReactionStoichiometryMatrix,2), [])';

            %% FBA reaction bounds
            m.fbaReactionBounds = zeros(numReactions, 2);

            %chemical reactions
            m.fbaReactionBounds(m.reactionIndexs_chemical, :) = repmat([-Inf Inf],length(m.reactionIndexs_chemical),1);
            m.fbaReactionBounds(strcmp(m.reactionWholeCellModelIDs, 'PTS,glk'),1) = 0;
            m.fbaReactionBounds(strcmp(m.reactionWholeCellModelIDs, 'pfk'),    1) = 0;
            m.fbaReactionBounds(strcmp(m.reactionWholeCellModelIDs, 'fbp'),    1) = 0;
            m.fbaReactionBounds(strcmp(m.reactionWholeCellModelIDs, 'pyk'),    1) = 0;
            m.fbaReactionBounds(strcmp(m.reactionWholeCellModelIDs, 'zwf'),    1) = 0;
            m.fbaReactionBounds(strcmp(m.reactionWholeCellModelIDs, 'pgl'),    1) = 0;
            m.fbaReactionBounds(strcmp(m.reactionWholeCellModelIDs, 'gnd'),    1) = 0;

            %transport reactions
            m.fbaReactionBounds(m.reactionIndexs_transport, :) = repmat([-Inf Inf], length(m.reactionIndexs_transport),1);

            %exchange reactions
            m.fbaReactionBounds(m.fbaReactionIndexs_metaboliteExternalExchange(strcmp(m.substrateWholeCellModelIDs, 'glucose')),2) =  Inf;
            m.fbaReactionBounds(m.fbaReactionIndexs_metaboliteExternalExchange(strcmp(m.substrateWholeCellModelIDs, 'ATP')),    :) = [-Inf Inf];
            m.fbaReactionBounds(m.fbaReactionIndexs_metaboliteExternalExchange(strcmp(m.substrateWholeCellModelIDs, 'ADP')),    :) = [-Inf Inf];
            m.fbaReactionBounds(m.fbaReactionIndexs_metaboliteExternalExchange(strcmp(m.substrateWholeCellModelIDs, 'PI')),     2) =  Inf;
            m.fbaReactionBounds(m.fbaReactionIndexs_metaboliteExternalExchange(strcmp(m.substrateWholeCellModelIDs, 'NAD')),    2) =  Inf;
            m.fbaReactionBounds(m.fbaReactionIndexs_metaboliteExternalExchange(strcmp(m.substrateWholeCellModelIDs, 'NADH')),   1) = -Inf;
            m.fbaReactionBounds(m.fbaReactionIndexs_metaboliteExternalExchange(strcmp(m.substrateWholeCellModelIDs, 'H2O')),    :) = [-Inf Inf];
            m.fbaReactionBounds(m.fbaReactionIndexs_metaboliteExternalExchange(strcmp(m.substrateWholeCellModelIDs, 'CO2')),    1) = -Inf;
            m.fbaReactionBounds(m.fbaReactionIndexs_metaboliteExternalExchange(strcmp(m.substrateWholeCellModelIDs, 'H')),      2) =  Inf;

            %biomass production, efflux
            m.fbaReactionBounds(m.fbaReactionIndexs_biomassProduction, :) = [0 Inf];
            m.fbaReactionBounds(m.fbaReactionIndexs_biomassExchange, :) = [-Inf Inf];

            %% other FBA properties
            m.fbaEnzymeBounds                = NaN(numReactions,2);
            m.fbaObjective                   = [zeros(numReactions-1,1); 1];
            m.fbaRightHandSide               = zeros(numSubstrates*numCompartments,1);
            m.fbaReactionCatalysisMatrix     = zeros(numReactions,0);
            m.metabolismProduction           = -permute(fbaReactionStoichiometryMatrix(1:24, m.fbaReactionIndexs_biomassProduction, :), [1 3 2]);
            m.metabolismNewProduction        = m.metabolismProduction;
            m.metabolismRecyclingProduction  = zeros(size(m.metabolismProduction));
            m.unaccountedEnergyConsumption   = 0;

            %% reaction m properties
            m.reactionStoichiometryMatrix = fbaReactionStoichiometryMatrix(:,m.fbaReactionIndexs_metabolicConversion,:);
            m.reactionCatalysisMatrix     = m.fbaReactionCatalysisMatrix(m.fbaReactionIndexs_metabolicConversion,:);
            m.reactionBounds              = m.fbaReactionBounds(m.fbaReactionIndexs_metabolicConversion,:);
            m.enzymeBounds                = m.fbaEnzymeBounds(m.fbaReactionIndexs_metabolicConversion,:);

            %% amounts of components
            m.stimuli    = zeros(0, numCompartments);
            m.substrates = zeros(length(m.substrateWholeCellModelIDs), numCompartments);
            m.enzymes    = zeros(0,1);
            m.metabolicReaction.fluxs = zeros(length(m.reactionWholeCellModelIDs),1);            
        end
    end

    %tests
    methods
        function testConstants(this)
            m = this.process;
            
            assertEqual(3, size(m.reactionStoichiometryMatrix, 3));
            assertAllEqual(false, isnan(m.reactionBounds));
            assertAllEqual(false, isnan(m.fbaReactionBounds));
            
            assertAllEqual(true, isfinite(m.substrateExternalExchangeBounds));
            assertAllEqual(true, isfinite(m.fbaReactionBounds(m.fbaReactionIndexs_metaboliteExternalExchange, :)));
            
            assertIn(m.metabolismRecyclingProduction(m.substrateIndexs_internalExchangedLimitedMetabolites), ...
                repmat([-Inf 0], size(m.substrateIndexs_internalExchangedLimitedMetabolites)));
            assertIn(m.metabolismRecyclingProduction(setdiff(m.substrateIndexs_internalExchangedMetabolites, m.substrateIndexs_internalExchangedLimitedMetabolites)), ...
                repmat([0 Inf], size(setdiff(m.substrateIndexs_internalExchangedMetabolites, m.substrateIndexs_internalExchangedLimitedMetabolites))));
            assertElementsAlmostEqual(m.metabolismNewProduction + m.metabolismRecyclingProduction, ...
                m.metabolismProduction, 'relative', 1e-6, 1e-6);
        end
        
        function testWildTypeGrowth(this)
            m = this.process;
            s = m.metabolicReaction;
            m.evolveState();

            %assert growth rate
            assertElementsAlmostEqual(this.expected_growthRate, s.growth, 'relative', 0.50);

            %substrates all positive
            assertFalse(any(any((m.substrates .* m.metabolismProduction > 0) < 0)));
            assertFalse(any(any(m.substrates(setdiff(1:size(m.substrates, 1), m.substrateMetaboliteLocalIndexs), :) < 0)));
        end

        %no substrates
        function testNoSubstrates(this)
            m = this.process;
            s = m.metabolicReaction;
            m.substrates(:) = 0;
            m.evolveState();
            assertElementsAlmostEqual(0, s.growth, 'absolute', 1e-8);
        end

        %no extracellular substrates
        function testNoMedia(this)
            m = this.process;
            s = m.metabolicReaction;
            m.substrates(:,m.compartmentIndexs_extracellular) = 0;
            m.evolveState();
            assertElementsAlmostEqual(0, s.growth, 'absolute', 1e-8);
        end

        %no proteins substrates
        function testNoProteinSubstrates(this)
            m = this.process;
            s = m.metabolicReaction;
            m.substrates(m.substrateIndexs_limitableProteins, :) = 0;
            m.evolveState();

            assertElementsAlmostEqual(0, s.growth, 'absolute', 1e-8);
            reactionIdxs = find(sum(sum(abs(m.reactionStoichiometryMatrix(m.substrateIndexs_limitableProteins,:,:)),3),1));
            assertElementsAlmostEqual(zeros(length(reactionIdxs),1), ...
                s.fluxs(reactionIdxs), 'absolute', 1e-5);
        end

        %enzymes all zero
        function testNoEnzymes(this)
            m = this.process;
            s = m.metabolicReaction;
            m.enzymes(:) = 0;
            m.evolveState();

            %growth
            assertElementsAlmostEqual(0, s.growth, 'absolute', 1e-8);

            %reaction fluxes
            reactionIndexs=~isinf(m.enzymeBounds(:,1));
            assertFalse(any(s.fluxs(reactionIndexs)<-1e-4));

            reactionIndexs=~isinf(m.enzymeBounds(:,2));
            assertFalse(any(s.fluxs(reactionIndexs)>1e-4));
        end

        %set all enzyme kinetic rates to zero
        function testNoCatalysis(this)
            m = this.process;
            s = m.metabolicReaction;
            m.fbaEnzymeBounds(~isnan(m.fbaEnzymeBounds))=0;
            m.evolveState();

            %growth
            assertElementsAlmostEqual(0, s.growth, 'absolute', 1e-8);

            %reaction fluxes
            reactionIndexs=~isinf(m.enzymeBounds(:,1));
            assertFalse(any(s.fluxs(reactionIndexs)<-1e-4));

            reactionIndexs=~isinf(m.enzymeBounds(:,2));
            assertFalse(any(s.fluxs(reactionIndexs)>1e-4));
        end

        function testNoReactions(this)
            m = this.process;
            s = m.metabolicReaction;
            m.fbaReactionBounds(:)=0;
            m.evolveState();

            %growth
            assertElementsAlmostEqual(0, s.growth, 'absolute', 1e-8);

            %fluxes
            assertElementsAlmostEqual(zeros(size(s.fluxs)),...
                s.fluxs,...
                'absolute',1e-4);
        end

        %turn off exchange reactions
        function testNoExternalExchangeReactions(this)
            m = this.process;
            s = m.metabolicReaction;
            m.fbaReactionBounds(m.fbaReactionIndexs_metaboliteExternalExchange, :) = 0;
            m.evolveState();

            %growth
            assertElementsAlmostEqual(0, s.growth, 'absolute', 1e-8);
        end

        %turn off transport reactions
        function testNoTransportReactions(this)
            m = this.process;
            s = m.metabolicReaction;
            m.fbaReactionBounds(m.fbaReactionIndexs_metabolicConversion(ismember(m.reactionIndexs_fba, m.reactionIndexs_transport)), :) = 0;
            m.evolveState();

            %growth
            assertElementsAlmostEqual(0, s.growth, 'absolute', 1e-8);

            %fluxes
            assertElementsAlmostEqual(zeros(size(m.reactionIndexs_transport)),...
                s.fluxs(m.reactionIndexs_transport),...
                'absolute',1e-4);
        end

        %turn off chemical reactions
        function testNoChemicalReactions(this)
            m = this.process;
            s = m.metabolicReaction;
            m.fbaReactionBounds(m.reactionIndexs_chemical,:)=0;
            m.evolveState();

            %growth
            assertElementsAlmostEqual(0, s.growth, 'absolute', 1e-8);

            %fluxes
            assertElementsAlmostEqual(zeros(size(m.reactionIndexs_chemical)),...
                s.fluxs(m.reactionIndexs_chemical),...
                'absolute',1e-4);
        end       

        %check that all metabolites in objective can be made individually
        function testObjectiveMetabolitesCanBeIndividuallyProduced(this)
            m = this.process;

            substrateIndexs = find(m.fbaReactionStoichiometryMatrix(m.fbaSubstrateIndexs_substrates, m.fbaReactionIndexs_biomassProduction) < 0);
            [j, ~] = ind2sub(size(m.substrates), m.substrateIndexs_fba(substrateIndexs));
            tfs = ~ismember(j, m.substrateIndexs({'MTHF', 'METTHF', 'FTHF5', 'FTHF10'}));
            substrateIndexs = substrateIndexs(tfs);
            j = j(tfs);
            
            for i = 1:numel(substrateIndexs)
                m.fbaReactionStoichiometryMatrix(:, m.fbaReactionIndexs_biomassProduction) = 0;
                m.fbaReactionStoichiometryMatrix(m.fbaSubstrateIndexs_substrates(substrateIndexs(i)), m.fbaReactionIndexs_biomassProduction) = -1;
                m.fbaReactionStoichiometryMatrix(m.fbaSubstrateIndexs_biomass, m.fbaReactionIndexs_biomassProduction) = 1;
                
                growth = m.calcGrowthRate(m.calcFluxBounds(m.substrates, m.enzymes, m.fbaReactionBounds, m.fbaEnzymeBounds, true, true, true, true, false, true));
                assertTrue(growth > 0, ...
                    sprintf('Cannot produce substrate %s, growth: %f', ...
                    m.substrateWholeCellModelIDs{j(i)}, ...
                    growth));
            end
        end
        
        %check media components uptaken
        function testMetaboliteUptake(this)
            m = this.process;
            cmIdx = [m.compartmentIndexs_cytosol; m.compartmentIndexs_membrane];
            eIdx = m.compartmentIndexs_extracellular;
            
            m.substrates(:, cmIdx) = 0;
            substrates = m.substrates;
            
            m.evolveState();
            
            tmp = substrates + m.metabolismNewProduction * m.metabolicReaction.growth * m.stepSizeSec;
            tmp(m.substrateIndexs_atpHydrolysis) = ...
                tmp(m.substrateIndexs_atpHydrolysis) + ...
                [-1; -1; 1; 1; 1] * ...
                m.unaccountedEnergyConsumption * m.metabolicReaction.growth * m.stepSizeSec;
            
            assertElementsAlmostEqual(tmp(:, cmIdx), m.substrates(:, cmIdx), ...
                'relative', 1e-3, 2);
            assertEqual(substrates(setdiff(1:end, m.substrateIndexs_externalExchangedMetabolites), eIdx), ...
                m.substrates(setdiff(1:end, m.substrateIndexs_externalExchangedMetabolites), eIdx));
            assertEqual(substrates(setdiff(1:end, m.substrateMetaboliteLocalIndexs), :), ...
                m.substrates(setdiff(1:end, m.substrateMetaboliteLocalIndexs), :));
        end
        
        function testInternalExchangeBounds(this)
            m = this.process;
            substrates = m.substrates;
            
            assertIn(numel(m.substrateIndexs_internalExchangedMetabolites), [1 Inf]);
            
            %plenty internal metabolites
            substrates(m.substrateIndexs_internalExchangedMetabolites) = 1e8;
            [growth, ~, fbaReactionFluxs] = m.calcGrowthRate(m.calcFluxBounds(substrates, m.enzymes, m.fbaReactionBounds, m.fbaEnzymeBounds, true, true, true, true, true, true));
            assertElementsAlmostEqual(this.expected_growthRate, growth, 'relative', 0.50);
            assertTrue(any(fbaReactionFluxs(m.fbaReactionIndexs_metaboliteInternalExchange)));
            
            %missing one internal metabolite
            idx = find(m.substrateIndexs_internalExchangedLimitedMetabolites(1) == m.substrateIndexs_internalExchangedMetabolites);
            substrates(m.substrateIndexs_internalExchangedMetabolites(idx)) = 0;
            [growth, ~, fbaReactionFluxs] = m.calcGrowthRate(m.calcFluxBounds(substrates, m.enzymes, m.fbaReactionBounds, m.fbaEnzymeBounds, true, true, true, true, true, true));
            assertElementsAlmostEqual(this.expected_growthRate, growth, 'relative', 0.50);
            assertEqual(0, fbaReactionFluxs(m.fbaReactionIndexs_metaboliteInternalExchange(idx)));
            
            %no internal metabolites
            substrates(m.substrateIndexs_internalExchangedLimitedMetabolites) = 0;
            [growth, ~, fbaReactionFluxs] = m.calcGrowthRate(m.calcFluxBounds(substrates, m.enzymes, m.fbaReactionBounds, m.fbaEnzymeBounds, true, true, true, true, true, true));
            assertElementsAlmostEqual(this.expected_growthRate, growth, 'relative', 0.50);
            assertAllEqual(0, fbaReactionFluxs(m.fbaReactionIndexs_metaboliteInternalExchange));
        end
        
        function testNTPRecycling(this)
            m = this.process;
            mr = m.metabolicReaction;
            %sim.applyOptions('verbosity', 1);
            %this.seedSimulation(sim, 1);
            
            iterMax = 30000;
            initEnzymes = m.enzymes;
            ntps = zeros(4, iterMax);
            ndps = zeros(4, iterMax);
            nmps = zeros(4, iterMax);
            ntpIdxs = m.substrateIndexs({'ATP'; 'CTP'; 'GTP'; 'UTP'});
            ndpIdxs = m.substrateIndexs({'ADP'; 'CDP'; 'GDP'; 'UDP'});
            nmpIdxs = m.substrateIndexs({'AMP'; 'CMP'; 'GMP'; 'UMP'});
            warning('off', 'WholeCell:warning:negativeMetabolites');
            for i = 1:iterMax
                %mock protein synthesis
                m.enzymes = ...
                    + m.enzymes ...
                    + m.randStream.stochasticRound(initEnzymes * log(2) / m.cellCycleLength * exp(log(2) * i / m.cellCycleLength));
                
                %simulate
                lastwarn('');
                m.evolveState();
                if i > 50 && ~isequal('', lastwarn())
                    break;
                end
                
                %mock other processes
                m.substrates = ...
                    + m.substrates ...
                    - m.randStream.stochasticRound(2 * mr.growth * m.randStream.rand() * m.metabolismRecyclingProduction);
                
                %track metabolic counts
                ntps(:, i) = m.substrates(ntpIdxs, m.compartmentIndexs_cytosol);
                nmps(:, i) = m.substrates(ndpIdxs, m.compartmentIndexs_cytosol);
                ndps(:, i) = m.substrates(nmpIdxs, m.compartmentIndexs_cytosol);
            end
            
            setWarnings;
            assertEqual('', lastwarn());

            assertIn(min(ndps(:, end)), [0 200]);
            assertIn(min(nmps(:, end)), [0 200]);
            
            %% plot
            %figHandle = figure();
            %clf(figHandle);
            %time = (1:iterMax)/3600;
            %
            %subplot(3, 1, 1);
            %h = plot(time, ntps);
            %legend(h, {'ATP'; 'CTP'; 'GTP'; 'UTP'}, 'Location', 'EastOutside');
            %
            %subplot(3, 1, 2);
            %h = plot(time, ndps);
            %legend(h, {'ADP'; 'CDP'; 'GDP'; 'UDP'}, 'Location', 'EastOutside');
            %
            %subplot(3, 1, 3);
            %h = plot(time, nmps);
            %legend(h, {'AMP'; 'CMP'; 'GMP'; 'UMP'}, 'Location', 'EastOutside');
        end
        
        function testGeneEssentiality(this)
            this.helpTestGeneEssentiality({
                'MG_006'; 'MG_013'; 'MG_023'; 'MG_034'; 'MG_037'; 'MG_038'; 'MG_041'; 
                'MG_042'; 'MG_043'; 'MG_044'; 'MG_045'; 'MG_047';
                'MG_053'; 'MG_058'; 'MG_066'; 'MG_069'; 'MG_071'; 
                'MG_077'; 'MG_078'; 'MG_079'; 'MG_080'; 'MG_102'; 'MG_107'; 
                'MG_111'; 'MG_112'; 'MG_114'; 'MG_118'; 'MG_124'; 'MG_128'; 
                'MG_137'; 'MG_145'; 'MG_171'; 'MG_179'; 'MG_180'; 'MG_181'; 
                'MG_212'; 'MG_215'; 'MG_216'; 'MG_228'; 'MG_229'; 'MG_230'; 
                'MG_231'; 'MG_245'; 'MG_270'; 'MG_271'; 'MG_272'; 'MG_273'; 
                'MG_274'; 'MG_275'; 'MG_276'; 'MG_278'; 'MG_287'; 'MG_299'; 
                'MG_300'; 'MG_301'; 'MG_302'; 'MG_303'; 'MG_304'; 'MG_321'; 
                'MG_322'; 'MG_323'; 'MG_330'; 'MG_342'; 'MG_351'; 'MG_357'; 
                'MG_368'; 'MG_382'; 'MG_383'; 'MG_394'; 'MG_396'; 'MG_407'; 
                'MG_429'; 'MG_430'; 'MG_431'; 'MG_434'; 'MG_437'; 'MG_453'; 
                'MG_458'; 'MG_517';}, ...
                @(m,i) m.metabolicReaction.growth > this.expected_growthRate / 5);
        end
    end
end
