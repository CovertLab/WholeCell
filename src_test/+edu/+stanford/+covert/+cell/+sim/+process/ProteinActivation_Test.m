%Protein activation process test case
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef ProteinActivation_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = ProteinActivation_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end

        function testActivationRuleTranscoding(this)
            m = this.process;

            A = 1; %#ok<NASGU>
            B = 2; %#ok<NASGU>

            %database syntax, transcoded syntax
            rules = {
                'A+B','A+B';
                'A-B','A-B';
                'A|B','A|B';
                'A&B','A&B';
                '!(A+B)','~(A+B)';
                'A>B','A>B';
                'A<B','A<B';
                'A>=B','A>=B';
                'A<=B','A<=B';
                'A==B','A==B'};

            for i = 1:size(rules,1)
                r = m.transcodeActivationRules(rules(i,1));
                assertEqual(rules{i,2}, r{1});
                [~] = eval(r{1});
            end
        end

        %tests scaling stimuli, substrates, and enzymes to concentrations in mM
        function testScaleState(this)
            m = this.process;

            nCmp = length(m.compartment.wholeCellModelIDs);

            %stimuli
            m.stimuliWholeCellModelIDs   = cellstr(('ABCDEFGHIJ')');
            m.substrateWholeCellModelIDs = cellstr(('ABCDEFGHIJ')');
            m.enzymeWholeCellModelIDs    = cellstr(('ABCDEFGHIJ')');

            m.stimulusStimulusLocalIndexs               = (1:2)';
            m.stimulusMetaboliteLocalIndexs             = (3:4)';
            m.stimulusRNALocalIndexs                    = (5:6)';
            m.stimulusMonomerLocalIndexs                = (7:8)';
            m.stimulusComplexLocalIndexs                = (9:10)';

            m.substrateStimulusLocalIndexs              = (1:2)';
            m.substrateMetaboliteLocalIndexs            = (3:4)';
            m.substrateRNALocalIndexs                   = (5:6)';
            m.substrateMonomerLocalIndexs               = (7:8)';
            m.substrateComplexLocalIndexs               = (9:10)';

            m.enzymeStimulusLocalIndexs                 = (1:2)';
            m.enzymeMetaboliteLocalIndexs               = (3:4)';
            m.enzymeRNALocalIndexs                      = (5:6)';
            m.enzymeMonomerLocalIndexs                  = (7:8)';
            m.enzymeComplexLocalIndexs                  = (9:10)';

            m.stimulusStimulusGlobalIndexs              = repmat((1:2)', 1, nCmp);
            m.stimulusMetaboliteGlobalIndexs            = repmat((1:2)', 1, nCmp);
            m.stimulusRNAGlobalIndexs                   = repmat((1:2)', 1, nCmp);
            m.stimulusMonomerGlobalIndexs               = repmat((1:2)', 1, nCmp);
            m.stimulusComplexGlobalIndexs               = repmat((1:2)', 1, nCmp);

            m.substrateStimulusGlobalIndexs             = repmat((1:2)', 1, nCmp);
            m.substrateMetaboliteGlobalIndexs           = repmat((1:2)', 1, nCmp);
            m.substrateRNAGlobalIndexs                  = repmat((1:2)', 1, nCmp);
            m.substrateMonomerGlobalIndexs              = repmat((1:2)', 1, nCmp);
            m.substrateComplexGlobalIndexs              = repmat((1:2)', 1, nCmp);

            m.enzymeStimulusGlobalIndexs                = repmat((1:2)', 1, nCmp);
            m.enzymeMetaboliteGlobalIndexs              = repmat((1:2)', 1, nCmp);
            m.enzymeRNAGlobalIndexs                     = repmat((1:2)', 1, nCmp);
            m.enzymeMonomerGlobalIndexs                 = repmat((1:2)', 1, nCmp);
            m.enzymeComplexGlobalIndexs                 = repmat((1:2)', 1, nCmp);

            m.stimulusStimulusCompartmentIndexs         = repmat(1:nCmp, 2, 1);
            m.stimulusMetaboliteCompartmentIndexs       = repmat(1:nCmp, 2, 1);
            m.stimulusRNACompartmentIndexs              = repmat(1:nCmp, 2, 1);
            m.stimulusMonomerCompartmentIndexs          = repmat(1:nCmp, 2, 1);
            m.stimulusComplexCompartmentIndexs          = repmat(1:nCmp, 2, 1);

            m.substrateStimulusCompartmentIndexs        = repmat(1:nCmp, 2, 1);
            m.substrateMetaboliteCompartmentIndexs      = repmat(1:nCmp, 2, 1);
            m.substrateRNACompartmentIndexs             = repmat(1:nCmp, 2, 1);
            m.substrateMonomerCompartmentIndexs         = repmat(1:nCmp, 2, 1);
            m.substrateComplexCompartmentIndexs         = repmat(1:nCmp, 2, 1);

            m.enzymeStimulusCompartmentIndexs           = repmat(1:nCmp, 2, 1);
            m.enzymeMetaboliteCompartmentIndexs         = repmat(1:nCmp, 2, 1);
            m.enzymeRNACompartmentIndexs                = repmat(1:nCmp, 2, 1);
            m.enzymeMonomerCompartmentIndexs            = repmat(1:nCmp, 2, 1);
            m.enzymeComplexCompartmentIndexs            = repmat(1:nCmp, 2, 1);

            m.geometry.volume = 2e-21;

            m.stimuli = 1e6 * rand(length(m.stimuliWholeCellModelIDs), nCmp);
            m.substrates = 1e6 * rand(length(m.stimuliWholeCellModelIDs), nCmp);
            m.enzymes = 1e6 * rand(length(m.stimuliWholeCellModelIDs), nCmp);

            nAvogadro = edu.stanford.covert.util.ConstantUtil.nAvogadro;

            scaledStimuli = m.stimuli / m.geometry.volume / nAvogadro * 1e3;
            scaledStimuli(m.stimulusStimulusLocalIndexs,:) = m.stimuli(m.stimulusStimulusLocalIndexs,:);

            scaledSubstrates = m.substrates / m.geometry.volume / nAvogadro * 1e3;
            scaledSubstrates(m.substrateStimulusLocalIndexs,:) = m.substrates(m.substrateStimulusLocalIndexs,:);

            scaledEnzymes = m.enzymes / m.geometry.volume / nAvogadro * 1e3;
            scaledEnzymes(m.enzymeStimulusLocalIndexs,:) = m.enzymes(m.enzymeStimulusLocalIndexs,:);

            %scale stimuli
            m.scaleState('stimulus',  'stimuli',    'mM');
            m.scaleState('substrate', 'substrates', 'mM');
            m.scaleState('enzyme',    'enzymes',    'mM');

            %assert stimuli were scaled properly
            assertElementsAlmostEqual(scaledStimuli,    m.stimuli,    'relative', 1e12, 'Stimuli incorrectly scaled');
            assertElementsAlmostEqual(scaledSubstrates, m.substrates, 'relative', 1e12, 'Substrates incorrectly scaled');
            assertElementsAlmostEqual(scaledEnzymes,    m.enzymes,    'relative', 1e12, 'Enzymes incorrectly scaled');
        end

        function testActivationRuleEvaluation(this)
            this.loadSimpleTestFixture();

            nAvogadro = edu.stanford.covert.util.ConstantUtil.nAvogadro;            
            m = this.process;
            
            %% example 1
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs,'G6P'),                  :) = 6 * m.geometry.volume * nAvogadro*1e-3;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs,'PI'),                   :) = 22 * m.geometry.volume * nAvogadro*1e-3;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs,'stimulus_gluconate'),   :) = 1;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs,'stimulus_ironStress'),  :) = 0;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs,'stimulus_thiolStress'), :) = 1;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs,'temperature'),          :) = 37;

            m.substrates(:) = 1;
            m.inactivatedSubstrates(:) = 1;

            %assert the evaluation of the activation rules
            ruleEvaluations = logical([
                0 1 0 1 0 1;
                0 1 0 1 0 1])';
            assertEqual(ruleEvaluations, m.evaluateActivationRules());
            
            %% example 2
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs,'G6P'),                  :) = 4 * m.geometry.volume * nAvogadro*1e-3;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs,'PI'),                   :) = 18 * m.geometry.volume * nAvogadro*1e-3;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs,'stimulus_gluconate'),   :) = 0;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs,'stimulus_ironStress'),  :) = 1;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs,'stimulus_thiolStress'), :) = 1;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs,'temperature'),          :) = 43;

            m.substrates(:) = 1;
            m.inactivatedSubstrates(:) = 1;

            %assert the evaluation of the activation rules
            ruleEvaluations = logical([
                1 1 1 0 1 0
                1 1 1 0 1 0])';
            assertEqual(ruleEvaluations, m.evaluateActivationRules());
        end

        %test that proteins can be activated and inactivated
        function testProteinRegulation(this)
            nAvogadro = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            m = this.process;
            sIdxs = m.substrateIndexs({
                'MG_101_MONOMER'
                'MG_127_MONOMER'
                'MG_236_MONOMER'
                'MG_085_HEXAMER'
                'MG_205_DIMER'
                'MG_409_DIMER'
                });
            
            %example 1
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs, 'G6P'),                  :) = 4 * m.geometry.volume * nAvogadro*1e-3;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs, 'PI'),                   :) = 18 * m.geometry.volume * nAvogadro*1e-3;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs, 'stimulus_gluconate'),   :) = 0;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs, 'stimulus_ironStress'),  :) = 1;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs, 'stimulus_thiolStress'), :) = 1;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs, 'temperature'),          :) = 43;

            m.substrates(:, 1:2) = 1;
            m.inactivatedSubstrates(:, 1:2) = 1;

            initial_stimuli = m.stimuli;

            m.evolveState();
            ruleEvaluations = logical([
                1 1 1 0 1 0
                1 1 1 0 1 0])';
            assertEqual(initial_stimuli, m.stimuli); %assert stimuli values haven't changed
            assertEqual(2 * ruleEvaluations, m.substrates(sIdxs, 1:2));
            assertEqual(2 * (1 - ruleEvaluations), m.inactivatedSubstrates(sIdxs, 1:2));
            
            %example 2
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs, 'G6P'),                  :) = 6 * m.geometry.volume * nAvogadro*1e-3;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs, 'PI'),                   :) = 22 * m.geometry.volume * nAvogadro*1e-3;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs, 'stimulus_gluconate'),   :) = 1;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs, 'stimulus_ironStress'),  :) = 0;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs, 'stimulus_thiolStress'), :) = 1;
            m.stimuli(strcmp(m.stimuliWholeCellModelIDs, 'temperature'),          :) = 37;

            m.substrates(:, 1:2) = 1;
            m.inactivatedSubstrates(:, 1:2) = 1;

            initial_stimuli = m.stimuli;

            m.evolveState();
            ruleEvaluations = logical([
                0 1 0 1 0 1;
                0 1 0 1 0 1])';
            assertEqual(initial_stimuli, m.stimuli); %assert stimuli values haven't changed
            assertEqual(2 * ruleEvaluations, m.substrates(sIdxs, 1:2));
            assertEqual(2 * (1 - ruleEvaluations), m.inactivatedSubstrates(sIdxs, 1:2));
        end

        function testGeneEssentiality(this)
            m = this.process;
            m.geometry.volume = 1e-21;
            m.substrates(:) = 0;
            m.inactivatedSubstrates(:) = 1;

            this.helpTestGeneEssentiality({}, @(m,i) any(m.substrates(:)));
        end
    end
    
    methods (Access = private)
        function loadSimpleTestFixture(this)
            m = this.process;

            %whole cell model IDs
            m.compartment.wholeCellModelIDs = {'c';'e'};
            m.compartment.cytosolIndexs = 1;
            m.compartment.extracellularIndexs = 2;
            m.stimuliWholeCellModelIDs = {
                'G6P';
                'PI';
                'stimulus_gluconate';
                'stimulus_ironStress';
                'stimulus_thiolStress';
                'temperature'};
            m.substrateWholeCellModelIDs = {
                'MG_101_MONOMER';
                'MG_127_MONOMER';
                'MG_236_MONOMER';
                'MG_085_HEXAMER';
                'MG_205_DIMER';
                'MG_409_DIMER'};

            %indices
            nCmp = length(m.compartment.wholeCellModelIDs);

            m.stimulusStimulusLocalIndexs          = (3:6)';
            m.stimulusMetaboliteLocalIndexs        = (1:2)';
            m.stimulusRNALocalIndexs               = zeros(0,1);
            m.stimulusMonomerLocalIndexs           = zeros(0,1);
            m.stimulusComplexLocalIndexs           = zeros(0,1);

            m.substrateStimulusLocalIndexs         = zeros(0,1);
            m.substrateMetaboliteLocalIndexs       = zeros(0,1);
            m.substrateRNALocalIndexs              = zeros(0,1);
            m.substrateMonomerLocalIndexs          = (1:3)';
            m.substrateComplexLocalIndexs          = (4:6)';

            m.enzymeStimulusLocalIndexs            = zeros(0,1);
            m.enzymeMetaboliteLocalIndexs          = zeros(0,1);
            m.enzymeRNALocalIndexs                 = zeros(0,1);
            m.enzymeMonomerLocalIndexs             = zeros(0,1);
            m.enzymeComplexLocalIndexs             = zeros(0,1);

            m.stimulusStimulusGlobalIndexs         = repmat((1:4)', 1, nCmp);
            m.stimulusMetaboliteGlobalIndexs       = repmat((1:2)', 1, nCmp);
            m.stimulusRNAGlobalIndexs              = zeros(0, nCmp);
            m.stimulusMonomerGlobalIndexs          = zeros(0, nCmp);
            m.stimulusComplexGlobalIndexs          = zeros(0, nCmp);

            m.substrateStimulusGlobalIndexs        = zeros(0, nCmp);
            m.substrateMetaboliteGlobalIndexs      = zeros(0, nCmp);
            m.substrateRNAGlobalIndexs             = zeros(0, nCmp);
            m.substrateMonomerGlobalIndexs         = repmat((1:3)', 1, nCmp);
            m.substrateComplexGlobalIndexs         = repmat((1:3)', 1, nCmp);

            m.enzymeStimulusGlobalIndexs           = zeros(0, nCmp);
            m.enzymeMetaboliteGlobalIndexs         = zeros(0, nCmp);
            m.enzymeRNAGlobalIndexs                = zeros(0, nCmp);
            m.enzymeMonomerGlobalIndexs            = zeros(0, nCmp);
            m.enzymeComplexGlobalIndexs            = zeros(0, nCmp);
            
            m.stimulusStimulusCompartmentIndexs    = repmat(1:nCmp, 6, 1);
            m.stimulusMetaboliteCompartmentIndexs  = repmat(1:nCmp, 6, 1);
            m.stimulusRNACompartmentIndexs         = repmat(1:nCmp, 6, 1);
            m.stimulusMonomerCompartmentIndexs     = repmat(1:nCmp, 6, 1);
            m.stimulusComplexCompartmentIndexs     = repmat(1:nCmp, 6, 1);
            
            m.substrateStimulusCompartmentIndexs   = repmat(1:nCmp, 6, 1);
            m.substrateMetaboliteCompartmentIndexs = repmat(1:nCmp, 6, 1);
            m.substrateRNACompartmentIndexs        = repmat(1:nCmp, 6, 1);
            m.substrateMonomerCompartmentIndexs    = repmat(1:nCmp, 6, 1);
            m.substrateComplexCompartmentIndexs    = repmat(1:nCmp, 6, 1);
            
            m.enzymeStimulusCompartmentIndexs      = zeros(0, nCmp);
            m.enzymeMetaboliteCompartmentIndexs    = zeros(0, nCmp);
            m.enzymeRNACompartmentIndexs           = zeros(0, nCmp);
            m.enzymeMonomerCompartmentIndexs       = zeros(0, nCmp);
            m.enzymeComplexCompartmentIndexs       = zeros(0, nCmp);
            
            %activation rules
            m.activationRules = {
                '~stimulus_gluconate';
                'stimulus_thiolStress';
                'stimulus_ironStress';
                'G6P>5';
                'temperature>=43';
                'PI>20'};

            %cell volume (L)
            m.geometry.volume = 2e-021;

            %initial state
            m.stimuli               = zeros(length(m.stimuliWholeCellModelIDs),   nCmp);
            m.substrates            = zeros(length(m.substrateWholeCellModelIDs), nCmp);
            m.inactivatedSubstrates = zeros(length(m.substrateWholeCellModelIDs), nCmp);
        end
    end
end
