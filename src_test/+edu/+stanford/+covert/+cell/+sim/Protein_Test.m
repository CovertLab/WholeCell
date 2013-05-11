%Protein medium test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 6/6/2011
classdef Protein_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = Protein_Test(name)
            this = this@TestCase(name);
        end
        
        function setUp(this)
            this.simulation = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'Metabolism'
                'tRNAAminoacylation'
                ...
                'RibosomeAssembly'
                'Translation'
                'ProteinProcessingI'
                'ProteinTranslocation'
                'ProteinProcessingII'
                'ProteinFolding'
                'ProteinModification'
                'MacromolecularComplexation'
                'TerminalOrganelleAssembly'
                'ProteinActivation'
                'ProteinDecay'});
        end
        
        function testSynthesisPathway(this)
            s = this.simulation;
            m = s.state('Metabolite');
            rna = s.state('Rna');
            pm = s.state('ProteinMonomer');
            pc = s.state('ProteinComplex');
            stim = s.state('Stimulus');
            rib = s.process('RibosomeAssembly');
            trl = s.process('Translation');
            pr1 = s.process('ProteinProcessingI');
            tlc = s.process('ProteinTranslocation');
            pr2 = s.process('ProteinProcessingII');
            fld = s.process('ProteinFolding');
            mod = s.process('ProteinModification');
            %cpx = s.process('MacromolecularComplexation');
            org = s.process('TerminalOrganelleAssembly');
            act = s.process('ProteinActivation');
            dcy = s.process('ProteinDecay');
            
            %% initialize protein assembly line
            rna.counts(rna.matureIndexs(rna.matureMRNAIndexs)) = 1;
            rna.counts(rna.matureIndexs(rna.matureRRNAIndexs(rib.rnaGlobalIndexs))) = 500;  %ribosome parts
            rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) = 1e6;  %aminoacylated tRNAs
            
            pm.counts(:) = 0;
            pm.counts(pm.matureIndexs(rib.monomerGlobalIndexs)) = 500;  %ribosome parts
            pm.counts(pm.matureIndexs(pm.getIndexs(rib.enzymeWholeCellModelIDs))) = 1;  %ribosome assembly enzymes
            pm.counts(pm.matureIndexs(trl.enzymeToMonomer([
                trl.enzymeIndexs_translationFactors]))) = 500;
            pm.counts(pm.matureIndexs(pr1.enzymeToMonomer([
                pr1.enzymeIndexs_methionineAminoPeptidase]))) = 200;
            pm.counts(pm.matureIndexs(tlc.enzymeToMonomer([
                tlc.enzymeIndexs_signalRecognitionParticleReceptor]))) = 500;
            pm.counts(pm.matureIndexs(pr2.enzymeToMonomer([
                pr2.enzymeIndexs_signalPeptidase;
                pr2.enzymeIndexs_diacylglycerylTransferase])),...
                s.compartment.membraneIndexs) = 500;
            pm.counts(pm.matureIndexs(fld.enzymeToMonomer([
                fld.enzymeIndexs_triggerFactor;
                fld.enzymeIndexs_dnaK]))) = 1;
            pm.counts(pm.matureIndexs(mod.enzymeToMonomer([
                mod.enzymeIndexs_lipoylTransferase;
                mod.enzymeIndexs_glutamateLigase]))) = 6000;
            pm.counts(pm.matureIndexs(org.enzymeToMonomer([
                org.enzymeIndexs_HMW1;
                org.enzymeIndexs_HMW2;
                org.enzymeIndexs_HMW3])),...
                s.compartment.terminalOrganelleCytosolIndexs) = 1;
            pm.counts(pm.matureIndexs(org.enzymeToMonomer([
                org.enzymeIndexs_P32])),...
                s.compartment.terminalOrganelleMembraneIndexs) = 1;
            
            pc.counts(:) = 0;
            pc.counts(pc.matureIndexs(trl.enzymeToComplex([
                trl.enzymeIndexs_translationFactors]))) = 500;
            pc.counts(pc.matureIndexs(pr1.enzymeToComplex([
                pr1.enzymeIndexs_deformylase]))) = 200;
            pc.counts(pc.matureIndexs(tlc.enzymeToComplex([
                tlc.enzymeIndexs_signalRecognitionParticle;
                tlc.enzymeIndexs_translocaseATPase]))) = 500;
            pc.counts(pc.matureIndexs(tlc.enzymeToComplex([
                tlc.enzymeIndexs_translocasePore])),...
                s.compartment.membraneIndexs) = 500;
            pc.counts(pc.matureIndexs(fld.enzymeToComplex([
                fld.enzymeIndexs_dnaJ;
                fld.enzymeIndexs_grpE;
                fld.enzymeIndexs_groELES]))) = 1;
            pc.counts(pc.matureIndexs(mod.enzymeToComplex([
                mod.enzymeIndexs_serineThreonineKinase])),...
                s.compartment.membraneIndexs) = 900;
            
            metCounts = m.counts;
            m.counts(:) = 1e9;
            m.counts(act.stimulusMetaboliteGlobalCompartmentIndexs) = ...
                metCounts(act.stimulusMetaboliteGlobalCompartmentIndexs);
            
            trl.ribosomeElongationRate = 200;
            
            %% assemble ribosome subunits
            ribosomeSubunitIdxs = pc.getIndexs({'RIBOSOME_30S' 'RIBOSOME_50S' 'RIBOSOME_30S_IF3'});
            s.evolveState();
            s.evolveState();
            assertEqual([500; 500], ...
                + pc.counts(pc.matureIndexs(ribosomeSubunitIdxs(1:2))) ...
                + [pc.counts(pc.matureIndexs(ribosomeSubunitIdxs(3))); 0] ...
                + pc.counts(pc.matureIndexs(pc.ribosome70SIndexs)) ...
                + pc.counts(pc.boundIndexs(pc.ribosome70SIndexs)) ...
                )
            
            %% translate, process, translocate, fold, modify, complex
            initMonomers = nnz(pm.counts(pm.matureIndexs,:)+pm.counts(pm.boundIndexs,:));
            initComplexs = nnz(pc.counts(pc.matureIndexs,:)+pc.counts(pc.boundIndexs,:));
            proteinSynthesized = false;
            for i = 1:20
                s.evolveState();
                
                if ...
                        nnz(pm.counts(pm.signalSequenceIndexs, :)) > 0 && ...
                        nnz(pm.counts(pm.matureIndexs, :)) > initMonomers && ...
                        nnz(pc.counts(pc.matureIndexs, :) + pc.counts(pc.boundIndexs, :)) > initComplexs
                    proteinSynthesized = true;
                end
            end
            assertTrue(proteinSynthesized);
            
            %% fold any new complexs (complexation runs after folding)
            fld.copyFromState();
            fld.evolveState();
            fld.copyToState();
            assertEqual(0, nnz(pc.counts(pc.nascentIndexs, :)));
            
            %% terminal organelle assembly
            assertAllEqual(0, nnz(pm.counts(...
                pm.matureIndexs(org.substrateMonomerGlobalIndexs), [
                s.compartment.cytosolIndexs;
                s.compartment.chromosomeIndexs;
                s.compartment.extracellularIndexs;
                s.compartment.membraneIndexs])));
            assertTrue(4 <= ...
                nnz(pm.counts(...
                pm.matureIndexs(org.substrateMonomerGlobalIndexs),...
                s.compartment.terminalOrganelleCytosolIndexs)));
            assertTrue(1 <= ...
                nnz(pm.counts(...
                pm.matureIndexs(org.substrateMonomerGlobalIndexs),...
                s.compartment.terminalOrganelleMembraneIndexs)));
            
            %% deactivation and reactivation
            [~,actMonomerIdxs] = ismember({
                'MG_101_MONOMER';
                'MG_127_MONOMER';
                'MG_236_MONOMER'}, pm.wholeCellModelIDs(pm.matureIndexs));
            [~,actComplexIdxs] = ismember({
                'MG_085_HEXAMER';
                'MG_205_DIMER';
                'MG_409_DIMER'}, pc.wholeCellModelIDs(pc.matureIndexs));
            pm.counts(pm.matureIndexs(actMonomerIdxs)) = 1;
            pc.counts(pc.matureIndexs(actComplexIdxs)) = 1;
            pm.counts(pm.inactivatedIndexs(actMonomerIdxs)) = 0;
            pc.counts(pc.inactivatedIndexs(actComplexIdxs)) = 0;
            
            j = edu.stanford.covert.cell.sim.constant.Condition.valueIndexs;
            stim.setValues(strcmp(stim.wholeCellModelIDs, 'temperature'), j) = 40;
            stim.setValues(strcmp(stim.wholeCellModelIDs, 'stimulus_gluconate'), j) = false;
            stim.setValues(strcmp(stim.wholeCellModelIDs, 'stimulus_thiolStress'), j) = false;
            stim.setValues(strcmp(stim.wholeCellModelIDs, 'stimulus_ironStress'), j) = false;
            m.counts(strcmp(m.wholeCellModelIDs, 'G6P'), :) = 1e12;
            m.counts(strcmp(m.wholeCellModelIDs, 'PI'), :) = 1e12;
            
            s.evolveState();
            assertEqual([false true  true ], any(pm.counts(pm.inactivatedIndexs(actMonomerIdxs), :), 2)');
            assertEqual([true  false false], any(pm.counts(pm.matureIndexs(actMonomerIdxs), :), 2)');
            assertEqual([false true  false], any(pc.counts(pc.inactivatedIndexs(actComplexIdxs), :), 2)');
            assertEqual([true  false true ], any(pc.counts(pc.matureIndexs(actComplexIdxs), :), 2)');
            
            %change conditions to flip the activation state of all the proteins
            j = edu.stanford.covert.cell.sim.constant.Condition.valueIndexs;
            stim.setValues(strcmp(stim.wholeCellModelIDs, 'temperature'), j) = 50;
            stim.setValues(strcmp(stim.wholeCellModelIDs, 'stimulus_gluconate'), j) = true;
            stim.setValues(strcmp(stim.wholeCellModelIDs, 'stimulus_thiolStress'), j) = true;
            stim.setValues(strcmp(stim.wholeCellModelIDs, 'stimulus_ironStress'), j) = true;
            m.counts(:) = 0;
            
            s.evolveState();
            assertEqual([true  false false], any(pm.counts(pm.inactivatedIndexs(actMonomerIdxs)), 2)');
            assertEqual([false true  true ], any(pm.counts(pm.matureIndexs(actMonomerIdxs)), 2)');
            assertEqual([true  false true ], any(pc.counts(pc.inactivatedIndexs(actComplexIdxs)), 2)');
            assertEqual([false true  false], any(pc.counts(pc.matureIndexs(actComplexIdxs)), 2)');            
            
            %% decay
            pm.counts(pm.matureIndexs(dcy.enzymeToMonomer(1:length(dcy.enzymes)))) = 1;
            pc.counts(pc.matureIndexs(dcy.enzymeToComplex(1:length(dcy.enzymes)))) = 1;
            m.counts(:) = 1e9;
            m.counts(strcmp(m.wholeCellModelIDs, 'PI'), :) = 0;
            m.counts(dcy.substrateToMetabolite(dcy.substrateIndexs_aminoAcids)) = 0;
            dcy.lonProteaseSpecificRate = inf;
            dcy.ftsHProteaseSpecificRate = inf;
            dcy.oligoendopeptidaseFSpecificRate = inf;
            pm.decayRates(:) = inf;
            pc.decayRates(:) = inf;
            
            s.evolveState();
            assertEqual([0 0], sum(pc.counts(:,[
                s.compartment.cytosolIndexs;
                s.compartment.terminalOrganelleCytosolIndexs])));
            assertEqual([0 0], sum(pm.counts(setdiff(1:end, pm.signalSequenceIndexs), [
                s.compartment.cytosolIndexs;
                s.compartment.terminalOrganelleCytosolIndexs])));
        end
        
        function testSynthesisMetabolismInterface1MRNAAtATime(this)
            s = this.simulation;
            s.applyOptions('verbosity', 0);
            
            comp = s.compartment;
            g = s.gene;
            r = s.state('Rna');
            m = s.state('Metabolite');
            pm = s.state('ProteinMonomer');
            pc = s.state('ProteinComplex');
            met = s.process('Metabolism');
            trl = s.process('Translation');
            pmod = s.process('ProteinModification');
            
            pcComp = sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3);
            
            pm.decayRates(:) = 0;
            pc.decayRates(:) = 0;
            
            pm.counts(setdiff(1:end, [pm.matureIndexs; pm.boundIndexs; pm.damagedIndexs; pm.misfoldedIndexs; pm.inactivatedIndexs]), :) = 0;
            pc.counts(setdiff(1:end, [pc.matureIndexs; pc.boundIndexs; pc.damagedIndexs; pc.misfoldedIndexs; pc.inactivatedIndexs]), :) = 0;
            
            pmod.enzymes(pmod.enzymeMonomerLocalIndexs) = 1e3;
            pmod.modifiedMonomers(pmod.enzymeMonomerGlobalIndexs) = 1e3;
            pmod.copyToState();
            
            met.substrates = ...
                met.substrates + ...
                max(0, met.randStream.stochasticRound(met.metabolismProduction));
            met.copyToState();
            
            initMetabolites = m.counts;
            initMonomers = pm.counts;
            initComplexs = pc.counts;
            
            initTotMonomers = sum(...
                + initMonomers(pm.matureIndexs, :)...
                + initMonomers(pm.boundIndexs, :) ...
                + initMonomers(pm.misfoldedIndexs, :) ...
                + initMonomers(pm.damagedIndexs, :) ...
                + initMonomers(pm.inactivatedIndexs, :) + pcComp * (...
                + initComplexs(pc.nascentIndexs, :) ...
                + initComplexs(pc.matureIndexs, :) ...
                + initComplexs(pc.boundIndexs, :) ...
                + initComplexs(pc.misfoldedIndexs, :) ...
                + initComplexs(pc.damagedIndexs, :) ...
                + initComplexs(pc.inactivatedIndexs, :)), 2);
            
            for i = 1:numel(r.matureMRNAIndexs)
                if s.verbosity >= 1
                    fprintf('RNA %d\n', i);
                end
                
                iMonomer = find(r.matureRNAGeneComposition(g.mRNAIndexs, r.matureMRNAIndexs(i)));
                
                m.counts = initMetabolites;
                pm.counts = initMonomers;
                pc.counts = initComplexs;
                
                r.counts(r.matureIndexs(r.matureMRNAIndexs), comp.cytosolIndexs) = 0;
                r.counts(r.matureIndexs(r.matureMRNAIndexs(i)), comp.cytosolIndexs) = 1;
                pm.counts(pm.nascentIndexs(iMonomer), comp.cytosolIndexs) = 1;
                
                trl.copyFromState();
                trl.initializeState();
                trl.copyToState();
               
                for j = 1:200
                    if mod(j, 10) == 1 && s.verbosity >= 2
                        fprintf('\ttime %d\n', j);
                    end
                    s.evolveState();
                    
                    totMonomers = sum(...
                        + pm.counts(pm.matureIndexs(iMonomer), :)...
                        + pm.counts(pm.boundIndexs(iMonomer), :) ...
                        + pm.counts(pm.misfoldedIndexs(iMonomer), :) ...
                        + pm.counts(pm.damagedIndexs(iMonomer), :) ...
                        + pm.counts(pm.inactivatedIndexs(iMonomer), :) + pcComp(iMonomer, :) * (...
                        + pc.counts(pc.nascentIndexs, :) ...
                        + pc.counts(pc.matureIndexs, :)...
                        + pc.counts(pc.boundIndexs, :) ...
                        + pc.counts(pc.misfoldedIndexs, :) ...
                        + pc.counts(pc.damagedIndexs, :) ...
                        + pc.counts(pc.inactivatedIndexs, :)), 2);
                    
                    monomerFlag = all(totMonomers - initTotMonomers(iMonomer) > 0);
                    
                    if monomerFlag
                        break;
                    end
                end
                
                assertTrue(monomerFlag, sprintf('Couldn''t produce RNA: %d, monomers: %s', i, strjoin(', ', pm.wholeCellModelIDs{pm.matureIndexs(iMonomer)})));
            end
        end
        
        function testComplexSynthesisMetabolismInterface(this)
            s = this.simulation;
            s.applyOptions('verbosity', 0);
            
            comp = s.compartment;
            r = s.state('Rna');
            pm = s.state('ProteinMonomer');
            pc = s.state('ProteinComplex');
            met = s.process('Metabolism');
            mc = s.process('MacromolecularComplexation');
            
            pm.decayRates(:) = 0;
            pc.decayRates(:) = 0;
            
            r.counts(r.matureIndexs(r.matureSRNAIndexs), comp.cytosolIndexs) = ...
                r.counts(r.matureIndexs(r.matureSRNAIndexs), comp.cytosolIndexs) + 200;
            pm.counts(sub2ind(size(pm.counts), pm.matureIndexs, pm.compartments(pm.matureIndexs))) = ...
                pm.counts(sub2ind(size(pm.counts), pm.matureIndexs, pm.compartments(pm.matureIndexs))) + 200;
            
            met.substrates = ...
                met.substrates + ...
                max(0, met.randStream.stochasticRound(met.metabolismProduction));
            met.copyToState();
            
            assertAllEqual(0, pc.counts(pc.nascentIndexs, :));
            
            initTotComplexs = sum(...
                + pc.counts(pc.matureIndexs, :) ...
                + pc.counts(pc.boundIndexs, :) ...
                + pc.counts(pc.misfoldedIndexs, :) ...
                + pc.counts(pc.damagedIndexs, :) ...
                + pc.counts(pc.inactivatedIndexs, :), 2);
            
            iComplex = mc.complexGlobalIndexs;
            
            for j = 1:50
                if mod(j, 10) == 1 && s.verbosity >= 2
                    fprintf('time %d\n', j);
                end
                s.evolveState();
                
                totComplexs = sum(...
                    + pc.counts(pc.matureIndexs, :) ...
                    + pc.counts(pc.boundIndexs, :) ...
                    + pc.counts(pc.misfoldedIndexs, :) ...
                    + pc.counts(pc.damagedIndexs, :) ...
                    + pc.counts(pc.inactivatedIndexs, :), 2);
                
                complexFlag = all(totComplexs(iComplex) - initTotComplexs(iComplex) > 0);
                
                if complexFlag
                    break;
                end
            end
            
            assertTrue(complexFlag, sprintf('Couldn''t make complexs %s', strjoin(', ', pc.wholeCellModelIDs{iComplex(totComplexs(iComplex) - initTotComplexs(iComplex) == 0)})));
        end
        
        function testPostTranslationSynthesisMetabolismInterface(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            s = this.simulation;
            s.applyOptions('verbosity', 0);
            
            comp = s.compartment;
            g = s.gene;
            mass = s.state('Mass');
            r = s.state('Rna');
            pm = s.state('ProteinMonomer');
            pc = s.state('ProteinComplex');
            met = s.process('Metabolism');
            mc = s.process('MacromolecularComplexation');
            trl = s.process('Translation');
            
            pm.decayRates(:) = 0;
            pc.decayRates(:) = 0;
            
            trl.ribosomeElongationRate = 0;
            
            rnaExp = r.expression(r.matureIndexs);
            mRNAExp = r.matureRNAGeneComposition(g.mRNAIndexs, :) * r.expression(r.matureIndexs);
            initTotRNAs = r.randStream.stochasticRound(...
                rnaExp * mass.cellInitialDryWeight * mass.dryWeightFractionRNA * mass.initialFractionNTPsInRNAs * ...
                ConstantUtil.nAvogadro / (rnaExp' * r.molecularWeights(r.matureIndexs)));
            initTotMonomers = pm.randStream.stochasticRound(mRNAExp * mass.cellInitialDryWeight * mass.dryWeightFractionProtein * mass.initialFractionAAsInMonomers * ...
                ConstantUtil.nAvogadro / (mRNAExp' * pm.molecularWeights(pm.matureIndexs)));
            
            initTotComplexs = sum(...
                + pc.counts(pc.matureIndexs, :) ...
                + pc.counts(pc.boundIndexs, :) ...
                + pc.counts(pc.inactivatedIndexs, :), 2);
            
            r.counts(r.matureIndexs, comp.cytosolIndexs) = ...
                r.counts(r.matureIndexs, comp.cytosolIndexs) + ...
                initTotRNAs;
            pm.counts(pm.nascentIndexs, comp.cytosolIndexs) = ...
                pm.counts(pm.nascentIndexs, comp.cytosolIndexs) + ...
                initTotMonomers;
            
            met.substrates = ...
                met.substrates + ...
                max(0, met.randStream.stochasticRound(met.metabolismProduction));
            met.copyToState();
            
            iComplex = mc.complexGlobalIndexs;
            
            s.process('ProteinTranslocation').copyFromState();
            s.process('ProteinTranslocation').translocaseSpecificRate = 1e6;
            s.process('ProteinTranslocation').enzymes([
                s.process('ProteinTranslocation').enzymeIndexs_signalRecognitionParticle
                s.process('ProteinTranslocation').enzymeIndexs_signalRecognitionParticleReceptor
                ]) = 1e6;
            s.process('ProteinTranslocation').copyToState();
            s.process('ProteinProcessingII').lipoproteinSignalPeptidaseSpecificRate = 1e6;
            s.process('ProteinProcessingII').lipoproteinDiacylglycerylTransferaseSpecificRate = 1e6;
            s.process('ProteinModification').enzymeBounds(:, 2) = 1e6;
            s.process('ProteinModification').initializeSpeciesNetwork();
            
            for j = 1:200
                if mod(j, 10) == 1 && s.verbosity >= 2
                    fprintf('time %d\n', j);
                end
                s.evolveState();
                if ...
                        ~nnz(pm.counts([pm.nascentIndexs; pm.processedIIndexs; pm.processedIIIndexs; pm.foldedIndexs], :)) && ...
                        ~nnz(pc.counts(pc.nascentIndexs, :))
                    break;
                end
            end
            
            assertAllEqual(0, pm.counts(pm.nascentIndexs, :));
            assertAllEqual(0, pm.counts(pm.processedIIndexs, :));
            assertAllEqual(0, pm.counts(pm.processedIIIndexs, :));
            assertAllEqual(0, pm.counts(pm.foldedIndexs, :));
            assertAllEqual(0, pc.counts(pc.nascentIndexs, :));
            
            totComplexs = sum(...
                + pc.counts(pc.matureIndexs, :) ...
                + pc.counts(pc.boundIndexs, :) ...
                + pc.counts(pc.inactivatedIndexs, :), 2);
            
            assertElementsAlmostEqual(...
                2 * initTotComplexs(setdiff(iComplex, pc.getIndexs({'MG_0001_048'; 'RNA_POLYMERASE'; 'MG_213_214_298_6MER'}))), ...
                totComplexs(setdiff(iComplex, pc.getIndexs({'MG_0001_048'; 'RNA_POLYMERASE'; 'MG_213_214_298_6MER'}))), ...
                'relative', 0.5, 10);
            assertEqual(...
                initTotComplexs(setdiff(1:end, [iComplex; pc.getIndexs({'RNA_POLYMERASE_HOLOENZYME'; 'RIBOSOME_30S'; 'RIBOSOME_30S_IF3'; 'RIBOSOME_50S'; 'RIBOSOME_70S'; 'MG_213_214_298_6MER_ATP'})])), ...
                totComplexs(setdiff(1:end, [iComplex; pc.getIndexs({'RNA_POLYMERASE_HOLOENZYME'; 'RIBOSOME_30S'; 'RIBOSOME_30S_IF3'; 'RIBOSOME_50S'; 'RIBOSOME_70S'; 'MG_213_214_298_6MER_ATP'})])) );
        end
        
        function testTranslationExpectation(this)
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            
            s = this.simulation;
            m = s.process('Translation');
            g = m.gene;
            met = m.metabolite;
            r = m.rna;
            pm = m.monomer;
            rib = m.ribosome;
            pol = m.polypeptide;
            
            iterMax = 5000;
            
            monomerTRNACounts = zeros(size(pol.monomerTRNASequences, 1), 36);
            for i = 1:numel(pol.monomerTRNASequences)
                monomerTRNACounts(i, :) = histc(pol.monomerTRNASequences{i}, 1:36);
            end
            
            m.substrates(:) = 1e7;
            m.aminoacylatedTRNAs(:) = 1e6;
            m.monomers(:) = 0;
            m.enzymes = m.enzymes + m.boundEnzymes;
            m.boundEnzymes(:) = 0;
            m.mRNAs = floor(3 * m.randStream.rand(size(m.mRNAs)));
            
            numRibs = min([
                m.enzymes(m.enzymeIndexs_ribosome30S)
                m.enzymes(m.enzymeIndexs_ribosome50S)
                m.enzymes(m.enzymeIndexs_elongationFactors)]);
            
            rib.boundMRNAs = zeros(numRibs, 1 );
            rib.mRNAPositions = zeros(numRibs, 1 );
            rib.tmRNAPositions = zeros(numRibs, 1 );
            rib.states = repmat(rib.notExistValue, numRibs, 1 );
            
            pol.boundMRNAs = zeros(numRibs, 1 );
            pol.nascentMonomerLengths = zeros(numRibs, 1 );
            pol.proteolysisTagLengths = zeros(numRibs, 1 );
            pol.abortedPolypeptides = zeros(0, 3);
            
            initMonomers = m.monomers;
            initBoundMRNAs = pol.boundMRNAs;
            initNascentMonomerLengths = pol.nascentMonomerLengths;
            initProteolysisTagLengths = pol.proteolysisTagLengths;
            for i = 1:iterMax
                m.evolveState();
            end
            
            %all ribosomes in active/free state
            assertFalse(any(rib.states == rib.stalledValue));
            
            %amino acids used and accounted
            usedAAs = pm.baseCounts(pm.nascentIndexs, met.aminoAcidIndexs)' * (m.monomers - initMonomers);
            usedTRNAs = monomerTRNACounts' * (m.monomers - initMonomers);
            
            for i = 1:size(initBoundMRNAs, 1)
                sequence = [];
                tRNASequence = [];
                if initNascentMonomerLengths(i) > 0
                    sequence = [sequence pol.monomerAASequences{initBoundMRNAs(i)}(1:initNascentMonomerLengths(i))]; %#ok<*AGROW>
                    tRNASequence = [tRNASequence pol.monomerTRNASequences{initBoundMRNAs(i)}(1:initNascentMonomerLengths(i))];
                end
                if initProteolysisTagLengths(i) > 0
                    sequence = [sequence pol.proteolysisTagAASequence(1:initProteolysisTagLengths(i))]; %#ok<*AGROW>
                    tRNASequence = [tRNASequence pol.proteolysisTagTRNASequence(1:initProteolysisTagLengths(i))]; %#ok<*AGROW>
                end
                if ~isempty(sequence)
                    usedAAs = usedAAs - ProteinMonomer.computeBaseCount(sequence, 21, 1:21, true)';
                    usedTRNAs = usedTRNAs - histc(tRNASequence, 1:36);
                end
            end
            
            for i = 1:size(pol.boundMRNAs, 1)
                sequence = [];
                tRNASequence = [];
                if pol.nascentMonomerLengths(i) > 0
                    sequence = [sequence pol.monomerAASequences{pol.boundMRNAs(i)}(1:pol.nascentMonomerLengths(i))]; %#ok<*AGROW>
                    tRNASequence = [tRNASequence pol.monomerTRNASequences{pol.boundMRNAs(i)}(1:pol.nascentMonomerLengths(i))]; %#ok<*AGROW>
                end
                if pol.proteolysisTagLengths(i) > 0
                    sequence = [sequence pol.proteolysisTagAASequence(1:pol.proteolysisTagLengths(i))]; %#ok<*AGROW>
                    tRNASequence = [tRNASequence pol.proteolysisTagTRNASequence(1:pol.proteolysisTagLengths(i))]; %#ok<*AGROW>
                end
                if ~isempty(sequence)
                    usedAAs = usedAAs + ProteinMonomer.computeBaseCount(sequence, 21, 1:21, true)';
                    usedTRNAs = usedTRNAs + histc(tRNASequence, 1:36);
                end
            end
            
            assertEqual(1e6 - m.aminoacylatedTRNAs, usedTRNAs);
            assertEqual(sum(usedTRNAs), sum(usedAAs));
            assertElementsAlmostEqual(sum(usedAAs), numRibs * m.ribosomeElongationRate * iterMax, 'relative', 1e-1);
            aaProd = pm.baseCounts(pm.nascentIndexs, met.aminoAcidIndexs)' * r.matureRNAGeneComposition(g.mRNAIndexs, :) * r.expression(r.matureIndexs);
            trnaUsage = monomerTRNACounts' * r.matureRNAGeneComposition(g.mRNAIndexs, :) * r.expression(r.matureIndexs);
            assertIn(corr(usedAAs, aaProd), [0.99 1]);
            assertIn(180 / pi * acos(usedAAs' * aaProd / (sqrt(usedAAs' * usedAAs) * sqrt(aaProd' * aaProd))), [0 5]);
            assertIn(corr(usedTRNAs, trnaUsage), [0.99 1]);
            assertIn(180 / pi * acos(usedTRNAs' * trnaUsage / (sqrt(usedTRNAs' * usedTRNAs) * sqrt(trnaUsage' * trnaUsage))), [0 5]);
            
            %monomer produced according to expectations
            assertIn(corr(m.monomers, m.mRNAs), [0.60 1]);
            assertIn(180 / pi * acos((m.monomers' * m.mRNAs) / (sqrt(m.monomers' * m.monomers) * sqrt(m.mRNAs' * m.mRNAs))), ...
                [0 40]);
        end
    end
    
    %amino acid production, aminoacylation, polymerization
    methods
        function testAminoacylation1(~)
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            import edu.stanford.covert.util.ConstantUtil;
            
            s = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'tRNAAminoacylation'
                });
            
            %% references
            comp = s.compartment;
            
            time = s.state('Time');
            m = s.state('Metabolite');
            rna = s.state('Rna');
            mr = s.state('MetabolicReaction');
            
            ta = s.process('tRNAAminoacylation');
            
            %% constants
            [~, idxs1] = ismember(m.aminoAcidIndexs(1:20), ta.substrateGlobalIndexs);
            [~, idxs2] = ismember(rna.matureTRNAIndexs, ta.aminoacylatedRNAGlobalIndexs);
            aaToTRNA = ...
                -ta.reactionStoichiometryMatrix(idxs1, :) * ...
                ta.reactionModificationMatrix(:, idxs2);
            
            %% initialize
            rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = ...
                + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs);
            rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = 0;
            
            m.counts(m.ntpIndexs([1 3]), comp.cytosolIndexs) = 1e12;
            m.counts(m.waterIndexs, comp.cytosolIndexs) = 1e12;
            m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = 0;
            
            %% simulation
            initAAs = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initRNAs = rna.counts;
            
            iterMax = 500;
            
            aas = zeros(21, iterMax);
            aaIncorporated = zeros(21, iterMax);
            aaProd = zeros(21, iterMax);
            for i = 1:iterMax
                %fake metabolism
                mr.growth = log(2) / time.cellCycleLength * exp(log(2) * i / time.cellCycleLength);
                newAAs = max(0, m.randStream.stochasticRound(mr.growth * (m.biomassProduction(m.aminoAcidIndexs, comp.cytosolIndexs) - m.byproducts(m.aminoAcidIndexs, comp.cytosolIndexs))));
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + newAAs;
                aaProd(:, i) = newAAs;
                
                %fake translation
                f = rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs));
                rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) = ...
                    + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) ...
                    -  f;
                rna.counts(rna.matureIndexs(rna.matureTRNAIndexs)) = ...
                    + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs)) ...
                    + f;
                aaIncorporated(1:20, i) = aaToTRNA * f;
                
                %simulate
                s.evolveState();
                
                %keep track
                aas(:, i) = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            end
            
            %% assertions
            % account for all amino acids
            initCounts = ...
                + initAAs ...
                + rna.baseCounts(:, m.aminoAcidIndexs)' * sum(initRNAs - rna.counts, 2) ...
                + sum(aaProd - aaIncorporated, 2);
            currCounts = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initCounts(13) = initCounts(13) + initCounts(21);
            currCounts(13) = currCounts(13) + currCounts(21);
            assertEqual(initCounts(1:20)', currCounts(1:20)');
            
            %amino acids are limiting
            assertIn(min(m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)), [0 1]);
        end
        
        function testAminoacylation2(~)
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            import edu.stanford.covert.util.ConstantUtil;
            
            s = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'tRNAAminoacylation'
                });
            
            %% references
            comp = s.compartment;
            g = s.gene;
            
            time = s.state('Time');
            m = s.state('Metabolite');
            rna = s.state('Rna');
            mr = s.state('MetabolicReaction');
            
            ta = s.process('tRNAAminoacylation');
            
            %% constants
            [~, idxs1] = ismember(m.aminoAcidIndexs(1:20), ta.substrateGlobalIndexs);
            [~, idxs2] = ismember(rna.matureTRNAIndexs, ta.aminoacylatedRNAGlobalIndexs);
            aaToTRNA = zeros(21, 36);
            aaToTRNA(1:20, :) = ...
                -ta.reactionStoichiometryMatrix(idxs1, :) * ...
                ta.reactionModificationMatrix(:, idxs2);
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            mRNAExp = mRNAExp / sum(mRNAExp);
            tRNAUsage = ta.monomerTRNACounts' * mRNAExp;
            tRNAUsage = tRNAUsage / sum(tRNAUsage);
            
            relAAProd = aaToTRNA * tRNAUsage;
            
            %% initialize
            rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = ...
                + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs);
            rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = 0;
            
            m.counts(m.ntpIndexs([1 3]), comp.cytosolIndexs) = 1e12;
            m.counts(m.waterIndexs, comp.cytosolIndexs) = 1e12;
            m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = 0;
            
            %% simulation
            initAAs = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initRNAs = rna.counts;
            
            iterMax = 500;
            
            aas = zeros(21, iterMax);
            aaIncorporated = zeros(21, iterMax);
            aaProd = zeros(21, iterMax);
            for i = 1:iterMax
                %fake metabolism
                mr.growth = log(2) / time.cellCycleLength * exp(log(2) * i / time.cellCycleLength);
                newAAs = m.randStream.stochasticRound(mr.growth * relAAProd);
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + newAAs;
                aaProd(:, i) = newAAs;
                
                %fake translation
                f = rna.randStream.stochasticRound(min(rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) ./ tRNAUsage) * tRNAUsage);
                rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) = ...
                    + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) ...
                    -  f;
                rna.counts(rna.matureIndexs(rna.matureTRNAIndexs)) = ...
                    + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs)) ...
                    + f;
                aaIncorporated(:, i) = aaToTRNA * f;
                
                %simulate
                s.evolveState();
                
                %keep track
                aas(:, i) = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            end
            
            %% assertions
            % account for all amino acids
            initCounts = ...
                + initAAs ...
                + rna.baseCounts(:, m.aminoAcidIndexs)' * sum(initRNAs - rna.counts, 2) ...
                + sum(aaProd - aaIncorporated, 2);
            currCounts = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initCounts(13) = initCounts(13) + initCounts(21);
            currCounts(13) = currCounts(13) + currCounts(21);
            assertEqual(initCounts(1:20)', currCounts(1:20)');
            
            %amino acids are limiting
            assertAllEqual(0, m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)');
        end
        
        function testAminoacylation3(~)
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            import edu.stanford.covert.util.ConstantUtil;
            
            s = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'tRNAAminoacylation'
                });
            
            %% references
            comp = s.compartment;
            g = s.gene;
            
            time = s.state('Time');
            mass = s.state('Mass');
            m = s.state('Metabolite');
            rna = s.state('Rna');
            mr = s.state('MetabolicReaction');
            
            ta = s.process('tRNAAminoacylation');
            
            %% constants
            [~, idxs1] = ismember(m.aminoAcidIndexs(1:20), ta.substrateGlobalIndexs);
            [~, idxs2] = ismember(rna.matureTRNAIndexs, ta.aminoacylatedRNAGlobalIndexs);
            aaToTRNA = zeros(21, 36);
            aaToTRNA(1:20, :) = ...
                -ta.reactionStoichiometryMatrix(idxs1, :) * ...
                ta.reactionModificationMatrix(:, idxs2);
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            mRNAExp = mRNAExp / sum(mRNAExp);
            tRNAUsage = ta.monomerTRNACounts' * mRNAExp;
            tRNAUsage = tRNAUsage / sum(tRNAUsage);
            
            relAAProd = aaToTRNA * tRNAUsage;
            
            initRnas = rna.expression * mass.cellInitialDryWeight * mass.dryWeightFractionRNA / ...
                (rna.expression' * rna.molecularWeights / edu.stanford.covert.util.ConstantUtil.nAvogadro);
            
            %% initialize
            rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = ...
                + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs);
            rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = 0;
            
            m.counts(m.ntpIndexs([1 3]), comp.cytosolIndexs) = 1e12;
            m.counts(m.waterIndexs, comp.cytosolIndexs) = 1e12;
            m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = 0;
            
            %% simulation
            initAAs = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initRNAs = rna.counts;
            
            iterMax = 500;
            
            aas = zeros(21, iterMax);
            aaIncorporated = zeros(21, iterMax);
            aaProd = zeros(21, iterMax);
            for i = 1:iterMax
                %fake RNA synthesis
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    + rna.randStream.stochasticRound((initRnas .* (log(2) / time.cellCycleLength + rna.decayRates)) * exp(log(2) * i / time.cellCycleLength));
                
                %fake RNA decay
                decayedRNAs = rna.randStream.stochasticRound(rna.decayRates .* rna.counts(:, comp.cytosolIndexs));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    - decayedRNAs;
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + rna.baseCounts(:, m.aminoAcidIndexs)' * decayedRNAs;
                m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) ...
                    + m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs);
                m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs) = 0;
                
                %fake metabolism
                mr.growth = log(2) / time.cellCycleLength * exp(log(2) * i / time.cellCycleLength);
                newAAs = m.randStream.stochasticRound(mr.growth * relAAProd);
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + newAAs;
                aaProd(:, i) = newAAs;
                
                %fake translation
                f = rna.randStream.stochasticRound(min(rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) ./ tRNAUsage) * tRNAUsage);
                rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) = ...
                    + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) ...
                    -  f;
                rna.counts(rna.matureIndexs(rna.matureTRNAIndexs)) = ...
                    + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs)) ...
                    + f;
                aaIncorporated(:, i) = aaToTRNA * f;
                
                %simulate
                s.evolveState();
                
                %keep track
                aas(:, i) = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            end
            
            %% assertions
            % account for all amino acids
            initCounts = ...
                + initAAs ...
                + rna.baseCounts(:, m.aminoAcidIndexs)' * sum(initRNAs - rna.counts, 2) ...
                + sum(aaProd - aaIncorporated, 2);
            currCounts = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initCounts(13) = initCounts(13) + initCounts(21);
            currCounts(13) = currCounts(13) + currCounts(21);
            assertEqual(initCounts(1:20)', currCounts(1:20)');
            
            %amino acids are limiting
            assertEqual(0, min(m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)));
        end
        
        function testAminoacylation4(~)
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            import edu.stanford.covert.util.ConstantUtil;
            
            s = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'tRNAAminoacylation'
                });
            
            %% references
            comp = s.compartment;
            g = s.gene;
            
            time = s.state('Time');
            mass = s.state('Mass');
            m = s.state('Metabolite');
            rna = s.state('Rna');
            pm = s.state('ProteinMonomer');
            pc = s.state('ProteinComplex');
            mr = s.state('MetabolicReaction');
            
            ta = s.process('tRNAAminoacylation');
            
            %% constants
            [~, idxs1] = ismember(m.aminoAcidIndexs(1:20), ta.substrateGlobalIndexs);
            [~, idxs2] = ismember(rna.matureTRNAIndexs, ta.aminoacylatedRNAGlobalIndexs);
            aaToTRNA = zeros(21, 36);
            aaToTRNA(1:20, :) = ...
                -ta.reactionStoichiometryMatrix(idxs1, :) * ...
                ta.reactionModificationMatrix(:, idxs2);
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            mRNAExp = mRNAExp / sum(mRNAExp);
            tRNAUsage = ta.monomerTRNACounts' * mRNAExp;
            tRNAUsage = tRNAUsage / sum(tRNAUsage);
            
            relAAProd = aaToTRNA * tRNAUsage;
            
            initRnas = rna.expression * mass.cellInitialDryWeight * mass.dryWeightFractionRNA / ...
                (rna.expression' * rna.molecularWeights / edu.stanford.covert.util.ConstantUtil.nAvogadro);
            
            rna.decayRates = min(1, rna.decayRates);
            pm.decayRates = min(1, pm.decayRates);
            pc.decayRates = min(1, pc.decayRates);
            
            %% initialize
            rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = ...
                + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs);
            rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = 0;
            
            m.counts(m.ntpIndexs([1 3]), comp.cytosolIndexs) = 1e12;
            m.counts(m.waterIndexs, comp.cytosolIndexs) = 1e12;
            m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = 0;
            
            %% simulation
            initAAs = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initRNAs = rna.counts;
            
            iterMax = 5000;
            
            aas = zeros(21, iterMax);
            aaIncorporated = zeros(21, iterMax);
            aaProd = zeros(21, iterMax);
            decayedAAs = zeros(21, 1);
            for i = 1:iterMax
                %fake RNA synthesis and decay
                newRNAs = rna.randStream.stochasticRound((initRnas .* (log(2) / time.cellCycleLength + rna.decayRates)) * exp(log(2) * i / time.cellCycleLength));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    + newRNAs;
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    - rna.baseCounts(:, m.aminoAcidIndexs)' * newRNAs;
                
                decayedRNAs = rna.randStream.stochasticRound(rna.decayRates .* rna.counts(:, comp.cytosolIndexs));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    - decayedRNAs;
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + rna.baseCounts(:, m.aminoAcidIndexs)' * decayedRNAs;
                m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) ...
                    + m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs);
                m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs) = 0;
                
                decayedMonomers = pm.randStream.stochasticRound(pm.decayRates(:, ones(size(pm.counts, 2), 1)) .* pm.counts);
                decayedComplexs = pc.randStream.stochasticRound(pc.decayRates(:, ones(size(pc.counts, 2), 1)) .* pc.counts);
                pm.counts = pm.counts - decayedMonomers;
                pc.counts = pc.counts - decayedComplexs;
                newDecayedAAs = m.randStream.stochasticRound(...
                    + relAAProd * sum(pm.baseCounts(:, m.aminoAcidIndexs), 2)' * sum(decayedMonomers, 2) ...
                    + relAAProd * sum(pc.baseCounts(:, m.aminoAcidIndexs), 2)' * sum(decayedComplexs, 2) ...
                    );
                decayedAAs = decayedAAs + newDecayedAAs;
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + newDecayedAAs;
                m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) ...
                    + m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs);
                m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs) = 0;
                
                %fake metabolism
                mr.growth = log(2) / time.cellCycleLength * exp(log(2) * i / time.cellCycleLength);
                newAAs = m.randStream.stochasticRound(mr.growth * relAAProd);
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + newAAs;
                aaProd(:, i) = newAAs;
                
                %fake translation
                f = rna.randStream.stochasticRound(min(rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) ./ tRNAUsage) * tRNAUsage);
                rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) = ...
                    + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) ...
                    -  f;
                rna.counts(rna.matureIndexs(rna.matureTRNAIndexs)) = ...
                    + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs)) ...
                    + f;
                aaIncorporated(:, i) = aaToTRNA * f;
                
                %simulate
                s.evolveState();
                
                %keep track
                aas(:, i) = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            end
            
            %% assertions
            % account for all amino acids
            initCounts = ...
                + initAAs ...
                + rna.baseCounts(:, m.aminoAcidIndexs)' * sum(initRNAs - rna.counts, 2) ...
                + decayedAAs ...
                + sum(aaProd - aaIncorporated, 2);
            currCounts = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initCounts(13) = initCounts(13) + initCounts(21);
            currCounts(13) = currCounts(13) + currCounts(21);
            assertEqual(initCounts(1:20)', currCounts(1:20)');
            
            %amino acids are limiting
            assertEqual(0, min(min(aas(1:20, end-500+1:end))));
        end
        
        function testAminoacylation5(~)
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            import edu.stanford.covert.util.ConstantUtil;
            
            s = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'tRNAAminoacylation'
                });
            
            %% references
            comp = s.compartment;
            g = s.gene;
            
            time = s.state('Time');
            mass = s.state('Mass');
            m = s.state('Metabolite');
            rna = s.state('Rna');
            pm = s.state('ProteinMonomer');
            pc = s.state('ProteinComplex');
            mr = s.state('MetabolicReaction');
            
            ta = s.process('tRNAAminoacylation');
            
            %% constants
            [~, idxs1] = ismember(m.aminoAcidIndexs(1:20), ta.substrateGlobalIndexs);
            [~, idxs2] = ismember(rna.matureTRNAIndexs, ta.aminoacylatedRNAGlobalIndexs);
            aaToTRNA = zeros(21, 36);
            aaToTRNA(1:20, :) = ...
                -ta.reactionStoichiometryMatrix(idxs1, :) * ...
                ta.reactionModificationMatrix(:, idxs2);
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            mRNAExp = mRNAExp / sum(mRNAExp);
            tRNAUsage = ta.monomerTRNACounts' * mRNAExp;
            tRNAUsage = tRNAUsage / sum(tRNAUsage);
            
            relAAProd = aaToTRNA * tRNAUsage;
            
            initRnas = rna.expression * mass.cellInitialDryWeight * mass.dryWeightFractionRNA / ...
                (rna.expression' * rna.molecularWeights / edu.stanford.covert.util.ConstantUtil.nAvogadro);
            initExpMons = pm.randStream.stochasticRound(mRNAExp * mass.cellInitialDryWeight * mass.dryWeightFractionProtein / ...
                (mRNAExp' * pm.molecularWeights(pm.nascentIndexs) / edu.stanford.covert.util.ConstantUtil.nAvogadro));
            initMonAAs = pm.baseCounts(pm.nascentIndexs, m.aminoAcidIndexs)' * initExpMons;
            initMonAAs(13) = initMonAAs(13) + initMonAAs(21);
            initMonAAs(21) = 0;
            assertElementsAlmostEqual(initMonAAs' / sum(initMonAAs), relAAProd' / sum(relAAProd), 'relative', 5e-3)
            
            pm.decayRates(~isinf(pm.decayRates)) = log(2) / (20 * 3600);
            pc.decayRates(~isinf(pc.decayRates)) = log(2) / (20 * 3600);
            rna.decayRates = min(1, rna.decayRates);
            pm.decayRates = min(1, pm.decayRates);
            pc.decayRates = min(1, pc.decayRates);
            
            %% initialize
            rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = ...
                + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs);
            rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = 0;
            
            m.counts(m.ntpIndexs([1 3]), comp.cytosolIndexs) = 1e12;
            m.counts(m.waterIndexs, comp.cytosolIndexs) = 1e12;
            m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = 0;
            m.counts(m.getIndexs('FTHF10')) = 1e12;
            
            %% simulation
            initAAs = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initRNAs = rna.counts;
            
            iterMax = 5000;
            
            aas = zeros(21, iterMax);
            aaIncorporated = zeros(21, iterMax);
            aaProd = zeros(21, iterMax);
            decayedAAs = zeros(21, 1);
            for i = 1:iterMax
                %fake RNA synthesis and decay
                newRNAs = rna.randStream.stochasticRound((initRnas .* (log(2) / time.cellCycleLength + rna.decayRates)) * exp(log(2) * i / time.cellCycleLength));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    + newRNAs;
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    - rna.baseCounts(:, m.aminoAcidIndexs)' * newRNAs;
                
                decayedRNAs = rna.randStream.stochasticRound(rna.decayRates .* rna.counts(:, comp.cytosolIndexs));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    - decayedRNAs;
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + rna.baseCounts(:, m.aminoAcidIndexs)' * decayedRNAs;
                m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) ...
                    + m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs);
                m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs) = 0;
                
                %fake protein decay
                decayedMonomers = pm.randStream.stochasticRound(pm.decayRates(pm.nascentIndexs) .* initExpMons);
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + pm.baseCounts(pm.nascentIndexs, m.aminoAcidIndexs)' * sum(decayedMonomers, 2);
                m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) ...
                    + m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs);
                m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs) = 0;
                decayedAAs = decayedAAs + pm.baseCounts(pm.nascentIndexs, m.aminoAcidIndexs)' * sum(decayedMonomers, 2);
                
                %fake metabolism
                mr.growth = log(2) / time.cellCycleLength * exp(log(2) * i / time.cellCycleLength);
                newAAs = m.randStream.stochasticRound(mr.growth * relAAProd);
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + newAAs;
                aaProd(:, i) = newAAs;
                
                %fake translation
                f = rna.randStream.stochasticRound(min(rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) ./ tRNAUsage) * tRNAUsage);
                rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) = ...
                    + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs)) ...
                    -  f;
                rna.counts(rna.matureIndexs(rna.matureTRNAIndexs)) = ...
                    + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs)) ...
                    + f;
                aaIncorporated(:, i) = aaToTRNA * f;
                
                %simulate
                s.evolveState();
                
                %keep track
                aas(:, i) = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            end
            
            %% assertions
            % account for all amino acids
            initCounts = ...
                + initAAs ...
                + rna.baseCounts(:, m.aminoAcidIndexs)' * sum(initRNAs - rna.counts, 2) ...
                + decayedAAs ...
                + sum(aaProd - aaIncorporated, 2);
            currCounts = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initCounts(13) = initCounts(13) + initCounts(21);
            currCounts(13) = currCounts(13) + currCounts(21);
            assertEqual(initCounts(1:20)', currCounts(1:20)');
            
            %amino acids are limiting
            assertEqual(0, min(min(aas(1:20, end-500+1:end))));
        end
        
        function testAminoacylation6(~)
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            import edu.stanford.covert.util.ConstantUtil;
            
            s = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'tRNAAminoacylation'
                'Translation'
                });
            
            %% references
            comp = s.compartment;
            g = s.gene;
            
            time = s.state('Time');
            mass = s.state('Mass');
            m = s.state('Metabolite');
            rna = s.state('Rna');
            pm = s.state('ProteinMonomer');
            pc = s.state('ProteinComplex');
            mr = s.state('MetabolicReaction');
            pol = s.state('Polypeptide');
            
            ta = s.process('tRNAAminoacylation');
            trl = s.process('Translation');
            
            %% constants
            [~, idxs1] = ismember(m.aminoAcidIndexs(1:20), ta.substrateGlobalIndexs);
            [~, idxs2] = ismember(rna.matureTRNAIndexs, ta.aminoacylatedRNAGlobalIndexs);
            aaToTRNA = zeros(21, 36);
            aaToTRNA(1:20, :) = ...
                -ta.reactionStoichiometryMatrix(idxs1, :) * ...
                ta.reactionModificationMatrix(:, idxs2);
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            mRNAExp = mRNAExp / sum(mRNAExp);
            tRNAUsage = ta.monomerTRNACounts' * mRNAExp;
            tRNAUsage = tRNAUsage / sum(tRNAUsage);
            
            relAAProd = aaToTRNA * tRNAUsage;
            
            initRnas = rna.expression * mass.cellInitialDryWeight * mass.dryWeightFractionRNA / ...
                (rna.expression' * rna.molecularWeights / edu.stanford.covert.util.ConstantUtil.nAvogadro);
            initExpMons = pm.randStream.stochasticRound(mRNAExp * mass.cellInitialDryWeight * mass.dryWeightFractionProtein / ...
                (mRNAExp' * pm.molecularWeights(pm.nascentIndexs) / edu.stanford.covert.util.ConstantUtil.nAvogadro));
            initMonAAs = pm.baseCounts(pm.nascentIndexs, m.aminoAcidIndexs)' * initExpMons;
            initMonAAs(13) = initMonAAs(13) + initMonAAs(21);
            initMonAAs(21) = 0;
            assertElementsAlmostEqual(initMonAAs' / sum(initMonAAs), relAAProd' / sum(relAAProd), 'relative', 5e-3)
            
            pm.decayRates(~isinf(pm.decayRates)) = log(2) / (20 * 3600);
            pc.decayRates(~isinf(pc.decayRates)) = log(2) / (20 * 3600);
            rna.decayRates = min(1, rna.decayRates);
            pm.decayRates = min(1, pm.decayRates);
            pc.decayRates = min(1, pc.decayRates);
            
            if ~isempty(trl)
                trl.tmRNABindingProbability = 0;
            end
            
            %% initialize
            rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = ...
                + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs);
            rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = 0;
            initRnas(rna.matureIndexs(rna.matureMRNAIndexs)) = 1e6 * initRnas(rna.matureIndexs(rna.matureMRNAIndexs));
            rna.counts(rna.matureIndexs(rna.matureMRNAIndexs), comp.cytosolIndexs) = rna.randStream.stochasticRound(initRnas(rna.matureIndexs(rna.matureMRNAIndexs)));
            
            m.counts(m.ntpIndexs([1 3]), comp.cytosolIndexs) = 1e12;
            m.counts(m.waterIndexs, comp.cytosolIndexs) = 1e12;
            m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = 0;
            m.counts(m.getIndexs('FTHF10')) = 1e12;
            
            %% simulation
            initAAs = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initRNAs = rna.counts;
            initMons = pm.counts;
            initCpxs = pc.counts;
            initAAsInPolypeptides = pol.totalBaseCounts;
            
            iterMax = 2000;
            
            aas = zeros(21, iterMax);
            aaProd = zeros(21, iterMax);
            decayedAAs = zeros(21, 1);
            for i = 1:iterMax
                %fake RNA synthesis and decay
                newRNAs = rna.randStream.stochasticRound((initRnas .* (log(2) / time.cellCycleLength + rna.decayRates)) * exp(log(2) * i / time.cellCycleLength));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    + newRNAs;
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    - rna.baseCounts(:, m.aminoAcidIndexs)' * newRNAs;
                
                decayedRNAs = rna.randStream.stochasticRound(rna.decayRates .* rna.counts(:, comp.cytosolIndexs));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    - decayedRNAs;
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + rna.baseCounts(:, m.aminoAcidIndexs)' * decayedRNAs;
                m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) ...
                    + m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs);
                m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs) = 0;
                
                %fake protein decay
                decayedMonomers = pm.randStream.stochasticRound(pm.decayRates(pm.nascentIndexs) .* initExpMons);
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + pm.baseCounts(pm.nascentIndexs, m.aminoAcidIndexs)' * sum(decayedMonomers, 2);
                m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) ...
                    + m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs);
                m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs) = 0;
                decayedAAs = decayedAAs + pm.baseCounts(pm.nascentIndexs, m.aminoAcidIndexs)' * sum(decayedMonomers, 2);
                
                %fake metabolism
                mr.growth = log(2) / time.cellCycleLength * exp(log(2) * i / time.cellCycleLength);
                newAAs = m.randStream.stochasticRound(mr.growth * relAAProd);
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + newAAs;
                aaProd(:, i) = newAAs;
                
                %simulate
                s.evolveState();
                
                %keep track
                aas(:, i) = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            end
            
            %% assertions
            % account for all amino acids
            initCounts = ...
                + initAAs ...
                + rna.baseCounts(:, m.aminoAcidIndexs)' * sum(initRNAs - rna.counts, 2) ...
                + pm.baseCounts(:, m.aminoAcidIndexs)' * sum(initMons - pm.counts, 2) ...
                + pc.baseCounts(:, m.aminoAcidIndexs)' * sum(initCpxs - pc.counts, 2) ...
                + decayedAAs ...
                + initAAsInPolypeptides' - pol.totalBaseCounts' ...
                + sum(aaProd, 2);
            currCounts = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initCounts(13) = initCounts(13) + initCounts(21);
            currCounts(13) = currCounts(13) + currCounts(21);
            assertEqual(initCounts(1:20)', currCounts(1:20)');
            
            %amino acids are limiting
            assertEqual(0, min(min(aas(1:20, end-500+1:end))));
        end
        
        function testAminoacylation7(~)
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            import edu.stanford.covert.util.ConstantUtil;
            
            s = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'tRNAAminoacylation'
                'Translation'
                });
            
            %% references
            comp = s.compartment;
            g = s.gene;
            
            time = s.state('Time');
            mass = s.state('Mass');
            m = s.state('Metabolite');
            rna = s.state('Rna');
            pm = s.state('ProteinMonomer');
            pc = s.state('ProteinComplex');
            mr = s.state('MetabolicReaction');
            pol = s.state('Polypeptide');
            
            ta = s.process('tRNAAminoacylation');
            trl = s.process('Translation');
            
            %% constants
            [~, idxs1] = ismember(m.aminoAcidIndexs(1:20), ta.substrateGlobalIndexs);
            [~, idxs2] = ismember(rna.matureTRNAIndexs, ta.aminoacylatedRNAGlobalIndexs);
            aaToTRNA = zeros(21, 36);
            aaToTRNA(1:20, :) = ...
                -ta.reactionStoichiometryMatrix(idxs1, :) * ...
                ta.reactionModificationMatrix(:, idxs2);
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            mRNAExp = mRNAExp / sum(mRNAExp);
            tRNAUsage = ta.monomerTRNACounts' * mRNAExp;
            tRNAUsage = tRNAUsage / sum(tRNAUsage);
            
            relAAProd = aaToTRNA * tRNAUsage;
            
            initRnas = rna.expression * mass.cellInitialDryWeight * mass.dryWeightFractionRNA / ...
                (rna.expression' * rna.molecularWeights / edu.stanford.covert.util.ConstantUtil.nAvogadro);
            initExpMons = pm.randStream.stochasticRound(mRNAExp * mass.cellInitialDryWeight * mass.dryWeightFractionProtein / ...
                (mRNAExp' * pm.molecularWeights(pm.nascentIndexs) / edu.stanford.covert.util.ConstantUtil.nAvogadro));
            initMonAAs = pm.baseCounts(pm.nascentIndexs, m.aminoAcidIndexs)' * initExpMons;
            initMonAAs(13) = initMonAAs(13) + initMonAAs(21);
            initMonAAs(21) = 0;
            assertElementsAlmostEqual(initMonAAs' / sum(initMonAAs), relAAProd' / sum(relAAProd), 'relative', 5e-3)
            
            pm.decayRates(~isinf(pm.decayRates)) = log(2) / (20 * 3600);
            pc.decayRates(~isinf(pc.decayRates)) = log(2) / (20 * 3600);
            rna.decayRates = min(1, rna.decayRates);
            pm.decayRates = min(1, pm.decayRates);
            pc.decayRates = min(1, pc.decayRates);
            
            if ~isempty(trl)
                trl.tmRNABindingProbability = 0;
            end
            
            %% initialize
            rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = ...
                + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs);
            rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = 0;
            initRnas(rna.matureIndexs(rna.matureMRNAIndexs)) = 1e6 * initRnas(rna.matureIndexs(rna.matureMRNAIndexs));
            rna.counts(rna.matureIndexs(rna.matureMRNAIndexs)) = rna.randStream.stochasticRound(initRnas(rna.matureIndexs(rna.matureMRNAIndexs)));
            
            m.counts(m.ntpIndexs([1 3]), comp.cytosolIndexs) = 1e12;
            m.counts(m.waterIndexs, comp.cytosolIndexs) = 1e12;
            m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = 0;
            m.counts(m.getIndexs('FTHF10')) = 1e12;
            
            %% simulation
            initAAs = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initRNAs = rna.counts;
            initMons = pm.counts;
            initCpxs = pc.counts;
            initAAsInPolypeptides = pol.totalBaseCounts;
            
            iterMax = 20000;
            
            aas = zeros(21, iterMax);
            aaProd = zeros(21, iterMax);
            for i = 1:iterMax
                %fake RNA synthesis and decay
                newRNAs = rna.randStream.stochasticRound((initRnas .* (log(2) / time.cellCycleLength + rna.decayRates)) * exp(log(2) * i / time.cellCycleLength));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    + newRNAs;
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    - rna.baseCounts(:, m.aminoAcidIndexs)' * newRNAs;
                
                decayedRNAs = rna.randStream.stochasticRound(rna.decayRates .* rna.counts(:, comp.cytosolIndexs));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    - decayedRNAs;
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + rna.baseCounts(:, m.aminoAcidIndexs)' * decayedRNAs;
                m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) ...
                    + m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs);
                m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs) = 0;
                
                %fake protein decay
                decayedMonomers = min(pm.counts, pm.randStream.stochasticRound(pm.decayRates(:, ones(6, 1)) .* pm.counts));
                decayedComplexs = min(pc.counts, pc.randStream.stochasticRound(pc.decayRates(:, ones(6, 1)) .* pc.counts));
                pm.counts = pm.counts - decayedMonomers;
                pc.counts = pc.counts - decayedComplexs;
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + pm.baseCounts(:, m.aminoAcidIndexs)' * sum(decayedMonomers, 2) ...
                    + pc.baseCounts(:, m.aminoAcidIndexs)' * sum(decayedComplexs, 2);
                m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) ...
                    + m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs);
                m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs) = 0;
                
                %fake metabolism
                mr.growth = log(2) / time.cellCycleLength * exp(log(2) * i / time.cellCycleLength);
                newAAs = m.randStream.stochasticRound(mr.growth * relAAProd);
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + newAAs;
                aaProd(:, i) = newAAs;
                
                %simulate
                s.evolveState();
                
                %keep track
                aas(:, i) = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            end
            
            %% assertions
            % account for all amino acids
            initCounts = ...
                + initAAs ...
                + rna.baseCounts(:, m.aminoAcidIndexs)' * sum(initRNAs - rna.counts, 2) ...
                + pm.baseCounts(:, m.aminoAcidIndexs)' * sum(initMons - pm.counts, 2) ...
                + pc.baseCounts(:, m.aminoAcidIndexs)' * sum(initCpxs - pc.counts, 2) ...
                + initAAsInPolypeptides' - pol.totalBaseCounts' ...
                + sum(aaProd, 2);
            currCounts = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initCounts(13) = initCounts(13) + initCounts(21);
            currCounts(13) = currCounts(13) + currCounts(21);
            assertEqual(initCounts(1:20)', currCounts(1:20)');
            
            %amino acids are limiting
            assertEqual(0, min(min(aas(1:20, end-500+1:end))));
        end
        
        function testAminoacylation8(~)
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            import edu.stanford.covert.util.ConstantUtil;
            
            s = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'tRNAAminoacylation'
                'Translation'
                'MacromolecularComplexation'
                });
            
            %% references
            comp = s.compartment;
            g = s.gene;
            
            time = s.state('Time');
            mass = s.state('Mass');
            m = s.state('Metabolite');
            rna = s.state('Rna');
            pm = s.state('ProteinMonomer');
            pc = s.state('ProteinComplex');
            mr = s.state('MetabolicReaction');
            pol = s.state('Polypeptide');
            
            ta = s.process('tRNAAminoacylation');
            trl = s.process('Translation');
            
            %% constants
            [~, idxs1] = ismember(m.aminoAcidIndexs(1:20), ta.substrateGlobalIndexs);
            [~, idxs2] = ismember(rna.matureTRNAIndexs, ta.aminoacylatedRNAGlobalIndexs);
            aaToTRNA = zeros(21, 36);
            aaToTRNA(1:20, :) = ...
                -ta.reactionStoichiometryMatrix(idxs1, :) * ...
                ta.reactionModificationMatrix(:, idxs2);
            
            mRNAExp = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);
            mRNAExp = mRNAExp / sum(mRNAExp);
            tRNAUsage = ta.monomerTRNACounts' * mRNAExp;
            tRNAUsage = tRNAUsage / sum(tRNAUsage);
            
            relAAProd = aaToTRNA * tRNAUsage;
            
            initRnas = rna.expression * mass.cellInitialDryWeight * mass.dryWeightFractionRNA / ...
                (rna.expression' * rna.molecularWeights / edu.stanford.covert.util.ConstantUtil.nAvogadro);
            initExpMons = pm.randStream.stochasticRound(mRNAExp * mass.cellInitialDryWeight * mass.dryWeightFractionProtein / ...
                (mRNAExp' * pm.molecularWeights(pm.nascentIndexs) / edu.stanford.covert.util.ConstantUtil.nAvogadro));
            initMonAAs = pm.baseCounts(pm.nascentIndexs, m.aminoAcidIndexs)' * initExpMons;
            initMonAAs(13) = initMonAAs(13) + initMonAAs(21);
            initMonAAs(21) = 0;
            assertElementsAlmostEqual(initMonAAs' / sum(initMonAAs), relAAProd' / sum(relAAProd), 'relative', 5e-3)
            
            pm.decayRates(~isinf(pm.decayRates)) = log(2) / (20 * 3600);
            pc.decayRates(~isinf(pc.decayRates)) = log(2) / (20 * 3600);
            rna.decayRates = min(1, rna.decayRates);
            pm.decayRates = min(1, pm.decayRates);
            pc.decayRates = min(1, pc.decayRates);
            
            if ~isempty(trl)
                trl.tmRNABindingProbability = 0;
            end
            
            %% initialize
            rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = ...
                + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs);
            rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = 0;
            initRnas(rna.matureIndexs(rna.matureMRNAIndexs)) = 1e6 * initRnas(rna.matureIndexs(rna.matureMRNAIndexs));
            rna.counts(rna.matureIndexs(rna.matureMRNAIndexs)) = rna.randStream.stochasticRound(initRnas(rna.matureIndexs(rna.matureMRNAIndexs)));
            
            m.counts(m.ntpIndexs([1 3]), comp.cytosolIndexs) = 1e12;
            m.counts(m.waterIndexs, comp.cytosolIndexs) = 1e12;
            m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = 0;
            m.counts(m.getIndexs('FTHF10')) = 1e12;
            
            %% simulation
            initAAs = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initRNAs = rna.counts;
            initMons = pm.counts;
            initCpxs = pc.counts;
            initAAsInPolypeptides = pol.totalBaseCounts;
            
            iterMax = 2000;
            
            aas = zeros(21, iterMax);
            aaProd = zeros(21, iterMax);
            for i = 1:iterMax
                %fake RNA synthesis and decay
                newRNAs = rna.randStream.stochasticRound((initRnas .* (log(2) / time.cellCycleLength + rna.decayRates)) * exp(log(2) * i / time.cellCycleLength));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    + newRNAs;
                m.counts(:, comp.cytosolIndexs) = ...
                    + m.counts(:, comp.cytosolIndexs) ...
                    - rna.baseCounts' * newRNAs;
                
                decayedRNAs = rna.randStream.stochasticRound(rna.decayRates .* rna.counts(:, comp.cytosolIndexs));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    - decayedRNAs;
                m.counts(:, comp.cytosolIndexs) = ...
                    + m.counts(:, comp.cytosolIndexs) ...
                    + rna.baseCounts' * decayedRNAs;
                
                %fake protein modification
                m.counts(:, comp.cytosolIndexs) = ...
                    + m.counts(:, comp.cytosolIndexs) ...
                    + (pm.baseCounts(pm.nascentIndexs, :) - pm.baseCounts(pm.matureIndexs, :))' * pm.counts(pm.nascentIndexs, comp.cytosolIndexs);
                pm.counts(sub2ind(size(pm.counts), pm.matureIndexs, pm.compartments(pm.matureIndexs))) = ...
                    + pm.counts(sub2ind(size(pm.counts), pm.matureIndexs, pm.compartments(pm.matureIndexs))) ...
                    + pm.counts(pm.nascentIndexs, comp.cytosolIndexs);
                pm.counts(pm.nascentIndexs, comp.cytosolIndexs) = 0;
                
                m.counts(:, comp.cytosolIndexs) = ...
                    + m.counts(:, comp.cytosolIndexs) ...
                    + (pc.baseCounts(pc.nascentIndexs, :) - pc.baseCounts(pc.matureIndexs, :))' * sum(pc.counts(pc.nascentIndexs, :), 2);
                pc.counts(pc.matureIndexs, :) = ...
                    + pc.counts(pc.matureIndexs, :) ...
                    + pc.counts(pc.nascentIndexs, :);
                pc.counts(pc.nascentIndexs, :) = 0;
                
                %fake protein decay
                decayedMonomers = min(pm.counts, pm.randStream.stochasticRound(pm.decayRates(:, ones(6, 1)) .* pm.counts));
                decayedComplexs = min(pc.counts, pc.randStream.stochasticRound(pc.decayRates(:, ones(6, 1)) .* pc.counts));
                pm.counts = pm.counts - decayedMonomers;
                pc.counts = pc.counts - decayedComplexs;
                m.counts(:, comp.cytosolIndexs) = ...
                    + m.counts(:, comp.cytosolIndexs) ...
                    + pm.baseCounts' * sum(decayedMonomers, 2) ...
                    + pc.baseCounts' * sum(decayedComplexs, 2);
                
                %fake metabolism
                m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) ...
                    + m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs);
                m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs) = 0;
                
                mr.growth = log(2) / time.cellCycleLength * exp(log(2) * i / time.cellCycleLength);
                newAAs = m.randStream.stochasticRound(mr.growth * relAAProd);
                m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) ...
                    + newAAs;
                aaProd(:, i) = newAAs;
                
                %simulate
                s.evolveState();
                
                %keep track
                aas(:, i) = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            end
            
            %% assertions
            % account for all amino acids
            initCounts = ...
                + initAAs ...
                + rna.baseCounts(:, m.aminoAcidIndexs)' * sum(initRNAs - rna.counts, 2) ...
                + pm.baseCounts(:, m.aminoAcidIndexs)' * sum(initMons - pm.counts, 2) ...
                + pc.baseCounts(:, m.aminoAcidIndexs)' * sum(initCpxs - pc.counts, 2) ...
                + initAAsInPolypeptides' - pol.totalBaseCounts' ...
                + sum(aaProd, 2);
            currCounts = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initCounts(13) = initCounts(13) + initCounts(21);
            currCounts(13) = currCounts(13) + currCounts(21);
            assertEqual(initCounts(1:20)', currCounts(1:20)');
            
            %amino acids are limiting
            assertIn(min(min(aas(1:20, end-500+1:end))), [-Inf 0]); %negative counts possible because mocking protein modification isn't limited by amino acids
        end
        
        function testAminoacylation9(~)
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            import edu.stanford.covert.util.ConstantUtil;
            
            s = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'tRNAAminoacylation'
                'RibosomeAssembly'
                'Translation'
                'ProteinProcessingI'
                'ProteinTranslocation'
                'ProteinProcessingII'
                'ProteinFolding'
                'ProteinModification'
                'MacromolecularComplexation'
                'TerminalOrganelleAssembly'
                'ProteinActivation'
                'ProteinDecay'
                'RNADecay'
                });
            
            %% references
            comp = s.compartment;
            
            time = s.state('Time');
            mass = s.state('Mass');
            m = s.state('Metabolite');
            rna = s.state('Rna');
            pm = s.state('ProteinMonomer');
            pc = s.state('ProteinComplex');
            mr = s.state('MetabolicReaction');
            pol = s.state('Polypeptide');
            
            trl = s.process('Translation');
            
            %% constants
            initRnas = rna.expression * mass.cellInitialDryWeight * mass.dryWeightFractionRNA / ...
                (rna.expression' * rna.molecularWeights / edu.stanford.covert.util.ConstantUtil.nAvogadro);
            
            rna.decayRates = min(1, rna.decayRates);
            pm.decayRates = min(1, pm.decayRates);
            pc.decayRates = min(1, pc.decayRates);
            
            trl.tmRNABindingProbability = 0;
            
            %% initialize
            rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = ...
                + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) ...
                + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs);
            rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) = 0;
            
            m.counts(m.ntpIndexs([1 3]), comp.cytosolIndexs) = 1e12;
            m.counts(m.waterIndexs, comp.cytosolIndexs) = 1e12;
            m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = 0;
            m.counts(m.getIndexs('FTHF10')) = 1e12;
            
            %% simulation
            initAAs = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initRNAs = rna.counts;
            initMons = pm.counts;
            initCpxs = pc.counts;
            initAAsInPolypeptides = pol.totalBaseCounts;
            
            iterMax = 2000;
            
            aas = zeros(21, iterMax);
            aaProd = zeros(21, iterMax);
            for i = 1:iterMax
                %fake RNA synthesis and decay
                newRNAs = rna.randStream.stochasticRound((initRnas .* (log(2) / time.cellCycleLength + rna.decayRates)) * exp(log(2) * i / time.cellCycleLength));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    + newRNAs;
                m.counts(:, comp.cytosolIndexs) = ...
                    + m.counts(:, comp.cytosolIndexs) ...
                    - rna.baseCounts' * newRNAs;
                
                %fake metabolism
                m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) = ...
                    + m.counts(m.aminoAcidIndexs(13), comp.cytosolIndexs) ...
                    + m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs);
                m.counts(m.aminoAcidIndexs(21), comp.cytosolIndexs) = 0;
                
                mr.growth = log(2) / time.cellCycleLength * exp(log(2) * i / time.cellCycleLength);
                newMets = max(0, m.randStream.stochasticRound(mr.growth * (m.biomassProduction - m.byproducts)));
                m.counts = m.counts + newMets;
                aaProd(:, i) = newMets(m.aminoAcidIndexs, comp.cytosolIndexs);
                
                %simulate
                s.evolveState();
                
                %keep track
                aas(:, i) = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            end
            
            %% assertions
            % account for all amino acids
            initCounts = ...
                + initAAs ...
                + rna.baseCounts(:, m.aminoAcidIndexs)' * sum(initRNAs - rna.counts, 2) ...
                + pm.baseCounts(:, m.aminoAcidIndexs)' * sum(initMons - pm.counts, 2) ...
                + pc.baseCounts(:, m.aminoAcidIndexs)' * sum(initCpxs - pc.counts, 2) ...
                + initAAsInPolypeptides' - pol.totalBaseCounts' ...
                + sum(aaProd, 2);
            currCounts = m.counts(m.aminoAcidIndexs, comp.cytosolIndexs);
            initCounts(13) = initCounts(13) + initCounts(21);
            currCounts(13) = currCounts(13) + currCounts(21);
            assertElementsAlmostEqual(initCounts(1:20)', currCounts(1:20)', 'relative', 10e-2, 1000);
            
            %amino acids are limiting
            assertIn(min(min(aas(1:20, end-500+1:end))), [0 1]);
        end
        
        function testSynthesis(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            %% references
            s = this.simulation;
            
            comp = s.compartment;
            g = s.gene;
            
            time = s.state('Time');
            mass = s.state('Mass');
            m = s.state('Metabolite');
            rna = s.state('Rna');
            pm = s.state('ProteinMonomer');
            pc = s.state('ProteinComplex');
            mr = s.state('MetabolicReaction');
            
            trl = s.process('Translation');
            ta = s.process('tRNAAminoacylation');
            
            %% constants
            rnaExp = rna.expression;
            rnaMWs = rna.molecularWeights;
            monMWs = pm.molecularWeights(pm.matureIndexs);
            
            pm.decayRates(setdiff(1:end, pm.signalSequenceIndexs)) = 0;
            pc.decayRates(:) = 0;
            
            %% simulation
            m.counts(m.aminoAcidIndexs, comp.cytosolIndexs) = 0;
            m.counts(m.atpIndexs, 1) = 1e6;
            m.counts(m.getIndexs('GTP'), 1) = 1e6;
            m.counts(m.waterIndexs, 1) = 1e6;
            initRNAs = rnaExp * (mass.cellInitialDryWeight * mass.dryWeightFractionRNA) / (rnaExp' * rnaMWs / ConstantUtil.nAvogadro);
            
            initMonomers = pm.counts;
            initComplexs = pc.counts;
            
            iterMax = 1000;
            
            growth = 0;
            simMrnaExp = zeros(size(rna.matureMRNAIndexs));
            aas = zeros(20, iterMax);
            for i = 1:iterMax
                %fake RNA synthesis and decay
                rna.counts(:, comp.cytosolIndexs) = ...
                    rna.counts(:, comp.cytosolIndexs) + ...
                    rna.randStream.stochasticRound(...
                    initRNAs * log(2) / time.cellCycleLength * exp(log(2) * i / time.cellCycleLength) + ...
                    initRNAs .* rna.decayRates * exp(log(2) * i / time.cellCycleLength));
                rna.counts(:, comp.cytosolIndexs) = ...
                    rna.counts(:, comp.cytosolIndexs) - ...
                    rna.randStream.stochasticRound(min(1, rna.decayRates) .* rna.counts(:, comp.cytosolIndexs));
                rna.counts(rna.matureIndexs(rna.matureMRNAIndexs), comp.cytosolIndexs) = round(1e3 * initRNAs(rna.matureIndexs(rna.matureMRNAIndexs)));
                
                %simulate
                s.evolveState();
                
                growth = growth + mr.growth;
                simMrnaExp = simMrnaExp + rna.counts(rna.matureIndexs(rna.matureMRNAIndexs), comp.cytosolIndexs);
                aas(:, i) = ta.substrates(ta.substrateIndexs_aminoAcids);
            end
            
            newMonomers = pm.counts - initMonomers;
            newComplexs = pc.counts - initComplexs;
            
            totNewMonomers = sum(...
                + newMonomers(pm.nascentIndexs, :) ...
                + newMonomers(pm.processedIIndexs, :) ...
                + newMonomers(pm.processedIIIndexs, :) ...
                + newMonomers(pm.foldedIndexs, :) ...
                + newMonomers(pm.matureIndexs, :) ...
                + newMonomers(pm.boundIndexs, :) ...
                + newMonomers(pm.misfoldedIndexs, :) ...
                + newMonomers(pm.damagedIndexs, :) ...
                + newMonomers(pm.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.nascentIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.matureIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.boundIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.misfoldedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.damagedIndexs, :), ...
                2);
            
            expProd = rna.matureRNAGeneComposition(g.mRNAIndexs, rna.matureMRNAIndexs) * simMrnaExp;
            
            %% assertions
            %amino acids are limiting
            assertEqual(0, min(trl.aminoacylatedTRNAs));
            assertIn(min(min(aas(:, end-100:end))), [0 0]);
            
            %amino acids are mostly used -- not too many are free
            assertIn(...
                (m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)' * m.molecularWeights(m.aminoAcidIndexs)) / ...
                sum(pm.counts' * pm.molecularWeights + pc.counts' * pc.molecularWeights), ...
                [0 0.01]);
            
            %mature/bound monomer counts increased
            assertIn(min(totNewMonomers), [0 Inf]);
            if ~all(expProd)
                assertAllEqual(0, totNewMonomers(expProd == 0));
            end
            assertElementsAlmostEqual(totNewMonomers' * monMWs / ConstantUtil.nAvogadro, ...
                growth * mass.cellInitialDryWeight * mass.dryWeightFractionProtein, ...
                'relative', 20e-2, 0);
            assertElementsAlmostEqual(...
                (m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)' * m.molecularWeights(m.aminoAcidIndexs) + ...
                totNewMonomers' * monMWs) / ConstantUtil.nAvogadro, ...
                growth * mass.cellInitialDryWeight * mass.dryWeightFractionProtein, ...
                'relative', 10e-2, 0);
            
            assertIn(sum(sum(pm.counts(pm.nascentIndexs, :))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.processedIIndexs, :))), [0 20]);
            assertIn(sum(sum(pm.counts(pm.processedIIIndexs, :))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.signalSequenceIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.foldedIndexs, :))), [0 20]);
            assertIn(sum(sum(pm.counts(pm.misfoldedIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.damagedIndexs, :))), [0 10]);
            
            assertIn(sum(sum(pc.counts(pc.nascentIndexs, :))), [0 10]);
            assertIn(sum(sum(pc.counts(pc.misfoldedIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pc.counts(pc.damagedIndexs, :))), [0 10]);
            
            %RNA production / expression matches expectations
            assertIn(corr(expProd, totNewMonomers), [0.75 1]);
            assertIn(180 / pi * acos((expProd' * totNewMonomers) / (sqrt(expProd' * expProd) * sqrt(totNewMonomers' * totNewMonomers))), ...
                [0 35]);
        end
        
        function testSynthesisAndDecay(this)
            import edu.stanford.covert.util.ConstantUtil;
            
            %% references
            s = this.simulation;
            %s.applyOptions('verbosity', 1);
            %this.seedSimulation(round(mod(now, 1) * 1e7));
            
            comp = s.compartment;
            g = s.gene;
            
            time = s.state('Time');
            mass = s.state('Mass');
            m = s.state('Metabolite');
            rna = s.state('Rna');
            pm = s.state('ProteinMonomer');
            pc = s.state('ProteinComplex');
            mr = s.state('MetabolicReaction');
            pol = s.state('Polypeptide');
            
            dcy = s.process('ProteinDecay');
            ta = s.process('tRNAAminoacylation');
            
            %% constants
            rnaExp = rna.expression;
            rnaMWs = rna.molecularWeights;
            monMWs = pm.molecularWeights(pm.matureIndexs);
            
            %% simulation
            initRNAs = rnaExp * (mass.cellInitialDryWeight * mass.dryWeightFractionRNA) / (rnaExp' * rnaMWs / ConstantUtil.nAvogadro);
            initRNAs(rna.matureIndexs(rna.matureMRNAIndexs)) = ...
                10 * initRNAs(rna.matureIndexs(rna.matureMRNAIndexs));
            rna.counts(rna.matureIndexs(rna.matureMRNAIndexs)) = ...
                rna.counts(rna.matureIndexs(rna.matureMRNAIndexs)) + ...
                rna.randStream.stochasticRound((10-1) * initRNAs(rna.matureIndexs(rna.matureMRNAIndexs)));
            
            initGrowth = mr.growth;
            initMetabolites = m.counts;
            initMonomers = pm.counts;
            initComplexs = pc.counts;
            iterMax = 10000;
            
            growth = 0;
            simRNAExp = zeros(size(rna.counts, 1), 1);
            simComplexExp = zeros(size(pc.counts, 1), 1);
            totRnas = zeros(numel(rna.matureMRNAIndexs), iterMax);
            tRNAs = zeros(numel(rna.matureTRNAIndexs), iterMax);
            freedAAs = zeros(numel(dcy.substrateIndexs_aminoAcids), iterMax);
            aaCnts = zeros(numel(ta.substrateIndexs_aminoAcids), iterMax);
            frac = zeros(iterMax, 1);
            nAbortedPolypeptides = zeros(iterMax, 1);
            oldNAbortedPolypeptides = 0;
            if s.verbosity > 0
                fprintf('%5s %8s %8s\n', 'Iter ', ' Min AA ', ' Sum AA ');
                fprintf('%5s %8s %8s\n', '=====', '========', '========');
            end
            for i = 1:iterMax
                if s.verbosity > 0 && mod(i, 100) == 1
                    [val, idx] = min(m.counts(m.aminoAcidIndexs(1:20), comp.cytosolIndexs));
                    fprintf('%5d %3s %4d %8d\n', i, ...
                        m.wholeCellModelIDs{m.aminoAcidIndexs(idx)}, val, ...
                        sum(m.counts(m.aminoAcidIndexs(1:20), comp.cytosolIndexs)));
                end
                
                %fake RNA synthesis and decay
                newRNAs = rna.randStream.stochasticRound(...
                    initRNAs * log(2) / time.cellCycleLength * exp(log(2) * i / time.cellCycleLength) + ...
                    initRNAs .* rna.decayRates * exp(log(2) * i / time.cellCycleLength));
                rnaDecays = rna.randStream.stochasticRound(min(1, rna.decayRates) .* rna.counts(:, comp.cytosolIndexs));
                rna.counts(:, comp.cytosolIndexs) = ...
                    + rna.counts(:, comp.cytosolIndexs) ...
                    + newRNAs ...
                    - rnaDecays;
                m.counts(m.aminoAcidIndexs(1:20), comp.cytosolIndexs) = ...
                    m.counts(m.aminoAcidIndexs(1:20), comp.cytosolIndexs) + ...
                    max(0, -ta.reactionStoichiometryMatrix(ta.substrateIndexs_aminoAcids, :)) * ta.reactionModificationMatrix * ...
                    rnaDecays(rna.aminoacylatedIndexs(ta.aminoacylatedRNAGlobalIndexs));
                
                %simulate
                s.evolveState();
                
                growth = growth + mr.growth;
                simRNAExp = simRNAExp + rna.counts(:, comp.cytosolIndexs);
                simComplexExp = simComplexExp + pc.counts(:, comp.cytosolIndexs);
                totRnas(i) = sum(rna.counts, 2)' * rna.molecularWeights - ...
                    sum(rna.counts(rna.matureIndexs(rna.matureMRNAIndexs), :), 2)' * rna.molecularWeights(rna.matureIndexs(rna.matureMRNAIndexs));
                tRNAs(:, i) = ...
                    + rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs) ...
                    + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), comp.cytosolIndexs);
                freedAAs(:, i) = dcy.substrates(dcy.substrateIndexs_aminoAcids);
                aaCnts(:, i) = ta.substrates(ta.substrateIndexs_aminoAcids);
                frac(i) = (m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)' * m.molecularWeights(m.aminoAcidIndexs)) / ...
                    (sum(pm.dryWeight + pc.dryWeight + pol.dryWeight) * ConstantUtil.nAvogadro);
                nAbortedPolypeptides(i) = max(0, numel(pol.abortedSequences) - oldNAbortedPolypeptides);
                oldNAbortedPolypeptides = numel(pol.abortedSequences);
                
                this.assertStateIsValid();
            end
            
            simMatureRNAExp = ...
                + simRNAExp(rna.matureIndexs) ...
                + simRNAExp(rna.aminoacylatedIndexs) ...
                + simRNAExp(rna.boundIndexs);
            simMatureRNAExp(setdiff(1:end, rna.matureMRNAIndexs)) = ...
                + simMatureRNAExp(setdiff(1:end, rna.matureMRNAIndexs)) ...
                + sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * (...
                + simComplexExp(pc.matureIndexs) ...
                + simComplexExp(pc.boundIndexs));
            
            %% assertions - growth
            assertElementsAlmostEqual(log(2) / time.cellCycleLength, initGrowth, 'relative', 0.3);
            assertElementsAlmostEqual(initGrowth * exp(log(2) * i / time.cellCycleLength), ...
                mr.growth, 'relative', 0.5, 0);
            assertElementsAlmostEqual(initGrowth * time.cellCycleLength / log(2) * (exp(log(2) * i / time.cellCycleLength) - 1), ...
                growth, 'relative', 0.3, 0);
            assertElementsAlmostEqual(mass.cellInitialDryWeight * initGrowth * time.cellCycleLength / log(2) * exp(log(2) * i / time.cellCycleLength), ...
                sum(mass.cellDry) + (...
                - rna.counts(rna.matureIndexs(rna.matureMRNAIndexs), comp.cytosolIndexs) ...
                + initRNAs(rna.matureIndexs(rna.matureMRNAIndexs)) * exp(log(2) * i / time.cellCycleLength) ...
                )' * rna.molecularWeights(rna.matureIndexs(rna.matureMRNAIndexs)) / ConstantUtil.nAvogadro, ...
                'relative', 0.3, 0);
            
            %% assertions -- RNA
            assertIn(corr(rna.expression(rna.matureIndexs), simMatureRNAExp), [0.95 1]);
            assertIn(180 / pi * acos(rna.expression(rna.matureIndexs)' * simMatureRNAExp / ...
                (sqrt(rna.expression(rna.matureIndexs)' * rna.expression(rna.matureIndexs)) * ...
                sqrt(simMatureRNAExp' * simMatureRNAExp))), ...
                [0 15]);
            
            %% assertions -- mRNA
            assertElementsAlmostEqual(...
                rna.molecularWeights(rna.matureIndexs(rna.matureMRNAIndexs))' * initRNAs(rna.matureIndexs(rna.matureMRNAIndexs)) * exp(i * log(2) / time.cellCycleLength), ...
                sum(rna.counts(rna.matureIndexs(rna.matureMRNAIndexs), :), 2)' * ...
                rna.molecularWeights(rna.matureIndexs(rna.matureMRNAIndexs)),...
                'relative', 1e-1);
            
            assertIn(corr(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs)), simMatureRNAExp(rna.matureMRNAIndexs)), [0.75 1])
            assertIn(180 / pi * acos(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))' * simMatureRNAExp(rna.matureMRNAIndexs) / ...
                (sqrt(rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))' * rna.expression(rna.matureIndexs(rna.matureMRNAIndexs))) * ...
                sqrt(simMatureRNAExp(rna.matureMRNAIndexs)' * simMatureRNAExp(rna.matureMRNAIndexs)))), ...
                [0 30]);
            
            %% assertions -- tRNA
            assertIn(min(tRNAs(:)), [4 Inf]);
            assertElementsAlmostEqual(...
                rna.molecularWeights(rna.matureIndexs(rna.matureTRNAIndexs))' * initRNAs(rna.matureIndexs(rna.matureTRNAIndexs)) * exp(i * log(2) / time.cellCycleLength), ...
                sum(rna.counts(rna.matureIndexs(rna.matureTRNAIndexs), :) + rna.counts(rna.aminoacylatedIndexs(rna.matureTRNAIndexs), :), 2)' * ...
                rna.molecularWeights(rna.matureIndexs(rna.matureTRNAIndexs)),...
                'relative', 5e-2);
            
            %% assertions -- proteins
            expProd = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * simMatureRNAExp;
            
            newMonomers = pm.counts - initMonomers;
            newComplexs = pc.counts - initComplexs;
            
            totNewMonomers = sum(...
                + newMonomers(pm.nascentIndexs, :) ...
                + newMonomers(pm.processedIIndexs, :) ...
                + newMonomers(pm.processedIIIndexs, :) ...
                + newMonomers(pm.foldedIndexs, :) ...
                + newMonomers(pm.matureIndexs, :) ...
                + newMonomers(pm.boundIndexs, :) ...
                + newMonomers(pm.misfoldedIndexs, :) ...
                + newMonomers(pm.damagedIndexs, :) ...
                + newMonomers(pm.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.nascentIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.matureIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.inactivatedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.boundIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.misfoldedIndexs, :) ...
                + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * newComplexs(pc.damagedIndexs, :), ...
                2);
            
            %amino acids are mostly used -- not too many are free
            assertIn(min(min(aaCnts(:, end-200+1:end))), [0 2 * min(initMetabolites(m.aminoAcidIndexs(1:20), comp.cytosolIndexs))]);
            assertIn(...
                (m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)' * m.molecularWeights(m.aminoAcidIndexs)) / ...
                (sum(pm.dryWeight + pc.dryWeight + pol.dryWeight) * ConstantUtil.nAvogadro), ...
                [0 0.02]);
            
            %mature/bound monomer counts increased
            assertIn(min(totNewMonomers), [-100 Inf]);
            if ~all(expProd)
                assertIn(max(totNewMonomers(expProd == 0)), [-Inf 0]);
            end
            assertElementsAlmostEqual(...
                totNewMonomers' * monMWs / ConstantUtil.nAvogadro, ...
                growth * mass.cellInitialDryWeight * mass.dryWeightFractionProtein, ...
                'relative', 10e-2, 0);
            assertElementsAlmostEqual(...
                (m.counts(m.aminoAcidIndexs, comp.cytosolIndexs)' * m.molecularWeights(m.aminoAcidIndexs) + ...
                totNewMonomers' * monMWs) / ConstantUtil.nAvogadro, ...
                growth * mass.cellInitialDryWeight * mass.dryWeightFractionProtein, ...
                'relative', 10e-2, 0);
            
            %proteins progress through entire synthesis pathway
            assertIn(sum(sum(pm.counts(pm.nascentIndexs, :))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.processedIIndexs, :))), [0 20]);
            assertIn(sum(sum(pm.counts(pm.processedIIIndexs, :))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.signalSequenceIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.foldedIndexs, :))), [0 50]);
            assertIn(sum(sum(pm.counts(pm.misfoldedIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pm.counts(pm.damagedIndexs, :))), [0 50]);
            
            assertIn(sum(sum(pc.counts(pc.nascentIndexs, :))), [0 10]);
            assertIn(sum(sum(pc.counts(pc.misfoldedIndexs, setdiff(1:end, comp.extracellularIndexs)))), [0 10]);
            assertIn(sum(sum(pc.counts(pc.damagedIndexs, :))), [0 10]);
            
            assertIn(numel(pol.abortedSequences), [0 2]);
            
            %RNA production / expression matches expectations
            assertIn(corr(expProd, totNewMonomers), [0.90 1]);
            assertIn(180 / pi * acos((expProd' * totNewMonomers) / (sqrt(expProd' * expProd) * sqrt(totNewMonomers' * totNewMonomers))), ...
                [0 20]);
            
            %protein decay
            assertElementsAlmostEqual((...
                + pm.lengths' * max(0, (pm.decayRates .* sum(initMonomers, 2))) ...
                + pm.lengths(pm.matureIndexs)' * sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * ...
                max(0, pc.decayRates(pc.matureIndexs) .* sum(initComplexs(pc.matureIndexs, :) + initComplexs(pc.boundIndexs, :), 2)) ...
                ) * time.cellCycleLength / log(2) * (exp(log(2) * iterMax / time.cellCycleLength)-1), ...
                sum(freedAAs(:)), ...
                'relative', 2.5e-1, 0);
        end
        
        function testSynthesisAndDecay_FewModules(this)
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'ProteinDecay'
                'tRNAAminoacylation'
                'Translation'});
            this.simulation = sim;
            
            %references
            trl = sim.process('Translation');
            ta = sim.process('tRNAAminoacylation');
            comp = trl.compartment;
            g = trl.gene;
            r = trl.ribosome;
            p = trl.polypeptide;
            met = trl.metabolite;
            rna = trl.rna;
            pm = trl.monomer;
            pc = trl.complex;
            
            %initialize at beginning of cell cycle
            nRibs = ...
                + trl.enzymes(trl.enzymeIndexs_ribosome70S) ...
                + min(trl.enzymes(trl.enzymeIndexs_ribosome30S), trl.boundEnzymes(trl.enzymeIndexs_ribosome50S)) ...
                + trl.boundEnzymes(trl.enzymeIndexs_ribosome70S);
            
            trl.substrates(:) = 1e12;
            ta.substrates(:) = 1e12;
            ta.substrates(ta.substrateIndexs_aminoAcids) = 0;
            ta.substrates(ta.substrateIndexs_fmethionine) = 0;
            aaProd = log(2) / (3600 * 8.5) * (...
                + met.biomassProduction(met.aminoAcidIndexs, comp.cytosolIndexs) ...
                - met.byproducts(met.aminoAcidIndexs, comp.cytosolIndexs));
            
            tmpProd = pm.baseCounts(pm.nascentIndexs, met.aminoAcidIndexs)' * ...
                rna.matureRNAGeneComposition(g.mRNAIndexs, :) * ...
                rna.expression(rna.matureIndexs);
            tmpProd(13) = tmpProd(13) + tmpProd(21);
            tmpProd(21) = 0;
            tmpProd = 397.8727 * tmpProd / sum(tmpProd);
            assertElementsAlmostEqual(aaProd, tmpProd, 'relative', 0.35, 1); %loose because amino acids are used for multiple processes and some are free
            assertIn(sum(tmpProd) / (nRibs * trl.ribosomeElongationRate), [0 0.9]);
            
            trl.enzymes = trl.enzymes + trl.boundEnzymes;
            trl.boundEnzymes(:) = 0;
            trl.enzymes(trl.enzymeIndexs_initiationFactor3) = ...
                trl.enzymes(trl.enzymeIndexs_initiationFactor3) + ...
                trl.enzymes(trl.enzymeIndexs_ribosome30SIF3);
            trl.enzymes(trl.enzymeIndexs_ribosome30S) = ...
                trl.enzymes(trl.enzymeIndexs_ribosome30S) + ...
                trl.enzymes(trl.enzymeIndexs_ribosome30SIF3) + ...
                trl.enzymes(trl.enzymeIndexs_ribosome70S);
            trl.enzymes(trl.enzymeIndexs_ribosome50S) = ...
                trl.enzymes(trl.enzymeIndexs_ribosome50S) + ...
                trl.enzymes(trl.enzymeIndexs_ribosome70S);
            trl.enzymes(trl.enzymeIndexs_ribosome30SIF3) = 0;
            trl.enzymes(trl.enzymeIndexs_ribosome70S) = 0;
            
            trl.monomers(:) = 0;
            totEnzymes = trl.enzymes;
            
            r.states = repmat(r.notExistValue, nRibs, 1);
            r.boundMRNAs = zeros(nRibs, 1);
            r.mRNAPositions = zeros(nRibs, 1);
            r.tmRNAPositions = zeros(nRibs, 1);
            p.boundMRNAs = r.boundMRNAs;
            p.nascentMonomerLengths = r.mRNAPositions;
            p.proteolysisTagLengths = r.tmRNAPositions;
            
            ta.freeRNAs = ta.freeRNAs + ta.aminoacylatedRNAs;
            ta.aminoacylatedRNAs(:) = 0;
            
            rnaExp = rna.expression;
            rnaMWs = rna.molecularWeights;
            initRNAs = rnaExp * sum(rna.dryWeight) / (rnaExp' * rnaMWs / edu.stanford.covert.util.ConstantUtil.nAvogadro);
            rna.counts(rna.matureIndexs(rna.matureMRNAIndexs), comp.cytosolIndexs) = ...
                round(1e3 * initRNAs(rna.matureIndexs(rna.matureMRNAIndexs))); %lots of RNA so that there isn't much variance in amino acid usage
            
            trl.copyToState();
            ta.copyToState();
            
            ta.copyFromState();
            ta.initializeState();
            ta.copyToState();
            trl.copyFromState();
            trl.initializeState();
            trl.copyToState();
            
            assertAllEqual(0, trl.monomers);
            
            iRib = find(r.states == r.activeValue);
            initNonFreeAACounts = zeros(21, 1);
            for i = 1:numel(iRib)
                seq = p.monomerAASequences{p.boundMRNAs(iRib(i))}(1:p.nascentMonomerLengths(iRib(i)));
                initNonFreeAACounts = initNonFreeAACounts + ProteinMonomer.computeBaseCount(seq, 21, 1:21, true)';
            end
            initNonFreeAACounts(1:end-1) = initNonFreeAACounts(1:end-1) + ...
                max(0, -ta.reactionStoichiometryMatrix(ta.substrateIndexs_aminoAcids, :) * ...
                ta.reactionModificationMatrix * rna.counts(rna.aminoacylatedIndexs(ta.aminoacylatedRNAGlobalIndexs), comp.cytosolIndexs));
            
            pm.decayRates(setdiff(1:end, pm.nascentIndexs)) = 0;
            pc.decayRates(:) = 0;
            trl.tmRNABindingProbability = 0;
            
            %evolve
            this.assertStateIsValid();
            iterMax = 2000;
            cumAAProd = zeros(size(aaProd));
            simMrnaExp = zeros(size(trl.mRNAs));
            for i = 1:iterMax
                tmpAAProd = trl.randStream.stochasticRound(aaProd);
                cumAAProd = cumAAProd + tmpAAProd;
                met.counts(met.aminoAcidIndexs, 1) = met.counts(met.aminoAcidIndexs, 1) + tmpAAProd;
                
                sim.evolveState();
                
                this.assertStateIsValid();
                
                simMrnaExp = simMrnaExp + trl.mRNAs;
            end
            
            ta.copyFromState();
            trl.copyFromState();
            
            %assert mRNA expression and synthesis are 1:1 correlated
            assertIn(corr(simMrnaExp, trl.monomers), [0.80 1]);
            assertIn(180 / pi * acos((simMrnaExp' * trl.monomers) / (sqrt(simMrnaExp' * simMrnaExp) * sqrt(trl.monomers' * trl.monomers))), ...
                [0 30]);
            
            %assert AAs correctly accounted
            iRib = find(r.states == r.activeValue);
            nonFreeAACounts = zeros(21, 1);
            for i = 1:numel(iRib)
                seq = p.monomerAASequences{p.boundMRNAs(iRib(i))}(1:p.nascentMonomerLengths(iRib(i)));
                nonFreeAACounts = nonFreeAACounts + ProteinMonomer.computeBaseCount(seq, 21, 1:21, true)';
            end
            nonFreeAACounts(1:end-1) = nonFreeAACounts(1:end-1) + ...
                max(0, -ta.reactionStoichiometryMatrix(ta.substrateIndexs_aminoAcids, :) * ...
                ta.reactionModificationMatrix * rna.counts(rna.aminoacylatedIndexs(ta.aminoacylatedRNAGlobalIndexs), comp.cytosolIndexs));
            
            usedAAs1 = cumAAProd - met.counts(met.aminoAcidIndexs, comp.cytosolIndexs);
            usedAAs1(13) = usedAAs1(13) + usedAAs1(21);
            usedAAs1(21) = 0;
            usedAAs2 = pm.baseCounts(pm.nascentIndexs, met.aminoAcidIndexs)' * sum(pm.counts(pm.nascentIndexs, :), 2) + ...
                nonFreeAACounts - initNonFreeAACounts;
            usedAAs2(13) = usedAAs2(13) + usedAAs2(21);
            usedAAs2(21) = 0;
            assertEqual(usedAAs1', usedAAs2');
            
            %assert AAs used
            assertEqual(0, min(ta.substrates(ta.substrateIndexs_aminoAcids)));
            assertIn(sum(ta.substrates(ta.substrateIndexs_aminoAcids)) / sum(cumAAProd), [0 0.10]);
            assertIn((ta.substrates(ta.substrateIndexs_aminoAcids)' * ta.substrateMolecularWeights(ta.substrateIndexs_aminoAcids)) / ...
                (sum(pm.dryWeight + pc.dryWeight)*edu.stanford.covert.util.ConstantUtil.nAvogadro), ...
                [0 0.01]);
        end
    end
    
    %RNA misfolding and polymerase decay
    methods
        function testDecayRNAPolymerases(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            
            pc.decayRates(:) = 0;
            pc.decayRates(pc.matureIndexs(pc.rnaPolymeraseIndexs)) = Inf;
            pc.decayRates(pc.boundIndexs(pc.rnaPolymeraseIndexs)) = Inf;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            nAbortedSequences = sum(t.boundTranscriptProgress > 1);
            
            %simulate protein decay
            for i = 1:10
                sim.evolveState();
            end
            
            %assert protein, chromosome, rna polymerase, transcript states
            %updated correctly
            assertTrue(all(pm.counts(:) >= 0));
            assertTrue(all(pc.counts(:) >= 0));
            assertAllEqual(0, pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(0, pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :));
            assertEqual(0, nnz(c.complexBoundSites == pc.rnaPolymeraseIndexs(1)));
            assertEqual(0, nnz(c.complexBoundSites == pc.rnaPolymeraseIndexs(2)));
            assertAllEqual(r.notExistValue, r.states);
            assertAllEqual(0, r.positionStrands);
            assertAllEqual(t.nullTranscriptValue, t.boundTranscriptProgress);
            assertAllEqual(t.nullTranscriptValue, t.boundTranscriptChromosome);
            assertAllEqual(t.nullTranscriptValue, t.boundTranscriptionUnits);
            assertEqual(nAbortedSequences, numel(t.abortedSequences));
        end
        
        function testATRNAPolymeraseDecay(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            pc = sim.state('ProteinComplex');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            
            pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs) = 1;
            
            iTU = find(t.transcriptionUnitDirections == 1 & t.transcriptionUnitLengths > 1000, 1, 'first');
            r.states = r.activelyTranscribingValue + 1;
            r.positionStrands = [t.transcriptionUnitFivePrimeCoordinates(iTU) + 1  1];
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptChromosome = 1;
            t.boundTranscriptionUnits = iTU;
            
            c.initialize();
            c.setSiteProteinBound(r.positionStrands, 1, 1, [], pc.rnaPolymeraseIndexs(1), [], [], true, true, 1, false, []);
            
            pc.decayRates(:) = 0;
            pc.decayRates(pc.boundIndexs(pc.rnaPolymeraseIndexs)) = Inf;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            %simulate protein decay
            sim.evolveState();
            
            %assert protein, chromosome, rna polymerase, transcript states
            %updated correctly
            assertAllEqual(0, pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(0, pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :));
            assertEqual(0, nnz(c.complexBoundSites));
            assertEqual(r.notExistValue, r.states);
            assertEqual([0 0], r.positionStrands);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptProgress);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptChromosome);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptionUnits);
            assertEqual(1, numel(t.abortedSequences));
            assertEqual(1, numel(t.abortedSequences{1}));
        end
        
        function testInitiatingRNAPolymeraseDecay(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            pc = sim.state('ProteinComplex');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            
            pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs) = 1;
            
            iTU = find(t.transcriptionUnitDirections == 1 & t.transcriptionUnitLengths > 1000, 1, 'first');
            r.states = r.activelyTranscribingValue;
            r.positionStrands = [t.transcriptionUnitFivePrimeCoordinates(iTU)  1];
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptChromosome = 1;
            t.boundTranscriptionUnits = iTU;
            
            c.initialize();
            c.setSiteProteinBound(r.positionStrands, 1, 1, [], pc.rnaPolymeraseIndexs(1), [], [], true, true, 1, false, []);
            
            pc.decayRates(:) = 0;
            pc.decayRates(pc.boundIndexs(pc.rnaPolymeraseIndexs)) = Inf;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            %simulate protein decay
            sim.evolveState();
            
            %assert protein, chromosome, rna polymerase, transcript states
            %updated correctly
            assertAllEqual(0, pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(0, pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :));
            assertEqual(0, nnz(c.complexBoundSites));
            assertEqual(r.notExistValue, r.states);
            assertEqual([0 0], r.positionStrands);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptProgress);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptChromosome);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptionUnits);
            assertEqual(0, numel(t.abortedSequences));
        end
        
        function testSBRNAPolymeraseDecay(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            pc = sim.state('ProteinComplex');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            
            pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs) = 1;
            
            iTU = find(t.transcriptionUnitDirections == 1 & t.transcriptionUnitLengths > 1000, 1, 'first');
            r.states = r.specificallyBoundValue;
            r.positionStrands = [t.transcriptionUnitFivePrimeCoordinates(iTU)  1];
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptChromosome = 1;
            t.boundTranscriptionUnits = iTU;
            
            c.initialize();
            c.setSiteProteinBound(r.positionStrands, 1, 1, [], pc.rnaPolymeraseIndexs(1), [], [], true, true, 1, false, []);
            
            pc.decayRates(:) = 0;
            pc.decayRates(pc.boundIndexs(pc.rnaPolymeraseIndexs)) = Inf;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            %simulate protein decay
            sim.evolveState();
            
            %assert protein, chromosome, rna polymerase, transcript states
            %updated correctly
            assertAllEqual(0, pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(0, pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :));
            assertEqual(0, nnz(c.complexBoundSites));
            assertEqual(r.notExistValue, r.states);
            assertEqual([0 0], r.positionStrands);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptProgress);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptChromosome);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptionUnits);
            assertEqual(0, numel(t.abortedSequences));
        end
        
        function testNSBRNAPolymeraseDecay(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            pc = sim.state('ProteinComplex');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            
            pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs) = 1;
            
            r.states = r.nonSpecificallyBoundValue;
            r.positionStrands = [1000 1];
            t.boundTranscriptProgress = t.nullTranscriptValue;
            t.boundTranscriptChromosome = t.nullTranscriptValue;
            t.boundTranscriptionUnits = t.nullTranscriptValue;
            
            c.initialize();
            c.setSiteProteinBound(r.positionStrands, 1, 1, [], pc.rnaPolymeraseIndexs(1), [], [], true, true, 1, false, []);
            
            pc.decayRates(:) = 0;
            pc.decayRates(pc.boundIndexs(pc.rnaPolymeraseIndexs)) = Inf;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            %simulate protein decay
            sim.evolveState();
            
            %assert protein, chromosome, rna polymerase, transcript states
            %updated correctly
            assertAllEqual(0, pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(0, pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :));
            assertEqual(0, nnz(c.complexBoundSites));
            assertEqual(r.notExistValue, r.states);
            assertEqual([0 0], r.positionStrands);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptProgress);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptChromosome);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptionUnits);
            assertEqual(0, numel(t.abortedSequences));
        end
        
        function testFreeRNAPolymeraseDecay(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            pc = sim.state('ProteinComplex');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            
            pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs) = 1;
            
            r.states = r.freeValue;
            r.positionStrands = [0 0];
            t.boundTranscriptProgress = t.nullTranscriptValue;
            t.boundTranscriptChromosome = t.nullTranscriptValue;
            t.boundTranscriptionUnits = t.nullTranscriptValue;
            
            c.initialize();
            
            pc.decayRates(:) = 0;
            pc.decayRates(pc.matureIndexs(pc.rnaPolymeraseIndexs)) = Inf;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            %simulate protein decay
            sim.evolveState();
            
            %assert protein, chromosome, rna polymerase, transcript states
            %updated correctly
            assertAllEqual(0, pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(0, pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :));
            assertEqual(0, nnz(c.complexBoundSites));
            assertEqual(r.notExistValue, r.states);
            assertEqual([0 0], r.positionStrands);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptProgress);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptChromosome);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptionUnits);
            assertEqual(0, numel(t.abortedSequences));
        end
        
        function testMisfoldingRNAPolymerases(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            dcy = sim.process('ProteinDecay');
            
            pm.counts(pm.boundIndexs, :) = 0;
            pc.counts(pc.matureIndexs(setdiff(1:end, pc.rnaPolymeraseIndexs)), :) = 0;
            pc.counts(pc.boundIndexs(setdiff(1:end, pc.rnaPolymeraseIndexs)), :) = 0;
            
            c.monomerBoundSites(find(c.monomerBoundSites)) = 0; %#ok<FNDSB>
            [subs, vals] = find(c.complexBoundSites);
            c.complexBoundSites(subs(~ismember(vals, pc.rnaPolymeraseIndexs), :)) = 0;
            
            dcy.proteinMisfoldingRate = 1;
            pc.counts(pc.matureIndexs(dcy.enzymeGlobalIndexs(dcy.enzymeIndexs_clpBProtease)), :) = 0;
            pc.decayRates(:) = 0;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            nRNAPols = [
                sum(r.states == r.freeValue | r.states == r.nonSpecificallyBoundValue | r.states > r.activelyTranscribingValue)
                sum(r.states == r.specificallyBoundValue | r.states == r.activelyTranscribingValue)];
            nAbortedSequences = sum(t.boundTranscriptProgress > 1);
            
            %simulate protein decay
            for i = 1:10
                sim.evolveState();
            end
            
            %assert protein, chromosome, rna polymerase, transcript states
            %updated correctly
            assertTrue(all(pm.counts(:) >= 0));
            assertTrue(all(pc.counts(:) >= 0));
            assertAllEqual(0, pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(0, pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :));
            assertEqual(nRNAPols, pc.counts(pc.misfoldedIndexs(pc.rnaPolymeraseIndexs), comp.cytosolIndexs));
            assertEqual(0, nnz(c.complexBoundSites == pc.rnaPolymeraseIndexs(1)));
            assertEqual(0, nnz(c.complexBoundSites == pc.rnaPolymeraseIndexs(2)));
            assertAllEqual(r.notExistValue, r.states);
            assertAllEqual(0, r.positionStrands);
            assertAllEqual(t.nullTranscriptValue, t.boundTranscriptProgress);
            assertAllEqual(t.nullTranscriptValue, t.boundTranscriptChromosome);
            assertAllEqual(t.nullTranscriptValue, t.boundTranscriptionUnits);
            assertEqual(nAbortedSequences, numel(t.abortedSequences));
        end
        
        function testATRNAPolymeraseMisfolding(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            rib = sim.state('Ribosome');
            pol = sim.state('Polypeptide');
            dcy = sim.process('ProteinDecay');
            
            pm.counts(pm.boundIndexs, :) = 0;
            pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs, :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs) = 1;
            rib.states(:) = rib.notExistValue;
            rib.boundMRNAs(:) = 0;
            rib.mRNAPositions(:) = 0;
            rib.tmRNAPositions(:) = 0;
            pol.boundMRNAs = rib.boundMRNAs;
            pol.nascentMonomerLengths = rib.mRNAPositions;
            pol.proteolysisTagLengths = rib.tmRNAPositions;
            
            iTU = find(t.transcriptionUnitDirections == 1 & t.transcriptionUnitLengths > 1000, 1, 'first');
            r.states = r.activelyTranscribingValue + 1;
            r.positionStrands = [t.transcriptionUnitFivePrimeCoordinates(iTU) + 1  1];
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptChromosome = 1;
            t.boundTranscriptionUnits = iTU;
            
            c.initialize();
            c.setSiteProteinBound(r.positionStrands, 1, 1, [], pc.rnaPolymeraseIndexs(1), [], [], true, true, 1, false, []);
            
            dcy.proteinMisfoldingRate = 1;
            pc.counts(pc.matureIndexs(dcy.enzymeGlobalIndexs(dcy.enzymeIndexs_clpBProtease)), :) = 0;
            pc.decayRates(:) = 0;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            %simulate protein decay
            sim.evolveState();
            
            %assert protein, chromosome, rna polymerase, transcript states
            %updated correctly
            assertAllEqual(0, pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(0, pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(1, pc.counts(pc.misfoldedIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs));
            assertEqual(0, nnz(c.complexBoundSites));
            assertEqual(r.notExistValue, r.states);
            assertEqual([0 0], r.positionStrands);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptProgress);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptChromosome);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptionUnits);
            assertEqual(1, numel(t.abortedSequences));
            assertEqual(1, numel(t.abortedSequences{1}));
        end
        
        function testInitiatingRNAPolymeraseMisfolding(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            dcy = sim.process('ProteinDecay');
            
            pm.counts(pm.boundIndexs, :) = 0;
            pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs, :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs) = 1;
            
            iTU = find(t.transcriptionUnitDirections == 1 & t.transcriptionUnitLengths > 1000, 1, 'first');
            r.states = r.activelyTranscribingValue;
            r.positionStrands = [t.transcriptionUnitFivePrimeCoordinates(iTU)  1];
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptChromosome = 1;
            t.boundTranscriptionUnits = iTU;
            
            c.initialize();
            c.setSiteProteinBound(r.positionStrands, 1, 1, [], pc.rnaPolymeraseIndexs(1), [], [], true, true, 1, false, []);
            
            dcy.proteinMisfoldingRate = 1;
            pc.counts(pc.matureIndexs(dcy.enzymeGlobalIndexs(dcy.enzymeIndexs_clpBProtease)), :) = 0;
            pc.decayRates(:) = 0;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            %simulate protein decay
            sim.evolveState();
            
            %assert protein, chromosome, rna polymerase, transcript states
            %updated correctly
            assertAllEqual(0, pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(0, pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(1, pc.counts(pc.misfoldedIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs));
            assertEqual(0, nnz(c.complexBoundSites));
            assertEqual(r.notExistValue, r.states);
            assertEqual([0 0], r.positionStrands);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptProgress);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptChromosome);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptionUnits);
            assertEqual(0, numel(t.abortedSequences));
        end
        
        function testSBRNAPolymeraseMisfolding(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            dcy = sim.process('ProteinDecay');
            
            pm.counts(pm.boundIndexs, :) = 0;
            pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs, :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs) = 1;
            
            iTU = find(t.transcriptionUnitDirections == 1 & t.transcriptionUnitLengths > 1000, 1, 'first');
            r.states = r.specificallyBoundValue;
            r.positionStrands = [t.transcriptionUnitFivePrimeCoordinates(iTU)  1];
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptChromosome = 1;
            t.boundTranscriptionUnits = iTU;
            
            c.initialize();
            c.setSiteProteinBound(r.positionStrands, 1, 1, [], pc.rnaPolymeraseIndexs(1), [], [], true, true, 1, false, []);
            
            dcy.proteinMisfoldingRate = 1;
            pc.counts(pc.matureIndexs(dcy.enzymeGlobalIndexs(dcy.enzymeIndexs_clpBProtease)), :) = 0;
            pc.decayRates(:) = 0;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            %simulate protein decay
            sim.evolveState();
            
            %assert protein, chromosome, rna polymerase, transcript states
            %updated correctly
            assertAllEqual(0, pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(0, pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(1, pc.counts(pc.misfoldedIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs));
            assertEqual(0, nnz(c.complexBoundSites));
            assertEqual(r.notExistValue, r.states);
            assertEqual([0 0], r.positionStrands);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptProgress);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptChromosome);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptionUnits);
            assertEqual(0, numel(t.abortedSequences));
        end
        
        function testNSBRNAPolymeraseMisfolding(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            dcy = sim.process('ProteinDecay');
            
            pm.counts(pm.boundIndexs, :) = 0;
            pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs, :) = 0;
            pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs) = 1;
            
            r.states = r.nonSpecificallyBoundValue;
            r.positionStrands = [1000 1];
            t.boundTranscriptProgress = t.nullTranscriptValue;
            t.boundTranscriptChromosome = t.nullTranscriptValue;
            t.boundTranscriptionUnits = t.nullTranscriptValue;
            
            c.initialize();
            c.setSiteProteinBound(r.positionStrands, 1, 1, [], pc.rnaPolymeraseIndexs(1), [], [], true, true, 1, false, []);
            
            dcy.proteinMisfoldingRate = 1;
            pc.counts(pc.matureIndexs(dcy.enzymeGlobalIndexs(dcy.enzymeIndexs_clpBProtease)), :) = 0;
            pc.decayRates(:) = 0;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            %simulate protein decay
            sim.evolveState();
            
            %assert protein, chromosome, rna polymerase, transcript states
            %updated correctly
            assertAllEqual(0, pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(0, pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(1, pc.counts(pc.misfoldedIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs));
            assertEqual(0, nnz(c.complexBoundSites));
            assertEqual(r.notExistValue, r.states);
            assertEqual([0 0], r.positionStrands);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptProgress);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptChromosome);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptionUnits);
            assertEqual(0, numel(t.abortedSequences));
        end
        
        function testFreeRNAPolymeraseMisfolding(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            dcy = sim.process('ProteinDecay');
            
            pm.counts(pm.boundIndexs, :) = 0;
            pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :) = 0;
            pc.counts(pc.boundIndexs, :) = 0;
            pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs) = 1;
            
            r.states = r.freeValue;
            r.positionStrands = [0 0];
            t.boundTranscriptProgress = t.nullTranscriptValue;
            t.boundTranscriptChromosome = t.nullTranscriptValue;
            t.boundTranscriptionUnits = t.nullTranscriptValue;
            
            c.initialize();
            
            dcy.proteinMisfoldingRate = 1;
            pc.counts(pc.matureIndexs(dcy.enzymeGlobalIndexs(dcy.enzymeIndexs_clpBProtease)), :) = 0;
            pc.decayRates(:) = 0;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            %simulate protein decay
            sim.evolveState();
            
            %assert protein, chromosome, rna polymerase, transcript states
            %updated correctly
            assertAllEqual(0, pc.counts(pc.matureIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(0, pc.counts(pc.boundIndexs(pc.rnaPolymeraseIndexs), :));
            assertAllEqual(1, pc.counts(pc.misfoldedIndexs(pc.rnaPolymeraseIndexs(1)), comp.cytosolIndexs));
            assertEqual(0, nnz(c.complexBoundSites));
            assertEqual(r.notExistValue, r.states);
            assertEqual([0 0], r.positionStrands);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptProgress);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptChromosome);
            assertEqual(t.nullTranscriptValue, t.boundTranscriptionUnits);
            assertEqual(0, numel(t.abortedSequences));
        end
    end
    
    methods
        function assertStateIsValid(this)
            sim = this.simulation;
            comp = sim.compartment;
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            rib = sim.state('Ribosome');
            pol = sim.state('Polypeptide');
            
            assertEqual(rib.boundMRNAs, pol.boundMRNAs);
            assertEqual(rib.mRNAPositions, pol.nascentMonomerLengths);
            assertEqual(rib.tmRNAPositions, pol.proteolysisTagLengths);
            
            if any(rib.states == 0)
                assertAllEqual(0, rib.boundMRNAs(rib.states == 0));
                assertAllEqual(0, rib.mRNAPositions(rib.states == 0));
                assertAllEqual(0, rib.tmRNAPositions(rib.states == 0));
            end
            assertTrue(all(rib.mRNAPositions(rib.states ~= 0) <= pol.monomerLengths(rib.boundMRNAs(rib.states ~= 0))));
            assertTrue(all(rib.mRNAPositions >= 0));
            if any(rib.states == rib.activeValue)
                assertAllEqual(0, rib.tmRNAPositions(rib.states == rib.activeValue));
            end
            assertTrue(all(rib.tmRNAPositions >= 0));
            assertTrue(all(rib.tmRNAPositions <= pol.proteolysisTagLength));
            
            assertEqual(pc.counts(pc.boundIndexs(pc.ribosome70SIndexs), comp.cytosolIndexs), sum(rib.states ~= 0));
            assertEqual(rna.counts(rna.boundIndexs(rna.matureTMRNAIndexs), comp.cytosolIndexs), sum(rib.states == rib.stalledValue));
            
            assertTrue(all(rna.counts(:) >= 0));
            assertTrue(all(pm.counts(:) >= 0));
            assertTrue(all(pc.counts(:) >= 0));
        end
        
        function seedSimulation(this, seed)
            s = this.simulation;
            
            if s.verbosity >= 1
                fprintf('Seed %d\n', seed);
            end
            
            s.applyOptions('seed', seed);
            s.seedRandStream();
            for i = 1:numel(s.states)
                o = s.states{i};
                o.seed = seed;
                o.seedRandStream();
            end
            for i = 1:numel(s.processes)
                o = s.processes{i};
                o.seed = seed;
                o.seedRandStream();
            end
        end
    end
end
