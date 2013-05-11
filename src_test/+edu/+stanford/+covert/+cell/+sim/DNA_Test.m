%DNA medium test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 4/19/2011
classdef DNA_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = DNA_Test(name)
            this = this@TestCase(name);
        end
        
        function setUp(this)
            sim = edu.stanford.covert.cell.sim.SimulationFixture.load([], {
                'ChromosomeCondensation'
                'ChromosomeSegregation'
                'DNADamage'
                'DNARepair'
                'DNASupercoiling'
                'Replication'
                'ReplicationInitiation'
                'Transcription'
                'TranscriptionalRegulation'
                });
            sim.applyOptions('verbosity', 0);
            
            this.simulation = sim;
        end
    end
    
    methods
        function testProteinDisplacementReactions(this)
            sim = this.simulation;
            c = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            cc = sim.process('ChromosomeCondensation');
            dr = sim.process('DNARepair');
            ds = sim.process('DNASupercoiling');
            rp = sim.process('Replication');
            ri = sim.process('ReplicationInitiation');
            ts = sim.process('Transcription');
            tr = sim.process('TranscriptionalRegulation');
            
            %- DNA repair enzymes can release everything, including DNA
            %  polymerase
            bindableEnzymes = [
                cc.enzymeWholeCellModelIDs(cc.enzymeIndexs_SMC_ADP)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_formamidopyrimidineGlycosylase)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_uracilGlycosylase)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_apurinicEndonuclease)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_incisionComplex)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_helicase35)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_recombinationStrandExchange)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_hollidayJunctionEndonuclease)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_hollidayJunctionHelicase)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_exonuclease53)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_polymerase)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_ligase)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_DisA)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_RM_typeI)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_RM_typeII)
                ds.enzymeWholeCellModelIDs
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_2coreBetaClampGammaComplexPrimase)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_coreBetaClampGammaComplex)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_coreBetaClampPrimase)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_betaClamp)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_helicase)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_ligase)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_ssb8mer)
                ri.enzymeWholeCellModelIDs(ri.enzymeIndexs_DnaA_Nmer_ATP)
                ri.enzymeWholeCellModelIDs(ri.enzymeIndexs_DnaA_Nmer_ADP)
                ts.enzymeWholeCellModelIDs(ts.enzymeIndexs_rnaPolymerase)
                ts.enzymeWholeCellModelIDs(ts.enzymeIndexs_rnaPolymeraseHoloenzyme)
                tr.enzymeWholeCellModelIDs
                ];
            
            %only bindable proteins can displace and be displaced
            assertTrue(all(ismember(c.reactionBoundMonomer, pm.getIndexs(bindableEnzymes))));
            assertTrue(all(ismember(c.reactionBoundComplex, pc.getIndexs(bindableEnzymes))));
            idxs = setdiff(1:numel(pm.matureIndexs), pm.getIndexs(bindableEnzymes));            
            assertFalse(any(any(c.reactionMonomerCatalysisMatrix(:, idxs))), ...
                sprintf('The following reactions must be removed\n- %s', ...
                strjoin(sprintf('\n- '), c.reactionWholeCellModelIDs{any(c.reactionMonomerCatalysisMatrix(:, idxs), 2)})));
            idxs = setdiff(1:numel(pc.matureIndexs), pc.getIndexs(bindableEnzymes));
            assertFalse(any(any(c.reactionComplexCatalysisMatrix(:, idxs))), ...
                sprintf('The following reactions must be removed\n- %s', ...
                strjoin(sprintf('\n- '), c.reactionWholeCellModelIDs{any(c.reactionComplexCatalysisMatrix(:, idxs), 2)})));           
            
            allReleaseReactions = cell(0, 2);            
            
            %DNA polymerase can release everything, except DNA repair
            %enzymes
            expReleasedProteins = [
                cc.enzymeWholeCellModelIDs(cc.enzymeIndexs_SMC_ADP)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_RM_typeI)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_RM_typeII)
                ds.enzymeWholeCellModelIDs
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_ligase)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_ssb8mer)
                ri.enzymeWholeCellModelIDs(ri.enzymeIndexs_DnaA_Nmer_ATP)
                ri.enzymeWholeCellModelIDs(ri.enzymeIndexs_DnaA_Nmer_ADP)
                tr.enzymeWholeCellModelIDs
                ts.enzymeWholeCellModelIDs(ts.enzymeIndexs_rnaPolymerase)
                ts.enzymeWholeCellModelIDs(ts.enzymeIndexs_rnaPolymeraseHoloenzyme)
                ];
            this.assertCorrectDisplacementReactions(rp.enzymeWholeCellModelIDs{rp.enzymeIndexs_2coreBetaClampGammaComplexPrimase}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(rp.enzymeWholeCellModelIDs{rp.enzymeIndexs_coreBetaClampGammaComplex}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(rp.enzymeWholeCellModelIDs{rp.enzymeIndexs_coreBetaClampPrimase}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(rp.enzymeWholeCellModelIDs{rp.enzymeIndexs_helicase}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(rp.enzymeWholeCellModelIDs{rp.enzymeIndexs_betaClamp}, expReleasedProteins);
            allReleaseReactions = [
                allReleaseReactions
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_2coreBetaClampGammaComplexPrimase(ones(size(expReleasedProteins)))) expReleasedProteins
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_coreBetaClampGammaComplex(ones(size(expReleasedProteins)))) expReleasedProteins
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_coreBetaClampPrimase(ones(size(expReleasedProteins)))) expReleasedProteins
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_helicase(ones(size(expReleasedProteins)))) expReleasedProteins
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_betaClamp(ones(size(expReleasedProteins)))) expReleasedProteins
                ];
            
            %RNA polymerase can release everything, except DNA polymerase
            %and DNA repair enzymes and including DnaA 1mers
            expReleasedProteins = [
                cc.enzymeWholeCellModelIDs(cc.enzymeIndexs_SMC_ADP)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_RM_typeI)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_RM_typeII)
                ds.enzymeWholeCellModelIDs
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_ssb8mer)
                ri.enzymeWholeCellModelIDs(ri.enzymeIndexs_DnaA_1mer_ATP)
                ri.enzymeWholeCellModelIDs(ri.enzymeIndexs_DnaA_1mer_ADP)
                tr.enzymeWholeCellModelIDs
                ];
            this.assertCorrectDisplacementReactions(ts.enzymeWholeCellModelIDs{ts.enzymeIndexs_rnaPolymerase}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(ts.enzymeWholeCellModelIDs{ts.enzymeIndexs_rnaPolymeraseHoloenzyme}, expReleasedProteins);
            allReleaseReactions = [
                allReleaseReactions
                ts.enzymeWholeCellModelIDs(ts.enzymeIndexs_rnaPolymerase(ones(size(expReleasedProteins)))) expReleasedProteins
                ts.enzymeWholeCellModelIDs(ts.enzymeIndexs_rnaPolymeraseHoloenzyme(ones(size(expReleasedProteins)))) expReleasedProteins
                ];
            
            %DNA repair enzymes can release everything, including DNA
            %polymerase. All repair enzymes releases DisA
            expReleasedProteins = [
                cc.enzymeWholeCellModelIDs(cc.enzymeIndexs_SMC_ADP)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_DisA)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_RM_typeI)
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_RM_typeII)
                ds.enzymeWholeCellModelIDs
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_2coreBetaClampGammaComplexPrimase)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_coreBetaClampGammaComplex)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_coreBetaClampPrimase)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_betaClamp)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_helicase)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_ligase)
                rp.enzymeWholeCellModelIDs(rp.enzymeIndexs_ssb8mer)
                ri.enzymeWholeCellModelIDs(ri.enzymeIndexs_DnaA_Nmer_ATP)
                ri.enzymeWholeCellModelIDs(ri.enzymeIndexs_DnaA_Nmer_ADP)
                tr.enzymeWholeCellModelIDs
                ts.enzymeWholeCellModelIDs(ts.enzymeIndexs_rnaPolymerase)
                ts.enzymeWholeCellModelIDs(ts.enzymeIndexs_rnaPolymeraseHoloenzyme)
                ];
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_formamidopyrimidineGlycosylase}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_uracilGlycosylase}, expReleasedProteins);            
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_apurinicEndonuclease}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_incisionComplex}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_helicase35}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_recombinationStrandExchange}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_hollidayJunctionEndonuclease}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_hollidayJunctionHelicase}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_exonuclease53}, expReleasedProteins);
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_polymerase}, [expReleasedProteins; dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_recombinationStrandExchange)]);
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_ligase}, setdiff(expReleasedProteins, {dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_ligase}, ...
                'DNA_GYRASE', 'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE', ...
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX', 'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE', 'MG_001_DIMER', 'MG_094_HEXAMER'}));
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_DisA}, setdiff(expReleasedProteins, {dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_DisA}, ...
                'DNA_GYRASE', 'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE', ...
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX', 'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE', 'MG_001_DIMER', 'MG_094_HEXAMER'}));
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_RM_typeI}, setdiff(expReleasedProteins, {dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_RM_typeI}, ...
                'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE', ...
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX', 'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE', 'MG_001_DIMER', 'MG_094_HEXAMER'}));
            this.assertCorrectDisplacementReactions(dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_RM_typeII}, setdiff(expReleasedProteins, {dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_RM_typeII}, ...
                'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE', ...
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX', 'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE', 'MG_001_DIMER', 'MG_094_HEXAMER'}));
            allReleaseReactions = [
                allReleaseReactions
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_formamidopyrimidineGlycosylase(ones(size(expReleasedProteins)))) expReleasedProteins
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_uracilGlycosylase(ones(size(expReleasedProteins)))) expReleasedProteins
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_apurinicEndonuclease(ones(size(expReleasedProteins)))) expReleasedProteins
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_incisionComplex(ones(size(expReleasedProteins)))) expReleasedProteins
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_helicase35(ones(size(expReleasedProteins)))) expReleasedProteins
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_recombinationStrandExchange(ones(size(expReleasedProteins)))) expReleasedProteins
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_hollidayJunctionEndonuclease(ones(size(expReleasedProteins)))) expReleasedProteins
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_hollidayJunctionHelicase(ones(size(expReleasedProteins)))) expReleasedProteins
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_exonuclease53(ones(size(expReleasedProteins)))) expReleasedProteins
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_polymerase(ones(size(expReleasedProteins)+[1 0]))) [expReleasedProteins; dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_recombinationStrandExchange)]
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_ligase(ones(size(expReleasedProteins)-[7 0]))) setdiff(expReleasedProteins, {dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_ligase}, ...
                'DNA_GYRASE', 'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE', ...
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX', 'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE', 'MG_001_DIMER', 'MG_094_HEXAMER'}')
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_DisA(ones(size(expReleasedProteins)-[7 0]))) setdiff(expReleasedProteins, {dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_DisA}, ...
                'DNA_GYRASE', 'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE', ...
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX', 'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE', 'MG_001_DIMER', 'MG_094_HEXAMER'}')
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_RM_typeI(ones(size(expReleasedProteins)-[6 0]))) setdiff(expReleasedProteins, {dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_RM_typeI}, ...
                'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE', ...
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX', 'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE', 'MG_001_DIMER', 'MG_094_HEXAMER'}')
                dr.enzymeWholeCellModelIDs(dr.enzymeIndexs_RM_typeII(ones(size(expReleasedProteins)-[6 0]))) setdiff(expReleasedProteins, {dr.enzymeWholeCellModelIDs{dr.enzymeIndexs_RM_typeII}, ...
                'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE', ...
                'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX', 'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE', 'MG_001_DIMER', 'MG_094_HEXAMER'}')
                ];
              
            %transcription factors releases SMCs, DnaA 1mers, topos
            expReleasedProteins = [
                cc.enzymeWholeCellModelIDs(cc.enzymeIndexs_SMC_ADP)
                ri.enzymeWholeCellModelIDs(ri.enzymeIndexs_DnaA_1mer_ATP)
                ri.enzymeWholeCellModelIDs(ri.enzymeIndexs_DnaA_1mer_ADP)
                ds.enzymeWholeCellModelIDs
                ];
            for i = 1:numel(tr.enzymeWholeCellModelIDs)
                this.assertCorrectDisplacementReactions(tr.enzymeWholeCellModelIDs{i}, expReleasedProteins);
                allReleaseReactions = [
                    allReleaseReactions
                    tr.enzymeWholeCellModelIDs(i(ones(size(expReleasedProteins)))) expReleasedProteins
                    ]; %#ok<AGROW>
            end            
            
            %SSBs, DnaA release nothing
            expReleasedProteins = cell(0, 1);
            this.assertCorrectDisplacementReactions(rp.enzymeWholeCellModelIDs{rp.enzymeIndexs_ssb8mer}, expReleasedProteins);
            for i = 1: numel(ri.enzymeIndexs_DnaA_Nmer_ATP)
                this.assertCorrectDisplacementReactions(ri.enzymeWholeCellModelIDs{ri.enzymeIndexs_DnaA_Nmer_ATP(i)}, expReleasedProteins);
                this.assertCorrectDisplacementReactions(ri.enzymeWholeCellModelIDs{ri.enzymeIndexs_DnaA_Nmer_ADP(i)}, expReleasedProteins);
            end
            
            %SMCs, gyrase release DnaA 1mers
            expReleasedProteins = [
                ri.enzymeWholeCellModelIDs(ri.enzymeIndexs_DnaA_1mer_ATP)
                ri.enzymeWholeCellModelIDs(ri.enzymeIndexs_DnaA_1mer_ADP)
                ];
            this.assertCorrectDisplacementReactions(cc.enzymeWholeCellModelIDs{cc.enzymeIndexs_SMC_ADP}, expReleasedProteins);
            allReleaseReactions = [
                allReleaseReactions
                cc.enzymeWholeCellModelIDs(cc.enzymeIndexs_SMC_ADP(ones(size(expReleasedProteins)))) expReleasedProteins
                ];
            for i = 1:numel(ds.enzymeWholeCellModelIDs)
                this.assertCorrectDisplacementReactions(ds.enzymeWholeCellModelIDs{i}, expReleasedProteins);
                allReleaseReactions = [
                    allReleaseReactions
                    ds.enzymeWholeCellModelIDs(i(ones(size(expReleasedProteins)))) expReleasedProteins
                    ]; %#ok<AGROW>
            end
            
            %no duplicates
            [tmp, i] = unique([c.reactionBoundMonomer c.reactionBoundComplex c.reactionMonomerCatalysisMatrix c.reactionComplexCatalysisMatrix], 'rows');            
            assertEqual(size(c.reactionBoundMonomer, 1), size(tmp, 1), ...
                sprintf('Duplicate reactions must be removed\n- %s', strjoin(sprintf('\n- '), c.reactionWholeCellModelIDs{setdiff(1:end, i)}))); 
            
            %no self displacement
            tmp = (1:numel(c.reactionBoundMonomer))';
            
            tmpIdxs = find(c.reactionBoundMonomer~=0);
            tmp2 = c.reactionMonomerCatalysisMatrix(sub2ind(...
                size(c.reactionMonomerCatalysisMatrix), ...
                tmp(tmpIdxs), c.reactionBoundMonomer(tmpIdxs)));
            assertFalse(any(tmp2), sprintf('Reactions must be removed\n- %s', strjoin(sprintf('\n- '), c.reactionWholeCellModelIDs{tmpIdxs(tmp2~=0)})));
            
            tmpIdxs = find(c.reactionBoundComplex~=0);
            tmp2 = c.reactionComplexCatalysisMatrix(sub2ind(...
                size(c.reactionComplexCatalysisMatrix), ...
                tmp(tmpIdxs), c.reactionBoundComplex(tmpIdxs)));
            assertFalse(any(tmp2), sprintf('Reactions must be removed\n- %s', strjoin(sprintf('\n- '), c.reactionWholeCellModelIDs{tmpIdxs(tmp2~=0)})));
            
            %no extraneous reactions
            extraneous = false(size(c.reactionWholeCellModelIDs));
            for i = 1:numel(c.reactionBoundMonomer)
                if c.reactionBoundMonomer(i) ~= 0
                    displacedID = pm.wholeCellModelIDs{c.reactionBoundMonomer(i)};
                else
                    displacedID = pc.wholeCellModelIDs{c.reactionBoundComplex(i)};
                end
                idx = find(c.reactionMonomerCatalysisMatrix(i, :));
                if ~isempty(idx)
                    releaserID = pm.wholeCellModelIDs{idx};
                else
                    idx = find(c.reactionComplexCatalysisMatrix(i, :));
                    releaserID = pc.wholeCellModelIDs{idx}; %#ok<FNDSB>
                end
                extraneous(i) = ~any(strcmp(allReleaseReactions(:, 1), releaserID) & strcmp(allReleaseReactions(:, 2), displacedID));                
            end
            assertFalse(any(extraneous), ...
                sprintf('%d Reactions must be removed\n- %s', sum(extraneous), strjoin(sprintf('\n- '), c.reactionWholeCellModelIDs{extraneous})));
            
            for i = 1:size(allReleaseReactions, 1)
                tfs = ...
                    c.reactionBoundMonomer == pm.getIndexs(allReleaseReactions{i, 2}) & ...
                    c.reactionBoundComplex == pc.getIndexs(allReleaseReactions{i, 2});
                if pm.getIndexs(allReleaseReactions{i, 1})
                    tfs = tfs & c.reactionMonomerCatalysisMatrix(:, pm.getIndexs(allReleaseReactions{i, 1}));
                else
                    tfs = tfs & c.reactionComplexCatalysisMatrix(:, pc.getIndexs(allReleaseReactions{i, 1}));
                end
                assertEqual(1, sum(tfs));
            end
        end
        
        function assertCorrectDisplacementReactions(this, enyzmeID, expReleasedProteins)
            sim = this.simulation;
            c = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            enyzmeIdx = pm.getIndexs(enyzmeID);
            if enyzmeIdx ~= 0
                releasedProteins = [
                    pm.wholeCellModelIDs(pm.boundIndexs(setdiff(c.reactionBoundMonomer(c.reactionMonomerCatalysisMatrix(:, enyzmeIdx) ~= 0), 0)))
                    pc.wholeCellModelIDs(pc.boundIndexs(setdiff(c.reactionBoundComplex(c.reactionMonomerCatalysisMatrix(:, enyzmeIdx) ~= 0), 0)))
                    ];
            else
                enyzmeIdx = pc.getIndexs(enyzmeID);
                if enyzmeIdx == 0
                    throw(MException('DNA_Test:error', 'undefined protein %s', enzymeId));
                end
                releasedProteins = [
                    pm.wholeCellModelIDs(pm.boundIndexs(setdiff(c.reactionBoundMonomer(c.reactionComplexCatalysisMatrix(:, enyzmeIdx) ~= 0), 0)))
                    pc.wholeCellModelIDs(pc.boundIndexs(setdiff(c.reactionBoundComplex(c.reactionComplexCatalysisMatrix(:, enyzmeIdx) ~= 0), 0)))
                    ];
            end
            
            diff = setdiff(expReleasedProteins, releasedProteins);
            assertEqual(cell(0, 1), diff, sprintf('Reactions must be added for release by %s of these proteins\n- %s\n', enyzmeID, strjoin(sprintf('\n- '), diff{:})));
            diff = setdiff(releasedProteins, expReleasedProteins);
            assertEqual(cell(0, 1), diff, sprintf('Reactions must be removed for release by %s of these proteins\n- %s\n', enyzmeID, strjoin(sprintf('\n- '), diff{:})));
        end
        
        function testTranscriptionDNADamageRepairInteraction_PosStrnd(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            ts = sim.process('Transcription');
            dd = sim.process('DNADamage');
            dr = sim.process('DNARepair');
            rp = sim.process('Replication');
            
            ts.enzymes = ts.enzymes + ts.boundEnzymes;
            ts.boundEnzymes(:) = 0;
            ts.enzymes(ts.enzymeIndexs_rnaPolymerase) = 1;
            ts.enzymes(ts.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0; 
            ts.RNAs(:) = 0;
            ts.copyToState();
            
            iTU = find(t.transcriptionUnitDirections == 1 & t.transcriptionUnitLengths > 1000, 1, 'first');
            r.states = r.activelyTranscribingValue + 1;
            r.positionStrands = [t.transcriptionUnitFivePrimeCoordinates(iTU) + 1  1];
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptChromosome = 1;
            t.boundTranscriptionUnits = iTU;
            
            c.initialize();
            dr.initializeState();
            ts.bindProteinToChromosome(r.positionStrands, ts.enzymeIndexs_rnaPolymerase);
            ts.transcriptionUnitBindingProbabilities(iTU) = 0;
            
            c.strandBreaks(t.transcriptionUnitFivePrimeCoordinates(iTU) + 800, 1) = 1;
            
            repairRates = dr.enzymeBounds;
            dd.reactionBounds(:) = 0;
            dr.enzymeBounds(:) = 0;
            rp.ligaseRate = 0;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            %assert reached damage
            for i = 1:50
                sim.evolveState();
                assertIn(r.positionStrands(1), ...
                    [t.transcriptionUnitFivePrimeCoordinates(iTU)+1  ...
                    t.transcriptionUnitFivePrimeCoordinates(iTU) + 800 - ts.enzymeDNAFootprints3Prime(ts.enzymeIndexs_rnaPolymerase)-1]);
                assertAllEqual(0, ts.RNAs);
                assertEqual(1, nnz(c.strandBreaks));
            end
            
            %assert damage repaired, RNA polymerase released and transcript aborted
            dr.enzymeBounds = repairRates;
            for i = 1:2
                sim.evolveState();
            end
            assertAllEqual(0, ts.RNAs);
            assertEqual(1, numel(t.abortedSequences));
            assertFalse(t.boundTranscriptionUnits == iTU);
            assertEqual(0, nnz(c.strandBreaks));
            assertEqual(1, ...
                + ts.enzymes(ts.enzymeIndexs_rnaPolymerase) ...
                + ts.enzymes(ts.enzymeIndexs_rnaPolymeraseHoloenzyme) ...
                + ts.boundEnzymes(ts.enzymeIndexs_rnaPolymerase) ...
                + ts.boundEnzymes(ts.enzymeIndexs_rnaPolymeraseHoloenzyme));
        end
        
        function testTranscriptionDNADamageRepairInteraction_NegStrnd(this)
            sim = this.simulation;
            comp = sim.compartment;
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            r = sim.state('RNAPolymerase');
            t = sim.state('Transcript');
            ts = sim.process('Transcription');
            dd = sim.process('DNADamage');
            dr = sim.process('DNARepair');
            rp = sim.process('Replication');
            
            c.rnaPolymerase = r;
            c.transcript = t;
            c.complexIndexs_rnaPolymerase = [200; 201];
                        
            ts.enzymes = ts.enzymes + ts.boundEnzymes;
            ts.boundEnzymes(:) = 0;
            ts.enzymes(ts.enzymeIndexs_rnaPolymerase) = 1;
            ts.enzymes(ts.enzymeIndexs_rnaPolymeraseHoloenzyme) = 0; 
            ts.RNAs(:) = 0;
            ts.copyToState();
            
            iTU = find(t.transcriptionUnitDirections == 0 & t.transcriptionUnitLengths > 1000, 1, 'first');
            r.states = r.activelyTranscribingValue + 1;
            r.positionStrands = [t.transcriptionUnitFivePrimeCoordinates(iTU)+t.transcriptionUnitLengths(iTU)-2  2];
            t.boundTranscriptProgress = r.states;
            t.boundTranscriptChromosome = 1;
            t.boundTranscriptionUnits = iTU;
            
            c.initialize();
            dr.initializeState();
            ts.bindProteinToChromosome(r.positionStrands, ts.enzymeIndexs_rnaPolymerase);
            ts.transcriptionUnitBindingProbabilities(iTU) = 0;
            
            c.strandBreaks(t.transcriptionUnitFivePrimeCoordinates(iTU) + t.transcriptionUnitLengths(iTU) - 800, 2) = 1;
            
            repairRates = dr.enzymeBounds;
            dd.reactionBounds(:) = 0;
            dr.enzymeBounds(:) = 0;
            rp.ligaseRate = 0;
            
            m.counts(:, comp.cytosolIndexs) = 1e6;
            
            %assert reached damage
            for i = 1:50
                sim.evolveState();
                assertIn(r.positionStrands(1), ...
                    [t.transcriptionUnitFivePrimeCoordinates(iTU)+t.transcriptionUnitLengths(iTU)-800+ts.enzymeDNAFootprints5Prime(ts.enzymeIndexs_rnaPolymerase)+1+1 ...
                    t.transcriptionUnitFivePrimeCoordinates(iTU)+t.transcriptionUnitLengths(iTU)-2]);
                assertAllEqual(0, ts.RNAs);
                assertEqual(1, nnz(c.strandBreaks));
            end
            
            %assert damage repaired, RNA polymerase released and transcript aborted
            dr.enzymeBounds = repairRates;
            for i = 1:2
                sim.evolveState();
            end
            assertAllEqual(0, ts.RNAs);
            assertEqual(1, numel(t.abortedSequences));
            assertFalse(t.boundTranscriptionUnits == iTU);
            assertEqual(0, nnz(c.strandBreaks));
            assertEqual(1, ...
                + ts.enzymes(ts.enzymeIndexs_rnaPolymerase) ...
                + ts.enzymes(ts.enzymeIndexs_rnaPolymeraseHoloenzyme) ...
                + ts.boundEnzymes(ts.enzymeIndexs_rnaPolymerase) ...
                + ts.boundEnzymes(ts.enzymeIndexs_rnaPolymeraseHoloenzyme));
        end
    end
    
    %test that chromosome modules all work when there are 1-2
    %chromosomes
    %- SMCs bind both chromosomes
    %- both chromosomes damaged, repaired
    %- both chromosomes transcribed
    %- both chromosomes transcriptionally regulated
    %- both chromosomes wound/unwound
    %- no DNA replication
    %- DnaA binds both chromosomes
    methods
        %assert
        %- both chromosomes damaged, repaired
        %- both chromosomes transcribed
        %- both chromosomes transcriptionally regulated
        %- both chromosomes wound/unwound
        %- no DNA replication
        %- DnaA binds both chromosomes
        function testMultipleChromosomes_IndividualProcesses(this)
            %2 full chromosomes
            this.helpTestMultipleChromosomes_individualProcess('DNADamage', 2);
            this.helpTestMultipleChromosomes_individualProcess('DNARepair', 2);
            this.helpTestMultipleChromosomes_individualProcess('Transcription', 2);
            this.helpTestMultipleChromosomes_individualProcess('ChromosomeCondensation', 2);
            this.helpTestMultipleChromosomes_individualProcess('DNASupercoiling', 2);
            this.helpTestMultipleChromosomes_individualProcess('ReplicationInitiation', 2);
            this.helpTestMultipleChromosomes_individualProcess('TranscriptionalRegulation', 2);
            
            %1.5 chromosomes
            this.helpTestMultipleChromosomes_individualProcess('DNADamage', 1.5);
            this.helpTestMultipleChromosomes_individualProcess('DNARepair', 1.5);
            this.helpTestMultipleChromosomes_individualProcess('ChromosomeCondensation', 1.5);
            this.helpTestMultipleChromosomes_individualProcess('DNASupercoiling', 1.5);
            this.helpTestMultipleChromosomes_individualProcess('ReplicationInitiation', 1.5);
            this.helpTestMultipleChromosomes_individualProcess('Transcription', 1.5);
            this.helpTestMultipleChromosomes_individualProcess('TranscriptionalRegulation', 1.5);
        end
        
        function helpTestMultipleChromosomes_individualProcess(this, processId, nChromosomes)
            if nChromosomes == 2
                this.helpTestTwoChromosomes();
            else
                this.helpTestOneAndAHalfChromosomes();
            end
            
            %references
            sim = this.simulation;
            c = sim.state('Chromosome');
            pol = sim.state('RNAPolymerase');
            trnscpt = sim.state('Transcript');
            m = sim.process(processId);
            dd = sim.process('DNADamage');
            dr = sim.process('DNARepair');
            
            dd.reactionBounds(:, 2) = 1e3 * dd.reactionBounds(:, 2);
            dr.enzymeBounds(:, 2) = 1e3 * dr.enzymeBounds(:, 2);
            
            polymerizedRegions = c.polymerizedRegions;
            linkingNumbers = c.linkingNumbers;
            complexBoundSites = c.complexBoundSites;
            monomerBoundSites = c.monomerBoundSites;
            gapSites = c.gapSites;
            abasicSites = c.abasicSites;
            damagedSugarPhosphates = c.damagedSugarPhosphates;
            damagedBases = c.damagedBases;
            intrastrandCrossLinks = c.intrastrandCrossLinks;
            strandBreaks = c.strandBreaks;
            hollidayJunctions = c.hollidayJunctions;
            
            positionStrands = pol.positionStrands;
            states = pol.states;
            boundTranscriptionUnits = trnscpt.boundTranscriptionUnits;
            boundTranscriptProgress = trnscpt.boundTranscriptProgress;
            boundTranscriptChromosome = trnscpt.boundTranscriptChromosome;
            
            cumLinkingNumbersInc = zeros(2, 1);
            cumLinkingNumbersDec = zeros(2, 1);
            cumProteinCnts = zeros(2, 1);
            cumDamageCnts = zeros(2, 1);
            for j = 1:5
                this.seedRandStream(j);
                
                c.polymerizedRegions = polymerizedRegions;
                c.linkingNumbers = linkingNumbers;
                c.complexBoundSites = complexBoundSites;
                c.monomerBoundSites = monomerBoundSites;
                c.gapSites = gapSites;
                c.abasicSites = abasicSites;
                c.damagedSugarPhosphates = damagedSugarPhosphates;
                c.damagedBases = damagedBases;
                c.intrastrandCrossLinks = intrastrandCrossLinks;
                c.strandBreaks = strandBreaks;
                c.hollidayJunctions = hollidayJunctions;
                c.invalidate();
                
                pol.positionStrands = positionStrands;
                pol.states = states;
                trnscpt.boundTranscriptionUnits = boundTranscriptionUnits;
                trnscpt.boundTranscriptProgress = boundTranscriptProgress;
                trnscpt.boundTranscriptChromosome = boundTranscriptChromosome;
                
                m.copyFromState();
                for i = 1:numel(m.enzymes)
                    m.releaseProteinFromChromosome(i, Inf, zeros(0, 2), zeros(0, 1));
                end
                [~, vals] = find(c.monomerBoundSites);
                assertFalse(any(ismember(vals, m.enzymeMonomerGlobalIndexs)));
                [~, vals] = find(c.complexBoundSites);
                assertFalse(any(ismember(vals, m.enzymeComplexGlobalIndexs)));
                m.boundEnzymes(:) = 0;
                
                % simulate
                m.substrates(m.substrateMetaboliteLocalIndexs) = 1e6;
                dd.substrates(dd.substrateMetaboliteLocalIndexs) = 1e6;
                for i = 1:100
                    if isequal(processId, 'DNARepair')
                        dd.evolveState();
                    end
                    m.evolveState();
                    
                    cumLinkingNumbersInc(1) = cumLinkingNumbersInc(1) + max(0, collapse(c.linkingNumbers(:, 1:2) - linkingNumbers(:, 1:2)));
                    cumLinkingNumbersInc(2) = cumLinkingNumbersInc(2) + max(0, collapse(c.linkingNumbers(:, 3:4) - linkingNumbers(:, 3:4)));
                    cumLinkingNumbersDec(1) = cumLinkingNumbersDec(1) + max(0, -collapse(c.linkingNumbers(:, 1:2) - linkingNumbers(:, 1:2)));
                    cumLinkingNumbersDec(2) = cumLinkingNumbersDec(2) + max(0, -collapse(c.linkingNumbers(:, 3:4) - linkingNumbers(:, 3:4)));
                    
                    [posStrnds, vals] = find(c.monomerBoundSites);
                    cumProteinCnts = cumProteinCnts + [
                        sum(ismember(vals, m.enzymeMonomerGlobalIndexs) & posStrnds(:, 2) <=2)
                        sum(ismember(vals, m.enzymeMonomerGlobalIndexs) & posStrnds(:, 2)  >2)];
                    [posStrnds, vals] = find(c.complexBoundSites);
                    cumProteinCnts = cumProteinCnts + [
                        sum(ismember(vals, m.enzymeComplexGlobalIndexs) & posStrnds(:, 2) <=2)
                        sum(ismember(vals, m.enzymeComplexGlobalIndexs) & posStrnds(:, 2)  >2)];
                    
                    cumDamageCnts(1) = cumDamageCnts(1) + nnz(c.damagedSites(:, 1:2));
                    cumDamageCnts(2) = cumDamageCnts(2) + nnz(c.damagedSites(:, 3:4));
                end
            end
            
            %assert linking number equally affected
            assertElementsAlmostEqual(cumLinkingNumbersInc(1), 1/(nChromosomes-1) * cumLinkingNumbersInc(2), 'relative', 100e-2, 100);
            assertElementsAlmostEqual(cumLinkingNumbersDec(1), 1/(nChromosomes-1) * cumLinkingNumbersDec(2), 'relative', 100e-2, 100);
            
            %assert proteins bind both chromosomes equally
            assertElementsAlmostEqual(cumProteinCnts(1), 1/(nChromosomes-1) * cumProteinCnts(2), 'relative', 25e-2);
            
            %assert damage affects chromosomes equally
            assertElementsAlmostEqual(cumDamageCnts(1), 1/(nChromosomes-1) * cumDamageCnts(2), 'relative', 100e-2);
        end
        
        function helpTestTwoChromosomes(this)
            this.setUp();
            
            %% references
            sim = this.simulation;
            c = sim.state('Chromosome');
            trnscpt = sim.state('Transcript');
            met = sim.state('Metabolite');
            pol = sim.state('RNAPolymerase');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            ring = sim.state('FtsZRing');
            
            %% initialize
            %2X RNA, protein
            met.counts = 2 * met.counts;
            rna.counts = 2 * rna.counts;
            pm.counts = 2 * pm.counts;
            pc.counts = 2 * pc.counts;
            
            %2 chromosomes
            c.polymerizedRegions(:, 3:4) = c.polymerizedRegions(:, 1:2);
            c.linkingNumbers(:, 3:4) = c.linkingNumbers(:, 1:2);
            c.complexBoundSites(:, 3:4) = c.complexBoundSites(:, 1:2);
            c.monomerBoundSites(:, 3:4) = c.monomerBoundSites(:, 1:2);
            c.gapSites(:, 3:4) = c.gapSites(:, 1:2);
            c.abasicSites(:, 3:4) = c.abasicSites(:, 1:2);
            c.damagedSugarPhosphates(:, 3:4) = c.damagedSugarPhosphates(:, 1:2);
            c.damagedBases(:, 3:4) = c.damagedBases(:, 1:2);
            c.intrastrandCrossLinks(:, 3:4) = c.intrastrandCrossLinks(:, 1:2);
            c.strandBreaks(:, 3:4) = c.strandBreaks(:, 1:2);
            c.hollidayJunctions(:, 3:4) = c.hollidayJunctions(:, 1:2);
            
            %2X active RNA polymerases
            pol.positionStrands = [
                pol.positionStrands
                pol.positionStrands + [zeros(size(pol.positionStrands, 1), 1)   2 * (pol.positionStrands(:, 2) ~= 0)]];
            pol.states = [pol.states; pol.states];
            trnscpt.boundTranscriptionUnits = [
                trnscpt.boundTranscriptionUnits
                trnscpt.boundTranscriptionUnits];
            trnscpt.boundTranscriptProgress = [
                trnscpt.boundTranscriptProgress
                trnscpt.boundTranscriptProgress];
            trnscpt.boundTranscriptChromosome = [
                trnscpt.boundTranscriptChromosome
                2 * trnscpt.boundTranscriptChromosome];
            
            %transcription fold changes
            pol.transcriptionFactorBindingProbFoldChange(:, 2) = pol.transcriptionFactorBindingProbFoldChange(:, 1);
            pol.supercoilingBindingProbFoldChange(:, 2) = pol.supercoilingBindingProbFoldChange(:, 1);
            
            %aborted sequences
            trnscpt.abortedTranscripts = zeros(0, 2);
            
            %FtsZ ring
            ring.numEdgesOneStraight = 2 * ring.numEdgesOneStraight;
            ring.numEdgesTwoStraight = 2 * ring.numEdgesTwoStraight;
            ring.numEdgesTwoBent = 2 * ring.numEdgesTwoBent;
            ring.numResidualBent = 2 * ring.numResidualBent;
        end
        
        function helpTestOneAndAHalfChromosomes(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            this.setUp();
            
            %% references
            sim = this.simulation;
            comp = sim.compartment;
            c = sim.state('Chromosome');
            trnscpt = sim.state('Transcript');
            met = sim.state('Metabolite');
            pol = sim.state('RNAPolymerase');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            rib = sim.state('Ribosome');
            poly = sim.state('Polypeptide');
            ring = sim.state('FtsZRing');
            tr = sim.process('TranscriptionalRegulation');
            sc = sim.process('DNASupercoiling');
            dr = sim.process('DNARepair');
            
            %% initialize
            %2 chromosomes
            posStrnds = [
                1 1
                1 2
                c.sequenceLen*1/4+1  2
                c.sequenceLen*3/4+1+1000 2
                1 3
                c.sequenceLen*3/4+1+1000 3
                1 4
                c.sequenceLen*3/4+1 4
                ];
            lengths = [
                c.sequenceLen
                c.sequenceLen/4-1000
                c.sequenceLen/2
                c.sequenceLen/4-1000
                c.sequenceLen/4-1000
                c.sequenceLen/4-1000
                c.sequenceLen/4
                c.sequenceLen/4
                ];
            c.polymerizedRegions = CircularSparseMat(posStrnds, lengths, [c.sequenceLen c.nCompartments], 1);
            
            posStrnds = [
                1 1
                1 2
                1 3
                1 4
                c.sequenceLen*3/4+1+1000 1
                c.sequenceLen*3/4+1+1000 2
                c.sequenceLen*3/4+1+1000 3
                c.sequenceLen*3/4+1+1000 4
                c.sequenceLen*1/4+1  1
                c.sequenceLen*1/4+1  2
                ];
            lengths = [
                c.sequenceLen/4-1000
                c.sequenceLen/4-1000
                c.sequenceLen/4-1000
                c.sequenceLen/4-1000
                c.sequenceLen/4-1000
                c.sequenceLen/4-1000
                c.sequenceLen/4-1000
                c.sequenceLen/4-1000
                c.sequenceLen/2
                c.sequenceLen/2
                ];
            c.linkingNumbers = CircularSparseMat(posStrnds, lengths / c.relaxedBasesPerTurn * (1 + c.equilibriumSuperhelicalDensity), [c.sequenceLen c.nCompartments], 1);
            
            [posStrnds, vals] = find(c.monomerBoundSites);
            posStrnds = [posStrnds; posStrnds(:, 1) posStrnds(:, 2) + 2];
            vals = [vals; vals];
            tfs = ...
                (posStrnds(:, 1) >= c.sequenceLen/4-2000 & ...
                posStrnds(:, 1) <= c.sequenceLen/4) | ...
                (posStrnds(:, 1) >= 3/4*c.sequenceLen/4-1000 & ...
                posStrnds(:, 1) <= 3/4*c.sequenceLen/4+1000) | ...
                (posStrnds(:, 1) >= 1/4*c.sequenceLen/4-2000 & ...
                posStrnds(:, 1) <= 3/4*c.sequenceLen/4+1000 & ...
                posStrnds(:, 2) >= 3);
            c.monomerBoundSites = CircularSparseMat(posStrnds(~tfs, :), vals(~tfs), [c.sequenceLen c.nCompartments], 1);
            
            [posStrnds, vals] = find(c.complexBoundSites);
            posStrnds = [posStrnds; posStrnds(:, 1) posStrnds(:, 2) + 2];
            vals = [vals; vals];
            tfs = ...
                (posStrnds(:, 1) >= c.sequenceLen/4-2000 & ...
                posStrnds(:, 1) <= c.sequenceLen/4) | ...
                (posStrnds(:, 1) >= 3/4*c.sequenceLen/4-1000 & ...
                posStrnds(:, 1) <= 3/4*c.sequenceLen/4+1000) | ...
                (posStrnds(:, 1) >= 1/4*c.sequenceLen/4-2000 & ...
                posStrnds(:, 1) <= 3/4*c.sequenceLen/4+1000 & ...
                posStrnds(:, 2) >= 3);
            c.complexBoundSites = CircularSparseMat(posStrnds(~tfs, :), vals(~tfs), [c.sequenceLen c.nCompartments], 1);
            
            posStrnds1 = dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1));
            posStrnds2 = dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1));
            posStrnds3 = dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2));
            posStrnds4 = dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2));
            tfs2 = ...
                posStrnds2 <= c.sequenceLen/4-1000 | ...
                (posStrnds2 >= c.sequenceLen/4+1 & posStrnds2 <= c.sequenceLen*3/4) | ...
                posStrnds2 >= c.sequenceLen*3/4+1000+1;
            tfs3 = ...
                posStrnds2 <= c.sequenceLen/4-1000 | ...
                posStrnds2 >= c.sequenceLen*3/4+1000+1;
            tfs4 = ...
                posStrnds2 <= c.sequenceLen/4 | ...
                posStrnds2 >= c.sequenceLen*3/4+1;
            posStrnds = [
                posStrnds1 ones(size(posStrnds1))
                posStrnds2(tfs2) 2*ones(sum(tfs2), 1)
                posStrnds3(tfs3) 3*ones(sum(tfs3), 1)
                posStrnds4(tfs4) 4*ones(sum(tfs4), 1)
                ];
            c.damagedBases = CircularSparseMat(posStrnds, met.m6ADIndexs, [c.sequenceLen c.nCompartments], 1);
            
            c.gapSites = CircularSparseMat([], [], [c.sequenceLen c.nCompartments], 1);
            c.abasicSites = CircularSparseMat([], [], [c.sequenceLen c.nCompartments], 1);
            c.damagedSugarPhosphates = CircularSparseMat([], [], [c.sequenceLen c.nCompartments], 1);
            c.intrastrandCrossLinks = CircularSparseMat([], [], [c.sequenceLen c.nCompartments], 1);
            c.strandBreaks = CircularSparseMat([], [], [c.sequenceLen c.nCompartments], 1);
            c.hollidayJunctions = CircularSparseMat([], [], [c.sequenceLen c.nCompartments], 1);
            
            %2X active RNA polymerases
            pol.states = [pol.states; pol.states];
            pol.positionStrands = [
                pol.positionStrands
                pol.positionStrands + [zeros(size(pol.positionStrands, 1), 1)   2 * (pol.positionStrands(:, 2) ~= 0)]];
            trnscpt.boundTranscriptionUnits = [
                trnscpt.boundTranscriptionUnits
                trnscpt.boundTranscriptionUnits];
            trnscpt.boundTranscriptProgress = [
                trnscpt.boundTranscriptProgress
                trnscpt.boundTranscriptProgress];
            trnscpt.boundTranscriptChromosome = [
                trnscpt.boundTranscriptChromosome
                2 * trnscpt.boundTranscriptChromosome];
            
            tfs = ...
                (pol.positionStrands(:, 1) >= c.sequenceLen/4-2000 & ...
                pol.positionStrands(:, 1) <= c.sequenceLen/4) | ...
                (pol.positionStrands(:, 1) >= 3/4*c.sequenceLen/4-1000 & ...
                pol.positionStrands(:, 1) <= 3/4*c.sequenceLen/4+1000) | ...
                (pol.positionStrands(:, 1) >= 1/4*c.sequenceLen/4-2000 & ...
                pol.positionStrands(:, 1) <= 3/4*c.sequenceLen/4+1000 & ...
                pol.positionStrands(:, 2) >= 3);
            
            assertEqual(sum(pol.positionStrands(:,1) ~= 0) - sum(tfs), ...
                nnz(c.complexBoundSites == pc.getIndexs('RNA_POLYMERASE')) + ...
                nnz(c.complexBoundSites == pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')));
            
            pol.states(tfs) = pol.freeValue;
            pol.positionStrands(tfs, :) = 0;
            trnscpt.boundTranscriptionUnits(tfs) = 0;
            trnscpt.boundTranscriptProgress(tfs) = 0;
            trnscpt.boundTranscriptChromosome(tfs) = 0;
            
            %transcription fold changes
            pol.transcriptionFactorBindingProbFoldChange = ...
                tr.calcBindingProbabilityFoldChange(tr.tfBoundPromoters);
            pol.supercoilingBindingProbFoldChange = ...
                sc.calcRNAPolymeraseBindingProbFoldChange();
            
            %aborted sequences
            trnscpt.abortedTranscripts = zeros(0, 2);
            
            %1.5X RNA, protein except bound tmRNA, translation factors,
            %ribosomes, FtsZ
            met.counts = ceil(1.5 * met.counts);
            rna.counts = ceil(1.5 * rna.counts);
            pm.counts = ceil(1.5 * pm.counts);
            pc.counts = ceil(1.5 * pc.counts);
            
            pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE')), comp.cytosolIndexs) = sum(pol.states == pol.freeValue);
            
            rna.counts(rna.boundIndexs, comp.cytosolIndexs) = 0;
            [~, vals] = find(c.monomerBoundSites);
            pm.counts(pm.boundIndexs, comp.cytosolIndexs) = histc(vals, (1:numel(pm.boundIndexs))');
            [~, vals] = find(c.complexBoundSites);
            pc.counts(pc.boundIndexs, comp.cytosolIndexs) = histc(vals, (1:numel(pc.boundIndexs))');
            
            pm.counts(pm.matureIndexs(tr.enzymeMonomerGlobalIndexs), comp.cytosolIndexs) = 20;
            pc.counts(pc.matureIndexs(tr.enzymeComplexGlobalIndexs), comp.cytosolIndexs) = 20;
            
            %ribosome
            rib.states(:) = rib.notExistValue;
            rib.boundMRNAs(:) = 0;
            rib.mRNAPositions(:) = 0;
            rib.tmRNAPositions(:) = 0;
            poly.boundMRNAs = rib.boundMRNAs;
            poly.nascentMonomerLengths = rib.mRNAPositions;
            poly.proteolysisTagLengths = rib.tmRNAPositions;
            poly.abortedPolypeptides = zeros(0, 3);
            
            %FtsZ ring
            ring.numEdgesOneStraight = 0;
            ring.numEdgesTwoStraight = 0;
            ring.numEdgesTwoBent = 0;
            ring.numResidualBent = 0;
        end
    end
    
    methods
        function testReplicationInitiationInContextOfOtherChromosomeModules(this)
            %% references
            sim = this.simulation;
            %sim.applyOptions('verbosity', 1);
            %this.seedRandStream(round(mod(now, 1) * 1e7));
            
            time = sim.state('Time');
            met = sim.state('Metabolite');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            c = sim.state('Chromosome');
            cc = sim.process('ChromosomeCondensation');
            repInit = sim.process('ReplicationInitiation');
            rep = sim.process('Replication');
            dr = sim.process('DNARepair');
            tr = sim.process('TranscriptionalRegulation');
            
            %% constants
            nCompartments = size(pm.counts, 2);
            
            %% setup
            %mean DnaA expression so that test doesn't take too long
            repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ATP) = ...
                + repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ATP) ...
                + max(0, 54 - repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ATP) - ...
                repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_1mer_ATP));
            repInit.copyToState();
            
            %% store initial state
            initMonomers = pm.counts;
            initMonomers(pm.matureIndexs, :) = ...
                + initMonomers(pm.matureIndexs, :) ...
                + initMonomers(pm.boundIndexs, :);
            initMonomers(pm.boundIndexs, :) = 0;
            initComplexs = pc.counts;
            initComplexs(pc.matureIndexs, :) = ...
                + initComplexs(pc.matureIndexs, :) ...
                + initComplexs(pc.boundIndexs, :);
            initComplexs(pc.boundIndexs, :) = 0;
            
            [~, vals] = find(c.monomerBoundSites);
            initBoundMonomers = reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1);
            [~, vals] = find(c.complexBoundSites);
            initBoundComplexs = reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1);
            
            %% simulate
            iterMax = 1.5 * time.replicationInitiationDuration;
            if sim.verbosity > 0
                fprintf('%5s %3s %7s %9s\n', 'Iter ', 'Sup', ' DnaA  ', '  R1-5   ');
                fprintf('%5s %3s %7s %9s\n', '=====', '===', '=======', '=========');
            end
            for iter = 1:iterMax
                if sim.verbosity > 0 && mod(iter, 100) == 1
                    dnaAPol = repInit.calculateDnaAR1234Polymerization();
                    fprintf('%5d %3d %3d %3d %d %d %d %d %d\n', iter, ...
                        collapse(c.supercoiled)/2, ...
                        repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ATP), ...
                        repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_1mer_ATP), ...
                        dnaAPol(1, 1), dnaAPol(2, 1), dnaAPol(3, 1), dnaAPol(4, 1), ...
                        any(repInit.calcuateIsDnaAR5Occupied()));
                end
                
                %mock metabolism
                met.counts = met.counts + ...
                    log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength) * ...
                    (met.biomassProduction - met.byproducts);
                
                %mock protein synthesis and decay
                decayedMonomers = pm.randStream.stochasticRound(pm.decayRates(:, ones(nCompartments, 1)) .* pm.counts);
                decayedComplexs = pc.randStream.stochasticRound(pc.decayRates(:, ones(nCompartments, 1)) .* pc.counts);
                notUpdatingMonomers = pm.updateExternalState(-decayedMonomers, true);
                notUpdatingComplexs = pc.updateExternalState(-decayedComplexs, true);
                pm.counts = pm.counts + pm.randStream.stochasticRound(...
                    + initMonomers * log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength) ...
                    - (decayedMonomers - notUpdatingMonomers));
                pc.counts = pc.counts + pc.randStream.stochasticRound(...
                    + initComplexs * log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength) ...
                    - (decayedComplexs - notUpdatingComplexs));
                
                %simulate
                sim.evolveState();
                
                %stop once replication initiated
                if rep.isAnyHelicaseBound && ...
                        ~any(any(repInit.calculateDnaAR1234Polymerization())) && ...
                        ~any(repInit.calcuateIsDnaAR5Occupied()) && ...
                        ~any(repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ATP(2:end)))
                    break;
                end
            end
            
            %% assertions
           
            %replication initiated
            assertFalse(rep.isDnaAORIComplexAssembled());
            assertTrue(all(rep.helicasePosition));
            assertTrue(all(rep.leadingPolymerasePosition));
            assertAllEqual(0, repInit.calculateDnaAR1234Polymerization());
            assertAllEqual(false, repInit.calcuateIsDnaAR5Occupied());
            assertIn(collapse(c.polymerizedRegions), [2 * c.sequenceLen + 1 4 * c.sequenceLen]);
            
            %time of replication approximately correct
            assertElementsAlmostEqual(iter, time.replicationInitiationDuration, 'relative', 0.9);
            
            %chromosome supercoiled
            assertIn(collapse(c.supercoiled), [2 Inf]);
            
            %not too much damage, R/M sites fully methylated
            assertIn(nnz(c.damagedSites), [0 10]);
            
            posStrnds = find(c.restrictableMunIRMSites);
            assertIn(sum(c.isRegionDoubleStranded(posStrnds, 1, false)), [0 1]); %#ok<FNDSB>
            posStrnds = find(c.hemiunmethylatedMunIRMSites);
            assertIn(sum(c.isRegionDoubleStranded(posStrnds, 1, false)), [0 2]); %#ok<FNDSB>
            posStrnds = [
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))      ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2))  2 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)];
            [~, ~, posStrnds] = c.isRegionDoubleStranded(posStrnds, 1, false);
            assertAllEqual(met.m6ADIndexs, c.damagedBases(posStrnds));
            
            %appropriate distribution of bound proteins
            %- gyrase, topoisomerase
            %- DnaA
            %- SMCs
            %- RNA polymerase, transcription factors
            %- repair machinery
            [~, vals] = find(c.monomerBoundSites);
            boundMonomers = reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1);
            [~, vals] = find(c.complexBoundSites);
            boundComplexs = reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1);
            
            assertElementsAlmostEqual(initBoundMonomers, boundMonomers, 'relative', 0.5, 4);
            
            assertElementsAlmostEqual(...
                initBoundComplexs([cc.enzymeComplexGlobalIndexs; tr.enzymeComplexGlobalIndexs]), ...
                boundComplexs([cc.enzymeComplexGlobalIndexs; tr.enzymeComplexGlobalIndexs]), ...
                'relative', 0.5, 4);
            assertElementsAlmostEqual(...
                exp(log(2) * iter / time.replicationInitiationDuration) * ...
                initBoundComplexs(setdiff(1:end, [cc.enzymeComplexGlobalIndexs; tr.enzymeComplexGlobalIndexs; rep.enzymeComplexGlobalIndexs; repInit.enzymeComplexGlobalIndexs])), ...
                boundComplexs(setdiff(1:end, [cc.enzymeComplexGlobalIndexs; tr.enzymeComplexGlobalIndexs; rep.enzymeComplexGlobalIndexs; repInit.enzymeComplexGlobalIndexs])), ...
                'relative', 0.5, 4);
            
            %DnaA 7mer broken down, but not renegerated yet
            assertElementsAlmostEqual(4 * 7 + 1, ...
                + repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(1)) ...
                + repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(1)), ...
                'relative', 0.1, 1);
            assertAllEqual(0, repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ATP(2:end)));
            assertAllEqual(0, repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(2:end)));
            assertAllEqual(0, repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_Nmer_ATP(2:end)));
            assertAllEqual(0, repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(2:end)));
        end
        
        %- test that replication proceeds correctly in presence of SMCs,
        %  DnaA, RNA polymerase, transcription factors, topoisomerases and
        %  gyrases, damage, repair enzymes
        %- test that metabolism can support replication
        %- test that R/M sites get fully methylated
        %- test that equilibrium superhelicity restored
        %- test that DNA pol stalls on contact with RNA pol
        function testReplicationInContextOfOtherChromosomeModules(this)
            %% references
            sim = this.simulation;
            %sim.applyOptions('verbosity', 1);
            %this.seedRandStream(round(mod(now, 1) * 1e7));
            
            comp = sim.compartment;
            time = sim.state('Time');
            c = sim.state('Chromosome');
            met = sim.state('Metabolite');
            rna = sim.state('Rna');
            trnscpt = sim.state('Transcript');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            rnaPol = sim.state('RNAPolymerase');
            rep = sim.process('Replication');
            repInit = sim.process('ReplicationInitiation');
            dr = sim.process('DNARepair');
            transcription = sim.process('Transcription');
            cc = sim.process('ChromosomeCondensation');
            tr = sim.process('TranscriptionalRegulation');
            
            %% assert
            %metabolism can support replication
            %- sufficient dNTPs for polymerization
            %- sufficient NAD for ligation
            %- sufficient m6ad for R/M site methylation
            f = exp(log(2) * (1 + 10e-2) * (time.replicationInitiationDuration + time.replicationDuration) / time.cellCycleLength);
            metVals = met.counts + (f - 1) * (met.biomassProduction - met.byproducts);
            nOkazakiFragments = sum(cellfun(@numel, rep.primaseBindingLocations)) + 2;
            assertIn(metVals(met.dntpIndexs, comp.cytosolIndexs), [getBaseCounts(c.sequence) Inf(4, 1)]);
            assertIn(metVals(met.getIndexs('NAD'), comp.cytosolIndexs), [nOkazakiFragments Inf]);
            assertIn(metVals(met.atpIndexs, comp.cytosolIndexs), [c.sequenceLen + nOkazakiFragments, Inf]);
            assertIn(metVals(met.waterIndexs, comp.cytosolIndexs), [c.sequenceLen + nOkazakiFragments, Inf]);
            assertIn(metVals(met.getIndexs('AMET'), comp.cytosolIndexs), [2 * size(dr.RM_MunI_MethylatedPositions, 1), Inf]);
            
            %% setup
            f = exp(log(2) * (1 + 10e-2) * (time.replicationInitiationDuration + time.replicationDuration) / time.cellCycleLength);
            met.counts = ceil(met.counts + (f - 1) * (met.biomassProduction - met.byproducts));
            rna.counts(setdiff(1:end, rna.boundIndexs), :) = ceil(f * rna.counts(setdiff(1:end, rna.boundIndexs), :));
            pm.counts(setdiff(1:end, pm.boundIndexs), :) = ceil(f * pm.counts(setdiff(1:end, pm.boundIndexs), :));
            pc.counts(setdiff(1:end, pc.boundIndexs), :) = ceil(f * pc.counts(setdiff(1:end, pc.boundIndexs), :));
            rna.counts(rna.matureIndexs, :) = rna.counts(rna.matureIndexs, :) + rna.counts(rna.boundIndexs, :);
            pm.counts(pm.matureIndexs, :) = pm.counts(pm.matureIndexs, :) + pm.counts(pm.boundIndexs, :);
            pc.counts(pc.matureIndexs, :) = pc.counts(pc.matureIndexs, :) + pc.counts(pc.boundIndexs, :);
            
            %equilibrate chromosome binding new additional proteins
            for iter = 1:100
                sim.evolveState();
            end
            
            %% setup for replication
            rep.copyFromState();
            rep.substrates(setdiff(1:end, rep.substrateIndexs_ntp)) = 1e8;
            c.complexBoundSites(rep.dnaAFunctionalBoxStartPositions(rep.dnaAFunctionalBoxIndexs_R1234), 1) = rep.complexIndexs_DnaA_7mer_ATP;
            c.complexBoundSites(rep.dnaAFunctionalBoxStartPositions(rep.dnaAFunctionalBoxIndexs_R5), 1) = rep.complexIndexs_DnaA_1mer_ATP;
            pc.counts(pc.boundIndexs(rep.complexIndexs_DnaA_7mer_ATP), comp.cytosolIndexs) = ...
                nnz(c.complexBoundSites == rep.complexIndexs_DnaA_7mer_ATP);
            pc.counts(pc.boundIndexs(rep.complexIndexs_DnaA_1mer_ATP), comp.cytosolIndexs) = ...
                nnz(c.complexBoundSites == rep.complexIndexs_DnaA_1mer_ATP);
            rep.copyToState();
            
            %% store initial state for later comparison
            [~, vals] = find(c.monomerBoundSites);
            initBoundMonomers = reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1);
            [~, vals] = find(c.complexBoundSites);
            initBoundComplexs = reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1);
            
            %% simulate
            iterMax = ceil(1.35 * time.replicationDuration);
            replicationFinishedIter = NaN;
            
            ticHandle = tic;
            if sim.verbosity > 0
                fprintf('%5s %5s %13s %13s %13s\n', 'Iter ', 'Time ', 'Helicase Pos ', ' Leading Pol ', ' Lagging Pol ');
                fprintf('%5s %5s %13s %13s %13s\n', '=====', '=====', '=============', '=============', '=============');
            end
            for iter = 1:iterMax
                if sim.verbosity > 0 && mod(iter, 100) == 1
                    fprintf('%5d %5.2f %6d %6d %6d %6d %6d %6d\n', ...
                        iter, toc(ticHandle), ...
                        rep.helicasePosition(1), rep.helicasePosition(2), ...
                        rep.leadingPosition(1), rep.leadingPosition(2), ...
                        rep.laggingPosition(1), rep.laggingPosition(2));
                    ticHandle = tic;
                end
                sim.evolveState();
                if all(rep.strandDuplicated)
                    replicationFinishedIter = min(replicationFinishedIter, iter);
                    if ~nnz(c.hemiunmethylatedMunIRMSites) && ~repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(1))
                        break;
                    end
                end
            end
            
            %% assert replication proceeeded correctly in presence of other proteins
            
            %DnaA not polymerized
            assertElementsAlmostEqual(zeros(4, 2), repInit.calculateDnaAR1234Polymerization(), 'absolute', 7);
            %assertAllEqual(false, repInit.calcuateIsDnaAR5Occupied());
            
            %replicated
            assertTrue(all(rep.strandDuplicated));
            
            %equilibrium superhelical density restored
            assertEqual([ones(4, 1) (1:4)'], find(c.linkingNumbers));
            assertTrue(all(full(c.supercoiled(1, :))));
            
            %R/M sites fully methylated
            assertEqual(0, nnz(c.restrictableMunIRMSites));
            assertEqual(0, nnz(c.hemiunmethylatedMunIRMSites));
            assertAllEqual(met.m6ADIndexs, c.damagedBases([
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))      ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2))  2 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))  3 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2))  4 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)]));
            
            %DNA polymerase stalled, and RNA polymerase aborted resulting
            %in incomplete transcripts
            nActRNAPol = rnaPol.stateExpectations(rnaPol.activelyTranscribingIndex) * (...
                + pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE')), comp.cytosolIndexs) ...
                + pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')), comp.cytosolIndexs) ...
                + pc.counts(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE')), comp.cytosolIndexs) ...
                + pc.counts(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')), comp.cytosolIndexs));
            assertElementsAlmostEqual(nActRNAPol, ...
                numel(trnscpt.abortedSequences), ...
                'relative', 1);
            probs = transcription.computeRNAPolymeraseTUBindingProbabilities();
            nStalls = nActRNAPol * (...
                + max(...
                    sum(probs(trnscpt.transcriptionUnitDirections == 1 & trnscpt.transcriptionUnitFivePrimeCoordinates <= c.sequenceLen/2, 1)), ...
                    sum(probs(trnscpt.transcriptionUnitDirections == 0 & trnscpt.transcriptionUnitFivePrimeCoordinates <= c.sequenceLen/2, 1))) ...
                + max(...
                    sum(probs(trnscpt.transcriptionUnitDirections == 1 & trnscpt.transcriptionUnitFivePrimeCoordinates > c.sequenceLen/2, 1)), ...
                    sum(probs(trnscpt.transcriptionUnitDirections == 0 & trnscpt.transcriptionUnitFivePrimeCoordinates > c.sequenceLen/2, 1)))) ...
                / sum(probs(:, 1));
            assertIn(replicationFinishedIter, ...
                [1 + c.sequenceLen / 2 / rep.dnaPolymeraseElongationRate + nStalls * rep.rnaPolymeraseCollisionMeanDwellTime  Inf]);
            
            %distribution of bound proteins restored
            %- for most proteins this means similar number bound as before
            %- except for SMC which should have the bound density preserved
            [~, vals] = find(c.monomerBoundSites);
            boundMonomers = reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1);
            [~, vals] = find(c.complexBoundSites);
            boundComplexs = reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1);
            
            assertElementsAlmostEqual(initBoundMonomers, boundMonomers, ...
                'relative', 0.5, 4);
            smcIdxs = cc.enzymeGlobalIndexs(cc.enzymeIndexs_SMC_ADP);
            assertElementsAlmostEqual(...
                initBoundComplexs(setdiff(1:end, [smcIdxs; repInit.enzymeComplexGlobalIndexs; tr.enzymeComplexGlobalIndexs])), ...
                boundComplexs(setdiff(1:end, [smcIdxs; repInit.enzymeComplexGlobalIndexs; tr.enzymeComplexGlobalIndexs])), ...
                'relative', 0.75, 4);
            assertElementsAlmostEqual(...
                2 * initBoundComplexs([smcIdxs; tr.enzymeComplexGlobalIndexs]), ...
                boundComplexs([smcIdxs; tr.enzymeComplexGlobalIndexs]), ...
                'relative', 0.75, 4);
            
            %DnaA 7mer broken down and regenerated
            assertElementsAlmostEqual(.45^(replicationFinishedIter / 3600) * (4 * 7 + 1), ...
                repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ADP) + repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_1mer_ADP), ...
                'relative', 0.35, 10);
            assertAllEqual(0, repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ATP(2:end)));
            assertAllEqual(0, repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(2:end)));
            assertElementsAlmostEqual(0, repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_Nmer_ATP(end)), 'absolute', 8);
            assertAllEqual(0, repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(2:end)));
        end
        
        %- test that replication proceeds correctly in presence of SMCs,
        %  DnaA, RNA polymerase, transcription factors, topoisomerases and
        %  gyrases, damage, repair enzymes
        %- test that metabolism can support replication
        %- test that R/M sites get fully methylated
        %- test that equilibrium superhelicity restored
        %- test that DNA pol stalls on contact with RNA pol
        function testInitiationAndReplicationInContextOfOtherChromosomeModules(this)
            %% references
            sim = this.simulation;
            %sim.applyOptions('verbosity', 1);
            %this.seedRandStream(round(mod(now, 1) * 1e7));
            
            comp = sim.compartment;
            
            time = sim.state('Time');
            met = sim.state('Metabolite');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            c = sim.state('Chromosome');
            rna = sim.state('Rna');
            trnscpt = sim.state('Transcript');
            rnaPol = sim.state('RNAPolymerase');
            
            repInit = sim.process('ReplicationInitiation');
            rep = sim.process('Replication');
            dr = sim.process('DNARepair');
            
            %% constants
            nCompartments = size(pm.counts, 2);
            
            %% setup
            %mean DnaA expression so test doesn't take too long
            repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ATP) = ...
                + repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ATP) ...
                + max(0, 54 - repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ATP) - ...
                repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_1mer_ATP));
            repInit.copyToState();
            
            %% store initial state
            initRNAs = rna.counts;
            initRNAs(rna.matureIndexs, :) = ...
                + initRNAs(rna.matureIndexs, :) ...
                + initRNAs(rna.boundIndexs, :);
            initRNAs(rna.boundIndexs, :) = 0;
            initMonomers = pm.counts;
            initMonomers(pm.matureIndexs, :) = ...
                + initMonomers(pm.matureIndexs, :) ...
                + initMonomers(pm.boundIndexs, :);
            initMonomers(pm.boundIndexs, :) = 0;
            initComplexs = pc.counts;
            initComplexs(pc.matureIndexs, :) = ...
                + initComplexs(pc.matureIndexs, :) ...
                + initComplexs(pc.boundIndexs, :);
            initComplexs(pc.boundIndexs, :) = 0;
            
            [~, vals] = find(c.monomerBoundSites);
            initBoundMonomers = reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1);
            [~, vals] = find(c.complexBoundSites);
            initBoundComplexs = reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1);
            
            %% simulate
            iterMax = 1.5 * (...
                + time.replicationInitiationDuration ...
                + time.replicationDuration);
            replicationInitiationFinishedIter = NaN;
            replicationFinishedIter = NaN;
            if sim.verbosity > 0
                fprintf('%5s %7s %9s %13s %13s %13s\n', ...
                    'Iter', ' DnaA  ', '  R1-5   ', 'Helicase Pos ', ' Leading Pol ', ' Lagging Pol ');
                fprintf('%5s %7s %9s %13s %13s %13s\n', ...
                    '=====', '=======', '=========', '=============', '=============', '=============');
            end
            for iter = 1:iterMax
                if sim.verbosity > 0 && mod(iter, 100) == 1
                    dnaAPol = repInit.calculateDnaAR1234Polymerization();
                    fprintf('%5d %3d %3d %d %d %d %d %d %6d %6d %6d %6d %6d %6d\n', iter, ...
                        repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ATP), repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_1mer_ATP), ...
                        dnaAPol(1, 1), dnaAPol(2, 1), dnaAPol(3, 1), dnaAPol(4, 1), ...
                        any(repInit.calcuateIsDnaAR5Occupied()), ...
                        rep.helicasePosition(1), rep.helicasePosition(2), ...
                        rep.leadingPosition(1), rep.leadingPosition(2), ...
                        rep.laggingPosition(1), rep.laggingPosition(2));
                end
                
                %mock metabolism
                met.counts = met.counts + ...
                    log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength) * ...
                    (met.biomassProduction - met.byproducts);
                
                %mock rna synthesis and decay
                decayedRNAs = rna.randStream.stochasticRound(rna.decayRates(:, ones(nCompartments, 1)) .* rna.counts);
                notUpdatingRNAs = rna.updateExternalState(-decayedRNAs, true);
                rna.counts = rna.counts + rna.randStream.stochasticRound(...
                    + initRNAs * log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength) ...
                    - (decayedRNAs - notUpdatingRNAs));
                
                %mock protein synthesis and decay
                decayedMonomers = pm.randStream.stochasticRound(pm.decayRates(:, ones(nCompartments, 1)) .* pm.counts);
                decayedComplexs = pc.randStream.stochasticRound(pc.decayRates(:, ones(nCompartments, 1)) .* pc.counts);
                notUpdatingMonomers = pm.updateExternalState(-decayedMonomers, true);
                notUpdatingComplexs = pc.updateExternalState(-decayedComplexs, true);
                pm.counts = pm.counts + pm.randStream.stochasticRound(...
                    + initMonomers * log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength) ...
                    - (decayedMonomers - notUpdatingMonomers));
                pc.counts = pc.counts + pc.randStream.stochasticRound(...
                    + initComplexs * log(2) / time.cellCycleLength * exp(log(2) * iter / time.cellCycleLength) ...
                    - (decayedComplexs - notUpdatingComplexs));
                
                %simulate
                sim.evolveState();
                
                if isnan(replicationInitiationFinishedIter) && ...
                        rep.isAnyHelicaseBound && ...
                        ~any(any(repInit.calculateDnaAR1234Polymerization())) && ...
                        ~any(repInit.calcuateIsDnaAR5Occupied())
                    replicationInitiationFinishedIter = iter;
                end
                
                %stop once replication initiated
                if all(rep.strandDuplicated)
                    replicationFinishedIter = min(replicationFinishedIter, iter);
                    if ~nnz(c.hemiunmethylatedMunIRMSites) && ~repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(1))
                        break;
                    end
                end
            end
            
            %% assert replication proceeeded correctly in presence of other proteins
            
            %replicated
            assertTrue(all(rep.strandDuplicated));
            
            %equilibrium superhelical density restored
            assertEqual([ones(4, 1) (1:4)'], find(c.linkingNumbers));
            assertTrue(all(full(c.supercoiled(1, :))));
            
            %not too much damage
            assertIn(nnz(c.damagedSites), [0 10]);

            %R/M sites fully methylated
            assertIn(nnz(c.restrictableMunIRMSites), [0 1]);
            assertIn(nnz(c.hemiunmethylatedMunIRMSites), [0 1]);
            assertAllEqual(met.m6ADIndexs, c.damagedBases([
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))      ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2))  2 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))  3 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2))  4 * ones(size(dr.RM_MunI_RecognitionSites, 1), 1)]));
            
            %DNA polymerase stalled, and RNA polymerase aborted resulting
            %in incomplete transcripts
            nActRNAPol = rnaPol.stateExpectations(rnaPol.activelyTranscribingIndex) * (...
                + pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE')), comp.cytosolIndexs) ...
                + pc.counts(pc.matureIndexs(pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')), comp.cytosolIndexs) ...
                + pc.counts(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE')), comp.cytosolIndexs) ...
                + pc.counts(pc.boundIndexs(pc.getIndexs('RNA_POLYMERASE_HOLOENZYME')), comp.cytosolIndexs));
            assertElementsAlmostEqual(nActRNAPol, ...
                numel(trnscpt.abortedSequences), ...
                'relative', 1);
            assertElementsAlmostEqual(iter, time.replicationInitiationDuration + time.replicationDuration, ...
                'relative', 0.8);
            
            %distribution of bound proteins restored
            %- for most proteins this means similar number bound as before
            %- except for SMC which should have the bound density preserved
            [~, vals] = find(c.monomerBoundSites);
            boundMonomers = reshape(histc(vals, (1:numel(pm.boundIndexs))'), [], 1);
            [~, vals] = find(c.complexBoundSites);
            boundComplexs = reshape(histc(vals, (1:numel(pc.boundIndexs))'), [], 1);
            
            assertElementsAlmostEqual(...
                exp(log(2) * iter / time.replicationInitiationDuration) * initBoundMonomers, ...
                boundMonomers, ...
                'relative', 0.5, 4);
            assertElementsAlmostEqual(...
                exp(log(2) * iter / time.replicationInitiationDuration) * ...
                initBoundComplexs(setdiff(1:end, [repInit.enzymeGlobalIndexs; pc.rnaPolymeraseIndexs])), ...
                boundComplexs(setdiff(1:end, [repInit.enzymeGlobalIndexs; pc.rnaPolymeraseIndexs])), ...
                'relative', 0.75, 4);
            assertElementsAlmostEqual(...
                exp(log(2) * iter / time.replicationInitiationDuration) * ...
                sum(initBoundComplexs(pc.rnaPolymeraseIndexs)), ...
                sum(boundComplexs(pc.rnaPolymeraseIndexs)), ...
                'relative', 0.75, 4);
            
            %DnaA 7mer broken down and regenerated, and not polymerized
            %(too much)
            assertElementsAlmostEqual(.45^((replicationFinishedIter - replicationInitiationFinishedIter) / 3600) * (4 * 7 + 1), ...
                repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ADP) + repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_1mer_ADP), ...
                'relative', 0.3, 5);
            assertAllEqual(0, repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ATP(2:end)));
            assertAllEqual(0, repInit.enzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(2:end)));
            assertAllEqual(0, repInit.boundEnzymes(repInit.enzymeIndexs_DnaA_Nmer_ADP(2:end)));
        end
    end
    
    methods
        function seedRandStream(this, seed)
            sim = this.simulation;
            
            if sim.verbosity > 0
                fprintf('Seed %d\n', seed);
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
