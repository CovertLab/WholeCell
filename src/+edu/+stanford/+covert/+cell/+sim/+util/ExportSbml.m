%Requires libSBML: http://sourceforge.net/projects/sbml/files/libsbml/
classdef ExportSbml
    methods (Static)
        function run(sim, outDir)
            import edu.stanford.covert.cell.sim.util.ExportSbml;
            
            %create output directory
            if ~exist(outDir, 'dir')
                mkdir(outDir);
            end
            
            %convert and export models
            abstkineticlawObj = sbioabstractkineticlaw('FirstOrder', 'k*x');
            set (abstkineticlawObj, 'SpeciesVariables', {'x'});
            set (abstkineticlawObj, 'ParameterVariables', {'k'});
            sbioaddtolibrary(abstkineticlawObj);
            
            mdl = ExportSbml.constructTranscriptionSimBiologyModel(sim);
            sbmlexport(mdl, fullfile(outDir, 'Transcription.xml'));
            
            mdl = ExportSbml.constructTranslationSimBiologyModel(sim);
            sbmlexport(mdl, fullfile(outDir, 'Translation.xml'));
            
            [mdl, species, enzymes, reactions] = ExportSbml.constructRNADecaySimBiologyModel(sim);
            ExportSbml.exportModel(mdl, species, enzymes, reactions, fullfile(outDir, 'RNADecay.xml'), fullfile(outDir, 'RNADecay.xls'));            
            
            mdl = ExportSbml.constructProteinDecaySimBiologyModel(sim);
            sbmlexport(mdl, fullfile(outDir, 'ProteinDecay.xml'));
        end
        
        function exportModel(mdl, species, enzymes, reactions, outSbml, outXls)
            sbmlexport(mdl, outSbml);
            xlswrite(outXls, [
                {'ID' 'Name' 'Compartment' 'Copy number, typical' 'Type'}
                [{species.id}' {species.name}' {species.compartment}' {species.copyNumber}' {species.type}']
                ], 'Species');
            xlswrite(outXls, [
                {'ID' 'Name' 'Compartment' 'Copy number, typical' 'Type'}
                [{enzymes.id}' {enzymes.name}' {enzymes.compartment}' {enzymes.copyNumber}' {enzymes.type}']
                ], 'Enzymes');
            xlswrite(outXls, [
                {'ID' 'Name' 'Stoichiometry' 'Enzyme', 'Rate law' 'Rate parameters'}
                [{reactions.id}' {reactions.name}' {reactions.stoichiometry}' {reactions.enzyme}' {reactions.rateLaw}' {reactions.rateParameters}']
                ], 'Reactions');
        end
        
        function mdl = constructTranscriptionSimBiologyModel(sim)
        end
        
        function mdl = constructTranslationSimBiologyModel(sim)
        end
        
        %RNA decay sub-model
        %  Models decay of each RNA species (unmodified, aminoacylated) as
        %  a single reaction with kinetic rate that is linear in the RNA
        %  copy number. The unmodified reactions are catalyzed by
        %  ribonuclease R. The aminoacylated reactions are catalyzed by 
        %  peptidyl-tRNA hydrolase.
        %
        %TODO: factor enzymes into rate law
        function [mdl, species, enzymes, reactions] = constructRNADecaySimBiologyModel(sim)
            %% get handle to submodel
            proc = sim.process('RNADecay');
            rnaState = sim.state('Rna');
            metState = sim.state('Metabolite');
            iCytComp = sim.compartment.cytosolIndexs;
            
            %% create output tables
            species = [];
            enzymes = [];
            reactions = [];
            
            %% create SimBiology model
            mdl = sbiomodel('RNADecay');
            
            %% create initial conditions
            initCond = addvariant(mdl, 'TypicalInitialConditions');
            
            %% create cytosol compartment
            cmp = mdl.addcompartment('c', 'tag', 'cytosol');
            
            %% add metabolite species
            iMet = [
                metState.waterIndexs
                metState.hydrogenIndexs
                metState.nmpIndexs
                metState.aminoAcidIndexs
                ];
            for j = 1:numel(iMet)
                wid = metState.wholeCellModelIDs{iMet(j)};
                name = metState.names{iMet(j)};
                cnt = metState.counts(iMet(j), iCytComp);
                
                cmp.addspecies(wid, 'tag', name);
                addcontent(initCond, {'species', wid, 'InitialAmount', cnt});
                
                species = [species
                    struct('id', wid, 'name', name, 'compartment', 'c', 'copyNumber', cnt, 'type', 'metabolite')
                    ];
            end
            
            %% add RNA species
            %RNA
            for iRna = 1:numel(rnaState.matureIndexs)
                wid = rnaState.wholeCellModelIDs{rnaState.matureIndexs(iRna)};
                name = rnaState.names{rnaState.matureIndexs(iRna)};
                if isempty(name)
                    name = wid;
                end
                cnt = rnaState.counts(rnaState.matureIndexs(iRna), iCytComp);
                
                cmp.addspecies(wid, 'tag', name);
                addcontent(initCond, {'species', wid, 'InitialAmount', cnt});
                
                species = [species
                    struct('id', wid, 'name', name, 'compartment', 'c', 'copyNumber', cnt, 'type', 'rna')
                    ];
            end
            
            %RNA-AA
            for iRna = 1:numel(rnaState.matureTRNAIndexs)
                wid = sprintf('%s_aa', rnaState.wholeCellModelIDs{rnaState.matureIndexs(rnaState.matureTRNAIndexs(iRna))});
                name = sprintf('%s-AA', rnaState.wholeCellModelIDs{rnaState.matureIndexs(rnaState.matureTRNAIndexs(iRna))});
                cnt = rnaState.counts(rnaState.aminoacylatedIndexs(rnaState.matureTRNAIndexs(iRna)), iCytComp);
                
                cmp.addspecies(wid, 'tag', name);
                addcontent(initCond, {'species', wid, 'InitialAmount', cnt});
                
                species = [species
                    struct('id', wid, 'name', name, 'compartment', 'c', 'copyNumber', cnt, 'type', 'rna-aa')
                    ];
            end
            
            %% add enzymes species
            for iEnz = 1:numel(proc.enzymeWholeCellModelIDs)
                wid = proc.enzymeWholeCellModelIDs{iEnz};
                name = proc.enzymeNames{iEnz};
                cnt = proc.enzymes(iEnz);
                
                cmp.addspecies(wid, 'tag', name);
                addcontent(initCond, {'species', wid, 'InitialAmount', cnt});
                
                if ismember(iEnz, proc.enzymeMonomerLocalIndexs)
                    type = 'protein-monomer';
                else
                    type = 'protein-complex';
                end
                
                enzymes = [enzymes
                    struct('id', wid, 'name', name, 'compartment', 'c', 'copyNumber', cnt, 'type', type)
                    ];
            end
            
            %% add decay reactions
            %RNA
            iMet = metState.nmpIndexs;
            for iRna = 1:numel(rnaState.matureIndexs)
                rnaWid = rnaState.wholeCellModelIDs{rnaState.matureIndexs(iRna)};
                rnaName = rnaState.names{rnaState.matureIndexs(iRna)};
                if isempty(rnaName)
                    rnaName = rnaWid;
                end
                
                rxnWid = sprintf('RNA_decay_%s', rnaWid);
                rxnName = sprintf('RNA decay (%s)', rnaName);
                
                baseCnts = rnaState.baseCounts(rnaState.processedIndexs(iRna), iMet);
                len = sum(baseCnts);
                stoichiometry = sprintf('1 %s + %d H2O -> %d AMP + %d CMP + %d GMP + %d UMP + %d H', rnaWid, len - 1, baseCnts, len - 1);
                rxn = mdl.addreaction(stoichiometry, 'tag', rxnWid);
                
                kineticLaw = rxn.addkineticlaw('FirstOrder');
                rateParamId = sprintf('k_dcy_%s', rnaWid);
                rateParamVal = rnaState.decayRates(rnaState.matureIndexs(iRna));
                kineticLaw.addparameter(rateParamId, 'Value', rateParamVal); %1/s
                kineticLaw.ParameterVariableNames = {rateParamId};
                kineticLaw.SpeciesVariableNames = {rnaWid};
                
                reactions = [reactions; struct(...
                    'id', rxnWid, ...
                    'name', rxnName, ...
                    'stoichiometry', stoichiometry, ...
                    'enzyme', 'MG_104_MONOMER', ...
                    'rateLaw', sprintf('%s * %s', rateParamId, rnaWid), ...
                    'rateParameters', sprintf('%s = %f', rateParamId, rateParamVal))
                    ];
            end
            
            %RNA-AA
            iMet = metState.nmpIndexs;
            for iRna = 1:numel(rnaState.matureTRNAIndexs)
                rnaWid = rnaState.wholeCellModelIDs{rnaState.matureIndexs(rnaState.matureTRNAIndexs(iRna))};
                rnaName = rnaState.wholeCellModelIDs{rnaState.matureIndexs(rnaState.matureTRNAIndexs(iRna))};
                rnaAaWid = sprintf('%s_aa', rnaWid);
                rnaAAName = sprintf('%s-AA', rnaName);
                
                rxnWid = sprintf('RNA_decay_%s', rnaAaWid);
                rxnName = sprintf('RNA decay, aminoacylated (%s)', rnaName);
                
                iAA = find(rnaState.baseCounts(rnaState.aminoacylatedIndexs(rnaState.matureTRNAIndexs(iRna)), metState.aminoAcidIndexs));
                aaWid = metState.wholeCellModelIDs{metState.aminoAcidIndexs(iAA)};
                stoichiometry = sprintf('1 %s + 1 H2O -> 1 %s + 1 %s + 1 H', rnaAaWid, rnaWid, aaWid);
                rxn = mdl.addreaction(stoichiometry, 'tag', rxnWid);
                
                kineticLaw = rxn.addkineticlaw('FirstOrder');
                rateParamId = sprintf('k_dcy_%s', rnaAaWid);
                rateParamVal = rnaState.decayRates(rnaState.matureIndexs(rnaState.matureTRNAIndexs(iRna)));
                kineticLaw.addparameter(rateParamId, 'Value', rateParamVal); %1/s
                kineticLaw.ParameterVariableNames = {rateParamId};
                kineticLaw.SpeciesVariableNames = {rnaAaWid};
                
                reactions = [reactions; struct(...
                    'id', rxnWid, ...
                    'name', rxnName, ...
                    'stoichiometry', stoichiometry, ...
                    'enzyme', 'MG_083_MONOMER', ...
                    'rateLaw', sprintf('%s * %s', rateParamId, rnaWid), ...
                    'rateParameters', sprintf('%s = %f', rateParamId, rateParamVal))
                    ];
            end
        end
        
        %RNA decay sub-model
        %  Models decay of each RNA species (unmodified, aminoacylated) as
        %  a single reaction with kinetic rate that is linear in the RNA
        %  copy number. The unmodified reactions are catalyzed by
        %  ribonuclease R. The aminoacylated reactions are catalyzed by 
        %  peptidyl-tRNA hydrolase.
        %
        %TODO: factor enzymes into rate law
        function [mdl, species, enzymes, reactions] = constructProteinDecaySimBiologyModel(sim)
            %% get handle to submodel
            proc = sim.process('ProteinDecay');
            gene = sim.gene;
            rnaState = sim.state('Rna');
            pmState = sim.state('ProteinMonomer');
            pcState = sim.state('ProteinComplex');
            metState = sim.state('Metabolite');            
            iCytComp = sim.compartment.cytosolIndexs;
            iMemComp = sim.compartment.membraneIndexs;
            iExtComp = sim.compartment.extracellularIndexs;
            iTermOrgCytComp = sim.compartment.terminalOrganelleCytosolIndexs;
            iTermOrgMemComp = sim.compartment.terminalOrganelleMembraneIndexs;
            
            %% create output tables
            species = [];
            enzymes = [];
            reactions = [];
            
            %% create SimBiology model
            mdl = sbiomodel('ProteinDecay');
            
            %% create initial conditions
            initCond = addvariant(mdl, 'TypicalInitialConditions');
            
            %% create cytosol compartment
            cytCmp = mdl.addcompartment('c', 'tag', 'cytosol');
            memCmp = mdl.addcompartment('m', 'tag', 'membrane');
            extCmp = mdl.addcompartment('e', 'tag', 'extracellular space');
            
            %% add metabolite species
            iMet = [
                metState.waterIndexs
                metState.aminoAcidIndexs
                metState.ntpIndexs(3)
                metState.ndpIndexs(3)
                metState.hydrogenIndexs
                metState.phosphateIndexs
                ];
            for j = 1:numel(iMet)
                wid = metState.wholeCellModelIDs{iMet(j)};
                name = metState.names{iMet(j)};
                cnt = metState.counts(iMet(j), iCytComp);
                
                cytCmp.addspecies(wid, 'tag', name);
                addcontent(initCond, {'species', wid, 'InitialAmount', cnt});
                
                species = [species
                    struct('id', wid, 'name', name, 'compartment', 'c', 'copyNumber', cnt, 'type', 'metabolite')
                    ];
            end
            
            %% add protein species
            %monomers
            monCnts = sum(...
                + pmState.counts(pmState.matureIndexs, :) ...
                + pmState.counts(pmState.inactivatedIndexs, :) ...
                + pmState.counts(pmState.boundIndexs, :) ...
                , 2);            
            for iMon = 1:numel(pm.matureIndexs)
                wid = pmState.wholeCellModelIDs{pmState.matureIndexs(iMon)};
                name = pmState.names{pmState.matureIndexs(iMon)};
                cnt = monCnts(iMon);
                
                switch pmState.compartments(pmState.matureIndexs(iMon))
                    case iCytComp, cmp = cytCmp;
                    case iMemComp, cmp = memCmp;
                    case iExtComp, cmp = extCmp;
                    case iTermOrgCytComp, cmp = cytCmp;
                    case iTermOrgMemComp, cmp = memCmp;
                end
                
                cmp.addspecies(wid, 'tag', name);
                addcontent(initCond, {'species', wid, 'InitialAmount', cnt});
                
                species = [species
                    struct('id', wid, 'name', name, 'compartment', cmp.Name, 'copyNumber', cnt, 'type', 'protein-monomer')
                    ];
            end
            
            %complexes
            cpxCnts = sum(...
                + pcState.counts(pcState.matureIndexs, :) ...
                + pcState.counts(pcState.inactivatedIndexs, :) ...
                + pcState.counts(pcState.boundIndexs, :) ...
                , 2);
            for iCpx = 1:numel(pcState.matureIndexs)
                wid = pcState.wholeCellModelIDs{pcState.matureIndexs(iCpx)};
                name = pcState.names{pcState.matureIndexs(iCpx)};
                cnt = cpxCnts(iCpx);
                
                switch pcState.compartments(pcState.matureIndexs(iCpx))
                    case iCytComp, cmp = cytCmp;
                    case iMemComp, cmp = memCmp;
                    case iExtComp, cmp = extCmp;
                    case iTermOrgCytComp, cmp = cytCmp;
                    case iTermOrgMemComp, cmp = memCmp;
                end
                
                cmp.addspecies(wid, 'tag', name);
                addcontent(initCond, {'species', wid, 'InitialAmount', cnt});
                
                species = [species
                    struct('id', wid, 'name', name, 'compartment', cmp.Name, 'copyNumber', cnt, 'type', 'protein-complex')
                    ];
            end
            
            %RNA
            iGene = setdiff(find(any(pcState.proteinComplexComposition(:, :, iCytComp), 2)), gene.mRNAIndexs);
            geneWids = gene.wholeCellModelIDs(iGene);
            [~, iRna] = ismember(geneWids, rnaState.wholeCellModelIDs(rnaState.matureIndexs));
            rnaCnts = rnaState.counts(rnaState.matureIndexs, iCytComp);
            for j = 1:numel(iRna)
                wid = rnaState.wholeCellModelIDs{rnaState.matureIndexs(iRna(j))};
                name = rnaState.names{rnaState.matureIndexs(iRna(j))};
                cnt = rnaCnts(iRna(j));
                
                cytCmp.addspecies(wid, 'tag', name);
                addcontent(initCond, {'species', wid, 'InitialAmount', cnt});
                
                species = [species
                    struct('id', wid, 'name', name, 'compartment', 'c', 'copyNumber', cnt, 'type', 'rna')
                    ];
            end
            
            %% add enzymes species
            for iEnz = 1:numel(proc.enzymeWholeCellModelIDs)
                wid = proc.enzymeWholeCellModelIDs{iEnz};
                name = proc.enzymeNames{iEnz};
                cnt = proc.enzymes(iEnz);
                
                cmp.addspecies(wid, 'tag', name);
                addcontent(initCond, {'species', wid, 'InitialAmount', cnt});
                
                if ismember(iEnz, proc.enzymeMonomerLocalIndexs)
                    type = 'protein-monomer';
                else
                    type = 'protein-complex';
                end
                
                enzymes = [enzymes
                    struct('id', wid, 'name', name, 'compartment', 'c', 'copyNumber', cnt, 'type', type)
                    ];
            end
            
            %% add decay reactions
            %monomers
            iMet = metState.nmpIndexs;
            for iRna = 1:numel(rnaState.matureIndexs)
                rnaWid = rnaState.wholeCellModelIDs{rnaState.matureIndexs(iRna)};
                rnaName = rnaState.names{rnaState.matureIndexs(iRna)};
                if isempty(rnaName)
                    rnaName = rnaWid;
                end
                
                rxnWid = sprintf('RNA_decay_%s', rnaWid);
                rxnName = sprintf('RNA decay (%s)', rnaName);
                
                baseCnts = rnaState.baseCounts(rnaState.processedIndexs(iRna), iMet);
                len = sum(baseCnts);
                stoichiometry = sprintf('1 %s + %d H2O -> %d AMP + %d CMP + %d GMP + %d UMP + %d H', rnaWid, len - 1, baseCnts, len - 1);
                rxn = mdl.addreaction(stoichiometry, 'tag', rxnWid);
                
                kineticLaw = rxn.addkineticlaw('FirstOrder');
                rateParamId = sprintf('k_dcy_%s', rnaWid);
                rateParamVal = rnaState.decayRates(rnaState.matureIndexs(iRna));
                kineticLaw.addparameter(rateParamId, 'Value', rateParamVal); %1/s
                kineticLaw.ParameterVariableNames = {rateParamId};
                kineticLaw.SpeciesVariableNames = {rnaWid};
                
                reactions = [reactions; struct(...
                    'id', rxnWid, ...
                    'name', rxnName, ...
                    'stoichiometry', stoichiometry, ...
                    'enzyme', 'MG_104_MONOMER', ...
                    'rateLaw', sprintf('%s * %s', rateParamId, rnaWid), ...
                    'rateParameters', sprintf('%s = %f', rateParamId, rateParamVal))
                    ];
            end
            
            %complexes
            iMet = metState.nmpIndexs;
            for iRna = 1:numel(rnaState.matureTRNAIndexs)
                rnaWid = rnaState.wholeCellModelIDs{rnaState.matureIndexs(rnaState.matureTRNAIndexs(iRna))};
                rnaName = rnaState.wholeCellModelIDs{rnaState.matureIndexs(rnaState.matureTRNAIndexs(iRna))};
                rnaAaWid = sprintf('%s_aa', rnaWid);
                rnaAAName = sprintf('%s-AA', rnaName);
                
                rxnWid = sprintf('RNA_decay_%s', rnaAaWid);
                rxnName = sprintf('RNA decay, aminoacylated (%s)', rnaName);
                
                iAA = find(rnaState.baseCounts(rnaState.aminoacylatedIndexs(rnaState.matureTRNAIndexs(iRna)), metState.aminoAcidIndexs));
                aaWid = metState.wholeCellModelIDs{metState.aminoAcidIndexs(iAA)};
                stoichiometry = sprintf('1 %s + 1 H2O -> 1 %s + 1 %s + 1 H', rnaAaWid, rnaWid, aaWid);
                rxn = mdl.addreaction(stoichiometry, 'tag', rxnWid);
                
                kineticLaw = rxn.addkineticlaw('FirstOrder');
                rateParamId = sprintf('k_dcy_%s', rnaAaWid);
                rateParamVal = rnaState.decayRates(rnaState.matureIndexs(rnaState.matureTRNAIndexs(iRna)));
                kineticLaw.addparameter(rateParamId, 'Value', rateParamVal); %1/s
                kineticLaw.ParameterVariableNames = {rateParamId};
                kineticLaw.SpeciesVariableNames = {rnaAaWid};
                
                reactions = [reactions; struct(...
                    'id', rxnWid, ...
                    'name', rxnName, ...
                    'stoichiometry', stoichiometry, ...
                    'enzyme', 'MG_083_MONOMER', ...
                    'rateLaw', sprintf('%s * %s', rateParamId, rnaWid), ...
                    'rateParameters', sprintf('%s = %f', rateParamId, rateParamVal))
                    ];
            end
        end
    end
end