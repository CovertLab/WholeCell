% Calculates minimal gene complete required for life
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Created: 10/24/2012
classdef MinimalGenes
    properties (Constant = true)
        UNCHARACTERIZED = -1;
        ALWAYS_ESSENTIAL = 0;
        SYNTHETICALLY_ESSENTIAL = 1;
        NEVER_ESSENTIAL = -2;
    end
    
    methods (Static = true)        
        function run(outputDirectory)
            import edu.stanford.covert.cell.sim.analysis.MinimalGenes;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.FitConstants;
            
            %% process arguments
            %default output directory
            if nargin < 1
                outputDirectory = 'output/runAnalysisTests/MinimalGenes';
            end
            
            %create output directory if neccessary
            if ~exist(outputDirectory, 'dir')
                mkdir(outputDirectory)
            end
            
            %% load simulation
            sim = CachedSimulationObjectUtil.load();
            g = sim.gene;
            pm = sim.state('ProteinMonomer');
            metabolism = sim.process('Metabolism');
            
            isGeneMetabolicEnzyme = ...
                any(metabolism.substrateGeneComposition(), 2) | ...
                any(metabolism.enzymeGeneComposition(), 2);
            
            sim.initializeState();
            
            %% single-gene analysis from Karr et al., 2012 Table S2G
            [~, ~, tmp] = xlsread('documentation/paper/tableS2.xls', 'S2G-Gene disruption strains');
            tmp = tmp(cellfun(@(x) ~isempty(x) && ~all(isnan(x)), tmp(:, 1)), :);
            [tfs, idxs] = ismember(tmp(:, 1), g.wholeCellModelIDs);
            geneData = cell(numel(g.wholeCellModelIDs), size(tmp, 2));
            geneData(idxs(tfs), :) = tmp(tfs, :);
            isGeneEssential = strcmp(geneData(:, 5), 'Y');
            isGeneEssential(1:300) = true;
            isGeneEssential(g.getIndexs({'MG_023', 'MG_041', 'MG_049', 'MG_051', 'MG_069', 'MG_111', 'MG_215', 'MG_429'})) = false;
            isGeneCharacterized = cellfun(@(x) ~all(isnan(x)), geneData(:, 3));
            
            geneClasses = NaN(size(g.wholeCellModelIDs, 1), 1);
            geneClasses(~isGeneCharacterized) = MinimalGenes.UNCHARACTERIZED;
            geneClasses(isGeneEssential) = MinimalGenes.ALWAYS_ESSENTIAL;
            geneClasses(isGeneCharacterized & ~isGeneEssential & ~isGeneMetabolicEnzyme) = MinimalGenes.NEVER_ESSENTIAL;
            
            %% calculate metabolic synthetic gene lethality
            %setup problem
            fbaCandidateGeneTfs = isGeneCharacterized & ~isGeneEssential & isGeneMetabolicEnzyme;
            fbaCandidateGeneIdxs = find(fbaCandidateGeneTfs);
            
            m = metabolism;
            [~, S] = m.calcEffectiveFBANetwork(0); %TODO
            S(:, m.fbaReactionIndexs_biomassProduction) = S(:, m.fbaReactionIndexs_biomassProduction) * 1e-5;
            S(m.fbaSubstrateIndexs_biomass, m.fbaReactionIndexs_biomassProduction) = 1;
            tmp = m.substrateGeneComposition();
            [i, ~] = ind2sub(size(m.substrates), m.substrateIndexs_fba);
            subGeneComp = zeros(size(tmp, 1), size(S, 1));
            subGeneComp(:, m.fbaSubstrateIndexs_substrates) = tmp(:, i);
            enzGeneComp = m.enzymeGeneComposition();
            
            nFbaGene = numel(fbaCandidateGeneIdxs);
            nRxn = size(m.reactionStoichiometryMatrix, 2);
            nEnz = size(m.fbaReactionCatalysisMatrix, 2);
            
            fitter = FitConstants(sim);
            [rnas, mons, cpxs, ~, ~, ~, totMons] = fitter.calcMacromolecularCounts(fitter.constructParameterVectorFromSimulation());
            totProtWt = totMons' * pm.molecularWeights(pm.matureIndexs);
            monTfs = isGeneEssential(g.mRNAIndexs);
            essProtWt = totMons(monTfs, 1)' * pm.molecularWeights(pm.matureIndexs(monTfs));
            monTfs = ismember(geneClasses(g.mRNAIndexs), [MinimalGenes.UNCHARACTERIZED; MinimalGenes.NEVER_ESSENTIAL]);
            nonEssProtWts = [
                totMons(monTfs, 1) .* pm.molecularWeights(pm.matureIndexs(monTfs))
                zeros(sum(ismember(geneClasses, [MinimalGenes.UNCHARACTERIZED; MinimalGenes.NEVER_ESSENTIAL])) - sum(monTfs), 1)
                ];
            
            enzCnts = zeros(nEnz, 1); %m.enzymes;
            enzCnts(m.enzymeRNALocalIndexs) = rnas(m.enzymeRNAGlobalIndexs);
            enzCnts(m.enzymeMonomerLocalIndexs) = mons(m.enzymeMonomerGlobalIndexs);
            enzCnts(m.enzymeComplexLocalIndexs) = cpxs(m.enzymeComplexGlobalIndexs);
            metCnts = m.substrates;
            fbaFluxBounds = m.calcFluxBounds(metCnts, enzCnts, ...
                m.fbaReactionBounds, m.fbaEnzymeBounds, ...
                true, true, true, true, true, true);
            fbaFluxBounds(:, 1) = max(fbaFluxBounds(:, 1), -m.realmax);
            fbaFluxBounds(:, 2) = min(fbaFluxBounds(:, 2),  m.realmax);
            
            %find gene lethality sets
            monTfs = fbaCandidateGeneTfs(g.mRNAIndexs);
            metNonEssGeneWts = totMons(monTfs) .* pm.molecularWeights(pm.matureIndexs(monTfs));
            
            opts = m.linearProgrammingOptions;
            [growths, incGenes, fbaFluxs, status] = MinimalGenes.calcMaxGeneSetGrowthRates(...
                S, fbaFluxBounds, ...
                subGeneComp(fbaCandidateGeneIdxs, :), ...
                enzGeneComp(fbaCandidateGeneIdxs, :), ...
                m.fbaReactionCatalysisMatrix, ...
                m.fbaReactionIndexs_biomassProduction, ...
                metNonEssGeneWts, ...
                opts);
            fluxs = zeros(nFbaGene+1, nRxn);
            fluxs(:, m.reactionIndexs_fba) = fbaFluxs(:, m.fbaReactionIndexs_metabolicConversion);
            
%             %enumerate all set
%             rxnGeneMat = m.fbaReactionCatalysisMatrix * enzGeneComp(fbaCandidateGeneIdxs, :)';
%             opts = m.linearProgrammingOptions;
%             opts.solver = 'cplex';
%             [growths, geneSets] = MinimalGenes.enumerateViableGeneSets(...
%                 true(1, numel(fbaCandidateGeneIdxs)), rxnGeneMat, ...
%                 fbaObj, S, m.fbaRightHandSide, ...
%                 fbaFluxBounds(:, 1), fbaFluxBounds(:, 2), ...
%                 opts);
            
            %% analyze
            for i = 1:nFbaGene
                geneClasses(fbaCandidateGeneIdxs(i)) = find(growths > 0 & incGenes(:, i), 1, 'last') - 1;
            end
            
            nGene = numel(g.wholeCellModelIDs);
            nUchar = sum(geneClasses == MinimalGenes.UNCHARACTERIZED);
            nEss = sum(geneClasses == MinimalGenes.ALWAYS_ESSENTIAL);
            nNonEss = sum(geneClasses == MinimalGenes.NEVER_ESSENTIAL);
            nSynEss = nGene - nUchar - nEss - nNonEss;
            
            corrGrowths = zeros(size(growths));
            for i = 1:numel(growths)
                corrGrowths(i) = growths(i) * totProtWt / (incGenes(i, :) * metNonEssGeneWts + essProtWt);
            end
            
            xData = [
                0
                nEss+(0:nSynEss)'
                nEss+nSynEss+(1:nNonEss+nUchar)'
                ];
            yData = [
                0
                corrGrowths
                growths(end) * totProtWt ./ (sum(metNonEssGeneWts) + sort(nonEssProtWts) + essProtWt);
                ];
            plot(xData, yData);
            
            xlim([-0.5 nGene + 0.5]);
            xlim([-0.5+nEss nEss+nSynEss+0.5]);
            ylim([0 max(corrGrowths)*1.05]);
            
            line(nEss * [1 1], ylim(), 'LineStyle', ':', 'Color', 0.5 * [1 1 1]);
            line((nEss + nSynEss) * [1 1], ylim(),  'LineStyle', ':', 'Color', 0.5 * [1 1 1]);
            line((nEss + nSynEss + nNonEss) * [1 1], ylim(),  'LineStyle', ':', 'Color', 0.5 * [1 1 1]);
            
            y = max(ylim());
            text(nEss/2,                              y, 'Ess',       'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            text(nEss + nSynEss/2,                    y, 'Synth Ess', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            text(nEss + nSynEss + nNonEss/2,          y, 'Non-ess',   'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            text(nEss + nSynEss + nNonEss + nUchar/2, y, 'Unchar',    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            
            xlabel('No. genes', 'FontSize', 12);
            ylabel('Growth (fg h^{-1})', 'FontSize', 12);
            
            box('on');
            
            saveas(gcf, [outputDirectory filesep 'Summary.pdf']);
        end
        
        function [growths, incGenes, fluxs, status] = calcMaxGeneSetGrowthRates(...
                S, fluxBounds, subGeneComp, enzGeneComp, rxnEnzCatMat, growthRxnIdx, geneProtMolWts, opts)
            import edu.stanford.covert.util.ComputationUtil;
            
            %setup
            nMet = size(S, 1);
            nRxn = size(S, 2);
            nGene = size(enzGeneComp, 1);
            nEnz = size(enzGeneComp, 2);
            subGeneComp = min(1, subGeneComp);
            enzGeneComp = min(1, enzGeneComp);
            catalyzedRxnIdxs = find(any(rxnEnzCatMat(:, any(enzGeneComp, 1)), 2));
            
            fluxBounds(:, 1) = max(fluxBounds(:, 1), -1e6);
            fluxBounds(:, 2) = min(fluxBounds(:, 2),  1e6);
            
            if isempty(geneProtMolWts)
                geneProtMolWts = ones(nGene, 1);
            end
            
            f = zeros(nRxn + nGene, 1);
            f(growthRxnIdx) = 1;
            g = [zeros(nRxn, 1); geneProtMolWts];
            
            lb = [fluxBounds(:, 1); zeros(nGene, 1)];
            ub = [fluxBounds(:, 2); ones(nGene, 1)];
            %lb(catalyzedRxnIdxs) = -1e6;
            %ub(catalyzedRxnIdxs) =  1e6;
            varTypes = [repmat('C', nRxn, 1); repmat('B', nGene, 1)];
            
            rxnGeneCatMat = (S~=0)' * subGeneComp' | rxnEnzCatMat * enzGeneComp';
            [iRxn, iGene] = find(rxnGeneCatMat);
            nRxnGeneConstraints = numel(iRxn);
            Aineq = [
                zeros(2 * nRxnGeneConstraints, nRxn)  zeros(2 * nRxnGeneConstraints, nGene)
                ];
            Aineq(sub2ind(size(Aineq), (1:nRxnGeneConstraints)', iRxn)) = -1;
            Aineq(sub2ind(size(Aineq), (1:nRxnGeneConstraints)', nRxn + iGene)) = fluxBounds(iRxn, 1);
            Aineq(sub2ind(size(Aineq), (1:nRxnGeneConstraints)' + nRxnGeneConstraints, iRxn)) = 1;
            Aineq(sub2ind(size(Aineq), (1:nRxnGeneConstraints)' + nRxnGeneConstraints, nRxn + iGene)) = -fluxBounds(iRxn, 2);
            bineq = zeros(2 * nRxnGeneConstraints, 1);
            
            Aeq = [
                S               zeros(nMet, nGene)
                zeros(1, nRxn)  ones(1, nGene)
                ];
            beq = zeros(nMet + 1, 1);
            constraintTypes = [repmat('U', 2 * nRxnGeneConstraints, 1);  repmat('S', nMet + 1, 1)];
            
            A = [Aineq; Aeq];
            b = [bineq; beq];
            
            %solve
            growths = NaN(nGene + 1, 1);
            incGenes = zeros(nGene + 1, nGene);
            fluxs = zeros(nGene + 1, nRxn);
            status = NaN(nGene + 1, 2);
            for N = 0:nGene
                beq(end) = N;
                b(end) = N;
                
                %find maximal growth for gene set of size N
                opts.solver = 'cplex';
                [x, ~, ~, errFlag, ~, extra] = ComputationUtil.linearProgramming(...
                    'maximize', f, A, b, lb, ub, constraintTypes, varTypes, opts);
                status(N + 1, 1) = extra.status;
                
                if errFlag || extra.status == 102
                    opts.solver = 'glpk';
                    opts.solverOptions.glpk.tmlim = 60;
                    [y, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                        'maximize', f, A, b, lb, ub, constraintTypes, varTypes, opts);
                    if errFlag
                        warning('WholeCell:warning', 'Linear programming error at N=%d: %s', N, errMsg);
                        continue;
                    end
                    x = y;
                end
                
                growth = x(growthRxnIdx);
                
                %find gene set of size N with growth = maxima and minimal
                %protein cost
                tmpLb = lb;
                tmpUb = ub;
                tmpLb(growthRxnIdx) = growth;
                tmpUb(growthRxnIdx) = growth;
                opts.solver = 'cplex';
                [x, ~, ~, errFlag, ~, extra] = ComputationUtil.linearProgramming(...
                    'minimize', g, A, b, tmpLb, tmpUb, constraintTypes, varTypes, opts);
                status(N + 1, 2) = extra.status;
                
                if errFlag || extra.status == 102
                    opts.solver = 'glpk';
                    opts.solverOptions.glpk.tmlim = 60;
                    [y, ~, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                        'maximize', f, A, b, lb, ub, constraintTypes, varTypes, opts);
                    if errFlag
                        warning('WholeCell:warning', 'Linear programming error at N=%d: %s', N, errMsg);
                        continue;
                    end
                    x = y;
                end
                
                %record output
                growths(N + 1) = growth;
                incGenes(N + 1, :) = x(nRxn+1:end);
                fluxs(N + 1, :) = x(1:nRxn);
            end
        end
        
        function [growths, geneSets] = enumerateViableGeneSets(geneSet, rxnGeneMat, f, S, b, lb, ub, opts, depth)
            import edu.stanford.covert.cell.sim.analysis.MinimalGenes;
            import edu.stanford.covert.util.ComputationUtil;
            
            if nargin < 9
                depth = [];
            end
            
            geneIdxs = find(geneSet);
            nGenes = numel(geneIdxs);
            
            growths = zeros(0, 1);
            geneSets = false(0, size(rxnGeneMat, 2));            
            for i = 1:nGenes
                fprintf('%2d ', [depth geneIdxs(i)]);
                fprintf('\n');
                
                tmpGeneSet = geneSet;
                tmpGeneSet(1, geneIdxs(i)) = false;
                
                tmpLb = lb .* ~any(rxnGeneMat(:, ~tmpGeneSet), 2);
                tmpUb = ub .* ~any(rxnGeneMat(:, ~tmpGeneSet), 2);
                
                if strcmp(opts.solver, 'cplex')
                    [~, fopt, exitflag, output] =  cplexlp(...
                        -f, ...
                        [], [], ...
                        S, b, ...
                        tmpLb, tmpUb);
                    if exitflag ~= 1
                        warning('WholeCell:warning', 'Linear programming error at N=%d: %s', N, output.message);
                        continue;
                    end
                else
                    [~, fopt, ~, errFlag, errMsg] = ComputationUtil.linearProgramming(...
                        'maximize', f, S, b, tmpLb, tmpUb, 'S', 'C', opts);
                    if errFlag
                        warning('WholeCell:warning', 'Linear programming error at N=%d: %s', N, errMsg);
                        continue;
                    end
                end
                growth = -fopt;
                
                if growth > 0
                    [tmp1, tmp2] = MinimalGenes.enumerateViableGeneSets(tmpGeneSet, rxnGeneMat, f, S, b, lb, ub, opts, [depth geneIdxs(i)]);
                    growths = [growths; growth; tmp1]; %#ok<AGROW>
                    geneSets = [geneSets; tmpGeneSet; tmp2]; %#ok<AGROW>
                end
            end
        end
    end
end