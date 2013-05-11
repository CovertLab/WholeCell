%SingleGeneDeletions
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 3/23/2011
classdef SingleGeneDeletions
    properties (Constant = true)
        %A. Non-essential (NON_ESSENTIAL)
        %B. Semi-essential
        %   1. Slow growing (SLOW_GROWING)
        %C. Essential
        %   1. No growth
        %      a. Macromolecule maintenance (NON_GROWING)
        %      b. No macromolecule maintenance, toxic metabolite levels (DECOMPOSING)
        %   2. Decaying growth
        %      a. Non-RNA, Non-protein synthesizing (DECAYING_GROWTH_NON_RNA_PROTEIN)
        %         i.  Non-transcribing
        %         ii. Non-maturing
        %      b. Non-protein synthesizing (DECAYING_GROWTH_NON_PROTEIN)
        %   3. Non-dividing
        %      a. Non-replicative (NON_REPLICATIVE)
        %      b. Non-fissive (NON_FISSIVE)
        %   4. Unable to sustain division over many generations (NON_PERPETUATING)
        %   5. No terminal organelle (NO_TERMINAL_ORGANELLE)
        %   6. Toxin accumulation (TOXIN_ACCUMULATION)
        DELETION_STRAIN_CLASSES = {
            'Non-essential'             'Doesn''t meet criterion of any of the other categories'
            'Decomposing'               'No growth, no energy production, no macromolecule maintenance, toxic metabolite levels'
            'Non-growing'               'No growth, energy production, macromolecule maintenance'
            'Non-RNA synthesizing'      'Derivative of growth is negative, time to mass doubling is long, time to end of replication is long, cell cycle is elongated'
            'Non-protein synthesizing'  'Derivative of growth is negative, time to mass doubling is long, time to end of replication is long, cell cycle is elongated'
            'Non-replicative'           'End ploidy = 1, high dNTPs, no simulation finishes replication (rep end time = NaN)'
            'Non-fissive'               'End ploidy = 2, but pinchedDiamter = diameter'
            'Slow growing'              'Grows and divides, but slower than wild-type'
            'Toxin accumulation'        'Accumulates toxic concentrations of metabolites'
            'Non-perpetuating'          'Gene is essential, but cells will survive for several generations until gene products are sufficiently diluted among progeny'
            'No terminal organelle'     'Lack terminal organelle'
            'Unobserved'                'Simulation not yet run'
            'Wild-type'                 'Wild-type simulation'
            };
        NON_ESSENTIAL                   = 1
        DECOMPOSING                     = 2
        NON_GROWING                     = 3
        DECAYING_GROWTH_NON_RNA_PROTEIN = 4
        DECAYING_GROWTH_NON_PROTEIN     = 5
        NON_REPLICATIVE                 = 6
        NON_FISSIVE                     = 7
        SLOW_GROWING                    = 8
        TOXIN_ACCUMULATION              = 9
        NON_PERPETUATING                = 10
        NO_TERMINAL_ORGANELLE           = 11
        UNOBSERVED                      = 12
        WILD_TYPE                       = 13
        
        COLORS = [
            0 0 1 %blue
            1 0 0 %red
            0 1 0 %green
            0 1 1 %cyan
            1 0.5 0 %orange
            0.5 0 0.5 %purple
            1 0.84 0 %gold
            0.25 0 0.5 %indigo
            0 1 1 %aqua
            0.5 1 0.8314 %aquamarine
            0.329 0.525 0.043 %dark gold
            1 0.753 0.795 %pink
            0.25 0.25 0.25 %dark grey
            ];
        
        nonPerpetuatingEssentialGenes = {
            'MG_0001'
            'MG_012'
            'MG_019'
            'MG_048'
            'MG_072'
            'MG_106'
            'MG_109'
            'MG_110'
            'MG_139'
            'MG_143'
            'MG_170'
            'MG_172'
            'MG_184'
            'MG_210'
            'MG_277'
            'MG_297'
            'MG_305'
            'MG_329'
            'MG_335'
            'MG_384'
            'MG_387'
            'MG_392'
            'MG_393'
            'MG_425'
            'MG_442'
            'MG_464'
            'MG_476'
            };
    end
    
    %printing
    methods (Static)
        function run(simBatchDir, useCachedData, iJob, nJobs)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 2
                useCachedData = false;
            end
            if nargin < 4
                iJob = 1;
                nJobs = 1;
            end
            
            %% constants
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            sim = CachedSimulationObjectUtil.load();
            
            %% plot overview of each single gene deletion strain
            SingleGeneDeletions.plotIndividualDeletionStrains(simBatchDir, useCachedData, iJob, nJobs);
            
            %% analyze specific single gene deletion strains
            SingleGeneDeletions.analyzeSpecificDeletionStrains(useCachedData);
            
            %% calculate strain statistics and classify strains
            SingleGeneDeletions.calcGeneDeletionStrainStats(useCachedData);
            SingleGeneDeletions.classifyIndividualDeletionStrains(useCachedData);
            
            %% plot overview of all single gene deletion strains
            SingleGeneDeletions.plotOverview();
            SingleGeneDeletions.plotAllDeletionStrains();
            
            %% plot signature of each deletion strain class
            SingleGeneDeletions.plotDeletionStrainClasses();
            SingleGeneDeletions.plotDeletionStrainClassification();
            SingleGeneDeletions.plotDeletionGrowthRateDistribution();
            
            %% Metabolic essentiality
            [deletionMetabolicOutputs, wtMetabolicOutput] = SingleGeneDeletions.calculateSingleGeneDeletionMetabolicOutputs(sim);
            
            %% plot summary of model accuracy
            load([baseDir 'wtStats.mat']);
            load([baseDir 'deletionStats.mat']);
            load([baseDir 'geneClasses.mat']);
            
            SingleGeneDeletions.plot(sim, deletionMetabolicOutputs, wtMetabolicOutput, deletionStats, wtStats, geneClasses, [baseDir 'summaryGrid.svg']);
            [~, ~] = system(sprintf('inkscape "%ssummaryGrid.svg" --export-pdf="%ssummaryGrid.pdf" --export-area-page', baseDir, baseDir));
            
            %% tables
            if ispc && exist([baseDir 'summary.xls'], 'file')
                delete([baseDir 'summary.xls']);
            end
            
            %classes
            [classContent, classColLabels] = SingleGeneDeletions.printSingleGeneDeletionStrainClasses(sim, geneClasses);
            PrintUtil.printToFile(classContent, classColLabels, [baseDir 'summary.xls'], 'Strains Classes');
            
            %strains
            [strainContent, strainColLabels] = SingleGeneDeletions.printSingleGeneDeletionStrains(sim, deletionMetabolicOutputs, wtMetabolicOutput, deletionStats, wtStats, geneClasses);
            PrintUtil.printToFile(strainContent, strainColLabels, [baseDir 'summary.xls'], 'Strains');
        end
        
        function plotIndividualDeletionStrains(simBatchDir, useCachedData, iJob, nJobs)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.PdfUtil;
            
            if nargin < 4
                iJob = 1;
                nJobs = 1;
            end
            
            %% constants
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            
            simDir = SimulationDiskUtil.getSimulation([simBatchDir filesep '1']);
            idxs = find(simDir == '/' | simDir == '\');
            simBatch = simDir(idxs(end-1)+1:idxs(end)-1);
            
            sim = CachedSimulationObjectUtil.load();
            g = sim.gene;
            
            %% overview of each single gene deletion strain
            
            %find simulation with single gene deletion data
            [metaDataAll, options] = SimulationDiskUtil.getSimulations();
            
            tfKO = false(size(options));
            geneticKnockouts = zeros(size(options));
            for i = 1:numel(options)
                tfKO(i) = ...
                    ~isempty(options(i).geneticKnockouts) && ...
                    isscalar(options(i).geneticKnockouts) && ...
                    isempty(options(i).media) && ...
                    isempty(options(i).stimulus);
                if tfKO(i)
                    geneticKnockouts(i) = find(strcmp(options(i).geneticKnockouts{1}, g.wholeCellModelIDs), 1);
                end
            end
            
            %find simulations with wild-type data and get data
            tfWt = ~tfKO;
            metaData = metaDataAll(tfWt);
            
            %get wild-type data
            if (nargin >= 2 && useCachedData) || ~any(strcmp({metaData.simGroup}, simBatch))
                %load cached data
                try %#ok<TRYNC>
                    wtData = load([baseDir 'WT.mat']);
                    metaData = wtData.metaData;
                end
            end
            
            if ~exist('wtData', 'var')
                %get wild-type data
                [wtData, metaData] = SingleGeneDeletions.getSimulationData(sim, metaData);
                try %#ok<TRYNC>
                    save([baseDir 'WT.mat'], '-struct', 'wtData');
                    save([baseDir 'WT.mat'], '-append', 'metaData');
                end
            end
            
            %plot wild-type data
            if any(strcmp({metaData.simGroup}, simBatch)) || ~exist([baseDir 'WT.pdf'], 'file')
                [~, figHandle] = PlotUtil.newAxesHandle();
                SingleGeneDeletions.plotSingleGeneDeletion({wtData}, 1, {'Wild-type'}, false, figHandle);
                try %#ok<TRYNC>
                    saveas(figHandle, [baseDir 'WT.pdf']);
                end
                close(figHandle);
                PdfUtil.rotatePages([baseDir 'WT.pdf']);
            end
            
            %get and plot single-gene deletion data
            geneIdxs = setdiff(unique(geneticKnockouts), 0);
            for i = iJob:nJobs:numel(geneIdxs)
                geneIdx = geneIdxs(i);
                geneID = g.wholeCellModelIDs{geneIdx};
                
                if nargin >= 2 && useCachedData
                    %load cached data
                    try %#ok<TRYNC>
                        deletionData = load([baseDir geneID '.mat']);
                        metaData = deletionData.metaData;
                        if ~any(strcmp({metaData.simGroup}, simBatch))
                            continue;
                        end
                    end
                end
                if ~exist('deletionData', 'var')
                    %find simulations
                    simTfs = geneticKnockouts == geneIdx;
                    metaData = metaDataAll(simTfs);
                    if ~any(strcmp({metaData.simGroup}, simBatch))
                        continue;
                    end
                    
                    %load data
                    [deletionData, metaData] = SingleGeneDeletions.getSimulationData(sim, metaData); %#ok<NASGU>
                    try %#ok<TRYNC>
                        save([baseDir geneID '.mat'], '-struct', 'deletionData');
                        save([baseDir geneID '.mat'], '-append', 'metaData');
                    end
                end
                
                [~, figHandle] = PlotUtil.newAxesHandle();
                SingleGeneDeletions.plotSingleGeneDeletion({deletionData; wtData}, [1; 2], {['{\Delta}' strrep(geneID, '_', '\_')]; 'Wild-Type'}, [false; true], figHandle);
                try
                    saveas(figHandle, [baseDir geneID '.pdf']);
                    PdfUtil.rotatePages([baseDir geneID '.pdf']);
                catch exception
                    warning('WholeCell:warning', 'Unable to save %s plot:\n%s', geneID, exception.getReport());
                end
                close(figHandle);
                
                clear deletionData;
            end
        end
        
        function [deletionStats, wtStats] = calcGeneDeletionStrainStats(useCachedData, iJob, nJobs)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 3
                iJob = 1;
                nJobs = 1;
            end
            
            sim = CachedSimulationObjectUtil.load();
            g = sim.gene;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            
            tmp = dir([baseDir filesep 'MG*.mat']);
            [~, geneticKnockouts] = ismember(cellfun(@(str) str(1:end-4), {tmp.name}', 'UniformOutput', false), g.wholeCellModelIDs);
            geneIdxs = setdiff(unique(geneticKnockouts), 0);
            
            if nargin >= 1 && useCachedData && exist([baseDir 'deletionStats.mat'], 'file')
                load([baseDir 'deletionStats.mat']);
                geneIdxs = intersect(geneIdxs, find(all(all(isnan(deletionStats), 2), 3))); %#ok<NODEF>
            else
                deletionStats = NaN(numel(g.wholeCellModelIDs), 19, 2);
            end
            for i = iJob:nJobs:numel(geneIdxs)
                geneID = g.wholeCellModelIDs{geneIdxs(i)};
                if useCachedData && exist([baseDir geneID '-stats.mat'], 'file')
                    load([baseDir geneID '-stats.mat']);
                else
                    tmp = SingleGeneDeletions.calcStrainStats(geneID);
                    save([baseDir geneID '-stats.mat'], 'tmp');
                end
                deletionStats(geneIdxs(i), :, :) = tmp;
            end
            save([baseDir 'deletionStats.mat'], 'deletionStats');
            
            if nargin < 0 || ~useCachedData || ~exist([baseDir 'wtStats.mat'], 'file')
                wtStats = SingleGeneDeletions.calcStrainStats('WT');
                save([baseDir 'wtStats.mat'], 'wtStats');
            else
                load([baseDir 'wtStats.mat']);
            end
        end
        
        function stats = calcStrainStats(geneID)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            data = load([baseDir geneID '.mat']);
            
            inds = sub2ind(size(data.growth), data.simEndTimes + 1, (1:numel(data.simEndTimes))');
            dt = diff(data.time(find(~isnan(data.growth), 2, 'first')));
            
            stats = NaN(1, 19, 2);
            stats(1,  1, :) = [nanmean(data.growth(inds)) nanstd(data.growth(inds))];
            stats(1,  2, :) = [nanmean(data.atpProduction(inds)) nanstd(data.atpProduction(inds))];
            stats(1,  3, :) = [nanmean(data.gtpProduction(inds)) nanstd(data.gtpProduction(inds))];
            stats(1,  4, :) = [nanmean(data.repInitDuration) nanstd(data.repInitDuration)];
            stats(1,  5, :) = [nanmean(data.repDuration) nanstd(data.repDuration)];
            stats(1,  6, :) = [nanmean(data.cytokinesisDuration) nanstd(data.cytokinesisDuration)];
            stats(1,  7, :) = [nanmean(data.cellCycleDuration) nanstd(data.cellCycleDuration)];
            stats(1,  8, :) = [nanmean(data.massDoublingDuration) nanstd(data.massDoublingDuration)];
            stats(1,  9, :) = [nanmean(nansum(data.atpUsage, 1) * dt, 2) nanstd(nansum(data.atpUsage, 1) * dt, [], 2)]; %atp usage
            stats(1, 10, :) = [nanmean(nansum(data.gtpUsage, 1) * dt, 2) nanstd(nansum(data.gtpUsage, 1) * dt, [], 2)]; %gtp usage
            stats(1, 11, :) = [nanmean(nanmean(data.dnaWt ./ data.cellWt, 1), 2) nanstd(nanmean(data.dnaWt ./ data.cellWt, 1), [], 2)]; %DNA content
            stats(1, 12, :) = [nanmean(nanmean(data.rnaWt ./ data.cellWt, 1), 2) nanstd(nanmean(data.rnaWt ./ data.cellWt, 1), [], 2)]; %RNA content
            stats(1, 13, :) = [nanmean(nanmean(data.proteinWt ./ data.cellWt, 1), 2) nanstd(nanmean(data.proteinWt ./ data.cellWt, 1), [], 2)]; %Protein content
            stats(1, 14, :) = [nanmean(nanmean(data.superhelicalDensity, 1), 2) nanstd(nanmean(data.superhelicalDensity, 1), [], 2)]; %superhelicity
            
            tau_rna = zeros(size(data.simEndTimes));
            tau_protein = zeros(size(data.simEndTimes));
            tau_damagedRna = zeros(size(data.simEndTimes));
            tau_damagedProtein = zeros(size(data.simEndTimes));
            adjrsquare_rna = zeros(size(data.simEndTimes));
            adjrsquare_protein = zeros(size(data.simEndTimes));
            adjrsquare_damagedRna = zeros(size(data.simEndTimes));
            adjrsquare_damagedProtein = zeros(size(data.simEndTimes));
            expFittype = fittype('a * exp(log(2) * t / tau)', ...
                'independent', 't', ...
                'coefficients', {'a', 'tau'});
            for i = 1:numel(data.simEndTimes)
                expFitOptions = fitoptions(...
                    'Method', 'NonlinearLeastSquares', ...
                    'Startpoint', [data.rnaWt(2, i) * 1e15  8.6], ...
                    'Lower', [0 -1e3 * (diff(data.rnaWt([2 data.simEndTimes(i)+1], i)) < 0)], ...
                    'Upper', [max(data.rnaWt(:, i)) * 1e15 1e3]);
                [fitResult, gof] = fit(data.time(2:data.simEndTimes(i)+1, 1), data.rnaWt(2:data.simEndTimes(i)+1, i) * 1e15, expFittype, expFitOptions);
                tau_rna(i) = fitResult.tau;
                adjrsquare_rna(i) = gof.adjrsquare;
                if gof.adjrsquare < 0.5
                    tau_rna(i) = NaN;
                end
                
                expFitOptions = fitoptions(...
                    'Method', 'NonlinearLeastSquares', ...
                    'Startpoint', [data.proteinWt(2, i) * 1e15  8.6], ...
                    'Lower', [0 -1e3 * (diff(data.proteinWt([2 data.simEndTimes(i)+1], i)) < 0)], ...
                    'Upper', [max(data.proteinWt(:, i)) * 1e15 1e3]);
                [fitResult, gof] = fit(data.time(2:data.simEndTimes(i)+1, 1), data.proteinWt(2:data.simEndTimes(i)+1, i) * 1e15, expFittype, expFitOptions);
                tau_protein(i) = fitResult.tau;
                adjrsquare_protein(i) = gof.adjrsquare;
                if gof.adjrsquare < 0.5
                    tau_protein(i) = NaN;
                end
                
                xdata = data.time(2:data.simEndTimes(i)+1, 1) / 3600;
                ydata = data.damagedRnas(2:data.simEndTimes(i)+1, i);
                xdata = xdata(~isnan(ydata));
                ydata = ydata(~isnan(ydata));
                expFitOptions = fitoptions(...
                    'Method', 'NonlinearLeastSquares', ...
                    'Startpoint', [ydata(2, 1)  8.6], ...
                    'Lower', [0 -1e3 * (diff(data.damagedRnas([2 data.simEndTimes(i)+1], i)) < 0)], ...
                    'Upper', [max(max(data.damagedRnas(:, i)), 1) 1e3]);
                [fitResult, gof] = fit(xdata, ydata, expFittype, expFitOptions);
                tau_damagedRna(i) = fitResult.tau;
                adjrsquare_damagedRna(i) = gof.adjrsquare;
                if gof.adjrsquare < 0.5
                    tau_damagedRna(i) = NaN;
                end
                
                xdata = data.time(2:data.simEndTimes(i)+1, 1) / 3600;
                ydata = data.damagedProteins(2:data.simEndTimes(i)+1, i);
                xdata = xdata(~isnan(ydata));
                ydata = ydata(~isnan(ydata));
                expFitOptions = fitoptions(...
                    'Method', 'NonlinearLeastSquares', ...
                    'Startpoint', [ydata(2, 1)  8.6], ...
                    'Lower', [0 -1e3 * (diff(data.damagedProteins([2 data.simEndTimes(i)+1], i)) < 0)], ...
                    'Upper', [max(max(data.damagedProteins(:, i)), 1) 1e3]);
                [fitResult, gof] = fit(xdata, ydata, expFittype, expFitOptions);
                tau_damagedProtein(i) = fitResult.tau;
                adjrsquare_damagedProtein(i) = gof.adjrsquare;
                if gof.adjrsquare < 0.5
                    tau_damagedProtein(i) = NaN;
                end
            end
            stats(1, 15, :) = [nanmean(tau_rna) nanstd(tau_rna)]; %rate of change of RNA content
            stats(1, 16, :) = [nanmean(tau_protein) nanstd(tau_protein)]; %rate of change of protein content
            stats(1, 17, :) = [nanmean(tau_damagedRna) nanstd(tau_damagedRna)]; %rate of change of damaged RNA content
            stats(1, 18, :) = [nanmean(tau_damagedProtein) nanstd(tau_damagedProtein)]; %rate of change of damaged protein content
            
            stats(1, 19, :) = [nanmean(nanmean(data.antibiotics, 1), 2) nanstd(nanmean(data.antibiotics, 1), [], 2)]; %antibiotics count
        end
        
        function geneClasses = classifyIndividualDeletionStrains(useCachedData)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            sim = CachedSimulationObjectUtil.load();
            g = sim.gene;
            
            wtData = load([baseDir 'WT.mat']);
            tmp = dir([baseDir filesep 'MG*.mat']);
            [~, geneticKnockouts] = ismember(cellfun(@(str) str(1:end-4), {tmp.name}', 'UniformOutput', false), g.wholeCellModelIDs);
            geneIdxs = setdiff(unique(geneticKnockouts), 0);
            
            if nargin >= 1 && useCachedData && exist([baseDir 'geneClasses.mat'], 'file')
                load([baseDir 'geneClasses.mat']);
                geneIdxs = intersect(geneIdxs, find(geneClasses == SingleGeneDeletions.UNOBSERVED)); %#ok<NODEF>
            else
                geneClasses = repmat(SingleGeneDeletions.UNOBSERVED, size(g.wholeCellModelIDs));
            end
            
            for i = 1:numel(geneIdxs)
                geneIdx = geneIdxs(i);
                geneID = g.wholeCellModelIDs{geneIdx};
                geneClasses(geneIdx) = SingleGeneDeletions.classifyDeletionStrain(geneID, wtData);
                clear deletionData;
            end
            
            try %#ok<TRYNC>
                save([baseDir 'geneClasses.mat'], 'geneClasses');
            end
        end
        
        function analyzeSpecificDeletionStrains(useCachedData)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            
            wtData = load([baseDir 'WT.mat']);
            
            md = meta.class.fromName('edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions');
            for i = 1:numel(md.Methods)
                if length(md.Methods{i}.Name) >= length('analyzeSingleGeneDeletionStrain_') && ...
                        isequal(md.Methods{i}.Name(1:length('analyzeSingleGeneDeletionStrain_')), 'analyzeSingleGeneDeletionStrain_') && ...
                        isequal(md.Methods{i}.InputNames, {'figHandle'; 'geneID'; 'deletionData'; 'wtData'; 'figData'})
                    geneID = md.Methods{i}.Name(length('analyzeSingleGeneDeletionStrain_')+1:end);
                    if ~exist([baseDir geneID '.mat'], 'file')
                        continue;
                    end
                    
                    deletionData = load([baseDir geneID '.mat']);
                    [~, figHandle] = PlotUtil.newAxesHandle();
                    if useCachedData && exist([baseDir geneID '-analysis.mat'], 'file')
                        figData = load([baseDir geneID '-analysis.mat']);
                        SingleGeneDeletions.(md.Methods{i}.Name)(figHandle, geneID, deletionData, wtData, figData);
                    else
                        figData = SingleGeneDeletions.(md.Methods{i}.Name)(figHandle, geneID, deletionData, wtData); %#ok<NASGU>
                        try %#ok<TRYNC>
                            save([baseDir geneID '-analysis.mat'], '-struct', 'figData');
                        end
                    end
                    clear figData;
                    
                    try %#ok<TRYNC>
                        saveas(figHandle, [baseDir geneID '-analysis.pdf']);
                    end
                    close(figHandle);
                end
            end
        end
        
        function plotAllDeletionStrains(iJob, nJobs)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.PdfUtil;
            
            if nargin < 2
                iJob = 1;
                nJobs = 1;
            end
            
            %% overview of all single gene deletion strains
            sim = CachedSimulationObjectUtil.load();
            g = sim.gene;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            
            load([baseDir 'geneClasses.mat']);
            tmp = dir([baseDir filesep 'MG*.mat']);
            [~, geneticKnockouts] = ismember(cellfun(@(str) str(1:end-4), {tmp.name}', 'UniformOutput', false), g.wholeCellModelIDs);
            geneIdxs = setdiff(unique(geneticKnockouts), 0);
            
            data = cell(numel(geneIdxs)+1, 1);
            for i = 1:numel(geneIdxs)
                geneID = g.wholeCellModelIDs{geneIdxs(i)};
                tmp = load([baseDir geneID '.mat']);
                if numel(tmp.simEndTimes) > 1
                    fields = setdiff(fieldnames(tmp), {'metaData', 'time', 'simEndTimes'});
                    simIdx = randperm(numel(tmp.simEndTimes));
                    simIdx = simIdx(1);
                    tmp.metaData = tmp.metaData(simIdx);
                    tmp.simEndTimes = tmp.simEndTimes(simIdx);
                    for k = 1:numel(fields)
                        tmp.(fields{k}) = tmp.(fields{k})(:, simIdx);
                    end
                end
                data{i} = tmp;
            end
            data{end} = load([baseDir 'WT.mat']);
            
            classLabels = [
                [SingleGeneDeletions.DELETION_STRAIN_CLASSES(geneClasses(geneIdxs)) ...
                cellfun(@(geneIdx) ['{\Delta}' strrep(g.wholeCellModelIDs{geneIdx}, '_', '\_')], num2cell(geneIdxs), 'UniformOutput', false)];
                {'Wild-Type' []}];
            
            [~, figHandle] = PlotUtil.newAxesHandle();
            nPlots = SingleGeneDeletions.plotSingleGeneDeletion(...
                data, ...
                [geneClasses(geneIdxs); SingleGeneDeletions.WILD_TYPE], ...
                classLabels, ...
                [false(numel(geneIdxs), 1); true], ...
                figHandle);
            
            try %#ok<TRYNC>
                saveas(figHandle, [baseDir 'all.pdf']);
                PdfUtil.rotatePages([baseDir 'all.pdf']);
            end
            close(figHandle);
            
            for i = iJob:nJobs:nPlots
                [~, figHandle] = PlotUtil.newAxesHandle();
                SingleGeneDeletions.plotSingleGeneDeletion(...
                    data, ...
                    [geneClasses(geneIdxs); SingleGeneDeletions.WILD_TYPE], ...
                    classLabels, ...
                    [false(numel(geneIdxs), 1); true], ...
                    figHandle, i);
                try %#ok<TRYNC>
                    saveas(figHandle, [baseDir sprintf('all-%02d.pdf', i)]);
                end
                close(figHandle);
            end
        end
        
        function plotDeletionStrainClasses(classIdxs)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.PdfUtil;
            
            if nargin < 1
                classIdxs = 1:size(SingleGeneDeletions.DELETION_STRAIN_CLASSES, 1) - 2;
            end
            
            sim = CachedSimulationObjectUtil.load();
            g = sim.gene;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            load([baseDir 'geneClasses.mat']);
            
            for i = 1:numel(classIdxs)
                geneIdxs = find(geneClasses == classIdxs(i));
                
                data = cell(numel(geneIdxs) + 1, 1);
                for j = 1:numel(geneIdxs)
                    tmp = load([baseDir g.wholeCellModelIDs{geneIdxs(j)} '.mat']);
                    if numel(tmp.simEndTimes) > 1
                        fields = setdiff(fieldnames(tmp), {'metaData', 'time', 'simEndTimes'});
                        simIdx = randperm(numel(tmp.simEndTimes));
                        simIdx = simIdx(1);
                        tmp.metaData = tmp.metaData(simIdx);
                        tmp.simEndTimes = tmp.simEndTimes(simIdx);
                        for k = 1:numel(fields)
                            tmp.(fields{k}) = tmp.(fields{k})(:, simIdx);
                        end
                    end
                    data{j} = tmp;
                end
                data{end} = load([baseDir 'WT.mat']);
                
                classLabels = [
                    cellfun(@(geneIdx) ['{\Delta}' strrep(g.wholeCellModelIDs{geneIdx}, '_', '\_')], ...
                    num2cell(geneIdxs), 'UniformOutput', false);
                    {'Wild-Type'}];
                
                [~, figHandle] = PlotUtil.newAxesHandle();
                SingleGeneDeletions.plotSingleGeneDeletion(...
                    data, ...
                    1:numel(data), ...
                    classLabels, ...
                    [false(numel(geneIdxs), 1); true], ...
                    figHandle, ...
                    [], [PlotUtil.getRedGreenColorOrder(1:numel(geneIdxs)); 0 1 1]);
                if nargin == 0
                    try %#ok<TRYNC>
                        saveas(figHandle, sprintf('%sclass-%d.pdf', baseDir, classIdxs(i)));
                        PdfUtil.rotatePages(sprintf('%sclass-%d.pdf', baseDir, classIdxs(i)));
                    end
                    close(figHandle);
                end
                
                clear data;
            end
        end
        
        function plotDeletionStrainClassification()
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            load([baseDir 'wtStats.mat']);
            load([baseDir 'deletionStats.mat']);
            load([baseDir 'geneClasses.mat']);
            
            %Decomposing, non-growing, slow growing, non-essential
            SingleGeneDeletions.plotDeletionStrainClassificationHelper(...
                deletionStats(:, 1, 1), deletionStats(:, 9, 1) + deletionStats(:, 10, 1), ...
                'Growth (fg h^-1)', 'Energy Usage', ...
                geneClasses, [
                SingleGeneDeletions.DECOMPOSING
                SingleGeneDeletions.NON_GROWING
                SingleGeneDeletions.NON_ESSENTIAL
                SingleGeneDeletions.SLOW_GROWING
                ], 1, 'SouthEast');
            
            %decaying growth
            SingleGeneDeletions.plotDeletionStrainClassificationHelper(...
                deletionStats(:, 1, 1), deletionStats(:, 15, 1), ...
                'Growth (fg h^-1)', 'RNA Mass Decay Constant (h)', ...
                geneClasses, [
                SingleGeneDeletions.DECAYING_GROWTH_NON_RNA_PROTEIN
                SingleGeneDeletions.DECAYING_GROWTH_NON_PROTEIN
                ], 2, 'NorthEast');
            
            %non-replicative
            SingleGeneDeletions.plotDeletionStrainClassificationHelper(...
                deletionStats(:, 1, 1), deletionStats(:, 11, 1), ...
                'Growth (fg h^-1)', 'DNA Mass (fg)', ...
                geneClasses, [
                SingleGeneDeletions.NON_REPLICATIVE
                ], 3, 'SouthWest');
            
            %non-fissive
            ydata = deletionStats(:, 7, 1);
            ydata(isnan(ydata)) = max(ydata);
            SingleGeneDeletions.plotDeletionStrainClassificationHelper(...
                deletionStats(:, 5, 1), ydata, ...
                'Replication Time (h)', 'Cell Cycle Length (h)', ...
                geneClasses, [
                SingleGeneDeletions.NON_FISSIVE
                ], 4, 'SouthEast');
            
            %toxin accumulation
            SingleGeneDeletions.plotDeletionStrainClassificationHelper(...
                deletionStats(:, 1, 1), deletionStats(:, 18, 1), ...
                'Growth (fg h^-1)', 'Damaged Protein Mass Constant (h)', ...
                geneClasses, [
                SingleGeneDeletions.TOXIN_ACCUMULATION
                ], 5, 'NorthEast');
        end
        
        function plotDeletionGrowthRateDistribution(wtGrowth, koGrowth, geneClasses, outputFile)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            sim = CachedSimulationObjectUtil.load();
            
            if nargin < 4
                outputFile = [baseDir 'deletionGrowthRateDistribution.pdf'];
            end
            
            %% get data
            if nargin < 1
                wtData = load([baseDir 'WT.mat'], 'simEndTimes', 'growth');
                wtGrowth = wtData.growth(sub2ind(size(wtData.growth), wtData.simEndTimes, (1:numel(wtData.simEndTimes))'));
            end
            if nargin < 2
                load([baseDir 'deletionStats.mat']);
                koGrowth = deletionStats(:, 1, 1); %#ok<NODEF>
            end
            if nargin < 3
                load([baseDir 'geneClasses.mat']);
            end
            
            koNonEssGrowth = koGrowth(geneClasses == SingleGeneDeletions.NON_ESSENTIAL);
            koSlowGrowth = koGrowth(geneClasses == SingleGeneDeletions.SLOW_GROWING);
            koGrowthEssGrowth = koGrowth(ismember(geneClasses, [
                SingleGeneDeletions.DECOMPOSING;
                SingleGeneDeletions.NON_GROWING
                SingleGeneDeletions.DECAYING_GROWTH_NON_RNA_PROTEIN;
                SingleGeneDeletions.DECAYING_GROWTH_NON_PROTEIN
                ]));
            koOtherEssGrowth = koGrowth(ismember(geneClasses, [
                SingleGeneDeletions.NON_REPLICATIVE;
                SingleGeneDeletions.NON_FISSIVE;
                SingleGeneDeletions.NON_PERPETUATING;
                SingleGeneDeletions.NO_TERMINAL_ORGANELLE
                SingleGeneDeletions.TOXIN_ACCUMULATION
                ]));
            
            wtGrowth =  ((wtGrowth*1e-15)/2) / sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight);
            koNonEssGrowth = ((koNonEssGrowth*1e-15)/2) / sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight);
            koSlowGrowth = ((koSlowGrowth*1e-15)/2) / sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight);
            koGrowthEssGrowth = ((koGrowthEssGrowth*1e-15)/2) / sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight);
            koOtherEssGrowth = ((koOtherEssGrowth*1e-15)/2) / sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight);
            
            edges = linspace(0, max([wtGrowth; koGrowth]), 50);
            freqWtGrowth = histc(wtGrowth, edges) / numel(wtGrowth);
            freqKONonEssGrowth = histc(koNonEssGrowth, edges) / numel(koNonEssGrowth);
            freqKOSlowGrowth = histc(koSlowGrowth, edges) / numel(koSlowGrowth);
            freqKOGrowthEssGrowth = histc(koGrowthEssGrowth, edges) / numel(koGrowthEssGrowth);
            freqKOOtherEssGrowth = histc(koOtherEssGrowth, edges) / numel(koOtherEssGrowth);
            
            freqWtGrowth(freqWtGrowth == 0) = NaN;
            freqKONonEssGrowth(freqKONonEssGrowth == 0) = NaN;
            freqKOSlowGrowth(freqKOSlowGrowth == 0) = NaN;
            freqKOGrowthEssGrowth(freqKOGrowthEssGrowth == 0) = NaN;
            freqKOOtherEssGrowth(freqKOOtherEssGrowth == 0) = NaN;
            
            %% plot
            [~, figHandle] = PlotUtil.newAxesHandle();
            clf(figHandle);
            
            w = 11.4;
            h = 6;
            set(figHandle, 'PaperUnits', 'centimeters');
            set(figHandle, 'PaperSize', [w h]); %6.5 in
            set(figHandle, 'PaperPositionMode', 'manual');
            set(figHandle, 'PaperPosition', [0 0 get(figHandle, 'PaperSize')]);
            figPos = get(figHandle, 'Position');
            set(figHandle, 'Position', [figPos(1:2) figPos(3) figPos(3) * h / w]);
            
            axesHandle = subplot('Position', [0.075 0.13 0.925 0.87], 'Parent', figHandle);
            hold(axesHandle, 'on');
            
            barHandles = bar(axesHandle, edges, [freqKOGrowthEssGrowth freqKOSlowGrowth freqWtGrowth] * 100);
            set(barHandles(1), 'FaceColor', [0 1 0])
            set(barHandles(2), 'FaceColor', [1 0 0])
            set(barHandles(3), 'FaceColor', [0 0 1])
            set(barHandles, 'EdgeColor', 'none', 'BarWidth', 2);
            xlim(axesHandle, [-0.05 1.2]);
            ylim(axesHandle, [0 0.62] * 100);
            set(axesHandle, 'XTick', 0:0.25:1.5, 'XTickLabel', {'0' '' '1' '' '2' '' '3'});
            set(axesHandle, 'YTick', (0:0.2:0.6) * 100);
            set(axesHandle, 'FontSize', 5, 'TickDir', 'out');
            
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.05, 'ytickoffset', 0.028);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, 'FontSize', 5);
            set(yTicks, 'FontSize', 5);
            delete(xTicks(2:2:end));
            
            text(0.87, 0.60 * 100, 'Wild type', ...
                'Color', [0 0 1], 'HorizontalAlign', 'center', 'FontSize', 7, 'Parent', axesHandle)
            text(0.48, 0.60 * 100, 'Quasi-essential', ...
                'Color', [1 0 0], 'HorizontalAlign', 'center', 'FontSize', 7, 'Parent', axesHandle)
            text(0.09, 0.60 * 100, 'Essential', ...
                'Color', [0 1 0], 'HorizontalAlign', 'center', 'FontSize', 7, 'Parent', axesHandle)
            
            xlabel(axesHandle, 'Growth rate constant (h^{-1})', 'FontSize', 7);
            ylabel(axesHandle, '% cells', 'FontSize', 7);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            ylabelPos = get(get(axesHandle, 'ylabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1) -0.05*100 xlabelPos(end)]);
            set(get(axesHandle, 'ylabel'), 'position', [-0.120 ylabelPos(2:end)]);
            
            line(xlim(axesHandle), [0 0], 'Parent', axesHandle, 'Color', 'k', 'LineWidth', 0.5)
            
            print(figHandle, outputFile, '-dpdf', '-rgb');
            close(figHandle);
        end
        
        function plotDeletionStrainClassificationHelper(xdata, ydata, xlabelStr, ylabelStr, geneClasses, classIdxs, plotIdx, legendLocation)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            hold(axesHandle, 'on');
            h = zeros(size(SingleGeneDeletions.DELETION_STRAIN_CLASSES, 1), 1);
            for i = 1:size(SingleGeneDeletions.DELETION_STRAIN_CLASSES, 1) - 2
                geneTfs = geneClasses == i;
                if ~any(geneTfs)
                    continue;
                end
                h(i) = plot(axesHandle, xdata(geneTfs), ydata(geneTfs), ...
                    'Marker', '.', ...
                    'LineStyle', 'none', ...
                    'Color', SingleGeneDeletions.COLORS(i, :));
            end
            classIdxs = intersect(classIdxs, find(h));
            uistack(h(classIdxs), 'top')
            set(h(classIdxs), 'MarkerSize', 12)
            legend(h(classIdxs), SingleGeneDeletions.DELETION_STRAIN_CLASSES(classIdxs), ...
                'Location', legendLocation);
            if ~isnan(range(xdata)) && range(xdata)
                xlim(axesHandle, [min(xdata) max(xdata)] + range(xdata) * [-0.05 0.05]);
            end
            if ~isnan(range(ydata)) && range(ydata)
                ylim(axesHandle, [min(ydata) max(ydata)] + range(ydata) * [-0.05 0.05]);
            end
            xlabel(axesHandle, xlabelStr, 'FontSize', 14);
            ylabel(axesHandle, ylabelStr, 'FontSize', 14);
            saveas(figHandle, [baseDir 'classification-' num2str(plotIdx) '.pdf']);
            close(figHandle);
        end
        
        function plotOverview(outFileName_Figure, outFileName_Table, ...
                sim, geneClasses, gridData, gridMetaData)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            if nargin < 1 || isempty(outFileName_Figure)
                outFileName_Figure = [baseDir filesep 'overview.pdf'];
            end
            if nargin < 2 || isempty(outFileName_Table)
                outFileName_Table = [baseDir 'overview.xls'];
            end
            
            %% initialize
            figW = 17.4;
            figH = 7.3;
            
            [~, figHandle] = PlotUtil.newAxesHandle();
            clf(figHandle);
            set(figHandle, 'PaperUnits', 'centimeters');
            set(figHandle, 'PaperSize', [figW figH]); %cell two column
            set(figHandle, 'PaperPositionMode', 'manual');
            set(figHandle, 'PaperPosition', [0 0 get(figHandle, 'PaperSize')]);
            
            %% get data
            if nargin < 3
                sim = CachedSimulationObjectUtil.load();
            end
            g = sim.gene;
            isGeneImplemented = SingleGeneDeletions.getGeneImplementations(sim);
            expEss = ismember(g.essential, {'Y', ''});
            if nargin < 4
                load([baseDir 'geneClasses.mat']);
            end
            modelEss = geneClasses ~= SingleGeneDeletions.NON_ESSENTIAL;
            
            showGeneSymbols = true;
            geneSymbols = struct;
            geneSymbols.MG006 = 'tmk';
            geneSymbols.MG022 = 'rpoE';
            geneSymbols.MG113 = 'asnS';
            geneSymbols.MG048 = 'ffh';
            geneSymbols.MG001 = 'dnaN';
            geneSymbols.MG204 = 'parC';
            geneSymbols.MG084 = 'tilS';
            
            %% sizes            
            fontSizeSubfigLabel = 8;
            fontSizeLarge = 7;
            fontSizeMed = 6;
            fontSizeSmall = 5;
            
            leftColX = 0.05;
            leftColLabelX = 0.051;
            leftColW = 0.15;
            leftColH1 = .23;
            leftColH2 = 0.23;
            leftColH3 = leftColH2;
            leftColY2 = 0.375;
            leftColY3 = 0.0622;
             
            %% part A -- model vs. observed gene essentiality
            PlotUtil.labelSubFigure('A', [leftColLabelX 1.01 -0.0286 -0.0690], figHandle);
            SingleGeneDeletions.plotGeneDeletionPrediction(sim, geneClasses, [leftColX 0.71 leftColW leftColH1]);
            
            %% part B -- deletion strain classes spark line matrix
            catLabels = {
                'WT'             SingleGeneDeletions.WILD_TYPE
                'Energy'         SingleGeneDeletions.DECOMPOSING
                'Metabolic'      SingleGeneDeletions.NON_GROWING
                'RNA'            SingleGeneDeletions.DECAYING_GROWTH_NON_RNA_PROTEIN
                'Protein'        SingleGeneDeletions.DECAYING_GROWTH_NON_PROTEIN
                'Other'          SingleGeneDeletions.NON_PERPETUATING
                'DNA'            SingleGeneDeletions.NON_REPLICATIVE
                'Cytokinesis'    SingleGeneDeletions.NON_FISSIVE
                'Term Org'       SingleGeneDeletions.NO_TERMINAL_ORGANELLE
                'Damaged'        SingleGeneDeletions.TOXIN_ACCUMULATION
                'Quasi-Ess'      SingleGeneDeletions.SLOW_GROWING
                'Non-Ess'        SingleGeneDeletions.NON_ESSENTIAL
                };
            propLabels = {
                {'NTP'}                 'ntps'
                {'Growth (fg h^{-1})'}  'growth'
                {'Protein (fg)'}        'proteinWt'
                {'RNA (fg)'}            'rnaWt'
                {'DNA (fg)'}            'dnaWt'
                {'Septum (nm)'}         'pinchedDiameter'
                {'Term Org (fg)'}       'terminalOrganelleWt'
                {'Damaged' 'Prot'}      'damagedProteins'
                };
            highlightPlot = false(size(propLabels, 1), size(catLabels, 1));
            highlightPlot(1:8, 2) = true;%energy
            highlightPlot(2:7, 3) = true;%metabolic
            highlightPlot([1:4 6:7], 4) = true;%RNA
            highlightPlot([2 3 6:7], 5) = true;%Protein
            highlightPlot(7, 6) = true;%other synthetic
            highlightPlot(5:6, 7) = true;%DNA
            highlightPlot(6, 8) = true;%cytokinesis
            highlightPlot(8, 10) = true;%damaged
            highlightPlot([2 3 7], 11) = true;%Slow
            
            groupCats = (1:size(catLabels, 1))';
            groupCats([2 9 10 12]) = [3 6 11 0];
            showProps = (1:size(propLabels, 1))';
            showProps([1 7 8]) = false;
            
            catIdxs = find(groupCats == (1:size(catLabels, 1))');
            propIdxs = find(showProps);
            
            nCats = numel(catIdxs);
            nProps = numel(propIdxs);
            
            time = (0:50000)' / 3600;
            
            if nargin >= 8
            elseif exist([baseDir 'overview.mat'], 'file')
                load([baseDir 'overview.mat']);
            else
                [gridData, gridMetaData] = SingleGeneDeletions.calcOverviewData();
                save([baseDir 'overview.mat'], 'gridData', 'gridMetaData');
            end
            
            gridData(:, :, 1) = NaN;
            
            %layout
            x = leftColX + leftColW + 0.13;
            y = 0.845;
            W = 1 - x - 0.005;
            H = 0.773;
            w = W / nCats;
            h = H / nProps;
            y1 = y + 0.12;
            y2 = y + 0.08;
            y3 = y + 0.06;
            yMargin = 0.13 / 2;
            xMargin = (h * figH) / (w * figW) * yMargin;
            xSpan = 1 - 2 * xMargin;
            ySpan = 1 - 2 * yMargin;            
            highlightColor = [255 230 230] / 255; %light red
            nonhighlightColor = [1 1 1]; %white
            
            %label
            PlotUtil.labelSubFigure('B', [x-0.0650 1.01 -0.0286 -0.0690], figHandle, fontSizeSubfigLabel);
            
            %data
            xlims = [0 time(end)] + [-1 1] * (1/xSpan - 1)/2 * time(end);
            for i = 1:nProps
                for j = 1:nCats
                    axesHandle = subplot('Position', [x+(j-1)*w  y-(i)*h  w  h]);
                    
                    switch propLabels{propIdxs(i), 2}
                        case 'ntps'
                            yticks = [0 1e6];
                            ytickLabels = {num2str(yticks(1))  '10^6'};
                        case 'growth'
                            yticks = [0 2.5];
                            ytickLabels = {num2str(yticks(1))  num2str(yticks(2))};
                        case 'rnaWt'
                            yticks = [0 0.4] * 1e-15;
                            ytickLabels = {num2str(yticks(1)*1e15)  num2str(yticks(2)*1e15)};
                        case 'proteinWt'
                            yticks = [2 7] * 1e-15;
                            ytickLabels = {num2str(yticks(1)*1e15)  num2str(yticks(2)*1e15)};
                        case 'dnaWt'
                            yticks = [0.6 1.2] * 1e-15;
                            ytickLabels = {num2str(yticks(1)*1e15)  num2str(yticks(2)*1e15)};
                        case 'pinchedDiameter'
                            yticks = [0 250e-9];
                            ytickLabels = {num2str(yticks(1))  num2str(yticks(2) * 1e9)};
                        case 'terminalOrganelleWt'
                            yticks = [0 0.1] * 1e-15;
                            ytickLabels = {num2str(yticks(1) * 1e15)  num2str(yticks(2) * 1e15)};
                        case 'damagedProteins'
                            yticks = [0 6000];
                            ytickLabels = {num2str(yticks(1))  num2str(yticks(2))};
                        otherwise
                            throw(MException('SingleGeneDeletions:error', 'undefined property %s', propLabels{propIdxs(i), 2}))
                    end
                    ylims0 = [min(yticks(1), min(min(gridData(propIdxs(i), :, :))))  max(yticks(2), max(max(gridData(propIdxs(i), :, :))))];
                    ylims = ylims0 + [-1 1] * (1/ySpan - 1)/2 * (ylims0(2) - ylims0(1));
                    
                    if highlightPlot(propIdxs(i), catIdxs(j))
                        faceColor = highlightColor;
                    else
                        faceColor = nonhighlightColor;
                    end
                    if ~all(faceColor == 1)
                        hold(axesHandle, 'on');
                        rectangle('Parent', axesHandle, 'Position', [xlims(1) ylims(1) range(xlims) range(ylims)], 'EdgeColor', 'none', 'FaceColor', faceColor);
                    end
                    
                    plot(axesHandle, time, squeeze(gridData(propIdxs(i), catIdxs(j), :)), 'Color', 'k');
                    
                    set(axesHandle, 'XTick', [], 'YTick', []);
                    set(axesHandle, 'Color', 'none', 'Box', 'off', 'Visible', 'off');
                    xlim(axesHandle, xlims);
                    ylim(axesHandle, ylims);
                    
                    if j == 1
                        ylims = ylim(axesHandle);
                        for k = 1:numel(yticks)
                            annotation(figHandle, 'textbox', [x-0.006 - xMargin*w  (y-i*h + h*(yticks(k)-ylims(1)) / (ylims(end)-ylims(1)))-0.005 0  0], ...
                                'String', ytickLabels{k}, ...
                                'FontSize', fontSizeSmall, ...
                                'EdgeColor', 'none', ...
                                'Margin', 0, ...
                                'HorizontalAlignment', 'Right', ...
                                'VerticalAlignment', 'Middle');
                            annotation(figHandle, 'line', [x  x-0.004] - xMargin*w, (y-i*h + h*(yticks(k)-ylims(1)) / (ylims(end)-ylims(1))) * [1 1] , 'Color', 'k');
                        end
                        annotation(figHandle, 'line', [x; x] - xMargin*w, [
                            y - i*h + yMargin*h + ySpan*h
                            y - i*h + yMargin*h
                            ] , 'Color', 'k');
                    end
                end
            end
            
            %time axis
            axesHandle = subplot('Position', [x+xMargin*w  y-H-yMargin*h  w*xSpan  1e-6]);
            xlim(axesHandle, [0 time(end)]);
            set(axesHandle, 'XTick', [0 10], 'YTick', [], 'TickDir', 'out', 'TickLen', [0.04 0.04], 'FontSize', fontSizeSmall, 'YColor', get(figHandle, 'Color'));
            xlabel(axesHandle, 'Time (h)', 'FontSize', fontSizeLarge);
            xlabelPos = get(get(axesHandle, 'xlabel'), 'position');
            set(get(axesHandle, 'xlabel'), 'position', [xlabelPos(1) xlabelPos(2)-50 xlabelPos(3)]);
            tick2text(axesHandle, 'axis', 'x', 'xtickoffset', 7);
            xTicks = getappdata(axesHandle, 'XTickText');
            set(xTicks, 'FontSize', fontSizeSmall, 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
            
            %grid
            for i = 0:nProps
                if i > 0
                    txtHandle = annotation(figHandle, 'textbox', [x-0.1-0.008-0.025+0.02  y - (i - 1 / 2) * h .2 0], ...
                        'String', strjoin(' ', propLabels{propIdxs(i), 1}{:}), ...
                        'HorizontalAlignment', 'right', ...
                        'VerticalAlignment', 'Middle', ...
                        'FontSize', fontSizeLarge, ...
                        'EdgeColor', 'none', 'Margin', 0, ...
                        'Interpreter', 'tex');
                    set(txtHandle, 'Position', [x-0.2-0.008-0.018  y - (i - 1 / 2) * h 0.2 0])
                end
                annotation(figHandle, 'line', x + [0 W], y - i * h * [1 1], 'Color', [0.75 0.75 0.75]);
            end
            txtHandles = zeros(nCats, 1);
            for i = 0:nCats
                if i > 0
                    txtHandles(i) = annotation(figHandle, 'textbox', [x+(i-1/2)*w  y+0.03  0  0], ...
                        'String', catLabels{catIdxs(i), 1}, ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'bottom', ...
                        'FontSize', fontSizeMed, ...
                        'EdgeColor', 'none', ...
                        'Margin', 0, ...
                        'FitBoxToText', 'on');
                end
                if i > 0 && catLabels{catIdxs(i), 2} ~= SingleGeneDeletions.WILD_TYPE
                    nGenes = 0;
                    nCorrect = 0;
                    for j = 1:numel(groupCats)
                        if groupCats(j) == catIdxs(i)
                            nGenes = nGenes + ...
                                sum(isGeneImplemented & geneClasses == catLabels{j, 2});
                            nCorrect = nCorrect + ...
                                sum(expEss == modelEss & isGeneImplemented & geneClasses == catLabels{j, 2});
                        end
                    end
                    
                    tmp = strrep(gridMetaData(catIdxs(i)).gene, '_', '');
                    if showGeneSymbols
                        tmp = sprintf('{\\it{%s}}', geneSymbols.(tmp));
                    end
                    
                    txtHandles(i) = annotation(figHandle, 'textbox', [x+(i-1/2)*w  y+0.005  0  0], ...
                        'String', sprintf('(%d, %s)', nGenes, tmp), ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'bottom', ...
                        'FontSize', fontSizeSmall, ...
                        'EdgeColor', 'none', ...
                        'Margin', 0, ...
                        'FitBoxToText', 'on');
                end
                annotation(figHandle, 'line', x + i * w * [1 1], y - [0 H], 'Color', [0.75 0.75 0.75]);
            end
            
            %dendrogram
            annotation(figHandle, 'line', x + ( 1.5) * w * [1 1] + [0 6*w], [y1 y1],  'Color', 'k');
            annotation(figHandle, 'line', x + ( 1.5) * w * [1 1], [y1 y3],  'Color', 'k');
            annotation(figHandle, 'line', x + ( 3.5) * w * [1 1], [y1 y2+0.03],  'Color', 'k');
            annotation(figHandle, 'line', x + ( 6) * w * [1 1], [y1 y2+0.03],  'Color', 'k');
            annotation(figHandle, 'line', x + (7.5) * w * [1 1], [y3 y1],  'Color', 'k');
            
            annotation(figHandle, 'line', x + ( 2.5) * w * [1 1], [y3 y2],  'Color', 'k');
            annotation(figHandle, 'line', x + ( 2.5) * w * [1 1] + [0 2*w], [y2 y2],  'Color', 'k');
            annotation(figHandle, 'line', x + ( 4.5) * w * [1 1], [y2 y3],  'Color', 'k');
            
            annotation(figHandle, 'line', x + (5.5) * w * [1 1], [y3 y2],  'Color', 'k');
            annotation(figHandle, 'line', x + (5.5) * w * [1 1] + [0 w], [y2 y2],  'Color', 'k');
            annotation(figHandle, 'line', x + (6.5) * w * [1 1], [y2 y3],  'Color', 'k');
            
            txtBox = [
                annotation(figHandle, 'textbox', [x+4.5*w-0.15  y1+0.02 0.3 0], 'String', 'Essential')
                annotation(figHandle, 'textbox', [x+3.5*w-0.15  y2+0.02 0.3 0], 'String', 'Macromolecule synthesis')
                annotation(figHandle, 'textbox', [x+6*w-0.15  y2+0.02 0.3 0],   'String', 'Cell cycle')
                ];
            set(txtBox, 'FontSize', fontSizeLarge, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Margin', 0, 'EdgeColor', 'none', 'FitBoxToText', 'on');
            
            %% part C -- methionine aminopeptidase, type I (MG_172)
            
            nGenerations = 10;
            nSamples = 20;
            if exist([baseDir filesep 'MG_172-analysis2.mat'], 'file')
                load([baseDir filesep 'MG_172-analysis2.mat'])
            else
                g = sim.gene;
                r = sim.state('Rna');
                pm = sim.state('ProteinMonomer');
                pc = sim.state('ProteinComplex');
                mr = sim.state('MetabolicReaction');
                ppI = sim.process('ProteinProcessingI');
                met = sim.process('Metabolism');
                
                sim.setForTest('processesInInitOrder', sim.processesInInitOrder(~ismember(...
                    cellfun(@(p) p.name, sim.processesInInitOrder, 'UniformOutput', false), {...
                    'Transcription'
                    'FtsZPolymerization'
                    'ChromosomeCondensation'
                    'Translation'
                    'DNASupercoiling'
                    'ReplicationInitiation'
                    })));
                
                pcComp = sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3);
                monIdxs = ppI.nascentMonomerNTerminalMethionineCleavages;
                unprocessedMonIdxs = ~ppI.nascentMonomerNTerminalMethionineCleavages;
                cpxIdxs = find(any(pcComp(monIdxs, :), 1));
                
                rnaExp = r.expression;
                rnaDecayRates = r.decayRates;
                monDecayRates = pm.decayRates;
                cpxDecayRates = pc.decayRates;
                
                cleavedMonWt = zeros(nGenerations, 3);
                growths = zeros(nGenerations, 3);
                for i = 1:nGenerations
                    tmpCleavedMonWt = zeros(nSamples, 1);
                    tmpGrowths = zeros(nSamples, 1);
                    for j = 1:nSamples
                        r.expression = rnaExp;
                        r.decayRates = rnaDecayRates;
                        pm.decayRates = monDecayRates;
                        pc.decayRates = cpxDecayRates;
                        mr.initialGrowthFilterWidth = Inf;
                        
                        sim.applyOptions('seed', i * nGenerations + j);
                        sim.allocateMemoryForState(1);
                        sim.initializeState();
                        sim.applyPerturbations();
                        
                        mons = pm.counts;
                        cpxs = pc.counts;
                        
                        pm.counts(pm.processedIIIndexs(monIdxs),    :) = pm.randStream.random('binomial', pm.counts(pm.processedIIIndexs(monIdxs), :), 1 / (2.^(i - 1)));
                        pm.counts(pm.signalSequenceIndexs(monIdxs), :) = pm.randStream.random('binomial', pm.counts(pm.signalSequenceIndexs(monIdxs), :), 1 / (2.^(i - 1)));
                        pm.counts(pm.foldedIndexs(monIdxs),         :) = pm.randStream.random('binomial', pm.counts(pm.foldedIndexs(monIdxs), :), 1 / (2.^(i - 1)));
                        pm.counts(pm.matureIndexs(monIdxs),         :) = pm.randStream.random('binomial', pm.counts(pm.matureIndexs(monIdxs), :), 1 / (2.^(i - 1)));
                        pm.counts(pm.inactivatedIndexs(monIdxs),    :) = pm.randStream.random('binomial', pm.counts(pm.inactivatedIndexs(monIdxs), :), 1 / (2.^(i - 1)));
                        pm.counts(pm.boundIndexs(monIdxs),          :) = pm.randStream.random('binomial', pm.counts(pm.boundIndexs(monIdxs), :), 1 / (2.^(i - 1)));
                        pm.counts(pm.misfoldedIndexs(monIdxs),      :) = pm.randStream.random('binomial', pm.counts(pm.misfoldedIndexs(monIdxs), :), 1 / (2.^(i - 1)));
                        pm.counts(pm.damagedIndexs(monIdxs),        :) = pm.randStream.random('binomial', pm.counts(pm.damagedIndexs(monIdxs), :), 1 / (2.^(i - 1)));
                        
                        pc.counts(pc.nascentIndexs(cpxIdxs),     :) = pc.randStream.random('binomial', pc.counts(pc.nascentIndexs(cpxIdxs), :), 1 / (2.^(i - 1)));
                        pc.counts(pc.matureIndexs(cpxIdxs),      :) = pc.randStream.random('binomial', pc.counts(pc.matureIndexs(cpxIdxs), :), 1 / (2.^(i - 1)));
                        pc.counts(pc.inactivatedIndexs(cpxIdxs), :) = pc.randStream.random('binomial', pc.counts(pc.inactivatedIndexs(cpxIdxs), :), 1 / (2.^(i - 1)));
                        pc.counts(pc.boundIndexs(cpxIdxs),       :) = pc.randStream.random('binomial', pc.counts(pc.boundIndexs(cpxIdxs), :), 1 / (2.^(i - 1)));
                        pc.counts(pc.misfoldedIndexs(cpxIdxs),   :) = pc.randStream.random('binomial', pc.counts(pc.misfoldedIndexs(cpxIdxs), :), 1 / (2.^(i - 1)));
                        pc.counts(pc.damagedIndexs(cpxIdxs),     :) = pc.randStream.random('binomial', pc.counts(pc.damagedIndexs(cpxIdxs), :), 1 / (2.^(i - 1)));
                        
                        diffMons = ...
                            + mons(pm.processedIIIndexs,    :) - pm.counts(pm.processedIIIndexs, :) ...
                            + mons(pm.signalSequenceIndexs, :) - pm.counts(pm.signalSequenceIndexs, :) ...
                            + mons(pm.foldedIndexs,         :) - pm.counts(pm.foldedIndexs, :) ...
                            + mons(pm.matureIndexs,         :) - pm.counts(pm.matureIndexs, :) ...
                            + mons(pm.inactivatedIndexs,    :) - pm.counts(pm.inactivatedIndexs, :) ...
                            + mons(pm.boundIndexs,          :) - pm.counts(pm.boundIndexs, :) ...
                            + mons(pm.misfoldedIndexs,      :) - pm.counts(pm.misfoldedIndexs, :) ...
                            + mons(pm.damagedIndexs,        :) - pm.counts(pm.damagedIndexs, :);
                        diffCpxs = ...
                            + cpxs(pc.nascentIndexs,     :) - pc.counts(pc.nascentIndexs, :) ...
                            + cpxs(pc.matureIndexs,      :) - pc.counts(pc.matureIndexs, :) ...
                            + cpxs(pc.inactivatedIndexs, :) - pc.counts(pc.inactivatedIndexs, :) ...
                            + cpxs(pc.boundIndexs,       :) - pc.counts(pc.boundIndexs, :) ...
                            + cpxs(pc.misfoldedIndexs,   :) - pc.counts(pc.misfoldedIndexs, :) ...
                            + cpxs(pc.damagedIndexs,     :) - pc.counts(pc.damagedIndexs, :);
                        pm.counts(pm.processedIIndexs, :) = ...
                            + pm.counts(pm.processedIIndexs, :) ...
                            + diffMons;
                        pm.counts(pm.processedIIndexs(monIdxs), :) = ...
                            + pm.counts(pm.processedIIndexs(monIdxs), :) ...
                            + pcComp(monIdxs, :) * diffCpxs;
                        pm.counts(pm.matureIndexs(unprocessedMonIdxs), :) = ...
                            + pm.counts(pm.matureIndexs(unprocessedMonIdxs), :) ...
                            + pcComp(unprocessedMonIdxs, :) * diffCpxs;
                        
                        pmCnts = sum(pm.counts, 2);
                        pcCnts = sum(pc.counts, 2);
                        totMons = ...
                            + pmCnts(pm.matureIndexs) ...
                            + pmCnts(pm.boundIndexs) ...
                            + sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3) * (...
                                + pcCnts(pc.matureIndexs) ...
                                + pcCnts(pc.boundIndexs));
                        tmpCleavedMonWt(j) = (totMons .* ppI.nascentMonomerNTerminalMethionineCleavages)' * pm.molecularWeights(pm.matureIndexs) / ConstantUtil.nAvogadro * 1e15;
                        
                        processes = sim.processes;
                        sim.setForTest('processes', {met});
                        sim.evolveState();
                        sim.setForTest('processes', processes);
                        
                        tmpGrowths(j) = mr.growth;
                    end
                    
                    cleavedMonWt(i, 1) = mean(tmpCleavedMonWt);
                    cleavedMonWt(i, 2:3) = quantile(tmpCleavedMonWt, [0.341 0.682]) - mean(tmpCleavedMonWt);
                    
                    growths(i, 1) = mean(tmpGrowths);
                    growths(i, 2:3) = quantile(tmpGrowths, [0.341 0.682]) - mean(tmpGrowths);
                end
                
                save([baseDir filesep 'MG_172-analysis2.mat'], 'growths', 'cleavedMonWt');
            end
            
            growths = growths * 3600 * sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight) * 1e15;
                        
            PlotUtil.labelSubFigure('C', [leftColLabelX leftColY2+leftColH2+0.05 -0.0286 -0.0690], figHandle, fontSizeSubfigLabel);
            
            axesHandle = subplot('Position', [leftColX leftColY2 leftColW leftColH2]);
            hold(axesHandle, 'on');
            errorbar((0:nGenerations-1)', cleavedMonWt(ones(nGenerations, 1), 1), cleavedMonWt(ones(nGenerations, 1), 2), cleavedMonWt(ones(nGenerations, 1), 3), 'Color', 'k');
            errorbar((0:nGenerations-1)', cleavedMonWt(:, 1), cleavedMonWt(:, 2), cleavedMonWt(:, 3), 'Color', [0 0 255] / 255);
            xlabel(axesHandle, 'Generation', 'FontSize', fontSizeLarge);
            ylabel(axesHandle, {'N-term cleaved' 'protein (fg)'}, 'FontSize', fontSizeLarge);
            xlim(axesHandle, [-0.25 8.25])
            ylim(axesHandle, [-0.01 max(cleavedMonWt(:, 1) + cleavedMonWt(:, 3))])
            set(axesHandle, 'XTick', [0 5])
            set(axesHandle, 'YTick', 0:0.1:0.3)
            set(axesHandle, 'FontSize', fontSizeSmall)
            text(8.2, 0.315, 'WT', ...
                'Margin', 1e-6, ...
                'FontSize', fontSizeLarge, ...
                'Color', 'k', ...
                'Parent', axesHandle, ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'Top')
            text(8.2, 0, '{\Delta}map', ...
                'Margin', 1e-6, ...
                'FontSize', fontSizeLarge, ...
                'Color', [0 0 255] / 255, ...
                'Parent', axesHandle, ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'Bottom');
            set(axesHandle, 'FontSize', fontSizeSmall, 'TickLen', [0.01 0.01], 'TickDir', 'out');
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.07, 'ytickoffset', 0.02);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, ...
                'FontSize', fontSizeSmall, ...
                'HorizontalAlignment', 'Center', ...
                'VerticalAlignment', 'Middle');
            set(yTicks, ...
                'FontSize', fontSizeSmall, ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'Middle');
            delete(yTicks(2:end-1));
            
            xLabelPos = get(get(axesHandle, 'xlabel'), 'Position');
            set(get(axesHandle, 'xlabel'), 'Position', [xLabelPos(1) -0.055 xLabelPos(3)]);
            
            yLabelPos = get(get(axesHandle, 'ylabel'), 'Position');
            set(get(axesHandle, 'ylabel'), 'Position', [-1.3 yLabelPos(2:end)]);
            
            %% part D -- methionine aminopeptidase, type I (MG_172)
            PlotUtil.labelSubFigure('D', [leftColLabelX leftColY3+leftColH3+0.05 -0.0286 -0.0690], figHandle, fontSizeSubfigLabel);
            
            axesHandle = subplot('Position', [leftColX leftColY3 leftColW leftColH3]);
            hold(axesHandle, 'on');
            errorbar((0:nGenerations-1)', growths(ones(nGenerations, 1), 1), growths(ones(nGenerations, 1), 2), growths(ones(nGenerations, 1), 3), 'Color', 'k');
            errorbar((0:nGenerations-1)', growths(:, 1), growths(:, 2), growths(:, 3), 'Color', [0 0 255] / 255);
            xlabel(axesHandle, 'Generation', 'FontSize', fontSizeLarge);
            ylabel(axesHandle, 'Growth (fg h^{-1})', 'FontSize', fontSizeLarge);
            xlim(axesHandle, [-0.25 8.25])
            ylim(axesHandle, [-0.05 max(growths(:, 1) + growths(:, 3))])
            set(axesHandle, 'XTick', [0 5])
            set(axesHandle, 'YTick', 0:0.25:1)
            set(axesHandle, 'FontSize', fontSizeSmall)
            text(8.2, 0.92, 'WT', ...
                'Margin', 1e-6, ...
                'FontSize', fontSizeLarge, ...
                'Color', 'k', ...
                'Parent', axesHandle, ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'Top')
            text(8.2, 0, '{\Delta}map', ...
                'Margin', 1e-6, ...
                'FontSize', fontSizeLarge, ...
                'Color', [0 0 255] / 255, ...
                'Parent', axesHandle, ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'Bottom');
            set(axesHandle, 'FontSize', fontSizeSmall, 'TickLen', [0.01 0.01], 'TickDir', 'out');
            tick2text(axesHandle, 'axis', 'xy', 'xtickoffset', 0.07, 'ytickoffset', 0.02);
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, ...
                'FontSize', fontSizeSmall, ...
                'HorizontalAlignment', 'Center', ...
                'VerticalAlignment', 'Middle');
            set(yTicks, ...
                'FontSize', fontSizeSmall, ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'Middle');
            delete(yTicks(2:end-1));
            
            xLabelPos = get(get(axesHandle, 'xlabel'), 'Position');
            set(get(axesHandle, 'xlabel'), 'Position', [xLabelPos(1) -0.18 xLabelPos(3)]);
            
            yLabelPos = get(get(axesHandle, 'ylabel'), 'Position');
            set(get(axesHandle, 'ylabel'), 'Position', [-1.3 yLabelPos(2:end)]);
            
            %% standarize tick lengths
            paperSize = get(figHandle, 'PaperSize');
            figW = paperSize(1);
            figH = paperSize(2);
            axesHandles = findobj(figHandle, 'type', 'axes', '-and', 'tag', '');
            axesPos = cell2mat(get(axesHandles, 'position'));
            tickLen = 0.0075 * max(max(axesPos(:, 3:4), [], 1) .* [figW figH]);
            for i = 1:numel(axesHandles)
                maxVal = max(axesPos(i, 3:4) .* [figW figH]);
                set(axesHandles(i), 'TickLen', [tickLen / maxVal  0.0250]);
            end
            
            %% save figure
            print(figHandle, outFileName_Figure, '-dpdf', '-rgb');
            close(figHandle);
            
            %% save
            colLabels = {'Class' 'Name' 'No Genes' 'No Simulations' 'Representative Gene Locus Tag' 'Representative Gene Symbol' 'Representative Simulation Batch' 'Representative Simulation Index'};
            content = cell(0, numel(colLabels));
            for i = 1:numel(catIdxs)
                nGenes = 0;
                nSimulations = 0;
                for j = 1:numel(groupCats)
                    if groupCats(j) == catIdxs(i)
                        nGenes = nGenes + gridMetaData(j).nGenes;
                        nSimulations = nSimulations + gridMetaData(j).nSimulations;
                    end
                end
                
                locusTag = gridMetaData(catIdxs(i)).gene;
                if ~isempty(locusTag)
                    symbol = geneSymbols.(strrep(locusTag, '_', ''));
                else
                    symbol = [];
                end
                
                content = [content; {
                    catLabels{catIdxs(i), 1} ...
                    SingleGeneDeletions.DELETION_STRAIN_CLASSES{catLabels{catIdxs(i), 2}, 1} ...
                    nGenes ...
                    nSimulations ...
                    locusTag ...                    
                    symbol, ...                   
                    gridMetaData(catIdxs(i)).simGroup ...
                    gridMetaData(catIdxs(i)).simIdx ...
                    }];
            end

            if ispc && exist(outFileName_Table, 'file')
                delete(outFileName_Table)
            end
            PrintUtil.printToFile(content, colLabels, outFileName_Table, 'Single Gene Deletions');
        end
        
        function plotGeneDeletionPrediction(sim, geneClasses, position)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            g = sim.gene;
            isGeneImplemented = SingleGeneDeletions.getGeneImplementations(sim);
            expEss = ismember(g.essential, {'Y', ''});
            modelEss = geneClasses ~= SingleGeneDeletions.NON_ESSENTIAL;
            
            %% fontsizes
            fontSizeLarge = 7;
            fontSizeSmall = 5;
                        
            axesHandle = subplot('Position', position);
            axis(axesHandle, 'equal');
            xlim(axesHandle, [-0.5 1.5]);
            ylim(axesHandle, [-0.5 1.5]);
            set(axesHandle, 'Box', 'on');
            line([0.5 0.5], [-0.5 1.5], 'Parent', axesHandle, 'Color', 'k')
            line([-0.5 1.5], [0.5 0.5], 'Parent', axesHandle, 'Color', 'k')
            set(axesHandle, 'XAxisLocation', 'top');
            set(axesHandle, 'XDir','reverse');
            
            set(axesHandle, 'TickLength', [0 0]);
            set(axesHandle, 'XTick', [0 1]);
            set(axesHandle, 'YTick', [0 1]);
            tick2text(axesHandle, 'axis', 'xy', 'ytickoffset', 0.06)
            xTicks = getappdata(axesHandle, 'XTickText');
            yTicks = getappdata(axesHandle, 'YTickText');
            set(xTicks, ...
                'FontSize', fontSizeSmall, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
            set(xTicks(1), 'String', 'Non-ess', 'Position', [0 1.68 -1]);
            set(xTicks(2), 'String', 'Essential', 'Position', [1 1.68 -1]);
            set(yTicks, 'rotation', 90, ...
                'FontSize', fontSizeSmall, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
            set(yTicks(1), 'String', 'Non-ess', 'Position', [1.68 0 -1]);
            set(yTicks(2), 'String', 'Essential', 'Position', [1.68 1 -1]);
            
            xlabel(axesHandle', 'Model');
            ylabel(axesHandle, 'Experiment');
            xlab = get(axesHandle, 'xLabel');
            ylab = get(axesHandle, 'yLabel');
            set(xlab, ...
                'FontSize', fontSizeLarge, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Position', [0.5 1.9 1]);
            set(ylab, ...
                'FontSize', fontSizeLarge, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Position', [1.9 0.5 1]);
            
            data = zeros(2, 2);
            data(1, 1) = sum(~modelEss & ~expEss & isGeneImplemented);
            data(1, 2) = sum(~modelEss &  expEss & isGeneImplemented);
            data(2, 1) = sum( modelEss & ~expEss & isGeneImplemented);
            data(2, 2) = sum( modelEss &  expEss & isGeneImplemented);
            for i = 0:1
                for j = 0:1
                    text(i, j, num2str(data(i+1, j+1)), ...
                        'Parent', axesHandle, ...
                        'FontSize', fontSizeLarge, ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle');
                end
            end
            
            text(1.75, -0.68, sprintf('Correct:'), ...
                'Parent', axesHandle, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', fontSizeLarge);
            text(1.75, -0.98, sprintf('Incorrect:'), ...
                'Parent', axesHandle, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', fontSizeLarge);
            text(-0.75, -0.68, sprintf('%d (%d%%)', sum(diag(data)), round(100 * sum(diag(data)) / sum(data(:)))), ...
                'Parent', axesHandle, ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', fontSizeLarge);
            text(-0.75, -0.98, sprintf('%d (%d%%)', sum(data(:)) - sum(diag(data)), round(100 * (1 - sum(diag(data)) / sum(data(:))))), ...
                'Parent', axesHandle, ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', fontSizeLarge);
        end
        
        function [gridData, gridMetaData] = calcOverviewData(iJob, nJobs)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 2
                iJob = 1;
                nJobs = 1;
            end
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            load([baseDir 'geneClasses.mat']);
            
            sim = CachedSimulationObjectUtil.load();
            g = sim.gene;
            
            isGeneImplemented = SingleGeneDeletions.getGeneImplementations(sim);
            
            nCats = size(SingleGeneDeletions.DELETION_STRAIN_CLASSES, 1) - 1;
            nProps = 8;
            nTime = 50001;
            catLabels = {
                'WT'             SingleGeneDeletions.WILD_TYPE                        'growth'
                'Energy'         SingleGeneDeletions.DECOMPOSING                      'ntps'
                'Metabolic'      SingleGeneDeletions.NON_GROWING                      'growth'
                'RNA'            SingleGeneDeletions.DECAYING_GROWTH_NON_RNA_PROTEIN  'rnaWt'
                'Protein'        SingleGeneDeletions.DECAYING_GROWTH_NON_PROTEIN      'proteinWt'
                'Other'          SingleGeneDeletions.NON_PERPETUATING                 'proteinWt'
                'DNA'            SingleGeneDeletions.NON_REPLICATIVE                  'dnaWt'
                'Cytokinesis'    SingleGeneDeletions.NON_FISSIVE                      'growth'
                'Term Org'       SingleGeneDeletions.NO_TERMINAL_ORGANELLE            'terminalOrganelleWt'
                'Damaged'        SingleGeneDeletions.TOXIN_ACCUMULATION               'damagedProteins'
                'Slow'           SingleGeneDeletions.SLOW_GROWING                     'rnaWt'
                'Non-Ess'        SingleGeneDeletions.NON_ESSENTIAL                    'growth'
                };
            propLabels = {
                {'NTP'}                 'ntps'
                {'Growth (fg h^{-1})'}  'growth'
                {'Protein (fg)'}        'proteinWt'
                {'RNA (fg)'}            'rnaWt'
                {'DNA (fg)'}            'dnaWt'
                {'Septum (nm)'}         'pinchedDiameter'
                {'Term Org (fg)'}       'terminalOrganelleWt'
                {'Damaged' 'Prot'}      'damagedProteins'
                };
            
            gridData = NaN(nProps, nCats, nTime);
            gridMetaData = repmat(struct('nSimulations', [], 'nGenes', [], 'gene', [], 'simGroup', [], 'simIdx', []), nCats, 1);
            for i = iJob:nJobs:nCats
                if exist([baseDir 'overview-' catLabels{i, 1} '.mat'], 'file')
                    load([baseDir 'overview-' catLabels{i, 1} '.mat']);
                else
                    %get data for all simulations in category
                    if catLabels{i, 2} == SingleGeneDeletions.WILD_TYPE
                        tmp2 = load([baseDir 'WT.mat'], 'metaData', 'simEndTimes', 'ntps', 'growth', 'rnaWt', 'proteinWt', 'dnaWt', 'pinchedDiameter', 'terminalOrganelleWt', 'damagedProteins');
                        tmp = tmp2.growth(sub2ind(size(tmp2.growth), tmp2.simEndTimes, (1:numel(tmp2.simEndTimes))'));
                        for k = 1:nProps
                            for l = 1:numel(tmp2.simEndTimes)
                                if isnan(tmp(l)) || ~isfield(tmp2, propLabels{k, 2}) || any(isnan(tmp2.(propLabels{k, 2})(2:tmp2.simEndTimes(l), l)))
                                    tmp(l) = NaN;
                                end
                            end
                        end
                        [~, idx] = min(abs(tmp - nanmedian(tmp)));
                        
                        tmpMetaData = struct(...
                            'nSimulations', numel(tmp), ...
                            'nGenes', [], ...
                            'gene', [], ...
                            'simGroup', tmp2.metaData(idx).simGroup, ...
                            'simIdx', tmp2.metaData(idx).simIdx);
                    else
                        geneIdxs = find(geneClasses == catLabels{i, 2} & isGeneImplemented);
                        
                        tmp = zeros(0, 1);
                        genes = zeros(0, 1);
                        simIdxs = zeros(0, 1);
                        for j = 1:numel(geneIdxs)
                            tmp2 = load([baseDir g.wholeCellModelIDs{geneIdxs(j)} '.mat'], 'simEndTimes', catLabels{i, 3});
                            if ~isfield(tmp2, catLabels{i, 3})
                                continue;
                            end
                            tmp = [tmp; tmp2.(catLabels{i, 3})(sub2ind(size(tmp2.(catLabels{i, 3})), tmp2.simEndTimes, (1:numel(tmp2.simEndTimes))'))];
                            genes = [genes; repmat(geneIdxs(j), numel(tmp2.simEndTimes), 1)];
                            simIdxs = [simIdxs; (1:numel(tmp2.simEndTimes))'];
                        end
                        
                        while true
                            switch catLabels{i, 2}
                                case SingleGeneDeletions.NO_TERMINAL_ORGANELLE
                                    [~, idx] = min(tmp);                                
                                otherwise
                                    [~, idx] = min(abs(tmp - nanmedian(tmp)));
                            end
                            tmp2 = load([baseDir g.wholeCellModelIDs{genes(idx)} '.mat'], 'metaData', 'simEndTimes', propLabels{:, 2});
                            
                            if ~all(ismember(propLabels(:, 2), fieldnames(tmp2)))
                                tmp(idx) = NaN;
                            end
                            for k = 1:nProps
                                if ~isfield(tmp2, propLabels{k, 2}) || any(isnan(tmp2.(propLabels{k, 2})(2:tmp2.simEndTimes(simIdxs(idx)), simIdxs(idx))))
                                    tmp(idx) = NaN;
                                end
                            end
                            
                            if catLabels{i, 2} == SingleGeneDeletions.NON_FISSIVE && (...
                                    isnan(range(tmp2.dnaWt(:, simIdxs(idx)))) || ~range(tmp2.dnaWt(:, simIdxs(idx))) || ...
                                    isnan(tmp2.rnaWt(end-4*3600, simIdxs(idx))) || tmp2.rnaWt(end, simIdxs(idx)) < tmp2.rnaWt(end-4*3600, simIdxs(idx)))
                                tmp(idx) = NaN;
                            end
                            
                            if ~isnan(tmp(idx))
                                break;
                            end
                        end
                        
                        tmpMetaData = struct(...
                            'nSimulations', numel(tmp), ...
                            'nGenes', numel(geneIdxs), ...
                            'gene', g.wholeCellModelIDs{genes(idx)}, ...
                            'simGroup', tmp2.metaData(simIdxs(idx)).simGroup, ...
                            'simIdx', tmp2.metaData(simIdxs(idx)).simIdx);
                        
                        idx = simIdxs(idx);
                    end
                    
                    %get data for representative cell
                    tmpData = struct;
                    for k = 1:nProps
                        tmpData.(propLabels{k, 2}) = tmp2.(propLabels{k, 2})(:, idx);
                    end
                    
                    %save
                    save([baseDir 'overview-' catLabels{i, 1} '.mat'], 'tmpMetaData', 'tmpData');
                    
                    %cleanup
                    clear tmp2;
                end
                
                gridMetaData(i) = tmpMetaData;
                for j = 1:nProps
                    gridData(j, i, 1:numel(tmpData.(propLabels{j, 2}))) = tmpData.(propLabels{j, 2});
                end
            end
        end
    end
    
    %metabolic simulations of single gene deletions
    methods (Static)
        function [deletionMetabolicOutputs, wtMetabolicOutput] = calculateSingleGeneDeletionMetabolicOutputs(sim)
            %constants
            g = sim.gene;
            comp = sim.compartment;
            met = sim.process('Metabolism');
            nGenes = size(g.wholeCellModelIDs, 1);
            ntpIdxs = met.substrateIndexs({'ATP'; 'GTP'});
            energyProductionRxnIdxs = met.fbaReactionIndexs_metaboliteInternalExchange(...
                ismember(...
                met.substrateIndexs_internalExchangedMetabolites, ...
                sub2ind(size(met.substrates), ntpIdxs, repmat(comp.cytosolIndexs, size(ntpIdxs))) ...
                ) ...
                );
            
            %setup substrates so that NMPs, NDPs are available to be
            %recycled and energy can be produced
            sim.initializeState();
            
            %calculate wild-type metabolic output: growth rate and energy
            %production
            wtMetabolicOutput = zeros(1, 3);
            [wtMetabolicOutput(1, 1), ~, fbaReactionFluxs] = met.calcGrowthRate(met.calcFluxBounds(met.substrates, met.enzymes, met.fbaReactionBounds, met.fbaEnzymeBounds));
            wtMetabolicOutput(1, 2:3) = fbaReactionFluxs(energyProductionRxnIdxs);
            
            %calculate single gene deletion metabolic output: growth rate
            %and energy production
            deletionMetabolicOutputs = zeros(nGenes,3);
            deletionMetabolicOutputs(~any(met.enzymeGeneComposition, 2) & ~any(met.substrateGeneComposition, 2), :) = NaN;
            enzGeneComp = met.enzymeGeneComposition();
            subGeneComp = met.substrateGeneComposition();
            for i = 1:nGenes
                enzIdxs = find(enzGeneComp(i, :));
                subIdxs = find(subGeneComp(i, :));
                if isempty(enzIdxs) && isempty(subIdxs)
                    continue;
                end
                
                enzymes = met.enzymes;
                substrates = met.substrates;
                enzymes(enzIdxs, :) = 0;
                substrates(subIdxs, :) = 0;
                [deletionMetabolicOutputs(i, 1), ~, fbaReactionFluxs] = ...
                    met.calcGrowthRate(met.calcFluxBounds(substrates, enzymes, met.fbaReactionBounds, met.fbaEnzymeBounds));
                deletionMetabolicOutputs(i, 2:3) = fbaReactionFluxs(energyProductionRxnIdxs);
            end
        end
    end
    
    %summary of all single gene deletion strains and simulations
    methods (Static)
        function [content, colLabels] = printSingleGeneDeletionStrains(sim, deletionMetabolicOutputs, wtMetabolicOutput, deletionStats, wtStats, geneClasses)
            %import classes
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.Simulation;
            import edu.stanford.covert.util.ConstantUtil;
            
            g = sim.gene;
            rna = sim.state('Rna');
            [~, geneImplementedMolecules, geneImplementedProcesses] = SingleGeneDeletions.getGeneImplementations(sim);
            geneImplementedMolecules = cellfun(@(molecules) strjoin(', ', molecules{:}), geneImplementedMolecules, 'UniformOutput', false);
            geneImplementedProcesses = cellfun(@(processes) strjoin(', ', processes{:}), geneImplementedProcesses, 'UniformOutput', false);
            
            colLabels = {
                'Gene'
                'Name'
                'Transcription Unit'
                'Transcription Unit Name'
                'Implementation Molecule'
                'Implementation Process'
                'Model Class'
                'Model Essentiality'
                'Mean Model Growth (fg h^-1)'
                'Std Model Growth (fg h^-1)'
                'Mean Model ATP Production (molecules s^-1)'
                'Std Model ATP Production (molecules s^-1)'
                'Mean Model GTP Production (molecules s^-1)'
                'Std Model GTP Production (molecules s^-1)'
                'Mean Model Replication Initiation Duration (h)'
                'Std Model Replication Initiation Duration (h)'
                'Mean Model Replication Duration (h)'
                'Std Model Replication Duration (h)'
                'Mean Model Cytokinesis Duration (h)'
                'Std Model Cytokinesis Duration (h)'
                'Mean Model Cell Cycle Duration (h)'
                'Std Model Cell Cycle Duration (h)'
                'Mean Model Mass Doubling Duration (h)'
                'Std Model Mass Doubling Duration (h)'
                'Mean Cumulative ATP Usage (molecules)'
                'Std Cumulative ATP Usage (molecules)'
                'Mean Cumulative GTP Usage (molecules)'
                'Std Cumulative GTP Usage (molecules)'
                'Mean Cell Weight Fraction DNA (%)'
                'Std Cell Weight Fraction DNA (%)'
                'Mean Cell Weight Fraction RNA (%)'
                'Std Cell Weight Fraction RNA (%)'
                'Mean Cell Weight Fraction Protein (%)'
                'Std Cell Weight Fraction Protein (%)'
                'Mean Superhelicity'
                'Std Superhelicity'
                'Mean Exponential Rate RNA growth (h^-1)'
                'Std Exponential Rate RNA growth (h^-1)'
                'Mean Exponential Rate Protein growth (h^-1)'
                'Std Exponential Rate Protein growth (h^-1)'
                'Mean Exponential Rate Damaged RNA growth (h^-1)'
                'Std Exponential Rate Damaged RNA growth (h^-1)'
                'Mean Exponential Rate Damaged Protein growth (h^-1)'
                'Std Exponential Rate Damaged Protein growth (h^-1)'
                'Mean Antibiotics (molecules)'
                'Std Antibiotics (molecules)'
                'FBA Growth (fg h^-1)'
                'FBA ATP Production (molecules s^-1)'
                'FBA GTP Production (molecules s^-1)'
                'Mean Experimental Mass Doubling Time (h)'
                'Std Experimental Mass Doubling Time (h)'
                'Experimental No. Replicates'
                'FBA Essentiality'
                'Experiment Essentiality'
                'Model Agree'
                'FBA Agree'
                }';
            content = cell(0, numel(colLabels));
            
            %content
            deletionGrowthRates = deletionMetabolicOutputs(:, 1);
            wtGrowthRate = wtMetabolicOutput(:, 1);
            deletionEnergyProductions = deletionMetabolicOutputs(:, 2:3);
            wtEnergyProduction = wtMetabolicOutput(:, 2:3);
            
            [deletionDoublingTimesExperimental, wtDoublingTimeExperimental] = SingleGeneDeletions.getExperimentalDoublingTimes();
            wtGrowthRateExperimental = log(2) ./ wtDoublingTimeExperimental(:, 1);
            deletionGrowthRatesExperimental = wtGrowthRateExperimental * (~ismember(g.essential, {'Y', ''}));
            essThresh = 0.2 * wtGrowthRateExperimental;
            
            ynChar = {'N'; 'Y'; ''};
            agStr = {'Correct'; 'False Essential'; 'False Non-Essential'; 'Not Predicted'};
            
            modelEss = 2 * ones(size(deletionGrowthRates));
            modelEss(geneClasses == SingleGeneDeletions.NON_ESSENTIAL) = 1;
            modelEss(geneClasses == SingleGeneDeletions.UNOBSERVED) = 3;
            
            fbaEss = ones(size(deletionGrowthRates));
            fbaEss(deletionGrowthRates >= essThresh) = 1;
            fbaEss(deletionGrowthRates < essThresh) = 2;
            fbaEss(isnan(deletionGrowthRates)) = 3;
            
            expEss = ones(size(deletionGrowthRatesExperimental));
            expEss(deletionGrowthRatesExperimental >= essThresh) = 1;
            expEss(deletionGrowthRatesExperimental < essThresh) = 2;
            expEss(isnan(deletionGrowthRatesExperimental)) = 3;
            
            modelAgrees = ones(size(deletionGrowthRates));
            modelAgrees(((geneClasses == SingleGeneDeletions.NON_ESSENTIAL) == (deletionGrowthRatesExperimental > essThresh))) = 1;
            modelAgrees(((geneClasses ~= SingleGeneDeletions.NON_ESSENTIAL) & (deletionGrowthRatesExperimental > essThresh))) = 2;
            modelAgrees(((geneClasses == SingleGeneDeletions.NON_ESSENTIAL) & (deletionGrowthRatesExperimental < essThresh))) = 3;
            modelAgrees(geneClasses == SingleGeneDeletions.UNOBSERVED) = 4;
            
            fbaAgrees = ones(size(deletionGrowthRates));
            fbaAgrees(((deletionGrowthRates > essThresh) == (deletionGrowthRatesExperimental > essThresh))) = 1;
            fbaAgrees(((deletionGrowthRates < essThresh) & (deletionGrowthRatesExperimental > essThresh))) = 2;
            fbaAgrees(((deletionGrowthRates > essThresh) & (deletionGrowthRatesExperimental < essThresh))) = 3;
            fbaAgrees(isnan(deletionGrowthRates)) = 4;
            
            geneClassNames = cell(size(geneClasses));
            for i = 1:size(SingleGeneDeletions.DELETION_STRAIN_CLASSES, 1)
                geneClassNames(geneClasses == i) = SingleGeneDeletions.DELETION_STRAIN_CLASSES(i);
            end
            wtGrowthRate = wtGrowthRate * 3600 * sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight) * 1e15;
            deletionGrowthRates = deletionGrowthRates * 3600 * sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight) * 1e15;
            
            content = [content; {
                'WT' ...
                'Wild-type' ...
                '' ...
                '' ...
                '' ...
                '' ...
                'Wild-type' ...
                '' ...
                } ...
                num2cell(reshape(permute(wtStats, [3 2 1]), [], size(wtStats, 1))') ...
                num2cell(wtGrowthRate) ...
                num2cell(wtEnergyProduction) ...
                num2cell(wtDoublingTimeExperimental) ...
                { ...
                '' ...
                '' ...
                '' ...
                '' ...
                }];
            
            [i, j] = find(rna.nascentRNAGeneComposition);
            tuIDs = cell(size(g.wholeCellModelIDs));
            tuNames = cell(size(g.wholeCellModelIDs));
            tuIDs(i) = rna.wholeCellModelIDs(rna.nascentIndexs(j));
            tuNames(i) = rna.names(rna.nascentIndexs(j));
            
            content = [content;
                g.wholeCellModelIDs ...
                g.names ...
                tuIDs ...
                tuNames ...
                geneImplementedMolecules ...
                geneImplementedProcesses ...
                geneClassNames ...
                ynChar(modelEss) ...
                num2cell(reshape(permute(deletionStats, [3 2 1]), [], size(deletionStats, 1))') ...
                num2cell(deletionGrowthRates) ...
                num2cell(deletionEnergyProductions) ...
                num2cell(deletionDoublingTimesExperimental) ...
                ynChar(fbaEss)  ...
                ynChar(expEss) ...
                agStr(modelAgrees) ...
                agStr(fbaAgrees)
                ];
            content(cellfun(@(x) ~isempty(x) && all(isnumeric(x)) && all(isnan(x)), content)) = {[]};
        end
        
        function [content, colLabels] = printSingleGeneDeletionStrainClasses(sim, geneClasses)
            %import classes
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.Simulation;
            import edu.stanford.covert.util.ConstantUtil;
            
            %constants
            g = sim.gene;
            isGeneImplemented = SingleGeneDeletions.getGeneImplementations(sim);
            
            %content
            expEss = ismember(g.essential, {'Y', ''});
            modelEss = geneClasses ~= SingleGeneDeletions.NON_ESSENTIAL;
            
            %overall p-value versus random model of implemented genes. that
            %is p-value that model isn't random model which
            %- Randomly guesses that implemented genes are essential with
            %  probability = sum(tfModelEss & tfImplemented) / nImplemented
            %- Deterministically reports that all non-implemented genes are
            %  non-essential
            nImplemented = sum(isGeneImplemented);
            nCorrectImplemented = sum(modelEss == expEss & isGeneImplemented);
            pCorrectGivenImplemented = ...
                + sum( modelEss & isGeneImplemented) / nImplemented * sum( expEss & isGeneImplemented) / nImplemented  ...
                + sum(~modelEss & isGeneImplemented) / nImplemented * sum(~expEss & isGeneImplemented) / nImplemented;
            totalPValue = 1 - binocdf(nCorrectImplemented - 1, nImplemented, pCorrectGivenImplemented);
            
            %ROC analysis
            rocAnalysis = zeros(size(SingleGeneDeletions.DELETION_STRAIN_CLASSES, 1), 15);
            pValues = zeros(size(SingleGeneDeletions.DELETION_STRAIN_CLASSES, 1), 1);
            for i = 1:size(SingleGeneDeletions.DELETION_STRAIN_CLASSES, 1)
                tfs = geneClasses == i;
                rocAnalysis(i, 1) = sum(tfs);
                rocAnalysis(i, 2:5) = sum([
                    modelEss(tfs) &  expEss(tfs) ...
                    ~modelEss(tfs) & ~expEss(tfs) ...
                    modelEss(tfs) & ~expEss(tfs) ...
                    ~modelEss(tfs) &  expEss(tfs) ...
                    ], 1);
                
                tfs = geneClasses == i  &  isGeneImplemented;
                rocAnalysis(i, 6) = sum(tfs);
                rocAnalysis(i, 7:10) = sum([
                    modelEss(tfs) &  expEss(tfs) ...
                    ~modelEss(tfs) & ~expEss(tfs) ...
                    modelEss(tfs) & ~expEss(tfs) ...
                    ~modelEss(tfs) &  expEss(tfs) ...
                    ], 1);
                
                tfs = geneClasses == i  &  ~isGeneImplemented;
                rocAnalysis(i, 11) = sum(tfs);
                rocAnalysis(i, 12:15) = sum([
                    modelEss(tfs) &  expEss(tfs) ...
                    ~modelEss(tfs) & ~expEss(tfs) ...
                    modelEss(tfs) & ~expEss(tfs) ...
                    ~modelEss(tfs) &  expEss(tfs) ...
                    ], 1);
                
                % p-value versus random model of implemented genes
                tfs = geneClasses == i  &  isGeneImplemented;
                nImplemented = sum(isGeneImplemented(tfs));
                nCorrectImplemented = sum(modelEss(tfs) == expEss(tfs) & isGeneImplemented(tfs));
                if nImplemented == 0
                    pValues(i) = NaN;
                else
                    pValues(i) = 1 - binocdf(nCorrectImplemented - 1, nImplemented, pCorrectGivenImplemented);
                end
            end
            
            %format output
            colLabels = {
                'Class'  ...
                'Total-Count'  'Total-True Essential'  'Total-True Non-Essential'  'Total-False Essential'  'Total-False Non-Essential' ...
                'Implemented-Count'  'Implemented-True Essential'  'Implemented-True Non-Essential'  'Implemented-False Essential'  'Implemented-False Non-Essential' ...
                'Unimplemented-Count'  'Unimplemented-True Essential'  'Unimplemented-True Non-Essential'  'Unimplemented-False Essential'  'Unimplemented-False Non-Essential' ...
                'p-Value' ...
                };
            content = [
                SingleGeneDeletions.DELETION_STRAIN_CLASSES(1:end-1, 1)  num2cell([rocAnalysis(1:end-1, 1)     rocAnalysis(1:end-1, 2:5) ./ rocAnalysis(1:end-1, ones(4, 1))   rocAnalysis(1:end-1, 6)     rocAnalysis(1:end-1, 7:10) ./ rocAnalysis(1:end-1, 6*ones(4, 1))    rocAnalysis(1:end-1, 11)     rocAnalysis(1:end-1, 12:15) ./ rocAnalysis(1:end-1, 11*ones(4, 1))  pValues(1:end-1, :)])
                'Total'                                                  num2cell([sum(rocAnalysis(:, 1), 1)   sum(rocAnalysis(:, 2:5), 1) / sum(rocAnalysis(:, 1), 1)         sum(rocAnalysis(:, 6), 1)   sum(rocAnalysis(:, 7:10), 1) / sum(rocAnalysis(:, 6), 1)            sum(rocAnalysis(:, 11), 1)   sum(rocAnalysis(:, 12:15), 1) / sum(rocAnalysis(:, 11), 1)          totalPValue])
                ];
            content(cellfun(@(el) any(isnan(el)), content)) = {[]};
        end
        
        function plot(sim, deletionMetabolicOutputs, wtMetabolicOutput, deletionStats, wtStats, geneClasses, fileName)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            
            %categories, IDs
            [content, colLabels] = SingleGeneDeletions.printSingleGeneDeletionStrains(sim, deletionMetabolicOutputs, wtMetabolicOutput, deletionStats, wtStats, geneClasses);
            geneIDs = content(2:end, ismember(colLabels, 'Gene'));
            geneCats = content(2:end, ismember(colLabels, 'Model Agree'));
            geneProcesses = content(2:end, ismember(colLabels, 'Implementation Process'));
            wtGrowthRateModel = cell2mat(content(1, ismember(colLabels, 'FBA Growth (fg h^-1)')));
            wtDoublingTimeModel = cell2mat(content(1, ismember(colLabels, 'Mean Model Mass Doubling Duration (h)')));
            wtCellCycleTimeModel = cell2mat(content(1, ismember(colLabels, 'Mean Model Cell Cycle Duration (h)')));
            wtDoublingTimeExp = cell2mat(content(1, ismember(colLabels, 'Mean Experimental Mass Doubling Time (h)')));
            geneGrowthRatesModel = content(2:end, ismember(colLabels, 'FBA Growth (fg h^-1)'));
            geneDoublingTimesModel = content(2:end, ismember(colLabels, 'Mean Model Mass Doubling Duration (h)'));
            geneCellCycleTimesModel = content(2:end, ismember(colLabels, 'Mean Model Cell Cycle Duration (h)'));
            geneDoublingTimesExp = content(2:end, ismember(colLabels, 'Mean Experimental Mass Doubling Time (h)'));
            
            geneGrowthRatesModel(cellfun(@isempty, geneGrowthRatesModel)) = {NaN};
            geneDoublingTimesModel(cellfun(@isempty, geneDoublingTimesModel)) = {NaN};
            geneCellCycleTimesModel(cellfun(@isempty, geneCellCycleTimesModel)) = {NaN};
            geneDoublingTimesExp(cellfun(@isempty, geneDoublingTimesExp)) = {NaN};
            
            geneGrowthRatesModel = cell2mat(geneGrowthRatesModel);
            geneDoublingTimesModel = cell2mat(geneDoublingTimesModel);
            geneCellCycleTimesModel = cell2mat(geneCellCycleTimesModel);
            geneDoublingTimesExp = cell2mat(geneDoublingTimesExp);
            
            geneTfsModelImputed = isnan(geneDoublingTimesModel) & ~isnan(geneGrowthRatesModel);
            geneDoublingTimesModel(isnan(geneDoublingTimesModel)) = wtDoublingTimeModel * ...
                wtGrowthRateModel ./ geneGrowthRatesModel(isnan(geneDoublingTimesModel));
            geneCellCycleTimesModel(isnan(geneCellCycleTimesModel)) = wtCellCycleTimeModel * ...
                wtGrowthRateModel ./ geneGrowthRatesModel(isnan(geneCellCycleTimesModel));
            
            %sort
            geneCats(cellfun(@isempty, geneProcesses)) = {'Uncharacterized'};
            [~, idxs] = sort(geneCats);
            geneIDs = geneIDs(idxs);
            geneCats = geneCats(idxs);
            
            %write svg grid
            categories = unique([geneCats; {'Correct'; 'False Essential'; 'False Non-Essential'; 'Uncharacterized'}]);
            colors = [
                149 179 215
                217 150 148
                195 214 155
                255 255 255];
            fid = fopen(fileName, 'w');
            fwrite(fid, SingleGeneDeletions.plotGrid(categories, geneIDs, geneCats, 720, 450, 11, 16, colors));
            fclose(fid);
            
            %% plot experimental observations vs model predictions
            
            geneTfs = geneDoublingTimesModel > 0 | geneDoublingTimesExp > 0;
            geneTfsExpImputed = geneTfs & geneDoublingTimesExp == wtDoublingTimeExp;
            geneTfsExpObserved = geneTfs & ~geneTfsExpImputed;
            
            labels = {'Imputed', 'Model Imputed', 'Experiment Imputed', 'Observed'};
            
            %doubling time vs doubling time
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            hold(axesHandle, 'on');
            h = zeros(4, 1);
            if any(geneTfsModelImputed & geneTfsExpImputed)
                h(1) = plot(axesHandle, geneDoublingTimesExp(geneTfsModelImputed & geneTfsExpImputed), geneDoublingTimesModel(geneTfsModelImputed & geneTfsExpImputed), '.', 'MarkerSize', 8, 'Color', 'r');
            end
            if any(geneTfsModelImputed & geneTfsExpObserved)
                h(2) = plot(axesHandle, geneDoublingTimesExp(geneTfsModelImputed & geneTfsExpObserved), geneDoublingTimesModel(geneTfsModelImputed & geneTfsExpObserved), '.', 'MarkerSize', 8, 'Color', [0 0.5 0.5]);
            end
            if any(~geneTfsModelImputed & geneTfsExpImputed)
                h(3) = plot(axesHandle, geneDoublingTimesExp(~geneTfsModelImputed & geneTfsExpImputed), geneDoublingTimesModel(~geneTfsModelImputed & geneTfsExpImputed), '.', 'MarkerSize', 8, 'Color', [0.5 0 0.5]);
            end
            if any(~geneTfsModelImputed & geneTfsExpObserved)
                h(4) = plot(axesHandle, geneDoublingTimesExp(~geneTfsModelImputed & geneTfsExpObserved), geneDoublingTimesModel(~geneTfsModelImputed & geneTfsExpObserved), '.', 'MarkerSize', 8, 'Color', 'b');
            end
            labels = labels(h ~= 0);
            h = h(h ~= 0);
            legend(h(end:-1:1), labels(end:-1:1), 'Location', 'SouthEast');
            box(axesHandle, 'on');
            xlim(axesHandle, [6 20]);
            ylim(axesHandle, [6 20]);
            xlabel(axesHandle, 'Observed Mass Doubling Time (h)', 'FontSize', 14);
            ylabel(axesHandle, 'Predicted Mass Doubling Time (h)', 'FontSize', 14);
            saveas(figHandle, [baseDir 'ModelDoublingTimeVsExperimentDoublingTime.pdf']);
            close(figHandle);
            
            %cell cycle length vs doubling time
            [axesHandle, figHandle] = PlotUtil.newAxesHandle();
            hold(axesHandle, 'on');
            h = zeros(4, 1);
            if any(geneTfsModelImputed & geneTfsExpImputed)
                plot(axesHandle, geneDoublingTimesExp(geneTfsModelImputed & geneTfsExpImputed), geneCellCycleTimesModel(geneTfsModelImputed & geneTfsExpImputed), '.', 'MarkerSize', 8, 'Color', 'r')
            end
            if any(geneTfsModelImputed & geneTfsExpObserved)
                plot(axesHandle, geneDoublingTimesExp(geneTfsModelImputed & geneTfsExpObserved), geneCellCycleTimesModel(geneTfsModelImputed & geneTfsExpObserved), '.', 'MarkerSize', 8, 'Color', [0 0.5 0.5])
            end
            if any(~geneTfsModelImputed & geneTfsExpImputed)
                plot(axesHandle, geneDoublingTimesExp(~geneTfsModelImputed & geneTfsExpImputed), geneCellCycleTimesModel(~geneTfsModelImputed & geneTfsExpImputed), '.', 'MarkerSize', 8, 'Color', [0.5 0 0.5])
            end
            if any(~geneTfsModelImputed & geneTfsExpObserved)
                plot(axesHandle, geneDoublingTimesExp(~geneTfsModelImputed & geneTfsExpObserved), geneCellCycleTimesModel(~geneTfsModelImputed & geneTfsExpObserved), '.', 'MarkerSize', 8, 'Color', 'b')
            end
            labels = labels(h ~= 0);
            h = h(h ~= 0);
            legend(h(end:-1:1), labels(end:-1:1), 'Location', 'SouthEast');
            box(axesHandle, 'on');
            xlim(axesHandle, [6 20]);
            ylim(axesHandle, [6 20]);
            xlabel(axesHandle, 'Observed Mass Doubling Time (h)', 'FontSize', 14);
            ylabel(axesHandle, 'Predicted Cell Cycle Length (h)', 'FontSize', 14);
            saveas(figHandle, [baseDir 'ModelCellCycleLengthVsExperimentDoublingTime.pdf']);
            close(figHandle);
        end
        
        function svg = plotGrid(categories, cellIDs, cellCats, width, height, gridFontSize, legendFontSize, colors)
            %layout grid
            gridWidth = width - 220;
            area = gridWidth * height / numel(cellIDs);
            cols = ceil(gridWidth / sqrt(area));
            cellWidth = gridWidth / cols;
            cellHeight = cellWidth;
            rows = min(floor(height / cellHeight), ceil(numel(cellIDs) / cols));
            
            while rows * cols < numel(cellIDs)
                cols = cols + 1;
                cellWidth = gridWidth / cols;
                cellHeight = cellWidth;
                rows = min(floor(height / cellHeight), ceil(numel(cellIDs) / cols));
            end
            
            %grid
            grid = [];
            grid = [grid sprintf('\t<g tansform="translate(0,0)">\n')];
            for i = 1:numel(cellIDs)
                j = floor((i-1) / cols);
                k = mod((i-1), cols);
                category = regexprep(cellCats{i}, '\W', '_');
                cellID = regexprep(regexprep(cellIDs{i}, 'MG_*', ''), 'rrnA', '');
                grid = [grid sprintf('\t\t<g transform="translate(%f,%f)"><rect transform="translate(%f,%f)" width="%f" height="%f" class="%s"/><text transform="translate(%f,%f)" class="grid %s">%s</text></g>\n', ...
                    k * cellWidth, j * cellHeight, 0, 0, cellWidth, cellHeight, category, cellWidth / 2, cellHeight / 2 + 0.4 * gridFontSize, category, cellID)]; %#ok<*AGROW>
            end
            grid = [grid sprintf('\t</g>')];
            
            %styles, legend
            styles = [];
            legend = sprintf('\t<g transform="translate(%f,%f)">\n', cellWidth * cols + 10, 0);
            for i = 1:numel(categories)
                category = regexprep(categories{i}, '\W', '_');
                r = colors(mod(i-1, size(colors, 1))+1, 1);
                g = colors(mod(i-1, size(colors, 1))+1, 2);
                b = colors(mod(i-1, size(colors, 1))+1, 3);
                
                cellCount = sum(strcmp(cellCats, categories{i}));
                
                if mean(colors(i, :)) < 255 / 2
                    textAlpha = 255;
                else
                    textAlpha = 0;
                end
                styles = [styles sprintf('\t\trect.%s{\n\t\t\tfill:rgb(%d,%d,%d);\n\t\t}\n', category, r, g, b)];
                styles = [styles sprintf('\t\ttext.%s{\n\t\t\tfill:rgb(%d,%d,%d);\n\t\t}\n', category, textAlpha, textAlpha, textAlpha)];
                
                legend = [legend sprintf('\t\t<g transform="translate(%f,%f)">\n', ...
                    0, (i-1) * (legendFontSize * 1.5))];
                legend = [legend sprintf('\t\t\t<rect transform="translate(%f,%f)" width="%f" height="%f" class="%s grid legend" />\n', ...
                    0, 0, legendFontSize, legendFontSize, category)];
                legend = [legend sprintf('\t\t\t<text transform="translate(%f,%f)" class="legend">%s (%d)</text>\n', ...
                    legendFontSize * 1.25, legendFontSize / 2 + 0.4 * legendFontSize, categories{i}, cellCount)];
                legend = [legend sprintf('\t\t</g>\n')];
            end
            legend = [legend sprintf('\t</g>\n')];
            
            %grid lines
            gridLines = [];
            gridLines = [gridLines sprintf('\t<g>\n')];
            gridLines = [gridLines sprintf('\t\t<g>\n')];
            for i = 0:rows
                if i == rows
                    c = numel(cellIDs) - (rows - 1) * cols;
                else
                    c = cols;
                end
                if i == 0 || i == rows
                    style = 'outter';
                else
                    style = 'inner';
                end
                gridLines = [gridLines sprintf('\t\t\t<line x1="%f" y1="%f" x2="%f" y2="%f" class="grid_%s"/>\n', 0, i * cellHeight, cellWidth * c, i * cellHeight, style)];
                if i == rows - 1 && rows * cols ~= numel(cellIDs)
                    c = mod(numel(cellIDs), cols);
                    gridLines = [gridLines sprintf('\t\t\t<line x1="%f" y1="%f" x2="%f" y2="%f" class="grid_outter"/>\n', cellWidth * c, i * cellHeight, cellWidth * cols, i * cellHeight)];
                end
            end
            
            gridLines = [gridLines sprintf('\t\t</g>\n')];
            gridLines = [gridLines sprintf('\t\t<g>\n')];
            for i = 0:cols
                if (rows - 1) * cols + i > numel(cellIDs)
                    r = rows - 1;
                else
                    r = rows;
                end
                if i == 0 || i == cols
                    style = 'outter';
                else
                    style = 'inner';
                end
                gridLines = [gridLines sprintf('\t\t\t<line x1="%f" y1="%f" x2="%f" y2="%f" class="grid_%s"/>\n', i * cellWidth, 0, i * cellWidth, cellHeight * r, style)];
                if (rows - 1) * cols + i == numel(cellIDs) && rows * cols ~= numel(cellIDs)
                    gridLines = [gridLines sprintf('\t\t\t<line x1="%f" y1="%f" x2="%f" y2="%f" class="grid_outter"/>\n', i * cellWidth, cellHeight * (rows-1), i * cellWidth, cellHeight * rows)];
                end
            end
            gridLines = [gridLines sprintf('\t\t</g>\n')];
            gridLines = [gridLines sprintf('\t</g>\n')];
            
            %svg
            svg = [];
            svg = [svg sprintf('<?xml version="1.0" encoding="UTF-8"?>\n')];
            svg = [svg sprintf('<svg xmlns="http://www.w3.org/2000/svg"\n')];
            svg = [svg sprintf('	viewBox="-1 -1 %d %d" width="%d" height="%d">\n', width-1, height-1, width, height)];
            svg = [svg sprintf('	<style type="text/css">\n')];
            svg = [svg sprintf('		text{\n')];
            svg = [svg sprintf('			font-family:Helvetica,Helvetica;\n')];
            svg = [svg sprintf('			fill:rgb(0,0,0);\n')];
            svg = [svg sprintf('		}\n')];
            svg = [svg sprintf('		line.grid_outter{\n')];
            svg = [svg sprintf('			stroke-width:2px;\n')];
            svg = [svg sprintf('			stroke:rgb(0,0,0);\n')];
            svg = [svg sprintf('		}\n')];
            svg = [svg sprintf('		line.grid_inner{\n')];
            svg = [svg sprintf('			stroke-width:0.5px;\n')];
            svg = [svg sprintf('			stroke:rgb(0,0,0);\n')];
            svg = [svg sprintf('		}\n')];
            svg = [svg sprintf('		text.legend{\n')];
            svg = [svg sprintf('			font-size:%dpx;\n', legendFontSize)];
            svg = [svg sprintf('			text-anchor:start;\n')];
            svg = [svg sprintf('		}\n')];
            svg = [svg sprintf('		text.grid{\n')];
            svg = [svg sprintf('			font-size:%dpx;\n', gridFontSize)];
            svg = [svg sprintf('			text-anchor:middle;\n')];
            svg = [svg sprintf('		}\n')];
            svg = [svg sprintf('		rect{\n')];
            svg = [svg sprintf('			fill:none;\n')];
            svg = [svg sprintf('			stroke-width:0px;\n')];
            svg = [svg sprintf('			stroke:none;\n')];
            svg = [svg sprintf('		}\n')];
            svg = [svg sprintf('		rect.legend{\n')];
            svg = [svg sprintf('			stroke-width:2px;\n')];
            svg = [svg sprintf('			stroke:rgb(0,0,0);\n')];
            svg = [svg sprintf('		}\n')];
            svg = [svg sprintf('%s\n', styles)];
            svg = [svg sprintf('	</style>\n\n')];
            svg = [svg sprintf('<!-- Grid -->\n')];
            svg = [svg sprintf('%s\n\n', grid)];
            svg = [svg sprintf('<!-- Grid lines -->\n')];
            svg = [svg sprintf('%s\n\n', gridLines)];
            svg = [svg sprintf('<!-- Legend -->\n')];
            svg = [svg sprintf('%s\n\n', legend)];
            svg = [svg sprintf('</svg>\n')];
        end
    end
    
    %overview of each single gene deletion strain
    methods (Static)
        function [data, simMetaData, validSims] = getSimulationData(sim, simMetaData)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            nSims = numel(simMetaData);
            simEndTimes = [simMetaData.lengthSec]';
            simBatchIDs = {simMetaData.simGroup}';
            simBatchIdxs = [simMetaData.simIdx]';
            
            termOrgIdxs = [
                sim.compartment.terminalOrganelleCytosolIndexs
                sim.compartment.terminalOrganelleMembraneIndexs
                ];
            
            %% get data
            nTimePoints = max(simEndTimes)+1;
            time = (0:nTimePoints-1)' / 3600;
            
            ensembleStateNames = {
                'mass'
                'growth_rate'
                'ploidy'
                'superhelicity'
                'dntps'
                'ntps'
                'amino_acids'
                'antibiotics'
                'mrnas'
                'rrnas'
                'srnas'
                'trnas'
                'matureMonomers'
                'matureComplexs'
                'damagedRnas'
                'damagedProteins'
                'atpUsage'
                'gtpUsage'
                'atpProduction'
                'gtpProduction'
                'pinchedDiameter'
                };
            
            stateNames = {
                'Mass' 'cell'       ':'  ':'
                'Mass' 'dnaWt'      ':'  '-sum'
                'Mass' 'rnaWt'      ':'  '-sum'
                'Mass' 'proteinWt'  ':'  '-sum'
                'Geometry' 'totalLength' ':' ':'
                'Host' 'isBacteriumAdherent' ':' ':'
                };
            
            mass = NaN(nTimePoints, nSims);
            growth = NaN(nTimePoints, nSims);
            atpUsage = NaN(nTimePoints, nSims);
            gtpUsage = NaN(nTimePoints, nSims);
            atpProduction = NaN(nTimePoints, nSims);
            gtpProduction = NaN(nTimePoints, nSims);
            ploidy = NaN(nTimePoints, nSims);
            superhelicalDensity = NaN(nTimePoints, nSims);
            rnas = NaN(nTimePoints, nSims);
            prots = NaN(nTimePoints, nSims);
            dntps = NaN(nTimePoints, nSims);
            ntps = NaN(nTimePoints, nSims);
            aminoAcids = NaN(nTimePoints, nSims);
            damagedRnas = NaN(nTimePoints, nSims);
            damagedProteins = NaN(nTimePoints, nSims);
            antibiotics = NaN(nTimePoints, nSims);
            repInitDuration = NaN(1, nSims);
            repDuration = NaN(1, nSims);
            cytokinesisDuration = NaN(1, nSims);
            massDoublingDuration = NaN(1, nSims);
            cellCycleDuration = NaN(1, nSims);
            cellWt = NaN(nTimePoints, nSims);
            dnaWt = NaN(nTimePoints, nSims);
            rnaWt = NaN(nTimePoints, nSims);
            proteinWt = NaN(nTimePoints, nSims);
            terminalOrganelleWt = NaN(nTimePoints, nSims);
            pinchedDiameter = NaN(nTimePoints, nSims);
            totalLength = NaN(nTimePoints, nSims);
            isBacteriumAdherent = NaN(nTimePoints, nSims);
            
            validSims = false(nSims, 1);
            for i = 1:nSims
                try
                    %summary statistics
                    simStats = SummaryLogger.getSimulationStatistics([SimulationDiskUtil.getBaseDir() filesep simBatchIDs{i}], simBatchIdxs(i)) / 3600;
                    repInitDuration(i) = simStats(SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME);
                    repDuration(i) = simStats(SummaryLogger.SIM_STATUS_INDEX_REP_TIME) - simStats(SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME);
                    cytokinesisDuration(i) = simStats(SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME) - simStats(SummaryLogger.SIM_STATUS_INDEX_REP_TIME);
                    massDoublingDuration(i) = simStats(SummaryLogger.SIM_STATUS_INDEX_MASS_DOUBLING_TIME);
                    cellCycleDuration(i) = simStats(SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME);
                    
                    clear simStats;
                    
                    %summary data
                    ensemble = SimulationEnsemble(simBatchIDs{i}, ensembleStateNames, [], simBatchIdxs(i));
                    dt = ensemble.stateData.downsampleStepSec;
                    timeIdxs = 1:dt:dt * floor(simEndTimes(i) / dt) + 1;
                    if timeIdxs(end) ~= simEndTimes(i)+1
                        timeIdxs = [timeIdxs  simEndTimes(i)+1];
                    end
                    if ensemble.stateData.time(end) > simEndTimes(i)
                        time = (0:ensemble.stateData.time(end))' / 3600;
                    end
                    mass(timeIdxs, i) = permute(ensemble.stateData.values(ensemble.getPropertyIndices('mass'), :, :), [3 1 2]);
                    growth(timeIdxs, i) = permute(ensemble.stateData.values(ensemble.getPropertyIndices('growth_rate'), :, :), [3 1 2]);
                    ploidy(timeIdxs, i) = permute(ensemble.stateData.values(ensemble.getPropertyIndices('ploidy'), :, :), [3 1 2]);
                    superhelicalDensity(timeIdxs, i) = permute(ensemble.stateData.values(ensemble.getPropertyIndices('superhelicity'), :, :), [3 1 2]);
                    dntps(timeIdxs, i) = permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('dntps'), :, :), 1), [3 1 2]);
                    ntps(timeIdxs, i) = permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('ntps'), :, :), 1), [3 1 2]);
                    aminoAcids(timeIdxs, i) = permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('amino_acids'), :, :), 1), [3 1 2]);
                    antibiotics(timeIdxs, i) = permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('antibiotics'), :, :), 1), [3 1 2]);
                    damagedRnas(timeIdxs, i) = permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('damagedRnas'), :, :), 1), [3 1 2]);
                    damagedProteins(timeIdxs, i) = permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('damagedProteins'), :, :), 1), [3 1 2]);
                    rnas(timeIdxs, i) = ...
                        + permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('mrnas'), :, :), 1), [3 1 2]) ...
                        + permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('rrnas'), :, :), 1), [3 1 2]) ...
                        + permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('srnas'), :, :), 1), [3 1 2]) ...
                        + permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('trnas'), :, :), 1), [3 1 2]);
                    prots(timeIdxs, i) = ...
                        + permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('matureMonomers'), :, :), 1), [3 1 2]) ...
                        + permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('matureComplexs'), :, :), 1), [3 1 2]);
                    atpUsage(timeIdxs, i) = permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('atpUsage'), :, :), 1), [3 1 2]);
                    gtpUsage(timeIdxs, i) = permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('gtpUsage'), :, :), 1), [3 1 2]);
                    atpProduction(timeIdxs, i) = permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('atpProduction'), :, :), 1), [3 1 2]);
                    gtpProduction(timeIdxs, i) = permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('gtpProduction'), :, :), 1), [3 1 2]);
                    pinchedDiameter(timeIdxs, i) = permute(sum(ensemble.stateData.values(ensemble.getPropertyIndices('pinchedDiameter'), :, :), 1), [3 1 2]);
                    
                    clear ensemble;
                    
                    states = SimulationEnsemble.load(simBatchIDs{i}, stateNames, [], [], 1, 'extract', simBatchIdxs(i));
                    
                    cellWt(2:simEndTimes(i)+1, i) = permute(sum(states.Mass.cell, 2), [3 1 2]);
                    dnaWt(2:simEndTimes(i)+1, i) = permute(states.Mass.dnaWt, [3 1 2]);
                    rnaWt(2:simEndTimes(i)+1, i) = permute(states.Mass.rnaWt, [3 1 2]);
                    proteinWt(2:simEndTimes(i)+1, i) = permute(states.Mass.proteinWt, [3 1 2]);
                    terminalOrganelleWt(2:simEndTimes(i)+1, i) = permute(sum(states.Mass.cell(:, termOrgIdxs, :), 2), [3 1 2]);
                    totalLength(2:simEndTimes(i)+1, i) = permute(states.Geometry.totalLength, [3 1 2]);
                    isBacteriumAdherent(2:simEndTimes(i)+1, i) = permute(states.Host.isBacteriumAdherent, [3 1 2]);
                    
                    clear states;
                    
                    validSims(i) = true;
                catch exception
                    if exist('simStats', 'var')
                        clear simStats;
                    end
                    if exist('ensemble', 'var')
                        clear ensemble;
                    end
                    warning('WholeCell:warning:loadSingleGeneDeletion', 'Unable to read simulation %s:%d\n%s', simBatchIDs{i}, simBatchIdxs(i), exception.getReport());
                end
            end
            
            if nSims > 0 && ~any(validSims)
                exception.addCause(MException('SingleGeneDeletions:error', 'No valid Simulations')).rethrow();
            end
            
            simEndTimes = simEndTimes(validSims, :);
            mass = mass(:, validSims);
            growth = growth(:, validSims);
            atpUsage = atpUsage(:, validSims);
            gtpUsage = gtpUsage(:, validSims);
            atpProduction = atpProduction(:, validSims);
            gtpProduction = gtpProduction(:, validSims);
            ploidy = ploidy(:, validSims);
            superhelicalDensity = superhelicalDensity(:, validSims);
            rnas = rnas(:, validSims);
            prots = prots(:, validSims);
            dntps = dntps(:, validSims);
            ntps = ntps(:, validSims);
            aminoAcids = aminoAcids(:, validSims);
            damagedRnas = damagedRnas(:, validSims);
            damagedProteins = damagedProteins(:, validSims);
            antibiotics = antibiotics(:, validSims);
            repInitDuration = repInitDuration(:, validSims);
            repDuration = repDuration(:, validSims);
            cytokinesisDuration = cytokinesisDuration(:, validSims);
            massDoublingDuration = massDoublingDuration(:, validSims);
            cellCycleDuration = cellCycleDuration(:, validSims);
            cellWt = cellWt(:, validSims);
            dnaWt = dnaWt(:, validSims);
            rnaWt = rnaWt(:, validSims);
            proteinWt = proteinWt(:, validSims);
            terminalOrganelleWt = terminalOrganelleWt(:, validSims);
            pinchedDiameter = pinchedDiameter(:, validSims);
            totalLength = totalLength(:, validSims);
            isBacteriumAdherent = isBacteriumAdherent(:, validSims);
            
            mass = mass * sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight) * 1e15;
            growth = growth * 3600 * sim.state('Mass').cellInitialDryWeight / (1 - sim.state('Mass').fractionWetWeight) * 1e15;
            
            %% figure data
            data = struct;
            data.simEndTimes = simEndTimes;
            data.time = time;
            data.mass = mass;
            data.growth = growth;
            data.atpUsage = atpUsage;
            data.gtpUsage = gtpUsage;
            data.atpProduction = atpProduction;
            data.gtpProduction = gtpProduction;
            data.ploidy = ploidy;
            data.superhelicalDensity = superhelicalDensity;
            data.rnas = rnas;
            data.prots = prots;
            data.dntps = dntps;
            data.ntps = ntps;
            data.aminoAcids = aminoAcids;
            data.damagedRnas = damagedRnas;
            data.damagedProteins = damagedProteins;
            data.antibiotics = antibiotics;
            data.repInitDuration = repInitDuration;
            data.repDuration = repDuration;
            data.cytokinesisDuration = cytokinesisDuration;
            data.massDoublingDuration = massDoublingDuration;
            data.cellCycleDuration = cellCycleDuration;
            data.cellWt = cellWt;
            data.dnaWt = dnaWt;
            data.rnaWt = rnaWt;
            data.proteinWt = proteinWt;
            data.terminalOrganelleWt = terminalOrganelleWt;
            data.pinchedDiameter = pinchedDiameter;
            data.totalLength = totalLength;
            data.isBacteriumAdherent = isBacteriumAdherent;
        end
        
        function nPlots = plotSingleGeneDeletion(data, classes, classLabels, plotRepresentativeCell, figHandle, plotIdx, colors)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if nargin < 7
                colors = SingleGeneDeletions.COLORS;
            end
            
            maxTime = NaN;
            for i = 1:numel(data)
                maxTime = max(maxTime, data{i}.time(end));
            end
            
            growthDataNamesFields = {
                'mass'                 {'Mass' '(fg)'}
                'growth'               {'Growth' '(fg h^{-1})'}
                'atpUsage'             {'ATP' 'Usage'}
                'gtpUsage'             {'GTP' 'Usage'}
                'ploidy'               {'Chr' 'Copy No'}
                'superhelicalDensity'  {'Super-' 'helicity'}
                'rnas'                 'RNA'
                'prots'                'Protein'
                };
            toxicDataNamesFields = {
                'dntps'            'dNTP'
                'ntps'             'NTP'
                'aminoAcids'       {'Amino' 'Acid'}
                'damagedRnas'      {'Damaged' 'RNA'}
                'damagedProteins'  {'Damaged' 'Protein'}
                'antibiotics'      'Antibiotics'
                };
            phaseDataNamesFields = {
                'repInitDuration'       {'Replication' 'Initiation'}
                'repDuration'           'Replication'
                'cytokinesisDuration'   'Cytokinesis'
                'massDoublingDuration'  'Growth'
                'cellCycleDuration'     'Cell Cycle'
                };
            
            if nargin >= 6 && ~isempty(plotIdx)
                for i = 1:size(growthDataNamesFields, 1)
                    if iscell(growthDataNamesFields{i, 2})
                        growthDataNamesFields{i, 2} = strjoin(' ', growthDataNamesFields{i, 2}{:});
                    end
                end
                for i = 1:size(toxicDataNamesFields, 1)
                    if iscell(toxicDataNamesFields{i, 2})
                        toxicDataNamesFields{i, 2} = strjoin(' ', toxicDataNamesFields{i, 2}{:});
                    end
                end
                for i = 1:size(phaseDataNamesFields, 1)
                    if iscell(phaseDataNamesFields{i, 2})
                        phaseDataNamesFields{i, 2} = strjoin(' ', phaseDataNamesFields{i, 2}{:});
                    end
                end
            end
            
            nPlots = ...
                + size(growthDataNamesFields, 1) ...
                + size(toxicDataNamesFields, 1) ...
                + size(phaseDataNamesFields, 1);
            
            
            %% plot data
            clf(figHandle);
            if nargin >= 6 && ~isempty(plotIdx)
                if plotIdx <= size(growthDataNamesFields, 1) + size(toxicDataNamesFields, 1)
                    [axesHandles, xAxesHandles] = PlotUtil.multiElementPlot(figHandle, {50; 50}, [0 maxTime], struct(...
                        'colWidths', [2.4 0.5], 'yAxesLabelWidths', [0.14 0.01], 'position', [0.02 0.32 0.95 0.6]));
                    delete(xAxesHandles{2});
                else
                    [axesHandles, xAxesHandles] = PlotUtil.multiElementPlot(figHandle, {50}, [0 maxTime], struct(...
                        'yAxesLabelWidths', 0.09, 'position', [0.02 0.32 0.95 0.6]));
                end
            else
                axesSizes = {
                    5 * ones(8, 1);   5 * ones(8, 1);
                    5 * ones(6, 1);   5 * ones(6, 1);
                    5 * ones(5, 1);
                    };
                [axesHandles, xAxesHandles] = PlotUtil.multiElementPlot(figHandle, axesSizes, [0 maxTime], struct(...
                    'colWidths', [2.4 0.5 2.4 0.5 2], 'yAxesLabelWidths', [0.14 0.01 0.14 0.01 0.09], 'position', [0.02 0.32 0.95 0.6]));
                delete(xAxesHandles{2});
                delete(xAxesHandles{4});
            end
            
            %% growth
            if nargin >= 6 && ~isempty(plotIdx) && plotIdx <= size(growthDataNamesFields, 1)
                SingleGeneDeletions.plotSingleGeneDeletionTemporalDynamics(...
                    axesHandles{1}, axesHandles{2}, data, classes, classLabels, ...
                    plotRepresentativeCell, growthDataNamesFields{plotIdx, :}, colors);
            elseif nargin < 6 || isempty(plotIdx)
                for i = 1:size(growthDataNamesFields, 1)
                    SingleGeneDeletions.plotSingleGeneDeletionTemporalDynamics(...
                        axesHandles{1}(i), axesHandles{2}(i), data, classes, classLabels, ...
                        plotRepresentativeCell, growthDataNamesFields{i, :}, colors);
                end
            end
            
            %% toxic build up
            if nargin >= 6 && ~isempty(plotIdx) && plotIdx > size(growthDataNamesFields, 1) && plotIdx <= size(growthDataNamesFields, 1) + size(toxicDataNamesFields, 1)
                tmp = plotIdx - size(growthDataNamesFields, 1);
                SingleGeneDeletions.plotSingleGeneDeletionTemporalDynamics(...
                    axesHandles{1}, axesHandles{2}, data, classes, classLabels, ...
                    plotRepresentativeCell, toxicDataNamesFields{tmp, :}, colors);
            elseif nargin < 6 || isempty(plotIdx)
                for i = 1:size(toxicDataNamesFields, 1)
                    SingleGeneDeletions.plotSingleGeneDeletionTemporalDynamics(...
                        axesHandles{3}(i), axesHandles{4}(i), data, classes, classLabels, ...
                        plotRepresentativeCell, toxicDataNamesFields{i, :}, colors);
                end
            end
            
            %% cell cycle phase durations
            if nargin >= 6 && ~isempty(plotIdx) && plotIdx > (size(growthDataNamesFields, 1) + size(toxicDataNamesFields, 1))
                tmp = plotIdx - size(growthDataNamesFields, 1) - size(toxicDataNamesFields, 1);
                SingleGeneDeletions.plotSingleGeneDeletionCellCyclePhaseHistogram(...
                    axesHandles{1}, data, classes, classLabels, phaseDataNamesFields{tmp, :}, colors);
            elseif nargin < 6 || isempty(plotIdx)
                for i = 1:size(phaseDataNamesFields, 1)
                    SingleGeneDeletions.plotSingleGeneDeletionCellCyclePhaseHistogram(...
                        axesHandles{5}(i), data, classes, classLabels, phaseDataNamesFields{i, :}, colors);
                end
            end
            
            %% formatting
            tmp = [];
            tmp2 = [];
            for i = 1:2:numel(axesHandles)
                tmp = [tmp; axesHandles{i}];
                tmp2 = [tmp2; xAxesHandles{i}];
            end
            
            if nargin >= 6 && ~isempty(plotIdx)
                set(tmp, 'FontSize', 10);
                set(tmp2, 'FontSize', 10);
                set(get(tmp2, 'xlabel'), 'FontSize', 12);
                set(get(tmp, 'ylabel'), 'FontSize', 12, 'Color', [0 0 0]);
                set(get(tmp, 'title'), 'FontSize', 12);
                PlotUtil.alignYAxesLabels(axesHandles);
                axesHandles2 = PlotUtil.offsetYAxes(tmp, 0.02);
                h = [
                    getappdata(tmp, 'YTickText')
                    getappdata(axesHandles2(1), 'YTickText')
                    ];
                set(h, 'FontSize', 10);
            else
                set(tmp, 'FontSize', 6);
                set(tmp2, 'FontSize', 6);
                set(cell2mat(get(tmp2, 'xlabel')), 'FontSize', 8);
                set(cell2mat(get(tmp, 'ylabel')), 'FontSize', 8, 'Color', [0 0 0]);
                set(cell2mat(get(tmp, 'title')), 'FontSize', 8);
                PlotUtil.alignYAxesLabels(axesHandles);
                PlotUtil.offsetYAxes(axesHandles(1:2:end-1), 0.05);
            end
            
            %% legend
            [uClasses, uClassIdxs] = unique(classes);
            h = zeros(size(uClasses));
            nCols = max(5, ceil(numel(uClasses) / 5));
            nRows = 5;
            for i = 1:numel(uClasses)
                x = 0.11 + 0.78 / nCols * (ceil(i / nRows) - 1);
                y = 0.32 - 0.10 - 0.03 * (mod(i - 1, nRows) + 1);
                h(i) = annotation(figHandle, 'textarrow', x + [0 -0.025], y(1, [1 1]), ...
                    'HeadStyle', 'none', ...
                    'String', classLabels{uClassIdxs(i), 1}, ...
                    'FontSize', 6, ...
                    'LineWidth', 4, ...
                    'Color', colors(uClasses(i), :), ...
                    'TextColor', 'k');
                setappdata(h(i), 'class', uClasses(i));
                setappdata(h(i), 'label', classLabels(uClassIdxs(i), 1));
                set(h(i), 'ButtonDownFcn', @SingleGeneDeletions.plotSingleGeneDeletionCallback);
            end
            
            setappdata(figHandle, 'legendEntries', h);
        end
        
        function plotSingleGeneDeletionTemporalDynamics(dynamicsAxesHandle, histogramAxesHandle, data, classes, classLabels, plotRepresentativeCell, dataField, dataLabel, colors)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            %options
            if nargin < 9
                colors = SingleGeneDeletions.COLORS;
            end
            
            %axes
            cla(dynamicsAxesHandle);
            cla(histogramAxesHandle);
            hold(dynamicsAxesHandle, 'on');
            hold(histogramAxesHandle, 'on');
            
            yRange = [NaN NaN];
            maxTime = NaN;
            for i = 1:numel(data)
                yRange(1) = min(yRange(1), min(data{i}.(dataField)(:)));
                yRange(2) = max(yRange(2), max(data{i}.(dataField)(:)));
                maxTime = max(maxTime, data{i}.time(end));
            end
            maxTime = min(maxTime, 20);
            if yRange(1) == yRange(2)
                yRange(1) = yRange(1) - sqrt(eps);
                yRange(2) = yRange(2) + sqrt(eps);
            end
            
            simEndTimes = cell(size(data));
            for i = 1:numel(data)
                simEndTimes{i} = zeros(size(data{i}.simEndTimes));
                for j = 1:numel(data{i}.simEndTimes);
                    simEndTimes{i}(j) = find(~isnan(data{i}.growth(:, j)), 1, 'last');
                end
            end
            
            %dynamics
            for i = 1:numel(data)
                if plotRepresentativeCell(i)
                    medianTime = median(simEndTimes{i});
                    [~, simIdx] = min(abs(medianTime - simEndTimes{i}));
                else
                    simIdx = 1:numel(simEndTimes{i});
                end
                xdata = data{i}.time;
                ydata = data{i}.(dataField)(:, simIdx);
                xdata = xdata(any(~isnan(ydata), 2), :);
                ydata = ydata(any(~isnan(ydata), 2), :);
                
                h = plot(dynamicsAxesHandle, xdata, ydata);
                set(h, 'Color', colors(classes(i), :));
                for j = 1:numel(h)
                    setappdata(h(j), 'class', i);
                    setappdata(h(j), 'label', classLabels(i, :));
                end
                set(h, 'ButtonDownFcn', @SingleGeneDeletions.plotSingleGeneDeletionCallback);
            end
            
            %marginal
            edges = linspace(yRange(1), yRange(2), 20);
            uClasses = unique(classes);
            freq = zeros(numel(edges), numel(uClasses));
            for i = 1:numel(data)
                tmp = data{i}.(dataField)(sub2ind(size(data{i}.(dataField)), simEndTimes{i}, (1:numel(simEndTimes{i}))'));
                freq(:, classes(i) == uClasses) = ...
                    + freq(:, classes(i) == uClasses) ...
                    + reshape(histc(tmp, edges), [], 1);
            end
            freq = freq ./ repmat(sum(freq, 1), [numel(edges) 1]);
            edges = [edges(1);  edges(:) + diff(edges(1:2)) / 2;  edges(end) + diff(edges(1:2))];
            h = plot(histogramAxesHandle, [zeros(1, numel(uClasses)); freq; zeros(1, numel(uClasses))], edges, 'b');
            
            for i = 1:numel(uClasses)
                set(h(i), 'Color', colors(uClasses(i), :));
                setappdata(h(i), 'class', uClasses(i));
                setappdata(h(i), 'label', classLabels(uClasses(i), 1));
                set(h(i), 'ButtonDownFcn', @SingleGeneDeletions.plotSingleGeneDeletionCallback);
            end
            
            %axes
            xlim(dynamicsAxesHandle, [0 maxTime]);
            xlim(histogramAxesHandle, [0 max(1, max(freq(:)))]);
            ylabel(dynamicsAxesHandle, dataLabel, 'FontSize', 8);
            
            if any(yRange) && range(yRange)
                ylim(dynamicsAxesHandle, yRange);
                ylim(histogramAxesHandle, yRange);
                set(dynamicsAxesHandle, 'YTick', yRange);
                PlotUtil.formatYAxes(dynamicsAxesHandle, yRange(1), yRange(end));
            end
            
            figHandle = get(dynamicsAxesHandle, 'Parent');
            
            box(dynamicsAxesHandle, 'off');
            set(dynamicsAxesHandle, ...
                'FontSize', 7, ...
                'TickDir', 'out', ...
                'YMinorTick', 'off', ...
                'XTick', [], ...
                'XColor', get(figHandle, 'Color'));
            
            box(histogramAxesHandle, 'off');
            set(histogramAxesHandle, ...
                'Visible', 'off', ...
                'FontSize', 7, ...
                'TickDir', 'out', ...
                'XMinorTick', 'off', ...
                'YMinorTick', 'off', ...
                'XTick', [], ...
                'YTick', [], ...
                'XColor', get(figHandle, 'Color'), ...
                'YColor', get(figHandle, 'Color'));
        end
        
        function plotSingleGeneDeletionCellCyclePhaseHistogram(axesHandle, data, classes, classLabels, dataField, dataLabel, colors)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if nargin < 7
                colors = SingleGeneDeletions.COLORS;
            end
            
            cla(axesHandle);
            
            maxTime = NaN;
            for i = 1:numel(data)
                maxTime = max(maxTime, data{i}.time(end));
            end
            
            edges = linspace(0, maxTime, 20);
            
            uClasses = unique(classes);
            freq = zeros(numel(edges), numel(uClasses));
            for i = 1:numel(data)
                freq(:, classes(i) == uClasses) = ...
                    + freq(:, classes(i) == uClasses) ...
                    + reshape(histc(data{i}.(dataField), edges), [], 1);
            end
            freq = freq ./ repmat(sum(freq, 1), [numel(edges) 1]);
            
            edges = [edges(1);  edges(:) + diff(edges(1:2)) / 2;  edges(end) + diff(edges(1:2))];
            h = plot(axesHandle, edges, [zeros(1, numel(uClasses)); freq; zeros(1, numel(uClasses))], 'b');
            for i = 1:numel(uClasses)
                set(h(i), 'Color', colors(uClasses(i), :));
                setappdata(h(i), 'class', uClasses(i));
                setappdata(h(i), 'label', classLabels(uClasses(i), 1));
                set(h(i), 'ButtonDownFcn', @SingleGeneDeletions.plotSingleGeneDeletionCallback);
            end
            
            xlim(axesHandle, [0 maxTime]);
            box(axesHandle, 'off');
            set(axesHandle, ...
                'FontSize', 7, ...
                'TickDir', 'out', ...
                'XMinorTick', 'off', ...
                'YMinorTick', 'off', ...
                'XTick', [], ...
                'YTick', [], ...
                'XColor', get(get(axesHandle, 'Parent'), 'Color'), ...
                'YColor', get(get(axesHandle, 'Parent'), 'Color'));
            ylabel(axesHandle, [dataLabel {'Freq'}], 'FontSize', 8);
        end
        
        function plotSingleGeneDeletionCallback(hObject, ~)
            class = getappdata(hObject, 'class');
            label = getappdata(hObject, 'label');
            
            figHandle = get(get(hObject, 'Parent'), 'Parent');
            axesHandles = findobj(get(figHandle, 'Children'), 'type', 'axes');
            for i = 1:numel(axesHandles)
                lineHandles = findobj(get(axesHandles(i), 'Children'), 'type', 'line');
                for j = 1:numel(lineHandles)
                    tmpClass = getappdata(lineHandles(j), 'class');
                    if ~isempty(tmpClass)
                        if tmpClass == class
                            if ~isequal(get(lineHandles(j), 'MarkerEdgeColor'), 'auto')
                                set(lineHandles(j), 'Color', get(lineHandles(j), 'MarkerEdgeColor'));
                            end
                            set(lineHandles(j), 'LineWidth', 2);
                            uistack(lineHandles(j), 'top');
                        else
                            if isequal(get(lineHandles(j), 'MarkerEdgeColor'), 'auto')
                                set(lineHandles(j), 'MarkerEdgeColor', get(lineHandles(j), 'Color'));
                            end
                            set(lineHandles(j), 'LineWidth', 1, 'Color', [0.9 0.9 0.9]);
                        end
                    end
                end
            end
            
            legendEntries = getappdata(figHandle, 'legendEntries');
            for i = 1:numel(legendEntries)
                tmpClass = getappdata(legendEntries(i), 'class');
                if tmpClass == class
                    if ~isequal(get(legendEntries(i), 'TextEdgeColor'), 'none')
                        set(legendEntries(i), 'Color', get(legendEntries(i), 'TextEdgeColor'));
                    end
                    uistack(legendEntries(i), 'top');
                else
                    if isequal(get(legendEntries(i), 'TextEdgeColor'), 'none')
                        set(legendEntries(i), 'TextEdgeColor', get(legendEntries(i), 'Color'), 'TextLineWidth', 0);
                    end
                    set(legendEntries(i), 'Color', [0.9 0.9 0.9]);
                end
                set(legendEntries(i), 'TextColor', 'k');
            end
            
            label = strjoin(' - ', label{:});
            label = strrep(label, '{\Delta}', char(916));
            label = strrep(label, '\_', '_');
            disp(['Selected: ', label]);
        end
        
        function [class, classLabel] = classifyDeletionStrain(geneID, wtData)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            baseDir = [SimulationDiskUtil.getBaseDir() filesep 'singleGeneDeletions' filesep];
            deletionData = load([baseDir geneID '.mat']);
            if nargin < 2
                wtData = load([baseDir 'WT.mat']);
            end
            
            nWtSims = numel(wtData.simEndTimes);
            nDeletionSims = numel(deletionData.simEndTimes);
            
            cellCycleLength = min([wtData.simEndTimes; deletionData.simEndTimes]);
            
            wtStartTimes = 2 * ones(size(wtData.simEndTimes));
            wtTauTimes = zeros(size(wtData.growth, 2), 1);
            wtFinTimes = wtData.simEndTimes + 1;
            deletionStartTimes = 101 * ones(size(deletionData.simEndTimes));
            deletionTauTimes = zeros(size(deletionData.growth, 2), 1);
            deletionFinTimes = deletionData.simEndTimes + 1;
            for i = 1:nWtSims
                wtTauTimes(i) = find(~isnan(wtData.growth(1:cellCycleLength, i)), 1, 'last');
            end
            for i = 1:nDeletionSims
                deletionTauTimes(i) = find(~isnan(deletionData.growth(1:cellCycleLength, i)), 1, 'last');
            end
            
            wtStartInds = sub2ind(size(wtData.growth), wtStartTimes', 1:nWtSims);
            wtFinInds = sub2ind(size(wtData.growth), wtFinTimes', 1:nWtSims); %#ok<NASGU>
            wtTauInds = sub2ind(size(wtData.growth), wtTauTimes', 1:nWtSims);
            deletionStartInds = sub2ind(size(deletionData.growth), deletionStartTimes', 1:nDeletionSims);
            deletionFinInds = sub2ind(size(deletionData.growth), deletionFinTimes', 1:nDeletionSims);
            deletionTauInds = sub2ind(size(deletionData.growth), deletionTauTimes', 1:nDeletionSims);
            
            wtTauGrowths = wtData.growth(wtTauInds);
            deletionFinGrowths = deletionData.growth(deletionFinInds);
            deletionTauGrowths = deletionData.growth(deletionTauInds);
            deletionInitGrowths = deletionData.growth(deletionStartInds);
            
            if ...
                    median(deletionInitGrowths) < 0.1 * median(wtData.growth(wtStartInds)) && (...
                    median(deletionData.damagedProteins(deletionFinInds)) > max(100, 100 * median(wtData.damagedProteins(wtTauInds))) || ...
                    median(deletionData.damagedRnas(deletionFinInds))     > max(100, 100 * median(wtData.damagedRnas(wtTauInds)))     || ...
                    median(deletionData.aminoAcids(deletionFinInds))      > max(100, 100 * median(wtData.aminoAcids(wtTauInds)))      || ...
                    median(deletionData.ntps(deletionFinInds))            > max(100, 100 * median(wtData.ntps(wtTauInds)))               ...
                    )
                class = SingleGeneDeletions.DECOMPOSING;
            elseif median(deletionInitGrowths) < 0.1 * median(wtData.growth(wtStartInds))
                class = SingleGeneDeletions.NON_GROWING;
            elseif all(deletionFinGrowths < deletionInitGrowths) && all(deletionData.rnas(deletionTauInds) < deletionData.rnas(deletionStartInds))
                class = SingleGeneDeletions.DECAYING_GROWTH_NON_RNA_PROTEIN;
            elseif all(deletionFinGrowths < deletionInitGrowths) && all(deletionData.prots(deletionTauInds) < deletionData.prots(deletionStartInds))
                class = SingleGeneDeletions.DECAYING_GROWTH_NON_PROTEIN;
            elseif ...
                    median(deletionData.damagedProteins(deletionFinInds)) > max(100, 100 * median(wtData.damagedProteins(wtTauInds))) || ...
                    median(deletionData.damagedRnas(deletionFinInds))     > max(100, 100 * median(wtData.damagedRnas(wtTauInds)))     || ...
                    median(deletionData.aminoAcids(deletionFinInds))      > max(100, 100 * median(wtData.aminoAcids(wtTauInds)))      || ...
                    median(deletionData.ntps(deletionFinInds))            > max(100, 100 * median(wtData.ntps(wtTauInds)))            || ...
                    median(deletionData.antibiotics(deletionFinInds))     > max(100, 100 * median(wtData.antibiotics(wtTauInds)))
                class = SingleGeneDeletions.TOXIN_ACCUMULATION;
            elseif ~any(deletionData.ploidy(:) > 1)
                class = SingleGeneDeletions.NON_REPLICATIVE;
            elseif median(deletionTauGrowths) < 0.80 * median(wtTauGrowths) && median(deletionTauGrowths) > median(deletionInitGrowths)
                class = SingleGeneDeletions.SLOW_GROWING;
            elseif ~any(deletionData.cytokinesisDuration(:))
                class = SingleGeneDeletions.NON_FISSIVE;
            elseif ~any(deletionData.isBacteriumAdherent(2, :))
                class = SingleGeneDeletions.NO_TERMINAL_ORGANELLE;
            elseif ismember(geneID, SingleGeneDeletions.nonPerpetuatingEssentialGenes)
                class = SingleGeneDeletions.NON_PERPETUATING;
            else
                class = SingleGeneDeletions.NON_ESSENTIAL;
            end
            
            classLabel = SingleGeneDeletions.DELETION_STRAIN_CLASSES{class};
        end
    end
    
    %Analyze specific single gene deletion strains
    methods (Static)
        function figData = analyzeSingleGeneDeletionStrain_MG_210(figHandle, geneID, deletionData, wtData, figData)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %% constants
            perturbations = struct;
            perturbations.geneticKnockouts = {geneID};
            
            sim = CachedSimulationObjectUtil.load();
            sim.applyOptions('lengthSec', 0);
            sim.applyOptions(perturbations);
            
            g = sim.gene;
            mass = sim.state('Mass');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            mr = sim.state('MetabolicReaction');
            ppII = sim.process('ProteinProcessingII');
            
            pcComp = sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3);
            monIdxs = [ppII.lipoproteinMonomerIndexs; ppII.secretedMonomerIndexs];
            unprocessedMonIdxs = ppII.unprocessedMonomerIndexs;
            cpxIdxs = find(any(pcComp(monIdxs, :), 1));
            
            cellInitialDryWeight = mass.cellInitialDryWeight;
            rnaExp = r.expression;
            rnaDecayRates = r.decayRates;
            monDecayRates = pm.decayRates;
            cpxDecayRates = pc.decayRates;
            
            %% data
            maxTime = NaN;
            data = {deletionData, wtData};
            for i = 1:numel(data)
                maxTime = max(maxTime, data{i}.time(end));
            end
            
            if nargin >= 5 && ~isempty(figData)
                nSamples = figData.nSamples;
                nGenerations = figData.nGenerations;
                generations = figData.generations;
            else
                nSamples = 20;
                nGenerations = 10;
                generations = 1:nGenerations;
                
                figData = struct;
                figData.data = cell(nGenerations, 1);
                figData.nSamples = nSamples;
                figData.nGenerations = nGenerations;
                figData.generations = generations;
                
                for i = 1:nGenerations
                    growths = zeros(4, nSamples);
                    for j = 1:5
                        for k = 1:nSamples
                            r.expression = rnaExp;
                            r.decayRates = rnaDecayRates;
                            pm.decayRates = monDecayRates;
                            pc.decayRates = cpxDecayRates;
                            mass.cellInitialDryWeight = cellInitialDryWeight * exp(log(2) * (j - 1) / 4);
                            mr.initialGrowthFilterWidth = Inf;
                            
                            sim.applyOptions('seed', i * nGenerations + j * 5 +  k);
                            sim.allocateMemoryForState(1);
                            sim.initializeState();
                            sim.applyPerturbations();
                            
                            mons = pm.counts;
                            cpxs = pc.counts;
                            
                            pm.counts(pm.processedIIIndexs(monIdxs),    :) = pm.randStream.random('binomial', pm.counts(pm.processedIIIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.signalSequenceIndexs(monIdxs), :) = pm.randStream.random('binomial', pm.counts(pm.signalSequenceIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.foldedIndexs(monIdxs),         :) = pm.randStream.random('binomial', pm.counts(pm.foldedIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.matureIndexs(monIdxs),         :) = pm.randStream.random('binomial', pm.counts(pm.matureIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.inactivatedIndexs(monIdxs),    :) = pm.randStream.random('binomial', pm.counts(pm.inactivatedIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.boundIndexs(monIdxs),          :) = pm.randStream.random('binomial', pm.counts(pm.boundIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.misfoldedIndexs(monIdxs),      :) = pm.randStream.random('binomial', pm.counts(pm.misfoldedIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.damagedIndexs(monIdxs),        :) = pm.randStream.random('binomial', pm.counts(pm.damagedIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            
                            pc.counts(pc.nascentIndexs(cpxIdxs),     :) = pc.randStream.random('binomial', pc.counts(pc.nascentIndexs(cpxIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pc.counts(pc.matureIndexs(cpxIdxs),      :) = pc.randStream.random('binomial', pc.counts(pc.matureIndexs(cpxIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pc.counts(pc.inactivatedIndexs(cpxIdxs), :) = pc.randStream.random('binomial', pc.counts(pc.inactivatedIndexs(cpxIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pc.counts(pc.boundIndexs(cpxIdxs),       :) = pc.randStream.random('binomial', pc.counts(pc.boundIndexs(cpxIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pc.counts(pc.misfoldedIndexs(cpxIdxs),   :) = pc.randStream.random('binomial', pc.counts(pc.misfoldedIndexs(cpxIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pc.counts(pc.damagedIndexs(cpxIdxs),     :) = pc.randStream.random('binomial', pc.counts(pc.damagedIndexs(cpxIdxs), :), 1 / (2.^(generations(i) - 1)));
                            
                            diffMons = ...
                                + mons(pm.processedIIIndexs,    :) - pm.counts(pm.processedIIIndexs, :) ...
                                + mons(pm.signalSequenceIndexs, :) - pm.counts(pm.signalSequenceIndexs, :) ...
                                + mons(pm.foldedIndexs,         :) - pm.counts(pm.foldedIndexs, :) ...
                                + mons(pm.matureIndexs,         :) - pm.counts(pm.matureIndexs, :) ...
                                + mons(pm.inactivatedIndexs,    :) - pm.counts(pm.inactivatedIndexs, :) ...
                                + mons(pm.boundIndexs,          :) - pm.counts(pm.boundIndexs, :) ...
                                + mons(pm.misfoldedIndexs,      :) - pm.counts(pm.misfoldedIndexs, :) ...
                                + mons(pm.damagedIndexs,        :) - pm.counts(pm.damagedIndexs, :);
                            diffCpxs = ...
                                + cpxs(pc.nascentIndexs,     :) - pc.counts(pc.nascentIndexs, :) ...
                                + cpxs(pc.matureIndexs,      :) - pc.counts(pc.matureIndexs, :) ...
                                + cpxs(pc.inactivatedIndexs, :) - pc.counts(pc.inactivatedIndexs, :) ...
                                + cpxs(pc.boundIndexs,       :) - pc.counts(pc.boundIndexs, :) ...
                                + cpxs(pc.misfoldedIndexs,   :) - pc.counts(pc.misfoldedIndexs, :) ...
                                + cpxs(pc.damagedIndexs,     :) - pc.counts(pc.damagedIndexs, :);
                            pm.counts(pm.processedIIndexs, :) = ...
                                + pm.counts(pm.processedIIndexs, :) ...
                                + diffMons;
                            pm.counts(pm.processedIIndexs(monIdxs), :) = ...
                                + pm.counts(pm.processedIIndexs(monIdxs), :) ...
                                + pcComp(monIdxs, :) * diffCpxs;
                            pm.counts(pm.matureIndexs(unprocessedMonIdxs), :) = ...
                                + pm.counts(pm.matureIndexs(unprocessedMonIdxs), :) ...
                                + pcComp(unprocessedMonIdxs, :) * diffCpxs;
                            
                            sim.evolveState();
                            growths(j, k) = mr.growth;
                        end
                    end
                    figData.data{i} = struct;
                    figData.data{i}.time = cumsum([ 0; log(2) ./ mean(growths(1:end-1, :), 2) .* exp(log(2) * (0:j-2)'/4) / 4]) / 3600;
                    figData.data{i}.growth = mean(growths, 2) * 3600 * cellInitialDryWeight / (1 - mass.fractionWetWeight) * 1e15;
                    figData.data{i}.growths = growths * 3600 * cellInitialDryWeight / (1 - mass.fractionWetWeight) * 1e15;
                    figData.data{i}.simEndTimes = j;
                end
            end
            
            %% plot data
            clf(figHandle);
            
            options = struct;
            options.colWidths = [2.4 0.5 2.4];
            options.yAxesLabelWidths = [0.14 0.01 0.01];
            options.xlabelStr = {'Time (h)', [], 'Generation'};
            [axesHandles, xAxesHandles] = PlotUtil.multiElementPlot(figHandle, {50; 50; 50}, {[0 maxTime] [0 maxTime] generations([1 end])}, options);
            delete(xAxesHandles{2});
            
            colors = PlotUtil.getRedGreenColorOrder(1:nGenerations);
            colors = [
                colors(1, :)
                0 1 1
                colors(1:end, :)
                ];
            classLabels = [
                {'Generation 1 (Full Model)'}
                cellfun(@(i) ['Generation ' num2str(i)], num2cell((1:nGenerations)'), 'UniformOutput', false)
                {'Wild-type'}
                ];
            SingleGeneDeletions.plotSingleGeneDeletionTemporalDynamics(axesHandles{1}(1), axesHandles{2}(1), ...
                [data figData.data'], 1:nGenerations+2, classLabels, [false; true; false(size(figData.data))], 'growth', 'Growth (fg h^{-1})', colors);
            
            axesHandle = axesHandles{3}(1);
            growths = [];
            groups = [];
            for i = 1:nGenerations
                growths = [growths figData.data{i}.growths(1, :)];
                groups = [groups repmat(generations(i), 1, nSamples)];
            end
            boxplot(growths, groups, 'Parent', axesHandle);
            ylim(axesHandle, ylim(axesHandles{1}(1)));
            set(axesHandle, 'YTick', []);
        end
        
        %methionine aminopeptidase
        function figData = analyzeSingleGeneDeletionStrain_MG_172(figHandle, geneID, deletionData, wtData, figData)
            import edu.stanford.covert.cell.sim.analysis.SingleGeneDeletions;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            %% constants
            perturbations = struct;
            perturbations.geneticKnockouts = {geneID};
            
            sim = CachedSimulationObjectUtil.load();
            sim.applyOptions('lengthSec', 0);
            sim.applyOptions(perturbations);
            
            g = sim.gene;
            mass = sim.state('Mass');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            mr = sim.state('MetabolicReaction');
            ppI = sim.process('ProteinProcessingI');
            met = sim.process('Metabolism');
            
            sim.setForTest('processesInInitOrder', sim.processesInInitOrder(~ismember(...
                cellfun(@(p) p.name, sim.processesInInitOrder, 'UniformOutput', false), {...
                'Transcription'
                'FtsZPolymerization'
                'ChromosomeCondensation'
                'Translation'
                'DNASupercoiling'
                'ReplicationInitiation'
                })));
            
            pcComp = sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3);
            monIdxs = ppI.nascentMonomerNTerminalMethionineCleavages;
            unprocessedMonIdxs = ~ppI.nascentMonomerNTerminalMethionineCleavages;
            cpxIdxs = find(any(pcComp(monIdxs, :), 1));
            
            cellInitialDryWeight = mass.cellInitialDryWeight;
            rnaExp = r.expression;
            rnaDecayRates = r.decayRates;
            monDecayRates = pm.decayRates;
            cpxDecayRates = pc.decayRates;
            
            %% data
            maxTime = NaN;
            data = {deletionData, wtData};
            for i = 1:numel(data)
                maxTime = max(maxTime, data{i}.time(end));
            end
            
            if nargin >= 5 && ~isempty(figData)
                nSamples = figData.nSamples;
                nGenerations = figData.nGenerations;
                generations = figData.generations;
            else
                nSamples = 20;
                nGenerations = 10;
                generations = 1:nGenerations;
                
                figData = struct;
                figData.data = cell(nGenerations, 1);
                figData.nSamples = nSamples;
                figData.nGenerations = nGenerations;
                figData.generations = generations;
                
                for i = 1:nGenerations
                    growths = zeros(4, nSamples);
                    for j = 1:5
                        for k = 1:nSamples
                            r.expression = rnaExp;
                            r.decayRates = rnaDecayRates;
                            pm.decayRates = monDecayRates;
                            pc.decayRates = cpxDecayRates;
                            mass.cellInitialDryWeight = cellInitialDryWeight * exp(log(2) * (j - 1) / 4);
                            mr.initialGrowthFilterWidth = Inf;
                                                        
                            sim.applyOptions('seed', i * nGenerations + j * 5 +  k);
                            sim.allocateMemoryForState(1);
                            sim.initializeState();
                            sim.applyPerturbations();
                            
                            mons = pm.counts;
                            cpxs = pc.counts;
                            
                            pm.counts(pm.processedIIIndexs(monIdxs),    :) = pm.randStream.random('binomial', pm.counts(pm.processedIIIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.signalSequenceIndexs(monIdxs), :) = pm.randStream.random('binomial', pm.counts(pm.signalSequenceIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.foldedIndexs(monIdxs),         :) = pm.randStream.random('binomial', pm.counts(pm.foldedIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.matureIndexs(monIdxs),         :) = pm.randStream.random('binomial', pm.counts(pm.matureIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.inactivatedIndexs(monIdxs),    :) = pm.randStream.random('binomial', pm.counts(pm.inactivatedIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.boundIndexs(monIdxs),          :) = pm.randStream.random('binomial', pm.counts(pm.boundIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.misfoldedIndexs(monIdxs),      :) = pm.randStream.random('binomial', pm.counts(pm.misfoldedIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pm.counts(pm.damagedIndexs(monIdxs),        :) = pm.randStream.random('binomial', pm.counts(pm.damagedIndexs(monIdxs), :), 1 / (2.^(generations(i) - 1)));
                            
                            pc.counts(pc.nascentIndexs(cpxIdxs),     :) = pc.randStream.random('binomial', pc.counts(pc.nascentIndexs(cpxIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pc.counts(pc.matureIndexs(cpxIdxs),      :) = pc.randStream.random('binomial', pc.counts(pc.matureIndexs(cpxIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pc.counts(pc.inactivatedIndexs(cpxIdxs), :) = pc.randStream.random('binomial', pc.counts(pc.inactivatedIndexs(cpxIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pc.counts(pc.boundIndexs(cpxIdxs),       :) = pc.randStream.random('binomial', pc.counts(pc.boundIndexs(cpxIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pc.counts(pc.misfoldedIndexs(cpxIdxs),   :) = pc.randStream.random('binomial', pc.counts(pc.misfoldedIndexs(cpxIdxs), :), 1 / (2.^(generations(i) - 1)));
                            pc.counts(pc.damagedIndexs(cpxIdxs),     :) = pc.randStream.random('binomial', pc.counts(pc.damagedIndexs(cpxIdxs), :), 1 / (2.^(generations(i) - 1)));
                            
                            diffMons = ...
                                + mons(pm.processedIIIndexs,    :) - pm.counts(pm.processedIIIndexs, :) ...
                                + mons(pm.signalSequenceIndexs, :) - pm.counts(pm.signalSequenceIndexs, :) ...
                                + mons(pm.foldedIndexs,         :) - pm.counts(pm.foldedIndexs, :) ...
                                + mons(pm.matureIndexs,         :) - pm.counts(pm.matureIndexs, :) ...
                                + mons(pm.inactivatedIndexs,    :) - pm.counts(pm.inactivatedIndexs, :) ...
                                + mons(pm.boundIndexs,          :) - pm.counts(pm.boundIndexs, :) ...
                                + mons(pm.misfoldedIndexs,      :) - pm.counts(pm.misfoldedIndexs, :) ...
                                + mons(pm.damagedIndexs,        :) - pm.counts(pm.damagedIndexs, :);
                            diffCpxs = ...
                                + cpxs(pc.nascentIndexs,     :) - pc.counts(pc.nascentIndexs, :) ...
                                + cpxs(pc.matureIndexs,      :) - pc.counts(pc.matureIndexs, :) ...
                                + cpxs(pc.inactivatedIndexs, :) - pc.counts(pc.inactivatedIndexs, :) ...
                                + cpxs(pc.boundIndexs,       :) - pc.counts(pc.boundIndexs, :) ...
                                + cpxs(pc.misfoldedIndexs,   :) - pc.counts(pc.misfoldedIndexs, :) ...
                                + cpxs(pc.damagedIndexs,     :) - pc.counts(pc.damagedIndexs, :);
                            pm.counts(pm.processedIIndexs, :) = ...
                                + pm.counts(pm.processedIIndexs, :) ...
                                + diffMons;
                            pm.counts(pm.processedIIndexs(monIdxs), :) = ...
                                + pm.counts(pm.processedIIndexs(monIdxs), :) ...
                                + pcComp(monIdxs, :) * diffCpxs;
                            pm.counts(pm.matureIndexs(unprocessedMonIdxs), :) = ...
                                + pm.counts(pm.matureIndexs(unprocessedMonIdxs), :) ...
                                + pcComp(unprocessedMonIdxs, :) * diffCpxs;
                            
                            processes = sim.processes;
                            sim.setForTest('processes', {met});
                            sim.evolveState();
                            sim.setForTest('processes', processes);
                            
                            growths(j, k) = mr.growth;
                        end
                    end
                    figData.data{i} = struct;
                    figData.data{i}.time = cumsum([ 0; log(2) ./ mean(growths(1:end-1, :), 2) .* exp(log(2) * (0:j-2)'/4) / 4]) / 3600;
                    figData.data{i}.growth = mean(growths, 2) * 3600 * cellInitialDryWeight / (1 - mass.fractionWetWeight) * 1e15;
                    figData.data{i}.growths = growths * 3600 * cellInitialDryWeight / (1 - mass.fractionWetWeight) * 1e15;
                    figData.data{i}.simEndTimes = j;
                end
            end
            
            %% plot data
            clf(figHandle);
            
            options = struct;
            options.colWidths = [2.4 0.5 2.4];
            options.yAxesLabelWidths = [0.14 0.01 0.01];
            options.xlabelStr = {'Time (h)', [], 'Generation'};
            [axesHandles, xAxesHandles] = PlotUtil.multiElementPlot(figHandle, {50; 50; 50}, {[0 maxTime] [0 maxTime] generations([1 end])}, options);
            delete(xAxesHandles{2});
            
            colors = PlotUtil.getRedGreenColorOrder(1:nGenerations);
            colors = [
                colors(1, :)
                0 1 1
                colors(1:end, :)
                ];
            classLabels = [
                {'Generation 1 (Full Model)'}
                cellfun(@(i) ['Generation ' num2str(i)], num2cell((1:nGenerations)'), 'UniformOutput', false)
                {'Wild-type'}
                ];
            SingleGeneDeletions.plotSingleGeneDeletionTemporalDynamics(axesHandles{1}(1), axesHandles{2}(1), ...
                [data figData.data'], 1:nGenerations+2, classLabels, [false; true; false(size(figData.data))], 'growth', 'Growth (fg h^{-1})', colors);
            
            axesHandle = axesHandles{3}(1);
            growths = [];
            groups = [];
            for i = 1:nGenerations
                growths = [growths figData.data{i}.growths(1, :)];
                groups = [groups repmat(generations(i), 1, nSamples)];
            end
            boxplot(growths, groups, 'Parent', axesHandle);
            ylim(axesHandle, ylim(axesHandles{1}(1)));
            set(axesHandle, 'YTick', []);
        end
    end
    
    %helper methods
    methods (Static = true)
        function [deletionDoublingTimes, wtDoublingTime] = getExperimentalDoublingTimes(sim)
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            if nargin < 1
                sim = CachedSimulationObjectUtil.load();
            end
            g = sim.gene;
            
            %read data
            [~, ~, raw] = xlsread('data/singleGeneDeletions.xls', 'Doubling Time', 'A1:AA100', 'basic');
            colLabels = raw(2, :);
            strainIDs = raw(3:end, strcmp(colLabels, 'ID'));
            strainMeanDoublingTimes = cell2mat(raw(3:end, strcmp(colLabels, 'Mean (h)')));
            strainStdDoublingTimes = cell2mat(raw(3:end, strcmp(colLabels, 'Standard Deviation (h)')));
            strainNumReplicateDoublingTimes = cell2mat(raw(3:end, strcmp(colLabels, 'No. Replicates')));
            
            %wild-type growth rate
            wtDoublingTime = NaN(1, 3);
            wtDoublingTime(1, 1) = strainMeanDoublingTimes(ismember(strainIDs, 'WT'));
            wtDoublingTime(1, 2) = strainStdDoublingTimes(ismember(strainIDs, 'WT'));
            wtDoublingTime(1, 3) = strainNumReplicateDoublingTimes(ismember(strainIDs, 'WT'));
            
            %deletion growth rates
            deletionDoublingTimes = NaN(size(g.wholeCellModelIDs, 1), 3);
            
            deletionDoublingTimes(strcmp(g.essential, 'Y'), 1) = 0;
            deletionDoublingTimes(strcmp(g.essential, 'Y'), 2) = NaN;
            deletionDoublingTimes(strcmp(g.essential, 'Y'), 3) = 0;
            
            [tfs, idxs] = ismember(strainIDs, g.wholeCellModelIDs);
            deletionDoublingTimes(idxs(tfs), 1) = strainMeanDoublingTimes(tfs);
            deletionDoublingTimes(idxs(tfs), 2) = strainStdDoublingTimes(tfs);
            deletionDoublingTimes(idxs(tfs), 3) = strainNumReplicateDoublingTimes(tfs);
            deletionDoublingTimes(isnan(deletionDoublingTimes(:, 1)), 2) = wtDoublingTime(1, 2);
            deletionDoublingTimes(isnan(deletionDoublingTimes(:, 1)), 3) = 0;
            deletionDoublingTimes(isnan(deletionDoublingTimes(:, 1)), 1) = wtDoublingTime(1, 1);
        end
        
        function [tfs, molecules, processes] = getGeneImplementations(sim)
            g = sim.gene;
            
            tfs = false(size(g.wholeCellModelIDs));
            molecules = repmat({{}}, size(g.wholeCellModelIDs));
            processes = repmat({{}}, size(g.wholeCellModelIDs));
            
            for i = 1:length(sim.processes)
                m = sim.processes{i};
                
                %stimuli
                stimuliGeneComposition = m.stimuliGeneComposition();
                for j = 1:length(sim.gene.wholeCellModelIDs)
                    for k = 1:length(m.stimuliWholeCellModelIDs)
                        if ~stimuliGeneComposition(j, k)
                            continue;
                        end
                        
                        tfs(j) = true;
                        molecules{j} = [molecules{j}; m.stimuliWholeCellModelIDs(k)];
                        processes{j} = [processes{j}; m.wholeCellModelID];
                    end
                end
                
                %substrates
                substrateGeneComposition = m.substrateGeneComposition();
                for j = 1:length(sim.gene.wholeCellModelIDs)
                    for k = 1:length(m.substrateWholeCellModelIDs)
                        if ~substrateGeneComposition(j, k)
                            continue;
                        end
                        
                        tfs(j) = true;
                        molecules{j} = [molecules{j}; m.substrateWholeCellModelIDs(k)];
                        processes{j} = [processes{j}; m.wholeCellModelID];
                    end
                end
                
                %enzymes
                enzymeGeneComposition = m.enzymeGeneComposition();
                for j = 1:length(sim.gene.wholeCellModelIDs)
                    for k = 1:length(m.enzymeWholeCellModelIDs)
                        if ~enzymeGeneComposition(j, k)
                            continue;
                        end
                        
                        tfs(j) = true;
                        molecules{j} = [molecules{j}; m.enzymeWholeCellModelIDs(k)];
                        processes{j} = [processes{j}; m.wholeCellModelID];
                    end
                end
            end
        end
    end
end