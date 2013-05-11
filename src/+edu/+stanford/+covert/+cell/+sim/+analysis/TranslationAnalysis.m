% TranslationAnalysis
%
% Author: Derek Macklin, macklin@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Created: 8/31/2011
classdef TranslationAnalysis
    methods (Static = true)
        function run(simDir, fileName)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.analysis.ChromosomePositionHistogram;
            import edu.stanford.covert.cell.sim.analysis.TranslationAnalysis;
            
            %% options
            if nargin < 1 || isempty(simDir)
                simDir = [SimulationDiskUtil.getLatestSimulationGroup() filesep '1'];
            end
            
            %% constants
            [simDir, ~, sim] = SimulationDiskUtil.getSimulation(simDir);
            simBatchDir = SimulationDiskUtil.getSimulationBatchDir(simDir);
            simIdx = num2str(SimulationDiskUtil.getSimulationIndex(simDir));
            simTimeStamp = SimulationDiskUtil.getSimulationTimeStamp(simDir);
            
            g = sim.gene;
            met = sim.state('Metabolite');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            pp = sim.state('Polypeptide');
            
            atpIdx = met.getIndexs('ATP');
            gtpIdx = met.getIndexs('GTP');
            
            efGIdx = pc.boundIndexs(pc.getIndexs('MG_089_DIMER'));
            efTsIdx = pc.boundIndexs(pc.getIndexs('MG_433_DIMER'));
            efTuIdx = pc.boundIndexs(pc.getIndexs('MG_451_DIMER'));
            
            efPIdx = pm.boundIndexs(pm.getIndexs('MG_026_MONOMER'));
            
            %% get data
            stateNames = {...
                'Metabolite'        'counts'         [atpIdx; gtpIdx; met.aminoAcidIndexs] 1
                'ProteinComplex'    'counts'         [efGIdx; efTsIdx; efTuIdx] 1
                'ProteinMonomer'    'counts'         efPIdx 1
                'Ribosome'          'boundMRNAs'	 ':' ':'
                'Ribosome'          'mRNAPositions'  ':' ':'
                'Ribosome'          'tmRNAPositions' ':' ':'
                'Ribosome'          'states'         ':' ':'
                'Rna'               'counts'         [rna.matureIndexs(rna.matureTRNAIndexs); rna.aminoacylatedIndexs(rna.matureTRNAIndexs)] 1
                'Polypeptide'       'boundMRNAs'	 ':' ':'
                'Time'              'values'         ':' ':'
                };
            
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], [], 'extract', str2double(simIdx));
            
            time = squeeze(states.Time.values) / 3600;
            
            aminoacylatedtRNACounts = permute(states.Rna.counts(1+numel(rna.matureIndexs(rna.matureTRNAIndexs)):end, 1, :), [1 3 2]);
            atpCounts = permute(states.Metabolite.counts(1, 1, :), [1 3 2]);
            gtpCounts = permute(states.Metabolite.counts(2, 1, :), [1 3 2]);
            aminoAcidCounts = permute(states.Metabolite.counts(3:end, 1, :), [1 3 2]);
            efGCounts = permute(states.ProteinComplex.counts(1, 1, :), [1 3 2]);
            
            ribStates = permute(states.Ribosome.states, [1 3 2]);
            pos = permute(states.Ribosome.mRNAPositions, [1 3 2]);
            stallPos = permute(states.Ribosome.tmRNAPositions, [1 3 2]);
            
            [starts, ends, pauses] = TranslationAnalysis.getStartsEndsPauses(sim, pos, stallPos, ribStates);            
            [startRibs, endRibs, startTimes, endTimes] = TranslationAnalysis.verifyStartsEnds(starts, ends);
           
            %% synthesis
            [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            TranslationAnalysis.plotMonomerSynthesis(axesHandle, states, sim, starts, ends);
            
            if exist('fileName', 'var')
                saveas(figHandle, [fileName '-MonomerSynthesis.pdf'])
                close(figHandle);
            end
            
            %% pauses
            pausesInfo = TranslationAnalysis.getPausesInfo(states, pauses, false);
            [~, tmp] = max(pausesInfo(:, 5));
            maxPausedRibosomeIdx = pausesInfo(tmp, 1);
            
            ribLP = pausesInfo(1, 1);
            bmRNALP = pausesInfo(1, 2);
            positionLP = pausesInfo(1, 3);
            startLP = pausesInfo(1, 4);
            timeLP = pausesInfo(1, 5);
            endLP = startLP + timeLP - 1;
            
            startInfo = [startRibs, startTimes];
            endInfo = [endRibs, endTimes];
            translateStartTime = startInfo(startInfo(:, 1) == ribLP & startInfo(:, 2) <= startLP, :);
            translateEndTime = endInfo(startInfo(:, 1) == maxPausedRibosomeIdx & endInfo(:, 2) >= endLP, :);
            translateStartTime = translateStartTime(end, 2);
            translateEndTime = translateEndTime(1, 2);
            
            [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            clf(figHandle);
            
            nRows = 7;
            
            posL = 0.10;
            posB = 0.10;
            posW = 0.80;
            posH = 0.75;
            
            timeSliceIdxs = translateStartTime:translateEndTime;
            figTitle = sprintf('Longest Ribosome Stall\nSimulation: %s #%s', simTimeStamp, simIdx);
            annotation(figHandle, 'TextBox', [0.25 0.925 0.5 0.025], ...
                'String', figTitle, 'HorizontalAlignment', 'Center', ...
                'VerticalAlignment', 'Middle', 'EdgeColor', 'None', ...
                'FontSize', 14, 'FontWeight', 'Normal');
            axesHandles = PlotUtil.multiElementPlot(...
                figHandle, 3.5 * ones(nRows, 1), time([translateStartTime translateEndTime]), struct(...
                'position', [posL posB posW posH] ...
                ));
                    
            geneWID = g.wholeCellModelIDs(g.mRNAIndexs(bmRNALP));
            PlotUtil.plotLine(axesHandles(1), time(timeSliceIdxs), pos(ribLP, timeSliceIdxs), true, true, false);
            title(axesHandles(1), sprintf('Ribosomal Progress on %s', geneWID{1}), 'Interpreter', 'none');
            ylabel(axesHandles(1), 'AAs Inc.');
            
            pausedTRNAIdx = pp.monomerTRNASequences{bmRNALP}(positionLP + 1);
            tRNAWIDs = rna.wholeCellModelIDs(rna.matureIndexs(rna.matureTRNAIndexs));
            PlotUtil.plotLine(axesHandles(2), time(timeSliceIdxs), aminoacylatedtRNACounts(pausedTRNAIdx, timeSliceIdxs), true, true, false);
            title(axesHandles(2), sprintf('Aminoacylated tRNA %s', tRNAWIDs{pausedTRNAIdx}), 'Interpreter', 'none');
            ylabel(axesHandles(2), 'Counts');
            
            aminoAcidAbbr = pp.monomerAASequences{bmRNALP}(positionLP + 1);
            [~, aminoAcidRelIdx] = ismember(aminoAcidAbbr, edu.stanford.covert.cell.kb.ProteinMonomer.bases);
            aminoAcidWID = met.wholeCellModelIDs(met.aminoAcidIndexs(aminoAcidRelIdx));
            PlotUtil.plotLine(axesHandles(3), time(timeSliceIdxs), aminoAcidCounts(aminoAcidRelIdx, timeSliceIdxs), true, true, false);
            title(axesHandles(3), sprintf('Amino Acid %s', aminoAcidWID{1}), 'Interpreter', 'none');
            ylabel(axesHandles(3), 'Counts');
            
            PlotUtil.plotLine(axesHandles(4), time(timeSliceIdxs), atpCounts(timeSliceIdxs), true, true, false);
            title(axesHandles(4), 'ATP', 'Interpreter', 'none');
            ylabel(axesHandles(4), 'Counts');
            
            PlotUtil.plotLine(axesHandles(5), time(timeSliceIdxs), gtpCounts(timeSliceIdxs), true, true, false);
            title(axesHandles(5), 'GTP', 'Interpreter', 'none');
            ylabel(axesHandles(5), 'Counts');
            
            PlotUtil.plotLine(axesHandles(6), time(timeSliceIdxs), efGCounts(timeSliceIdxs), true, true, false);
            title(axesHandles(6), 'Bound Elongation Factors', 'Interpreter', 'none');
            ylabel(axesHandles(6), 'Counts');
            
            PlotUtil.plotLine(axesHandles(7), time(timeSliceIdxs), sum(ribStates(:, timeSliceIdxs), 1), true, true, false);
            title(axesHandles(7), 'Number of Active Ribosomes', 'Interpreter', 'none');
            ylabel(axesHandles(7), 'Counts');
            
            if exist('fileName', 'var') && ~isempty(fileName)
                saveas(figHandle, [fileName '-Pauses.pdf']);
                close(figHandle);
            end
            
            %% histogram of pause lengths
            [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            
            hist(axesHandle, pausesInfo(:, 5));
            title(axesHandle, 'Duration Ribosome Pauses', 'FontSize', 12)
            xlabel(axesHandle, 'Pause Duration (s)', 'FontSize', 10);
            ylabel(axesHandle, 'Frequency', 'FontSize', 10);
            
            if exist('fileName', 'var') && ~isempty(fileName)
                saveas(figHandle, [fileName '-PauseDistribution.pdf']);
                close(figHandle);
            end
        end
    end
    
    methods (Static = true)
        function [starts, ends, pauses] = getStartsEndsPauses(sim, pos, stallPos, ribStates)
            rib = sim.state('Ribosome');
                
            dpos = diff(pos, 1, 2);
            dStallPos = diff(stallPos, 1, 2);
            
            ends = false(size(pos));
            ends(:, 1:end-1) = ...
                ((dpos < 0 | dStallPos < 0) & ribStates(:, 1:end-1) ~= rib.notExistValue) | ...
                (ribStates(:, 1:end-1) ~= rib.notExistValue & ribStates(:, 2:end) == rib.notExistValue);
            ends(:, end) = ribStates(:, end) ~= rib.notExistValue;
            
            starts = false(size(pos));
            starts(:, 2:end) = ...
                ((dpos < 0 | dStallPos < 0) & ribStates(:, 2:end) ~= rib.notExistValue) | ...
                (ribStates(:, 1:end-1) == rib.notExistValue & ribStates(:, 2:end) ~= rib.notExistValue);
            starts(:, 1) = ribStates(:, 1) ~= rib.notExistValue;
            
            pauses = false(size(pos));
            pauses(dpos == 0 & ribStates(:, 1:(end - 1)) == rib.activeValue) = true;
        end
        
        function [startRibs, endRibs, startTimes, endTimes] = verifyStartsEnds(starts, ends)
            [startTimes, startRibs] = find(starts');
            [endTimes, endRibs] = find(ends');

            assert(isequal(startRibs, endRibs));
            assert(~any(endTimes < startTimes));
        end
        
        function plotMonomerSynthesis(axesHandle, states, sim, starts, ~)
            import edu.stanford.covert.cell.sim.util.PlotUtil;

            g = sim.gene;
            pm = sim.state('ProteinMonomer');
            rna = sim.state('Rna');
            
            monSynthesis = histc(states.Ribosome.boundMRNAs(starts), (1:numel(pm.matureIndexs)));
            monSynthesisExpected = rna.matureRNAGeneComposition(g.mRNAIndexs, :) * rna.expression(rna.matureIndexs);

            plot(axesHandle, monSynthesisExpected, monSynthesis, '.')
            xlabel(axesHandle, 'Expected');
            ylabel(axesHandle, 'Simulated');            
        end
        
        function pausesInfo = getPausesInfo(states, pauses, verbose)
            bmRNA = states.Ribosome.boundMRNAs;
            pos = states.Ribosome.mRNAPositions;
            pausesTFs = false(size(pos, 1), 1);
            time = squeeze(states.Time.values);
            
            for i = 1:size(pos, 1)
                if ~isempty(pos(pauses(i, :, :)))
                    pausesTFs(i) = true;
                    if verbose
                        fprintf('%d\n', i);
                    end
                end
            end
            
            pausesInds = 1:size(pos, 1);
            pausesInds = pausesInds(pausesTFs);
            
            pausesInfo = zeros(sum(pauses(:)), 5);
            
            j = 1;
            for i = 1:numel(pausesInds)
                elemThisIter = sum(squeeze(pauses(pausesInds(i), :, :)));
                pausesInfo(j:(j + elemThisIter - 1), 1:4) = [...
                    pausesInds(i) * ones(elemThisIter, 1), ...
                    squeeze(bmRNA(pausesInds(i), :, pauses(pausesInds(i), :, :))), ...
                    squeeze(pos(pausesInds(i), :, pauses(pausesInds(i), :, :))), ...
                    time(squeeze(pauses(pausesInds(i), :, :)))...
                    ];
                j = j + elemThisIter;
            end
            
            i = 1;
            while i < size(pausesInfo, 1)
                j = 1;
                while ((i + j) < size(pausesInfo, 1)) && isequal(pausesInfo(i+j, 1:3), pausesInfo(i, 1:3))
                    j = j + 1;
                end
                pausesInfo(i, 5) = j;
                i = i + j + 1;
            end
            pausesInfo(end, 5) = 1;
            
            pausesInfo(pausesInfo(:, 5) == 0, :) = [];
            
            pausesInfo = flipud(sortrows(pausesInfo, 5));
        end
    end
end
