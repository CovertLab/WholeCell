% ChromosomePositionHistogram
% NOTE: Histograms for RNA Polymerase
% Author: Derek Macklin, macklin@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Created: 8/18/2011
classdef ChromosomePositionHistogram
    properties (Constant)
        boundIndices = struct( ...
            'complex', {{...
                82, ...                                             % SMC
                50, ...                                             % Helicase
                [180, 181, 183, 185, 187, 189, 191, 193], ...       % DnaA
                1, ...                                              % Gyrase
                78, ...                                             % Topo IV
                165, ...                                            % LuxR
                }}, ...
            'monomer', {{
                101, ...                                            % HTH regulator
                242, ...                                            % Ferric uptake repressor
                }}...
            );
        boundIndicesNames = struct( ...
            'complex', {{...
                'SMC', ...
                'Helicase', ...
                'DnaA', ...
                'Gyrase', ...
                'Topo IV', ...
                'LuxR'
                }}, ...
            'monomer', {{
                'HTH regulator', ...
                'Ferric uptake\nrepressor', ...
                }}...
            );
    end
    
    methods (Static = true)
        function run(simDir, fileName)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.analysis.ChromosomePositionHistogram;
            
            %% options
            if nargin < 1 || isempty(simDir)
                simDir = SimulationDiskUtil.getLatestSimulation();
            end
            
            [simDir, ~, sim] = SimulationDiskUtil.getSimulation(simDir);
            c = sim.state('Chromosome');
            
            %% get data
            stateNames = {
                'ProteinComplex' 'counts'
                'Chromosome' 'complexBoundSites'
                'Chromosome' 'monomerBoundSites'
                'RNAPolymerase' 'positionStrands'
                'RNAPolymerase' 'states'
                'Time' 'values'};
            states = DiskLogger.load(simDir, stateNames, [], [], [], 'extract');
            
            cDensityMatrix = ChromosomePositionHistogram.makeDensityMatrix(states, ChromosomePositionHistogram.boundIndices, sim);
            
            %% figure 1
            [~, figHandle1] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            clf(figHandle1);
            nPlots = numel(ChromosomePositionHistogram.boundIndices.complex);
            [axesHandles, xAxesHandle] = PlotUtil.multiElementPlot(...
                figHandle1, 3.5 * ones(nPlots, 1), [1 c.sequenceLen], ...
                struct(...
                'titleStr', 'Protein Density on Chromosome 1', ...
                'xlabelStr', 'Position on Chromosome 1'));
            set(xAxesHandle, 'XTick', [1 c.terCPosition c.sequenceLen], 'XTickLabel', {'oriC', 'terC', 'oriC'});
            
            for i = 1:size(cDensityMatrix, 2)
                PlotUtil.plotLine(axesHandles(i), 1:size(cDensityMatrix, 1), cDensityMatrix(:, i), false, true, false);
                ylabel(axesHandles(i), {ChromosomePositionHistogram.boundIndicesNames.complex{i} '(Counts)'});
            end
            if exist('fileName', 'var') && ~isempty(fileName)
                saveas(figHandle1, sprintf('%s-1.pdf', fileName));
                close(figHandle1);
            end
            
            %% figure 2
            [~, figHandle2] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
            ChromosomePositionHistogram.plotRnaPolymeraseDensity(sim, states, figHandle2);
            
            if exist('fileName', 'var') && ~isempty(fileName)
                saveas(figHandle2, sprintf('%s-2.pdf', fileName));
                close(figHandle2);
            end
        end
    end
    
    %multiple simulations
    methods (Static = true)
        
        % TODO: Handle cases where we don't have one of monomers or complexes
        function [cDensityMatrix, mDensityMatrix] = makeDensityMatrix(states, boundIndices, simulation)
            complexIdxs = boundIndices.complex;
            monomerIdxs = boundIndices.monomer;
            
            cbs = states.Chromosome.complexBoundSites;
            mbs = states.Chromosome.monomerBoundSites;
            
            ch = simulation.state('Chromosome');
            
            cDensityMatrix = zeros(size(cbs, 1), numel(complexIdxs));
            for i = 1:numel(complexIdxs)
                ftpt = ones(ch.complexDNAFootprints(complexIdxs{i}(1)), 1);
                for j = 1:numel(complexIdxs{i})
                    complexIdx = complexIdxs{i}(j);
                    cDensityMatrix(:, i) = cDensityMatrix(:, i) + ...
                        cconv(...
                        full(sum(sum(cbs(:, 1:2, :) == complexIdx, 3), 2)), ...
                        ftpt, ...
                        size(cDensityMatrix, 1) ...
                        );
                end
            end
            cDensityMatrix = cDensityMatrix / size(cbs, 3);
            
            mDensityMatrix = zeros(size(mbs, 1), numel(monomerIdxs));
            for i = 1:numel(monomerIdxs)
                ftpt = ones(ch.monomerDNAFootprints(monomerIdxs{i}(1)), 1);
                for j = 1:numel(monomerIdxs{i})
                    monomerIdx = monomerIdxs{i}(j);
                    mDensityMatrix(:, i) = mDensityMatrix(:, i) + ...
                        cconv(...
                        full(sum(sum(mbs(:, 1:2, :) == monomerIdx, 3), 2)), ...
                        ftpt, ...
                        size(mDensityMatrix, 1) ...
                        );
                end
            end
            mDensityMatrix = mDensityMatrix / size(mbs, 3);
            
            % Hack to deal with rounding issues that lead to small
            % negative numbers
            cDensityMatrix = max(0, cDensityMatrix);
            mDensityMatrix = max(0, mDensityMatrix);
        end
        
        function plotRnaPolymeraseDensity(sim, states, figHandle)
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.util.CircularSparseMat;
            
            c = sim.state('Chromosome');
            pol = sim.state('RNAPolymerase');
            pc = sim.state('ProteinComplex');
            
            nTimePoints = size(states.RNAPolymerase.positionStrands, 3);
            nRnaPol = size(states.RNAPolymerase.positionStrands, 1);
            polStates = reshape(states.RNAPolymerase.states, [nRnaPol * nTimePoints 1]);
            posStrndTimes = reshape(permute([states.RNAPolymerase.positionStrands repmat(permute(1:nTimePoints, [1 3 2]), [nRnaPol 1])], [1 3 2]), [nRnaPol * nTimePoints 3]);
            
            polStates = polStates(posStrndTimes(:, 1) ~= 0, :);
            posStrndTimes = posStrndTimes(posStrndTimes(:, 1) ~= 0, :);
            
            % Get rid of strands 3 and 4 (i.e., only look at Chromosome 1)
            polStates = polStates(posStrndTimes(:, 2) <= 2, :);
            posStrndTimes = posStrndTimes(posStrndTimes(:, 2) <= 2, :);
            
            boundSites = CircularSparseMat(posStrndTimes, polStates, [c.sequenceLen 4 nTimePoints], 1);
            
            activeDensity = cconv(...
                full(sum(sum(boundSites >= pol.activelyTranscribingValue, 3), 2)), ...
                ones(c.complexDNAFootprints(pc.getIndexs('RNA_POLYMERASE')), 1), ...
                c.sequenceLen) / nTimePoints;
            nonSpecBoundDensity = cconv(...
                full(sum(sum(boundSites == pol.nonSpecificallyBoundValue, 3), 2)), ...
                ones(c.complexDNAFootprints(pc.getIndexs('RNA_POLYMERASE')), 1), ...
                c.sequenceLen) / nTimePoints;
            specBoundDensity = cconv(...
                full(sum(sum(boundSites == pol.specificallyBoundValue, 3), 2)), ...
                ones(c.complexDNAFootprints(pc.getIndexs('RNA_POLYMERASE')), 1), ...
                c.sequenceLen) / nTimePoints;
            
            % Hack to deal with rounding issues that lead to small
            % negative numbers
            activeDensity = max(0, activeDensity);
            nonSpecBoundDensity = max(0, nonSpecBoundDensity);
            specBoundDensity = max(0, specBoundDensity);
            
            %layout axes
            clf(figHandle);
            [axesHandles, xAxesHandle] = PlotUtil.multiElementPlot(figHandle, 8 * ones(3, 1), [1 c.sequenceLen], struct( ...
                'titleStr', 'RNA Polymerase Density on Chromosome 1', ...
                'xlabelStr', 'Position on Chromosome 1'));
            set(xAxesHandle, 'XTick', [1 c.terCPosition c.sequenceLen], 'XTickLabel', {'oriC', 'terC', 'oriC'});
            
            %plot data
            PlotUtil.plotLine(axesHandles(1), 1:c.sequenceLen, activeDensity, false, true, false);
            ylabel(axesHandles(1), 'Active');
            
            PlotUtil.plotLine(axesHandles(2), 1:c.sequenceLen, specBoundDensity, false, true, false);
            ylabel(axesHandles(2), {'Specifically' 'Bound'});
            
            PlotUtil.plotLine(axesHandles(3), 1:c.sequenceLen, nonSpecBoundDensity, false, true, false);
            ylabel(axesHandles(3), {'Non-' 'Specifically' 'Bound'});
            
            % Format Y-axis
            PlotUtil.alignYAxesLabels(axesHandles);
            PlotUtil.offsetYAxes(axesHandles, 0.015);
        end
    end
end
