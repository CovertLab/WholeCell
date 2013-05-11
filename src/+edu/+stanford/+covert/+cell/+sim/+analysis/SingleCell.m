%SingleCell
%
% Author: Derek Macklin, macklin@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 7/6/2011
classdef SingleCell
    methods (Static = true)
        function run(simBatchDir, simNum, fileName, methodsToRun)
            import edu.stanford.covert.cell.sim.analysis.SingleCell;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            %options
            if nargin < 1 || isempty(simBatchDir)
                simBatchDir = SimulationDiskUtil.getLatestSimulationGroup();
            end
            
            if nargin < 2 || isempty(simNum)
                simNum = 1;
            end
            
            %simulation
            [~, ~, sim] = SimulationDiskUtil.getSimulation([simBatchDir filesep num2str(simNum)]);
            
            %evaluate each method and optionally save at the result as
            %pdf files
            metaData = meta.class.fromName('edu.stanford.covert.cell.sim.analysis.SingleCell');
            for i = 1:numel(metaData.Methods)
                %only plot methods matching the signature
                %  plotFcn(sim, simBatchDir, simNum, figHandle)
                if ...
                        numel(metaData.Methods{i}.InputNames) < 4 || ~isequal(metaData.Methods{i}.InputNames(1:4), {'sim'; 'simBatchDir'; 'simNum'; 'figHandle'}) || ...
                        numel(metaData.Methods{i}.OutputNames) < 1 || ~isequal(metaData.Methods{i}.OutputNames(1), {'figData'})
                    continue;
                end
                
                %get method name
                methodName = metaData.Methods{i}.Name;
                MethodName = [upper(methodName(1)) methodName(2:end)];
                if nargin >= 4 && ~ismember(metaData.Methods{i}.Name, methodsToRun)
                    continue;
                end
                
                %run plot method
                [~, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle();
                try
                    argout = cell(size(metaData.Methods{i}.OutputNames));
                    [argout{:}] = SingleCell.(methodName)(sim, simBatchDir, simNum, figHandle);
                    
                    % save plot
                    if nargin >= 3
                        figData = argout{1}; %#ok<NASGU>
                        save([fileName '-' MethodName '.mat'], '-struct', 'figData');
                        
                        orient(figHandle, 'portrait');
                        saveas(figHandle, [fileName '-' MethodName '.pdf']);
                    end
                    
                    % save table
                    if numel(argout) > 1
                        if nargin >= 3
                            PrintUtil.printToFile(argout{2}, argout{3}, [fileName '-' MethodName '.xls'], MethodName);
                        else
                            PrintUtil.printToStdIO(argout{2}, argout{3});
                        end
                    end
                catch exception
                    warning('WholeCell:warning', 'Unable to make %s plot:\n%s', methodName, exception.getReport());
                end
                
                if nargin >= 2
                    close(figHandle);
                    clear figHandle;
                end
                clear argout figData;
            end
        end
        
        function figData = growth(sim, simBatchDir, simNum, figHandle)
            import edu.stanford.covert.cell.sim.analysis.ChromosomeSpaceTimePlot;
            import edu.stanford.covert.cell.sim.analysis.SingleCell;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            %% get data
            [time, mass, rnas, proteins, ...
                dntps, ntps, aminoAcids, ...
                growth, reactionFlux, ...
                enzymeComplex, enzymeMonomer, enzymeRNA, enzymeGene, ...
                ploidy, superhelicalDensity] = ...
                SingleCell.calcGrowthData(sim, simBatchDir, simNum);
            
            %% layout figure, plot data
            clf(figHandle);
            
            options = struct();
            options.titleStr = {'Composition'; 'Growth'; 'Organization'};
            options.position = [0.12 0.52 0.87 0.4];
            options.xdata = time;
            options.ydata = {{
                permute(mass, [1 3 2 4])
                permute(rnas, [1 3 2 4])
                permute(proteins, [1 3 2 4])
                permute(dntps, [1 3 2 4])
                permute(ntps, [1 3 2 4])
                permute(aminoAcids, [1 3 2 4])
                }
                {
                permute(growth, [1 3 2 4])
                permute(reactionFlux, [1 3 2 4])
                permute(enzymeComplex, [1 3 2 4])
                permute(enzymeMonomer, [1 3 2 4])
                permute(enzymeRNA, [1 3 2 4])
                permute(enzymeGene, [1 3 2 4])
                }
                {
                permute(ploidy, [1 3 2 4])
                permute(superhelicalDensity, [1 3 2 4])
                []
                }};
            options.ylabelStr = {{
                'Mass (fg)'
                'RNAs'
                'Proteins'
                'dNTPs'
                'NTPs'
                {'Amino' 'acids'}
                }
                {
                {'Growth' '(fg h^{-1})'}
                {sprintf('%s', rxnId) 'flux' '(MHz)'}
                {sprintf('%s', rxnId) 'complex'}
                {sprintf('%s', rxnId) 'monomer'}
                {sprintf('%s', rxnId) 'mRNA'}
                sprintf('{\\it %s}', rxnId)
                }
                {
                {'Copy', 'number'}
                {'Super-' 'helicity'}
                {'Polymerase' 'position'}
                }
                };
            axesHandles = PlotUtil.multiElementPlot(figHandle, {3 * ones(6, 1); 3 * ones(6, 1); [3 3 4*3+(4-1)*1]}, time([1 end])', ...
                options);
            ChromosomeSpaceTimePlot.plotSpaceTime(axesHandles{3}(3), simBatchDir, simNum, [], false);
            
            clear options;
            
            %% Format Y-axis
            PlotUtil.labelSubplots([axesHandles{1}(1); axesHandles{2}(1); axesHandles{3}(1)]);
            
            %% collect data
            figData = struct;
            figData.time = time;
            figData.mass = mass;
            figData.rnas = rnas;
            figData.proteins = proteins;
            figData.dntps = dntps;
            figData.ntps = ntps;
            figData.aminoAcids = aminoAcids;
            figData.growth = growth;
            figData.reactionFlux = reactionFlux;
            figData.enzymeComplex = enzymeComplex;
            figData.enzymeMonomer = enzymeMonomer;
            figData.enzymeRNA = enzymeRNA;
            figData.enzymeGene = enzymeGene;
            figData.ploidy = ploidy;
            figData.superhelicalDensity = superhelicalDensity;
        end
        
        function [time, mass, rnas, proteins, ...
                dntps, ntps, aminoAcids, ...
                growth, reactionFlux, ...
                enzymeComplex, enzymeMonomer, enzymeRNA, enzymeGene, ...
                ploidy, superhelicalDensity] = ...
                calcGrowthData(sim, simBatchDir, simNum, rxnId)
            
            import edu.stanford.covert.cell.sim.util.MacromoleculeUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.analysis.CellOverview;
            import edu.stanford.covert.util.ConstantUtil;
            
            comp = sim.compartment;
            g = sim.gene;
            
            c = sim.state('Chromosome');
            m = sim.state('Metabolite');
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            sc = sim.process('DNASupercoiling');
            met = sim.process('Metabolism');
            
            %% get data
            stateNames = {
                'mass'              'Mass'
                'ploidy'            {'Chromosome' 'Copy Number'}
                'rnas'              'RNA'
                'proteins'          'Proteins'
                'lipids'            'Lipids'
                'growth_rate'       'Growth'
                'amino_acids'       {'Amino' 'Acids'}
                'ntps'              'NTPs'
                };
            ensemble = SimulationEnsemble(simBatchDir, stateNames, [], simNum);
            
            time = ensemble.stateData.time / 3600;
            mass = ensemble.stateData.values(ensemble.getPropertyIndices('mass'), :, :, :) * sim.state('Mass').cellInitialDryWeight / (1-sim.state('Mass').fractionWetWeight) * 1e15;
            ploidy = ensemble.stateData.values(ensemble.getPropertyIndices('ploidy'), :, :, :);
            rnas = ensemble.stateData.values(ensemble.getPropertyIndices('rnas'), :, :, :);
            proteins = ensemble.stateData.values(ensemble.getPropertyIndices('proteins'), :, :, :);
            %lipids = ensemble.stateData.values(ensemble.getPropertyIndices('lipids'), :, :, :);
            growth = ensemble.stateData.values(ensemble.getPropertyIndices('growth_rate'), :, :, :) * 3600 * sim.state('Mass').cellInitialDryWeight / (1-sim.state('Mass').fractionWetWeight)  * 1e15;
            
            clear ensemble;
            
            if nargin < 4
                rxnId = 'AckA';
            end
            rxnIdx = met.reactionIndexs(rxnId);
            enzId = met.enzymeWholeCellModelIDs{met.reactionCatalysisMatrix(rxnIdx, :) ~= 0};
            [~, ~, ~, ~, ~, ...
                geneIdxs, ~, ~, ~, ...
                ~, ~, ~, ~, ...
                matureRnaIdxs, ~, ~, ~, ...
                monomerIdxs, ~, ~, ~, ...
                complexIdxs, ~, ~, ~] = ...
                MacromoleculeUtil.getMacromoleculeIndexsIDsNames(enzId, sim);
            stateNames = {
                'Metabolite'        'counts'              ':'  ':'
                'MetabolicReaction' 'fluxs'               rxnIdx ':'
                'Chromosome'        'polymerizedRegions'  ':'  ':'
                'Chromosome'        'linkingNumbers'       ':'  ':'
                'Rna'               'counts' rna.matureIndexs(matureRnaIdxs)         comp.cytosolIndexs
                'ProteinMonomer'    'counts' [pm.matureIndexs([monomerIdxs; sc.enzymeGlobalIndexs(sc.enzymeIndexs_topoI)]); pm.boundIndexs([monomerIdxs; sc.enzymeGlobalIndexs(sc.enzymeIndexs_topoI)])] comp.cytosolIndexs
                'ProteinComplex'    'counts' [pc.matureIndexs([complexIdxs; sc.enzymeGlobalIndexs([sc.enzymeIndexs_gyrase; sc.enzymeIndexs_topoIV])]); pc.boundIndexs([complexIdxs; sc.enzymeGlobalIndexs([sc.enzymeIndexs_gyrase; sc.enzymeIndexs_topoIV])])] comp.cytosolIndexs
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', simNum);
            tmp = full(states.Metabolite.counts);
            %smallMolecules = sum(sum(tmp(setdiff(1:end, m.lipidIndexs), [comp.cytosolIndexs; comp.membraneIndexs], :), 2), 1);
            dntps = cat(3, NaN(4, 1), tmp(m.dntpIndexs, 1, :));
            ntps = cat(3, NaN(4, 1), tmp(m.ntpIndexs, 1, :));
            %nmps = tmp(m.nmpIndexs(1:20), comp.cytosolIndexs, :);
            aminoAcids = cat(3, NaN(20, 1), tmp(m.aminoAcidIndexs(1:20), 1, :));
                        
            enzymeGene = cat(3, NaN(size(geneIdxs)), max(1, permute(MacromoleculeUtil.extractCopyNumberTimeCourse(c, states.Chromosome.polymerizedRegions, ...
                g.startCoordinates(geneIdxs), g.strands(geneIdxs), g.lengths(geneIdxs)), [1 3 2 4])));
            enzymeRNA = cat(3, NaN(size(matureRnaIdxs)), states.Rna.counts);
            enzymeMonomer = cat(3, NaN(size(monomerIdxs)), states.ProteinMonomer.counts(1:numel(monomerIdxs), 1, :, :));
            enzymeComplex = cat(3, NaN(size(complexIdxs)), states.ProteinComplex.counts(1:numel(complexIdxs), 1, :, :));
            reactionFlux = cat(3, NaN, states.MetabolicReaction.fluxs * 3600 * 1e-6);
            
            lk0 = full(sum(sum(states.Chromosome.polymerizedRegions, 2), 1)) / c.relaxedBasesPerTurn;
            superhelicalDensity = cat(3, NaN, (full(sum(sum(states.Chromosome.linkingNumbers, 2), 1)) - lk0) ./ lk0);
            %gyrase = ...
            %    + states.ProteinComplex.counts(end/2-1, 1, :, :) ...
            %    + states.ProteinComplex.counts(end-1, 1, :, :);
            %topoI = ...
            %    + states.ProteinMonomer.counts(end/2, 1, :, :) ...
            %    + states.ProteinMonomer.counts(end, 1, :, :);
            %topoIV = ...
            %    + states.ProteinComplex.counts(end/2, 1, :, :) ...
            %    + states.ProteinComplex.counts(end, 1, :, :);
            
            clear tmp states lk0;
            
            %% collect data
            figData = struct;
            figData.time = time;
            figData.mass = mass;
            figData.rnas = rnas;
            figData.proteins = proteins;
            figData.dntps = dntps;
            figData.ntps = ntps;
            figData.aminoAcids = aminoAcids;
            figData.growth = growth;
            figData.reactionFlux = reactionFlux;
            figData.enzymeComplex = enzymeComplex;
            figData.enzymeMonomer = enzymeMonomer;
            figData.enzymeRNA = enzymeRNA;
            figData.enzymeGene = enzymeGene;
            figData.ploidy = ploidy;
            figData.superhelicalDensity = superhelicalDensity;
        end
        
        function figData = mass(sim, simBatchDir, simNum, figHandle, normalized)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            comp = sim.compartment;
            
            %% get data
            stateNames = {
                'Time' 'values'       ':' ':'
                'Mass' 'waterWt'      ':' comp.cellularIndexs(:)'
                'Mass' 'dnaWt'        ':' '-sum'
                'Mass' 'rnaWt'        ':' '-sum'
                'Mass' 'proteinWt'    ':' '-sum'
                'Mass' 'metaboliteWt' ':' [comp.cytosolIndexs; comp.membraneIndexs]
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', simNum);
            
            %% layout figure, plot data
            clf(figHandle);
            
            options = struct();
            options.titleStr = 'Cell Mass';
                
            options.plotFunc = @PlotUtil.plotStackedArea;
            axesHandle = PlotUtil.multiElementPlot(figHandle, 50, states.Time.values([1 end]), options);
            
            ydata = [
                permute(states.Mass.waterWt, [1 3 2])
                permute(states.Mass.proteinWt, [1 3 2])
                permute(states.Mass.dnaWt, [1 3 2])
                permute(states.Mass.metaboliteWt(:, 2, :), [1 3 2])
                permute(states.Mass.metaboliteWt(:, 1, :), [1 3 2])
                permute(states.Mass.rnaWt, [1 3 2])
                ] * 1e15;
            if nargin < 5 || ~normalized
                h = PlotUtil.plotStackedArea(axesHandle, permute(states.Time.values, [1 3 2]), ydata);
                ylabel(axesHandle, 'Mass (fg)');
                ylim(axesHandle, [0 max(sum(ydata, 1))]);
            else
                h = PlotUtil.plotNormalizedStackedArea(axesHandle, permute(states.Time.values, [1 3 2]), ydata);
                ylabel(axesHandle, 'Percent Mass');
                ylim(axesHandle, [0 100]);
            end
            legend(h, {'Water', 'Protein', 'DNA', 'Membrane', 'Metabolite', 'RNA'}, 'Location', 'SouthEast');
            
            PlotUtil.offsetYAxes(axesHandle, 0.02);
            
            clear ydata;
            
            %% collect data
            figData = struct;
            figData.time = states.Time.values;
            figData.waterWt = states.Mass.waterWt;
            figData.proteinWt = states.Mass.proteinWt;
            figData.dnaWt = states.Mass.dnaWt;
            figData.metaboliteWt = states.Mass.metaboliteWt;
            figData.rnaWt = states.Mass.rnaWt;
        end
        
        function figData = massNormalized(sim, simBatchDir, simNum, figHandle)
            import edu.stanford.covert.cell.sim.analysis.SingleCell;
            
            figData = SingleCell.mass(sim, simBatchDir, simNum, figHandle, true);
        end
        
        function figData = metabolites(sim, simBatchDir, simNum, figHandle, normalized)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            import edu.stanford.covert.util.ConstantUtil;
            
            comp = sim.compartment;
            m = sim.state('Metabolite');
            metIdxs = [
                m.ntpIndexs
                m.ndpIndexs
                m.nmpIndexs
                m.dntpIndexs
                m.aminoAcidIndexs
                m.phosphateIndexs
                m.diphosphateIndexs
                m.hydrogenIndexs
                ];
            
            %% get data
            stateNames = {
                'Time'       'values'  ':'      ':'
                'Geometry'   'volume'  ':'      ':'
                'Metabolite' 'counts'  metIdxs  comp.cytosolIndexs
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', simNum);
            time = permute(states.Time.values, [1 3 2]) / 3600;
            
            states.Metabolite.counts = ...
                states.Metabolite.counts ./ ...
                states.Geometry.volume(ones(size(metIdxs)), :, :) / ...
                ConstantUtil.nAvogadro * 1000;
            
            ntps = permute(states.Metabolite.counts(1:4, :, :), [1 3 2]);
            ndps = permute(states.Metabolite.counts(5:8, :, :), [1 3 2]);
            nmps = permute(states.Metabolite.counts(9:12, :, :), [1 3 2]);
            dntps = permute(states.Metabolite.counts(13:16, :, :), [1 3 2]);
            aas = permute(states.Metabolite.counts(17:37, :, :), [1 3 2]);
            others = permute(states.Metabolite.counts(38:40, :, :), [1 3 2]);
           
            %% layout figure
            clf(figHandle);
            
            axesHandles = PlotUtil.multiElementPlot(figHandle, {10 * ones(3, 1) 10 * ones(3, 1)}, time([1 end]));
            
            %% plot data
            if nargin >= 5 && normalized
                plotFcn = @PlotUtil.plotNormalizedStackedArea;
            else
                plotFcn = @PlotUtil.plotStackedArea;
            end
            
            %ntp
            axesHandle = axesHandles{1}(1);
            h = plotFcn(axesHandle, time, ntps);
            legend(h, m.wholeCellModelIDs(m.ntpIndexs), 'Location', 'SouthEast');
            
            %ndp
            axesHandle = axesHandles{1}(2);
            h = plotFcn(axesHandle, time, ndps);
            legend(h, m.wholeCellModelIDs(m.ndpIndexs), 'Location', 'SouthEast');
            
            %nmp
            axesHandle = axesHandles{1}(3);
            h = plotFcn(axesHandle, time, nmps);
            legend(h, m.wholeCellModelIDs(m.nmpIndexs), 'Location', 'SouthEast');
            
            %dNTP
            axesHandle = axesHandles{2}(1);
            h = plotFcn(axesHandle, time, dntps);
            legend(h, m.wholeCellModelIDs(m.dntpIndexs), 'Location', 'SouthEast');
            
            %AA
            axesHandle = axesHandles{2}(2);
            plotFcn(axesHandle, time, aas);
            
            %other
            axesHandle = axesHandles{2}(3);
            h = plotFcn(axesHandle, time, others);
            legend(h, m.wholeCellModelIDs([m.phosphateIndexs; m.diphosphateIndexs; m.hydrogenIndexs]), 'Location', 'SouthEast');
            
            %labels
            for i = 1:numel(axesHandles)
                for j = 1:numel(axesHandles{i})
                    if nargin >= 5 && normalized
                        ylabel(axesHandles{i}(j), {'Relative' 'Copy Number'});
                        ylim(axesHandles{i}(j), [0 100]);
                    else
                        ylabel(axesHandles{i}(j), 'Concentration (mM)');
                    end
                end
            end
            
            %% ylabels
            PlotUtil.alignYAxesLabels(axesHandles);
            PlotUtil.offsetYAxes(axesHandles, 0.02);
            
            %% collect data
            figData = struct;
            figData.time = time;
            figData.ntps = ntps;
            figData.ndps = ndps;
            figData.nmps = nmps;
            figData.dntps = dntps;
            figData.aas = aas;
            figData.others = others;
        end
        
        function figData = metabolitesNormalized(sim, simBatchDir, simNum, figHandle)
            import edu.stanford.covert.cell.sim.analysis.SingleCell;
            
            figData = SingleCell.metabolites(sim, simBatchDir, simNum, figHandle, true);
        end
        
        function figData = ribosomes(sim, simBatchDir, simNum, figHandle, useGenomicCoordinates)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            %get constants
            g = sim.gene;
            c = sim.state('Chromosome');
            rib = sim.state('Ribosome');
            
            if nargin >= 5 && useGenomicCoordinates
                startCoordinates = g.startCoordinates(g.mRNAIndexs);
                lengths = g.lengths(g.mRNAIndexs);
                strands = g.strands(g.mRNAIndexs);
            else
                startCoordinates = cumsum(g.lengths(g.mRNAIndexs)+50);
                lengths = g.lengths(g.mRNAIndexs);
                strands = ones(size(g.mRNAIndexs));
            end
            
            %get data
            stateNames = {
                'Time'      'values'
                'Ribosome'  'states' 
                'Ribosome'  'boundMRNAs'
                'Ribosome'  'mRNAPositions'
                };
            states = SimulationEnsemble.load(simBatchDir, stateNames, [], [], 1, 'extract', simNum);
            time = permute(states.Time.values, [1 3 2]) / 3600;
            ribStates = permute(states.Ribosome.states, [1 3 2]);
            boundMRNAs = permute(states.Ribosome.boundMRNAs, [1 3 2]);
            mRNAPos = permute(states.Ribosome.mRNAPositions, [1 3 2]);
            
            %layout figure
            clf(figHandle);
            options = struct;
            axesHandle = PlotUtil.multiElementPlot(figHandle, 30, time([1 end]), options);
            
            %plot data
            colorOrder = PlotUtil.getRedGreenColorOrder(randperm(size(ribStates, 1)));
            for i = 1:size(ribStates, 1)
                pos = NaN(size(time));
                tfs = ribStates(i, :) == rib.activeValue & mRNAPos(i, :) > 0;
                pos(tfs) = ...
                    + startCoordinates(boundMRNAs(i, tfs)) ...
                    + (strands(boundMRNAs(i, tfs)) == 2) .* (lengths(boundMRNAs(i, tfs)) - 1) ...
                    + (2 * (2 - strands(boundMRNAs(i, tfs))) - 1) .* (3 * mRNAPos(i, tfs)' - 1);
                h = plot(axesHandle, time, pos);
                set(h, 'Color', colorOrder(i, :));
            end
            
            if nargin >= 5 && useGenomicCoordinates
                ylabel(axesHandle, 'Ribosome Position');
                ylim(axesHandle, [1 c.sequenceLen]);
                set(axesHandle, ...
                    'YTick', [1 c.terCPosition c.sequenceLen], ...
                    'YTickLabel', {'OriC' 'TerC' 'Oric'});
                PlotUtil.offsetYAxes(axesHandle, 0.02);
            else
                tmp = xlim(axesHandle);
                x = tmp(1) - 0.02 * range(tmp);
                xlim(axesHandle, [x tmp(2)]);
                ylim(axesHandle, [min(startCoordinates) max(startCoordinates + lengths - 1)]);
                line(x * ones(3 * numel(startCoordinates), 1), ...
                    reshape([startCoordinates  startCoordinates + lengths - 1  NaN(size(startCoordinates))]', [], 1), ...
                    'Parent', axesHandle, 'Color', 'k')
                set(axesHandle, 'YTick', [], 'visible', 'off');
                text(tmp(1) - 0.03 * range(tmp), mean(ylim(axesHandle)), 'Ribosome Position', ...
                    'Parent', axesHandle, ...
                    'rotation', 90, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom');
            end
            
            %% collect data
            figData = struct;
            figData.time = time;
            figData.ribStates = ribStates;
            figData.boundMRNAs = boundMRNAs;
            figData.mRNAPos = mRNAPos;
        end
        
        function figData = ribosomesGenomicCoordinates(sim, simBatchDir, simNum, figHandle)
            import edu.stanford.covert.cell.sim.analysis.SingleCell;
            
            figData = SingleCell.ribosomes(sim, simBatchDir, simNum, figHandle, true);
        end
    end
end
