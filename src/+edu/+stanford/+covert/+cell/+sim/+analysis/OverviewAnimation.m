%OverviewAnimation
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 10/6/2011
classdef OverviewAnimation < edu.stanford.covert.cell.sim.analysis.Animation
    properties (SetAccess = protected)
        title = 'Whole-Cell Animation'
        description = 'Animation of Mycoplasma genitalium whole-cell simulations.'
        author = 'Jonathan Karr, Covert Lab, Department of Bioengineering, Stanford University'

        displayAuthorList = false;
		useTransitions = true;
    end
    
    methods
        function this = OverviewAnimation(varargin)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            
            this = this@edu.stanford.covert.cell.sim.analysis.Animation(varargin{:});
            
            this.movieWidth = 800;
            this.movieHeight = 600;
            this.frameRate = 30;
        end
        
        function [times, data] = loadSimulationData(this)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            %% constants
            sim = this.simulation;
            comp = sim.compartment;
            mass = sim.state('Mass');
            met = sim.state('Metabolite');
            mr = sim.state('MetabolicReaction');
            pc = sim.state('ProteinComplex');
            
            %% metabolic map
            dbConnectionParameters = getConfig();
            database = edu.stanford.covert.db.MySQLDatabase(dbConnectionParameters);
            database.setNullValue(0);
            
            kbWID = edu.stanford.covert.cell.kb.KnowledgeBaseUtil.selectLatestKnowledgeBase(database);
            
            database.prepareStatement('CALL get_metabolicmapmetabolites("{Si}", null)', kbWID);
            data = database.query();
            metabolicMapMetabolites = struct;
            [~, metabolicMapMetabolites.idxs] = ismember(data.Metabolite, met.wholeCellModelIDs);
            [~, metabolicMapMetabolites.compIdxs] = ismember(data.Compartment, comp.wholeCellModelIDs);
            [~, metabolicMapMetabolites.uidxs] = ismember(metabolicMapMetabolites.idxs, unique(metabolicMapMetabolites.idxs));
            metabolicMapMetabolites.x = data.X;
            metabolicMapMetabolites.y = cellfun(@str2double, data.Y);
            
            database.prepareStatement('CALL get_metabolicmapreactions("{Si}", null)', kbWID);
            data = database.query();
            metabolicMapReactions = struct;
            [~, metabolicMapReactions.idxs] = ismember(data.Reaction, mr.reactionWholeCellModelIDs);
            [~, metabolicMapReactions.uidxs] = ismember(metabolicMapReactions.idxs, unique(metabolicMapReactions.idxs));
            metabolicMapReactions.path = data.Path;
            metabolicMapReactions.labelX = data.LabelX;
            metabolicMapReactions.labelY = data.LabelY;
            metabolicMapReactions.valueX = data.ValueX;
            metabolicMapReactions.valueY = data.ValueY;
            
            %% data
            stateNames = {
                'mass'
                'growth_rate'
                'dnaA_box1'
                'dnaA_box2'
                'dnaA_box3'
                'dnaA_box4'
                'dnaA_box5'
                'ploidy'
                'superhelicity'
                'helicase1'
                'helicase2'
                'rnas'
                'mrnas'
                'rrnas'
                'srnas'
                'trnas'
                'immatureRnas'
                'matureMonomers'
                'immatureMonomers'
                'matureComplexs'
                'immatureComplexs'
                'lipids'
                'polysaccharides'
                'atp'
                'ntps'
                'dntps'
                'amino_acids'
                'atpUsage'
                'gtpUsage'
                'ftsZ'
                'ftsZRing1st'
                'ftsZRing2st'
                'ftsZRing2bt'
                'ftsZRingRbt'
                'pinchedDiameter'
                };
            ensemble = SimulationEnsemble(this.simGroupId, stateNames, [], this.simId);
            
            ensemble.stateData.values(ensemble.getPropertyIndices('mass'), :, :, :) = ...
                ensemble.stateData.values(ensemble.getPropertyIndices('mass'), :, :, :) * ...
                mass.cellInitialDryWeight / (1 - mass.fractionWetWeight) * 1e15;
            ensemble.stateData.values(ensemble.getPropertyIndices('growth_rate'), :, :, :) = ...
                ensemble.stateData.values(ensemble.getPropertyIndices('growth_rate'), :, :, :) * ...
                3600 * mass.cellInitialDryWeight / (1 - mass.fractionWetWeight) * 1e15;
            
            ensemble.stateData.values = ensemble.stateData.values(:, :, 2:ensemble.stateData.simulationEndTimes(1)+1, :);
            times = ensemble.stateData.time(2:ensemble.stateData.simulationEndTimes(1)+1);
            
            stateNames = {
                'Geometry'           'width'               ':'             ':'
                'Geometry'           'totalLength'         ':'             ':'
                'Geometry'           'pinchedDiameter'     ':'             ':'
                'Geometry'           'cylindricalLength'   ':'             ':'
                'Chromosome'         'polymerizedRegions'  ':'             ':'
                'Chromosome'         'monomerBoundSites'   ':'             ':'
                'Chromosome'         'complexBoundSites'   ':'             ':'
                'Ribosome'           'states'              ':'             ':'
                'Ribosome'           'boundMRNAs'          ':'             ':'
                'Ribosome'           'mRNAPositions'       ':'             ':'
                'Ribosome'           'tmRNAPositions'      ':'             ':'
                'Metabolite'         'processUsages'       met.ntpIndexs'  ':'
                'Metabolite'         'counts'              ':'             ':'
                'MetabolicReaction'  'fluxs'               ':'             ':'
                'ProteinComplex'     'counts'              pc.matureIndexs([pc.ribosome30SIndexs; pc.ribosome50SIndexs])  comp.cytosolIndexs
                };
            states = SimulationEnsemble.load(this.simGroupId, stateNames, [], [], 1, 'extract', 1);
            states.Metabolite.processUsages = full(states.Metabolite.processUsages);
            states.Metabolite.counts = full(states.Metabolite.counts);
            states.Geometry.septumLength = (states.Geometry.width - states.Geometry.pinchedDiameter) / 2;            
            states.ProteinComplex.ribosome30SCounts = states.ProteinComplex.counts(1, :, :, :);
            states.ProteinComplex.ribosome50SCounts = states.ProteinComplex.counts(2, :, :, :);
            states.ProteinComplex = rmfield(states.ProteinComplex, 'counts');
            
            %%
            data = struct;
            data.ensemble = ensemble;
            data.states = states;
            data.metabolicMapMetabolites = metabolicMapMetabolites;
            data.metabolicMapReactions = metabolicMapReactions;
        end
        
        function svg = drawFrame(this, timeIdx)
            import edu.stanford.covert.cell.sim.util.SimulationEnsemble;
            
            sim = this.simulation;
            chr = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            sequenceLen = chr.sequenceLen;
            
            ensemble = this.states.ensemble;
            states = this.states.states;
            
            t1 = find(ensemble.stateData.values(ensemble.getPropertyIndices('ploidy'), :, :) > 1, 1, 'first');
            t2 = find(states.Geometry.pinchedDiameter < states.Geometry.width, 1, 'first') + 30 * 60;
            transitionTimes = [
                t2
                t1 + 60 * 60
                t1 + 30 * 60
                t1 + (t2 - t1) * 1 / 3;
                t1 + (t2 - t1) * 2 / 3;
                ];
			if ~this.useTransitions
				transitionTimes(:) = Inf;
			end
            tickLen = 4;
            tAvg = 500;
            
            %% open svg
            svg = [];
            svg = [svg sprintf('<?xml version="1.0" standalone="no"?>\n')];
            svg = [svg sprintf('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n')];
            svg = [svg sprintf('<svg width="%d" height="%d" version="1.1" xmlns="http://www.w3.org/2000/svg" viewBox="%d %d %d %d">\n', 800, 600, 0, 0, 800, 600)];
            
            %definitions
            svg = [svg sprintf('<defs>\n')];
            svg = [svg sprintf('  <linearGradient id="cellBg" x1="0%%" y1="0%%" x2="0%%" y2="100%%">\n')];
            svg = [svg sprintf('    <stop style="stop-color:#97b0fb;stop-opacity:0.05" offset="0"/>\n')];
            svg = [svg sprintf('    <stop style="stop-color:#89b8d6;stop-opacity:1" offset="1"/>\n')];
            svg = [svg sprintf('  </linearGradient>')];
            svg = [svg sprintf('  <linearGradient id="barBg" x1="0%%" y1="0%%" x2="0%%" y2="100%%">\n')];
            svg = [svg sprintf('    <stop style="stop-color:#97b0fb;stop-opacity:0.25" offset="0"/>\n')];
            svg = [svg sprintf('    <stop style="stop-color:#89b8d6;stop-opacity:1" offset="1"/>\n')];
            svg = [svg sprintf('  </linearGradient>')];
            svg = [svg sprintf('  <radialGradient id="metabolite" x1="0%%" y1="0%%" x2="0%%" y2="100%%">\n')];
            svg = [svg sprintf('    <stop style="stop-color:#97b0fb;stop-opacity:0.25" offset="0"/>\n')];
            svg = [svg sprintf('    <stop style="stop-color:#89b8d6;stop-opacity:1" offset="1"/>\n')];
            svg = [svg sprintf('  </radialGradient>')];
            svg = [svg sprintf('  <radialGradient id="activeRibosome" x1="0%%" y1="0%%" x2="0%%" y2="100%%">\n')];
            svg = [svg sprintf('    <stop style="stop-color:#FFFFFF;stop-opacity:0.25" offset="0"/>\n')];
            svg = [svg sprintf('    <stop style="stop-color:#8BF763;stop-opacity:1" offset="1"/>\n')];
            svg = [svg sprintf('  </radialGradient>')];
            svg = [svg sprintf('  <radialGradient id="stalledRibosome" x1="0%%" y1="0%%" x2="0%%" y2="100%%">\n')];
            svg = [svg sprintf('    <stop style="stop-color:#FFFFFF;stop-opacity:0.25" offset="0"/>\n')];
            svg = [svg sprintf('    <stop style="stop-color:#F96464;stop-opacity:1" offset="1"/>\n')];
            svg = [svg sprintf('  </radialGradient>')];
            svg = [svg sprintf('  <radialGradient id="free30SRibosome" x1="0%%" y1="0%%" x2="0%%" y2="100%%">\n')];
            svg = [svg sprintf('    <stop style="stop-color:#FFFFFF;stop-opacity:0.25" offset="0"/>\n')];
            svg = [svg sprintf('    <stop style="stop-color:#CD64F9;stop-opacity:1" offset="1"/>\n')];
            svg = [svg sprintf('  </radialGradient>')];
            svg = [svg sprintf('  <radialGradient id="free50SRibosome" x1="0%%" y1="0%%" x2="0%%" y2="100%%">\n')];
            svg = [svg sprintf('    <stop style="stop-color:#FFFFFF;stop-opacity:0.25" offset="0"/>\n')];
            svg = [svg sprintf('    <stop style="stop-color:#F9B264;stop-opacity:1" offset="1"/>\n')];
            svg = [svg sprintf('  </radialGradient>')];
            svg = [svg sprintf('  <marker id="EndArrowPos"\n')];
            svg = [svg sprintf('    viewBox="0 0 10 10" refX="0" refY="5" \n')];
            svg = [svg sprintf('    markerUnits="strokeWidth"\n')];
            svg = [svg sprintf('    markerWidth="4" markerHeight="3"\n')];
            svg = [svg sprintf('    orient="auto">\n')];
            svg = [svg sprintf('    <path d="M 0 0 L 10 5 L 0 10 z" style="fill:green;" />\n')];
            svg = [svg sprintf('  </marker>\n')];
            svg = [svg sprintf('  <marker id="StartArrowNeg"\n')];
            svg = [svg sprintf('    viewBox="0 0 10 10" refX="0" refY="5" \n')];
            svg = [svg sprintf('    markerUnits="strokeWidth"\n')];
            svg = [svg sprintf('    markerWidth="4" markerHeight="3"\n')];
            svg = [svg sprintf('    orient="auto">\n')];
            svg = [svg sprintf('    <path d="M 10 0 L 0 5 L 10 10 z" style="fill:#ff241E;" />\n')];
            svg = [svg sprintf('  </marker>\n')];
            svg = [svg sprintf('</defs>\n')];
            
            %style
            svg = [svg sprintf('<style type="text/css">\n')];
            svg = [svg sprintf('#background{fill:white; stroke:none;}\n')];
            svg = [svg sprintf('text{font-family:arial;}\n')];
            svg = [svg sprintf('#title text{font-variant:small-caps; fill:white; font-weight:bold;}\n')];
            svg = [svg sprintf('#time text{font-variant:small-caps; fill:white; font-weight:bold;}\n')];
            svg = [svg sprintf('text.sectionLabel{font-size:20px; font-weight:bold; text-anchor:middle;}\n')];
            svg = [svg sprintf('text.axisLabel{font-size:15px; font-weight:bold;}\n')];
            svg = [svg sprintf('text.legendLabel{font-size:15px; font-weight:normal;}\n')];            
            svg = [svg sprintf('text.tickLabel{font-size:12px; font-weight:bold;}\n')];
            svg = [svg sprintf('line.axis{stroke:black; stroke-width:1;}\n')];
            svg = [svg sprintf('circle.activeRibosome, circle.stalledRibosome, circle.free30SRibosome, circle.free50SRibosome{stroke:#666666; stroke-width:0.5;}\n')];
            svg = [svg sprintf('circle.activeRibosome{fill:url(#activeRibosome);}\n')];
            svg = [svg sprintf('circle.stalledRibosome{fill:url(#stalledRibosome);}\n')];
            svg = [svg sprintf('circle.free30SRibosome{fill:url(#free30SRibosome);}\n')];
            svg = [svg sprintf('circle.free50SRibosome{fill:url(#free50SRibosome);}\n')];            
            svg = [svg sprintf('rect.activeRibosome, rect.stalledRibosome, rect.free30SRibosome, rect.free50SRibosome{stroke:none; stroke-width:0.5; height:2}\n')];
            svg = [svg sprintf('rect.activeRibosome{fill:#8BF763;}\n')];
            svg = [svg sprintf('rect.stalledRibosome{fill:#F96464;}\n')];
            svg = [svg sprintf('rect.free30SRibosome{fill:#CD64F9;}\n')];
            svg = [svg sprintf('rect.free50SRibosome{fill:#F9B264;}\n')];
            
            svg = [svg sprintf('</style>\n')];
            
            %white background
            svg = [svg sprintf('<rect x="%d" y="%d" width="%d" height="%d" id="background"/>\n', ...
                0, 0, 800, 600)];
            
            %% title
            svg = [svg sprintf('<g transform="translate(%d, %d)" id="title">\n', 0, 0)];
            svg = [svg sprintf('  <rect x="%d" y="%d" width="%d" height="%d" style="fill:#ff241E; stroke:none"/>\n', ...
                0, 0, 800, 50)];
            if this.displayAuthorList
                svg = [svg sprintf('  <text x="%d" y="%d" style="font-size:24px; text-anchor:start">%s</text>\n', ...
                    4, 22, 'A Whole-Cell Model Predicts Phenotype from Genotype')];
                svg = [svg sprintf('  <text x="%d" y="%d" style="font-size:18px; text-anchor:start">%s</text>\n', ...
                    4, 42, 'Karr JR, Sanghvi JC, Macklin DN, Gutschow MV, Jacobs JM, Bolival B, Assad-Garcia N, Glass JI, and Covert MW')];
            else
                svg = [svg sprintf('  <text x="%d" y="%d" style="font-size:24px; text-anchor:start">%s</text>\n', ...
                    4, 32, 'A Whole-Cell Model Predicts Phenotype from Genotype')];
            end
            svg = [svg sprintf('</g>\n')];
            
            %% footer
            svg = [svg sprintf('<g transform="translate(%d, %d)" id="time">\n', 0, 570)];
            svg = [svg sprintf('  <rect x="%d" y="%d" width="%d" height="%d" style="fill:#ff241E; stroke:none"/>\n', ...
                0, 0, 800, 30)];
            svg = [svg sprintf('  <text x="%d" y="%d" style="font-size:24px; text-anchor:start">%s</text>\n', ...
                664, 23, 'Time')];
            svg = [svg sprintf('  <text x="%d" y="%d" style="font-size:18px; text-anchor:end">%s</text>\n', ...
                796, 23, datestr(this.times(timeIdx) / (24 * 60 * 60), 'HH:MM:SS'))];            
            svg = [svg sprintf('  <text x="%d" y="%d" style="font-size:24px; text-anchor:start">%s</text>\n', ...
                4, 23, 'wholecell.stanford.edu')];
            svg = [svg sprintf('</g>\n')];
            
            %% top left
            H = 150;
            W = 230;
            
            ftsZRing = this.simulation.state('FtsZRing');
            timeIdx2 = min(timeIdx, numel(this.times) - 1);
            
            svg = [svg sprintf('<g transform="translate(%d, %d)">\n', 0, 100)];
            svg = [svg sprintf('  <text x="%d" y="%d" class="sectionLabel">%s</text>\n', 230/2, -20, 'Cell Growth')];
            
            rib = sim.state('Ribosome');
            ribStates = this.states.states.Ribosome.states(:, :, timeIdx);
            nActive = sum(ribStates == rib.activeValue);
            nStalled = sum(ribStates == rib.stalledValue);
            nFree30 = this.states.states.ProteinComplex.ribosome30SCounts(:, :, timeIdx);
            nFree50 = this.states.states.ProteinComplex.ribosome50SCounts(:, :, timeIdx);
            
            scale = min(W / max(states.Geometry.totalLength(:, :, 1:min(end, transitionTimes(1)))), H / max(states.Geometry.width(:, :, 1:min(end, transitionTimes(1)))));
            svg = [svg this.drawCell(...
                scale * states.Geometry.width(timeIdx), ...
                scale * states.Geometry.cylindricalLength(timeIdx), ...
                scale * states.Geometry.septumLength(timeIdx), ...
                W, H, nActive, nStalled, nFree30, nFree50, ...
                ftsZRing.calcNumEdges(states.Geometry.pinchedDiameter(timeIdx2), ftsZRing.filamentLengthInNm), ...
                ensemble.stateData.values(ensemble.getPropertyIndices('ftsZRing1st'), :, timeIdx2, 1), ...
                ensemble.stateData.values(ensemble.getPropertyIndices('ftsZRing2st'), :, timeIdx2, 1), ...
                ensemble.stateData.values(ensemble.getPropertyIndices('ftsZRing2bt'), :, timeIdx2, 1), ...
                ensemble.stateData.values(ensemble.getPropertyIndices('ftsZRingRbt'), :, timeIdx2, 1), ...
                this.times(timeIdx) >= transitionTimes(1))];
            
            svg = [svg sprintf('</g>\n')];            
            
            %% bottom left
            if this.times(timeIdx) < transitionTimes(2)
                H = 150;
                W = 230;
                r = H/2;
                
                svg = [svg sprintf('<g transform="translate(%d, %d)">\n', 0, 360)];
                svg = [svg sprintf('  <text x="%d" y="%d" class="sectionLabel">%s</text>\n', 230/2, -20, 'Chromosome')];
                
                %polymerized Regions
                [posStrnds, lens] = find(states.Chromosome.polymerizedRegions);
                lens = lens(posStrnds(:, 3) == timeIdx, :);
                posStrnds = posStrnds(posStrnds(:, 3) == timeIdx, 1:2);
                twoChrs = any(posStrnds(:, 2) > 2); 
                svg = [svg sprintf(this.drawArc(posStrnds, lens, W/2, H/2+5, r, 5, 'grey', 1, twoChrs, sequenceLen))];
                
                %bound proteins
                [monPosStrnds, monIdxs] = find(states.Chromosome.monomerBoundSites);
                [cpxPosStrnds, cpxIdxs] = find(states.Chromosome.complexBoundSites);
                posStrnds = [
                    monPosStrnds(monPosStrnds(:, 3) == timeIdx, 1:2)
                    cpxPosStrnds(cpxPosStrnds(:, 3) == timeIdx & ~ismembc(cpxIdxs, pc.rnaPolymeraseIndexs), 1:2)
                    ];
                lens = [
                    chr.monomerDNAFootprints(monIdxs(monPosStrnds(:, 3) == timeIdx, :))
                    chr.complexDNAFootprints(cpxIdxs(cpxPosStrnds(:, 3) == timeIdx, :))
                    ];
                svg = [svg sprintf(this.drawArc(posStrnds, 2 * lens, W/2, H/2+5, r, 10, '#ff241E', 2, twoChrs, sequenceLen))];
                
                posStrnds = cpxPosStrnds(cpxPosStrnds(:, 3) == timeIdx & ismembc(cpxIdxs, pc.rnaPolymeraseIndexs), 1:2);
                lens = 10000 * ones(size(posStrnds, 1), 1);
                svg = [svg sprintf(this.drawArc(posStrnds, lens, W/2, H/2+5, r, 10, 'darkblue', 2, twoChrs, sequenceLen))];
                
                %key
                svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:grey; stroke-width:2"/>\n', 20, H+36, 30, H+36)];
                svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:#ff241E; stroke-width:2"/>\n', 75, H+36, 85, H+36)];
                svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:darkblue; stroke-width:2"/>\n', 145, H+36, 155, H+36)];
                svg = [svg sprintf('  <text x="%d" y="%d" class="legendLabel" style="text-anchor:start">%s</text>\n', 35, H+41, 'DNA')];
                svg = [svg sprintf('  <text x="%d" y="%d" class="legendLabel" style="text-anchor:start">%s</text>\n', 90, H+41, 'Protein')];
                svg = [svg sprintf('  <text x="%d" y="%d" class="legendLabel" style="text-anchor:start">%s</text>\n', 160, H+41, 'RNA Pol')];
                svg = [svg sprintf('</g>\n')];
            else
                H = 150;
                W = 165;
                dx = 5;
                
                svg = [svg sprintf('<g transform="translate(%d, %d)">\n', 65, 360)];
                svg = [svg sprintf('  <text x="%d" y="%d" class="sectionLabel">%s</text>\n', 230/2-65, -20, 'DNA Polymerase')];
                
                %axes
                svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:black; stroke-width:1"/>\n', -dx, 0, -dx, H)];
                svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:black; stroke-width:1"/>\n', -dx, 0, -dx-tickLen, 0)];
                svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:black; stroke-width:1"/>\n', -dx, H/2, -dx-tickLen, H/2)];
                svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:black; stroke-width:1"/>\n', -dx, H, -dx-tickLen, H)];
                svg = [svg sprintf('  <text x="0" y="0" class="axisLabel" transform="translate(%d, %d) rotate(270)" style="text-anchor:middle">%s</text>\n', ...
                    -dx-36, H/2, 'Position')];
                svg = [svg sprintf('  <text x="%d" y="%d" class="tickLabel" style="text-anchor:end">%s</text>\n', ...
                    -dx-tickLen-2, 0+3, 'oriC')];
                svg = [svg sprintf('  <text x="%d" y="%d" class="tickLabel" style="text-anchor:end">%s</text>\n', ...
                    -dx-tickLen-2, H/2+3, 'terC')];
                svg = [svg sprintf('  <text x="%d" y="%d" class="tickLabel" style="text-anchor:end">%s</text>\n', ...
                    -dx-tickLen-2, H+3, 'oriC')];
                
                svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:black; stroke-width:1"/>\n', 0, H+dx, W, H+dx)];
                ticks = [0 5];
                if this.times(end) / 3600 >= 10
                    ticks(end) = 10;
                end
                for i = 1:numel(ticks)
                    svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:black; stroke-width:1"/>\n', ...
                        W * ticks(i) / (this.times(end) / 3600), H+dx, W * ticks(i) / (this.times(end) / 3600), H+dx+tickLen)];
                    svg = [svg sprintf('  <text x="%d" y="%d" class="tickLabel" style="text-anchor:middle">%d</text>\n', ...
                        W * ticks(i) / (this.times(end) / 3600), H+dx+tickLen+16, ticks(i))];
                end
                svg = [svg sprintf('  <text x="%d" y="%d" class="axisLabel" style="text-anchor:middle">%s</text>\n', W/2, H+dx+36, 'Time (h)')];
                
                %data
                for k = 1:2
                    data = permute(ensemble.stateData.values(ensemble.getPropertyIndices(['helicase' num2str(k)]), :, :, 1), [4 3 1 2]);
                    
                    minVal = 1;
                    maxVal = sequenceLen;
                    
					if ~any(~isnan(data(1:timeIdx)) & max(0, data(1:timeIdx)))
						continue;
					end
					
                    path = sprintf('M %f %f ', W * this.times(find(~isnan(data), 1, 'first')) / this.times(end), H * (1 - (data(find(~isnan(data), 1)) - minVal) / (maxVal - minVal)));
                    idxs = find(diff(data)~=0 & ~isnan(diff(data)) & ~isnan(data(:, 2:end)) & max(0, data(:, 2:end)));
                    if numel(idxs) > tAvg
                        idxs = idxs(1:floor(numel(idxs) / tAvg):end);
                    end
                    idxs = idxs(idxs + 1 <= timeIdx);
                    idxs(end+1) = find(~isnan(data(1:timeIdx)), 1, 'last') - 1;
                    for j = 1:numel(idxs)
                        if idxs(j) + 1 > timeIdx
                            break;
                        end
                        path = [path sprintf('L %f %f ', W * this.times(idxs(j)+1) / this.times(end), H * (1 - (data(idxs(j)+1) - minVal) / (maxVal - minVal)))];
                    end
                    
                    svg = [svg sprintf('  <path d="%s" style="stroke:darkblue; stroke-width:2px; fill:none;"/>', path)];
                end
                
                svg = [svg sprintf('</g>\n')];
            end
            
            %% top middle
            H = 200;
            W = 280;
            
            if this.times(timeIdx) < transitionTimes(3)
                svg = [svg sprintf('<g transform="translate(%d, %d)">\n', 255, 100)];
                svg = [svg sprintf('  <text x="%d" y="%d" class="sectionLabel">%s</text>\n', 230/2+25, -20, 'Metabolism')];
                
                mapMets = this.states.metabolicMapMetabolites;
                mapRxns = this.states.metabolicMapReactions;

                minX = min(mapMets.x);
                minY = min(mapMets.y);
                scaleX = W / range(mapMets.x);
                scaleY = H / range(mapMets.y);
                mapMets.x = scaleX * (mapMets.x - minX);
                mapMets.y = scaleY * (mapMets.y - minY);
                for i = 1:numel(mapRxns.path)
                    [starts, ends, tokens] = regexp(mapRxns.path{i}, '(-*\d*\.*\d*),(-*\d*\.*\d*)', 'start', 'end', 'tokens');
                    for j = 1:numel(starts)
                        valX = num2str(scaleX * (str2double(tokens{j}{1}) - minX));
                        valY = num2str(scaleY * (str2double(tokens{j}{2}) - minY));
                        mapRxns.path{i} = [...
                            mapRxns.path{i}(1:starts(j)-1) ...
                            valX ...
                            ',' ...
                            valY ...
                            mapRxns.path{i}(ends(j)+1:end) ...
                            ];
                        starts = starts + numel(valX) + numel(valY) - numel(tokens{j}{1}) - numel(tokens{j}{2});
                        ends = ends + numel(valX) + numel(valY) - numel(tokens{j}{1}) - numel(tokens{j}{2});
                    end
                end
                
                data = 2 + 7 * (states.Metabolite.counts(:, :, timeIdx) - min(states.Metabolite.counts, [], 3)) ./ ...
                    (max(states.Metabolite.counts, [], 3) -  min(states.Metabolite.counts, [], 3));
                data(isnan(data)) = 5;
                for i = 1:numel(mapMets.idxs)
                    svg = [svg sprintf('    <circle cx="%d" cy="%d" r="%d" style="stroke:none; fill:white"/>\n', ...
                        mapMets.x(i), mapMets.y(i), data(mapMets.idxs(i), mapMets.compIdxs(i)))];
                    svg = [svg sprintf('    <circle cx="%d" cy="%d" r="%d" style="stroke:darkblue; stroke-width:0.5px; fill:url(#metabolite)"/>\n', ...
                        mapMets.x(i), mapMets.y(i), data(mapMets.idxs(i), mapMets.compIdxs(i)))];
                end
                
                data1 = states.MetabolicReaction.fluxs(:, :, timeIdx);
                data2 = 0.25 + 0.75 * (abs(data1) - min(abs(states.MetabolicReaction.fluxs), [], 3)) ./ ...
                    (max(abs(states.MetabolicReaction.fluxs), [], 3) -  min(abs(states.MetabolicReaction.fluxs), [], 3));
                data2(isnan(data2)) = 0.25;
                for i = 1:numel(mapRxns.idxs)
                    if data1(mapRxns.idxs(i)) == 0
                        data2(mapRxns.idxs(i)) = 1;
                        stroke = 'lightgrey';
                        style = 'stroke-dasharray: 2, 2;';
                        markerStart = 'none';
                        markerEnd = 'none';
                    elseif data1(mapRxns.idxs(i)) > 0
                        stroke = 'green';
                        style = '';
                        markerStart = 'none';
                        markerEnd = 'url(#EndArrowPos)';
                    else
                        stroke = '#ff241E';
                        style = '';
                        markerStart = 'url(#StartArrowNeg)';
                        markerEnd = 'none';
                    end
                    svg = [svg sprintf('    <path d="%s" style="fill:none; stroke:%s; stroke-width:1px; opacity:%0.3f; %s; marker-start:%s; marker-end:%s"/>\n', ...
                        mapRxns.path{i}, stroke, data2(mapRxns.idxs(i)), style, markerStart, markerEnd)];
                end
                
                svg = [svg sprintf('</g>\n')];
            else
                H = 150;
                W = 230;
                dx = 5;
                
                svg = [svg sprintf('<g transform="translate(%d, %d)">\n', 280, 100)];
                svg = [svg sprintf('  <text x="%d" y="%d" class="sectionLabel">%s</text>\n', 230/2, -20, 'Energy Usage')];
                
                %axes
                svg = [svg sprintf('  <rect x="%d" y="%d" width="%d" height="%d" style="stroke:none; fill:rgb(210,210,210); stroke-width:1"/>\n', 0, 0, W, H)];
                
                %data
                data = permute(states.Metabolite.processUsages, [2 3 1]);
                nProcesses = size(data, 1);
                data = [data NaN(nProcesses, tAvg-mod(size(data, 2), tAvg))];
                data = permute(nanmean(reshape(data, [nProcesses tAvg size(data, 2) / tAvg]), 2), [1 3 2]);
                sums = sum(data, 2);
                data = data(sums > 0, :);
                sums = sums(sums > 0);
                nProcesses = numel(sums);
                [~, order] = sort(sums);
                data = data(order, :);
                
                data = log(max(1, data));
                data = data - min(data(:));
                maxVal = max(data(:));
                r = round(255 * max(0, (data - maxVal / 2) / (maxVal / 2)));
                g = round(255 * max(0, (maxVal / 2 - data) / (maxVal / 2)));
                rows = floor(H * (0:nProcesses) / nProcesses);
                rows(end) = H;
                rowHeights = diff(rows)+1;
                cols = floor(W * (0:ceil(numel(this.times) / tAvg)) / ceil(numel(this.times) / tAvg));
                cols(end) = W;
                colWidths = diff(cols) + 1;
                colWidths(ceil(timeIdx / tAvg)) = colWidths(ceil(timeIdx / tAvg)) * mod(timeIdx / tAvg, 1);
                if ceil(timeIdx / tAvg) > 1
                    colWidths(ceil(timeIdx / tAvg) - 1) = ...
                        + colWidths(ceil(timeIdx / tAvg) - 1) ...
                        - 1 ...
                        + min(1, colWidths(ceil(timeIdx / tAvg)));
                end
                for i = 1:ceil(timeIdx / tAvg)
                    for j = 1:nProcesses
                        svg = [svg sprintf('  <rect x="%d" y="%d" width="%d" height="%d" style="stroke:none; stroke-width:1px; fill:rgb(%d, %d, 0);"/>\n', ...
                            cols(i), rows(j), colWidths(i), rowHeights(j), ...
                            r(j, i), g(j, i))];
                    end
                end
                
                %axes
                svg = [svg sprintf('  <rect x="%d" y="%d" width="%d" height="%d" style="stroke:black; fill:none; stroke-width:1"/>\n', 0, 0, W, H)];
                svg = [svg sprintf('  <text x="0" y="0" class="axisLabel" transform="translate(%d, %d) rotate(270)" style="text-anchor:middle">%s</text>\n', ...
                    -6, H/2, 'Process')];
                
                svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:black; stroke-width:1"/>\n', 0, H+dx, W, H+dx)];
                ticks = [0 5];
                if this.times(end) / 3600 >= 10
                    ticks(end) = 10;
                end
                for i = 1:numel(ticks)
                    svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:black; stroke-width:1"/>\n', ...
                        W * ticks(i) / (this.times(end) / 3600), H+dx, W * ticks(i) / (this.times(end) / 3600), H+dx+tickLen)];
                    svg = [svg sprintf('  <text x="%d" y="%d" class="tickLabel" style="text-anchor:middle">%d</text>\n', ...
                        W * ticks(i) / (this.times(end) / 3600), H+dx+tickLen+16, ticks(i))];
                end
                svg = [svg sprintf('  <text x="%d" y="%d" class="axisLabel" style="text-anchor:middle">%s</text>\n', W/2, H+dx+36, 'Time (h)')];
                
                svg = [svg sprintf('</g>\n')];
            end
            
            %% bottom middle
            if this.times(timeIdx) < transitionTimes(4)
                H = 150;
                W = 230;
                dx = 5;
                w = (W - dx * (5 - 1)) / 5;
                svg = [svg sprintf('<g transform="translate(%d, %d)">\n', 280, 360)];
                svg = [svg sprintf('  <text x="%d" y="%d" class="sectionLabel">%s</text>\n', 230/2, -20, 'Replication Initiation')];
                svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="fill:none; stroke:black;"/>\n', ...
                    -dx, 0, -dx, H)];
                svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="fill:none; stroke:black;"/>\n', ...
                    -dx, 0, -dx-tickLen, 0)];
                svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="fill:none; stroke:black;"/>\n', ...
                    -dx, H, -dx-tickLen, H)];
                svg = [svg sprintf('  <text x="0" y="0" class="axisLabel" transform="translate(%d, %d) rotate(270)" style="text-anchor:middle">%s</text>\n', -dx-tickLen-8, H/2, 'DnaA')];
                svg = [svg sprintf('  <text x="%d" y="%d" class="tickLabel" style="text-anchor:end">%d</text>\n', -dx-tickLen-2, 3, 7)];
                svg = [svg sprintf('  <text x="%d" y="%d" class="tickLabel" style="text-anchor:end">%d</text>\n', -dx-tickLen-2, H+3, 0)];
                
                for i = 1:5
                    data = permute(ensemble.stateData.values(ensemble.getPropertyIndices(['dnaA_box' num2str(i)]), :, timeIdx, 1), [4 3 1 2]);
                    
                    svg = [svg sprintf('  <rect x="%d" y="%d" width="%d" height="%d" style="fill:url(#barBg); stroke:darkblue;"/>\n', ...
                        round((i - 1) * (w + dx)), H * (1 - data/7), round(w), H * data / 7)];
                    svg = [svg sprintf('  <line x1="%d" y1="%d" x2="%d" y2="%d" style="fill:none; stroke:black;"/>\n', ...
                        round((i - 1) * (w + dx)), H+dx, round((i - 1) * (w + dx) + w), H+dx)];
                    svg = [svg sprintf('  <text x="%d" y="%d" class="tickLabel" style="text-anchor:middle">R%d</text>\n', ...
                        round((i - 1) * (w + dx) + w/2), H+dx+16, i)];
                end
                svg = [svg sprintf('  <text x="%d" y="%d" class="axisLabel" style="text-anchor:middle">%s</text>\n', W/2, H+dx+36, 'Functional DnaA Boxes')];
                svg = [svg sprintf('</g>\n')];
            else
                H = 155;
                W = 230;
                dx = 5;
                
                lens = pm.lengths(pm.nascentIndexs);
                
                nRows = 20;
                nCols = ceil(numel(lens) / nRows);
                
                svg = [svg sprintf('<g transform="translate(%d, %d)">\n', 280, 360)];
                svg = [svg sprintf('  <text x="%d" y="%d" class="sectionLabel">%s</text>\n', 230/2, -20, 'Translation')];
                
                h = (H - dx * (nRows - 1)) / nRows;
                for i = 1:nRows
                    if i == nRows
                        norm = sum(lens(1:nCols));
                    else
                        norm = sum(lens((i - 1) * nCols + 1:min(i * nCols, numel(lens))));
                    end
                    w = (W - dx * (nCols - 1)) / norm;
                    for j = (i - 1) * nCols + 1:min(i * nCols, numel(lens))
                        svg = [svg sprintf('  <line x1="%d" x2="%d" y1="%d" y2="%d" style="fill:none; stroke:grey; stroke-width:1px;"/>\n', ...
                            dx * (j - (i - 1) * nCols - 1) + w * sum(lens((i - 1) * nCols + 1:j-1)), ...
                            dx * (j - (i - 1) * nCols - 1) + w * sum(lens((i - 1) * nCols + 1:j)), ...
                            (i - 1) * (h + dx) + h, (i - 1) * (h + dx) + h)];
                        
                        nascentLens = max(states.Ribosome.mRNAPositions(states.Ribosome.boundMRNAs(:, :, timeIdx) == j, :, timeIdx));
                        abortedLens = max(states.Ribosome.tmRNAPositions(states.Ribosome.boundMRNAs(:, :, timeIdx) == j, :, timeIdx));
                        for k = 1:numel(nascentLens)
                            if abortedLens(k) == 0
                                color = 'darkblue';
                            else
                                color = '#ff241E';
                            end
                            
                            svg = [svg sprintf('  <rect x="%d" width="%d" y="%d" height="%d" style="fill:%s; stroke:none; stroke-width:1px;"/>\n', ...
                                dx * (j - (i - 1) * nCols - 1) + w * sum(lens((i - 1) * nCols + 1:j-1)), ...
                                w * nascentLens(k), ...
                                (i - 1) * (h + dx), h, color)];
                        end
                    end
                end
                
                svg = [svg sprintf('</g>\n')];
            end
            
            %% right
            svg = [svg sprintf('<g transform="translate(%d, %d)" id="right">\n', 650, 70)];
            
            %data
            if this.times(timeIdx) < transitionTimes(5)
                simIdxs = 1;
                states = {
                    'growth_rate'       {'Growth' '(fg h<tspan dy="-8" style="font-size:60%">-1</tspan><tspan dy="8">)</tspan>'}
                    'matureComplexs'    {'Cpx'}
                    'matureMonomers'    {'Mon'}
                    'mrnas'             {'mRNA'}
                    'rrnas'             {'rRNA'}
                    'trnas'             {'tRNA'}
                    'amino_acids'       {'Amino' 'Acid'}
                    'ntps'              {'NTP'}
                    };
            else
                simIdxs = 1:numel(this.simId);
                states = {
                    'mass'              {'Mass' '(fg)'}
                    'growth_rate'       {'Growth' '(fg h<tspan dy="-8" style="font-size:60%">-1</tspan><tspan dy="8">)</tspan>'}
                    'matureComplexs'    {'Cpx'}
                    'matureMonomers'    {'Mon'}
                    'matureMonomers'    {'RNA'}
                    'ploidy'            {'Chr' 'Copy' 'No'}
                    'superhelicity'     {'Super-' 'helicity'}
                    'dntps'             {'dNTP'}
                    };
            end
            dy = 16;
            H = 430;
            h = (H - dy * (size(states, 1) - 1)) / size(states, 1);
            w = 132;
            for i = 1:size(states, 1)
                y = round((i-1) * (h + dy));
                %open
                svg = [svg sprintf('  <g transform="translate(%d, %d)">\n', 0, y)]; %#ok<*AGROW>
                
                %axis
                svg = [svg sprintf('    <g transform="translate(%d, %d)">\n', -8, 0)]; %#ok<*AGROW>
                svg = [svg sprintf('      <line x1="%d" y1="%d" x2="%d" y2="%d" class="axis"/>\n', 0, 0, 0, h)];
                svg = [svg sprintf('      <line x1="%d" y1="%d" x2="%d" y2="%d" class="axis"/>\n', 0, 0, -tickLen, 0)];
                svg = [svg sprintf('      <line x1="%d" y1="%d" x2="%d" y2="%d" class="axis"/>\n', 0, h, -tickLen, h)];
                for j = 1:numel(states{i,2})
                    svg = [svg sprintf('      <text x="0" y="0" transform="translate(%d, %d) rotate(270)" class="axisLabel" style="text-anchor:middle; alignment-baseline:bottom">%s</text>\n', ...
                        -18 * (numel(states{i,2}) - j + 1) - 36, round(h/2), states{i, 2}{j})];
                end
                data = permute(ensemble.stateData.values(ensemble.getPropertyIndices(states{i, 1}), :, :, simIdxs), [4 3 1 2]);
                minVal = min(data(:));
                maxVal = max(data(:));
                if minVal == ceil(minVal)
                    minStr = sprintf('%d', minVal);
                else
                    minStr = sprintf('%.2f', minVal);
                end
                if maxVal == ceil(maxVal)
                    maxStr = sprintf('%d', maxVal);
                else
                    maxStr = sprintf('%.2f', maxVal);
                end
                svg = [svg sprintf('      <text x="%d" y="%d" class="tickLabel" style="text-anchor:end">%s</text>\n', -tickLen-2, h+3, minStr)];
                svg = [svg sprintf('      <text x="%d" y="%d" class="tickLabel" style="text-anchor:end">%s</text>\n', -tickLen-2, 0+3, maxStr)];
                svg = [svg sprintf('    </g>\n')];
                
                %data
                for k = numel(simIdxs):-1:1
                    path = sprintf('M %f %f ', w * this.times(1) / this.times(end), h * (1 - (data(k, 1) - minVal) / (maxVal - minVal)));
                    idxs = find(diff(data(k, :)));
                    if numel(idxs) > tAvg
                        idxs = idxs(1:floor(numel(idxs)/tAvg):end);
                    end
                    idxs = idxs(idxs + 1 <= timeIdx);
                    idxs(end+1) = find(~isnan(data(k, 1:timeIdx)), 1, 'last') - 1;
                    for j = 1:numel(idxs)
                        if idxs(j) + 1 > timeIdx
                            break;
                        end
                        path = [path sprintf('L %f %f ', w * this.times(idxs(j)+1) / this.times(end), h * (1 - (data(k, idxs(j)+1) - minVal) / (maxVal - minVal)))];
                    end
                    if k == 1
                        svg = [svg sprintf('    <path d="%s" style="stroke:darkblue; stroke-width:2px; fill:none;"/>', path)];
                    else
                        svg = [svg sprintf('    <path d="%s" style="stroke:rgb(200,200,200); stroke-width:1px; fill:none;"/>', path)];
                    end
                end
                
                %close
                svg = [svg sprintf('  </g>\n')];
            end
            
            %time axes
            svg = [svg sprintf('  <g transform="translate(%d, %d)">\n', 0, H + dy)];
            svg = [svg sprintf('    <line x1="%d" y1="%d" x2="%d" y2="%d" class="axis"/>\n', 0, 0, w, 0)];
            svg = [svg sprintf('    <text x="%d" y="%d" class="axisLabel" style="text-anchor:middle">%s</text>\n', 75, 36, 'Time (h)')];
            ticks = [0 5];
            if this.times(end) / 3600 >= 10
                ticks(end) = 10;
            end
            for i = 1:numel(ticks)
                svg = [svg sprintf('    <line x1="%d" y1="%d" x2="%d" y2="%d" class="axis"/>\n', w * ticks(i) / (this.times(end) / 3600), 0, w * ticks(i) / (this.times(end) / 3600), tickLen)];
                svg = [svg sprintf('    <text x="%d" y="%d" class="tickLabel" style="text-anchor:middle">%d</text>\n', w * ticks(i) / (this.times(end) / 3600), 18, ticks(i))];
            end
            svg = [svg sprintf('  </g>\n')];
            
            %close
            svg = [svg sprintf('</g>\n')];
            
            %% close svg
            svg = [svg sprintf('</svg>\n')];
        end
        
        function svg = drawCell(~, width, cylindricalLength, septumLength, W, H, nActive, nStalled, nFree30, nFree50, numEdges, ftsZRing1st, ftsZRing2st, ftsZRing2bt, ftsZRingRbt, zoomIn)
            %initialize
            svg = [];
            
            %start group
            svg = [svg sprintf('<g>\n')];
                     
            %rectangle clip
            svg = [svg ...
                sprintf('  <clipPath id="RectClip">\n') ...
                sprintf('    <path d="M %d,0  L %d,0  L%d,%d  L%d,%d z" style="fill:white;"/>', 0.05 * W, 0.95 * W, 0.95 * W, H, 0.05 * W, H) ...
                sprintf('  </clipPath>\n') ...
                ];
            
            %open cell and contents
            if zoomIn
                svg = [svg sprintf('  <g transform="translate(%d, %d) scale(0.95)" clip-path="url(#RectClip)">\n', 0.03 * W, -0.07 * H)];
            else
                svg = [svg sprintf('  <g transform="translate(%d, %d) scale(0.95)">\n', 0.03 * W, -0.07 * H)];
            end
            
            %background
            svg = [svg ...
                sprintf('  <path d="\n') ...
                sprintf('    M%0.4f,%0.4f a%0.4f,%0.4f 0 1 0 0,%0.4f\n', W/2 - septumLength - cylindricalLength/2, H/2 - width/2, width/2, width/2, width), ...
                sprintf('    l%0.4f,0\n', cylindricalLength/2) ...
                sprintf('    a%0.4f,%0.4f 0 0 0 %0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('    a%0.4f,%0.4f 0 0 0 %0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('    l%0.4f,0\n', cylindricalLength/2) ...
                sprintf('    a%0.4f,%0.4f 0 1 0 0,-%0.4f\n', width/2,width/2, width) ...
                sprintf('    l-%0.4f,0\n', cylindricalLength/2) ...
                sprintf('    a%0.4f,%0.4f 0 0 0 -%0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('    a%0.4f,%0.4f 0 0 0 -%0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('    z\n') ...
                '    " style="stroke:none;fill:url(#cellBg)"/>\n'];
            
            %cell clip
            svg = [svg ...
                sprintf('  <clipPath id="CellClip">\n') ...
                sprintf('    <path d="\n') ...
                sprintf('      M%0.4f,%0.4f a%0.4f,%0.4f 0 1 0 0,%0.4f\n', W/2 - septumLength - cylindricalLength/2, H/2 - width/2, width/2, width/2, width), ...
                sprintf('      l%0.4f,0\n', cylindricalLength/2) ...
                sprintf('      a%0.4f,%0.4f 0 0 0 %0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('      a%0.4f,%0.4f 0 0 0 %0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('      l%0.4f,0\n', cylindricalLength/2) ...
                sprintf('      a%0.4f,%0.4f 0 1 0 0,-%0.4f\n', width/2, width/2, width) ...
                sprintf('      l-%0.4f,0\n', cylindricalLength/2) ...
                sprintf('      a%0.4f,%0.4f 0 0 0 -%0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('      a%0.4f,%0.4f 0 0 0 -%0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('      z\n') ...
                '      " style="stroke:none;fill:url(#cellBg)"/>\n' ...
                sprintf('  </clipPath>\n') ...
                ];
            
            %ribosomes
            svg = [svg sprintf('  <g clip-path="url(#CellClip)">\n')];            
            nRibs = nActive + nStalled + nFree30 + nFree50;
            cx = 0.05 * W + 0.9 * W * rand(nRibs, 1);
            cy = 0.05 * H + 0.9 * H * rand(nRibs, 1);
            r = [5 * ones(nActive, 1); 5 * ones(nStalled, 1); 3 * ones(nFree30, 1); 3 * ones(nFree50, 1)];
            style = [
                repmat({'activeRibosome'}, nActive, 1)
                repmat({'stalledRibosome'}, nStalled, 1)
                repmat({'free30SRibosome'}, nFree30, 1)
                repmat({'free50SRibosome'}, nFree50, 1)
                ];
            order = randperm(nRibs);
            for j = 1:nRibs
                i = order(j);
                svg = [svg sprintf('    <g transform="translate(%d,%d)">\n', cx(i), cy(i))];
                svg = [svg sprintf('      <circle cx="0" cy="0" r="%d" style="fill:white;"/>\n', r(i))];
                svg = [svg sprintf('      <circle cx="0" cy="0" r="%d" class="%s"/>\n', r(i), style{i})];
                svg = [svg sprintf('    </g>\n')];
            end
            svg = [svg sprintf('  </g>\n')];
            
            %FtsZ
            if ftsZRing1st || ftsZRing2st || ftsZRing2bt || ftsZRingRbt
                ftsZRing1st = ftsZRing1st / numEdges * 2 * (septumLength + 6);
                ftsZRing2st = ftsZRing2st / numEdges * 2 * (septumLength + 6);
                ftsZRing2bt = ftsZRing2bt / numEdges * 2 * (septumLength + 6);
                ftsZRingRbt = ftsZRingRbt / numEdges * 2 * (septumLength + 6);
                svg = [svg sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" style="stroke:#ff241E; stroke-width:2px; fill:none" clip-path="url(#CellClip)"/>', ...
                    W/2 - 6, W/2 - 6, 0, H/2 - septumLength + ftsZRing1st)];
                svg = [svg sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" style="stroke:#ff241E; stroke-width:2px; fill:none" clip-path="url(#CellClip)"/>', ...
                    W/2 - 2, W/2 - 2, 0, H/2 - septumLength + ftsZRing2st)];
                svg = [svg sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" style="stroke:green; stroke-width:2px; fill:none" clip-path="url(#CellClip)"/>', ...
                    W/2 + 2, W/2 + 2, H, H/2 + septumLength - ftsZRing2bt)];
                svg = [svg sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" style="stroke:green; stroke-width:2px; fill:none" clip-path="url(#CellClip)"/>', ...
                    W/2 + 6, W/2 + 6, H, H/2 + septumLength - ftsZRingRbt)];
            end
            
            %membrane
            svg = [svg ...
                sprintf('  <path d="\n') ...
                sprintf('    M%0.4f,%0.4f a%0.4f,%0.4f 0 1 0 0,%0.4f\n', W/2 - septumLength - cylindricalLength/2, H/2 - width/2, width/2, width/2, width), ...
                sprintf('    l%0.4f,0\n', cylindricalLength/2) ...
                sprintf('    a%0.4f,%0.4f 0 0 0 %0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('    a%0.4f,%0.4f 0 0 0 %0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('    l%0.4f,0\n', cylindricalLength/2) ...
                sprintf('    a%0.4f,%0.4f 0 1 0 0,-%0.4f\n', width/2,width/2, width) ...
                sprintf('    l-%0.4f,0\n', cylindricalLength/2) ...
                sprintf('    a%0.4f,%0.4f 0 0 0 -%0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('    a%0.4f,%0.4f 0 0 0 -%0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('    z\n') ...
                '    " style="stroke:darkblue;stroke-width:2;stroke-linecap:round;stroke-linejoin:miter;fill:none;"/>\n'];
            
            %close cell and contents
            svg = [svg sprintf('  </g>\n')];
            
            %legend
            svg = [svg sprintf('  <g transform="translate(6, 5)">\n')];
            svg = [svg sprintf('    <text x="%d" y="%d" class="%s">%s</text>\n', 33, H, 'axisLabel', 'Ribosomes')];
            svg = [svg sprintf('    <rect x="%d" width="%d" y="%d" height="%d" class="%s"/><text x="%d" y="%d" class="%s">%s</text>\n', 0, 10, 16 + H-5,    2, 'activeRibosome',  13, 16 + H,    'legendLabel', 'Active')];
            svg = [svg sprintf('    <rect x="%d" width="%d" y="%d" height="%d" class="%s"/><text x="%d" y="%d" class="%s">%s</text>\n', 0, 10, 16 + H-5+16, 2, 'stalledRibosome', 13, 16 + H+16, 'legendLabel', 'Stalled')];
            svg = [svg sprintf('    <rect x="%d" width="%d" y="%d" height="%d" class="%s"/><text x="%d" y="%d" class="%s">%s</text>\n', 68, 10, 16 + H-5,    2, 'free30SRibosome', 81, 16 + H,    'legendLabel', 'Free 30S')];
            svg = [svg sprintf('    <rect x="%d" width="%d" y="%d" height="%d" class="%s"/><text x="%d" y="%d" class="%s">%s</text>\n', 68, 10, 16 + H-5+16, 2, 'free50SRibosome', 81, 16 + H+16, 'legendLabel', 'Free 50S')];
            
            svg = [svg sprintf('    <text x="%d" y="%d" class="%s">%s</text>\n', 168, H, 'axisLabel', 'FtsZ')];
            svg = [svg sprintf('    <rect x="%d" width="%d" y="%d" height="%d" style="fill:%s"/><text x="%d" y="%d" class="%s">%s</text>\n', 153, 10, 16 + H-5,    2, '#ff241E', 166, 16 + H,    'legendLabel', 'Straight')];
            svg = [svg sprintf('    <rect x="%d" width="%d" y="%d" height="%d" style="fill:%s"/><text x="%d" y="%d" class="%s">%s</text>\n', 153, 10, 16 + H-5+16, 2, 'green',   166, 16 + H+16, 'legendLabel', 'Bent')];
            svg = [svg sprintf('  </g>\n')];
            
            %end group
            svg = [svg sprintf('</g>\n')];
        end
        
        function svg = drawCellDivision(~, width, numEdges, septumLength, W, H, ftsZRing1st, ftsZRing2st, ftsZRing2bt, ftsZRingRbt)
            svg = [];
            
            %background
            svg = [svg ...
                sprintf('<path d="\n') ...
                sprintf('  M%0.4f,%0.4f l0,%0.4f\n', W/2 - septumLength, H/2 - width/2, width), ...
                sprintf('  a%0.4f,%0.4f 0 0 0 %0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('  a%0.4f,%0.4f 0 0 0 %0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('  l0,%0.4f\n', -width) ...
                sprintf('  a%0.4f,%0.4f 0 0 0 -%0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('  a%0.4f,%0.4f 0 0 0 -%0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                '  " style="stroke:none;fill:url(#cellBg)"/>\n'];
            
            %FtsZ
            svg = [svg ...
                sprintf('<clipPath id="FtsZClip"><path d="\n') ...
                sprintf('  M%0.4f,%0.4f l0,%0.4f\n', W/2 - septumLength, H/2 - width/2, width), ...
                sprintf('  a%0.4f,%0.4f 0 0 0 %0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('  a%0.4f,%0.4f 0 0 0 %0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('  l0,%0.4f\n', -width) ...
                sprintf('  a%0.4f,%0.4f 0 0 0 -%0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('  a%0.4f,%0.4f 0 0 0 -%0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                '  " style="stroke:darkblue;stroke-width:2;stroke-linecap:round;stroke-linejoin:miter;fill:none"/></clipPath>\n'];
            
            ftsZRing1st = ftsZRing1st / numEdges * 2 * (septumLength + 6);
            ftsZRing2st = ftsZRing2st / numEdges * 2 * (septumLength + 6);
            ftsZRing2bt = ftsZRing2bt / numEdges * 2 * (septumLength + 6);
            ftsZRingRbt = ftsZRingRbt / numEdges * 2 * (septumLength + 6);
            svg = [svg sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" style="stroke:#ff241E; stroke-width:2px; fill:none" clip-path="url(#FtsZClip)"/>', ...
                W/2 - 6, W/2 - 6, 0, H/2 - septumLength + ftsZRing1st)];
            svg = [svg sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" style="stroke:#ff241E; stroke-width:2px; fill:none" clip-path="url(#FtsZClip)"/>', ...
                W/2 - 2, W/2 - 2, 0, H/2 - septumLength + ftsZRing2st)];
            svg = [svg sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" style="stroke:green; stroke-width:2px; fill:none" clip-path="url(#FtsZClip)"/>', ...
                W/2 + 2, W/2 + 2, H, H/2 + septumLength - ftsZRing2bt)];
            svg = [svg sprintf('<line x1="%d" x2="%d" y1="%d" y2="%d" style="stroke:green; stroke-width:2px; fill:none" clip-path="url(#FtsZClip)"/>', ...
                W/2 + 6, W/2 + 6, H, H/2 + septumLength - ftsZRingRbt)];
            
            %outline
            svg = [svg ...
                sprintf('<path d="\n') ...
                sprintf('  M%0.4f,%0.4f m0,%0.4f\n', W/2 - septumLength, H/2 - width/2, width), ...
                sprintf('  a%0.4f,%0.4f 0 0 0 %0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('  a%0.4f,%0.4f 0 0 0 %0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('  m0,%0.4f\n', -width) ...
                sprintf('  a%0.4f,%0.4f 0 0 0 -%0.4f,%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                sprintf('  a%0.4f,%0.4f 0 0 0 -%0.4f,-%0.4f\n', septumLength, septumLength, septumLength, septumLength) ...
                '  " style="stroke:rgb(0,0,0);stroke-width:2;stroke-linecap:round;stroke-linejoin:miter;fill:none"/>\n'];
        end
        
        function svg = drawArc(~, posStrnds, lens, cx, cy, r0, dr, strokeColor, strokeWidth, twoChrs, sequenceLen)
            svg = [];
            for j = 1:size(posStrnds, 1)
                theta = posStrnds(j, 1) / sequenceLen * 2 * pi - pi/2;
                dTheta = lens(j) / sequenceLen * 2 * pi;
                largeArc = dTheta >= pi;
                tmpcx = cx;
                if ~twoChrs
                    r = r0 + (2 * iseven(posStrnds(j, 2)) - 1) * dr;
                else
                    r = r0/2 + (2 * iseven(posStrnds(j, 2)) - 1) * dr;
                    if posStrnds(j, 2) > 2
                        tmpcx = 3/2 * cx;
                    else
                        tmpcx = cx / 2;
                    end
                end
                if lens(j) == sequenceLen
                    svg = [svg sprintf('<circle cx="%d" cy="%d" r="%d" stroke="%s" stroke-width="%f" fill="none"/>\n', ...
                        tmpcx, cy, r, ...
                        strokeColor, strokeWidth)]; %#ok<*AGROW>
                else
                    svg = [svg sprintf('<path d="M%f,%f A%f,%f %d %d,%d %f %f" stroke="%s" stroke-width="%d" fill="none"/>\n', ...
                        tmpcx + r*cos(theta), cy + r*sin(theta), ...
                        r, r, 0, largeArc, 1, tmpcx + r*cos(theta+dTheta), cy + r*sin(theta+dTheta), ...
                        strokeColor, strokeWidth)];
                end
            end
        end
    end
end
