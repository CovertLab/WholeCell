%MacromoleculeUtil
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Derek Macklin, macklin@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 6/28/2011
classdef MacromoleculeUtil
    methods (Static = true)
        function [idx, compIdx, id, name, type, ...
                geneIdxs, geneCompIdxs, geneIds, geneNames, ...
                nascentRnaIdxs, nascentRnaCompIdxs, nascentRnaIds, nascentRnaNames, ...
                matureRnaIdxs, matureRnaCompIdxs, matureRnaIds, matureRnaNames, ...
                monomerIdxs, monomerCompIdxs, monomerIds, monomerNames, ...
                complexIdxs, complexCompIdxs, complexIds, complexNames] = ...
                getMacromoleculeIndexsIDsNames(id, sim)
            
            validateattributes(id, {'char'}, {'row'})
            
            comp = sim.compartment;
            g = sim.gene;
            rna = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            geneIdxs = g.getIndexs(id);
            [~, nascentRnaIdxs] = ismember(id, rna.wholeCellModelIDs(rna.nascentIndexs));
            matureRnaIdxs = rna.getIndexs(id);
            monomerIdxs = pm.getIndexs(id);
            complexIdxs = pc.getIndexs(id);
            
            if geneIdxs
                type = 'gene';
                idx = geneIdxs;
                compIdx = comp.cytosolIndexs;
                name = g.names{geneIdxs};
                
                nascentRnaIdxs = find(any(rna.nascentRNAGeneComposition(geneIdxs, :), 1));
                matureRnaIdxs = find(any(rna.matureRNAGeneComposition(geneIdxs, :), 1));
                monomerIdxs = find(g.mRNAIndexs == geneIdxs);
                complexIdxs = find(any(any(pc.proteinComplexComposition(g.mRNAIndexs(monomerIdxs), :, :), 3), 1));
            elseif nascentRnaIdxs
                type = 'nascentRna';
                idx = nascentRnaIdxs;
                compIdx = comp.cytosolIndexs;
                name = rna.names{rna.nascentIndexs(nascentRnaIdxs)};
                
                geneIdxs = find(any(rna.nascentRNAGeneComposition(:, nascentRnaIdxs), 2));
                matureRnaIdxs = find(any(rna.nascentRNAMatureRNAComposition(:, nascentRnaIdxs), 2));
                monomerIdxs = find(any(rna.matureRNAGeneComposition(g.mRNAIndexs, matureRnaIdxs), 2));
                complexIdxs = find(any(any(pc.proteinComplexComposition(g.mRNAIndexs(monomerIdxs), :, :), 3), 1));
            elseif matureRnaIdxs
                type = 'matureRna';
                idx = matureRnaIdxs;
                compIdx = comp.cytosolIndexs;
                name = rna.names{rna.matureIndexs(matureRnaIdxs)};
                
                geneIdxs = find(any(rna.matureRNAGeneComposition(:, matureRnaIdxs), 2));
                nascentRnaIdxs = find(any(rna.nascentRNAMatureRNAComposition(matureRnaIdxs, :), 1));
                monomerIdxs = find(any(rna.matureRNAGeneComposition(g.mRNAIndexs, matureRnaIdxs), 2));
                complexIdxs = find(any(any(pc.proteinComplexComposition(g.mRNAIndexs(monomerIdxs), :, :), 3), 1));
            elseif monomerIdxs
                type = 'monomer';
                idx = monomerIdxs;
                compIdx = pm.compartments(pm.matureIndexs(monomerIdxs));
                name = pm.names{pm.matureIndexs(monomerIdxs)};
                
                geneIdxs = g.mRNAIndexs(monomerIdxs);
                matureRnaIdxs = find(any(rna.matureRNAGeneComposition(g.mRNAIndexs(monomerIdxs), :), 1));
                nascentRnaIdxs = find(any(rna.nascentRNAGeneComposition(g.mRNAIndexs(monomerIdxs), :), 1));
                complexIdxs = find(any(any(pc.proteinComplexComposition(g.mRNAIndexs(monomerIdxs), :, :), 3), 1));
            elseif complexIdxs
                type = 'complex';
                idx = complexIdxs;
                compIdx = pc.compartments(pc.matureIndexs(complexIdxs));
                name = pc.names{pc.matureIndexs(complexIdxs)};
                
                monomerIdxs = find(any(any(pc.proteinComplexComposition(g.mRNAIndexs, complexIdxs, :), 3), 2));
                geneIdxs = g.mRNAIndexs(monomerIdxs);
                matureRnaIdxs = find(any(rna.matureRNAGeneComposition(g.mRNAIndexs(monomerIdxs), :), 1));
                nascentRnaIdxs = find(any(rna.nascentRNAGeneComposition(g.mRNAIndexs(monomerIdxs), :), 1));
            else
                idx = []; compIdx = []; id = []; name = []; type = [];
                geneIdxs = []; geneCompIdxs = []; geneIds = []; geneNames = [];
                nascentRnaIdxs = []; nascentRnaCompIdxs = []; nascentRnaIds = []; nascentRnaNames = [];
                matureRnaIdxs = []; matureRnaCompIdxs = []; matureRnaIds = []; matureRnaNames = [];
                monomerIdxs = []; monomerCompIdxs = []; monomerIds = []; monomerNames = [];
                complexIdxs = []; complexCompIdxs = []; complexIds = []; complexNames = [];
                return;
            end
            
            geneCompIdxs = comp.cytosolIndexs(ones(size(geneIdxs)));
            nascentRnaCompIdxs = comp.cytosolIndexs(ones(size(nascentRnaIdxs)));
            matureRnaCompIdxs = comp.cytosolIndexs(ones(size(nascentRnaIdxs)));
            monomerCompIdxs = pm.compartments(pm.matureIndexs(monomerIdxs));
            complexCompIdxs = pc.compartments(pc.matureIndexs(complexIdxs));
            
            geneIds = g.wholeCellModelIDs(geneIdxs);
            nascentRnaIds = rna.wholeCellModelIDs(rna.nascentIndexs(nascentRnaIdxs));
            matureRnaIds = rna.wholeCellModelIDs(rna.matureIndexs(matureRnaIdxs));
            monomerIds = pm.wholeCellModelIDs(pm.matureIndexs(monomerIdxs));
            complexIds = pc.wholeCellModelIDs(pc.matureIndexs(complexIdxs));
            
            geneNames = g.names(geneIdxs);
            nascentRnaNames = rna.names(rna.nascentIndexs(nascentRnaIdxs));
            matureRnaNames = rna.names(rna.matureIndexs(matureRnaIdxs));
            monomerNames = pm.names(pm.matureIndexs(monomerIdxs));
            complexNames = pc.names(pc.matureIndexs(complexIdxs));
        end
        
        function [cnts, ...
                geneCnts, ...
                nascentRnaCnts, processedRnaCnts, matureRnaCnts, boundRnaCnts, misfoldedRnaCnts, damagedRnaCnts, aminoacylatedRnaCnts, ...
                nascentMonCnts, processedIMonCnts, translocatedMonCnts, processedIIMonCnts, foldedMonCnts, matureMonCnts, inactivatedMonCnts, boundMonCnts, misfoldedMonCnts, damagedMonCnts, ...
                nascentCpxCnts, matureCpxCnts, inactivatedCpxCnts, boundCpxCnts, misfoldedCpxCnts, damagedCpxCnts, ...
                rnaPolTUs, rnaPolPosStrnds, ribMRNAs, ribPos] = ...
                getMacromoleculeCounts(id, sim, simDir)
            import edu.stanford.covert.cell.sim.state.RNAPolymerase;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.MacromoleculeUtil;
            
            [idx, compIdx, ~, ~, type, ...
                geneIdxs, ~, ~, ~, ...
                nascentRnaIdxs, nascentRnaCompIdxs, ~, ~, ...
                matureRnaIdxs, matureRnaCompIdxs, ~, ~, ...
                monIdxs, monCompIdxs, ~, ~, ...
                cpxIdxs, cpxCompIdxs, ~, ~] = ...
                MacromoleculeUtil.getMacromoleculeIndexsIDsNames(id, sim);
            
            g = sim.gene;
            comp = sim.compartment;
            chr = sim.state('Chromosome');
            rna = sim.state('Rna');
            mon = sim.state('ProteinMonomer');
            cpx = sim.state('ProteinComplex');
            
            stateNames = {
                'Time'               'values'
                'Chromosome'         'polymerizedRegions'
                'Rna'                'counts'
                'ProteinMonomer'     'counts'
                'ProteinComplex'     'counts'
                };
            states = DiskLogger.load(simDir, stateNames, [], [], [], 'extract');
            times = permute(states.Time.values, [1 3 2]);
            polRgns = states.Chromosome.polymerizedRegions;
            rnaCnts = states.Rna.counts;
            monCnts = states.ProteinMonomer.counts;
            cpxCnts = states.ProteinComplex.counts;
            
            switch type
                case 'gene'
                    cnts = MacromoleculeUtil.extractCopyNumberTimeCourse(chr, polRgns, g.startCoordinates(idx), g.strands(idx), g.lengths(idx));
                case 'nascentRna'
                    cnts = MacromoleculeUtil.extractCountTimeCourse(rnaCnts, rna.processedIndexs(idx), compIdx);
                case 'matureRna'
                    cnts = MacromoleculeUtil.extractCountTimeCourse(rnaCnts, rna.matureIndexs(idx), compIdx);
                case 'monomer'
                    cnts = MacromoleculeUtil.extractCountTimeCourse(monCnts, mon.matureIndexs(idx), compIdx);
                case 'complex'
                    cnts = MacromoleculeUtil.extractCountTimeCourse(cpxCnts, cpx.matureIndexs(idx), compIdx);
            end
            
            if nargout == 1
                return;
            end
            
            geneCnts = MacromoleculeUtil.extractCopyNumberTimeCourse(chr, polRgns, g.startCoordinates(geneIdxs), g.strands(geneIdxs), g.lengths(geneIdxs));
            
            nascentRnaCnts       = MacromoleculeUtil.extractCountTimeCourse(rnaCnts, rna.nascentIndexs(nascentRnaIdxs),      nascentRnaCompIdxs);
            processedRnaCnts     = MacromoleculeUtil.extractCountTimeCourse(rnaCnts, rna.processedIndexs(matureRnaIdxs),     matureRnaCompIdxs);
            matureRnaCnts        = MacromoleculeUtil.extractCountTimeCourse(rnaCnts, rna.matureIndexs(matureRnaIdxs),        matureRnaCompIdxs);
            boundRnaCnts         = MacromoleculeUtil.extractCountTimeCourse(rnaCnts, rna.boundIndexs(matureRnaIdxs),         matureRnaCompIdxs);
            misfoldedRnaCnts     = MacromoleculeUtil.extractCountTimeCourse(rnaCnts, rna.misfoldedIndexs(matureRnaIdxs),     matureRnaCompIdxs);
            damagedRnaCnts       = MacromoleculeUtil.extractCountTimeCourse(rnaCnts, rna.damagedIndexs(matureRnaIdxs),       matureRnaCompIdxs);
            aminoacylatedRnaCnts = MacromoleculeUtil.extractCountTimeCourse(rnaCnts, rna.aminoacylatedIndexs(matureRnaIdxs), matureRnaCompIdxs);
            
            nascentMonCnts       = MacromoleculeUtil.extractCountTimeCourse(monCnts, mon.nascentIndexs(monIdxs),     comp.cytosolIndexs(ones(size(monIdxs))));
            processedIMonCnts    = MacromoleculeUtil.extractCountTimeCourse(monCnts, mon.processedIIndexs(monIdxs),  comp.cytosolIndexs(ones(size(monIdxs))));
            translocatedMonCnts  = MacromoleculeUtil.extractCountTimeCourse(monCnts, mon.processedIIndexs(monIdxs),  monCompIdxs);
            processedIIMonCnts   = MacromoleculeUtil.extractCountTimeCourse(monCnts, mon.processedIIIndexs(monIdxs), monCompIdxs);
            foldedMonCnts        = MacromoleculeUtil.extractCountTimeCourse(monCnts, mon.foldedIndexs(monIdxs),      monCompIdxs);
            matureMonCnts        = MacromoleculeUtil.extractCountTimeCourse(monCnts, mon.matureIndexs(monIdxs),      monCompIdxs);
            inactivatedMonCnts   = MacromoleculeUtil.extractCountTimeCourse(monCnts, mon.inactivatedIndexs(monIdxs), monCompIdxs);
            boundMonCnts         = MacromoleculeUtil.extractCountTimeCourse(monCnts, mon.boundIndexs(monIdxs),       monCompIdxs);
            misfoldedMonCnts     = MacromoleculeUtil.extractCountTimeCourse(monCnts, mon.misfoldedIndexs(monIdxs),   monCompIdxs);
            damagedMonCnts       = MacromoleculeUtil.extractCountTimeCourse(monCnts, mon.damagedIndexs(monIdxs),     monCompIdxs);
            
            nascentCpxCnts       = MacromoleculeUtil.extractCountTimeCourse(cpxCnts, cpx.nascentIndexs(cpxIdxs),     cpxCompIdxs);
            matureCpxCnts        = MacromoleculeUtil.extractCountTimeCourse(cpxCnts, cpx.matureIndexs(cpxIdxs),      cpxCompIdxs);
            inactivatedCpxCnts   = MacromoleculeUtil.extractCountTimeCourse(cpxCnts, cpx.inactivatedIndexs(cpxIdxs), cpxCompIdxs);
            boundCpxCnts         = MacromoleculeUtil.extractCountTimeCourse(cpxCnts, cpx.boundIndexs(cpxIdxs),       cpxCompIdxs);
            misfoldedCpxCnts     = MacromoleculeUtil.extractCountTimeCourse(cpxCnts, cpx.misfoldedIndexs(cpxIdxs),   cpxCompIdxs);
            damagedCpxCnts       = MacromoleculeUtil.extractCountTimeCourse(cpxCnts, cpx.damagedIndexs(cpxIdxs),     cpxCompIdxs);
            
            if nargout >= 26
                stateNames = {
                    'RNAPolymerase'      'states'
                    'RNAPolymerase'      'positionStrands'
                    'Transcript'         'boundTranscriptionUnits'
                    'Transcript'         'boundTranscriptProgress'
                    'Transcript'         'boundTranscriptChromosome'
                    };
                states = DiskLogger.load(simDir, stateNames, [], [], [], 'extract');                                
                polStates = permute(states.RNAPolymerase.states, [1 3 2]);
                pos = permute(states.RNAPolymerase.positionStrands(:, 1, :), [1 3 2]);
                strnds = permute(states.RNAPolymerase.positionStrands(:, 2, :), [1 3 2]);
                tus = permute(states.Transcript.boundTranscriptionUnits, [1 3 2]);
                lens = permute(states.Transcript.boundTranscriptProgress, [1 3 2]);
                chrs = permute(states.Transcript.boundTranscriptChromosome, [1 3 2]);
                
                rnaPolTUs = zeros(0, 3);
                rnaPolPosStrnds = cell(0, 1);
                for i = 1:numel(nascentRnaIdxs)
                    for j = 1:size(tus, 1)
                        for k = 1:2
                            startIdxs = find(...
                                tus(j, :) == nascentRnaIdxs(i) & ...
                                chrs(j, :) == k & ...
                                polStates(j, :) == RNAPolymerase.specificallyBoundValue & ...
                                [true ...
                                tus(j, 1:end-1)~=nascentRnaIdxs(i) | ...
                                chrs(j, 1:end-1)~=k | ...
                                polStates(j, 1:end-1) ~= RNAPolymerase.specificallyBoundValue...
                                ]);
                            endIdxs = find(...
                                tus(j, :) == nascentRnaIdxs(i) & ...
                                chrs(j, :) == k & ...
                                [tus(j, 2:end) ~= nascentRnaIdxs(i) | ...
                                chrs(j, 2:end) ~= k ...
                                true]);
                            
                            for l = 1:numel(startIdxs)
                                rnaPolTUs = [rnaPolTUs; nascentRnaIdxs(i) k strnds(j, startIdxs(l))]; %#ok<AGROW>
                                rnaPolPosStrnds = [rnaPolPosStrnds;
                                    {[times(startIdxs(l):endIdxs(l)) 
                                    lens(j, startIdxs(l):endIdxs(l))
                                    pos(j, startIdxs(l):endIdxs(l))]}]; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
            
            if nargout >= 28
                stateNames = {
                    'Polypeptide'        'boundMRNAs'
                    'Polypeptide'        'nascentMonomerLengths'
                    };
                states = DiskLogger.load(simDir, stateNames, [], [], [], 'extract');
                mrnas = permute(states.Polypeptide.boundMRNAs, [1 3 2]);
                lens = permute(states.Polypeptide.nascentMonomerLengths, [1 3 2]);
                
                ribMRNAs = zeros(0, 1);
                ribPos = cell(0, 1);
                for i = 1:numel(monIdxs)
                    for j = 1:size(mrnas, 1)
                        startIdxs = find(...
                            mrnas(j, :) == monIdxs(i) & ...
                            [true mrnas(j, 1:end-1) ~= monIdxs(i)]);
                        endIdxs = find(...
                            mrnas(j, :) == monIdxs(i) & ...
                            [mrnas(j, 2:end) ~= monIdxs(i) true]);
                        
                        for l = 1:numel(startIdxs)
                            ribMRNAs = [ribMRNAs; monIdxs(i)]; %#ok<AGROW>
                            ribPos = [ribPos;
                                {[times(startIdxs(l):endIdxs(l)); lens(j, startIdxs(l):endIdxs(l))]}]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
    
    methods (Static = true)
        function value = extractCopyNumberTimeCourse(chr, polRgns, startCoordinates, strands, lengths)
            numTimePoints = size(polRgns, 3);
            times = 1:numTimePoints;
            
            posStrndTimes = [
                reshape(startCoordinates(:, ones(numTimePoints, 1)), [], 1)  reshape(strands(:, ones(numTimePoints, 1)),   [], 1)  reshape(times(ones(size(strands)), :), [], 1)
                reshape(startCoordinates(:, ones(numTimePoints, 1)), [], 1)  reshape(strands(:, ones(numTimePoints, 1))+2, [], 1)  reshape(times(ones(size(strands)), :), [], 1)
                ];
            lens = [
                reshape(lengths(:, ones(numTimePoints, 1)), [], 1)
                reshape(lengths(:, ones(numTimePoints, 1)), [], 1)
                ];
            value = sum(reshape(reshape(chr.isRegionPolymerized(posStrndTimes, lens, false, false, false, polRgns), [], 2), [numel(strands) numTimePoints 2]), 3);
        end
        
        function value = extractCountTimeCourse(timeCourse, idxs, compIdxs)
            numTimePoints = size(timeCourse, 3);
            times = 1:numTimePoints;
            compIdxs = compIdxs(:);
            
            if isa(timeCourse, 'edu.stanford.covert.util.SparseMat')
                value = reshape(timeCourse([...
                    reshape(idxs(:, ones(numTimePoints, 1)), [], 1) ...
                    reshape(compIdxs(:, ones(numTimePoints, 1)), [], 1) ...
                    reshape(times(ones(size(idxs)), :), [], 1) ...
                    ]), [], numTimePoints);
            else
                value = timeCourse(sub2ind(...
                    size(timeCourse), ...
                    idxs(:, ones(numTimePoints, 1)), ...
                    compIdxs(:, ones(numTimePoints, 1)), ...
                    times(ones(size(idxs)), :) ...
                    ));
            end
        end
    end
end