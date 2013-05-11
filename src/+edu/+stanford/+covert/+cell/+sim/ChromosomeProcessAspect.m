%ChromosomeProcessAspect
% Aspect class which provides processes with a handle to the chromosome state, and
% which provides some convenience methods for the process to use to change or
% query chromosome state.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/9/2010
classdef ChromosomeProcessAspect < handle
    %constants
    properties
        enzymeDNAFootprints         %DNA footprints of enzymes (nt)
        enzymeDNAFootprints3Prime   %DNA footprints of enzymes, 3' to center base (nt)
        enzymeDNAFootprints5Prime   %DNA footprints of enzymes, 5' to center base (nt)
    end
    
    %global state (stored locally for convenience)
    properties
        chromosome %instance of Chromosome
    end
    
    methods
        %retrieve references to state objects from simulation
        function storeObjectReferences(this, simulation)
            this.chromosome = simulation.state('Chromosome');
            this.states = [this.states; {this.chromosome}];
        end
    end
    
    methods
        function initializeConstants(this, ~, ~)
            c = this.chromosome;
            
            %DNA footprints
            this.enzymeDNAFootprints = zeros(size(this.enzymeWholeCellModelIDs));
            this.enzymeDNAFootprints(this.enzymeMonomerLocalIndexs) = c.monomerDNAFootprints(this.enzymeMonomerGlobalIndexs);
            this.enzymeDNAFootprints(this.enzymeComplexLocalIndexs) = c.complexDNAFootprints(this.enzymeComplexGlobalIndexs);
            [this.enzymeDNAFootprints3Prime, this.enzymeDNAFootprints5Prime] = c.calculateFootprintOverhangs(this.enzymeDNAFootprints);
        end
    end
    
    methods
        %this.enzymeMonomerLocalIndexs must be sorted
        function [nBound, posStrnds] = bindProteinToChromosomeStochastically(this,...
                enzymeIndex, nProteins, positionsStrands, lengths,...
                calcRegionWeightsFun, calcBindingPositionFun, calcNewRegionsFun, ...
                checkRegionSupercoiled)
            %% initialize output
            nBound = 0;
            
            %% process arguments, supply optional arguments
            if nargin < 3 || isempty(nProteins)
                nProteins = this.enzymes(enzymeIndex);
            end
            if nProteins == 0
                posStrnds = zeros(0,2);
                return;
            end
            
            c = this.chromosome;
            footprint = this.enzymeDNAFootprints(enzymeIndex);
            if nargin < 4 || isscalar(positionsStrands) && isnan(positionsStrands)
                [positionsStrands, lengths] = find(c.polymerizedRegions);
            end
            if nargin < 6 || isempty(calcRegionWeightsFun)
                calcRegionWeightsFun = @calcRegionWeights;
            end
            if nargin < 7 || isempty(calcBindingPositionFun)
                calcBindingPositionFun = @calcBindingPosition;
            end
            if nargin < 8 || isempty(calcNewRegionsFun)
                calcNewRegionsFun = @calcNewRegions;
            end
            if nargin < 9
                checkRegionSupercoiled = false;
            end
            
            %% find accessible regions
            tf = ismembc(enzymeIndex, this.enzymeMonomerLocalIndexs);
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions(...
                this.enzymeGlobalIndexs(enzymeIndex( tf, 1), 1), ...
                this.enzymeGlobalIndexs(enzymeIndex(~tf, 1), 1));
            [rgnPosStrnds, rgnLens] = c.intersectRegions(...
                rgnPosStrnds, rgnLens, positionsStrands, lengths);
            
            %% compute probability of binding each region
            rgnProbs = calcRegionWeightsFun(rgnLens);
            
            %% randomly select regions to bind
            posStrnds = zeros(nProteins, 2);
            nBound = 0;
            for i = 1:nProteins
                if ~any(rgnProbs); break; end
                
                %pick a region to bind
                rgnIdx = this.randStream.randsample(numel(rgnProbs), 1, true, rgnProbs);
                
                %pick a position within region to bind
                offset = calcBindingPositionFun(this.randStream.rand, rgnLens(rgnIdx)) - 1;
                
                %store selected position and strand
                posStrnds(i, :) = rgnPosStrnds(rgnIdx,:) + [offset 0];
                
                %split region about new protein position
                [rgnPosStrnds, rgnLens, rgnProbs] = calcNewRegionsFun(rgnPosStrnds, rgnLens, rgnProbs, rgnIdx, offset);
                
                nBound = nBound + 1;
            end
            
            posStrnds = posStrnds(1:nBound, :);
            if nargin < 4
                posStrnds(:, 1) = mod(posStrnds(:, 1) - 1, c.sequenceLen) + 1;
            end
            
            if nBound == 0
                return;
            end
            
            %% bind to positions
            tmpMonomerBoundSites = c.monomerBoundSites;
            tmpComplexBoundSites = c.complexBoundSites;
            if ~all(this.bindProteinToChromosome(posStrnds, enzymeIndex, nBound, [], [], false, 1, false, [], checkRegionSupercoiled))
                if isscalar(lengths)
                    lengths = lengths(ones(size(positionsStrands, 1), 1));
                end
                
                msg = [];
                msg = [msg sprintf('Region should be accessible\n')];
                
                msg = [msg sprintf('\tEnzyme: %s, index = %d, nProteins = %d, nBound = %d\n\n', this.enzymeWholeCellModelIDs{enzymeIndex}, enzymeIndex, nProteins, nBound)];
                
                msg = [msg sprintf('\tpositionStrands\n')];
                msg = [msg sprintf('\t%8s %6s %6s\n', 'Position', 'Strand', 'Length')];
                msg = [msg sprintf('\t%8s %6s %6s\n', '========', '======', '======')];
                for i = 1:size(positionsStrands, 1)
                    msg = [msg sprintf('\t%8d %6d %6d\n', positionsStrands(i, 1), positionsStrands(i, 2), lengths(i))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                
                msg = [msg sprintf('\tposStrnds\n')];
                msg = [msg sprintf('\t%8s %6s\n', 'Position', 'Strand')];
                msg = [msg sprintf('\t%8s %6s\n', '========', '======')];
                for i = 1:size(posStrnds, 1)
                    msg = [msg sprintf('\t%8d %6d\n', posStrnds(i, 1), posStrnds(i, 2))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                
                msg = [msg sprintf('\tRegion Weights Function\n\t%s\n\n', func2str(calcRegionWeightsFun))];
                msg = [msg sprintf('\tBinding Position Function\n\t%s\n\n', func2str(calcBindingPositionFun))];
                msg = [msg sprintf('\tNew Region Function\n\t%s\n\n', func2str(calcNewRegionsFun))];
                msg = [msg sprintf('\tCheck region supercoiled: %d\n\n', checkRegionSupercoiled)];
                
                msg = [msg sprintf('\tpolymerizedRegions\n')];
                msg = [msg sprintf('\t%8s %6s %6s\n', 'Position', 'Strand', 'Length')];
                msg = [msg sprintf('\t%8s %6s %6s\n', '========', '======', '======')];
                [subs, vals] = find(c.polymerizedRegions);
                for i = 1:size(subs, 1)
                    msg = [msg sprintf('\t%8d %6d %6d\n', subs(i, 1), subs(i, 2), vals(i))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                
                msg = [msg sprintf('\tlinkingNumbers\n')];
                msg = [msg sprintf('\t%8s %6s %14s\n', 'Position', 'Strand', 'Linking Number')];
                msg = [msg sprintf('\t%8s %6s %14s\n', '========', '======', '==============')];
                [subs, vals] = find(c.linkingNumbers);
                for i = 1:size(subs, 1)
                    msg = [msg sprintf('\t%8d %6d %14.4f\n', subs(i, 1), subs(i, 2), vals(i))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                
                msg = [msg sprintf('\tmonomerBoundSites\n')];
                msg = [msg sprintf('\t%8s %6s %7s\n', 'Position', 'Strand', 'Monomer')];
                msg = [msg sprintf('\t%8s %6s %7s\n', '========', '======', '=======')];
                [subs, vals] = find(tmpMonomerBoundSites);
                for i = 1:size(subs, 1)
                    msg = [msg sprintf('\t%8d %6d %6d\n', subs(i, 1), subs(i, 2), vals(i))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                
                msg = [msg sprintf('\tcomplexBoundSites\n')];
                msg = [msg sprintf('\t%8s %6s %7s\n', 'Position', 'Strand', 'Complex')];
                msg = [msg sprintf('\t%8s %6s %7s\n', '========', '======', '=======')];
                [subs, vals] = find(tmpComplexBoundSites);
                for i = 1:size(subs, 1)
                    msg = [msg sprintf('\t%8d %6d %6d\n', subs(i, 1), subs(i, 2), vals(i))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                
                msg = [msg sprintf('\tdamagedSites\n')];
                msg = [msg sprintf('\t%8s %6s %6s\n', 'Position', 'Strand', 'Damage')];
                msg = [msg sprintf('\t%8s %6s %6s\n', '========', '======', '======')];
                [subs, vals] = find(c.damagedSites);
                for i = 1:size(subs, 1)
                    msg = [msg sprintf('\t%8d %6d %6d\n', subs(i, 1), subs(i, 2), vals(i))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                
                [subs, vals] = find(c.monomerBoundSites - tmpMonomerBoundSites);
                msg = [msg sprintf('\tBinding monomers\n')];
                msg = [msg sprintf('\t%8s %6s %7s\n', 'Position', 'Strand', 'Monomer')];
                msg = [msg sprintf('\t%8s %6s %7s\n', '========', '======', '=======')];
                for i = 1:size(subs, 1)
                    if vals(i) < 0
                        break;
                    end
                    msg = [msg sprintf('\t%8d %6d %7d\n', subs(i, 1), subs(i, 2), vals(i))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                if any(this.enzymeMonomerGlobalIndexs == enzymeIndex)
                    subs2 = posStrnds(c.monomerBoundSites(posStrnds) ~= this.enzymeGlobalIndexs(enzymeIndex), :);
                else
                    subs2 = zeros(0, 2);
                end
                msg = [msg sprintf('\tNot binding monomers\n')];
                msg = [msg sprintf('\t%8s %6s\n', 'Position', 'Strand')];
                msg = [msg sprintf('\t%8s %6s\n', '========', '======')];
                for i = 1:size(subs2, 1)
                    msg = [msg sprintf('\t%8d %6d\n', subs2(i, 1), subs2(i, 2))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                msg = [msg sprintf('\tUnbinding monomers\n')];
                msg = [msg sprintf('\t%8s %6s %7s\n', 'Position', 'Strand', 'Monomer')];
                msg = [msg sprintf('\t%8s %6s %7s\n', '========', '======', '=======')];
                for i = 1:size(subs, 1)
                    if vals(i) > 0
                        break;
                    end
                    msg = [msg sprintf('\t%8d %6d %7d\n', subs(i, 1), subs(i, 2), -vals(i))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                
                [subs, vals] = find(c.complexBoundSites - tmpComplexBoundSites);
                msg = [msg sprintf('\tBinding complexes\n')];
                msg = [msg sprintf('\t%8s %6s %7s\n', 'Position', 'Strand', 'Complex')];
                msg = [msg sprintf('\t%8s %6s %7s\n', '========', '======', '=======')];
                for i = 1:size(subs, 1)
                    if vals(i) < 0
                        break;
                    end
                    msg = [msg sprintf('\t%8d %6d %7d\n', subs(i, 1), subs(i, 2), vals(i))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                if any(this.enzymeComplexGlobalIndexs == enzymeIndex)
                    subs2 = posStrnds(c.complexBoundSites(posStrnds) ~= this.enzymeGlobalIndexs(enzymeIndex), :);
                else
                    subs2 = zeros(0, 2);
                end
                msg = [msg sprintf('\tNot binding complexes\n')];
                msg = [msg sprintf('\t%8s %6s\n', 'Position', 'Strand')];
                msg = [msg sprintf('\t%8s %6s\n', '========', '======')];
                for i = 1:size(subs2, 1)
                    msg = [msg sprintf('\t%8d %6d\n', subs2(i, 1), subs2(i, 2))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                msg = [msg sprintf('\tUnbinding complexes\n')];
                msg = [msg sprintf('\t%8s %6s %7s\n', 'Position', 'Strand', 'Complex')];
                msg = [msg sprintf('\t%8s %6s %7s\n', '========', '======', '=======')];
                for i = 1:size(subs, 1)
                    if vals(i) > 0
                        break;
                    end
                    msg = [msg sprintf('\t%8d %6d %7d\n', subs(i, 1), subs(i, 2), -vals(i))]; %#ok<AGROW>
                end
                msg = [msg sprintf('\n')];
                
                throw(MException('ChromosomeProcessAspect:bindProteinToChromosomeStochastically', msg));
            end
            
            function weights = calcRegionWeights(lens)
                weights = max(0, lens - footprint + 1);
            end
            
            function position = calcBindingPosition(randReal, len)
                position = ceil(randReal * (len - footprint + 1));
            end
            
            function [rgnPosStrnds, rgnLens, rgnProbs] = calcNewRegions(rgnPosStrnds, rgnLens, rgnProbs, rgnIdx, offset)
                rgnPosStrnds(end + 1, :) = [
                    rgnPosStrnds(rgnIdx, 1) + offset + footprint, rgnPosStrnds(rgnIdx, 2)];
                rgnLens(end + 1) = rgnLens(rgnIdx) - offset - footprint;
                rgnLens(rgnIdx) = offset;
                rgnProbs([rgnIdx end+1]) = calcRegionWeightsFun(rgnLens([rgnIdx end]));
            end
        end
        
        function [tfs, idxs, positionsStrands, maxBindings, processivityLengths] = bindProteinToChromosome(this, ...
                positionsStrands, proteinIndexs, maxBindings, weights, ...
                isBindingStable, isPositionsStrandFootprintCentroid, ...
                processivityLengths, isBindingProcessive, ...
                ignoreDamageFilter, checkRegionSupercoiled)
            
            %process options
            if nargin < 4 || isempty(maxBindings)
                maxBindings = size(positionsStrands, 1);
            else
                maxBindings = min(maxBindings, size(positionsStrands, 1));
            end
            if maxBindings == 0
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                positionsStrands = zeros(0, 2);
                processivityLengths = zeros(0, 1);
                return;
            end
            
            if islogical(proteinIndexs)
                proteinIndexs = find(proteinIndexs);
            end
            if nargin < 5 || isempty(weights)
                weights = ones(size(positionsStrands, 1), 1);
            end
            if nargin < 6 || isempty(isBindingStable)
                isBindingStable = true;
            end
            if nargin < 7 || isempty(isPositionsStrandFootprintCentroid)
                isPositionsStrandFootprintCentroid = true;
            end
            if nargin < 8 || isempty(processivityLengths)
                processivityLengths = 1;
            end
            if nargin < 9 || isempty(isBindingProcessive)
                isBindingProcessive = false;
            end
            if nargin < 10
                ignoreDamageFilter = [];
            end
            if nargin < 11
                checkRegionSupercoiled = false;
            end
            
            %get protein global indexs
            if isscalar(proteinIndexs)
                monomerGblIndexs = this.enzymeMonomerGlobalIndexs(this.enzymeMonomerLocalIndexs == proteinIndexs);
                complexGblIndexs = this.enzymeComplexGlobalIndexs(this.enzymeComplexLocalIndexs == proteinIndexs);
            else
                monomerGblIndexs = this.enzymeMonomerGlobalIndexs(ismember(this.enzymeMonomerLocalIndexs, proteinIndexs));
                complexGblIndexs = this.enzymeComplexGlobalIndexs(ismember(this.enzymeComplexLocalIndexs, proteinIndexs));
            end
            
            %bind protein to accessible sites
            [releasedMonomers, releasedComplexs, sideEffects, tfs, idxs, positionsStrands, processivityLengths] = ...
                this.chromosome.setSiteProteinBound(...
                positionsStrands, maxBindings, weights, ...
                monomerGblIndexs, complexGblIndexs, ...
                this.enzymeMonomerGlobalIndexs, this.enzymeComplexGlobalIndexs, ...
                isBindingStable, isPositionsStrandFootprintCentroid, processivityLengths, isBindingProcessive, ignoreDamageFilter, checkRegionSupercoiled);
            maxBindings = numel(idxs);
            
            this.boundEnzymes(this.enzymeMonomerLocalIndexs) = this.boundEnzymes(this.enzymeMonomerLocalIndexs) - releasedMonomers;
            this.boundEnzymes(this.enzymeComplexLocalIndexs) = this.boundEnzymes(this.enzymeComplexLocalIndexs) - releasedComplexs;
            this.enzymes(this.enzymeMonomerLocalIndexs)      = this.enzymes(this.enzymeMonomerLocalIndexs)      + releasedMonomers;
            this.enzymes(this.enzymeComplexLocalIndexs)      = this.enzymes(this.enzymeComplexLocalIndexs)      + releasedComplexs;
            
            %side effects
            this.simulationStateSideEffects = [this.simulationStateSideEffects; sideEffects];
        end
        
        function modifyProteinOnChromosome(this, positionsStrands, newProteinIndexs)
            if isempty(positionsStrands)
                return;
            end
            
            if isscalar(newProteinIndexs)
                newMonomerIndexs = this.enzymeMonomerGlobalIndexs(newProteinIndexs == this.enzymeMonomerLocalIndexs);
                if isempty(newMonomerIndexs)
                    newMonomerIndexs = 0;
                    newComplexIndexs = this.enzymeComplexGlobalIndexs(newProteinIndexs == this.enzymeComplexLocalIndexs);
                else
                    newComplexIndexs = 0;
                end
            else
                newMonomerIndexs = zeros(size(positionsStrands, 1), 1);
                [tfs, idxs] = ismember(newProteinIndexs, this.enzymeMonomerLocalIndexs);
                newMonomerIndexs(tfs) = this.enzymeMonomerGlobalIndexs(idxs(tfs));
                
                newComplexIndexs = zeros(size(positionsStrands, 1), 1);
                [tfs, idxs] = ismember(newProteinIndexs, this.enzymeComplexLocalIndexs);
                newComplexIndexs(tfs) = this.enzymeComplexGlobalIndexs(idxs(tfs));
            end
            
            [releasedMonomers, releasedComplexs, sideEffects] = ...
                this.chromosome.modifyBoundProtein(positionsStrands, ...
                newMonomerIndexs, newComplexIndexs, ...
                this.enzymeMonomerGlobalIndexs, this.enzymeComplexGlobalIndexs);
            this.simulationStateSideEffects = [this.simulationStateSideEffects; sideEffects];
            
            this.boundEnzymes(this.enzymeMonomerLocalIndexs) = this.boundEnzymes(this.enzymeMonomerLocalIndexs) - releasedMonomers;
            this.boundEnzymes(this.enzymeComplexLocalIndexs) = this.boundEnzymes(this.enzymeComplexLocalIndexs) - releasedComplexs;
            this.enzymes(this.enzymeMonomerLocalIndexs)      = this.enzymes(this.enzymeMonomerLocalIndexs)      + releasedMonomers;
            this.enzymes(this.enzymeComplexLocalIndexs)      = this.enzymes(this.enzymeComplexLocalIndexs)      + releasedComplexs;
        end
        
        function positionsStrands = releaseProteinFromChromosome(...
                this, proteinLclIdx, rate, protectedPositionsStrands, protectedLengths)
            %process arguments
            if isempty(proteinLclIdx) || rate == 0
                positionsStrands = zeros(0, 2);
                return;
            end
            
            %global index of protein
            if any(this.enzymeMonomerLocalIndexs == proteinLclIdx)
                monomerGblIdx = this.enzymeGlobalIndexs(proteinLclIdx);
                complexGblIdx = [];
            else
                monomerGblIdx = [];
                complexGblIdx = this.enzymeGlobalIndexs(proteinLclIdx);
            end
            
            %release protein from chromosome
            positionsStrands = this.chromosome.stochasticallySetProteinUnbound(monomerGblIdx, complexGblIdx, rate * this.stepSizeSec, ...
                protectedPositionsStrands, protectedLengths, false, false);
            nReleasedProteins = size(positionsStrands, 1);
            if nReleasedProteins == 0
                return;
            end
            
            this.enzymes(proteinLclIdx)      = this.enzymes(proteinLclIdx)      + nReleasedProteins;
            this.boundEnzymes(proteinLclIdx) = this.boundEnzymes(proteinLclIdx) - nReleasedProteins;
        end
        
        function releasedProteins = releaseProteinFromSites(this, positionsStrands, ...
                bothStrands, proteinLclIdx, isPositionsStrandFootprintCentroid, ...
                suspendExternalStateUpdating)
            releasedProteins = zeros(size(this.enzymes));
            if isempty(positionsStrands)
                return;
            end
            
            if nargin < 3 || isempty(bothStrands)
                bothStrands = true;
            end
            
            if nargin < 5
                isPositionsStrandFootprintCentroid = true;
            end
            
            if nargin >= 4 && isPositionsStrandFootprintCentroid
                positionsStrands(isodd(positionsStrands(:, 2)), 1) = ...
                    positionsStrands(isodd(positionsStrands(:, 2)), 1) - ...
                    this.enzymeDNAFootprints5Prime(proteinLclIdx);
                positionsStrands(iseven(positionsStrands(:, 2)), 1) = ...
                    positionsStrands(iseven(positionsStrands(:, 2)), 1) - ...
                    this.enzymeDNAFootprints3Prime(proteinLclIdx);
            end
            
            if nargin < 6
                suspendExternalStateUpdating = false;
            end
            
            [releasedMonomers, releasedComplexs, sideEffects] = ...
                this.chromosome.setRegionProteinUnbound(positionsStrands, 1, ...
                this.enzymeMonomerGlobalIndexs, this.enzymeComplexGlobalIndexs, ...
                bothStrands, bothStrands, suspendExternalStateUpdating, false);
            
            if any(releasedMonomers)
                this.boundEnzymes(this.enzymeMonomerLocalIndexs) = this.boundEnzymes(this.enzymeMonomerLocalIndexs) - releasedMonomers;
                this.enzymes(this.enzymeMonomerLocalIndexs)      = this.enzymes(this.enzymeMonomerLocalIndexs)      + releasedMonomers;
                releasedProteins(this.enzymeMonomerLocalIndexs)  = releasedMonomers;
            end
            
            if any(releasedComplexs)
                this.boundEnzymes(this.enzymeComplexLocalIndexs) = this.boundEnzymes(this.enzymeComplexLocalIndexs) - releasedComplexs;
                this.enzymes(this.enzymeComplexLocalIndexs)      = this.enzymes(this.enzymeComplexLocalIndexs)      + releasedComplexs;
                releasedProteins(this.enzymeComplexLocalIndexs)  = releasedComplexs;
            end
            
            %side effects
            if numel(sideEffects) > 0
                if numel(this.simulationStateSideEffects) == 0
                    this.simulationStateSideEffects = sideEffects;
                else
                    this.simulationStateSideEffects = [this.simulationStateSideEffects; sideEffects];
                end
            end
        end
        
        %Returns true/false whether or not each position and strands is bound by
        %the input enzyme.
        %
        %Note: this.enzymeMonomerLocalIndexs must be sorted
        function result = isDnaBound(this, positionStrands, enzymeIndexs)
            if isempty(positionStrands)
                result = false(0, 1);
                return;
            end
            
            c = this.chromosome;
            
            result = true(size(positionStrands, 1), 1);
            
            tmpTfs = ismembc(enzymeIndexs, this.enzymeMonomerLocalIndexs);
            
            [posStrnds, proteins] = find(c.monomerBoundSites);
            result(tmpTfs) = edu.stanford.covert.util.SparseMat.ismember_subs(...
                [positionStrands(tmpTfs, :) this.enzymeGlobalIndexs(enzymeIndexs(tmpTfs, 1), 1)], [posStrnds proteins], [c.sequenceLen c.nCompartments numel(c.monomerDNAFootprints)]);
            
            [posStrnds, proteins] = find(c.complexBoundSites);
            result(~tmpTfs) = edu.stanford.covert.util.SparseMat.ismember_subs(...
                [positionStrands(~tmpTfs, :) this.enzymeGlobalIndexs(enzymeIndexs(~tmpTfs, 1), 1)], [posStrnds proteins], [c.sequenceLen c.nCompartments numel(c.complexDNAFootprints)]);
        end
        
        function positionsStrands = findProteinInRegion(this, ...
                position, strand, regionLength, enzymeIndex)
            if regionLength == 0
                positionsStrands = zeros(0, 2);
                return;
            end
            
            if any(this.enzymeMonomerLocalIndexs == enzymeIndex)
                enzymeGblIdx = this.enzymeGlobalIndexs(enzymeIndex);
                [positionsStrands, proteins] = find(this.chromosome.monomerBoundSites);
            else
                enzymeGblIdx = this.enzymeGlobalIndexs(enzymeIndex);
                [positionsStrands, proteins] = find(this.chromosome.complexBoundSites);
            end
            
            filter = proteins == enzymeGblIdx;
            filter(filter) = positionsStrands(filter, 2) == strand;
            filter(filter) = positionsStrands(filter, 1) >= position;
            filter(filter) = positionsStrands(filter, 1) < position + regionLength;
            
            positionsStrands = positionsStrands(filter, :);
        end
    end
end
