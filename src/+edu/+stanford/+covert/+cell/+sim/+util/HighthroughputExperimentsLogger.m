%HighthroughputExperimentsLogger.
% Efficiently logs several in silico experiments:
% - Dynamics
%   - Growth (g/s)
%   - Mass (g)
%   - Volume (L)
% - Event times
%   - Replication initiation time (s)
%   - Replication terminatation time (s)
%   - Cell cycle length (s)
% - Averages
%   - Metabolite concentrations (M)
%   - DNA-seq (freq)
%   - RNA-seq (freq)
%   - ChIP-seq (freq)
%   - RNA expression array (M)
%   - Protein expression array (M)
%   - Metabolic reaction fluxes (rxn/s/gDCW)
%
% Input:
% - outPath [.mat file path]: file path where simulated data should be
%   saved.
%
% Output
% - .mat file containing a struct of the simulated data. .mat is saved at
%   the location specified by outPath. Row and column labels are returned
%   in the "labels" field of the returned struct.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/12/2013
classdef HighthroughputExperimentsLogger < edu.stanford.covert.cell.sim.util.Logger
    %config
    properties (Access = protected)
        outPath %output directory
    end
    
    %indices, labels
    properties
        sIdx      %indices of states
        labels    %row, col labels
        
        idxsDnaBoundMon %DNA-bound monomer mapping
        idxsDnaBoundCpx %DNA-bound monomer mapping
    end
    
    %dynamics
    properties (Access = protected)
        time         %s
        growth       %g/s
        mass         %g
        volume       %L
    end
    
    %event timing
    properties (Access = protected)
        repInitTime  %s
        repTermTime  %s
        cellCycleLen %s
    end
    
    %averages
    properties (Access = protected)
        metConcs     %M
        dnaSeq       %freq
        rnaSeq       %freq
        chipSeq      %freq
        rnaArray     %M
        protArray    %M
        rxnFluxes    %rxn/s/gDCW
    end
    
    %temporary values to help calculate averages
    properties
        tmpRnaSeqTranscripts
        tmpRnaSeqNascent
        tmpRnaSeqProcessed
        tmpRnaSeqIntergenic
        
        tmpChipSeqMon
        tmpChipSeqCpx
    end
    
    methods
        function this = HighthroughputExperimentsLogger(outPath)
            this.outPath = outPath;
        end
        
        %Allocates space for small summary stats gathered after each segment.
        function this = initialize(this, sim)
            %indices
            this.sIdx = struct();
            this.sIdx.time       = sim.stateIndex('Time');
            this.sIdx.mass       = sim.stateIndex('Mass');
            this.sIdx.geometry   = sim.stateIndex('Geometry');
            this.sIdx.metabolite = sim.stateIndex('Metabolite');
            this.sIdx.chromosome = sim.stateIndex('Chromosome');
            this.sIdx.rna        = sim.stateIndex('Rna');
            this.sIdx.protmon    = sim.stateIndex('ProteinMonomer');
            this.sIdx.protcpx    = sim.stateIndex('ProteinComplex');
            this.sIdx.metrxn     = sim.stateIndex('MetabolicReaction');
            this.sIdx.transcript = sim.stateIndex('Transcript');
            
            %handles
            g = sim.gene;
            met = sim.states{this.sIdx.metabolite};
            chr = sim.states{this.sIdx.chromosome};
            rna = sim.states{this.sIdx.rna};
            pm  = sim.states{this.sIdx.protmon};
            pc  = sim.states{this.sIdx.protcpx};
            mr  = sim.states{this.sIdx.metrxn};
            
            %indices
            tmp = setdiff(unique([chr.reactionBoundMonomer; find(any(chr.reactionMonomerCatalysisMatrix, 1))']), 0);
            this.idxsDnaBoundMon = zeros(size(pm.matureIndexs));
            this.idxsDnaBoundMon(tmp) = 1:numel(tmp);
            
            tmp = setdiff(unique([chr.reactionBoundComplex; find(any(chr.reactionComplexCatalysisMatrix, 1))']), 0);
            this.idxsDnaBoundCpx = zeros(size(pc.matureIndexs));
            this.idxsDnaBoundCpx(tmp) = 1:numel(tmp);
            
            %dynamics
            this.time   = zeros(1, sim.lengthSec);
            this.growth = zeros(1, sim.lengthSec);
            this.mass   = zeros(1, sim.lengthSec);
            this.volume = zeros(1, sim.lengthSec);
            
            %event timing
            this.repInitTime  = NaN;
            this.repTermTime  = NaN;
            this.cellCycleLen = NaN;
            
            %averages
            this.labels = struct('rows', struct(), 'cols', struct());
            this.labels.rows = struct();
            this.labels.rows.metConcs = met.wholeCellModelIDs;
            this.labels.rows.rnaArray = g.wholeCellModelIDs;
            this.labels.rows.protArray = g.wholeCellModelIDs(g.mRNAIndexs);
            this.labels.rows.rxnFluxes = mr.reactionWholeCellModelIDs;
            this.labels.cols.chipSeq = g.wholeCellModelIDs(g.mRNAIndexs);
            
            this.metConcs  = zeros(numel(this.labels.rows.metConcs), 1);
            this.dnaSeq    = zeros(chr.sequenceLen, 1);
            this.rnaSeq    = zeros(chr.sequenceLen, 1);
            this.chipSeq   = sparse(chr.sequenceLen, numel(this.labels.cols.chipSeq));
            this.rnaArray  = zeros(numel(this.labels.rows.rnaArray), 1);
            this.protArray = zeros(numel(this.labels.rows.protArray), 1);
            this.rxnFluxes = zeros(numel(this.labels.rows.rxnFluxes), 1);
            
            %temporary variables
            this.tmpRnaSeqTranscripts = zeros(chr.sequenceLen, 1);
            this.tmpRnaSeqNascent     = zeros(numel(rna.nascentIndexs), 1);
            this.tmpRnaSeqProcessed   = zeros(numel(rna.processedIndexs), 1);
            this.tmpRnaSeqIntergenic  = zeros(numel(rna.intergenicIndexs), 1);
            
            this.tmpChipSeqMon = sparse(chr.sequenceLen, nnz(this.idxsDnaBoundMon));
            this.tmpChipSeqCpx = sparse(chr.sequenceLen, nnz(this.idxsDnaBoundCpx));
        end
        
        %append
        function this = append(this, sim)
            import edu.stanford.covert.util.ConstantUtil;
            
            %% handles
            g     = sim.gene;
            comp  = sim.compartment;
            sTime = sim.states{this.sIdx.time};
            sMass = sim.states{this.sIdx.mass};
            sGeom = sim.states{this.sIdx.geometry};
            sMet  = sim.states{this.sIdx.metabolite};
            sChr  = sim.states{this.sIdx.chromosome};
            sRna  = sim.states{this.sIdx.rna};
            sPm   = sim.states{this.sIdx.protmon};
            sPc   = sim.states{this.sIdx.protcpx};
            sMr   = sim.states{this.sIdx.metrxn};
            sTrn  = sim.states{this.sIdx.transcript};
            
            %% dynamics
            idx = sTime.values + 1;
            this.time(1, idx)   = sTime.values;
            this.growth(1, idx) = sMr.growth * sMass.cellInitialDryWeight / (1 - sMass.fractionWetWeight);
            this.mass(1, idx)   = sum(sMass.cell);
            this.volume(1, idx) = sGeom.volume;
            
            %% event timing
            if isnan(this.repInitTime)
                if sChr.ploidy > 1
                    this.repInitTime = sTime.values;
                end
            elseif isnan(this.repTermTime)
                if sChr.ploidy >= 2
                    this.repTermTime = sTime.values;
                end
            elseif isnan(this.cellCycleLen)
                if sGeom.pinchedDiameter == 0
                    this.cellCycleLen = sTime.values;
                end
            end
            
            %% sum
            this.metConcs = this.metConcs + sum(sMet.counts(:, [comp.cytosolIndexs comp.membraneIndexs]), 2) / ConstantUtil.nAvogadro / sGeom.volume;
            
            %DNA seq
            [subs, vals] = find(sChr.polymerizedRegions);
            for i = 1:size(subs, 1)
                this.dnaSeq(subs(i, 1) + (0:vals(i) - 1)) = ...
                    + this.dnaSeq(subs(i, 1) + (0:vals(i) - 1)) ...
                    + 1;
            end
            
            %RNA seq
            for i = 1:size(sTrn.boundTranscriptionUnits, 1)
                tuIdx = sTrn.boundTranscriptionUnits(i);
                if tuIdx == 0
                    continue;
                end
                
                start = sTrn.transcriptionUnitFivePrimeCoordinates(tuIdx);
                dir = sTrn.transcriptionUnitDirections(tuIdx);
                len = sTrn.boundTranscriptProgress(i) - 1;
                
                this.tmpRnaSeqTranscripts(start + dir * ((1:len) - 1)) = ...
                    this.tmpRnaSeqTranscripts(start + dir * ((1:len) - 1)) + 1;
            end
            for i = 1:size(sTrn.abortedTranscripts, 1)
                tuIdx = sTrn.abortedTranscripts(i, 1);
                if tuIdx == 0
                    continue;
                end
                
                start = sTrn.transcriptionUnitFivePrimeCoordinates(tuIdx);
                dir = sTrn.transcriptionUnitDirections(tuIdx);
                len = sTrn.abortedTranscripts(i, 2);
                
                this.tmpRnaSeqTranscripts(start + dir * ((1:len) - 1)) = ...
                    this.tmpRnaSeqTranscripts(start + dir * ((1:len) - 1)) + 1;
            end
            rnas = sum(sRna.counts, 2);
            this.tmpRnaSeqNascent = ...
                + this.tmpRnaSeqNascent ...
                + rnas(sRna.nascentIndexs);
            this.tmpRnaSeqProcessed = ...
                + this.tmpRnaSeqProcessed ...
                + rnas(sRna.processedIndexs) ...
                + rnas(sRna.matureIndexs) ...
                + rnas(sRna.boundIndexs) ...
                + rnas(sRna.misfoldedIndexs) ...
                + rnas(sRna.damagedIndexs) ...
                + rnas(sRna.aminoacylatedIndexs);
            this.tmpRnaSeqIntergenic = ...
                + this.tmpRnaSeqIntergenic ...
                + rnas(sRna.intergenicIndexs);
            
            %ChIP-seq
            [monSubs, monVals] = find(sChr.monomerBoundSites);
            [cpxSubs, cpxVals] = find(sChr.complexBoundSites);
            for i = 1:4
                tfs = monSubs(:, 2) == i;
                idxs = sub2ind(size(this.tmpChipSeqMon), monSubs(tfs,1), this.idxsDnaBoundMon(monVals(tfs, 1)));
                this.tmpChipSeqMon(idxs) = this.tmpChipSeqMon(idxs) + 1;
                
                tfs = cpxSubs(:, 2) == i;
                idxs = sub2ind(size(this.tmpChipSeqCpx), cpxSubs(tfs,1), this.idxsDnaBoundCpx(cpxVals(tfs, 1)));
                this.tmpChipSeqCpx(idxs) = this.tmpChipSeqCpx(idxs) + 1;
            end
            
            %RNA, protein array
            rnas = sum(sRna.counts, 2);
            rnas = ...
                + sRna.nascentRNAMatureRNAComposition * rnas(sRna.nascentIndexs) ...
                + rnas(sRna.processedIndexs) ...
                + rnas(sRna.matureIndexs) ...
                + rnas(sRna.boundIndexs) ...
                + rnas(sRna.misfoldedIndexs) ...
                + rnas(sRna.damagedIndexs) ...
                + rnas(sRna.aminoacylatedIndexs);
            mons = sum(reshape(sum(sPm.counts, 2), numel(sPm.matureIndexs), []), 2);
            cpxs = sum(reshape(sum(sPc.counts, 2), numel(sPc.matureIndexs), []), 2);
            
            rnas = rnas / ConstantUtil.nAvogadro / sGeom.volume;
            mons = mons / ConstantUtil.nAvogadro / sGeom.volume;
            cpxs = cpxs / ConstantUtil.nAvogadro / sGeom.volume;
            
            protComp = sum(sPc.proteinComplexComposition, 3);
            
            this.rnaArray = ...
                + this.rnaArray ...
                + sRna.matureRNAGeneComposition * rnas;
            this.rnaArray([g.rRNAIndexs; g.sRNAIndexs; g.tRNAIndexs]) = ...
                + this.rnaArray([g.rRNAIndexs; g.sRNAIndexs; g.tRNAIndexs]) ...
                + protComp([g.rRNAIndexs; g.sRNAIndexs; g.tRNAIndexs], :) * cpxs;
            this.protArray = ...
                + this.protArray ...
                + mons ...
                + protComp(g.mRNAIndexs, :) * cpxs;
            
            %reaction fluxes
            this.rxnFluxes = this.rxnFluxes + sMr.fluxs / sum(sMass.cellDry);
        end
        
        %finalizes log
        function this = finalize(this, sim)
            %% transform
            %handles
            g = sim.gene;
            chr = sim.states{this.sIdx.chromosome};
            rna = sim.states{this.sIdx.rna};
            pc  = sim.states{this.sIdx.protcpx};
            trn = sim.states{this.sIdx.transcript};
            
            %RNA seq
            this.rnaSeq = this.tmpRnaSeqTranscripts;
            for i = size(this.tmpRnaSeqNascent)
                start = trn.transcriptionUnitFivePrimeCoordinates(i);
                dir = trn.transcriptionUnitDirections(i);
                len = trn.transcriptionUnitLengths(i);
                this.rnaSeq(start + dir * ((1:len) - 1)) = ...
                    + this.rnaSeq(start + dir * ((1:len) - 1)) ...
                    + this.tmpRnaSeqNascent(i);
                
                pIdxs = find(rna.nascentRNAMatureRNAComposition(:, i));
                iIdxs = find(rna.intergenicRNAMatrix(:, i));
                if isscalar(pIdxs)
                    this.rnaSeq(start + dir * ((1:len) - 1)) = ...
                        + this.rnaSeq(start + dir * ((1:len) - 1)) ...
                        + this.tmpRnaSeqProcessed(pIdxs);
                else
                    for j = 1:numel(pIdxs)
                        gIdx = g.getIndexs(rna.wholeCellModelIDs{rna.processedIndexs(pIdxs(j))});
                        validateattributes(gIdx, {'numeric'}, {'positive'});
                        
                        if j > 1
                            start = start + len;
                            len = g.startCoordinates(gIdx) - start;
                            validateattributes(len, {'numeric'}, {'nonnegative'});
                            
                            this.rnaSeq(start + ((1:len) - 1)) = ...
                                + this.rnaSeq(start + ((1:len) - 1)) ...
                                + this.tmpRnaSeqIntergenic(iIdxs(j - 1));
                        end
                        
                        start = g.startCoordinates(gIdx);
                        len = g.lengths(gIdx);
                        this.rnaSeq(start + ((1:len) - 1)) = ...
                            + this.rnaSeq(start + ((1:len) - 1)) ...
                            + this.tmpRnaSeqProcessed(pIdxs(j));
                    end
                end
            end
            
            %ChIP seq
            nMon = size(this.tmpChipSeqMon, 2);
            nCpx = size(this.tmpChipSeqCpx, 2);
            
            idxs = find(this.idxsDnaBoundMon);
            for i = 1:nMon
                if any(this.tmpChipSeqMon(:, i))
                    this.chipSeq(:, idxs(i)) = ...
                        + this.chipSeq(:, idxs(i)) ...
                        + sparse(cconv(full(this.tmpChipSeqMon(:, i)), ones(chr.monomerDNAFootprints(idxs(i)), 1), chr.sequenceLen));
                    this.tmpChipSeqMon(:, i) = 0;
                end
            end
            
            idxs = find(this.idxsDnaBoundCpx);
            protComp = sparse(sum(pc.proteinComplexComposition(g.mRNAIndexs, idxs, :), 3));
            for i = 1:nCpx
                if any(this.tmpChipSeqCpx(:, i))
                    tmp = sparse(cconv(full(this.tmpChipSeqCpx(:, i)), ones(chr.complexDNAFootprints(idxs(i)), 1), chr.sequenceLen));
                    this.tmpChipSeqCpx(:, i) = 0;
                    
                    monIdxs = find(protComp(:, i));
                    for j = 1:numel(monIdxs)
                        this.chipSeq(:, monIdxs(j)) = ...
                            + this.chipSeq(:, monIdxs(j)) ...
                            + tmp * protComp(monIdxs(j), i);
                    end
                end
            end
            
            %% average over time
            nTime = numel(this.time);
            this.metConcs  = 1 / nTime * this.metConcs;
            this.dnaSeq    = 1 / nTime * this.dnaSeq;
            this.rnaSeq    = 1 / nTime * this.rnaSeq;
            this.chipSeq   = 1 / nTime * this.chipSeq;
            this.rnaArray  = 1 / nTime * this.rnaArray;
            this.protArray = 1 / nTime * this.protArray;
            this.rxnFluxes = 1 / nTime * this.rxnFluxes;
            
            %% save
            data = struct(...
                'time',         this.time, ...
                'growth',       this.growth, ...
                'mass',         this.mass, ...
                'volume',       this.volume, ...
                'repInitTime',  this.repInitTime, ...
                'repTermTime',  this.repTermTime, ...
                'cellCycleLen', this.cellCycleLen, ...
                'metConcs',     this.metConcs, ...
                'dnaSeq',       this.dnaSeq, ...
                'rnaSeq',       this.rnaSeq, ...
                'chipSeq',      this.chipSeq, ...
                'rnaArray',     this.rnaArray, ...
                'protArray',    this.protArray, ...
                'rxnFluxes',    this.rxnFluxes, ...
                'labels',       this.labels ...
                ); %#ok<NASGU>
            save(this.outPath, '-struct', 'data');
        end
    end
    
    methods (Static = true)
        function [avgVals, labels] = average(inPathPattern)
            %get matching files
            [inPathBase, ~, ~] = fileparts(inPathPattern);
            files = dir(inPathPattern);

            %get longest simulation length
            nTimeMax = 0;
            isFilesSim = false(size(files));
            for i = 1:numel(files)
                tmpPath = fullfile(inPathBase, files(i).name);
                matObj = matfile(tmpPath);
                vars = whos(matObj);
                if ~isequal({vars.name}', {
                        'cellCycleLen'; 'chipSeq'; 'dnaSeq'; 'growth';
                        'labels'; 'mass'; 'metConcs'; 'protArray';
                        'repInitTime'; 'repTermTime'
                        'rnaArray'; 'rnaSeq'; 'rxnFluxes'; 'time'; 'volume'
                        })
                    continue;
                end
                if size(matObj.growth, 1) > 1
                    continue;
                end
                
                isFilesSim(i) = true;
                nTimeMax = max(nTimeMax, numel(matObj.time));
            end
            files = files(isFilesSim);
            fprintf('Averaging %d simulations ...', numel(files));
            
            %prep row labels
            sim = edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil.load();
            g = sim.gene;
            met = sim.state('Metabolite');
            chr = sim.state('Chromosome');
            mr  = sim.state('MetabolicReaction');
            
            labels = struct('rows', struct(), 'cols', struct());
            labels.rows = struct();
            labels.rows.metConcs  = met.wholeCellModelIDs;
            labels.rows.rnaArray  = g.wholeCellModelIDs;
            labels.rows.protArray = g.wholeCellModelIDs(g.mRNAIndexs);
            labels.rows.rxnFluxes = mr.reactionWholeCellModelIDs;
            labels.cols.chipSeq   = g.wholeCellModelIDs(g.mRNAIndexs);
            
            %initialize
            nSim = numel(files);
            
            avgVals = struct();
            
            avgVals.growth = NaN(nSim, nTimeMax);
            avgVals.mass   = NaN(nSim, nTimeMax);
            avgVals.volume = NaN(nSim, nTimeMax);
            
            avgVals.repInitTime  = NaN(nSim, 1);
            avgVals.repTermTime  = NaN(nSim, 1);
            avgVals.cellCycleLen = NaN(nSim, 1);
            
            avgVals.metConcs  = [];
            avgVals.dnaSeq    = [];
            avgVals.rnaSeq    = [];
            avgVals.chipSeq   = [];
            avgVals.rnaArray  = [];
            avgVals.protArray = [];
            avgVals.rxnFluxes = [];
            
            avgVals.metConcs  = zeros(numel(labels.rows.metConcs), 1);
            avgVals.dnaSeq    = zeros(chr.sequenceLen, 1);
            avgVals.rnaSeq    = zeros(chr.sequenceLen, 1);
            avgVals.chipSeq   = sparse(chr.sequenceLen, numel(labels.cols.chipSeq));
            avgVals.rnaArray  = zeros(numel(labels.rows.rnaArray), 1);
            avgVals.protArray = zeros(numel(labels.rows.protArray), 1);
            avgVals.rxnFluxes = zeros(numel(labels.rows.rxnFluxes), 1);
            
            %average
            for i = 1:nSim
                vals = load(fullfile(inPathBase, files(i).name));
                nTime = numel(vals.time);
                
                avgVals.growth(i, 1:nTime) = vals.growth;
                avgVals.mass(  i, 1:nTime) = vals.mass;
                avgVals.volume(i, 1:nTime) = vals.volume;
                
                avgVals.repInitTime( i, 1) = vals.repInitTime;
                avgVals.repTermTime( i, 1) = vals.repTermTime;
                avgVals.cellCycleLen(i, 1) = vals.cellCycleLen;
                
                avgVals.metConcs  = 1 / nSim * vals.metConcs;
                avgVals.dnaSeq    = 1 / nSim * vals.dnaSeq;
                avgVals.rnaSeq    = 1 / nSim * vals.rnaSeq;
                avgVals.chipSeq   = 1 / nSim * vals.chipSeq;
                avgVals.rnaArray  = 1 / nSim * vals.rnaArray;
                avgVals.protArray = 1 / nSim * vals.protArray;
                avgVals.rxnFluxes = 1 / nSim * vals.rxnFluxes;
            end

            %print status
            fprintf('done.\n');
        end
    end
end
