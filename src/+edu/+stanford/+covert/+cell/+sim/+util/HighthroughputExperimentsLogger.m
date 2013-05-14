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
    properties (Constant = true)
        ntResolution = 100; %nt resolution (downsampling) for dnaSeq, rnaSeq, chipSeq
    end
    
    properties (Access = protected)
        outPath      %output directory
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
        dnaSeq       %freq/nt
        rnaSeq       %freq/nt
        chipSeq      %freq/nt
        rnaArray     %M
        protArray    %M
        rxnFluxes    %rxn/s/gDCW
    end
    
    %temporary values to help calculate averages
    properties
        tmpRnaSeqTranscripts
        tmpRnaSeqNascentCnt
        tmpRnaSeqProcessedCnt
        tmpRnaSeqIntergenicCnt
        
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
            
            this.metConcs  = struct('mean', [], 'std', [], 'exp_x2', []);
            this.dnaSeq    = struct('mean', [], 'std', [], 'exp_x2', []);
            this.rnaSeq    = struct('mean', [], 'std', [], 'exp_x2', []);
            this.chipSeq   = struct('mean', [], 'std', [], 'exp_x2', []);
            this.rnaArray  = struct('mean', [], 'std', [], 'exp_x2', []);
            this.protArray = struct('mean', [], 'std', [], 'exp_x2', []);
            this.rxnFluxes = struct('mean', [], 'std', [], 'exp_x2', []);
            
            this.metConcs.mean  = zeros(numel(this.labels.rows.metConcs), 1);
            this.dnaSeq.mean    = zeros(ceil(chr.sequenceLen / this.ntResolution), 1);
            this.rnaSeq.mean    = zeros(ceil(chr.sequenceLen / this.ntResolution), 1);
            this.chipSeq.mean   = sparse(ceil(chr.sequenceLen / this.ntResolution), numel(this.labels.cols.chipSeq));
            this.rnaArray.mean  = zeros(numel(this.labels.rows.rnaArray), 1);
            this.protArray.mean = zeros(numel(this.labels.rows.protArray), 1);
            this.rxnFluxes.mean = zeros(numel(this.labels.rows.rxnFluxes), 1);
            
            this.metConcs.exp_x2  = this.metConcs.mean;
            this.dnaSeq.exp_x2    = this.dnaSeq.mean;
            this.rnaSeq.exp_x2    = this.rnaSeq.mean;
            this.chipSeq.exp_x2   = this.chipSeq.mean;
            this.rnaArray.exp_x2  = this.rnaArray.mean;
            this.protArray.exp_x2 = this.protArray.mean;
            this.rxnFluxes.exp_x2 = this.rxnFluxes.mean;
            
            this.metConcs.std  = this.metConcs.mean;
            this.dnaSeq.std    = this.dnaSeq.mean;
            this.rnaSeq.std    = this.rnaSeq.mean;
            this.chipSeq.std   = this.chipSeq.mean;
            this.rnaArray.std  = this.rnaArray.mean;
            this.protArray.std = this.protArray.mean;
            this.rxnFluxes.std = this.rxnFluxes.mean;
            
            %temporary variables
            this.tmpRnaSeqTranscripts   = struct('mean', [], 'exp_x2', []);
            this.tmpRnaSeqNascentCnt    = struct('mean', [], 'exp_x2', []);
            this.tmpRnaSeqProcessedCnt  = struct('mean', [], 'exp_x2', []);
            this.tmpRnaSeqIntergenicCnt = struct('mean', [], 'exp_x2', []);
            this.tmpChipSeqMon          = struct('mean', [], 'exp_x2', []);
            this.tmpChipSeqCpx          = struct('mean', [], 'exp_x2', []);
            
            this.tmpRnaSeqTranscripts.mean   = zeros(ceil(chr.sequenceLen / this.ntResolution), 1);
            this.tmpRnaSeqNascentCnt.mean    = zeros(numel(rna.nascentIndexs), 1);
            this.tmpRnaSeqProcessedCnt.mean  = zeros(numel(rna.processedIndexs), 1);
            this.tmpRnaSeqIntergenicCnt.mean = zeros(numel(rna.intergenicIndexs), 1);
            this.tmpChipSeqMon.mean          = sparse(chr.sequenceLen, nnz(this.idxsDnaBoundMon));
            this.tmpChipSeqCpx.mean          = sparse(chr.sequenceLen, nnz(this.idxsDnaBoundCpx));
            
            this.tmpRnaSeqTranscripts.exp_x2   = this.tmpRnaSeqTranscripts.mean;
            this.tmpRnaSeqNascentCnt.exp_x2    = this.tmpRnaSeqNascentCnt.mean;
            this.tmpRnaSeqProcessedCnt.exp_x2  = this.tmpRnaSeqProcessedCnt.mean;
            this.tmpRnaSeqIntergenicCnt.exp_x2 = this.tmpRnaSeqIntergenicCnt.mean;
            this.tmpChipSeqMon.exp_x2          = this.tmpChipSeqMon.mean;
            this.tmpChipSeqCpx.exp_x2          = this.tmpChipSeqCpx.mean;
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
            
            %metabolite concentrations
            tmpMet = sum(sMet.counts(:, [comp.cytosolIndexs comp.membraneIndexs]), 2) / ConstantUtil.nAvogadro / sGeom.volume;
            this.metConcs.mean = ...
                + this.metConcs.mean ...
                + tmpMet;
            this.metConcs.exp_x2 = ...
                + this.metConcs.exp_x2 ...
                + tmpMet .^ 2;
            clear tmpMet;
            
            %DNA seq
            [subs, vals] = find(sChr.polymerizedRegions);
            
            tmpDnaSeq = zeros(size(this.dnaSeq.mean));
            for i = 1:size(subs, 1)
                tmpDnaSeq = this.addDownsampled(tmpDnaSeq, subs(i, 1), 1, vals(i), 1);
            end
            
            tmpDnaSeq(1:end - 1) = tmpDnaSeq(1:end - 1) / this.ntResolution;
            tmpDnaSeq(end) = tmpDnaSeq(end) / mod(sChr.sequenceLen, this.ntResolution);
            
            this.dnaSeq.mean   = this.dnaSeq.mean   + tmpDnaSeq;
            this.dnaSeq.exp_x2 = this.dnaSeq.exp_x2 + tmpDnaSeq .^ 2;
            
            clear tmpDnaSeq subs vals;
            
            %RNA seq
            tmpMean = zeros(size(this.tmpRnaSeqTranscripts.mean));
            tus = [sTrn.boundTranscriptionUnits; sTrn.abortedTranscripts(:, 1)];
            lens = [sTrn.boundTranscriptProgress - 1; sTrn.abortedTranscripts(:, 2)];
            for i = 1:numel(tus)
                tuIdx = tus(i);
                if tuIdx == 0
                    continue;
                end
                
                fivePrime = sTrn.transcriptionUnitFivePrimeCoordinates(tuIdx);
                dir = sTrn.transcriptionUnitDirections(tuIdx);
                len = lens(i);
                tmpMean = this.addDownsampled(tmpMean, fivePrime, dir, len, 1);
            end
            tmpMean(1:end-1) = tmpMean(1:end-1) / this.ntResolution;
            tmpMean(1:end)   = tmpMean(1:end)   / mod(sChr.sequenceLen, this.ntResolution);
            this.tmpRnaSeqTranscripts.mean   = this.tmpRnaSeqTranscripts.mean   + tmpMean;
            this.tmpRnaSeqTranscripts.exp_x2 = this.tmpRnaSeqTranscripts.exp_x2 + tmpMean .^ 2;
            
            rnas = sum(sRna.counts, 2);
            
            this.tmpRnaSeqNascentCnt.mean = ...
                + this.tmpRnaSeqNascentCnt.mean ...
                + rnas(sRna.nascentIndexs);
            this.tmpRnaSeqNascentCnt.exp_2 = ...
                + this.tmpRnaSeqNascentCnt.exp_x2 ...
                + rnas(sRna.nascentIndexs) .^ 2;
            
            tmpProcessed = ...
                + rnas(sRna.processedIndexs) ...
                + rnas(sRna.matureIndexs) ...
                + rnas(sRna.boundIndexs) ...
                + rnas(sRna.misfoldedIndexs) ...
                + rnas(sRna.damagedIndexs) ...
                + rnas(sRna.aminoacylatedIndexs);
            this.tmpRnaSeqProcessedCnt.mean = ...
                + this.tmpRnaSeqProcessedCnt.mean ...
                + tmpProcessed;
            this.tmpRnaSeqProcessedCnt.exp_x2 = ...
                + this.tmpRnaSeqProcessedCnt.exp_x2 ...
                + tmpProcessed .^2;
            
            this.tmpRnaSeqIntergenicCnt.mean = ...
                + this.tmpRnaSeqIntergenicCnt.mean ...
                + rnas(sRna.intergenicIndexs);
            this.tmpRnaSeqIntergenicCnt.exp_x2 = ...
                + this.tmpRnaSeqIntergenicCnt.exp_x2 ...
                + rnas(sRna.intergenicIndexs) .^ 2;
            
            clear rnas start dir len;
            
            %ChIP-seq
            [subs, vals] = find(sChr.monomerBoundSites);
            tmp = sparse(size(this.tmpChipSeqMon.mean, 1), size(this.tmpChipSeqMon.mean, 2));
            for i = 1:4
                tfs = subs(:, 2) == i;
                idxs = sub2ind(size(this.tmpChipSeqMon.mean), subs(tfs, 1), this.idxsDnaBoundMon(vals(tfs, 1)));
                tmp(idxs) = tmp(idxs) + 1;
            end
            this.tmpChipSeqMon.mean   = this.tmpChipSeqMon.mean   + tmp;
            this.tmpChipSeqMon.exp_x2 = this.tmpChipSeqMon.exp_x2 + tmp .^ 2;
            
            [subs, vals] = find(sChr.complexBoundSites);
            tmp = sparse(size(this.tmpChipSeqCpx.mean, 1), size(this.tmpChipSeqCpx.mean, 2));
            for i = 1:4
                tfs = subs(:, 2) == i;
                idxs = sub2ind(size(this.tmpChipSeqCpx.mean), subs(tfs, 1), this.idxsDnaBoundCpx(vals(tfs, 1)));
                tmp(idxs) = tmp(idxs) + 1;
            end
            this.tmpChipSeqCpx.mean   = this.tmpChipSeqCpx.mean   + tmp;
            this.tmpChipSeqCpx.exp_x2 = this.tmpChipSeqCpx.exp_x2 + tmp .^ 2;
            
            clear subs vals tfs idxs tmp;
            
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
            
            tmpRna = sRna.matureRNAGeneComposition * rnas;
            tmpRna([g.rRNAIndexs; g.sRNAIndexs; g.tRNAIndexs]) = ...
                + tmpRna([g.rRNAIndexs; g.sRNAIndexs; g.tRNAIndexs]) ...
                + protComp([g.rRNAIndexs; g.sRNAIndexs; g.tRNAIndexs], :) * cpxs;
            this.rnaArray.mean = ...
                + this.rnaArray.mean ...
                + tmpRna;
            this.rnaArray.exp_x2 = ...
                + this.rnaArray.exp_x2 ...
                + tmpRna .^ 2;
            clear tmpRna rnas;
            
            tmpProt = mons + protComp(g.mRNAIndexs, :) * cpxs;
            this.protArray.mean = ...
                + this.protArray.mean ...
                + tmpProt;
            this.protArray.exp_x2 = ...
                + this.protArray.exp_x2 ...
                + tmpProt .^ 2;
            clear tmpProt mons cpxs protComp;
            
            %reaction fluxes
            this.rxnFluxes.mean = ...
                + this.rxnFluxes.mean ...
                + sMr.fluxs / sum(sMass.cellDry);
            this.rxnFluxes.exp_x2 = ...
                + this.rxnFluxes.exp_x2 ...
                + (sMr.fluxs / sum(sMass.cellDry)) .^ 2;
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
            tmpRnaSeqNascent    = struct('mean', [], 'exp_x2', []);
            tmpRnaSeqProcessed  = struct('mean', [], 'exp_x2', []);
            tmpRnaSeqIntergenic = struct('mean', [], 'exp_x2', []);
            
            tmpRnaSeqNascent.mean    = zeros(ceil(chr.sequenceLen / this.ntResolution), 1);
            tmpRnaSeqProcessed.mean  = zeros(ceil(chr.sequenceLen / this.ntResolution), 1);
            tmpRnaSeqIntergenic.mean = zeros(ceil(chr.sequenceLen / this.ntResolution), 1);
            
            tmpRnaSeqNascent.exp_x2    = tmpRnaSeqNascent.mean;
            tmpRnaSeqProcessed.exp_x2  = tmpRnaSeqProcessed.mean;
            tmpRnaSeqIntergenic.exp_x2 = tmpRnaSeqIntergenic.mean;
            
            for i = size(this.tmpRnaSeqNascentCnt)
                fivePrime = trn.transcriptionUnitFivePrimeCoordinates(i);
                dir = trn.transcriptionUnitDirections(i);
                len = trn.transcriptionUnitLengths(i);
                
                tmpRnaSeqNascent.mean = this.addDownsampled(tmpRnaSeqNascent.mean, ...
                    fivePrime, dir, len, this.tmpRnaSeqNascentCnt.mean(i));
                tmpRnaSeqNascent.exp_x2 = this.addDownsampled(tmpRnaSeqNascent.exp_x2, ...
                    fivePrime, dir, len, this.tmpRnaSeqNascentCnt.exp_x2(i));
                
                pIdxs = find(rna.nascentRNAMatureRNAComposition(:, i));
                iIdxs = find(rna.intergenicRNAMatrix(:, i));
                if isscalar(pIdxs)
                    tmpRnaSeqProcessed.mean = this.addDownsampled(tmpRnaSeqProcessed.mean, ...
                        fivePrime, dir, len, this.tmpRnaSeqProcessedCnt.mean(pIdxs));
                    tmpRnaSeqProcessed.exp_x2 = this.addDownsampled(tmpRnaSeqProcessed.exp_x2, ...
                        fivePrime, dir, len, this.tmpRnaSeqProcessedCnt.exp_x2(pIdxs));
                else
                    for j = 1:numel(pIdxs)
                        gIdx = g.getIndexs(rna.wholeCellModelIDs{rna.processedIndexs(pIdxs(j))});
                        validateattributes(gIdx, {'numeric'}, {'positive'});
                        
                        if j > 1
                            start = start + len;
                            len = g.startCoordinates(gIdx) - start;
                            validateattributes(len, {'numeric'}, {'nonnegative'});
                            
                            tmpRnaSeqIntergenic.mean = this.addDownsampled(tmpRnaSeqIntergenic.mean, ...
                                fivePrime, dir, len, this.tmpRnaSeqIntergenicCnt.mean(iIdxs(j - 1)));
                            tmpRnaSeqIntergenic.exp_x2 = this.addDownsampled(tmpRnaSeqIntergenic.exp_x2, ...
                                fivePrime, dir, len, this.tmpRnaSeqIntergenicCnt.exp_x2(iIdxs(j - 1)));
                        end
                        
                        start = g.startCoordinates(gIdx);
                        len = g.lengths(gIdx);
                        tmpRnaSeqProcessed.mean = this.addDownsampled(tmpRnaSeqProcessed.mean, ...
                            fivePrime, dir, len, this.tmpRnaSeqProcessedCnt.mean(pIdxs(j)));
                        tmpRnaSeqProcessed.exp_x2 = this.addDownsampled(tmpRnaSeqProcessed.exp_x2, ...
                            fivePrime, dir, len, this.tmpRnaSeqProcessedCnt.exp_x2(pIdxs(j)));
                    end
                end
            end
            
            tmpRnaSeqNascent.mean(   1:end-1) = tmpRnaSeqNascent.mean(   1:end-1) / this.ntResolution;
            tmpRnaSeqProcessed.mean( 1:end-1) = tmpRnaSeqProcessed.mean( 1:end-1) / this.ntResolution;
            tmpRnaSeqIntergenic.mean(1:end-1) = tmpRnaSeqIntergenic.mean(1:end-1) / this.ntResolution;
            tmpRnaSeqNascent.mean(   1:end)   = tmpRnaSeqNascent.mean(   1:end)   / mod(chr.sequenceLen, this.ntResolution);
            tmpRnaSeqProcessed.mean( 1:end)   = tmpRnaSeqProcessed.mean( 1:end)   / mod(chr.sequenceLen, this.ntResolution);
            tmpRnaSeqIntergenic.mean(1:end)   = tmpRnaSeqIntergenic.mean(1:end)   / mod(chr.sequenceLen, this.ntResolution);
            
            tmpRnaSeqNascent.exp_x2(   1:end-1) = tmpRnaSeqNascent.exp_x2(   1:end-1) / (this.ntResolution ^ 2);
            tmpRnaSeqProcessed.exp_x2( 1:end-1) = tmpRnaSeqProcessed.exp_x2( 1:end-1) / (this.ntResolution ^ 2);
            tmpRnaSeqIntergenic.exp_x2(1:end-1) = tmpRnaSeqIntergenic.exp_x2(1:end-1) / (this.ntResolution ^ 2);
            tmpRnaSeqNascent.exp_x2(   1:end)   = tmpRnaSeqNascent.exp_x2(   1:end)   / (mod(chr.sequenceLen, this.ntResolution) ^ 2);
            tmpRnaSeqProcessed.exp_x2( 1:end)   = tmpRnaSeqProcessed.exp_x2( 1:end)   / (mod(chr.sequenceLen, this.ntResolution) ^ 2);
            tmpRnaSeqIntergenic.exp_x2(1:end)   = tmpRnaSeqIntergenic.exp_x2(1:end)   / (mod(chr.sequenceLen, this.ntResolution) ^ 2);
            
            %ChIP seq
            nMon = size(this.tmpChipSeqMon.mean, 2);
            nCpx = size(this.tmpChipSeqCpx.mean, 2);
            
            idxs = find(this.idxsDnaBoundMon);
            for i = 1:nMon
                if any(this.tmpChipSeqMon.mean(:, i))
                    tmp = sum(reshape([
                        cconv(full(this.tmpChipSeqMon.mean(:, i)), ones(chr.monomerDNAFootprints(idxs(i)), 1), chr.sequenceLen)
                        zeros(-mod(chr.sequenceLen, -this.ntResolution), 1)
                        ], this.ntResolution, []), 1)';
                    tmp(1:end-1) = tmp(1:end-1) / this.ntResolution;
                    tmp(end)     = tmp(end)     / mod(chr.sequenceLen, this.ntResolution);
                    this.chipSeq.mean(:, idxs(i)) = ...
                        + this.chipSeq.mean(:, idxs(i)) ...
                        + tmp;
                    this.tmpChipSeqMon.mean(:, i) = 0;
                    
                    tmp = sum(reshape([
                        cconv(full(this.tmpChipSeqMon.exp_x2(:, i)), ones(chr.monomerDNAFootprints(idxs(i)), 1), chr.sequenceLen)
                        zeros(-mod(chr.sequenceLen, -this.ntResolution), 1)
                        ], this.ntResolution, []), 1)';
                    tmp(1:end-1) = tmp(1:end-1) / this.ntResolution;
                    tmp(end)     = tmp(end)     / mod(chr.sequenceLen, this.ntResolution);
                    this.chipSeq.exp_x2(:, idxs(i)) = ...
                        + this.chipSeq.exp_x2(:, idxs(i)) ...
                        + tmp;
                    this.tmpChipSeqMon.exp_x2(:, i) = 0;
                end
            end
            
            idxs = find(this.idxsDnaBoundCpx);
            protComp = sparse(sum(pc.proteinComplexComposition(g.mRNAIndexs, idxs, :), 3));
            for i = 1:nCpx
                if any(this.tmpChipSeqCpx.mean(:, i))
                    tmp1 = sum(reshape([
                        cconv(full(this.tmpChipSeqCpx.mean(:, i)), ones(chr.complexDNAFootprints(idxs(i)), 1), chr.sequenceLen)
                        zeros(-mod(chr.sequenceLen, -this.ntResolution), 1)
                        ], this.ntResolution, []), 1)';
                    tmp1(1:end-1) = tmp1(1:end-1) / this.ntResolution;
                    tmp1(end)     = tmp1(end)     / mod(chr.sequenceLen, this.ntResolution);
                    
                    tmp2 = sum(reshape([
                        cconv(full(this.tmpChipSeqCpx.exp_x2(:, i)), ones(chr.complexDNAFootprints(idxs(i)), 1), chr.sequenceLen)
                        zeros(-mod(chr.sequenceLen, -this.ntResolution), 1)
                        ], this.ntResolution, []), 1)';
                    tmp2(1:end-1) = tmp2(1:end-1) / this.ntResolution;
                    tmp2(end)     = tmp2(end)     / mod(chr.sequenceLen, this.ntResolution);
                    
                    this.tmpChipSeqCpx.mean(:, i) = 0;
                    this.tmpChipSeqCpx.exp_x2(:, i) = 0;
                    
                    monIdxs = find(protComp(:, i));
                    for j = 1:numel(monIdxs)
                        this.chipSeq.mean(:, monIdxs(j)) = ...
                            + this.chipSeq.mean(:, monIdxs(j)) ...
                            + tmp1 * protComp(monIdxs(j), i);
                        this.chipSeq.exp_x2(:, monIdxs(j)) = ...
                            + this.chipSeq.exp_x2(:, monIdxs(j)) ...
                            + tmp2 * protComp(monIdxs(j), i);
                    end
                end
            end
            
            %% average over time
            nTime = numel(this.time);
            
            this.metConcs.mean  = 1 / nTime * this.metConcs.mean;
            this.dnaSeq.mean    = 1 / nTime * this.dnaSeq.mean;
            this.chipSeq.mean   = 1 / nTime * this.chipSeq.mean;
            this.rnaArray.mean  = 1 / nTime * this.rnaArray.mean;
            this.protArray.mean = 1 / nTime * this.protArray.mean;
            this.rxnFluxes.mean = 1 / nTime * this.rxnFluxes.mean;
            
            this.metConcs.exp_x2  = 1 / nTime * this.metConcs.exp_x2;
            this.dnaSeq.exp_x2    = 1 / nTime * this.dnaSeq.exp_x2;
            this.chipSeq.exp_x2   = 1 / nTime * this.chipSeq.exp_x2;
            this.rnaArray.exp_x2  = 1 / nTime * this.rnaArray.exp_x2;
            this.protArray.exp_x2 = 1 / nTime * this.protArray.exp_x2;
            this.rxnFluxes.exp_x2 = 1 / nTime * this.rxnFluxes.exp_x2;
            
            this.metConcs.std  = sqrt(this.metConcs.exp_x2  - this.metConcs.mean  .^ 2);
            this.dnaSeq.std    = sqrt(this.dnaSeq.exp_x2    - this.dnaSeq.mean    .^ 2);
            this.chipSeq.std   = sqrt(this.chipSeq.exp_x2   - this.chipSeq.mean   .^ 2);
            this.rnaArray.std  = sqrt(this.rnaArray.exp_x2  - this.rnaArray.mean  .^ 2);
            this.protArray.std = sqrt(this.protArray.exp_x2 - this.protArray.mean .^ 2);
            this.rxnFluxes.std = sqrt(this.rxnFluxes.exp_x2 - this.rxnFluxes.mean .^ 2);
            
            this.metConcs  = rmfield(this.metConcs,  'exp_x2');
            this.dnaSeq    = rmfield(this.dnaSeq,    'exp_x2');
            this.chipSeq   = rmfield(this.chipSeq,   'exp_x2');
            this.rnaArray  = rmfield(this.rnaArray,  'exp_x2');
            this.protArray = rmfield(this.protArray, 'exp_x2');
            this.rxnFluxes = rmfield(this.rxnFluxes, 'exp_x2');
            
            %RNA seq
            this.tmpRnaSeqTranscripts.mean = 1 / nTime * this.tmpRnaSeqTranscripts.mean;
            tmpRnaSeqNascent.mean          = 1 / nTime * tmpRnaSeqNascent.mean;
            tmpRnaSeqProcessed.mean        = 1 / nTime * tmpRnaSeqProcessed.mean;
            tmpRnaSeqIntergenic.mean       = 1 / nTime * tmpRnaSeqIntergenic.mean;
            
            this.tmpRnaSeqTranscripts.exp_x2 = 1 / nTime * this.tmpRnaSeqTranscripts.exp_x2;
            tmpRnaSeqNascent.exp_x2          = 1 / nTime * tmpRnaSeqNascent.exp_x2;
            tmpRnaSeqProcessed.exp_x2        = 1 / nTime * tmpRnaSeqProcessed.exp_x2;
            tmpRnaSeqIntergenic.exp_x2       = 1 / nTime * tmpRnaSeqIntergenic.exp_x2;
            
            this.tmpRnaSeqTranscripts.var    = this.tmpRnaSeqTranscripts.exp_x2 - this.tmpRnaSeqTranscripts.mean .^ 2;
            tmpRnaSeqNascent.var             = tmpRnaSeqNascent.exp_x2          - tmpRnaSeqNascent.mean          .^ 2;
            tmpRnaSeqProcessed.var           = tmpRnaSeqProcessed.exp_x2        - tmpRnaSeqProcessed.mean        .^ 2;
            tmpRnaSeqIntergenic.var          = tmpRnaSeqIntergenic.exp_x2       - tmpRnaSeqIntergenic.mean       .^ 2;
            
            this.rnaSeq.mean = ...
                + this.tmpRnaSeqTranscripts.mean ...
                + tmpRnaSeqNascent.mean ...
                + tmpRnaSeqProcessed.mean ...
                + tmpRnaSeqIntergenic.mean;
            this.rnaSeq.std = sqrt(...
                + this.tmpRnaSeqTranscripts.var ...
                + tmpRnaSeqNascent.var ...
                + tmpRnaSeqProcessed.var ...
                + tmpRnaSeqIntergenic.var ...
                );
            
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
            import edu.stanford.covert.cell.sim.util.HighthroughputExperimentsLogger;
            
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
            
            avgVals.metConcs  = struct('mean', [], 'std', [], 'var', []);
            avgVals.dnaSeq    = struct('mean', [], 'std', [], 'var', []);
            avgVals.rnaSeq    = struct('mean', [], 'std', [], 'var', []);
            avgVals.chipSeq   = struct('mean', [], 'std', [], 'var', []);
            avgVals.rnaArray  = struct('mean', [], 'std', [], 'var', []);
            avgVals.protArray = struct('mean', [], 'std', [], 'var', []);
            avgVals.rxnFluxes = struct('mean', [], 'std', [], 'var', []);
            
            avgVals.metConcs.mean  = zeros(numel(labels.rows.metConcs), 1);
            avgVals.dnaSeq.mean    = zeros(ceil(chr.sequenceLen  / HighthroughputExperimentsLogger.ntResolution), 1);
            avgVals.rnaSeq.mean    = zeros(ceil(chr.sequenceLen  / HighthroughputExperimentsLogger.ntResolution), 1);
            avgVals.chipSeq.mean   = sparse(ceil(chr.sequenceLen / HighthroughputExperimentsLogger.ntResolution), numel(labels.cols.chipSeq));
            avgVals.rnaArray.mean  = zeros(numel(labels.rows.rnaArray), 1);
            avgVals.protArray.mean = zeros(numel(labels.rows.protArray), 1);
            avgVals.rxnFluxes.mean = zeros(numel(labels.rows.rxnFluxes), 1);
            
            avgVals.metConcs.var  = avgVals.metConcs.mean;
            avgVals.dnaSeq.var    = avgVals.dnaSeq.mean;
            avgVals.rnaSeq.var    = avgVals.rnaSeq.mean;
            avgVals.chipSeq.var   = avgVals.chipSeq.mean;
            avgVals.rnaArray.var  = avgVals.rnaArray.mean;
            avgVals.protArray.var = avgVals.protArray.mean;
            avgVals.rxnFluxes.var = avgVals.rxnFluxes.mean;
            
            avgVals.metConcs.std  = avgVals.metConcs.mean;
            avgVals.dnaSeq.std    = avgVals.dnaSeq.mean;
            avgVals.rnaSeq.std    = avgVals.rnaSeq.mean;
            avgVals.chipSeq.std   = avgVals.chipSeq.mean;
            avgVals.rnaArray.std  = avgVals.rnaArray.mean;
            avgVals.protArray.std = avgVals.protArray.mean;
            avgVals.rxnFluxes.std = avgVals.rxnFluxes.mean;
            
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
                
                avgVals.metConcs.mean  = avgVals.metConcs.mean  + 1 / nSim * vals.metConcs.mean;
                avgVals.dnaSeq.mean    = avgVals.dnaSeq.mean    + 1 / nSim * vals.dnaSeq.mean;
                avgVals.rnaSeq.mean    = avgVals.rnaSeq.mean    + 1 / nSim * vals.rnaSeq.mean;
                avgVals.chipSeq.mean   = avgVals.chipSeq.mean   + 1 / nSim * vals.chipSeq.mean;
                avgVals.rnaArray.mean  = avgVals.rnaArray.mean  + 1 / nSim * vals.rnaArray.mean;
                avgVals.protArray.mean = avgVals.protArray.mean + 1 / nSim * vals.protArray.mean;
                avgVals.rxnFluxes.mean = avgVals.rxnFluxes.mean + 1 / nSim * vals.rxnFluxes.mean;
                
                avgVals.metConcs.var  = avgVals.metConcs.var  + (1 / nSim * vals.metConcs.std ) .^ 2;
                avgVals.dnaSeq.var    = avgVals.dnaSeq.var    + (1 / nSim * vals.dnaSeq.std   ) .^ 2;
                avgVals.rnaSeq.var    = avgVals.rnaSeq.var    + (1 / nSim * vals.rnaSeq.std   ) .^ 2;
                avgVals.chipSeq.var   = avgVals.chipSeq.var   + (1 / nSim * vals.chipSeq.std  ) .^ 2;
                avgVals.rnaArray.var  = avgVals.rnaArray.var  + (1 / nSim * vals.rnaArray.std ) .^ 2;
                avgVals.protArray.var = avgVals.protArray.var + (1 / nSim * vals.protArray.std) .^ 2;
                avgVals.rxnFluxes.var = avgVals.rxnFluxes.var + (1 / nSim * vals.rxnFluxes.std) .^ 2;
            end
            
            avgVals.metConcs.std  = sqrt(avgVals.metConcs.var );
            avgVals.dnaSeq.std    = sqrt(avgVals.dnaSeq.var   );
            avgVals.rnaSeq.std    = sqrt(avgVals.rnaSeq.var   );
            avgVals.chipSeq.std   = sqrt(avgVals.chipSeq.var  );
            avgVals.rnaArray.std  = sqrt(avgVals.rnaArray.var );
            avgVals.protArray.std = sqrt(avgVals.protArray.var);
            avgVals.rxnFluxes.std = sqrt(avgVals.rxnFluxes.var);
            
            avgVals.metConcs  = rmfield(avgVals.metConcs,  'var');
            avgVals.dnaSeq    = rmfield(avgVals.dnaSeq,    'var');
            avgVals.rnaSeq    = rmfield(avgVals.rnaSeq,    'var');
            avgVals.chipSeq   = rmfield(avgVals.chipSeq,   'var');
            avgVals.rnaArray  = rmfield(avgVals.rnaArray,  'var');
            avgVals.protArray = rmfield(avgVals.protArray, 'var');
            avgVals.rxnFluxes = rmfield(avgVals.rxnFluxes, 'var');
            
            %print status
            fprintf('done.\n');
        end
    end
    
    %helpers
    methods (Static = true)
        function arr = addDownsampled(arr, fivePrime, dir, len, val)
            import edu.stanford.covert.cell.sim.util.HighthroughputExperimentsLogger;
            
            if len < 1
                return;
            end
            
            if dir
                startCoor = fivePrime;
                endCoor   = fivePrime + (len - 1);
            else
                startCoor = fivePrime - (len - 1);
                endCoor   = fivePrime;
            end
            
            startSeg = ceil(startCoor / HighthroughputExperimentsLogger.ntResolution);
            endSeg = ceil(endCoor / HighthroughputExperimentsLogger.ntResolution);
            
            arr(startSeg) = arr(startSeg) - mod(startCoor, -HighthroughputExperimentsLogger.ntResolution) * val;
            arr(endSeg)   = arr(endSeg)   + mod(endCoor,    HighthroughputExperimentsLogger.ntResolution) * val;
            arr(startSeg+1:endSeg-1) = ...
                + arr(startSeg+1:endSeg-1) ...
                + HighthroughputExperimentsLogger.ntResolution * val;
        end
    end
end
