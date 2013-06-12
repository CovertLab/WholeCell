%HighthroughputExperimentsLogger.
% Efficiently logs several in silico experiments:
% - Single-cell data
%   - Dynamics: rows correspond to individual cells, columns correspond to timepoints 0..N (s)
%     - Growth (g/s)
%     - Mass (g)
%     - Volume (L)
%   - Event times: rows correspond to individual cells
%     - Replication initiation time (s)
%     - Replication termination time (s)
%     - Cell cycle length (s)
% - Population and time averages: properties are each structs with two
%   fields containing the mean and standard deviation
%   - Metabolite concentrations (M)
%     Rows correspond to metabolite species. Rows are labeled by
%     Metabolite.wholeCellModelIDs
%   - DNA-seq (DNA molecules/nt)
%     Average DNA copy number of each 100 nt region of the chromosome. Row
%     1 corresponds to bases 1..100, Rows 2 corresponds to bases 101..200,
%     etc.
%   - RNA-seq (transcripts/nt)
%     Average number of mapped RNA transcripts of each 100 nt region of the
%     chromosome. Row 1 corresponds to bases 1..100, row 2 corresponds to
%     bases 101..200, etc. 
%   - ChIP-seq (protein molecules/nt)
%     Average DNA-bound protein density of each 100 nt region of the
%     chromosome. Row 1 corresponds to bases 1..100, row 2 corresponds 
%     to bases 101..200, etc. Columns corresponds to mRNA-coding genes and
%     are labeled by gene.wholeCellModelIDs(gene.mRNAIndexs)
%   - RNA expression array (M)
%     Concentrations of RNA by gene. Rows correspond to genes. Rows are
%     labeled by gene.wholeCellModelIDs.
%   - Protein expression array (M)
%     Concentrations of protein. Rows correspond to protein-coding genes.
%     Rows are labeled by gene.wholeCellModelIDs(gene.mRNAIndexs).
%   - Metabolic reaction fluxes (rxn/s/gDCW)
%     Rows are labeled by MetabolicReaction.reactionWholeCellModelIDs
%
% Input:
% - simPath [.mat file path]: file path where simulated data should be
%   saved.
%
% Output
% - .mat file containing a struct of the simulated data. .mat is saved at
%   the location specified by simPath. Row and column labels are returned
%   in the "labels" field of the returned struct.
%
% Notation
% - N = no. data points (time) to average over
% - mean = E[x] = sum / N
% - exp_x2 = E[x^2] = sum_x2 / N
% - var = E[x^2] - E[x]^2
% - std = sqrt(var)
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
        simPath      %output directory
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
        function this = HighthroughputExperimentsLogger(simPath)
            this.simPath = simPath;
        end
        
        %Allocates space for small summary stats gathered after each segment.
        function this = initialize(this, sim)
            %% initialize
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
            this.time   = zeros(1, sim.lengthSec + 1);
            this.growth = zeros(1, sim.lengthSec + 1);
            this.mass   = zeros(1, sim.lengthSec + 1);
            this.volume = zeros(1, sim.lengthSec + 1);
            
            %event timing
            this.repInitTime  = NaN;
            this.repTermTime  = NaN;
            this.cellCycleLen = NaN;
            
            %averages: sum, sum_x2 are temporary fields used to calculate
            %the variance online using little memory
            this.labels = struct('rows', struct(), 'cols', struct());
            this.labels.rows = struct();
            this.labels.rows.metConcs = met.wholeCellModelIDs;
            this.labels.rows.rnaArray = g.wholeCellModelIDs;
            this.labels.rows.protArray = g.wholeCellModelIDs(g.mRNAIndexs);
            this.labels.rows.rxnFluxes = mr.reactionWholeCellModelIDs;
            this.labels.cols.chipSeq = g.wholeCellModelIDs(g.mRNAIndexs);
            
            this.metConcs  = struct('mean', [], 'std', [], 'sum', [], 'sum_x2', []);
            this.dnaSeq    = struct('mean', [], 'std', [], 'sum', [], 'sum_x2', []);
            this.rnaSeq    = struct('mean', [], 'std', [], 'sum', [], 'sum_x2', []);
            this.chipSeq   = struct('mean', [], 'std', [], 'sum', [], 'sum_x2', []);
            this.rnaArray  = struct('mean', [], 'std', [], 'sum', [], 'sum_x2', []);
            this.protArray = struct('mean', [], 'std', [], 'sum', [], 'sum_x2', []);
            this.rxnFluxes = struct('mean', [], 'std', [], 'sum', [], 'sum_x2', []);
            
            this.metConcs.sum  = zeros(numel(this.labels.rows.metConcs), 1);
            this.dnaSeq.sum    = zeros(ceil(chr.sequenceLen / this.ntResolution), 1);
            this.rnaSeq.sum    = zeros(ceil(chr.sequenceLen / this.ntResolution), 1);
            this.chipSeq.sum   = zeros(ceil(chr.sequenceLen / this.ntResolution), numel(this.labels.cols.chipSeq));
            this.rnaArray.sum  = zeros(numel(this.labels.rows.rnaArray), 1);
            this.protArray.sum = zeros(numel(this.labels.rows.protArray), 1);
            this.rxnFluxes.sum = zeros(numel(this.labels.rows.rxnFluxes), 1);
            
            this.metConcs.sum_x2  = this.metConcs.sum;
            this.dnaSeq.sum_x2    = this.dnaSeq.sum;
            this.rnaSeq.sum_x2    = this.rnaSeq.sum;
            this.chipSeq.sum_x2   = this.chipSeq.sum;
            this.rnaArray.sum_x2  = this.rnaArray.sum;
            this.protArray.sum_x2 = this.protArray.sum;
            this.rxnFluxes.sum_x2 = this.rxnFluxes.sum;
            
            %temporary variables
            this.tmpRnaSeqTranscripts   = struct('sum', [], 'sum_x2', []);
            this.tmpRnaSeqNascentCnt    = struct('sum', [], 'sum_x2', []);
            this.tmpRnaSeqProcessedCnt  = struct('sum', [], 'sum_x2', []);
            this.tmpRnaSeqIntergenicCnt = struct('sum', [], 'sum_x2', []);
            this.tmpChipSeqMon          = struct('sum', [], 'sum_x2', []);
            this.tmpChipSeqCpx          = struct('sum', [], 'sum_x2', []);
            
            this.tmpRnaSeqTranscripts.sum   = zeros(ceil(chr.sequenceLen / this.ntResolution), 1);
            this.tmpRnaSeqNascentCnt.sum    = zeros(numel(rna.nascentIndexs), 1);
            this.tmpRnaSeqProcessedCnt.sum  = zeros(numel(rna.processedIndexs), 1);
            this.tmpRnaSeqIntergenicCnt.sum = zeros(numel(rna.intergenicIndexs), 1);
            this.tmpChipSeqMon.sum          = sparse(chr.sequenceLen, nnz(this.idxsDnaBoundMon));
            this.tmpChipSeqCpx.sum          = sparse(chr.sequenceLen, nnz(this.idxsDnaBoundCpx));
            
            this.tmpRnaSeqTranscripts.sum_x2   = this.tmpRnaSeqTranscripts.sum;
            this.tmpRnaSeqNascentCnt.sum_x2    = this.tmpRnaSeqNascentCnt.sum;
            this.tmpRnaSeqProcessedCnt.sum_x2  = this.tmpRnaSeqProcessedCnt.sum;
            this.tmpRnaSeqIntergenicCnt.sum_x2 = this.tmpRnaSeqIntergenicCnt.sum;
            this.tmpChipSeqMon.sum_x2          = this.tmpChipSeqMon.sum;
            this.tmpChipSeqCpx.sum_x2          = this.tmpChipSeqCpx.sum;
            
            %% log
            this.append(sim);
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
            this.metConcs.sum = ...
                + this.metConcs.sum ...
                + tmpMet;
            this.metConcs.sum_x2 = ...
                + this.metConcs.sum_x2 ...
                + tmpMet .^ 2;
            clear tmpMet;
            
            %DNA seq
            [subs, vals] = find(sChr.polymerizedRegions);
            
            tmpDnaSeq = zeros(size(this.dnaSeq.sum));
            for i = 1:size(subs, 1)
                tmpDnaSeq = this.addDownsampled(tmpDnaSeq, subs(i, 1), 1, vals(i), 1);
            end
            
            tmpDnaSeq(1:end - 1) = tmpDnaSeq(1:end - 1) / this.ntResolution;
            tmpDnaSeq(end)       = tmpDnaSeq(end)       / mod(sChr.sequenceLen, this.ntResolution);
            
            this.dnaSeq.sum    = this.dnaSeq.sum    + tmpDnaSeq;
            this.dnaSeq.sum_x2 = this.dnaSeq.sum_x2 + tmpDnaSeq .^ 2;
            
            clear tmpDnaSeq subs vals;
            
            %RNA seq
            tmpSum = zeros(size(this.tmpRnaSeqTranscripts.sum));
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
                tmpSum = this.addDownsampled(tmpSum, fivePrime, dir, len, 1);
            end
            tmpSum(1:end-1) = tmpSum(1:end-1) / this.ntResolution;
            tmpSum(1:end)   = tmpSum(1:end)   / mod(sChr.sequenceLen, this.ntResolution);
            this.tmpRnaSeqTranscripts.sum    = this.tmpRnaSeqTranscripts.sum    + tmpSum;
            this.tmpRnaSeqTranscripts.sum_x2 = this.tmpRnaSeqTranscripts.sum_x2 + tmpSum .^ 2;
            
            rnas = sum(sRna.counts, 2);
            
            this.tmpRnaSeqNascentCnt.sum = ...
                + this.tmpRnaSeqNascentCnt.sum ...
                + rnas(sRna.nascentIndexs);
            this.tmpRnaSeqNascentCnt.exp_2 = ...
                + this.tmpRnaSeqNascentCnt.sum_x2 ...
                + rnas(sRna.nascentIndexs) .^ 2;
            
            tmpProcessed = ...
                + rnas(sRna.processedIndexs) ...
                + rnas(sRna.matureIndexs) ...
                + rnas(sRna.boundIndexs) ...
                + rnas(sRna.misfoldedIndexs) ...
                + rnas(sRna.damagedIndexs) ...
                + rnas(sRna.aminoacylatedIndexs);
            this.tmpRnaSeqProcessedCnt.sum = ...
                + this.tmpRnaSeqProcessedCnt.sum ...
                + tmpProcessed;
            this.tmpRnaSeqProcessedCnt.sum_x2 = ...
                + this.tmpRnaSeqProcessedCnt.sum_x2 ...
                + tmpProcessed .^2;
            
            this.tmpRnaSeqIntergenicCnt.sum = ...
                + this.tmpRnaSeqIntergenicCnt.sum ...
                + rnas(sRna.intergenicIndexs);
            this.tmpRnaSeqIntergenicCnt.sum_x2 = ...
                + this.tmpRnaSeqIntergenicCnt.sum_x2 ...
                + rnas(sRna.intergenicIndexs) .^ 2;
            
            clear rnas start dir len;
            
            %ChIP-seq
            [subs, vals] = find(sChr.monomerBoundSites);
            tmp = sparse(size(this.tmpChipSeqMon.sum, 1), size(this.tmpChipSeqMon.sum, 2));
            for i = 1:4
                tfs = subs(:, 2) == i;
                idxs = sub2ind(size(this.tmpChipSeqMon.sum), subs(tfs, 1), this.idxsDnaBoundMon(vals(tfs, 1)));
                tmp(idxs) = tmp(idxs) + 1;
            end
            this.tmpChipSeqMon.sum    = this.tmpChipSeqMon.sum    + tmp;
            this.tmpChipSeqMon.sum_x2 = this.tmpChipSeqMon.sum_x2 + tmp .^ 2;
            
            [subs, vals] = find(sChr.complexBoundSites);
            tmp = sparse(size(this.tmpChipSeqCpx.sum, 1), size(this.tmpChipSeqCpx.sum, 2));
            for i = 1:4
                tfs = subs(:, 2) == i;
                idxs = sub2ind(size(this.tmpChipSeqCpx.sum), subs(tfs, 1), this.idxsDnaBoundCpx(vals(tfs, 1)));
                tmp(idxs) = tmp(idxs) + 1;
            end
            this.tmpChipSeqCpx.sum    = this.tmpChipSeqCpx.sum    + tmp;
            this.tmpChipSeqCpx.sum_x2 = this.tmpChipSeqCpx.sum_x2 + tmp .^ 2;
            
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
            this.rnaArray.sum = ...
                + this.rnaArray.sum ...
                + tmpRna;
            this.rnaArray.sum_x2 = ...
                + this.rnaArray.sum_x2 ...
                + tmpRna .^ 2;
            clear tmpRna rnas;
            
            tmpProt = mons + protComp(g.mRNAIndexs, :) * cpxs;
            this.protArray.sum = ...
                + this.protArray.sum ...
                + tmpProt;
            this.protArray.sum_x2 = ...
                + this.protArray.sum_x2 ...
                + tmpProt .^ 2;
            clear tmpProt mons cpxs protComp;
            
            %reaction fluxes
            tmpFluxes = sMr.fluxs / sum(sMass.cellDry);
            this.rxnFluxes.sum = ...
                + this.rxnFluxes.sum ...
                + tmpFluxes;
            this.rxnFluxes.sum_x2 = ...
                + this.rxnFluxes.sum_x2 ...
                + tmpFluxes .^ 2;
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
            
            %% average over time
            %calculate means and variances
            nTime = numel(this.time);
            
            this.metConcs.mean  = 1 / nTime * this.metConcs.sum;
            this.dnaSeq.mean    = 1 / nTime * this.dnaSeq.sum;
            this.rnaArray.mean  = 1 / nTime * this.rnaArray.sum;
            this.protArray.mean = 1 / nTime * this.protArray.sum;
            this.rxnFluxes.mean = 1 / nTime * this.rxnFluxes.sum;
            
            this.metConcs.exp_x2  = 1 / nTime * this.metConcs.sum_x2;
            this.dnaSeq.exp_x2    = 1 / nTime * this.dnaSeq.sum_x2;
            this.rnaArray.exp_x2  = 1 / nTime * this.rnaArray.sum_x2;
            this.protArray.exp_x2 = 1 / nTime * this.protArray.sum_x2;
            this.rxnFluxes.exp_x2 = 1 / nTime * this.rxnFluxes.sum_x2;
            
            this.metConcs.std  = sqrt(max(0, this.metConcs.exp_x2  - this.metConcs.mean  .^ 2));
            this.dnaSeq.std    = sqrt(max(0, this.dnaSeq.exp_x2    - this.dnaSeq.mean    .^ 2));
            this.rnaArray.std  = sqrt(max(0, this.rnaArray.exp_x2  - this.rnaArray.mean  .^ 2));
            this.protArray.std = sqrt(max(0, this.protArray.exp_x2 - this.protArray.mean .^ 2));
            this.rxnFluxes.std = sqrt(max(0, this.rxnFluxes.exp_x2 - this.rxnFluxes.mean .^ 2));
            
            this.metConcs  = rmfield(this.metConcs,  {'sum', 'sum_x2', 'exp_x2'});
            this.dnaSeq    = rmfield(this.dnaSeq,    {'sum', 'sum_x2', 'exp_x2'});
            this.rnaArray  = rmfield(this.rnaArray,  {'sum', 'sum_x2', 'exp_x2'});
            this.protArray = rmfield(this.protArray, {'sum', 'sum_x2', 'exp_x2'});
            this.rxnFluxes = rmfield(this.rxnFluxes, {'sum', 'sum_x2', 'exp_x2'});
            
            %RNA seq
            this.tmpRnaSeqTranscripts.mean   = 1 / nTime * this.tmpRnaSeqTranscripts.sum;
            this.tmpRnaSeqNascentCnt.mean    = 1 / nTime * this.tmpRnaSeqNascentCnt.sum;
            this.tmpRnaSeqProcessedCnt.mean  = 1 / nTime * this.tmpRnaSeqProcessedCnt.sum; 
            this.tmpRnaSeqIntergenicCnt.mean = 1 / nTime * this.tmpRnaSeqIntergenicCnt.sum; 
            
            this.tmpRnaSeqTranscripts.exp_x2   = 1 / nTime * this.tmpRnaSeqTranscripts.sum_x2; 
            this.tmpRnaSeqNascentCnt.exp_x2    = 1 / nTime * this.tmpRnaSeqNascentCnt.sum_x2; 
            this.tmpRnaSeqProcessedCnt.exp_x2  = 1 / nTime * this.tmpRnaSeqProcessedCnt.sum_x2; 
            this.tmpRnaSeqIntergenicCnt.exp_x2 = 1 / nTime * this.tmpRnaSeqIntergenicCnt.sum_x2; 
            
            this.tmpRnaSeqTranscripts.var   = max(0, this.tmpRnaSeqTranscripts.exp_x2 - this.tmpRnaSeqTranscripts.mean .^ 2);
            this.tmpRnaSeqNascentCnt.var    = max(0, this.tmpRnaSeqNascentCnt.exp_x2 - this.tmpRnaSeqNascentCnt.mean .^ 2);
            this.tmpRnaSeqProcessedCnt.var  = max(0, this.tmpRnaSeqProcessedCnt.exp_x2 - this.tmpRnaSeqProcessedCnt.mean .^ 2); 
            this.tmpRnaSeqIntergenicCnt.var = max(0, this.tmpRnaSeqIntergenicCnt.exp_x2 - this.tmpRnaSeqIntergenicCnt.mean .^ 2); 
                        
            tmpRnaSeqNascent    = struct('mean', [], 'var', []);
            tmpRnaSeqProcessed  = struct('mean', [], 'var', []);
            tmpRnaSeqIntergenic = struct('mean', [], 'var', []);
            
            tmpRnaSeqNascent.mean    = zeros(ceil(chr.sequenceLen / this.ntResolution), 1);
            tmpRnaSeqProcessed.mean  = zeros(ceil(chr.sequenceLen / this.ntResolution), 1);
            tmpRnaSeqIntergenic.mean = zeros(ceil(chr.sequenceLen / this.ntResolution), 1);
            
            tmpRnaSeqNascent.var    = tmpRnaSeqNascent.mean;
            tmpRnaSeqProcessed.var  = tmpRnaSeqProcessed.mean;
            tmpRnaSeqIntergenic.var = tmpRnaSeqIntergenic.mean;
            
            for i = size(this.tmpRnaSeqNascentCnt)
                fivePrime = trn.transcriptionUnitFivePrimeCoordinates(i);
                dir = trn.transcriptionUnitDirections(i);
                len = trn.transcriptionUnitLengths(i);
                
                tmpRnaSeqNascent.mean = this.addDownsampled(tmpRnaSeqNascent.mean, ...
                    fivePrime, dir, len, this.tmpRnaSeqNascentCnt.mean(i));
                tmpRnaSeqNascent.var = this.addDownsampled(tmpRnaSeqNascent.var, ...
                    fivePrime, dir, len, this.tmpRnaSeqNascentCnt.var(i));
                
                pIdxs = find(rna.nascentRNAMatureRNAComposition(:, i));
                iIdxs = find(rna.intergenicRNAMatrix(:, i));
                if isscalar(pIdxs)
                    tmpRnaSeqProcessed.mean = this.addDownsampled(tmpRnaSeqProcessed.mean, ...
                        fivePrime, dir, len, this.tmpRnaSeqProcessedCnt.mean(pIdxs));
                    tmpRnaSeqProcessed.var = this.addDownsampled(tmpRnaSeqProcessed.var, ...
                        fivePrime, dir, len, this.tmpRnaSeqProcessedCnt.var(pIdxs));
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
                            tmpRnaSeqIntergenic.var = this.addDownsampled(tmpRnaSeqIntergenic.var, ...
                                fivePrime, dir, len, this.tmpRnaSeqIntergenicCnt.var(iIdxs(j - 1)));
                        end
                        
                        start = g.startCoordinates(gIdx);
                        len = g.lengths(gIdx);
                        tmpRnaSeqProcessed.mean = this.addDownsampled(tmpRnaSeqProcessed.mean, ...
                            fivePrime, dir, len, this.tmpRnaSeqProcessedCnt.mean(pIdxs(j)));
                        tmpRnaSeqProcessed.var = this.addDownsampled(tmpRnaSeqProcessed.var, ...
                            fivePrime, dir, len, this.tmpRnaSeqProcessedCnt.var(pIdxs(j)));
                    end
                end
            end
            
            tmpRnaSeqNascent.mean(   1:end-1) = tmpRnaSeqNascent.mean(   1:end-1) / this.ntResolution;
            tmpRnaSeqProcessed.mean( 1:end-1) = tmpRnaSeqProcessed.mean( 1:end-1) / this.ntResolution;
            tmpRnaSeqIntergenic.mean(1:end-1) = tmpRnaSeqIntergenic.mean(1:end-1) / this.ntResolution;
            tmpRnaSeqNascent.mean(   1:end)   = tmpRnaSeqNascent.mean(   1:end)   / mod(chr.sequenceLen, this.ntResolution);
            tmpRnaSeqProcessed.mean( 1:end)   = tmpRnaSeqProcessed.mean( 1:end)   / mod(chr.sequenceLen, this.ntResolution);
            tmpRnaSeqIntergenic.mean(1:end)   = tmpRnaSeqIntergenic.mean(1:end)   / mod(chr.sequenceLen, this.ntResolution);
            
            tmpRnaSeqNascent.var(   1:end-1) = tmpRnaSeqNascent.var(   1:end-1) / (this.ntResolution ^ 2);
            tmpRnaSeqProcessed.var( 1:end-1) = tmpRnaSeqProcessed.var( 1:end-1) / (this.ntResolution ^ 2);
            tmpRnaSeqIntergenic.var(1:end-1) = tmpRnaSeqIntergenic.var(1:end-1) / (this.ntResolution ^ 2);
            tmpRnaSeqNascent.var(   1:end)   = tmpRnaSeqNascent.var(   1:end)   / (mod(chr.sequenceLen, this.ntResolution) ^ 2);
            tmpRnaSeqProcessed.var( 1:end)   = tmpRnaSeqProcessed.var( 1:end)   / (mod(chr.sequenceLen, this.ntResolution) ^ 2);
            tmpRnaSeqIntergenic.var(1:end)   = tmpRnaSeqIntergenic.var(1:end)   / (mod(chr.sequenceLen, this.ntResolution) ^ 2);

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
                        
            this.tmpRnaSeqTranscripts   = [];
            this.tmpRnaSeqNascentCnt    = [];
            this.tmpRnaSeqProcessedCnt  = [];
            this.tmpRnaSeqIntergenicCnt = [];
            this.rnaSeq = rmfield(this.rnaSeq, {'sum', 'sum_x2'});
            
            %ChIP seq
            this.tmpChipSeqMon.mean   = 1 / nTime * this.tmpChipSeqMon.sum;
            this.tmpChipSeqMon.exp_x2 = 1 / nTime * this.tmpChipSeqMon.sum_x2;
            this.tmpChipSeqMon.var    = max(0, this.tmpChipSeqMon.exp_x2 - this.tmpChipSeqMon.mean .^ 2);
            
            this.tmpChipSeqCpx.mean   = 1 / nTime * this.tmpChipSeqCpx.sum;
            this.tmpChipSeqCpx.exp_x2 = 1 / nTime * this.tmpChipSeqCpx.sum_x2;
            this.tmpChipSeqCpx.var    = max(0, this.tmpChipSeqCpx.exp_x2 - this.tmpChipSeqCpx.mean .^ 2);
            
            pcComp = sum(pc.proteinComplexComposition(g.mRNAIndexs, :, :), 3)';
            this.chipSeq.mean = ...
                + this.avgFootprintDownsampled(this.tmpChipSeqMon.mean, this.idxsDnaBoundMon, chr.monomerDNAFootprints, this.ntResolution) ...
                + this.avgFootprintDownsampled(this.tmpChipSeqCpx.mean, this.idxsDnaBoundCpx, chr.complexDNAFootprints, this.ntResolution) * pcComp;
            this.chipSeq.std  = sqrt(...
                + this.avgFootprintDownsampled(this.tmpChipSeqMon.var, this.idxsDnaBoundMon, chr.monomerDNAFootprints, this.ntResolution) ...
                + this.avgFootprintDownsampled(this.tmpChipSeqCpx.var, this.idxsDnaBoundCpx, chr.complexDNAFootprints, this.ntResolution) * pcComp ...
                );
            
            this.tmpChipSeqMon = [];
            this.tmpChipSeqCpx = [];
            this.chipSeq = rmfield(this.chipSeq, {'sum', 'sum_x2'});
            
            %% save
            data = struct(...
                'singleCell', struct(...
                    'time',         this.time, ...
                    'growth',       this.growth, ...
                    'mass',         this.mass, ...
                    'volume',       this.volume, ...
                    'repInitTime',  this.repInitTime, ...
                    'repTermTime',  this.repTermTime, ...
                    'cellCycleLen', this.cellCycleLen ...
                    ), ...
                'metConcs',     this.metConcs, ...
                'dnaSeq',       this.dnaSeq, ...
                'rnaSeq',       this.rnaSeq, ...
                'chipSeq',      this.chipSeq, ...
                'rnaArray',     this.rnaArray, ...
                'protArray',    this.protArray, ...
                'rxnFluxes',    this.rxnFluxes, ...
                'labels',       this.labels ...
                ); %#ok<NASGU>
            save(this.simPath, '-struct', 'data');
        end
    end
    
    methods (Static = true)
        function [avgVals, labels] = average(simPathPattern, verbosity)
            import edu.stanford.covert.cell.sim.util.HighthroughputExperimentsLogger;
            
            if nargin < 2
                verbosity = 1;
            end
            
            %get matching files
            [inPathBase, ~, ~] = fileparts(simPathPattern);
            files = dir(simPathPattern);
            
            %get longest simulation length
            nTimeMax = 0;
            isFilesSim = false(size(files));
            for i = 1:numel(files)
                tmpPath = fullfile(inPathBase, files(i).name);
                vars = whos('-file', tmpPath);
                if ~isequal({vars.name}', {
                        'chipSeq'; 'dnaSeq'; 'labels'; 'metConcs'; 'protArray';
                        'rnaArray'; 'rnaSeq'; 'rxnFluxes'; 'singleCell'
                        })
                    continue;
                end
                tmp = load(tmpPath, 'singleCell');
                if size(tmp.singleCell.growth, 1) > 1
                    continue;
                end
                
                isFilesSim(i) = true;
                nTimeMax = max(nTimeMax, numel(tmp.singleCell.time));
            end
            files = files(isFilesSim);
            if verbosity >= 1
                fprintf('Averaging %d simulations ...', numel(files));
            end
            
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
            
            avgVals.singleCell = struct();
            avgVals.singleCell.time   = NaN(1, nTimeMax);
            avgVals.singleCell.growth = NaN(nSim, nTimeMax);
            avgVals.singleCell.mass   = NaN(nSim, nTimeMax);
            avgVals.singleCell.volume = NaN(nSim, nTimeMax);
            avgVals.singleCell.repInitTime  = NaN(nSim, 1);
            avgVals.singleCell.repTermTime  = NaN(nSim, 1);
            avgVals.singleCell.cellCycleLen = NaN(nSim, 1);
            
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
                nTime = numel(vals.singleCell.time);
                
                avgVals.singleCell.time(1, 1:nTime)   = vals.singleCell.time;
                avgVals.singleCell.growth(i, 1:nTime) = vals.singleCell.growth;
                avgVals.singleCell.mass(  i, 1:nTime) = vals.singleCell.mass;
                avgVals.singleCell.volume(i, 1:nTime) = vals.singleCell.volume;
                
                avgVals.singleCell.repInitTime( i, 1) = vals.singleCell.repInitTime;
                avgVals.singleCell.repTermTime( i, 1) = vals.singleCell.repTermTime;
                avgVals.singleCell.cellCycleLen(i, 1) = vals.singleCell.cellCycleLen;
                
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
            if verbosity >= 1
                fprintf(' done.\n');
            end
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
            
            arr(startSeg) = arr(startSeg) + val * (HighthroughputExperimentsLogger.ntResolution - mod(startCoor, HighthroughputExperimentsLogger.ntResolution) + 1);
            arr(endSeg)   = arr(endSeg)   + val *                                                 mod(endCoor,   HighthroughputExperimentsLogger.ntResolution);
            arr(startSeg+1:endSeg-1) = ...
                + arr(startSeg+1:endSeg-1) ...
                + HighthroughputExperimentsLogger.ntResolution * val;
        end
        
        function downSampArr = avgFootprintDownsampled(fullSampArr, protToBoundProtIdxs, footprints, downSampleSize)
            L = size(fullSampArr, 1);
            idxs = find(protToBoundProtIdxs);
            downSampArr = zeros(ceil(L / downSampleSize), numel(protToBoundProtIdxs));
            for i = 1:numel(idxs)
                if any(fullSampArr(:, i))
                    tmpFull = full([
                        fullSampArr(end-1000+1:end, i)
                        fullSampArr(:, i)
                        fullSampArr(1:1000, i)
                        ]);
                    tmpFullConv = conv(tmpFull, ones(footprints(idxs(i)), 1), 'same');
                    
                    tmpDown = sum(reshape([
                        tmpFullConv(1001:end-1000)
                        zeros(-mod(L, -downSampleSize), 1)
                        ], downSampleSize, []), 1)';
                    tmpDown(1:end-1) = tmpDown(1:end-1) / downSampleSize;
                    tmpDown(end)     = tmpDown(end)     / mod(L, downSampleSize);
                    
                    downSampArr(:, idxs(i)) = ...
                        + downSampArr(:, idxs(i)) ...
                        + tmpDown;
                end
            end
        end
    end
end
