%CellState
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef CellState
    %printing
    methods (Static)
        function run(sim, fileName)
            import edu.stanford.covert.cell.sim.analysis.CellState;
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            import edu.stanford.covert.cell.sim.util.PrintUtil;            
            
            %excel file
            [content, colLabels, indentation] = CellState.state(sim);
            if nargin == 1
                PrintUtil.printToStdIO(content, colLabels, struct('indentation', indentation));
            else                
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'State', struct('indentation', indentation));
            end
            
            [content, colLabels, indentation] = CellState.processes(sim);
            if nargin == 1
                PrintUtil.printToStdIO(content, colLabels, struct('indentation', indentation));
            else                
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'Process', struct('indentation', indentation));
            end
            
            % plots
            if nargin == 1
                CellState.plotStateOnChromosome(sim, PlotUtil.newAxesHandle());
            else
                [axesHandle, figHandle] = PlotUtil.newAxesHandle();
                
                cla(axesHandle);
                CellState.plotStateOnChromosome(sim, axesHandle);
                saveas(figHandle, [fileName '-OnChromosome.pdf']);
                
                close(figHandle);
            end
        end
        
        function [content, colLabels, indentation] = state(sim)
            %import classes
            import edu.stanford.covert.util.ConstantUtil;
            import edu.stanford.covert.cell.sim.Simulation;
                        
            content = cell(0, 4);
            colLabels = {'Quantity', 'Value', 'Units'};
            
            %% time
            time = sim.state('Time');
            content = [content;{
                0 'Time' time.values(1)/ConstantUtil.secondsPerHour 'h';
                }];
            
            %% mass
            mass = sim.state('Mass');
            cellDryMass = sum(mass.cellDry);
            content = [content;{
                0 'Mass' sum(mass.cell) 'g';
                    1 'Dry' sum(mass.cellDry) 'g';
                    1 'Metabolite' sum(mass.metaboliteWt) 'g';
                        2 'Percent Dry' sum(mass.metaboliteWt) / cellDryMass * 100 '%';
                    1 'DNA' sum(mass.dnaWt) 'g';
                        2 'Percent Dry' sum(mass.dnaWt) / cellDryMass * 100 '%';
                    1 'RNA' sum(mass.rnaWt) 'g';
                        2 'Percent Dry' sum(mass.rnaWt) / cellDryMass * 100 '%';
                    1 'Protein' sum(mass.proteinWt) 'g';
                        2 'Percent Dry' sum(mass.proteinWt) / cellDryMass * 100 '%';
                    1 'No. cells' sum(mass.cellDry)/mass.cellInitialDryWeight 'cells';
                }];
            
            %% growth
            s = sim.state('MetabolicReaction');
            time = sim.state('Time');
            content = [content;{
                0 'Growth' s.growth 'cells/s';
                    1 'Relative to expected rate' s.growth / (1 / time.cellCycleLength * log(2) * exp(time.values / time.cellCycleLength * log(2))) [];
                    1 'Doubling Time', s.doublingTime / ConstantUtil.secondsPerHour 'h';
                }];
            
            %% shape
            s = sim.state('Geometry');
            content = [content;{
                0 'Shape' 'Cyclinder with hemispherical caps' []
                    1 'Width' s.width 'm'
                    1 'Pinched Diameter' s.pinchedDiameter 'm'
                    1 'Cylindrical length' s.cylindricalLength 'm'
                    1 'Total length' s.totalLength 'm'
                    1 'Pinched circumference' s.pinchedCircumference 'm'
                    1 'Volume' s.volume 'L'
                    1 'Surface area' s.surfaceArea' 'm^2'
                    1 'Pinched' s.pinched 'true/false'
                }];
            
            %% FtsZ ring
            s = sim.state('FtsZRing');
            content = [content;{
                0 'FtsZ ring' [] []
                    1 'No. edges' s.numEdges []
                    1 'No. edges bound by 1 straight FtsZ 9-mer' s.numEdgesOneStraight []
                    1 'No. edges bound by 2 straight FtsZ 9-mers' s.numEdgesTwoStraight []
                    1 'No. edges bound by 2 bent FtsZ 9-mers' s.numEdgesTwoBent []
                    1 'No. edges bound by 1 bent FtsZ 9-mer' s.numResidualBent []
                }];
            
            %% stimuli
            content = [content;{
                0 'Stimuli' '' []
                }];
            stim = sim.state('Stimulus');
            for i = 1:numel(stim.wholeCellModelIDs)
                value = stim.values(i, stim.values(i, :) ~= 0);
                if numel(value) == 0
                    value = 0;
                elseif numel(value) > 1
                    value = mean(value);
                end
                content = [content;{
                    1 stim.wholeCellModelIDs{i} value []
                    }]; %#ok<AGROW>
            end
            
            %% metabolites
            c = sim.state('Chromosome');
            s = sim.state('Metabolite');
            
            freedNTPs = s.counts(s.dntpIndexs, sim.compartment.cytosolIndexs);
            freeNTPs = s.counts(s.ntpIndexs, sim.compartment.cytosolIndexs);
            freeNDPs = s.counts(s.ndpIndexs, sim.compartment.cytosolIndexs);
            freeNMPs = s.counts(s.nmpIndexs, sim.compartment.cytosolIndexs);
            freeAAs = s.counts(s.aminoAcidIndexs, sim.compartment.cytosolIndexs);
            
            [mindNTPVal, mindNTPIdx] = min(freedNTPs);
            [minNTPVal, minNTPIdx] = min(freeNTPs);
            [minNDPVal, minNDPIdx] = min(freeNDPs);
            [minNMPVal, minNMPIdx] = min(freeNMPs);
            [minAAVal, minAAIdx] = min(freeAAs);
            
            freeNTPMass = freeNTPs' * s.molecularWeights(s.ntpIndexs) / ConstantUtil.nAvogadro;
            freeAAMass = freeAAs' * s.molecularWeights(s.aminoAcidIndexs) / ConstantUtil.nAvogadro;
            
            content = [content;{
                0 'Metabolites' [] []
                    1 'Free' sum(s.counts(:)) 'molecules'
                    1 'dNTPs' [] []
                        2 'Free' sum(freedNTPs) 'molecules'                        
                        2 'Min free' s.wholeCellModelIDs{s.dntpIndexs(mindNTPIdx)} []
                        2 'Min free' mindNTPVal 'molecules'
                        2 'Incorporated into DNA' (1 - sum(freedNTPs) / (sum(freedNTPs) + collapse(c.polymerizedRegions)))*100 '%'
                    1 'NTPs' [] []
                        2 'Free' sum(freeNTPs) 'molecules'                        
                        2 'Min free' s.wholeCellModelIDs{s.ntpIndexs(minNTPIdx)} []
                        2 'Min free' minNTPVal 'molecules'
                        2 'Incorporated into RNA' (1 - freeNTPMass / (freeNTPMass + sum(sim.state('Mass').rnaWt)))*100 '%'
                    1 'NDPs' [] []
                        2 'Free' sum(freeNDPs) 'molecules'                       
                        2 'Min free' s.wholeCellModelIDs{s.ndpIndexs(minNDPIdx)} []
                        2 'Min free' minNDPVal 'molecules'
                    1 'NMPs' [] []
                        2 'Free' sum(freeNMPs) 'molecules'                        
                        2 'Min free' s.wholeCellModelIDs{s.nmpIndexs(minNMPIdx)} []
                        2 'Min free' minNMPVal 'molecules'
                    1 'Amino acids' [] []
                        2 'Free' sum(freeAAs) 'molecules'                        
                        2 'Min free' s.wholeCellModelIDs{s.aminoAcidIndexs(minAAIdx)} []
                        2 'Min free' minAAVal 'molecules'
                        2 'Incorporated into Protein' (1 - freeAAMass / (freeAAMass + sum(sim.state('Mass').proteinWt)))*100 '%'
                    1 'Energy' [] []
                        2 'ATP' sum(s.counts(s.getIndexs({'ATP'}), :)) 'molecules'
                        2 'ADP' sum(s.counts(s.getIndexs({'ADP'}), :)) 'molecules'
                        2 'AMP' sum(s.counts(s.getIndexs({'AMP'}), :)) 'molecules'
                        2 'PPi' sum(s.counts(s.getIndexs({'PPI'}), :)) 'molecules'
                        2 'Pi'  sum(s.counts(s.getIndexs({'PI'}),  :)) 'molecules'
                        2 'H2O' sum(s.counts(s.getIndexs({'H2O'}), :)) 'molecules'
                        2 'H'   sum(s.counts(s.getIndexs({'H'}),   :)) 'molecules'
                }];
            
            %% DNA
            c = sim.state('Chromosome');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            ri = sim.process('ReplicationInitiation');
            dr = sim.process('DNARepair');
            dnaABoxBound = ...
                ri.calculateDnaABoxStatus() == ri.dnaABoxStatus_DnaAATPBound | ...
                ri.calculateDnaABoxStatus() == ri.dnaABoxStatus_DnaAADPBound;
            [polATPs, polADPs] = ri.calculateDnaAR1234Polymerization();
            oriCComplexSize = ri.calculateDnaAR1234ComplexSize(polATPs, polADPs);
            [~, boundMonomers] = find(c.monomerBoundSites);
            [~, boundComplexs] = find(c.complexBoundSites);
            occupiedBases = sum(c.monomerDNAFootprints(boundMonomers) + c.complexDNAFootprints(boundComplexs));
                        
            methylatedRMSites = nnz(c.damagedBases([
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1))   ones(size(dr.RM_MunI_RecognitionSites, 1), 1);
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2)) 2*ones(size(dr.RM_MunI_RecognitionSites, 1), 1);
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(1)) 3*ones(size(dr.RM_MunI_RecognitionSites, 1), 1);
                dr.RM_MunI_RecognitionSites(:, dr.RM_MunI_MethylatedPositions(2)) 4*ones(size(dr.RM_MunI_RecognitionSites, 1), 1);
                ])) / 2;
            polymerizedRMSites = sum(c.isRegionPolymerized([
                dr.RM_MunI_RecognitionSites(:, 1)   ones(size(dr.RM_MunI_RecognitionSites, 1), 1);
                dr.RM_MunI_RecognitionSites(:, 1) 2*ones(size(dr.RM_MunI_RecognitionSites, 1), 1);
                dr.RM_MunI_RecognitionSites(:, 1) 3*ones(size(dr.RM_MunI_RecognitionSites, 1), 1);
                dr.RM_MunI_RecognitionSites(:, 1) 4*ones(size(dr.RM_MunI_RecognitionSites, 1), 1);
                ], 6, false)) / 2;
            
            content = [content;{
                0 'No. chromosomes' collapse(c.polymerizedRegions)/c.sequenceLen [];
                    1 'Bases' collapse(c.polymerizedRegions) 'nt';
                    1 'Linking number' collapse(c.linkingNumbers)/2 'Links';
                    1 'Superhelical density', (collapse(c.linkingNumbers)  - collapse(c.linkingNumbers_minFreeEnergy)) / collapse(c.linkingNumbers_minFreeEnergy) [];
                    1 'Bound proteins' nnz(c.monomerBoundSites) + nnz(c.complexBoundSites) 'molecules'
                        2 'Monomers' nnz(c.monomerBoundSites) 'molecules'
                            3 'GntR' nnz(c.monomerBoundSites == pm.getIndexs({'MG_101_MONOMER'})) 'molecules'
                        2 'Complexes' nnz(c.monomerBoundSites) 'molecules'
                            3 'DnaA' nnz(dnaABoxBound) 'molecules'
                                4 '7-mer' nnz(dnaABoxBound(ri.dnaABoxIndexs_7mer, :)) 'molecules'
                                4 '8-mer' nnz(dnaABoxBound(ri.dnaABoxIndexs_8mer, :)) 'molecules'
                                4 '9-mer' nnz(dnaABoxBound(ri.dnaABoxIndexs_9mer, :)) 'molecules'
                                4 'OriC Complex' oriCComplexSize(1) 'molecules'
                                4 'R1' nnz(dnaABoxBound(ri.dnaABoxIndexs_R12345(1), :)) 'molecules'
                                4 'R2' nnz(dnaABoxBound(ri.dnaABoxIndexs_R12345(2), :)) 'molecules'
                                4 'R3' nnz(dnaABoxBound(ri.dnaABoxIndexs_R12345(3), :)) 'molecules'
                                4 'R4' nnz(dnaABoxBound(ri.dnaABoxIndexs_R12345(3), :)) 'molecules'
                                4 'R5' nnz(dnaABoxBound(ri.dnaABoxIndexs_R12345(5), :)) 'molecules'
                            3 'Gyrase' nnz(c.complexBoundSites == pc.getIndexs({'DNA_GYRASE'})) 'molecules'
                            3 'LuxR' nnz(c.complexBoundSites == pc.getIndexs({'MG_428_DIMER'})) 'molecules'
                            3 'SMC' nnz(c.complexBoundSites == pc.getIndexs({'MG_213_214_298_6MER_ADP'})) 'molecules'
                                4 'Average spacing' collapse(c.polymerizedRegions)/2 / nnz(c.complexBoundSites == pc.getIndexs({'MG_213_214_298_6MER_ADP'})) 'nt'
                        2 'Occupied' occupiedBases 'nt'
                        2 'Occupied' occupiedBases / collapse(c.polymerizedRegions) / 2 * 100 '%'
                    1 'Damage' [] []
                        2 'Gap sites' nnz(c.gapSites) 'nt'
                        2 'Abasic sites' nnz(c.abasicSites) 'nt'
                        2 'Damaged bases' nnz(c.damagedBases) 'nt'
                            3 'm6AD' nnz(c.damagedBases == c.metabolite.m6ADIndexs) 'nt'
                            3 'methylated R/M sites' methylatedRMSites 'nt'
                            3 'methylated R/M sites' methylatedRMSites/polymerizedRMSites*100 '%'
                        2 'Damaged sugar-phosphate' nnz(c.damagedSugarPhosphates) 'nt'
                        2 'Strand breaks' nnz(c.strandBreaks) 'nt'
                        2 'Intrastrand crosslinks' nnz(c.intrastrandCrossLinks) 'nt'
                        2 'Holliday junctions' nnz(c.hollidayJunctions) 'nt'
                }];
            
            %% RNA
            r = sim.state('Rna');
            
            matureRNAs = sum(r.counts(r.matureIndexs, :), 2);
            aminoacylatedRNAs = sum(r.counts(r.aminoacylatedIndexs, :), 2);
            content = [content;{
                0 'RNA' sum(r.counts(:)) 'molecules'
                    1 'mRNA' sum(matureRNAs(r.matureMRNAIndexs)) 'molecules'
                        2 'No. species expressed' nnz(r.counts(r.matureIndexs(r.matureMRNAIndexs), :)) []
                    1 'rRNA' sum(matureRNAs(r.matureRRNAIndexs)) 'molecules'
                        2 '5S'  matureRNAs(r.matureRibosomalRRNAIndexs(1)) 'molecules'
                        2 '15S' matureRNAs(r.matureRibosomalRRNAIndexs(2)) 'molecules'
                        2 '23S' matureRNAs(r.matureRibosomalRRNAIndexs(3)) 'molecules'
                    1 'sRNA' sum(matureRNAs(r.matureSRNAIndexs) + aminoacylatedRNAs(r.matureSRNAIndexs)) 'molecules'
                        2 'scRNA'   matureRNAs(r.matureSRNAIndexs(1))+aminoacylatedRNAs(r.matureSRNAIndexs(1)) 'molecules'
                        2 'MG170'   matureRNAs(r.matureSRNAIndexs(2))+aminoacylatedRNAs(r.matureSRNAIndexs(2)) 'molecules'
                        2 'RNase P' matureRNAs(r.matureSRNAIndexs(3))+aminoacylatedRNAs(r.matureSRNAIndexs(3)) 'molecules'
                        2 'tmRNA'   matureRNAs(r.matureSRNAIndexs(4))+aminoacylatedRNAs(r.matureSRNAIndexs(4)) 'molecules'
                            3 'Free' matureRNAs(r.matureSRNAIndexs(4)) 'molecules'
                            4 'Aminoacylated' aminoacylatedRNAs(r.matureSRNAIndexs(4)) 'molecules'
                    1 'tRNA' sum(matureRNAs(r.matureTRNAIndexs) + aminoacylatedRNAs(r.matureTRNAIndexs)) 'molecules'
                        2 'Free' sum(matureRNAs(r.matureTRNAIndexs)) 'molecules'
                        2 'Aminoacylated' sum(aminoacylatedRNAs(r.matureTRNAIndexs)) 'molecules'
                        2 'No. species expressed' nnz(matureRNAs(r.matureTRNAIndexs) + aminoacylatedRNAs(r.matureTRNAIndexs)) 'molecules'
                    1 'Immature' nnz(sim.process('Transcription').transcripts.boundTranscriptionUnits) 'molecules'
                    1 'Aborted' size(sim.state('Transcript').abortedTranscripts, 1) 'molecules'
                }];
            
            %% Protein
            m = sim.state('ProteinMonomer');
            c = sim.state('ProteinComplex');
            
            ribosomalMonomerIdxs = sum(c.proteinComplexComposition(sim.gene.mRNAIndexs, c.getIndexs({'RIBOSOME_70S'}), :),3);
            ribosomalMonomers = sum(m.counts(m.matureIndexs(ribosomalMonomerIdxs > 0), :), 2);
            [minRibosomalMonomers, minRibosomalMonomerIdx] = min(ribosomalMonomers);
            minRibosomalMonomerID = m.wholeCellModelIDs{m.matureIndexs(minRibosomalMonomerIdx)};
            
            [idxs, ~, counts] = find(sum(m.counts(m.inactivatedIndexs, :), 2));
            inactiveMonomers = cell(numel(idxs), 4);
            inactiveMonomers(:, 1) = {2};
            inactiveMonomers(:, 2) = m.wholeCellModelIDs(m.inactivatedIndexs(idxs));
            inactiveMonomers(:, 3) = num2cell(counts);

            [idxs, ~, counts] = find(sum(c.counts(c.inactivatedIndexs, :), 2));
            inactiveComplexs = cell(numel(idxs), 4);
            inactiveComplexs(:, 1) = {2};
            inactiveComplexs(:, 2) = c.wholeCellModelIDs(c.inactivatedIndexs(idxs));
            inactiveComplexs(:, 3) = reshape(num2cell(counts), [], 1);
            
            content = [content;{
                0 'Protein' sum(m.counts(:)) + sum(c.counts(:)) 'molecules'
                    1 'Monomers' sum(m.counts(:)) 'molecules'
                        2 'Ribosomal proteins' sum(ribosomalMonomers) 'molecules'
                            3 'Min expressed species' minRibosomalMonomers 'molecules'
                            3 'Min expressed species' minRibosomalMonomerID []
                    1 'Complexs' sum(c.counts(:)) 'molecules'
                    1 'Bound' sum(sum(m.counts(m.boundIndexs, :)))+sum(sum(c.counts(c.boundIndexs, :))) 'molecules'
                        2 'Chromosome' nnz(sim.state('Chromosome').monomerBoundSites)+nnz(sim.state('Chromosome').complexBoundSites) 'molecules'
                        2 'RNA polymerase' nnz(sim.process('Transcription').transcripts.boundTranscriptionUnits) 'molecules'
                        2 'Ribosome' nnz(sim.state('Ribosome').boundMRNAs) 'molecules'
                    1 'Immature' nnz(sim.state('Ribosome').boundMRNAs) 'molecules'
                    1 'Inactive' sum(sum(m.counts(m.inactivatedIndexs, :)))+sum(sum(c.counts(c.inactivatedIndexs, :))) 'molecules'
                        };
                        inactiveMonomers;
                        inactiveComplexs;
                        {
                    1 'Damaged' sum(sum(m.counts(m.damagedIndexs, :)))+sum(sum(c.counts(c.damagedIndexs, :))) 'molecules'
                    1 'Aborted' size(sim.state('Polypeptide').abortedPolypeptides, 1) 'molecules'
                    1 'Cytosol' sum(m.counts(:, sim.compartment.cytosolIndexs))+ sum(c.counts(:, sim.compartment.cytosolIndexs)) 'molecules'
                    1 'Membrane' sum(m.counts(:, sim.compartment.membraneIndexs))+ sum(c.counts(:, sim.compartment.membraneIndexs)) 'molecules'
                    1 'Term. Org. Cytosol' sum(m.counts(:, sim.compartment.terminalOrganelleCytosolIndexs))+ sum(c.counts(:, sim.compartment.terminalOrganelleCytosolIndexs)) 'molecules'
                    1 'Term. Org. Membrane' sum(m.counts(:, sim.compartment.terminalOrganelleMembraneIndexs))+ sum(c.counts(:, sim.compartment.terminalOrganelleMembraneIndexs)) 'molecules'
                    1 'Extracellular' sum(m.counts(:, sim.compartment.extracellularIndexs))+ sum(c.counts(:, sim.compartment.extracellularIndexs)) 'molecules'
                }];        
            
            %% DNA Polymerase
            r = sim.process('Replication');
            free = r.enzymes(r.enzymeIndexs_core) + r.enzymeComposition(r.enzymeIndexs_core, :) * r.enzymes;
            bound = r.boundEnzymes(r.enzymeIndexs_core) + r.enzymeComposition(r.enzymeIndexs_core, :) * r.boundEnzymes;            
            content = [content;{
                0 'DNA polymerase' free+bound 'molecules'
                    1 'Free' free 'molecules'
                    1 'Bound' bound 'molecules'
                }];
            
            %% Transcription enzymes
            t = sim.process('Transcription');
            s = sim.state('RNAPolymerase');
            rna = t.rnaPolymerases;
            rnaPolymeraseStateOccupancies = s.stateOccupancies*100;
            content = [content;{
                0 'Transcription enzymes' [] []                
                    1 'RNA polymerase' t.enzymes(t.enzymeIndexs_rnaPolymerase)+t.boundEnzymes(t.enzymeIndexs_rnaPolymerase) 'molecules'
                        2 'Active' rnaPolymeraseStateOccupancies(rna.activelyTranscribingIndex) ' %'
                        2 'Specifically bound' rnaPolymeraseStateOccupancies(rna.specificallyBoundIndex) ' %'
                        2 'Non-specifically bound' rnaPolymeraseStateOccupancies(rna.nonSpecificallyBoundIndex) ' %'
                        2 'Free' rnaPolymeraseStateOccupancies(rna.freeIndex) ' %'
                    1 'Initiation/Elongation factors' sum(t.enzymes(t.enzymeIndexs_transcriptionFactors)+t.boundEnzymes(t.enzymeIndexs_transcriptionFactors)) 'molecules'
                        2 'Free' sum(t.enzymes(t.enzymeIndexs_transcriptionFactors)) 'molecules'
                        2 'Bound' sum(t.boundEnzymes(t.enzymeIndexs_transcriptionFactors)) 'molecules'
                }];        
            
            %% tRNA aminoacylation enzymes
            t = sim.process('tRNAAminoacylation');
            totEnzymes = t.enzymes + t.boundEnzymes;
            content = [content;{
                0 'tRNA aminoacylation enzymes' sum(totEnzymes) 'molecules'
                    1 'tRNA synthetases' sum(totEnzymes(t.enzymeIndexs_tRNASynthetases)) 'molecules'
                    1 'tRNA transferases' sum(totEnzymes(t.enzymeIndexs_tRNATransferases)) 'molecules'
                }];
                                    
            %% Translation enzymes
            t = sim.process('Translation');
            r = sim.state('Ribosome');
            totEnzymes = t.enzymes + t.boundEnzymes;
            content = [content;{
                0 'Translation enzymes' [] []
                    1 'Ribosome' [] []
                        2 'Ribosome 30S' totEnzymes(t.enzymeIndexs_ribosome30S) 'molecules'
                        2 'Ribosome 30S-IF3' totEnzymes(t.enzymeIndexs_ribosome30SIF3) 'molecules'
                        2 'Ribosome 50S' totEnzymes(t.enzymeIndexs_ribosome50S) 'molecules'
                        2 'Ribosome 70S' totEnzymes(t.enzymeIndexs_ribosome70S) 'molecules'
                            3 'Free' t.enzymes(t.enzymeIndexs_ribosome70S) 'molecules'
                            3 'Bound' t.boundEnzymes(t.enzymeIndexs_ribosome70S) 'molecules'
                            3 'Active' nnz(r.boundMRNAs) 'molecules'
                    1 'Translation factors' sum(totEnzymes(t.enzymeIndexs_translationFactors)) 'molecules'
                        1 'Free' sum(t.enzymes(t.enzymeIndexs_translationFactors)) 'molecules'
                        1 'Bound' sum(t.boundEnzymes(t.enzymeIndexs_translationFactors)) 'molecules'
                }];                       
            
            %% FtsZ
            f = sim.process('FtsZPolymerization');
            totEnzymes = f.enzymes + f.boundEnzymes;
            content = [content;{
                0 'FtsZ' [1 1 1 2:9] * totEnzymes 'molecules'
                    1 'Free monomers' totEnzymes(f.enzymeIndexs_FtsZ) 'molecules'
                    1 'GDP 1-mer' totEnzymes(f.enzymeIndexs_FtsZ_GDP) 'molecules'
                    1 'GTP 1-mer' totEnzymes(f.enzymeIndexs_FtsZ_GTP) 'molecules'
                    1 'GTP 2-mer' totEnzymes(f.enzymeIndexs_FtsZ_activated(2)) 'molecules'
                    1 'GTP 3-mer' totEnzymes(f.enzymeIndexs_FtsZ_activated(3)) 'molecules'
                    1 'GTP 4-mer' totEnzymes(f.enzymeIndexs_FtsZ_activated(4)) 'molecules'
                    1 'GTP 5-mer' totEnzymes(f.enzymeIndexs_FtsZ_activated(5)) 'molecules'
                    1 'GTP 6-mer' totEnzymes(f.enzymeIndexs_FtsZ_activated(6)) 'molecules'
                    1 'GTP 7-mer' totEnzymes(f.enzymeIndexs_FtsZ_activated(7)) 'molecules'
                    1 'GTP 8-mer' totEnzymes(f.enzymeIndexs_FtsZ_activated(8)) 'molecules'
                    1 'GTP 9-mer' totEnzymes(f.enzymeIndexs_FtsZ_activated(9)) 'molecules'
                    1 'Bound 9-mer' f.boundEnzymes(f.enzymeIndexs_FtsZ_activated(9)) 'molecules'
                }];            
            
            %% format output
            indentation = cell2mat(content(:, 1));
            content = content(:, 2:end);
        end
            
        %diagnostics
        function [content, colLabels, indentation] = processes(sim)
            content = cell(0, 4);
            colLabels = {'Quantity', 'Value', 'Units'};
            
            %% Transcription
            t = sim.process('Transcription');
            s = sim.state('RNAPolymerase');
            
            %building blocks
            availableNTPs = t.substrates(t.substrateIndexs_ntp);
            rnaPolymeraseMaxNTPIncorporation = s.nActive * t.rnaPolymeraseElongationRate;
            transcriptionFactors = t.enzymes(t.enzymeIndexs_transcriptionFactors);
            
            %limits
            transcriptionNTPLimit = min(availableNTPs, ...
                availableNTPs / sum(availableNTPs) * rnaPolymeraseMaxNTPIncorporation);
            if abs(rnaPolymeraseMaxNTPIncorporation - sum(transcriptionNTPLimit)) < 1e-6
                transcriptionLimitingFractor = 'Active RNA Polymerases';
            else
                transcriptionLimitingFractor = 'Available NMPs, NTPs';
            end
            
            %limiting blocks
            [minNTPs, minNTP] = min(availableNTPs);
            minNTP = t.substrateWholeCellModelIDs{t.substrateIndexs_ntp(minNTP)};
            
            [minTranscriptionFactors, minTranscriptionFactor] = min(transcriptionFactors);
            minTranscriptionFactor = t.enzymeWholeCellModelIDs{t.enzymeIndexs_transcriptionFactors(minTranscriptionFactor)};
            
            %print
            content = [content;{
                0 'Transcription' [] []
                    1 'Allocated ATP' t.substrates(t.substrateIndexs_atp) 'molecules'
                    1 'Available NTPs' sum(availableNTPs) 'molecules'
                    1 'Min NTP' minNTPs 'molecules'
                    1 'Min NTP' minNTP []
                    1 'Active RNA Polymerases' s.nActive 'molecules'
                    1 'RNA Polymerase Elongation Rate' t.rnaPolymeraseElongationRate 'nt/s'
                    1 'Max RNA Polymerase NTP Incorporation'  rnaPolymeraseMaxNTPIncorporation 'nt/s'
                    1 'Transcription factors' sum(transcriptionFactors) 'molecules'
                    1 'Min transcription factors' minTranscriptionFactors 'molecules'
                    1 'Min transcription factors' minTranscriptionFactor []
                    1 'Max NTP Incorporation' sum(transcriptionNTPLimit) 'nt/s'
                    1 'Limiting Factor' transcriptionLimitingFractor []
             }];
            
            %% Translation
            t = sim.process('Translation');
            r = sim.state('Ribosome');
            
            %building blocks
            energy = t.substrates(t.substrateIndexs_gtp);
            
            aminoAcids = t.substrates(t.substrateIndexs_aminoAcids);
            ribosomeMaxAAIncorporation = r.nActive * t.ribosomeElongationRate;
            
            tRNAs = t.freeTRNAs + t.aminoacylatedTRNAs;
            translationFactors = t.enzymes(t.enzymeIndexs_translationFactors) + t.boundEnzymes(t.enzymeIndexs_translationFactors);
            
            %limiting blocks
            [minAminoAcids, minAminoAcid] = min(aminoAcids);
            minAminoAcid = t.substrateWholeCellModelIDs{t.substrateIndexs_aminoAcids(minAminoAcid)};
            
            [minTRNAs, minTRNA] = min(tRNAs);
            minTRNA = t.freeTRNAWholeCellModelIDs{minTRNA};
            
            [minTranslationFactors, minTranslationFactor]=min(translationFactors);
            minTranslationFactor = t.enzymeWholeCellModelIDs{t.enzymeIndexs_translationFactors(minTranslationFactor)};
            
            %limits            
            content = [content;{
                0 'Translation' [] []
                    1 'Allocated ATP' energy 'molecules'
                    1 'Available AAs'  sum(aminoAcids) 'molecules'
                    1 'Min Amino Acid' minAminoAcids 'molecules'
                    1 'Min Amino Acid' minAminoAcid []
                    1 'tRNAs' sum(tRNAs) 'molecules'
                    1 'Min tRNAs' minTRNAs 'molecules'
                    1 'Min tRNAs' minTRNA []
                    1 'Active Ribosomes' r.nActive 'molecules'
                    1 'Ribosome Elongation Rate' t.ribosomeElongationRate 'aa/s'
                    1 'Max Ribosome AA Incorporation' ribosomeMaxAAIncorporation 'aa/s'
                    1 'Translation factors' sum(translationFactors) 'molecules'
                    1 'Min translation factors' minTranslationFactors 'molecules'
                    1 'Min translation factors' minTranslationFactor []
                }];
            
            %% format output
            indentation = cell2mat(content(:, 1));
            content = content(:, 2:end);
        end
    end
    
    %plotting
    methods (Static)
        function plotStateOnChromosome(sim, axesHandle)
            c = sim.state('Chromosome');
            
            cla(axesHandle);
            hold on;
            legendHandles = [];
            
            %base chromosome
            for i = 1:4
                [r, xoffset] = calcRadiusXOffset(i);
                h = plot(axesHandle, xoffset + r*cos(pi/180 * (1:360)), r*sin(pi/180 * (1:360)), 'Color', [0.7 0.7 0.7], 'LineStyle', ':');
            end
            legendHandles = [legendHandles; h];
            
            %polymerized regions
            for i = 1:4
                [pos, lens] = find(c.polymerizedRegions(:, i));
                [r, xoffset] = calcRadiusXOffset(i);
                for j = 1:size(pos, 1)
                    theta = 2*pi * (pos(j,1) + (0:lens(j)/360:lens(j)))/c.sequenceLen;
                    h = plot(axesHandle, xoffset + r*cos(theta), r*sin(theta), 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
                end
            end
            legendHandles = [legendHandles; h];
            
            %polymerized transcription units
            for i = 1:2
                for j = 1:numel(c.transcriptionUnitWholeCellModelIDs)
                    if ~c.polymerizedTranscriptionUnits(j, i); continue; end;
                    [r, xoffset] = calcRadiusXOffset(2*(i-1) + c.transcriptionUnitStrands(j));
                    theta = 2*pi * (c.transcriptionUnitStartCoordinates(j) + (0:c.transcriptionUnitLengths(j)/360:c.transcriptionUnitLengths(j))) / c.sequenceLen;
                    h = plot(axesHandle, xoffset + r*cos(theta), r*sin(theta), 'Color', [0 1 0], 'LineWidth', 2);
                end                
            end
            legendHandles = [legendHandles; h];
            
            %polymerized genes
            for i = 1:2
                for j = 1:numel(sim.gene.wholeCellModelIDs)
                    if ~c.polymerizedGenes(j, i); continue; end;
                    [r, xoffset] = calcRadiusXOffset(2*(i-1) + sim.gene.strands(j));
                    theta = 2*pi * (sim.gene.startCoordinates(j) + (0:sim.gene.lengths(j)/360:sim.gene.lengths(j))) / c.sequenceLen;
                    h = plot(axesHandle, xoffset + r*cos(theta), r*sin(theta), 'Color', [0 0 1], 'LineWidth', 2);
                end
            end
            legendHandles = [legendHandles; h];
            
            %bound proteins
            [monomerPosStrnds, monomers] = find(c.monomerBoundSites);
            [complexPosStrnds, complexs] = find(c.complexBoundSites);
            posStrnds = [monomerPosStrnds; complexPosStrnds];
            lens = [c.monomerDNAFootprints(monomers); c.complexDNAFootprints(complexs)];
            ssStranded = [c.monomerDNAFootprintRegionStrandedness(monomers); c.complexDNAFootprintRegionStrandedness(complexs)] == c.dnaStrandedness_ssDNA;
            for i = 1:size(posStrnds, 1)
                [r, xoffset] = calcRadiusXOffset(posStrnds(i,2));
                theta = 2*pi * (posStrnds(i,1) + (0:lens(i)/360:lens(i))) / c.sequenceLen;
                h = plot(axesHandle, xoffset + r*cos(theta), r*sin(theta), 'Color', [1 1 0], 'LineWidth', 4);
                if ~ssStranded(i)
                    [r, xoffset] = calcRadiusXOffset(posStrnds(i, 2) -iseven(posStrnds(i, 2)) + isodd(posStrnds(i, 2)));
                    plot(axesHandle, xoffset + r*cos(theta), r*sin(theta), 'Color', [1 1 0], 'LineWidth', 4);
                end
            end
            legendHandles = [legendHandles; h];
            
            %damage
            posStrnds = find(c.getDamagedSites(1,1,1,1,1,1,1));
            for i = 1:size(posStrnds, 1)
                [r, xoffset] = calcRadiusXOffset(posStrnds(i,2));
                theta = 2*pi * posStrnds(i,1) / c.sequenceLen;
                h = plot(axesHandle, xoffset + r*cos(theta), r*sin(theta), 'Color', [1 0 0], 'LineWidth', 4);
            end
            legendHandles = [legendHandles; h];
            
            %expressed RNAs
            r = sim.state('Rna');
            genes = r.matureRNAGeneComposition * sum(r.counts(r.matureIndexs, :) + r.counts(r.aminoacylatedIndexs, :), 2);
            for i = 1:numel(sim.gene.wholeCellModelIDs)
                if genes(i) == 0; continue; end;
                
                r = 1 + [0.05;0.25] * iseven(sim.gene.strands(i)) - [0.05;0.25] * isodd(sim.gene.strands(i));
                xoffset = 0;
                
                theta = 2*pi * (sim.gene.startCoordinates(i) + sim.gene.lengths(i)/2) / c.sequenceLen;
                h = plot(axesHandle, xoffset + r*cos(theta), r*sin(theta), 'Color', [1 0 1], 'LineWidth', 2);
            end
            legendHandles = [legendHandles; h];
            
            %expressed proteins
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
             
            monomers = ...
                sum(pm.counts(pm.matureIndexs, :) + pm.counts(pm.boundIndexs, :) + pm.counts(pm.inactivatedIndexs, :), 2) + ...
                sum(pc.proteinComplexComposition(sim.gene.mRNAIndexs, :, :), 3) * sum(pc.counts(pc.matureIndexs, :) + pc.counts(pc.boundIndexs, :) + pc.counts(pc.inactivatedIndexs, :), 2);
            for i = 1:numel(monomers)
                if monomers(i) == 0; continue; end;
                
                r = 1 + [0.05;0.25] * iseven(sim.gene.strands(sim.gene.mRNAIndexs(i))) - [0.05;0.25] * isodd(sim.gene.strands(sim.gene.mRNAIndexs(i)));
                xoffset = 0;
                
                theta = 2*pi * (sim.gene.startCoordinates(sim.gene.mRNAIndexs(i)) + sim.gene.lengths(sim.gene.mRNAIndexs(i))/2) / c.sequenceLen;
                h = plot(axesHandle, xoffset + r*cos(theta), r*sin(theta), 'Color', [0 1 1], 'LineWidth', 2);
            end
            legendHandles = [legendHandles; h];
            
            %ribosomes
            ribosome = sim.state('Ribosome');
            for i = 1:numel(ribosome.boundMRNAs)
                if ribosome.boundMRNAs(i) == 0; continue; end;
                
                r = 1 + [0.05;0.25] * iseven(sim.gene.strands(sim.gene.mRNAIndexs(ribosome.boundMRNAs(i)))) - [0.05;0.25] * isodd(sim.gene.strands(sim.gene.mRNAIndexs(ribosome.boundMRNAs(i))));
                xoffset = 0;
                
                theta = 2*pi * (sim.gene.startCoordinates(sim.gene.mRNAIndexs(ribosome.boundMRNAs(i))) + ribosome.states(i)) / c.sequenceLen;
                plot(axesHandle, xoffset + r*cos(theta), r*sin(theta), 'Color', [1 1 0], 'LineWidth', 2);
            end

            %legend            
            legend(legendHandles, {
                'Chromosome'
                'Polymerized'
                '  Gene'
                '  TU'
                'Bnd Prot'
                'Damage'
                'Exp RNA'
                'Exp Prot'}, ...
                'Orientation', 'Horizontal', ...
                'Location', 'SouthOutside');
            
            %format            
            xlim([-1.4 2.5+1.4]);
            ylim([-1.4 1.4]);
            
            set(axesHandle, 'DataAspectRatio', [1 1 1])
            set(axesHandle, 'Visible', 'off')
            set(get(axesHandle, 'Parent'), 'Color', [1 1 1]);
            position = get(gca, 'Position');
            scale = 1.35;
            set(gca, 'Position', [position(1:2)-position(3:4)*(scale-1)/2 position(3:4)*scale]);
            
            function [r, xoffset] = calcRadiusXOffset(strnd)
                r = 1 + 0.05 * iseven(strnd) - 0.05 * isodd(strnd);
                xoffset = 2.5 * floor((strnd-1)/2);
            end
        end
    end
end