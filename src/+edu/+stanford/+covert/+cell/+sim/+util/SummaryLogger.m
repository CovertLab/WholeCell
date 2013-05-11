%Summary simulation logger. Implements simulation logger interface.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/10/2011
classdef SummaryLogger < edu.stanford.covert.cell.sim.util.Logger
    %options
    properties (Access = protected)
        stepSizeSec
        verbosity
        outputDirectory
    end
    
    %stored data
    properties (SetAccess = protected)
        log
        figureHandle
    end
    
    properties (Access = protected)
        timer
    end
    
    %indices into simulation states, processes
    properties (Access = protected)
        stateIndex_time
        stateIndex_mass
        stateIndex_rna
        stateIndex_metabolite
        stateIndex_monomer
        stateIndex_complex
        stateIndex_metabolicReaction
        stateIndex_chromosome
        stateIndex_geometry
        stateIndex_ftsZRing
        
		processIndex_transcription
		processIndex_translation
        processIndex_replicationInitiation
        processIndex_replication
        processIndex_ftsZPolymerization
    end
    
    properties(Constant = true)
        SIM_STATUS_DIDNT_START = -3;
        SIM_STATUS_STILL_RUNNING = -2;
        SIM_STATUS_TERMINATED_WITH_ERROR = -1;
        SIM_STATUS_COMPLETED_WITH_DIVISION = 0;
        SIM_STATUS_COMPLETED_WITHOUT_DIVISION = 1;
        
        SIM_STATUS_INDEX_STATUS = 1;
        SIM_STATUS_INDEX_REP_INIT_TIME = 2;
        SIM_STATUS_INDEX_INIT_REP_TIME = 3;
        SIM_STATUS_INDEX_REP_TIME = 4;
        SIM_STATUS_INDEX_CYTOKINESIS_TIME = 5;
        SIM_STATUS_INDEX_INIT_GROWTH = 6;
        SIM_STATUS_INDEX_FIN_GROWTH = 7;
        SIM_STATUS_INDEX_CUM_GROWTH = 8;
        SIM_STATUS_INDEX_FIN_MASS = 9;
        SIM_STATUS_INDEX_INIT_MASS = 10;
        SIM_STATUS_INDEX_MASS_DOUBLING_TIME = 11;
    end
    
    %ids
    properties
        ntpWholeCellModelIDs
        aaWholeCellModelIDs
    end
    
    methods
        function this = SummaryLogger(stepSizeSec, verbosity, outputDirectory)
            if nargin >= 1
                this.stepSizeSec = stepSizeSec;
            else
                this.stepSizeSec = 1;
            end
            
            if nargin >= 2
                this.verbosity = verbosity;
            else
                this.verbosity = 2;
            end
            
            if nargin >= 3
                this.outputDirectory = outputDirectory;
            else
                this.outputDirectory = [];
            end
        end
    end
    
    methods
        function setOptions(this, options)
            metaClass = metaclass(this);
            fields = intersect(fieldnames(options), cellfun(@(x) x.Name, metaClass.Properties, 'UniformOutput',false));
            for i = 1:numel(fields)
                this.(fields{i}) = options.(fields{i});
            end
        end
        
        %Allocates space for small summary stats gathered after each segment.
        function this = initialize(this, sim)
            %indices
            this.stateIndex_time = sim.stateIndex('Time');
            this.stateIndex_mass = sim.stateIndex('Mass');
            this.stateIndex_rna = sim.stateIndex('Rna');
            this.stateIndex_metabolite = sim.stateIndex('Metabolite');
            this.stateIndex_monomer = sim.stateIndex('ProteinMonomer');
            this.stateIndex_complex = sim.stateIndex('ProteinComplex');
            this.stateIndex_metabolicReaction = sim.stateIndex('MetabolicReaction');
            this.stateIndex_chromosome = sim.stateIndex('Chromosome');
            this.stateIndex_geometry = sim.stateIndex('Geometry');
            this.stateIndex_ftsZRing = sim.stateIndex('FtsZRing');
            
			this.processIndex_transcription = sim.processIndex('Transcription');
			this.processIndex_translation = sim.processIndex('Translation');
            this.processIndex_replicationInitiation = sim.processIndex('ReplicationInitiation');
            this.processIndex_replication = sim.processIndex('Replication');
            this.processIndex_ftsZPolymerization = sim.processIndex('FtsZPolymerization');
            
            %ids
            this.ntpWholeCellModelIDs = sim.state('Metabolite').wholeCellModelIDs(sim.state('Metabolite').ntpIndexs);
            this.aaWholeCellModelIDs = sim.state('Metabolite').wholeCellModelIDs(sim.state('Metabolite').aminoAcidIndexs);
            
            %initialize timer
            this.timer = tic;
            
            %allocate memory
            nSeg = sim.lengthSec / this.stepSizeSec + 1;
            validateattributes(nSeg, {'numeric'}, {'real', 'positive', 'integer'});
            
            this.log = struct();
            this.log.runTime               = zeros(1, nSeg, 'double');         %run time of simulation
            this.log.growth                = zeros(1, nSeg, 'double');         %growth rate
            this.log.mass                  = zeros(1, nSeg, 'double');         %cell mass
            this.log.time                  = 0:this.stepSizeSec:sim.lengthSec; %time (seconds)
            this.log.metabolites           = zeros(13, nSeg, 'int32');         %total numbers of ATP, ADP, AMP
            this.log.metaboliteUsages      = zeros(2, nSeg, 'int32');         %total numbers of ATP, ADP, AMP
            this.log.metaboliteProductions = zeros(2, nSeg, 'int32');         %total numbers of ATP, ADP, AMP
            this.log.rnas                  = zeros(7, nSeg, 'int32');          %total numbers of RNA; mature mRNA, rRNA, sRNA, tRNA
            this.log.proteins              = zeros(7, nSeg, 'int32');          %total number of protein monomers and complexes
            this.log.rnaPolymerases        = zeros(1, nSeg, 'int16');          %total number of RNA polymerases
            this.log.ribosomes             = zeros(1, nSeg, 'int16');          %total number of ribosomes
            this.log.supercoiled           = zeros(2, nSeg, 'double');         %number of supercoiled regions
            this.log.dnaA                  = zeros(2, nSeg, 'int16');          %numbers of free and total DnaA monomers
            this.log.dnaAFunctionalBoxes   = zeros(5, nSeg, 'int8');           %polymerization status of R1-5 functional DnaA boxes
            this.log.replisome             = zeros(6, nSeg, 'int32');          %positions of helicase and leading and lagging polymerases
            this.log.ftsZ                  = zeros(2, nSeg, 'int16');          %numbers of free and total FtsZ monomers
            this.log.ftsZRing              = zeros(4, nSeg, 'int32');          %size of FtsZ ring
            this.log.pinchedDiameter       = zeros(1, nSeg, 'double');         %pinched diameter
            this.log.ploidy                = zeros(1, nSeg, 'double');         %ploidy (genetic content--should range from 1 to 2)
            
            %initialize print, plot
            this.initializePrint();
            this.initializePlot();
            
            %append log
            this.append(sim);
        end
        
        %Appends a summary of the current state at a specific index.
        function this = append(this, sim)
            t = sim.states{this.stateIndex_time};
            i = t.values / this.stepSizeSec + 1;
            if mod(i, 1) ~= 0
                return;
            end
            
            %append log
            this.appendHelper(i, sim);
            
            %print and plot
            this.print(i);
            this.plot(i);
        end
        
        function appendHelper(this, i, sim)
            %references
            t = sim.states{this.stateIndex_time};
            mass = sim.states{this.stateIndex_mass};
            r = sim.states{this.stateIndex_rna};
            m = sim.states{this.stateIndex_metabolite};
            pm = sim.states{this.stateIndex_monomer};
            pc = sim.states{this.stateIndex_complex};
            mr = sim.states{this.stateIndex_metabolicReaction};
            c = sim.states{this.stateIndex_chromosome};
            geometry = sim.states{this.stateIndex_geometry};
            ftsZRing = sim.states{this.stateIndex_ftsZRing};
            
            %indices
            iCytosol = sim.compartment.cytosolIndexs;
            iMembrane = sim.compartment.membraneIndexs;
            
            %append
            this.log.runTime(1, i) = toc(this.timer); this.timer = tic;
            this.log.time(1, i) = t.values;
            
            this.log.mass(1, i) = sum(mass.cellDry) / mass.cellInitialDryWeight;
            this.log.growth(1, i) = mr.growth;
            
            this.log.metabolites(1, i) = sum(m.counts(m.atpIndexs, iCytosol));
            this.log.metabolites(2, i) = sum(m.counts(m.adpIndexs, iCytosol));
            this.log.metabolites(3, i) = sum(m.counts(m.ampIndexs, iCytosol));
            this.log.metabolites(4, i) = sum(m.counts(m.ntpIndexs, iCytosol));
            [this.log.metabolites(5, i), this.log.metabolites(6, i)] = min(m.counts(m.ntpIndexs, iCytosol));
            this.log.metabolites(7, i) = sum(m.counts(m.aminoAcidIndexs(1:20), iCytosol));
            [this.log.metabolites(8, i), this.log.metabolites(9, i)] = min(m.counts(m.aminoAcidIndexs(1:20), iCytosol));
            this.log.metabolites(10,i) = sum(m.counts(m.lipidIndexs, iMembrane));
            this.log.metabolites(11,i) = sum(m.counts(m.carbohydrateIndexs, iCytosol));
            this.log.metabolites(12,i) = sum(m.counts(m.dntpIndexs, iCytosol));
            this.log.metabolites(13,i) = sum(m.counts(m.antibioticIndexs, iCytosol));
            
            if ~isempty(m.processUsages)
                tmp1 = m.processUsages(sub2ind(size(m.counts), m.ntpIndexs(1), iCytosol), :);
                tmp2 = m.processUsages(sub2ind(size(m.counts), m.ntpIndexs(3), iCytosol), :);
                this.log.metaboliteUsages(1, i) = sum(max(0, -tmp1), 2);
                this.log.metaboliteUsages(2, i) = sum(max(0, -tmp2), 2);
                this.log.metaboliteProductions(1, i) = sum(max(0, tmp1), 2);
                this.log.metaboliteProductions(2, i) = sum(max(0, tmp2), 2);
            end
            
            this.log.rnas(1, i) = sum(r.counts(:));
            this.log.rnas(2, i) = sum(...
                + r.counts(r.matureIndexs(r.matureMRNAIndexs), iCytosol)...
                + r.counts(r.misfoldedIndexs(r.matureMRNAIndexs), iCytosol)...
                + r.counts(r.damagedIndexs(r.matureMRNAIndexs), iCytosol)...
                );
            this.log.rnas(3, i) = sum(...
                + r.counts(r.matureIndexs(r.matureRRNAIndexs), iCytosol)...
                + r.counts(r.misfoldedIndexs(r.matureRRNAIndexs), iCytosol)...
                + r.counts(r.damagedIndexs(r.matureRRNAIndexs), iCytosol)...
                );
            this.log.rnas(4, i) = sum(...
                + r.counts(r.matureIndexs(r.matureSRNAIndexs), iCytosol) ...
                + r.counts(r.boundIndexs(r.matureSRNAIndexs), iCytosol) ...
                + r.counts(r.aminoacylatedIndexs(r.matureSRNAIndexs), iCytosol) ...
                + r.counts(r.misfoldedIndexs(r.matureSRNAIndexs), iCytosol) ...
                + r.counts(r.damagedIndexs(r.matureSRNAIndexs), iCytosol) ...
                );
            this.log.rnas(5, i) = sum(...
                + r.counts(r.matureIndexs(r.matureTRNAIndexs), iCytosol) ...
                + r.counts(r.aminoacylatedIndexs(r.matureTRNAIndexs), iCytosol) ...
                + r.counts(r.misfoldedIndexs(r.matureTRNAIndexs), iCytosol) ...
                + r.counts(r.damagedIndexs(r.matureTRNAIndexs), iCytosol) ...
                );
            this.log.rnas(6, i) = this.log.rnas(1, i) - sum(this.log.rnas(2:5, i));
            this.log.rnas(7, i) = sum(r.counts(r.damagedIndexs, iCytosol));
            
            this.log.proteins(1, i) = sum(pm.counts(:));
            this.log.proteins(2, i) = sum(sum(...
                + pm.counts(pm.matureIndexs, :) ...
                + pm.counts(pm.boundIndexs, :) ...
                + pm.counts(pm.inactivatedIndexs, :) ...
                + pm.counts(pm.damagedIndexs, :) ...
                + pm.counts(pm.misfoldedIndexs, :) ...
                ));
            this.log.proteins(3, i) = this.log.proteins(1, i) - this.log.proteins(2, i);
            this.log.proteins(4, i) = sum(pc.counts(:));
            this.log.proteins(5, i) = sum(sum(...
                + pc.counts(pc.matureIndexs, :) ...
                + pc.counts(pc.boundIndexs, :) ...
                + pc.counts(pc.inactivatedIndexs, :) ...
                + pc.counts(pc.damagedIndexs, :) ...
                + pc.counts(pc.misfoldedIndexs, :) ...
                ));
            this.log.proteins(6, i) = this.log.proteins(4, i) - this.log.proteins(5, i);
            this.log.proteins(7, i) = ...
                + sum(sum(pm.counts(pm.damagedIndexs, :))) ...
                + sum(sum(pc.counts(pc.damagedIndexs, :)));
			
            if this.processIndex_transcription
                transcription = sim.processes{this.processIndex_transcription};               			
                this.log.rnaPolymerases(1, i) = ...
                    + transcription.enzymes(transcription.enzymeIndexs_rnaPolymerase) ...
                    + transcription.enzymes(transcription.enzymeIndexs_rnaPolymeraseHoloenzyme) ...
                    + transcription.boundEnzymes(transcription.enzymeIndexs_rnaPolymerase) ...
                    + transcription.boundEnzymes(transcription.enzymeIndexs_rnaPolymeraseHoloenzyme);
            end
            
            if this.processIndex_translation
                translation = sim.processes{this.processIndex_translation};
                this.log.ribosomes(1, i) = ...
                    + translation.enzymes(translation.enzymeIndexs_ribosome30S) ...
                    + translation.enzymes(translation.enzymeIndexs_ribosome30SIF3) ...
                    + translation.enzymes(translation.enzymeIndexs_ribosome50S) ...
                    + translation.enzymes(translation.enzymeIndexs_ribosome70S) ...
                    + translation.boundEnzymes(translation.enzymeIndexs_ribosome30S) ...
                    + translation.boundEnzymes(translation.enzymeIndexs_ribosome30SIF3) ...
                    + translation.boundEnzymes(translation.enzymeIndexs_ribosome50S) ...
                    + translation.boundEnzymes(translation.enzymeIndexs_ribosome70S);
            end
            
            this.log.supercoiled(1, i) = collapse(c.supercoiled) / 2;
            this.log.supercoiled(2, i) = (collapse(c.linkingNumbers) - collapse(c.polymerizedRegions)/c.relaxedBasesPerTurn) / (collapse(c.polymerizedRegions)/c.relaxedBasesPerTurn);
            
            if this.processIndex_replicationInitiation
                repInit = sim.processes{this.processIndex_replicationInitiation};
                
                this.log.dnaA(:, i) = [
                    repInit.enzymes(repInit.enzymeIndexs_DnaA_1mer_ATP)
                    [1 1:7 1:7] * (...
                    + repInit.enzymes([repInit.enzymeIndexs_DnaA; repInit.enzymeIndexs_DnaA_Nmer_ATP; repInit.enzymeIndexs_DnaA_Nmer_ADP]) ...
                    + repInit.boundEnzymes([repInit.enzymeIndexs_DnaA; repInit.enzymeIndexs_DnaA_Nmer_ATP; repInit.enzymeIndexs_DnaA_Nmer_ADP]) ...
                    )];
                dnaAPol = [
                    repInit.calculateDnaAR1234Polymerization()
                    repInit.calcuateIsDnaAR5Occupied()
                    ];
                this.log.dnaAFunctionalBoxes(:, i) = dnaAPol(:, 1);
            end
            if this.processIndex_replication
                rep = sim.processes{this.processIndex_replication};
                
                this.log.replisome(:, i) = [
                    rep.helicasePosition rep.leadingPosition rep.laggingPosition
                    ];
            end
            
            this.log.ploidy(1,i) = c.ploidy;
            
            if this.processIndex_ftsZPolymerization
                ftsZPol = sim.processes{this.processIndex_ftsZPolymerization};
                
                freeFtsZMonomers = ftsZPol.enzymes(ftsZPol.enzymeIndexs_FtsZ);
                totFtsZMonomers = ...
                    + ftsZPol.enzymes(ftsZPol.enzymeIndexs_FtsZ) ...
                    + ftsZPol.enzymes(ftsZPol.enzymeIndexs_FtsZ_GDP) ...
                    + ftsZPol.enzymes(ftsZPol.enzymeIndexs_FtsZ_activated)' * (1:9)' ...
                    + ftsZPol.boundEnzymes(ftsZPol.enzymeIndexs_FtsZ) ...
                    + ftsZPol.boundEnzymes(ftsZPol.enzymeIndexs_FtsZ_GDP) ...
                    + ftsZPol.boundEnzymes(ftsZPol.enzymeIndexs_FtsZ_activated)' * (1:9)';
                this.log.ftsZ(:, i) = [
                    freeFtsZMonomers
                    totFtsZMonomers
                    ];
            end
            this.log.ftsZRing(:, i) = [
                ftsZRing.numEdgesOneStraight
                ftsZRing.numEdgesTwoStraight
                ftsZRing.numEdgesTwoBent
                ftsZRing.numResidualBent
                ];
            this.log.pinchedDiameter(1, i) = geometry.pinchedDiameter;
        end
        
        function this = finalize(this, sim)
            t = sim.states{this.stateIndex_time};
            tmpi = t.values / this.stepSizeSec + 1;
            i = ceil(tmpi);
            
            %trim log
            fields = fieldnames(this.log);
            for j = 1:numel(fields)
                this.log.(fields{j}) = this.log.(fields{j})(:, 1:i);
            end
            
            %append
            this.appendHelper(i, sim);
            
            %print and plot
            if rem(tmpi, 1) ~= 0
                this.print(i);
            end
            this.finalizePrint();
            this.plot(i);
            
            %save
            this.save();
        end
        
        function this = save(this)
            if isempty(this.outputDirectory)
                return;
            end
            
            log = this.log; %#ok<PROP,NASGU>
            save([this.outputDirectory filesep 'summary.mat'], '-struct', 'log');
            
            if this.verbosity > 1
                try
                    saveas(this.figureHandle, [this.outputDirectory filesep 'summary.pdf']);
                catch exception
                    warning('WholeCell:warning', 'Summary figure could not be saved: %s.', exception.message);
                end
            end
        end
    end
    
    methods (Access = protected)
        function initializePrint(this)
            if this.verbosity <= 0
                return;
            end
            
            fprintf('%-5s %-3s %-4s %-10s %-89s %-29s %-11s %-11s %-3s %-7s %-9s %-41s %-6s %-7s %-11s %-10s\n', ...
                '', ...
                '', ...
                '', ...
                '', ...
                'Metabolites',...
                'RNAs', ...
                'Monomers', 'Complexes', ...
                '', '', '', ...
                'Replisome', ...
                '', ...
                '', '', '' ...
                );
            fprintf('%-5s %-3s %-4s %-10s %-89s %-29s %-11s %-11s %-3s %-7s %-9s %-41s %-6s %-7s %-11s %-10s\n', ...
                '', ...
                '', ...
                '', ...
                '', ...
                repmat('=', 1, 89),...
                repmat('=', 1, 29), ...
                repmat('=', 1, 11), repmat('=', 1, 11), ...
                '', '', '', ...
                repmat('=', 1, 41), ...
                '', ...
                '', '', '' ...
                );
            fprintf('%-5s %-3s %-4s %-10s %-8s %-8s %-8s %-8s %-9s %-8s %-9s %-8s %-15s %-4s %-4s %-4s %-4s %-4s %-4s %-5s %-5s %-5s %-5s %-3s %-7s %-9s %-13s %-13s %-13s %-6s %-7s %-11s %-10s\n',...
                'Time', ...
                'RT', ...
                'Mass', ...
                'Growth', ...
                'ATP', 'ADP', 'AMP', 'NTPs', 'Min NTP', 'AAs', 'Min AA', 'Lipids', 'Polysaccharides', ...
                'Tot', 'Nasc', 'mRNA', 'rRNA', 'sRNA', 'tRNA', ...
                'Nasct', 'Matur', 'Nasct', 'Matur', ...
                'Sup', 'DnaA', 'R1-5', ...
                'Helicase Pos', 'Leading Pos', 'Lagging Pos', ...
                'Ploidy', ...
                'FtsZ', 'FtsZ Ring', 'Pinch Diam' ...
                );
            fprintf('%-5s %-3s %-4s %-10s %-8s %-8s %-8s %-8s %-9s %-8s %-9s %-8s %-15s %-4s %-4s %-4s %-4s %-4s %-4s %-5s %-5s %-5s %-5s %-3s %-7s %-9s %-13s %-13s %-13s %-6s %-7s %-11s %-10s\n',...
                repmat('=', 1, 5), ...
                repmat('=', 1, 3), ...
                repmat('=', 1, 4), ...
                repmat('=', 1, 10), ...
                repmat('=', 1, 8), repmat('=', 1, 8), repmat('=', 1, 8), repmat('=', 1, 8), repmat('=', 1, 9), repmat('=', 1, 8), repmat('=', 1, 9), repmat('=', 1, 8), repmat('=', 1, 15), ...
                repmat('=', 1, 4), repmat('=', 1, 4), repmat('=', 1, 4), repmat('=', 1, 4), repmat('=', 1, 4), repmat('=', 1, 4), ...
                repmat('=', 1, 5), repmat('=', 1, 5), repmat('=', 1, 5), repmat('=', 1, 5), ...
                repmat('=', 1, 3), repmat('=', 1, 7), repmat('=', 1, 9), ...
                repmat('=', 1, 13), repmat('=', 1, 13), repmat('=', 1, 13), ...
                repmat('=', 1,6), ...
                repmat('=', 1, 7), repmat('=', 1, 11), repmat('=', 1, 10) ...
                );
        end
        
        function print(this, i)
            if this.verbosity <= 0
                return;
            end
            
            if ispc
                format = '%5d %3.0f %4.2f %5.3e %8d %8d %8d %8d %5d %3s %8d %5d %3s %8d %15d %4d %4d %4d %4d %4d %4d %5d %5d %5d %5d %3d %3d %3d %d %d %d %d %d %6d %6d %6d %6d %6d %6d %6.3f %3d %3d %2d %2d %2d %2d %5.3e\n';
            else
                format = '%5d %3.0f %4.2f %6.3e %8d %8d %8d %8d %5d %3s %8d %5d %3s %8d %15d %4d %4d %4d %4d %4d %4d %5d %5d %5d %5d %3d %3d %3d %d %d %d %d %d %6d %6d %6d %6d %6d %6d %6.3f %3d %3d %2d %2d %2d %2d %6.3e\n';
            end
            fprintf(format,...
                this.log.time(1, i),...
                this.log.runTime(1, i), ...
                this.log.mass(1, i), ...
                this.log.growth(1, i), ...
                this.log.metabolites(1, i),...
                this.log.metabolites(2, i),...
                this.log.metabolites(3, i),...
                this.log.metabolites(4, i),...
                this.log.metabolites(5, i),...
                this.ntpWholeCellModelIDs{this.log.metabolites(6, i)},...
                this.log.metabolites(7, i),...
                this.log.metabolites(8, i),...
                this.aaWholeCellModelIDs{this.log.metabolites(9, i)},...
                this.log.metabolites(10,i),...
                this.log.metabolites(11,i),...
                this.log.rnas(1, i),...
                this.log.rnas(6, i),...
                this.log.rnas(2, i),...
                this.log.rnas(3, i),...
                this.log.rnas(4, i),...
                this.log.rnas(5, i),...
                this.log.proteins(3, i),...
                this.log.proteins(2, i),...
                this.log.proteins(6, i),...
                this.log.proteins(5, i),...
                this.log.supercoiled(1, i), ...
                this.log.dnaA(1, i), ...
                this.log.dnaA(2, i), ...
                this.log.dnaAFunctionalBoxes(1, i), ...
                this.log.dnaAFunctionalBoxes(2, i), ...
                this.log.dnaAFunctionalBoxes(3, i), ...
                this.log.dnaAFunctionalBoxes(4, i), ...
                this.log.dnaAFunctionalBoxes(5, i), ...
                this.log.replisome(1, i), ...
                this.log.replisome(2, i), ...
                this.log.replisome(3, i), ...
                this.log.replisome(4, i), ...
                this.log.replisome(5, i), ...
                this.log.replisome(6, i), ...
                this.log.ploidy(1,i), ...
                this.log.ftsZ(1, i), ...
                this.log.ftsZ(2, i), ...
                this.log.ftsZRing(1, i), ...
                this.log.ftsZRing(2, i), ...
                this.log.ftsZRing(3, i), ...
                this.log.ftsZRing(4, i), ...
                this.log.pinchedDiameter(1, i) ...
                );
        end
        
        function finalizePrint(this)
            if this.verbosity <= 0
                return;
            end
            
            fprintf('Total simulation time: %d s\n', this.log.time(1, end));
            fprintf('Total runtime: %.3f s\n', sum(this.log.runTime));
            fprintf('\n');
        end
        
        function initializePlot(this)
            if this.verbosity <= 1
                return;
            end
            
            %create new figure
            this.figureHandle = figure();
        end
        
        %Plot progress of evolution of state of organism.
        function plot(this, ~)
            if this.verbosity <= 1
                return;
            end
            
            %run time
            axesHandle = subplot(2, 4, 1, 'Parent', this.figureHandle);
            plot(axesHandle, 0:numel(this.log.time)-1, this.log.runTime/60);
            title(axesHandle, 'Runtime')
            xlabel(axesHandle, 'Segment');
            ylabel(axesHandle, 'Run Time (m)');
            xlim(axesHandle, [0 numel(this.log.time)-1]);
            
            %metabolites
            axesHandle = subplot(2, 4, 2, 'Parent', this.figureHandle);
            h = plot(axesHandle, this.log.time / 3600, this.log.metabolites([1 2 3 4 5 7 8], :)');
            title(axesHandle, 'Metabolites')
            xlabel(axesHandle, 'Time (h)');
            ylabel(axesHandle, 'Metabolites');
            xlim(axesHandle, [0 this.log.time(end)] / 3600);
            legend(h, {'ATP', 'ADP', 'AMP', 'Tot NTPs', 'Min NTP', 'Tot AAs', 'Min AA'});
            
            %RNAs
            axesHandle = subplot(2, 4, 3, 'Parent', this.figureHandle);
            h = plot(axesHandle, this.log.time / 3600, this.log.rnas([1 6 2 3 4 5], :)');
            title(axesHandle, 'RNA')
            xlabel(axesHandle, 'Time (h)');
            ylabel(axesHandle, 'RNA');
            xlim(axesHandle, [0 this.log.time(end)] / 3600);
            legend(h, {'Total', 'Nascent', 'mRNA', 'rRNA', 'sRNA', 'tRNA'});
            
            %proteins
            axesHandle = subplot(2, 4, 4, 'Parent', this.figureHandle);
            h = plot(axesHandle, this.log.time/3600, this.log.proteins([1 3 2 4 6 5], :)');
            title(axesHandle, 'Protein')
            xlabel(axesHandle, 'Time (h)');
            ylabel(axesHandle, 'Protein');
            xlim(axesHandle, [0 this.log.time(end)]/3600);
            legend(h, {'Total Monomers', 'Nascent Monomers', 'Mature Monomers', 'Total Complexes', 'Nascent Complexes', 'Mature Complexes'});
            
            %Supercoiled
            axesHandle = subplot(2, 4, 5, 'Parent', this.figureHandle);
            plot(axesHandle, this.log.time/3600, this.log.supercoiled(1, :)');
            title(axesHandle, 'Supercoiled')
            xlabel(axesHandle, 'Time (h)');
            ylabel(axesHandle, 'No. Supercoiled Regions');
            xlim(axesHandle, [0 this.log.time(end)]/3600);
            
            %DnaA
            axesHandle = subplot(2, 4, 6, 'Parent', this.figureHandle);
            h = plot(axesHandle, this.log.time/3600, [this.log.dnaA/max(1, this.log.dnaA(1,1)); int16(this.log.dnaAFunctionalBoxes)]');
            title(axesHandle, 'DnaA')
            xlabel(axesHandle, 'Time (h)');
            ylabel(axesHandle, 'DnaA Count');
            xlim(axesHandle, [0 this.log.time(end)]/3600);
            legend(h, {'Free DnaA-ATP', 'Total DnaA', 'R1', 'R2', 'R3', 'R4', 'R5'});
            
            %Replisome
            axesHandle = subplot(2, 4, 7, 'Parent', this.figureHandle);
            h = plot(axesHandle, this.log.time/3600, this.log.replisome');
            title(axesHandle, 'Replisome')
            xlabel(axesHandle, 'Time (h)');
            ylabel(axesHandle, 'Position');
            xlim(axesHandle, [0 this.log.time(end)]/3600);
            legend(h, {'Helicase-1', 'Helicase-2', 'Leading-1', 'Leading-2', 'Lagging-1', 'Lagging-2'});
            
            %FtsZ
            axesHandle = subplot(2, 4, 8, 'Parent', this.figureHandle);
            h = plot(axesHandle, this.log.time/3600, [double(this.log.ftsZ/this.log.ftsZ(1, 1)); double(this.log.ftsZRing); this.log.pinchedDiameter/this.log.pinchedDiameter(1)]');
            title(axesHandle, 'FtsZ')
            xlabel(axesHandle, 'Time (h)');
            ylabel(axesHandle, 'FtsZ');
            xlim(axesHandle, [0 this.log.time(end)]/3600);
            legend(h, {'Free FtsZ-GTP', 'Total FtsZ', 'Ring-1strt', 'Ring-1bent', 'Ring-2strt', 'Ring-2 bent', 'Pinched Diameter'});
            
            %force display of figures
            drawnow;
        end
    end
    
    % Analyzes summary output of a batch of simulations. Stores analysis in mat file.
    methods (Static = true)
        function [simStats] = summarizeSimulations(baseDir)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            simStats = SummaryLogger.getSimulationStatistics(baseDir);
            
            nSimulations = size(simStats, 1);
            
            summaryFileName = [baseDir filesep 'summary.html'];
            stateFileName = [baseDir filesep 'initialAndFinalStates.html'];
            errorFileName = [baseDir filesep 'errors.html'];
            
            %% calc statistics
            stateFid = fopen(stateFileName, 'w');
            errorFid = fopen(errorFileName, 'w');
            fprintf(stateFid, '<html>\n<head>\n');
            fprintf(stateFid, '<style type="text/css">\n');
            fprintf(stateFid, '  table.states {width:12000px;}\n');
            fprintf(stateFid, '  table.states thead tr th, table.states tbody tr td{width:auto;padding-left:10px;}\n');
            fprintf(stateFid, '  table.states thead tr th:first-child, table.states tbody tr td:first-child{width:auto;padding-left:0px;}\n');
            fprintf(stateFid, '  table.states thead tr:first-child th {border-top:1px solid black;}\n');
            fprintf(stateFid, '  table.states thead tr:last-child th {border-bottom:1px solid black;}\n');
            fprintf(stateFid, '  table.states tbody tr:last-child td {border-bottom:1px solid black;}\n');
            fprintf(stateFid, '  table.states tr th:first-child, table.states tr td:first-child {border-right:1px solid black}\n');
            fprintf(stateFid, '  table.states thead tr:last-child th:nth-child(%d), table.states tbody tr td:nth-child(%d) {border-right:1px solid black}\n', 1+6, 1+6);
            fprintf(stateFid, '  table.states thead tr:last-child th:nth-child(%d), table.states tbody tr td:nth-child(%d) {border-right:1px solid black}\n', 1+6+48, 1+6+48);
            fprintf(stateFid, '  table.states thead tr th:first-child, table.states tbody tr td:first-child {border-left:1px solid black}\n');
            fprintf(stateFid, '  table.states thead tr th:last-child, table.states tbody tr td:last-child {border-right:1px solid black}\n');
            fprintf(stateFid, '  table.states thead tr th {background-color: #CCCCCC;}\n');
            fprintf(stateFid, '  table.states tbody tr:nth-child(even) td {background-color: #DDDDDD;}\n');
            fprintf(stateFid, '  table.states tbody tr:nth-child(odd) td {background-color: #EEEEEE;}\n');
            fprintf(stateFid, '</style>\n');
            fprintf(stateFid, '</head>\n');
            fprintf(stateFid, '<body>\n<table class="states" cellpadding="0" cellspacing="0">\n<thead>\n%s\n</thead>', SummaryLogger.formatColNames());
            fprintf(errorFid, '<html>\n<body>\n');
            for i = 1:nSimulations
                %calc simulation status
                if ~exist(sprintf('%s/%d/state-0.mat', baseDir, i), 'file')
                elseif ~exist(sprintf('%s/%d/err.mat', baseDir, i), 'file') && ~exist(sprintf('%s/%d/summary.mat', baseDir, i), 'file')
                elseif exist(sprintf('%s/%d/err.mat', baseDir, i), 'file')
                    tmpFid = fopen(sprintf('%s/%d/err.log', baseDir, i), 'r');
                    error = reshape(fread(tmpFid, '*char'), 1, []);
                    error = strrep(error, sprintf('\n\n\n\n'), sprintf('\n'));
                    error = strrep(error, sprintf('\n\n'), sprintf('\n'));
                    fclose(tmpFid);
                    fprintf(errorFid, '<p>\n<b>Simulation #%d</b><br/>\n<pre>\n%s\n</pre>\n</p>\n', i, error);
                else
                    summary = load(sprintf('%s/%d/summary.mat', baseDir, i));
                end
                
                if simStats(i, SummaryLogger.SIM_STATUS_INDEX_STATUS) >= 0
                    fprintf(stateFid, '<tr>\n<td>%d</td>\n%s\n%s\n%s\n</tr>\n', i, SummaryLogger.formatStats(simStats(i, :)), SummaryLogger.formatState(summary, 1), SummaryLogger.formatState(summary, numel(summary.time)));
                else
                    fprintf(stateFid, '<tr>\n<td>%d</td>\n%s\n%s\n%s\n</tr>\n', i, SummaryLogger.formatStats(simStats(i, :)), SummaryLogger.formatState(), SummaryLogger.formatState());
                end
            end
            fprintf(stateFid, '</tbody>\n</table>\n</body>\n</html>\n');
            fprintf(errorFid, '</body>\n</html>\n');
            fclose(stateFid);
            fclose(errorFid);
            
            if ~any(simStats(:, SummaryLogger.SIM_STATUS_INDEX_STATUS) == -1)
                errorFid = fopen(errorFileName, 'w');
                fprintf(errorFid, '<html>\n<body>\n');
                fprintf(errorFid, 'No run-time errors were encountered.\n');
                fprintf(errorFid, '</body>\n</html>\n');
                fclose(errorFid);
            end
            
            fid = fopen(summaryFileName, 'w');
            fwrite(fid, SummaryLogger.getSummaryReport(simStats));
            fclose(fid);
            
            %% plot
            figHandle = figure();
            
            %initial vs. cumulative growth
            plot(simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_GROWTH)'*3600, ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_CUM_GROWTH)', '.');
            xlabel('Initial growth rate (cell/h)');
            ylabel('Cumulative growth rate (cell)');
            saveas(gcf, [baseDir filesep 'initialVsCumulativeGrowth.pdf']);
            
            %replication initiation time vs replication, cytokinesis times
            h = plot(simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME)' / 3600, [...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME)]'/ 3600, '.');
            xlabel('Replication Initiation Time (h)');
            ylabel('Time (h)');
            legend(h, {'Replication', 'Cytokinesis'});
            saveas(gcf, [baseDir filesep 'cellCyclePhaseLengths.pdf']);
            
            %replication initiation duration vs replication, cytokinesis
            %durations
            h = plot(simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME)' / 3600, [...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME)-simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME)-simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME)]'/ 3600, '.');
            xlabel('Replication Initiation Time (h)');
            ylabel('Duration (h)');
            legend(h, {'Replication', 'Cytokinesis'});
            saveas(gcf, [baseDir filesep 'cellCyclePhaseDurations.pdf']);
            
            %initial growth rate vs replication initiation, replication,
            %cytokinesis times
            h = plot(simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_GROWTH)'*3600, [...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME)]'/ 3600, '.');
            xlabel('Initial growth rate (cell/h)');
            ylabel('Time (h)');
            legend(h, {'Replication initiation', 'Replication', 'Cytokinesis'});
            saveas(gcf, [baseDir filesep 'initialGrowthVsCellCyclePhaseTimes.pdf']);
            
            %cumulative growth vs replication initiation, replication,
            %cytokinesis times
            h = plot(simStats(:, SummaryLogger.SIM_STATUS_INDEX_CUM_GROWTH), [...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME) simStats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME)]'/ 3600, '.');
            xlabel('Cumulative growth rate (cell)');
            ylabel('Time (h)');
            legend(h, {'Replication initiation', 'Replication', 'Cytokinesis'});
            saveas(gcf, [baseDir filesep 'cumulativeGrowthVsCellCyclePhaseTimes.pdf']);
            
            %initial growth rate vs replication initiation, replication,
            %cytokinesis durations
            h = plot(simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_GROWTH)'*3600, [...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) ....
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME)-simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME)-simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME)]'/ 3600, '.');
            xlabel('Initial growth rate (cell/h)');
            ylabel('Duration (h)');
            legend(h, {'Replication initiation', 'Replication', 'Cytokinesis'});
            saveas(gcf, [baseDir filesep 'initialGrowthVsCellCyclePhaseDurations.pdf']);
            
            %cumulative growth vs replication initiation, replication,
            %cytokinesis durations
            h = plot(simStats(:, SummaryLogger.SIM_STATUS_INDEX_CUM_GROWTH), [...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME)-simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) ...
                simStats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME)-simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME)]'/ 3600, '.');
            xlabel('Cumulative growth rate (cell)');
            ylabel('Duration (h)');
            legend(h, {'Replication initiation', 'Replication', 'Cytokinesis'});
            saveas(gcf, [baseDir filesep 'cumulativeGrowthVsCellCyclePhaseDurations.pdf']);
            
            close(figHandle);
        end
        
        function simStats = getSimulationStatistics(baseDir, selectedSimulations)
            import edu.stanford.covert.cell.sim.util.SimulationDiskUtil;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            if nargin < 2
                selectedSimulations = 1:SimulationDiskUtil.getNumSimulations(baseDir);
            end
            
            nSimulations = numel(selectedSimulations);
            
            %% calc statistics
            simStats = NaN(nSimulations, 11);
            for j = 1:nSimulations
                i = selectedSimulations(j);
                
                %calc simulation status
                if ~exist(sprintf('%s/%d/state-0.mat', baseDir, i), 'file')
                    simStats(j, SummaryLogger.SIM_STATUS_INDEX_STATUS) = SummaryLogger.SIM_STATUS_DIDNT_START;
                elseif ~exist(sprintf('%s/%d/err.mat', baseDir, i), 'file') && ~exist(sprintf('%s/%d/summary.mat', baseDir, i), 'file')
                    simStats(j, SummaryLogger.SIM_STATUS_INDEX_STATUS) = SummaryLogger.SIM_STATUS_STILL_RUNNING;
                elseif exist(sprintf('%s/%d/err.mat', baseDir, i), 'file')
                    simStats(j, SummaryLogger.SIM_STATUS_INDEX_STATUS) = SummaryLogger.SIM_STATUS_TERMINATED_WITH_ERROR;
                else
                    summary = load(sprintf('%s/%d/summary.mat', baseDir, i), 'time', 'mass', 'growth', 'dnaAFunctionalBoxes', 'replisome', 'pinchedDiameter');
                    if summary.pinchedDiameter(end) == 0
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_STATUS) = SummaryLogger.SIM_STATUS_COMPLETED_WITH_DIVISION;
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME) = summary.time(end);
                    else
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_STATUS) = SummaryLogger.SIM_STATUS_COMPLETED_WITHOUT_DIVISION;
                    end
                    
                    %calc simulation statistics
                    idx = find(all(summary.dnaAFunctionalBoxes == repmat([7; 7; 7; 7; 1], 1, numel(summary.time)), 1), 1, 'first');
                    if ~isempty(idx)
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) = summary.time(idx);
                    end
                    
                    idx = find(any(summary.replisome, 1), 1, 'first');
                    if ~isempty(idx)
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_INIT_REP_TIME) = summary.time(idx);
                    end
                    simStats(j, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) = min(...
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME), ...
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_INIT_REP_TIME));
                    
                    idx = find(any(summary.replisome, 1), 1, 'last')+1;
                    if ~isempty(idx) && idx <= size(summary.time, 2)
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_REP_TIME) = summary.time(idx);
                    end
                    
                    simStats(j, SummaryLogger.SIM_STATUS_INDEX_INIT_GROWTH) = summary.growth(1);
                    cellCycleLength = 9.0 * 3600;
                    idx = find(summary.time > cellCycleLength, 1, 'first');
                    if ~isempty(idx)
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_FIN_GROWTH) = summary.growth(idx);
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_CUM_GROWTH) = diff(summary.time(1:2)) * sum(summary.growth(1:idx));
                    else
                        f = (2 - 1) / (exp(log(2) * summary.time(end) / cellCycleLength) - 1); %#ok<CPROP,PROP>
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_FIN_GROWTH) = f * summary.growth(end);
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_CUM_GROWTH) = f * diff(summary.time(1:2)) * sum(summary.growth);
                    end
                    
                    idx = find(summary.mass > 2 * summary.mass(1), 1, 'first');
                    if ~isempty(idx)
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_MASS_DOUBLING_TIME) = summary.time(idx-1);
                    else
                        simStats(j, SummaryLogger.SIM_STATUS_INDEX_MASS_DOUBLING_TIME) = ...
                            cellCycleLength / log(2) * log(2 / summary.mass(end)) + summary.time(end); %#ok<CPROP,PROP>
                    end
                    
                    simStats(j, SummaryLogger.SIM_STATUS_INDEX_FIN_MASS) = summary.mass(end);
                    simStats(j, SummaryLogger.SIM_STATUS_INDEX_INIT_MASS) = summary.mass(1);
                end
            end
        end
        
        function value = formatColNames()
            stateColNames = {
                'Time'
                'Growth'
                'Mass'
                'ATP'
                'ADP'
                'AMP'
                'Total NTPs'
                'Min NTP'
                'Min NTP'
                'Total AAs'
                'Min AA'
                'Min AA'
                'Total Lipids'
                'Total Carbohydrates'
                'Total RNA'
                'Total mature mRNA'
                'Total mature rRNA'
                'Total mature sRNA'
                'Total mature tRNA'
                'Total immature RNA'
                'Total Monomers'
                'Total Mature Monomers'
                'Total Immature Monomers'
                'Total Complexs'
                'Total Mature Complexs'
                'Total Immature Complexs'
                'Supercoiled'
                'Free'
                'Total'
                'R1'
                'R2'
                'R3'
                'R4'
                'R5'
                'Helicase Position - 1'
                'Helicase Position - 2'
                'Leading Polymerase Position - 1'
                'Leading Polymerase Position - 2'
                'Lagging Polymerase Position - 2'
                'Lagging Polymerase Position - 2'
                'Ploidy'
                'Free'
                'Total'
                'One straight'
                'One bent'
                'Two straight'
                'Residual'
                'Pinched Diameter'
                };
            colNames = [{
                'Simulation'
                'Status'
                'Replication initiation time'
                'Initiation of replication time'
                'Replication time'
                'Cytokinesis time'
                'Cumulative growth'
                };
                stateColNames; stateColNames];
            
            value = [];
            
            value = [value sprintf('<tr>\n')];
            value = [value sprintf('<th>&nbsp;</th>\n')];
            value = [value sprintf('<th colspan="6">%s</th>\n', 'Summary')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', numel(stateColNames), 'Initial State')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', numel(stateColNames), 'Final State')];
            value = [value sprintf('</tr>\n')];
            
            value = [value sprintf('<tr>\n')];
            value = [value sprintf('<th>&nbsp;</th>\n')];
            value = [value sprintf('<th colspan="6">%s</th>\n', '&nbsp;')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 1, '&nbsp;')]; %Time
            value = [value sprintf('<th colspan="%d">%s</th>\n', 1, '&nbsp;')]; %Growth
            value = [value sprintf('<th colspan="%d">%s</th>\n', 1, '&nbsp;')]; %Mass
            value = [value sprintf('<th colspan="%d">%s</th>\n', 11, 'Metabolite')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 6, 'RNA')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 6, 'Protein')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 1, '&nbsp;')]; %Supercoiled
            value = [value sprintf('<th colspan="%d">%s</th>\n', 2, 'DnaA')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 5, 'Functional DnaA Boxes')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 6, 'Replisome')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 1, '&nbsp;')]; %ploidy
            value = [value sprintf('<th colspan="%d">%s</th>\n', 2, 'FtsZ')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 4, 'FtsZ Ring')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 1, '&nbsp;')]; %Pinched Diameter
            value = [value sprintf('<th colspan="%d">%s</th>\n', 1, '&nbsp;')]; %Time
            value = [value sprintf('<th colspan="%d">%s</th>\n', 1, '&nbsp;')]; %Growth
            value = [value sprintf('<th colspan="%d">%s</th>\n', 1, '&nbsp;')]; %Mass
            value = [value sprintf('<th colspan="%d">%s</th>\n', 11, 'Metabolite')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 6, 'RNA')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 6, 'Protein')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 1, '&nbsp;')]; %Supercoiled
            value = [value sprintf('<th colspan="%d">%s</th>\n', 2, 'DnaA')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 5, 'Functional DnaA Boxes')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 6, 'Replisome')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 1, '&nbsp;')]; %ploidy
            value = [value sprintf('<th colspan="%d">%s</th>\n', 2, 'FtsZ')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 4, 'FtsZ Ring')];
            value = [value sprintf('<th colspan="%d">%s</th>\n', 1, '&nbsp;')]; %Pinched Diameter
            value = [value sprintf('</tr>\n')];
            
            value = [value sprintf('<tr>\n')];
            for i = 1:numel(colNames)
                value = [value sprintf('<th>%s</th>\n', colNames{i})]; %#ok<AGROW>
            end
            value = [value sprintf('</tr>\n')];
        end
        
        function value = formatStats(stats)
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            switch stats(SummaryLogger.SIM_STATUS_INDEX_STATUS)
                case -3, status = 'Not yet started';
                case -2, status = 'Not yet terminated';
                case -1, status = 'Terminated with error';
                case 0, status = 'Terminated with division';
                case 1, status = 'Terminated without division';
            end
            value = sprintf('<td>%s</td>\n', status);
            
            if ~isnan(stats(SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME))
                value = [value sprintf('<td>%d</td>\n', stats(SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME))];
            else
                value = [value sprintf('<td>&nbsp;</td>\n')];
            end
            
            if ~isnan(stats(SummaryLogger.SIM_STATUS_INDEX_INIT_REP_TIME))
                value = [value sprintf('<td>%d</td>\n', stats(SummaryLogger.SIM_STATUS_INDEX_INIT_REP_TIME))];
            else
                value = [value sprintf('<td>&nbsp;</td>\n')];
            end
            
            if ~isnan(stats(SummaryLogger.SIM_STATUS_INDEX_REP_TIME))
                value = [value sprintf('<td>%d</td>\n', stats(SummaryLogger.SIM_STATUS_INDEX_REP_TIME))];
            else
                value = [value sprintf('<td>&nbsp;</td>\n')];
            end
            
            if ~isnan(stats(SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME))
                value = [value sprintf('<td>%d</td>\n', stats(SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME))];
            else
                value = [value sprintf('<td>&nbsp;</td>\n')];
            end
            
            if ~isnan(stats(SummaryLogger.SIM_STATUS_INDEX_CUM_GROWTH))
                value = [value sprintf('<td>%e</td>\n', stats(SummaryLogger.SIM_STATUS_INDEX_CUM_GROWTH))];
            else
                value = [value sprintf('<td>&nbsp;</td>\n')];
            end
        end
                
        function value = formatState(summary, idx)
            fields = {
                'time' 1
                'growth' 1
                'mass' 1
                'metabolites' 11
                'rnas' 6
                'proteins' 6
                'supercoiled' 1
                'dnaA' 2
                'dnaAFunctionalBoxes' 5
                'replisome' 6
                'ploidy' 1
                'ftsZ' 2
                'ftsZRing' 4
                'pinchedDiameter' 1
                };
            value = [];
            for i = 1:size(fields, 1)                
                for j = 1:fields{i, 2}
                    if nargin == 0 || ~isfield(summary, fields{i, 1}) || j > size(summary.(fields{i, 1}), 1)
                        value = [value sprintf('<td>&nbsp;</td>\n')]; %#ok<AGROW>
                    elseif summary.(fields{i, 1})(j, idx) == ceil(summary.(fields{i, 1})(j, idx))
                        value = [value sprintf('<td>%d</td>\n', summary.(fields{i, 1})(j, idx))]; %#ok<AGROW>
                    else
                        value = [value sprintf('<td>%e</td>\n', summary.(fields{i, 1})(j, idx))]; %#ok<AGROW>
                    end
                end
            end
        end
        
        %% build summary report
        function summary = getSummaryReport(simStats)
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            
            nSimulations = size(simStats, 1);
            nTerminatedWithDivision = sum(simStats(:, SummaryLogger.SIM_STATUS_INDEX_STATUS) == 0);
            nTerminatedWithoutDivision = sum(simStats(:, SummaryLogger.SIM_STATUS_INDEX_STATUS) == 1);
            nTerminatedWithError = sum(simStats(:, SummaryLogger.SIM_STATUS_INDEX_STATUS) == -1);
            nNotYetTerminated = sum(simStats(:, SummaryLogger.SIM_STATUS_INDEX_STATUS) == -2);
            nNotYetStarted = sum(simStats(:, SummaryLogger.SIM_STATUS_INDEX_STATUS) == -3);
            repInitTimes = simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_INIT_TIME) / 3600;
            initRepTimes = simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_REP_TIME) / 3600;
            repTimes = simStats(:, SummaryLogger.SIM_STATUS_INDEX_REP_TIME) / 3600;
            cytTimes = simStats(:, SummaryLogger.SIM_STATUS_INDEX_CYTOKINESIS_TIME) / 3600;
            initGrowths = simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_GROWTH) * 3600;
            finGrowths = simStats(:, SummaryLogger.SIM_STATUS_INDEX_FIN_GROWTH) * 3600;
            cumGrowths = simStats(:, SummaryLogger.SIM_STATUS_INDEX_CUM_GROWTH);
            initMass = simStats(:, SummaryLogger.SIM_STATUS_INDEX_INIT_MASS);
            finMass = simStats(:, SummaryLogger.SIM_STATUS_INDEX_FIN_MASS);
            
            %styles
            summary = [];
            summary = [summary sprintf('<html>\n')];
            summary = [summary sprintf('<head>\n')];
            summary = [summary sprintf('<style type="text/css">\n')];
            summary = [summary sprintf('  table.summary tr th {text-align:left; width:400px;}\n')];
            summary = [summary sprintf('  table.summary tr td {text-align:right; padding-left:10px}\n')];
            summary = [summary sprintf('</style>\n')];
            summary = [summary sprintf('</head>\n')];
            summary = [summary sprintf('<body>\n')];
            
            %status
            summary = [summary sprintf('<table class="summary" cellpadding="0" cellspacing="0">\n')];
            summary = [summary sprintf('<tr><th>%s</th><td>%d</td><td>%5.1f%%</td></tr>\n', 'No. simulations', nSimulations, 100)];
            summary = [summary sprintf('<tr><th>%s</th><td>%d</td><td>%5.1f%%</td></tr>\n', 'Terminated with division', nTerminatedWithDivision, nTerminatedWithDivision/nSimulations*100)];
            summary = [summary sprintf('<tr><th>%s</th><td>%d</td><td>%5.1f%%</td></tr>\n', 'Terminated without division',nTerminatedWithoutDivision, nTerminatedWithoutDivision/nSimulations*100)];
            summary = [summary sprintf('<tr><th>%s</th><td>%d</td><td>%5.1f%%</td></tr>\n', 'Terminated with error', nTerminatedWithError, nTerminatedWithError/nSimulations*100)];
            summary = [summary sprintf('<tr><th>%s</th><td>%d</td><td>%5.1f%%</td></tr>\n', 'Not yet terminated', nNotYetTerminated, nNotYetTerminated/nSimulations*100)];
            summary = [summary sprintf('<tr><th>%s</th><td>%d</td><td>%5.1f%%</td></tr>\n', 'Not yet started', nNotYetStarted, nNotYetStarted/nSimulations*100)];
            summary = [summary sprintf('</table>\n<br/>\n')];
            
            %times
            meanMassDoublingTime = 9.0 / nanmean(cumGrowths);
            stdMassDoublingTime = sqrt(9.0 / nanmean(cumGrowths)) * nanstd(cumGrowths);
            summary = [summary sprintf('<table class="summary" cellpadding="0" cellspacing="0">\n')];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.1f &plusmn; %5.1f</td></tr>\n', 'Mass doubling time (h)', meanMassDoublingTime, stdMassDoublingTime)];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.1f &plusmn; %5.1f</td></tr>\n', 'Cell cycle length (h)', nanmean(cytTimes), nanstd(cytTimes))];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.1f &plusmn; %5.1f</td></tr>\n', 'Replication initiation time (h)', nanmean(repInitTimes), nanstd(repInitTimes))];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.1f &plusmn; %5.1f</td></tr>\n', 'Initiation of replication time (h)', nanmean(initRepTimes), nanstd(initRepTimes))];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.1f &plusmn; %5.1f</td></tr>\n', 'Replication time (h)', nanmean(repTimes), nanstd(repTimes))];
            summary = [summary sprintf('</table>\n<br/>\n')];
            
            %durations
            summary = [summary sprintf('<table class="summary" cellpadding="0" cellspacing="0">\n')];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.1f &plusmn; %5.1f</td></tr>\n', 'Replication initiation duration (h)', ...
                nanmean(repInitTimes), nanstd(repInitTimes))];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.1f &plusmn; %5.1f</td></tr>\n', 'Replication initiation-replisome formation duration (h)', ...
                nanmean(initRepTimes - repInitTimes), nanstd(initRepTimes - repInitTimes))];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.1f &plusmn; %5.1f</td></tr>\n', 'Replication duration (h)', ...
                nanmean(repTimes - initRepTimes), nanstd(repTimes - initRepTimes))];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.1f &plusmn; %5.1f</td></tr>\n', 'Cytokinesis duration (h)', ...
                nanmean(cytTimes - repTimes), nanstd(cytTimes - repTimes))];
            summary = [summary sprintf('</table>\n<br/>\n')];
            
            %growth
            summary = [summary sprintf('<table class="summary" cellpadding="0" cellspacing="0">\n')];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.3f &plusmn; %5.3f</td></tr>\n', 'Growth at t==0 (1/cells/h)', ...
                nanmean(initGrowths), nanstd(initGrowths))];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.3f &plusmn; %5.3f</td></tr>\n', 'Growth at t==&tau; (1/cells/h)', ...
                nanmean(finGrowths), nanstd(finGrowths))];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.3f &plusmn; %5.3f</td></tr>\n', 'Cumulative Growth over &tau; (1/cells)', ...
                nanmean(cumGrowths), nanstd(cumGrowths))];
            summary = [summary sprintf('</table>\n<br/>\n')];
            
            %mass
            summary = [summary sprintf('<table class="summary" cellpadding="0" cellspacing="0">\n')];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.3f &plusmn; %5.3f</td></tr>\n', 'Mean Mass at t==0 (cells)', ...
                nanmean(initMass), nanstd(initMass))];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.3f</td></tr>\n', 'Median Mass at t==0; (cells)', ...
                nanmedian(initMass))];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.3f &plusmn; %5.3f</td></tr>\n', 'Mass Mass at t==&tau; (cells)', ...
                nanmean(finMass), nanstd(finMass))];
            summary = [summary sprintf('<tr><th>%s</th><td>%5.3f</td></tr>\n', 'Median at t==&tau; (cells)', ...
                nanmedian(finMass))];
            summary = [summary sprintf('</table>\n')];
            
            summary = [summary sprintf('</body>\n')];
            summary = [summary sprintf('</html>\n')];
        end
    end
end
