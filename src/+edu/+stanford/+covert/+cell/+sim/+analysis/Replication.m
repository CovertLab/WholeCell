%Replication
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 3/2/2011
classdef Replication
    properties (Constant)
        verbose = false;
    end
    
    methods (Static)
        function run(outputFileName)
            import edu.stanford.covert.cell.sim.analysis.Replication;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
                        
            % calculate
            [time primases leadingPols laggingPols ssbs ligases segregation] = Replication.simulateReplication();
            durations = Replication.sampleDuration();
            
            % save
            save([outputFileName '.mat'], ...
                'time', 'primases', ...
                'leadingPols', 'laggingPols', ...
                'ssbs', 'ligases', 'segregation', ...
                'durations');            
            
            % print
            colLabels = {...
                'Time', 'Primase-1', 'Primase-2', ...
                'LeadingPol-1', 'LeadingPol-2', 'LaggingPol-1', 'LaggingPol-2', ...
                'SSBs-1', 'SSBs-2', 'Ligase-1', 'Ligase-2', 'Segregation'};
            if nargin == 0
                PrintUtil.printToStdIO(num2cell(durations), colLabels);
            else
                PrintUtil.printToFile(num2cell([time primases leadingPols laggingPols ssbs ligases segregation]), colLabels, [outputFileName '.xls'], 'Statistics');
            end
            
            colLabels = {'Duration'};
            PrintUtil.printToFile(num2cell(durations), colLabels, [outputFileName '.xls'], 'Duration');
                        
            % plot
            if nargin == 1, [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle(); end
            
            maxTime = find(segregation);
            
            if nargin == 0, [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle(); end
            plot(axesHandle, time/60, primases)
            xlim([0 maxTime]/60);
            xlabel('Time (m)', 'fontsize', 12)
            ylabel('Primase Activity', 'fontsize', 12);
            if nargin == 1, saveas(figHandle, [outputFileName '-primases.pdf']); end
            
            if nargin == 0, [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle(); end
            plot(axesHandle, time/60, leadingPols)
            xlim([0 maxTime]/60);
            xlabel('Time (m)', 'fontsize', 12)
            ylabel('Leading Position', 'fontsize', 12);
            if nargin == 1, saveas(figHandle, [outputFileName '-leadingPols.pdf']); end
            
            if nargin == 0, [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle(); end
            plot(axesHandle, time/60, laggingPols)
            xlim([0 maxTime]/60);
            xlabel('Time (m)', 'fontsize', 12)
            ylabel('Lagging Position', 'fontsize', 12);
            if nargin == 1, saveas(figHandle, [outputFileName '-laggingPols.pdf']); end
            
            if nargin == 0, [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle(); end
            plot(axesHandle, time/60, ssbs)
            xlim([0 maxTime]/60);
            xlabel('Time (m)', 'fontsize', 12)
            ylabel('Bound SSBs', 'fontsize', 12);
            if nargin == 1, saveas(figHandle, [outputFileName '-ssbs.pdf']); end
            
            if nargin == 0, [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle(); end
            plot(axesHandle, time/60, ligases)
            xlim([0 maxTime]/60);
            xlabel('Time (m)', 'fontsize', 12)
            ylabel('Ligase Activity', 'fontsize', 12);
            if nargin == 1, saveas(figHandle, [outputFileName '-ligases.pdf']); end
            
            if nargin == 0, [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle(); end
            plot(axesHandle, time/60, segregation)
            xlim([0 maxTime]/60);
            xlabel('Time (m)', 'fontsize', 12)
            ylabel('Segregation Activity', 'fontsize', 12);
            if nargin == 1, saveas(figHandle, [outputFileName '-segregation.pdf']); end
            
            if nargin == 0, [axesHandle, figHandle] = edu.stanford.covert.cell.sim.util.PlotUtil.newAxesHandle(); end
            hist(axesHandle, durations / 60);
            xlabel('Duration (m)', 'fontsize', 12);
            ylabel('Frequency', 'fontsize', 12);
            if nargin == 1, saveas(figHandle, [outputFileName '-segregation.pdf']); end
            
            if nargin == 1, close(figHandle); end
        end
    end
    
    methods (Static)
        function durations = sampleDuration()
            import edu.stanford.covert.cell.sim.analysis.Replication;
            
            nTrials = 0;
            durations = zeros(nTrials, 1);
            parfor i = 1:nTrials
                if Replication.verbose
                    fprintf('Trial = %d\n', i); 
                end
                
                [~, ~, ~, ~, ~, ~, segregation] = Replication.simulateReplication(i);
                durations(i) = find(segregation);
            end
        end
        
        function [time, primases, leadingPols, laggingPols, ssbs, ligases, segregation] = simulateReplication(seed)            
            import edu.stanford.covert.cell.sim.analysis.Replication;
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;
            import edu.stanford.covert.util.CircularSparseMat;
            
            m = edu.stanford.covert.cell.sim.ProcessFixture.load(edu.stanford.covert.cell.sim.process.Replication([], 'Replication'));
            m.verbosity = 0;
            if ~exist('seed', 'var')
                seed = 1;
            end
            m.seed = seed;
            m.seedRandStream();
            
            c = m.chromosome;           
            
            % initial state
            c.initialize();
            c.linkingNumbers(1, 1:2) = 0;
            
            m.substrates(:) = 0;
            m.substrates(m.substrateIndexs_dntp) = [
                sum(sum(c.sequence(:, 2:3) == 'A'));
                sum(sum(c.sequence(:, 2:3) == 'C'));
                sum(sum(c.sequence(:, 2:3) == 'G'));
                sum(sum(c.sequence(:, 2:3) == 'T'))];
            m.substrates(m.substrateIndexs_atp) = ...
                + size(c.sequence, 1) ...                              %helicase unwinding each base pair
                + sum(cellfun(@numel, m.primaseBindingLocations)) ...  %lagging strand beta-clamp formation
                + 2;                                                   %leading strand beta-clamp formation
            m.substrates(m.substrateIndexs_water) = m.substrates(m.substrateIndexs_atp);
            m.substrates(m.substrateIndexs_nad) = sum(cellfun(@numel,m.primaseBindingLocations)) + 2;
            
            m.enzymes = m.enzymeComposition(:, [m.enzymeIndexs_replisome; m.enzymeIndexs_betaClamp; m.enzymeIndexs_ssb8mer]) * [2; 2; 100];
            m.enzymes(m.enzymeIndexs_ligase) = 100;
            m.boundEnzymes(:) = 0;
            
            %set up OriC complex
            c.complexBoundSites(m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R1234), 1) = m.complexIndexs_DnaA_7mer_ATP;
            c.complexBoundSites(m.dnaAFunctionalBoxStartPositions(m.dnaAFunctionalBoxIndexs_R5), 1) = m.complexIndexs_DnaA_1mer_ATP;
            
            % evolve state
            nTime = 2*3600;
            time = (1:nTime)';
            primases = zeros(nTime, 2);
            leadingPols = zeros(nTime, 2);
            laggingPols = zeros(nTime, 2);
            ssbs = zeros(nTime, 2);
            ligases = zeros(nTime, 2);
            segregation = zeros(nTime, 1);
            
            for i = 1:nTime
                if Replication.verbose && mod(i, 100) == 1
                    fprintf('t = %d\n', i); 
                end

                oldStrandBreaks = c.strandBreaks;
                
                m.evolveState();
                                
                primases(i, :) = ...
                    + ((m.okazakiFragmentIndex > 0) & (m.okazakiFragmentProgress <= m.primerLength)) ...
                    + (m.leadingPosition <= m.primerLength);
                
                leadingPols(i, :) = m.leadingPosition;
                laggingPols(i, :) = m.laggingPosition;
                
                ssbs(i, :) = ...
                    + m.numLeadingTemplateBoundSSBs ...
                    + m.numLaggingTemplateBoundSSBs;
                
                if nnz(oldStrandBreaks)
                    subs1 = find(oldStrandBreaks);
                    subs2 = find(c.strandBreaks);
                    
                    ligases(i, 1) = sum(~ismembc(subs1(subs1(:, 2)==2, 1), subs2(subs2(:,2)==2, 1)));
                    ligases(i, 2) = sum(~ismembc(subs1(subs1(:, 2)==3, 1), subs2(subs2(:,2)==3, 1)));
                end
                
                if all(m.strandDuplicated)
                    break;
                end
            end            
            segregation(i+1, :) = 1;
            
            %% assertions
            assertEqual([true true], m.strandDuplicated);            
        end
    end
end