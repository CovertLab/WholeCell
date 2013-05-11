%Transcription
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/6/2011
classdef Transcription
    %plotting
    methods (Static = true)
        function plotRNA_Polymerase_States(sim, axesHandle, time, ~)
            import edu.stanford.covert.util.ConstantUtil;
            
            rnaPol = sim.state('RNAPolymerase');
            
            maxState = max(rnaPol.transcripts.transcriptionUnitLengths);
            if ndims(rnaPol.states) == 3
                colors = zeros([size(rnaPol.states) 3]);
            else
                colors = zeros([size(rnaPol.states) 1 3]);
            end
            for i = 1:size(rnaPol.states, 1)
                for k = 1:size(rnaPol.states, 3)
                    switch rnaPol.states(i, 1, k)
                        case rnaPol.notExistValue
                            colors(i,1,k,:) = [1 0 0];
                        case rnaPol.freeValue
                            colors(i,1,k,:) = [0 1 0];
                        case rnaPol.nonSpecificallyBoundValue
                            colors(i,1,k,:) = [0 0 1];
                        otherwise
                            colors(i,1,k,:) = repmat(max(0, rnaPol.states(i,1,k)) / maxState, 1, 3);
                    end
                end
            end
            
            colors = reshape(colors, [], size(rnaPol.states, 3), 3);
            image(colors, 'Parent', axesHandle);
            
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'RNA Polymerase', 'fontSize', 8);
            
            xTick = get(axesHandle, 'XTick');
            xTickLabel = round((xTick - 0.5) * max(time) / ConstantUtil.secondsPerHour / max(xTick));
            set(axesHandle, 'XTickLabel', xTickLabel);
        end
        
        function plotRNA_Polymerase_State_Occupancies(sim, axesHandle, time, ~)
            rnaPol = sim.state('RNAPolymerase');
            
            labels = cell(4,1);
            labels{rnaPol.activelyTranscribingIndex} = 'Actively Transcribing';
            labels{rnaPol.specificallyBoundIndex} = 'Specifically Bound';
            labels{rnaPol.nonSpecificallyBoundIndex} = 'Non-specifically Bound';
            labels{rnaPol.freeIndex} = 'Free';
            
            plotHandles = plot(axesHandle, time / 3600, permute(rnaPol.stateOccupancies, [3 1 2]));
            legend(plotHandles, labels);
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'Occupancy', 'fontSize', 8);
            xlim(axesHandle, [0 max(1, time(end))] / 3600);
        end
        
        function plotRNA_Polymerases(sim, axesHandle, time, ~)
            rnaPol = sim.state('RNAPolymerase');
            
            rnaPolymeraseGeneStates = rnaPol.states;
            rnaPolymeraseGenomeStates = zeros(size(rnaPolymeraseGeneStates));
            for i = 1:size(rnaPolymeraseGeneStates,1)
                for k = 1:size(rnaPolymeraseGeneStates,3)
                    switch rnaPolymeraseGeneStates(i, 1, k)
                        case rnaPol.specificallyBoundValue
                            rnaPolymeraseGenomeStates(i, 1, k) = rnaPol.transcripts.transcriptionUnitFivePrimeCoordinates(rnaPol.transcripts.boundTranscriptionUnits(i,1,k));
                        case rnaPol.activelyTranscribingValue
                            rnaPolymeraseGenomeStates(i, 1, k) = rnaPol.transcripts.transcriptionUnitFivePrimeCoordinates(rnaPol.transcripts.boundTranscriptionUnits(i,1,k))+...
                                (2*rnaPol.transcripts.transcriptionUnitDirections(rnaPol.transcripts.boundTranscriptionUnits(i,1,k))-1)*(rnaPolymeraseGeneStates(i,1,k)-1);
                    end
                end
            end
            
            plot(axesHandle, time/3600, permute(rnaPolymeraseGenomeStates, [3 1 2]), '.');
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'Genome Coordinate', 'fontSize', 8);
            xlim(axesHandle, [0 max(1, time(end))] / 3600);
            ylim(axesHandle, [rnaPol.freeValue-1 rnaPol.transcripts.genomeLength+1]);
        end
        
        function plotRNA_Polymerase_Bound_Transcription_Units(sim, axesHandle, time, ~)
            rnaPol = sim.state('RNAPolymerase');
            
            numTimePoints = size(rnaPol.states,3);
            
            polymeraseTUs = rnaPol.transcripts.boundTranscriptionUnits;
            boundTUs = zeros(length(rnaPol.transcripts.transcriptionUnitLengths), 1, numTimePoints);
            for i = 1:size(polymeraseTUs,1)
                for k = 1:size(polymeraseTUs,3)
                    if polymeraseTUs(i,1,k) > 0
                        boundTUs(polymeraseTUs(i,1,k),1,k) = ...
                            boundTUs(polymeraseTUs(i,1,k),1,k) + 1;
                    end
                end
            end
            
            boundTUs = permute(boundTUs, [1 3 2]);
            colors = repmat(boundTUs/(max(max(boundTUs))),[1 1 3]);
            
            image(colors,'Parent',axesHandle);
            xlabel(axesHandle,'Time (h)','fontSize',8);
            ylabel(axesHandle,'Transcription Unit','fontSize',8);
            
            xTick = get(axesHandle,'XTick');
            xTickLabel = round((xTick-0.5)*max(time)/3600/max(xTick));
            set(axesHandle,'XTickLabel',xTickLabel);
        end
        
        function plotNTPs(sim, axesHandle, time, ~)
            p = sim.process('Transcription');
            
            plotHandles = plot(axesHandle, time / 3600, permute(p.substrates(p.substrateIndexs_ntp, :, :), [3 1 2]));
            legend(plotHandles, p.substrateWholeCellModelIDs(p.substrateIndexs_ntp));
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'NTPs', 'fontSize', 8);
            xlim(axesHandle, [0 max(1, time(end))] / 3600);
        end
        
        function plotNMPs(sim, axesHandle, time, ~)
            p = sim.process('Transcription');
            
            plotHandles = plot(axesHandle, time / 3600, permute(p.substrates(p.substrateIndexs_nmp, :, :), [3 1 2]));
            legend(plotHandles, p.substrateWholeCellModelIDs(p.substrateIndexs_nmp));
            xlabel(axesHandle,'Time (h)', 'fontSize', 8);
            ylabel(axesHandle,'NMPs', 'fontSize', 8);
            xlim(axesHandle, [0 max(1, time(end))] / 3600);
        end
        
        %plot translation limits
        function plotTranscription_Limits(sim, axesHandle, time, ~)
            p = sim.process('Transcription');
            
            plotHandles = plot(axesHandle, time / 3600, permute(p.substrates(p.substrateIndexs_ntp), [3 1 2]));
            legend(plotHandles, p.substrateWholeCellModelIDs(p.substrateIndexs_ntp));
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'NTP', 'fontSize', 8);
            xlim(axesHandle, [0 max(1, time(end))] / 3600);
        end
        
        function plotBound_Transcription_Factors(sim, axesHandle, time, ~)
            p = sim.process('Transcription');
            
            transcriptionFactors      = p.enzymes(p.enzymeIndexs_transcriptionFactors, :, :);
            boundTranscriptionFactors = p.boundEnzymes(p.enzymeIndexs_transcriptionFactors, :, :);
            
            plotHandles = plot(axesHandle, time / 3600, permute(boundTranscriptionFactors./(transcriptionFactors+boundTranscriptionFactors), [3 1 2]));
            legend(plotHandles, p.enzymeNames(p.enzymeIndexs_transcriptionFactors));
            xlabel(axesHandle, 'Time (h)', 'fontSize', 8);
            ylabel(axesHandle, 'Fraction Bound', 'fontSize', 8);
            xlim(axesHandle, [0 max(1, time(end))] / 3600);
            ylim(axesHandle, [0 1]);
        end
    end
end