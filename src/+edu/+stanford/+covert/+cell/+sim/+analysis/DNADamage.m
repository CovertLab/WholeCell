%DNADamage
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 4/19/2011
classdef DNADamage
    methods (Static = true)
        function run(sim, fileName)
            import edu.stanford.covert.cell.sim.analysis.DNADamage;
            import edu.stanford.covert.cell.sim.util.PrintUtil;
            
            %excel file
            [content, colLabels, indentation] = DNADamage.calcExpectedRates(sim);
            if nargin == 1
                PrintUtil.printToStdIO(content, colLabels, struct('indentation', indentation));
            else
                PrintUtil.printToFile(content, colLabels, [fileName '.xls'], 'State', struct('indentation', indentation));
            end
        end
        
        function [content, colLabels, indentation] = calcExpectedRates(sim)
            d = sim.process('DNADamage');
                        
            reactionVulnerableMotifs = d.reactionVulnerableMotifs;
            tfs = ismember(d.reactionVulnerableMotifTypes, {'damagedBases','damagedSugarPhosphates','intrastrandCrossLinks'});
            reactionVulnerableMotifs(tfs) = ...
                d.metabolite.wholeCellModelIDs(cell2mat(d.reactionVulnerableMotifs(tfs)));
            tfs = ismember(d.reactionVulnerableMotifTypes, {'abasicSites','gapSites','hollidayJunctions', 'strandBreaks'});
            reactionVulnerableMotifs(tfs) = cell(sum(tfs), 1);
            
            rates = d.calcExpectedReactionRates();
            
            reactionRadiation = cell(size(d.reactionRadiation));
            reactionRadiationLevels = zeros(size(d.reactionRadiation));
            reactionRadiation(d.reactionRadiation~=0) = ...
                d.substrateWholeCellModelIDs(d.reactionRadiation(d.reactionRadiation~=0));
            reactionRadiationLevels(d.reactionRadiation~=0) = ...
                d.substrates(d.reactionRadiation(d.reactionRadiation~=0));
            
            colLabels = {'ID', 'Name', ...
                'Damage Type', 'Damage Motif', 'No. Vulnerable Site', ...
                'Specific Rate (eg 1/s/nt/gy)', 'Rate (1/s/cell)', 'Rate (1/cell cycle)', ...
                'Radiation ID', 'Radiation Value'};
            
            content = [...
                d.reactionWholeCellModelIDs d.reactionNames ...
                d.reactionVulnerableMotifTypes reactionVulnerableMotifs num2cell(d.calcNumberVulnerableSites()) ...
                num2cell(d.reactionBounds(:, 2)) num2cell(rates) num2cell(rates * ((9.0 + 2) * 3600)) ...
                reactionRadiation num2cell(reactionRadiationLevels) ...
                ];
            
            indentation = zeros(size(content, 1), 1);
        end
    end
end