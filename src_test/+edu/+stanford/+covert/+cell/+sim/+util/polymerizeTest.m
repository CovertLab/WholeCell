%polymerize test
%
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updaetd 2/3/2011
classdef polymerizeTest < TestCase
    methods
        function this = polymerizeTest(name)
            this = this@TestCase(name);
        end

        function testNoSequences(~)
            [progress, baseAmounts, baseCosts, energy, energyCost] = ...
                edu.stanford.covert.cell.sim.util.polymerize(...
                    [], [9 9 9 9]', 'ACGU', ' ', 0, 0, edu.stanford.covert.util.RandStream('mcg16807'));
            assertEqual(zeros(0, 1), progress);
            assertEqual([9 9 9 9]', baseAmounts);
            assertEqual([0 0 0 0]', baseCosts);
            assertEqual(0, energy);
            assertEqual(0, energyCost);
        end

        function testEnergyLimited(~)
            [progress, baseAmounts, baseCosts, energy, energyCost] = ...
                edu.stanford.covert.cell.sim.util.polymerize(...
                    ['ABAB';'BBBB'], [9 9]', 'AB', ' ', 6, 2,...
                    RandStream('mlfg6331_64', 'Seed', 0));
            assertEqual([2 1]', progress);
            assertEqual([8 7]', baseAmounts);
            assertEqual([1 2]', baseCosts);
            assertEqual(0, energy);
            assertEqual(6, energyCost);
        end

        function testFairness(~)
            [progress, baseAmounts, baseCosts, energy, energyCost] = ...
                edu.stanford.covert.cell.sim.util.polymerize([
                    'ABCDABCDABCD';
                    'DDDDDDDDDDDD';
                    'CCCCCCCCCCCC';
                    '            '],...
                    [11 11 11 11]', 'ABCD', ' ', 30, 1, edu.stanford.covert.util.RandStream('mcg16807'));
            assertEqual([10 9 9 0]', progress);
            assertEqual([8 8 0 0]', baseAmounts);
            assertEqual([3 3 11 11]', baseCosts);
            assertEqual(2, energy);
            assertEqual(28, energyCost);
        end

        function testAllAvailableBases(~)
            [progress, baseAmounts, baseCosts, energy, energyCost] = ...
                edu.stanford.covert.cell.sim.util.polymerize( ...
                    ['ABCDABCDABBB';...
                     'DDDDDDDDDDDD';...
                     'CCCCCCCCCCCC'],...
                     [11 11 11 11]', 'ABCD', ' ', 30, 1, edu.stanford.covert.util.RandStream('mcg16807'));
            assertEqual([12 9 9]', progress);
            assertEqual([8 6 0 0]', baseAmounts);
            assertEqual([3 5 11 11]', baseCosts);
            assertEqual(0, energy);
            assertEqual(30, energyCost);
        end
        
        function testNoElongationOfFullyPolymerizedSequences(~)
            %polymerize with saturating tRNA, limiting energy, and 1
            %sequence that has already been fully polymerized
            nTRNAs = 36;
            energyCostPerBase = 2;
            elngSeqs = [
                0     0
                24    21
                29    31
                ];
            energy = 3;
            aminoacylatedTRNAs = 50 * ones(36, 1);
            
            randStream = edu.stanford.covert.util.RandStream('mcg16807');
            randStream.reset(round(mod(now, 1) * 1e7));
           
            [translationProgress, newAminoacylatedTRNAs, aminoacylatedTRNACosts, newEnergy, energyCost] = ...
                edu.stanford.covert.cell.sim.util.polymerize(...
                elngSeqs, aminoacylatedTRNAs, 1:nTRNAs, 0, energy, energyCostPerBase, randStream);
            
            %no elongation of fully polymerized sequences
            assertAllEqual(0, translationProgress(elngSeqs(:, 1) == 0));
            
            %base costs correctly calculated
            tmp = zeros(size(aminoacylatedTRNACosts));
            for i = 1:numel(translationProgress)
                tmp = tmp + reshape(histc(elngSeqs(i, 1:translationProgress(i)), 1:nTRNAs), [], 1);
            end
            assertEqual(aminoacylatedTRNACosts, aminoacylatedTRNAs-newAminoacylatedTRNAs);
            assertEqual(aminoacylatedTRNACosts', tmp');
            
            %energy costs correctly calculated
            assertEqual(energyCost, energy-newEnergy);
            assertEqual(energyCost, energyCostPerBase * sum(aminoacylatedTRNACosts));
        end
    end
end
