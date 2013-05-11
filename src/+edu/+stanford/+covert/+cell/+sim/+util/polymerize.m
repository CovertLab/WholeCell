% Computes the maximum polymerization for the available amounts of bases
% for all active polymers using a greedy algorithm.
%
% 1. Finds the most limited bases and the position at which they become limiting.
% 2. Finds the energy limit of polymerization.
% 3. Polymerizes all polymers up to, but not including, the limiting position.
% 4. Polymerizes polymers at the limiting position according to base
%    availability.
% 5. Culls polymers that cannot make additional progress.
% 6. Repeats 1-5 until no additional bases can be polymerized.
%
% Author: Markus Covert, mcovert@stanford.edu
% Author: Jayodita Sanghvi, jayodita@stanfod.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/10/2010
function [progress, baseAmounts, baseCosts, energy, energyCost] = polymerize(...
    sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase, randStream)
    %% optional arguments
    if nargin < 5
        energy = 0;
        energyCostPerBase = 0;
    end
    if nargin < 7
        randStream = RandStream.getDefaultStream();
    end

    %% variables
    nSequences    = size(sequences, 1);               %number of sequences
    baseCosts     = zeros(length(bases), 1);          %usage of bases, in order of bases
    energyCost    = 0;                                %usage of energy
    progress      = zeros(nSequences, 1);             %number of bases polymerized in each sequence
    activeSeqIdxs = 1:nSequences;                     %indexs of sequences still being polymerized
    seqLengths    = lengths(sequences, basePadValue); %lengths of sequences
    
    %% stop early if any baseAmounts < 0
    if any(baseAmounts < 0) || any(isnan(baseAmounts))
        return;
    end

    %% elongate sequences
    while ~isempty(sequences) && energy >= energyCostPerBase
        %eliminate sequences with nothing to polymerize
        tf = progress(activeSeqIdxs) < seqLengths(activeSeqIdxs);
        activeSeqIdxs = activeSeqIdxs(tf);
        sequences = sequences(tf, :);
        if isempty(sequences)
            break; 
        end

        %calculate limit of elongation, the number of bases each of the active
        %sequences can be elongated
        [elongation, baseUsage, limitingBases] = calculateElongationLimits(...
            sequences, bases, baseAmounts, energy, energyCostPerBase);
        
        %elongate all active sequences up to elongation
        progress(activeSeqIdxs) = min([
            progress(activeSeqIdxs) + elongation...
            seqLengths(activeSeqIdxs)], [], 2);
        baseAmounts = baseAmounts - baseUsage;
        baseCosts   = baseCosts   + baseUsage;
        energy      = energy      - energyCostPerBase * sum(baseUsage);
        energyCost  = energyCost  + energyCostPerBase * sum(baseUsage);
        sequences   = sequences(:, elongation+1:end);
        if isempty(sequences)
            break; 
        end

        if ~isempty(limitingBases)
            %eliminate sequences for which bases are not available
            tf = true(size(sequences, 1), 1);
            for i = 1:length(limitingBases)
                seqIdxs = find(sequences(:,1) == bases(limitingBases(i)));
                randOrder = randStream.randperm(numel(seqIdxs));
                tf(seqIdxs(randOrder(baseAmounts(limitingBases(i))+1:end))) = false;
            end
            activeSeqIdxs = activeSeqIdxs(tf);
            sequences = sequences(tf, :);
        elseif elongation == 0
            %limited energy and no limiting bases;
            %randomly select sequences to receive one last base
            n = fix(energy / energyCostPerBase);
            idxs = randStream.randperm(numel(activeSeqIdxs));
            idxs = idxs(1:n);
            progress(activeSeqIdxs(idxs)) = progress(activeSeqIdxs(idxs)) + 1;
            baseUsage = countBases(sequences(idxs,1), bases);
            baseUsage = baseUsage(:,2);
            baseAmounts = baseAmounts - baseUsage;
            baseCosts   = baseCosts   + baseUsage;
            energy      = energy      - energyCostPerBase * n;
            energyCost  = energyCost  + energyCostPerBase * n;
            break;
        end
    end
end

function [elongation, baseUsage, limitingBases] = calculateElongationLimits(...
    sequences, bases, baseAmounts, energy, energyCostPerBase)

    cumBaseCounts = cumsum(countBases(sequences, bases), 2);
    nBases = size(cumBaseCounts, 1);

    %individual base limits
    baseLimits = zeros(nBases,1);
    for i = 1:nBases
        baseLimits(i) = find(cumBaseCounts(i,:) <= baseAmounts(i), 1, 'last');
    end

    %energy limit
    energyLimit = find(sum(cumBaseCounts) * energyCostPerBase <= energy, 1, 'last');

    %elongate up to minimum of individual base and energy limits
    elongation = min([baseLimits; energyLimit]) - 1;
    baseUsage = cumBaseCounts(:, elongation + 1);

    %if applicable, find bases which are setting the elongation limit
    if elongation + 1 < size(cumBaseCounts, 2)
        limitingBases = find(baseLimits == elongation + 1);
    else
        limitingBases = [];
    end
end

function counts = countBases(sequences, bases)
    counts = zeros(length(bases), size(sequences,2) + 1);
    for i = 1:length(bases)
        counts(i, 2:end) = sum(sequences == bases(i), 1);
    end
end

function result = lengths(sequences, padValue)
    present = sequences ~= padValue;
    n = size(sequences, 1);
    result = zeros(n, 1);
    for i = 1:n
        result(i) = max([0 find(present(i,:), 1, 'last')]);
    end
end
