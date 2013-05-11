% ChromosomeSequence
%   Object to represent the circular, complementary nature of chromosomal
%   sequences. Externally this object appears as a circular len X 4 character
%   array, where len is the length of the positive strand sequence, and the four
%   columns represent: (1) mother chromosome, positive strand, (2) mother
%   chromosome, negative strand, (3), daughter chromosome, positive strand,
%   (4) daughter chromosome, negative strand. Subscript referencing is
%   automatically wrapped over both dimensions so that
%     (a) x(-1,1) references to x(end-1,1) and
%     (b) x(1,5) references to x(1,1).
%
%   This object behaves differently that character arrays in several ways:
%   - Subscript assignment is prohibited
%   - Reshaping, squeezing, permuting, transposing, and concatenating return
%     results NOT of type ChromosomeSequence, but of type char
%
%   Internally this class uses one property positiveStrandSequence to represent
%   the mother and daughter chromosomes positive and negative sequences.
%   Negative sequences are computed using the seqcomplement method of the
%   bioinformatics toolbox.
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/6/2010
classdef ChromosomeSequence < handle
    properties (Access = protected)
        positiveStrandSequence
        gcContent
        baseCounts
    end

    properties (Access = protected, Dependent = true)
        negativeStrandSequence
        positiveAndNegtiveStrandSequences
        positiveAndNegtiveStrandMotherAndDaughterSequences
    end

    properties (Constant = true, GetAccess = protected)
        iPositiveStrand = 1;
        iNegativeStrand = 2;
        iMotherChromosome = 1;
        iDaughterChromosome = 2;
    end

    %constructor
    methods
        function this = ChromosomeSequence(positiveStrandSequence)
            if ...
                    ~ischar(positiveStrandSequence) || ...
                    ~all(ismembc(positiveStrandSequence, 'ACGT')) || ...
                    ~isvector(positiveStrandSequence)
                throw(MException('ChromosomeSequence:invalidSequence', 'positiveStrandSequence must be a vector containing only the characters A, C, G, T'));
            end

            this.positiveStrandSequence = reshape(positiveStrandSequence,[],1);
        end
    end
    
    %getters
    methods
        function value = getGCContent(this)
            if isempty(this.gcContent)
                this.gcContent = (...
                    sum(this.positiveStrandSequence == 'G') + ...
                    sum(this.positiveStrandSequence == 'C')) / ...
                    length(this.positiveStrandSequence);
            end
            
            value = this.gcContent;
        end
        
        function value = getBaseCounts(this, strand)
            if isempty(this.baseCounts)
                this.baseCounts = [
                    sum(this.positiveStrandSequence == 'A') sum(this.positiveStrandSequence == 'T')
                    sum(this.positiveStrandSequence == 'C') sum(this.positiveStrandSequence == 'G')
                    sum(this.positiveStrandSequence == 'G') sum(this.positiveStrandSequence == 'C')
                    sum(this.positiveStrandSequence == 'T') sum(this.positiveStrandSequence == 'A')];
            end
            
            if nargin == 1
                value = sum(this.baseCounts, 2);
            else
                value = this.baseCounts(:, mod(strand - 1, 2) + 1);
            end
        end
    end

    %getters
    methods
        function seq = get.negativeStrandSequence(this)
            seq = seqcomplement(this.positiveStrandSequence')';
        end

        function seq = get.positiveAndNegtiveStrandSequences(this)
            seq = [this.positiveStrandSequence this.negativeStrandSequence];
        end

        function seq = get.positiveAndNegtiveStrandMotherAndDaughterSequences(this)
            seq = repmat(this.positiveAndNegtiveStrandSequences, 1, 2);
        end
    end

    %size
    methods
        function val = isempty(this)
            val = isempty(this.positiveStrandSequence);
        end

        function val = numel(this)
            val = prod(size(this)); %#ok<PSIZE>
        end

        function val = length(this)
            val = max(size(this));
        end

        function val = size(this, dim)
            if exist('dim','var')
                if dim < 1 || ~isequal(dim, ceil(real(dim))) || numel(dim) ~= 1
                    throw(MException('ChromosomeSequence:invalidDimension', 'dim must be a positive real integer'));
                elseif dim == 1
                    val = numel(this.positiveStrandSequence);
                elseif dim == 2
                    val = 4;
                else
                    val = 1;
                end
            else
                val = [numel(this.positiveStrandSequence) 4];
            end
        end

        function val = ndims(this)
            val = numel(size(this));
        end
    end

    %size methods that return character arrays
    methods
        function seq = reshape(this, varargin)
            seq = reshape(this.positiveAndNegtiveStrandMotherAndDaughterSequences, varargin{:});
        end

        function seq = transpose(this)
            seq = transpose(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
        end

        function seq = ctranspose(this)
            seq = transpose(this);
        end

        function seq = permute(this, dims)
            seq = permute(this.positiveAndNegtiveStrandMotherAndDaughterSequences, dims);
        end

        function seq = squeeze(this)
            seq = squeeze(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
        end
    end

    %subscript reference, assignment
    %return character arrays
    methods
        function varargout = subsref(this, s)
            varargout = cell(max(1, nargout), 1);
            if strcmp(s(1).type, '()')
                if numel(s.subs) == 1
                    varargout{1} = subsref(this.positiveAndNegtiveStrandMotherAndDaughterSequences, s);
                else
                    tmp = s;
                    if isnumeric(s.subs{1})
                        s.subs{1} = mod(s.subs{1} - 1, size(this, 1)) + 1;
                        tmp.subs{1} = s.subs{1};
                    end
                    if isnumeric(s.subs{2})
                        s.subs{2} = mod(s.subs{2}, 2) == 0;
                        tmp.subs{2} = ones(size(s.subs{2}));
                    else
                        s.subs{2} = repmat([false; true], 2, 1);
                        tmp.subs{2} = ones(4, 1);
                    end

                    varargout{1} = subsref(this.positiveStrandSequence, tmp);
                    siz = size(varargout{1});
                    siz(2) = sum(s.subs{2});
                    varargout{1}(:, s.subs{2}, :) = reshape(seqcomplement(reshape(varargout{1}(:, s.subs{2}, :), 1, [])), siz);
                end
            else
                [varargout{:}] = builtin('subsref', this, s);
            end
        end

        function val = subsasgn(this, s, rhs)
            if strcmp(s(1).type, '()')
                throw(MException('ChromosomeSequence:invalidSyntax', 'undefined syntax'));
            else
                val = builtin('subsasgn', this, s, rhs);
            end
        end
        
        function e = end(this, k, n)
            if n == 1
                e = numel(this);
            else
                e = size(this, k);
            end
        end
    end

    %concatentation: return character arrays NOT ChromosomeSequences
    methods
        function seq = cat(dim, varargin)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            for i = 1:numel(varargin)
                if isa(varargin{i}, 'ChromosomeSequence')
                    varargin{i} = varargin{i}.positiveAndNegtiveStrandMotherAndDaughterSequences;
                end
            end

            seq = cat(dim, varargin{:});
        end

        function seq = horzcat(varargin)
            seq = cat(2, varargin{:});
        end

        function seq = vertcat(varargin)
            seq = cat(1, varargin{:});
        end
    end

    methods
        function val = isequal(A, B)
            val = ...
                isequal(class(A), class(B)) && ...
                isequal(A.positiveStrandSequence, B.positiveStrandSequence);
        end

        function C = eq(A, B)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            if isa(A, 'ChromosomeSequence')
                if isa(B, 'ChromosomeSequence')
                    C = eq(A.positiveAndNegtiveStrandMotherAndDaughterSequences, B.positiveAndNegtiveStrandMotherAndDaughterSequences);
                else
                    C = eq(A.positiveAndNegtiveStrandMotherAndDaughterSequences, B);
                end
            else
                C = eq(A, B.positiveAndNegtiveStrandMotherAndDaughterSequences);
            end
        end

        function C = ne(A, B)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            if isa(A, 'ChromosomeSequence')
                if isa(B, 'ChromosomeSequence')
                    C = ne(A.positiveAndNegtiveStrandMotherAndDaughterSequences, B.positiveAndNegtiveStrandMotherAndDaughterSequences);
                else
                    C = ne(A.positiveAndNegtiveStrandMotherAndDaughterSequences, B);
                end
            else
                C = ne(A, B.positiveAndNegtiveStrandMotherAndDaughterSequences);
            end
        end
    end

    %find
    methods
        function varargout = find(this)
            varargout = cell(nargout,1);
            switch nargout
                case 1,
                    if size(this, 1) == 1 && ndims(this) == 2
                        varargout{1} = 1:numel(this);
                    else
                        varargout{1} = (1:numel(this))';
                    end
                case 2
                    if size(this, 1) == 1 && ndims(this) == 2
                        varargout{1} = ones(1, size(this, 2));
                        varargout{2} = 1:size(this, 2);
                    else
                        varargout{1} = repmat((1:size(this, 1))', size(this, 2), 1);
                        varargout{2} = sort(repmat((1:size(this, 2))', size(this, 1), 1));
                    end
                case 3
                    if size(this, 1) == 1 && ndims(this) == 2
                        varargout{1} = ones(1, size(this, 2));
                        varargout{2} = 1:size(this, 2);
                        varargout{3} = double(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
                    else
                        varargout{1} = repmat((1:size(this,1))', size(this,2), 1);
                        varargout{2} = sort(repmat((1:size(this,2))', size(this,1), 1));
                        varargout{3} = double(reshape(this.positiveAndNegtiveStrandMotherAndDaughterSequences,[],1));
                    end
                otherwise
                    throw(MException('ChromosomeSequence:tooManyOutputs','find produces at most 3 outputs'));
            end
        end

        function positionStrands = findSubsequence(this, subsequence, subposition, bothStrands)
            if isempty(subposition) || isempty(this.positiveStrandSequence)
                positionStrands = [];
                return;
            end

            sequence = this.positiveStrandSequence;
            [~, positions1] = restrict([sequence; sequence(1:length(subsequence)-1)]', subsequence, subposition);
            if length(positions1) == 1 && positions1 == 0
                positions1 = zeros(0, 1);
            else
                positions1 = positions1(2:end);
            end

            if ~exist('bothStrands', 'var') || ~bothStrands
                positionStrands = positions1;
                return;
            end

            sequence = this.negativeStrandSequence;
            [~, positions2] = restrict([sequence(end-length(subsequence) + 2:end); sequence]', seqreverse(subsequence), subposition);
            if length(positions2) == 1 && positions2 == 0
                positions2 = zeros(0, 1);
            else
                positions2 = positions2(2:end);
            end
            
            positionStrands = [
                positions1     ones(size(positions1));
                positions2 2 * ones(size(positions2))];
        end

        %L = length of chromosome
        %S = number of strands
        %F = number of features
        %l = length of longest feature
        %
        %features:  L x S boolean indicating positions of some feature in chromosome(s)
        %positions: F x 1 integer indicating coordinates of features
        %strands:   F x 1 integer indicating strands of features
        %sequences: F x l char indicating sequences of features; vector is char(0) right padded
        function [allPositions, allLengths, allStrands, allSequences] = locateFeatures(this, features)
            allPositions = zeros(0, 1);
            allLengths   = zeros(0, 1);
            allStrands   = zeros(0, 1);
            allSequences = cell(0, 1);

            %loop over strands
            for i = 1:max(features(:, 2))
                positions = sort(features(features(:, 2) == i, 1));
                lengths   = ones(size(positions));
                seq = this.positiveStrandSequence(positions, 1);
                if mod(i, 2) == 0 && ~isempty(seq)
                    seq = seqcomplement(seq')';
                end
                sequences = num2cell(seq);

                if isempty(positions);
                    continue;
                end

                for j = length(positions):-1:2
                    if diff(positions(j-1:j)) == 1
                        lengths(j-1) = lengths(j-1) + lengths(j);
                        if i == 1
                            sequences{j-1} = [sequences{j-1} sequences{j}];
                        else
                            sequences{j-1} = [sequences{j} sequences{j-1}];
                        end
                        lengths(j) = 0;
                        sequences{j} = '';
                    end
                end

                positions = positions(lengths > 0);
                sequences = sequences(lengths > 0);
                lengths   = lengths(  lengths > 0);

                if positions(1) == 1 && positions(end) + lengths(end)-1 == size(this.positiveStrandSequence, 1)
                    if i == 1
                        lengths(end) = lengths(end) + lengths(1);
                        sequences{end} = [sequences{end} sequences{1}];
                        positions(1) = [];
                        sequences(1) = [];
                        lengths(1) = [];
                    else
                        lengths(1)   = lengths(1) + lengths(end);
                        sequences{1} = [sequences{1} sequences{end}];
                        positions(1) = positions(end);
                        positions(end) = [];
                        sequences(end) = [];
                        lengths(end) = [];
                    end
                end

                if i == 2
                    positions = mod(positions + lengths -1 -1, size(this.positiveStrandSequence,1))+1;
                end

                allPositions = [allPositions; positions];                %#ok<AGROW>
                allLengths   = [allLengths; lengths];                    %#ok<AGROW>
                allStrands   = [allStrands; repmat(i, size(positions))]; %#ok<AGROW>
                allSequences = [allSequences; sequences];                %#ok<AGROW>
            end
        end

        function subsequence = subsequence(this, positions, strands)
            if ~exist('strands', 'var')
                strands = positions(:, 2);
                positions = positions(:, 1);
            end

            if isempty(positions)
                subsequence = char(zeros(1, 0));
                return;
            end

            if ~isequal(size(strands), size(positions))
                strands = repmat(strands, size(positions) ./ size(strands));
            end

            positions = mod(positions - 1, size(this.positiveStrandSequence, 1)) + 1;
            subsequence = reshape(this.positiveStrandSequence(sub2ind(...
                size(this.positiveStrandSequence), positions, ones(size(positions)))), ...
                size(positions));
            subsequence(mod(strands, 2) == 0) = seqcomplement(reshape(subsequence(mod(strands, 2) == 0), 1, []));
        end

        function baseCounts = subsequenceBaseCounts(this, varargin)
            subsequence = reshape(this.subsequence(varargin{:}), 1, []);
            baseCounts = [
                sum(subsequence == 'A');
                sum(subsequence == 'C');
                sum(subsequence == 'G');
                sum(subsequence == 'T')];
        end
    end

    %cast
    methods
        function val = cast(this, newclass)
            val = cast(this.positiveAndNegtiveStrandMotherAndDaughterSequences, newclass);
        end

        function val = char(this)
            val = this.positiveAndNegtiveStrandMotherAndDaughterSequences;
        end

        function val = uint8(this)
            val = uint8(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
        end

        function val = uint16(this)
            val = uint16(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
        end

        function val = uint32(this)
            val = uint32(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
        end

        function val = uint64(this)
            val = uint64(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
        end

        function val = int8(this)
            val = int8(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
        end

        function val = int16(this)
            val = int16(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
        end

        function val = int32(this)
            val = int32(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
        end

        function val = int64(this)
            val = int64(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
        end

        function val = single(this)
            val = single(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
        end

        function val = double(this)
            val = double(this.positiveAndNegtiveStrandMotherAndDaughterSequences);
        end
    end

    %display
    methods
        function disp(this)
            display(this);
        end

        function display(this)
            fprintf('ChromosomeSequence of 2 chromosomes (mother and daughter), of 2 strands each, each of length %d\n', size(this, 1));
        end
    end
end