% Wrapper over random method of statisics toolbox. Changes default random
% stream to the current object prior to evaluating statistics toolbox
% method in order to draw random numbers from the current stream.
%
% Also adds methods
% - stochasticRound: rounding of floats to integers weighted by decimal
%   part
%
% Usage
% import edu.stanford.covert.util.RandStream;
% RandStream('mcg16807', 'Seed', randomSeed)
%
% Requirements:
% - MATLAB statistics toolbox
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/22/2012
classdef RandStream < handle
    properties
        defaultStream
        randStream
    end
    
    properties (Dependent = true)
        type
        seed
        numStreams
        streamIndex
        state
        substream
        randnAlg
        antithetic
        fullPrecision
    end
    
    %constructor
    methods
        function this = RandStream(type, varargin)
            this.randStream = RandStream(type, varargin{:});
        end
    end
    
    methods
        function reset(this, varargin)
            this.randStream.reset(varargin{:});
        end
    end
    
    methods
        function setDefault(this)
            if verLessThan('matlab', '7.13.0')
                %hold on to default random stream
                this.defaultStream = RandStream.getDefaultStream(); %#ok<GETRS>
                
                %make this stream the default random stream
                RandStream.setDefaultStream(this.randStream); %#ok<SETRS>
            else
                %hold on to default random stream
                this.defaultStream = RandStream.getGlobalStream();
                
                %make this stream the default random stream
                RandStream.setGlobalStream(this.randStream);
            end
        end
        
        function resetDefault(this)
            %reset default stream
            if verLessThan('matlab', '7.13.0')
                RandStream.setDefaultStream(this.defaultStream); %#ok<SETRS>
            else
                RandStream.setGlobalStream(this.defaultStream);
            end
        end
    end
    
    methods
        function value = rand(this, varargin)
            value = rand(this.randStream, varargin{:});
        end
        
        function value = randi(this, varargin)
            if numel(varargin{1})==1 && varargin{1} == 1
                siz = cell2mat(varargin(2:end));
                value = ones(siz);
                return;
            end
            value = randi(this.randStream, varargin{:});
        end
        
        function value = randn(this, varargin)
            value = randn(this.randStream, varargin{:});
        end
        
        function value = randperm(this, varargin)
            value = randperm(this.randStream, varargin{:});
        end
        
        %Wrapper over random method of statisics toolbox. Changes
        %default random stream to the current object prior to evaluating
        %statistics toolbox method in order to draw random numbers from the
        %current stream
        function value = random(this, name, varargin)
            %make this stream the default random stream
            this.setDefault();
            
            %call database toolbox method, which implicity will use this
            %stream
            value = random(name, varargin{:});
            
            %reset default stream
            this.resetDefault();
        end
        
        function value = poissrnd(this, varargin)
            %make this stream the default random stream
            this.setDefault();
            
            %call database toolbox method, which implicity will use this
            %stream
            value = poissrnd(varargin{:});
            
            %reset default stream
            this.resetDefault();
        end
        
        function value = mnrnd(this, varargin)
            %make this stream the default random stream
            this.setDefault();
            
            %call database toolbox method, which implicity will use this
            %stream
            value = mnrnd(varargin{:});
            
            %reset default stream
            this.resetDefault();
        end
        
        function value = binornd(this, varargin)
            %make this stream the default random stream
            this.setDefault();
            
            %call database toolbox method, which implicity will use this
            %stream
            value = binornd(varargin{:});
            
            %reset default stream
            this.resetDefault();
        end
        
        %weighted random selection of integers
        function integers = randsample(this, n, k, replacement, w)
            validateattributes(n, {'numeric'}, {'nonnegative', 'integer', 'finite', 'scalar'});
            validateattributes(k, {'numeric'}, {'nonnegative', 'integer', 'finite', 'scalar'});
            validateattributes(w, {'numeric'}, {'nonnegative', 'finite'});
            
            if ~replacement
                validateattributes(k, {'numeric'}, {'<=', n});
            end
            if n == 0 || k <= 0 || ~any(w)
                integers = zeros(0, 1);
                return;
            end
            
            if ~replacement && k > 1
                if all(w == w(1))
                    integers = randsample(this.randStream, n, k, replacement);
                    return;
                end
                
                integers = find(w >= realmax);
                if numel(integers) > k
                    order = this.randperm(numel(integers));
                    integers = integers(order(1:k));
                elseif numel(integers) == k
                    return;
                end
                
                w(integers) = 0;
                tfs = false(size(w));
                tfs(integers) = true;
                nMore = k - numel(integers);
                while nMore > 0 && any(w)
                    tmpIdxs = randsample(this.randStream, n, nMore, true, w);
                    tmpTfs = false(size(tmpIdxs));
                    for i = 1:numel(tmpIdxs)
                        if ~tfs(tmpIdxs(i))
                            tfs(tmpIdxs(i)) = true;
                            tmpTfs(i) = true;
                        end
                    end
                    w(tmpIdxs) = 0;
                    integers = [integers; tmpIdxs(tmpTfs)]; %#ok<AGROW>
                    nMore = k - numel(integers);
                end
                return;
            end
            
            if k == 1
                replacement = true;
            end
            integers = randsample(this.randStream, n, k, replacement, w);
        end
        
        %returns number of times of each of numel(counts) item is selected of N
        %total selections, where item selection is weighted by counts and
        %without replacement.
        function selectedCounts = randCounts(this, counts, N)
            validateattributes(counts, {'numeric'}, {'nonnegative', 'integer', 'vector'});
            validateattributes(N, {'numeric'}, {'nonnegative', 'integer', 'scalar'});
            
            cumsumCounts = cumsum(counts);
            positiveSelect = true;
            validateattributes(N, {'numeric'}, {'<=', cumsumCounts(end)});
            if N == cumsumCounts(end)
                selectedCounts = counts;
                return;
            elseif N > cumsumCounts(end) / 2
                positiveSelect = false;
                N = cumsumCounts(end) - N;
            end
            
            selectedCounts = zeros(size(counts));
            for i = 1:N
                idx = find(randi(this.randStream, cumsumCounts(end)) <= cumsumCounts, 1, 'first');
                selectedCounts(idx) = selectedCounts(idx) + 1;
                cumsumCounts(idx:end) = cumsumCounts(idx:end) - 1;
            end
            
            if ~positiveSelect
                selectedCounts = counts - selectedCounts;
            end
        end
        
        %rounding of floats to integers weighted by decimal part
        function value = stochasticRound(this, value)
            roundUp = rand(this.randStream, size(value)) < mod(value,1);
            value(roundUp) = ceil(value(roundUp));
            value(~roundUp) = floor(value(~roundUp));
        end
        
        %randomly selects rows of a matrix mat, each with probability prob
        function [mat, rndIdxs] = randomlySelectRows(this, mat, prob)
            nRndRows = this.stochasticRound(prob * size(mat,1));
            [mat, rndIdxs] = this.randomlySelectNRows(mat, nRndRows);
        end
        
        %randomly selects n rows of a matrix mat
        function [mat, rndIdxs] = randomlySelectNRows(this, mat, nRndRows)
            rndIdxs = sort(randsample(this.randStream, size(mat, 1), min(size(mat,1), nRndRows), false));
            mat = mat(rndIdxs, :);
        end
    end
    
    %getters, setters
    methods
        function value = get.type(this)
            value = this.randStream.Type;
        end
        
        function value = get.seed(this)
            value = this.randStream.Seed;
        end
        
        function value = get.numStreams(this)
            value = this.randStream.NumStreams;
        end
        
        function value = get.streamIndex(this)
            value = this.randStream.StreamIndex;
        end
        
        function value = get.state(this)
            value = this.randStream.State;
        end
        
        function set.state(this, value)
            this.randStream.State = value;
        end
        
        function value = get.substream(this)
            value = this.randStream.Substream;
        end
        
        function value = get.randnAlg(this)
            value = this.randStream.RandnAlg;
        end
        
        function set.randnAlg(this, value)
            this.randStream.RandnAlg = value;
        end
        
        function value = get.antithetic(this)
            value = this.randStream.Antithetic;
        end
        
        function set.antithetic(this, value)
            this.randStream.Antithetic = value;
        end
        
        function value = get.fullPrecision(this)
            value = this.randStream.FullPrecision;
        end
        
        function set.fullPrecision(this, value)
            this.randStream.FullPrecision = value;
        end
    end
    
    methods
        function result = isequal(this, otherStream)
            result = this == otherStream;
        end
        
        function result = eq(this, otherStream)
            result = ...
                isequal(class(this), class(otherStream)) && ...
                this.isEqual_RandStream(this.randStream, otherStream.randStream) && ...
                ((isempty(this.defaultStream) && isempty(otherStream.defaultStream)) || this.isEqual_RandStream(this.defaultStream, otherStream.defaultStream)) && ...
                ~this.isEqual_RandStream(this.randStream, this.defaultStream);
        end
        
        function result = ne(this, otherStream)
            result = ~eq(this, otherStream);
        end
    end
    
    methods (Static = true)
        function result = isEqual_RandStream(r1, r2)
            result = ...
                isequal('RandStream', class(r1)) && ...
                isequal('RandStream', class(r2)) && ...
                isequal(r1.Type, r2.Type) && ...
                isequal(r1.Seed, r2.Seed) && ...
                isequal(r1.NumStreams, r2.NumStreams) && ...
                isequal(r1.StreamIndex, r2.StreamIndex) && ...
                isequal(r1.State, r2.State) && ...
                isequal(r1.Substream, r2.Substream) && ...
                isequal(r1.RandnAlg, r2.RandnAlg) && ...
                isequal(r1.Antithetic, r2.Antithetic) && ...
                isequal(r1.FullPrecision, r2.FullPrecision);
        end
    end
end