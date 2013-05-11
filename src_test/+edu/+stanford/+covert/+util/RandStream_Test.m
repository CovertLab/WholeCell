% RandStream test cases
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/5/2010
classdef RandStream_Test < TestCase
    %properties
    properties
        randStream
    end
    
    %constructor
    methods
        function this = RandStream_Test(name)
            this = this@TestCase(name);
            
            this.randStream = edu.stanford.covert.util.RandStream('mcg16807', 'seed', 1);
        end
    end
    
    %fixture
    methods
        function setUp(this)
            this.randStream.reset(1);
        end
    end
    
    %tests
    methods
        function testRandsample(this)
            r = this.randStream;
            
            %% weights
            assertExceptionThrown(@nonVectorTest, 'MATLAB:noSuchMethodOrField');
            function nonVectorTest()
                r.randw([0 0;0 0]);
            end
            
            assertEqual(zeros(0,1), r.randsample(2, 1, true, [0 0]));
            
            assertExceptionThrown(@negativeTest, 'MATLAB:expectedNonnegative');
            function negativeTest()
                r.randsample(2, 1, true, [-1 0]);
            end
            
            assertExceptionThrown(@infTest, 'MATLAB:expectedFinite');
            function infTest()
                r.randsample(2, 1, true, [0 Inf]);
            end
            
            assertExceptionThrown(@nanTest, 'MATLAB:expectedFinite');
            function nanTest()
                r.randsample(2, 1, true, [0 NaN]);
            end
            
            %% numIntegers
            assertExceptionThrown(@negativeNTest, 'MATLAB:expectedNonnegative');
            function negativeNTest()
                r.randsample(2, -1, true, [0 1]);
            end
            
            assertExceptionThrown(@decimalNTest, 'MATLAB:expectedInteger');
            function decimalNTest()
                r.randsample(2, 1.1, true, [0 1]);
            end
            
            assertExceptionThrown(@infNTest, 'MATLAB:expectedInteger');
            function infNTest()
                r.randsample(2, Inf, true, [0 1]);
            end
            
            assertExceptionThrown(@nanNTest, 'MATLAB:expectedInteger');
            function nanNTest()
                r.randsample(2, NaN, true, [0 1]);
            end
            
            assertExceptionThrown(@matNTest, 'MATLAB:expectedScalar');
            function matNTest()
                r.randsample(2, [1 1], true, [0 1]);
            end
            
            %% sampling
            assertEqual(1, numel(r.randsample(7, 1, true, [0 1 1 1 1 0 2])));
            assertEqual(0, numel(r.randsample(7, 0, true, [0 1 1 1 1 0 2])));
            assertEqual(10, numel(r.randsample(7, 10, true, [0 1 1 1 1 0 2])));
            
            assertExceptionThrown(@replacementTest, 'MATLAB:notLessEqual');
            function replacementTest()
                r.randsample(7, 10, false, [0 1 1 1 1 0 2]);
            end
            
            assertEqual([4; 2; 7; 5; 3], r.randsample(7, 5, false, [0 1 1 1 1 0 2]));
            assertEqual([5; 4; 2; 3; 7], r.randsample(7, 5, false, [0 1 1 1 1 0 2]'));
            
            %% fairness
            counts = zeros(2, 1);
            for i = 1:200
                idx = r.randsample(2, 1, true, [1 1]);
                counts(idx) = counts(idx) + 1;
            end
            assertTrue(range(counts) < 0.1 * max(counts));
        end
        
        %tests testRandomlySelectRows and implicity testRandomlyNSelectRows
        function testRandomlySelectRows(this)
            r = this.randStream;
            
            mat = reshape(1:10,5,2);
            assertEqual(mat, r.randomlySelectRows(mat, 1));
            
            mat = reshape(1:4,2,2);
            counts = zeros(2,1);
            for i=1:1000
                [selMat, selIdx] = r.randomlySelectRows(mat, 0.5);
                assertEqual(selMat, mat(selIdx,:), 'Randomly selected rows should correspond to randomly selected indices');
                assertEqual(selIdx, unique(selIdx), 'Random selection should not be done with replacement');
                counts(selIdx)=counts(selIdx)+1;
            end
            assertTrue(max(counts)-min(counts) < 0.1 * max(counts), 'Random selection should be unbiased');
        end
        
        %tests testRandomlyNSelectRows
        function testRandomlyNSelectRows(this)
            r = this.randStream;
            
            mat = reshape(1:10,5,2);
            assertEqual(mat, r.randomlySelectNRows(mat, 10));
            assertEqual(5, size(r.randomlySelectNRows(mat, 5),1));
            
            mat = reshape(1:4,2,2);
            counts = zeros(2,1);
            for i=1:1000
                [selMat, selIdx] = r.randomlySelectNRows(mat, 1);
                assertEqual(selMat, mat(selIdx,:), 'Randomly selected rows should correspond to randomly selected indices');
                assertEqual(selIdx, unique(selIdx), 'Random selection should not be done with replacement');
                counts(selIdx)=counts(selIdx)+1;
            end
            assertTrue(max(counts)-min(counts) < 0.1 * max(counts), 'Random selection should be unbiased');
        end
        
        function testRandCounts(this)
            r = this.randStream;
            
            assertExceptionThrown(@() r.randCounts(2, -1), 'MATLAB:expectedNonnegative');
            assertEqual(0, r.randCounts(2, 0));
            assertEqual(1, r.randCounts(2, 1));
            assertEqual(2, r.randCounts(2, 2));
            
            assertEqual([2 2 3], r.randCounts([2 2 3], 7));
            assertExceptionThrown(@() r.randCounts([2 2 3], 8), 'MATLAB:notLessEqual');
            
            assertEqual(2, sum(r.randCounts([2 2 3], 2)));
            assertTrue(all([2 2 3] >= r.randCounts([2 2 3], 6)));
            
            counts = zeros(1,2);
            for i=1:250
                tmp = r.randCounts([2 2], 1);
                assertEqual(1, sum(tmp));
                counts = counts + tmp;
            end
            assertTrue(range(counts) < 0.1*max(counts));
        end
        
        function testGetters(this)
            assertEqual('mcg16807', this.randStream.type);
            assertEqual(uint32(1), this.randStream.seed);
            assertEqual(uint64(1), this.randStream.numStreams);
            assertEqual(uint64(1), this.randStream.streamIndex);
            assertEqual(1, this.randStream.state);
            assertEqual(1, this.randStream.substream);
            assertEqual('Polar', this.randStream.randnAlg);
            assertEqual(false, this.randStream.antithetic);
            assertEqual(true, this.randStream.fullPrecision);
        end
        
        function testSetters(this)
            this.randStream.state = 10;
            assertEqual(10, this.randStream.state);
            this.randStream.state = 2;
            assertEqual(2, this.randStream.state);
            
            this.randStream.randnAlg = 'Ziggurat';
            assertEqual('Ziggurat', this.randStream.randnAlg);
            this.randStream.randnAlg = 'Polar';
            assertEqual('Polar', this.randStream.randnAlg);
            
            this.randStream.antithetic = true;
            assertEqual(true, this.randStream.antithetic);
            this.randStream.antithetic = false;
            assertEqual(false, this.randStream.antithetic);
            
            this.randStream.fullPrecision = false;
            assertEqual(false, this.randStream.fullPrecision);
            this.randStream.fullPrecision = true;
            assertEqual(true, this.randStream.fullPrecision);
        end
        
        function testEquality(~)
            r1 = edu.stanford.covert.util.RandStream('mcg16807', 'seed', 1);
            r2 = edu.stanford.covert.util.RandStream('mcg16807', 'seed', 1);
            
            r1.rand(1238, 1);
            r2.rand(1238, 1);
            
            assertTrue(isequal(r1, r2));
            assertTrue(r1 == r2);
            assertFalse(r1 ~= r2);
        end
        
        function testReproducibility(this)
            %% hold state of default stream
            defaultStream = RandStream.getDefaultStream();
            defaultStreamType = defaultStream.Type;
            defaultStreamSeed = defaultStream.Seed;
            defaultStreamNumStreams = defaultStream.NumStreams;
            defaultStreamStreamIndex = defaultStream.StreamIndex;
            defaultStreamState = defaultStream.State;
            defaultStreamSubstram = defaultStream.Substream;
            defaultStreamRandnAlg = defaultStream.RandnAlg;
            defaultStreamAntithetic = defaultStream.Antithetic;
            defaultStreamFullPrecision = defaultStream.FullPrecision;
            
            %% use random streams
            [r1, val1] = this.helper_testReproducibility(112823);
            [r2, val2] = this.helper_testReproducibility(112823);
            
            %% assert that their states are equal
            assertTrue(isequal(r1, r2));
            assertTrue(r1 == r2);
            assertFalse(r1 ~= r2);
            assertEqual(val1, val2);
            
            %% assert state of default stream unchanged
            defaultStream = RandStream.getDefaultStream();
            assertEqual(defaultStreamType, defaultStream.Type);
            assertEqual(defaultStreamSeed, defaultStream.Seed);
            assertEqual(defaultStreamNumStreams, defaultStream.NumStreams);
            assertEqual(defaultStreamStreamIndex, defaultStream.StreamIndex);
            assertEqual(defaultStreamState, defaultStream.State);
            assertEqual(defaultStreamSubstram, defaultStream.Substream);
            assertEqual(defaultStreamRandnAlg, defaultStream.RandnAlg);
            assertEqual(defaultStreamAntithetic, defaultStream.Antithetic);
            assertEqual(defaultStreamFullPrecision, defaultStream.FullPrecision);
        end
    end
    
    methods
        function [stream, val] = helper_testReproducibility(~, seed)
            stream = edu.stanford.covert.util.RandStream('mcg16807', 'seed', seed);
            val.res1 = stream.rand(238, 2130);
            val.res2 = stream.randi(1239, 28, 2138);
            val.res3 = stream.randn(1230, 489);
            val.res4 = stream.randperm(178);
            val.res5 = stream.random('poisson', 17, 123, 1629);
            
            val.res6 = stream.rand(127, 1);
            probs = val.res6 / sum(val.res6);
            val.res7 = stream.mnrnd(12378, probs);
            
            val.res8 = stream.randsample(127, 61, true, ones(127, 1));
            val.res9 = stream.randsample(127, 61, false, ones(127, 1));
            val.res10 = stream.randsample(127, 61, true, probs);
            val.res11 = stream.randsample(127, 61, false, probs);
            
            val.res12 = stream.randCounts(ceil(1000 * probs), 101);
            val.res13 = stream.stochasticRound(probs);
            val.res14 = stream.randomlySelectRows(stream.rand(1238, 8), 0.20);
            val.res15 = stream.randomlySelectNRows(stream.rand(145, 8), 19);
        end
    end
end