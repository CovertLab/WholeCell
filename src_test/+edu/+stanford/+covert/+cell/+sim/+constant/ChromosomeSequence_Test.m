% ChromosomeSequence test cases
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/5/2010
classdef ChromosomeSequence_Test < TestCase
    %constructor
    methods
        function this = ChromosomeSequence_Test(name)
            this = this@TestCase(name);
        end
    end

    %tests
    methods
        function testConstructor(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = [
                'AGATCATTACGCGGGAGTTCTGCAATAGTAAAGAGATCACTCATATACGG' ...
                'ATCTCCACCTTATGGGGTGCGGATGATAGCGGATAGCGGATGTTCCTTGC' ...
                'GAAGTCGCGCAACTGTCTCTGAGTTGGCCTGATAGGGTAGCACCGCCTAT' ...
                'CCCTATTGCACAAGGTAACTTCAGTTATGAGGGCCACGTATCCCGCCAGT' ...
                'GTCGAGAACGACATGATGGGGGAACGGTTTTCTGTAACCTAGAGAACATT' ...
                'TTTGCCTAGCTAACCTCTATGGTCGATTGGCCATCTTAGGGTCCTTGGTC' ...
                'GCGTGATTTTGCGACGTGCCAATTCCTTACCGTGCCCCGTCCCCGAATAA' ...
                'GGTTAGAGTTCTCATCATTCGAAGTCGATATGGTTACAGGGCCTCCAACG' ...
                'TCCTGTACATACGGCGAAGAAACCACCGACTTGAGAGTCACTCAGCTAAT' ...
                'TCCGTTCTCGGAGCACTTCAACTGCGCGGTCACCCACGACAGATTTGGGA'];
            cseq = ChromosomeSequence(seq); %#ok<NASGU>
        end

        function testSize(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = [
                'AGATCATTACGCGGGAGTTCTGCAATAGTAAAGAGATCACTCATATACGG' ...
                'ATCTCCACCTTATGGGGTGCGGATGATAGCGGATAGCGGATGTTCCTTGC' ...
                'GAAGTCGCGCAACTGTCTCTGAGTTGGCCTGATAGGGTAGCACCGCCTAT' ...
                'CCCTATTGCACAAGGTAACTTCAGTTATGAGGGCCACGTATCCCGCCAGT' ...
                'GTCGAGAACGACATGATGGGGGAACGGTTTTCTGTAACCTAGAGAACATT' ...
                'TTTGCCTAGCTAACCTCTATGGTCGATTGGCCATCTTAGGGTCCTTGGTC' ...
                'GCGTGATTTTGCGACGTGCCAATTCCTTACCGTGCCCCGTCCCCGAATAA' ...
                'GGTTAGAGTTCTCATCATTCGAAGTCGATATGGTTACAGGGCCTCCAACG' ...
                'TCCTGTACATACGGCGAAGAAACCACCGACTTGAGAGTCACTCAGCTAAT' ...
                'TCCGTTCTCGGAGCACTTCAACTGCGCGGTCACCCACGACAGATTTGGGA'];
            cseq = ChromosomeSequence(seq);

            %isempty
            assertEqual(false, isempty(cseq));

            %numel
            assertEqual(4*numel(seq), numel(cseq));

            %length
            assertEqual(length(seq), length(cseq));

            %size
            assertEqual([numel(seq) 4], size(cseq));

            %ndims
            assertEqual(2, ndims(cseq));
        end

        function testResize(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = [
                'AGATCATTACGCGGGAGTTCTGCAATAGTAAAGAGATCACTCATATACGG' ...
                'ATCTCCACCTTATGGGGTGCGGATGATAGCGGATAGCGGATGTTCCTTGC' ...
                'GAAGTCGCGCAACTGTCTCTGAGTTGGCCTGATAGGGTAGCACCGCCTAT' ...
                'CCCTATTGCACAAGGTAACTTCAGTTATGAGGGCCACGTATCCCGCCAGT' ...
                'GTCGAGAACGACATGATGGGGGAACGGTTTTCTGTAACCTAGAGAACATT' ...
                'TTTGCCTAGCTAACCTCTATGGTCGATTGGCCATCTTAGGGTCCTTGGTC' ...
                'GCGTGATTTTGCGACGTGCCAATTCCTTACCGTGCCCCGTCCCCGAATAA' ...
                'GGTTAGAGTTCTCATCATTCGAAGTCGATATGGTTACAGGGCCTCCAACG' ...
                'TCCTGTACATACGGCGAAGAAACCACCGACTTGAGAGTCACTCAGCTAAT' ...
                'TCCGTTCTCGGAGCACTTCAACTGCGCGGTCACCCACGACAGATTTGGGA'];
            fullSeq = repmat([seq; seqcomplement(seq)],2,1)';
            cseq = ChromosomeSequence(seq);

            %permute
            assertEqual(permute(fullSeq,[2 4 1 3]), permute(cseq, [2 4 1 3]));
            assertEqual(permute(fullSeq,[2 5 4 1 3]), permute(cseq, [2 5 4 1 3]));

            %transpose
            assertEqual(fullSeq', cseq');
            assertEqual(transpose(fullSeq), transpose(cseq));

            %ctranspose
            assertEqual(ctranspose(fullSeq), ctranspose(cseq));

            %reshape
            assertEqual(reshape(fullSeq,4,5,[]), reshape(cseq,4,5,[]));

            %squeeze
            assertEqual(squeeze(fullSeq), squeeze(cseq));
        end

        function testSubscriptReferenceAssignment(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = [
                'AGATCATTACGCGGGAGTTCTGCAATAGTAAAGAGATCACTCATATACGG' ...
                'ATCTCCACCTTATGGGGTGCGGATGATAGCGGATAGCGGATGTTCCTTGC' ...
                'GAAGTCGCGCAACTGTCTCTGAGTTGGCCTGATAGGGTAGCACCGCCTAT' ...
                'CCCTATTGCACAAGGTAACTTCAGTTATGAGGGCCACGTATCCCGCCAGT' ...
                'GTCGAGAACGACATGATGGGGGAACGGTTTTCTGTAACCTAGAGAACATT' ...
                'TTTGCCTAGCTAACCTCTATGGTCGATTGGCCATCTTAGGGTCCTTGGTC' ...
                'GCGTGATTTTGCGACGTGCCAATTCCTTACCGTGCCCCGTCCCCGAATAA' ...
                'GGTTAGAGTTCTCATCATTCGAAGTCGATATGGTTACAGGGCCTCCAACG' ...
                'TCCTGTACATACGGCGAAGAAACCACCGACTTGAGAGTCACTCAGCTAAT' ...
                'TCCGTTCTCGGAGCACTTCAACTGCGCGGTCACCCACGACAGATTTGGGA'];
            fullSeq = repmat([seq; seqcomplement(seq)],2,1)';
            cseq = ChromosomeSequence(seq);

            %reference
            assertEqual(fullSeq(1:10), cseq(1:10));
            assertEqual(fullSeq(end-9:end), cseq(end-9:end));
            assertEqual(fullSeq(1:10,1), cseq(1:10,1));
            assertEqual(fullSeq(1:10,:), cseq(1:10,:));
            assertEqual(fullSeq(:,[1:4 1:2]), cseq(:,1:6));
            assertEqual(fullSeq(:,1:end), cseq(:,1:end));
            assertEqual(fullSeq(1:10,1:2), cseq(1:10,5:6));
            assertEqual(fullSeq(1:10,1:2), cseq(end+1:end+10,5:6));
            assertEqual(fullSeq([end-9:end 1:10],1:2), cseq(end-9:end+10,5:6));
            assertEqual(fullSeq([end-9:end 1:10],1:2), cseq(end-9:end+10,5:6));

            %assignment
            assertExceptionThrown(@assignmentTest1, 'ChromosomeSequence:invalidSyntax');
            function assignmentTest1()
                cseq(1)='A';
            end

            assertExceptionThrown(@assignmentTest2, 'ChromosomeSequence:invalidSyntax');
            function assignmentTest2()
                cseq(:,1)='A';
            end

            assertExceptionThrown(@assignmentTest3, 'ChromosomeSequence:invalidSyntax');
            function assignmentTest3()
                cseq(end+1:end+10,:)='A';
            end

            assertExceptionThrown(@assignmentTest4, 'ChromosomeSequence:invalidSyntax');
            function assignmentTest4()
                cseq(:)='A';
            end
        end

        function testEnd(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = [
                'AGATCATTACGCGGGAGTTCTGCAATAGTAAAGAGATCACTCATATACGG' ...
                'ATCTCCACCTTATGGGGTGCGGATGATAGCGGATAGCGGATGTTCCTTGC' ...
                'GAAGTCGCGCAACTGTCTCTGAGTTGGCCTGATAGGGTAGCACCGCCTAT' ...
                'CCCTATTGCACAAGGTAACTTCAGTTATGAGGGCCACGTATCCCGCCAGT' ...
                'GTCGAGAACGACATGATGGGGGAACGGTTTTCTGTAACCTAGAGAACATT' ...
                'TTTGCCTAGCTAACCTCTATGGTCGATTGGCCATCTTAGGGTCCTTGGTC' ...
                'GCGTGATTTTGCGACGTGCCAATTCCTTACCGTGCCCCGTCCCCGAATAA' ...
                'GGTTAGAGTTCTCATCATTCGAAGTCGATATGGTTACAGGGCCTCCAACG' ...
                'TCCTGTACATACGGCGAAGAAACCACCGACTTGAGAGTCACTCAGCTAAT' ...
                'TCCGTTCTCGGAGCACTTCAACTGCGCGGTCACCCACGACAGATTTGGGA'];
            fullSeq = repmat([seq; seqcomplement(seq)],2,1)';
            cseq = ChromosomeSequence(seq);

            assertEqual(fullSeq(end-2:end,1), cseq(end-2:end,1));
            assertEqual(fullSeq(:,[1 end-1]), cseq(:,[1 end-1]));
            assertEqual(fullSeq(1:end), cseq(1:end));
        end

        %concatenation
        function testConcatenation(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = [
                'AGATCATTACGCGGGAGTTCTGCAATAGTAAAGAGATCACTCATATACGG' ...
                'ATCTCCACCTTATGGGGTGCGGATGATAGCGGATAGCGGATGTTCCTTGC' ...
                'GAAGTCGCGCAACTGTCTCTGAGTTGGCCTGATAGGGTAGCACCGCCTAT' ...
                'CCCTATTGCACAAGGTAACTTCAGTTATGAGGGCCACGTATCCCGCCAGT' ...
                'GTCGAGAACGACATGATGGGGGAACGGTTTTCTGTAACCTAGAGAACATT' ...
                'TTTGCCTAGCTAACCTCTATGGTCGATTGGCCATCTTAGGGTCCTTGGTC' ...
                'GCGTGATTTTGCGACGTGCCAATTCCTTACCGTGCCCCGTCCCCGAATAA' ...
                'GGTTAGAGTTCTCATCATTCGAAGTCGATATGGTTACAGGGCCTCCAACG' ...
                'TCCTGTACATACGGCGAAGAAACCACCGACTTGAGAGTCACTCAGCTAAT' ...
                'TCCGTTCTCGGAGCACTTCAACTGCGCGGTCACCCACGACAGATTTGGGA'];
            fullSeq = repmat([seq; seqcomplement(seq)],2,1)';
            cseq = ChromosomeSequence(seq);

            seq2 = [
                'TCAATTGTTTTTACAGGTTAACATAGCCCATGCGACTACGAGAATATCTC' ...
                'ATTGAATAAATCGAAAACACCCGTGTACGTAACCAAGGATACTCTGCATT' ...
                'CGAGTATTATTCGGCTCTCTGGTGACAAATGGCAGGCAAGAACCTTCGGG' ...
                'AGCACGTCGGCTAAATTAAAGTCAATACGCCCTGTACCCTATCGAAACGT' ...
                'TGGAGATCTAAGAATATTAGTAAGACTGAAATAGTTCCCTGTTATAAAGG' ...
                'TGCATAGCTTGCTTCTAACGTTCCTATCGCTGTGATTATCGTTAGCAGCT' ...
                'TTTATCATCTCAATTTCGTGAGCTTGCGTTCCACCCTACGCCTAAAAAGC' ...
                'AGGAATGGACCATCGTAGACGTGCATTGCTTCCCGAGTCGCTGTCCGTTG' ...
                'CAATGTACATTTTGGATCGGGAGGCGTTCGACTGCAGGACATTTTCAGAA' ...
                'AGGATATGCAGGAAGTCAGCCCCGAAGAATCTCTTCGCTTTACGTCGTGG'];
            fullSeq2 = repmat([seq2; seqcomplement(seq2)],2,1)';
            cseq2 = ChromosomeSequence(seq2);


            assertEqual([fullSeq;fullSeq2], [cseq;cseq2]);
            assertEqual([fullSeq fullSeq2], [cseq cseq2]);
        end

        function testEquality(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = [
                'AGATCATTACGCGGGAGTTCTGCAATAGTAAAGAGATCACTCATATACGG' ...
                'ATCTCCACCTTATGGGGTGCGGATGATAGCGGATAGCGGATGTTCCTTGC' ...
                'GAAGTCGCGCAACTGTCTCTGAGTTGGCCTGATAGGGTAGCACCGCCTAT' ...
                'CCCTATTGCACAAGGTAACTTCAGTTATGAGGGCCACGTATCCCGCCAGT' ...
                'GTCGAGAACGACATGATGGGGGAACGGTTTTCTGTAACCTAGAGAACATT' ...
                'TTTGCCTAGCTAACCTCTATGGTCGATTGGCCATCTTAGGGTCCTTGGTC' ...
                'GCGTGATTTTGCGACGTGCCAATTCCTTACCGTGCCCCGTCCCCGAATAA' ...
                'GGTTAGAGTTCTCATCATTCGAAGTCGATATGGTTACAGGGCCTCCAACG' ...
                'TCCTGTACATACGGCGAAGAAACCACCGACTTGAGAGTCACTCAGCTAAT' ...
                'TCCGTTCTCGGAGCACTTCAACTGCGCGGTCACCCACGACAGATTTGGGA'];
            fullSeq = repmat([seq; seqcomplement(seq)],2,1)';
            cseq = ChromosomeSequence(seq);

            seq2 = [
                'TCAATTGTTTTTACAGGTTAACATAGCCCATGCGACTACGAGAATATCTC' ...
                'ATTGAATAAATCGAAAACACCCGTGTACGTAACCAAGGATACTCTGCATT' ...
                'CGAGTATTATTCGGCTCTCTGGTGACAAATGGCAGGCAAGAACCTTCGGG' ...
                'AGCACGTCGGCTAAATTAAAGTCAATACGCCCTGTACCCTATCGAAACGT' ...
                'TGGAGATCTAAGAATATTAGTAAGACTGAAATAGTTCCCTGTTATAAAGG' ...
                'TGCATAGCTTGCTTCTAACGTTCCTATCGCTGTGATTATCGTTAGCAGCT' ...
                'TTTATCATCTCAATTTCGTGAGCTTGCGTTCCACCCTACGCCTAAAAAGC' ...
                'AGGAATGGACCATCGTAGACGTGCATTGCTTCCCGAGTCGCTGTCCGTTG' ...
                'CAATGTACATTTTGGATCGGGAGGCGTTCGACTGCAGGACATTTTCAGAA' ...
                'AGGATATGCAGGAAGTCAGCCCCGAAGAATCTCTTCGCTTTACGTCGTGG'];
            fullSeq2 = repmat([seq2; seqcomplement(seq2)],2,1)';
            cseq2 = ChromosomeSequence(seq2);

            assertEqual(fullSeq == fullSeq2, cseq == cseq2);
            assertEqual(fullSeq == fullSeq2, fullSeq == cseq2);
            assertEqual(fullSeq == fullSeq2, cseq == fullSeq2);

            assertEqual(fullSeq ~= fullSeq2, cseq ~= cseq2);
            assertEqual(fullSeq ~= fullSeq2, fullSeq ~= cseq2);
            assertEqual(fullSeq ~= fullSeq2, cseq ~= fullSeq2);
        end

        function testFind(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = [
                'AGATCATTACGCGGGAGTTCTGCAATAGTAAAGAGATCACTCATATACGG' ...
                'ATCTCCACCTTATGGGGTGCGGATGATAGCGGATAGCGGATGTTCCTTGC' ...
                'GAAGTCGCGCAACTGTCTCTGAGTTGGCCTGATAGGGTAGCACCGCCTAT' ...
                'CCCTATTGCACAAGGTAACTTCAGTTATGAGGGCCACGTATCCCGCCAGT' ...
                'GTCGAGAACGACATGATGGGGGAACGGTTTTCTGTAACCTAGAGAACATT' ...
                'TTTGCCTAGCTAACCTCTATGGTCGATTGGCCATCTTAGGGTCCTTGGTC' ...
                'GCGTGATTTTGCGACGTGCCAATTCCTTACCGTGCCCCGTCCCCGAATAA' ...
                'GGTTAGAGTTCTCATCATTCGAAGTCGATATGGTTACAGGGCCTCCAACG' ...
                'TCCTGTACATACGGCGAAGAAACCACCGACTTGAGAGTCACTCAGCTAAT' ...
                'TCCGTTCTCGGAGCACTTCAACTGCGCGGTCACCCACGACAGATTTGGGA'];
            fullSeq = repmat([seq; seqcomplement(seq)],2,1)';
            cseq = ChromosomeSequence(seq);

            %1 output
            i1=find(fullSeq);
            i2=find(cseq);

            assertEqual(i1, i2);

            %2 outputs
            [i1,j1]=find(fullSeq);
            [i2,j2]=find(cseq);

            assertEqual(i1, i2);
            assertEqual(j1, j2);

            %3 outputs
            [i1,j1,k1]=find(fullSeq);
            [i2,j2,k2]=find(cseq);

            assertEqual(i1, i2);
            assertEqual(j1, j2);
            assertEqual(k1, k2);
        end

        function testFindSubsequence(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = 'AAAAATTAATGCTCAGTCCAAAAAA';
            cseq = ChromosomeSequence(seq);

            %example 1
            assertEqual([6 1], cseq.findSubsequence('TTANNNNNNNGTCY', 1, 1));

            %example 2
            assertEqual([9 1], cseq.findSubsequence('TTANNNNNNNGTCY',4,1));

            %example 3
            assertEqual(20, cseq.findSubsequence('AAAAAA',1,0));

            %example 4
            assertEqual(20, cseq.findSubsequence('AAAAAAAAAAATT',1,0));

            %example 5
            assertEqual([7 2], cseq.findSubsequence('AATTTTTTTTTTT',1,1));

            %example 6
            assertEqual([16 2], cseq.findSubsequence('CTGAGC',1,1));

            %example 7
            assertEqual([6 1 3 5 9 21 23 25;1 2 2 2 2 2 2 2]', cseq.findSubsequence('TT',1,1));

            %example 8
            assertEqual(zeros(0,1), cseq.findSubsequence('GG',1,0));

            %example 9
            assertEqual(zeros(0,2), cseq.findSubsequence('GGGGGG',1,1));
        end

        function testLocateFeatures(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = 'AAAAATTAATGCTCAGTCCAAAAAA';
            cseq = ChromosomeSequence(seq);

            %example 1
            features = [1 4:5 8:11 15:18 16:19 24:25; 2 1 2 1 1 1 1 2 2 2 2 1 1 1 1 2 2]';
            [allPositions, allLengths, allStrands, allSequences] = cseq.locateFeatures(features);
            assertEqual([4; 8; 16; 1; 5; 18], allPositions);
            assertEqual([1; 4; 4; 3; 1; 4], allLengths);
            assertEqual([1; 1; 1; 2; 2; 2], allStrands);
            assertEqual({'A';'AATG';'GTCC';'TTT';'T';'GACT'}, allSequences);

            %example 2
            features = [1 4:5 8:11 15:18 16:19 24:25; 1 1 2 1 1 1 1 2 2 2 2 1 1 1 1 1 1]';
            [allPositions, allLengths, allStrands, allSequences] = cseq.locateFeatures(features);
            assertEqual([4; 8; 16; 24; 5; 18], allPositions);
            assertEqual([1; 4; 4; 3; 1; 4], allLengths);
            assertEqual([1; 1; 1; 1; 2; 2], allStrands);
            assertEqual({'A';'AATG';'GTCC';'AAA';'T';'GACT'}, allSequences);
        end

        function testSubsequence(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = 'AAAAATTAATGCTCAGTCCAAAAAA';
            cseq = ChromosomeSequence(seq);

            assertEqual('AAAAAT',cseq.subsequence(1:6,1));
            assertEqual('TTTTTA',cseq.subsequence(1:6,2));
            assertEqual('AAATTA',cseq.subsequence(1:6,[1 1 1 2 2 2]));
            assertEqual('AT',cseq.subsequence([0 0],[1 2]));
            assertEqual(['AAAAAT';'TAATGC'],cseq.subsequence([1:6;7:12],1));
        end

        function testSubsequenceBaseCounts(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = 'AAAAATTAATGCTCAGTCCAAAAAA';
            cseq = ChromosomeSequence(seq);

            assertEqual([5;0;0;1],cseq.subsequenceBaseCounts(1:6,1));
            assertEqual([1;0;0;5],cseq.subsequenceBaseCounts(1:6,2));
        end

        %type casting
        function testTypeCasting(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = [
                'AGATCATTACGCGGGAGTTCTGCAATAGTAAAGAGATCACTCATATACGG' ...
                'ATCTCCACCTTATGGGGTGCGGATGATAGCGGATAGCGGATGTTCCTTGC' ...
                'GAAGTCGCGCAACTGTCTCTGAGTTGGCCTGATAGGGTAGCACCGCCTAT' ...
                'CCCTATTGCACAAGGTAACTTCAGTTATGAGGGCCACGTATCCCGCCAGT' ...
                'GTCGAGAACGACATGATGGGGGAACGGTTTTCTGTAACCTAGAGAACATT' ...
                'TTTGCCTAGCTAACCTCTATGGTCGATTGGCCATCTTAGGGTCCTTGGTC' ...
                'GCGTGATTTTGCGACGTGCCAATTCCTTACCGTGCCCCGTCCCCGAATAA' ...
                'GGTTAGAGTTCTCATCATTCGAAGTCGATATGGTTACAGGGCCTCCAACG' ...
                'TCCTGTACATACGGCGAAGAAACCACCGACTTGAGAGTCACTCAGCTAAT' ...
                'TCCGTTCTCGGAGCACTTCAACTGCGCGGTCACCCACGACAGATTTGGGA'];
            fullSeq = repmat([seq; seqcomplement(seq)],2,1)';
            cseq = ChromosomeSequence(seq);

            %cast
            assertEqual(fullSeq, cast(cseq,'char'));
            assertEqual(cast(fullSeq,'double'), cast(cseq,'double'));
            assertEqual(cast(fullSeq,'int32'), cast(cseq,'int32'));

            assertEqual(char(fullSeq), char(cseq));
            assertEqual(int8(fullSeq), int8(cseq));
            assertEqual(int16(fullSeq), int16(cseq));
            assertEqual(int32(fullSeq), int32(cseq));
            assertEqual(int64(fullSeq), int64(cseq));
            assertEqual(uint8(fullSeq), uint8(cseq));
            assertEqual(uint16(fullSeq), uint16(cseq));
            assertEqual(uint32(fullSeq), uint32(cseq));
            assertEqual(uint64(fullSeq), uint64(cseq));
            assertEqual(single(fullSeq), single(cseq));
            assertEqual(double(fullSeq), double(cseq));
        end
        
        function testIsEqual(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = [
                'AGATCATTACGCGGGAGTTCTGCAATAGTAAAGAGATCACTCATATACGG' ...
                'ATCTCCACCTTATGGGGTGCGGATGATAGCGGATAGCGGATGTTCCTTGC' ...
                'GAAGTCGCGCAACTGTCTCTGAGTTGGCCTGATAGGGTAGCACCGCCTAT' ...
                'CCCTATTGCACAAGGTAACTTCAGTTATGAGGGCCACGTATCCCGCCAGT' ...
                'GTCGAGAACGACATGATGGGGGAACGGTTTTCTGTAACCTAGAGAACATT' ...
                'TTTGCCTAGCTAACCTCTATGGTCGATTGGCCATCTTAGGGTCCTTGGTC' ...
                'GCGTGATTTTGCGACGTGCCAATTCCTTACCGTGCCCCGTCCCCGAATAA' ...
                'GGTTAGAGTTCTCATCATTCGAAGTCGATATGGTTACAGGGCCTCCAACG' ...
                'TCCTGTACATACGGCGAAGAAACCACCGACTTGAGAGTCACTCAGCTAAT' ...
                'TCCGTTCTCGGAGCACTTCAACTGCGCGGTCACCCACGACAGATTTGGGA'];
            cseq1 = ChromosomeSequence(seq);
            cseq2 = ChromosomeSequence([seq 'A']);
            
            assertTrue(isequal(cseq1, cseq1));
            assertFalse(isequal(cseq1, cseq2));
        end
        
        function testDisplay(~)
            import edu.stanford.covert.cell.sim.constant.ChromosomeSequence;

            seq = [
                'AGATCATTACGCGGGAGTTCTGCAATAGTAAAGAGATCACTCATATACGG' ...
                'ATCTCCACCTTATGGGGTGCGGATGATAGCGGATAGCGGATGTTCCTTGC' ...
                'GAAGTCGCGCAACTGTCTCTGAGTTGGCCTGATAGGGTAGCACCGCCTAT' ...
                'CCCTATTGCACAAGGTAACTTCAGTTATGAGGGCCACGTATCCCGCCAGT' ...
                'GTCGAGAACGACATGATGGGGGAACGGTTTTCTGTAACCTAGAGAACATT' ...
                'TTTGCCTAGCTAACCTCTATGGTCGATTGGCCATCTTAGGGTCCTTGGTC' ...
                'GCGTGATTTTGCGACGTGCCAATTCCTTACCGTGCCCCGTCCCCGAATAA' ...
                'GGTTAGAGTTCTCATCATTCGAAGTCGATATGGTTACAGGGCCTCCAACG' ...
                'TCCTGTACATACGGCGAAGAAACCACCGACTTGAGAGTCACTCAGCTAAT' ...
                'TCCGTTCTCGGAGCACTTCAACTGCGCGGTCACCCACGACAGATTTGGGA'];
            cseq = ChromosomeSequence(seq);
            
            disp(cseq);
            display(cseq);
        end
    end
end