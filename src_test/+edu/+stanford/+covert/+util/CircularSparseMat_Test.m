% CircularSparseMat test cases
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/5/2010
classdef CircularSparseMat_Test < TestCase
    methods
        function this = CircularSparseMat_Test(name)
            this = this@TestCase(name);
        end
        
        function testSubscriptReferenceAssignment(~)
            import edu.stanford.covert.util.CircularSparseMat;

            %reference
            spmat = CircularSparseMat([2 2],1,[100 2],1);
            spmat2 = CircularSparseMat([2 2],1,[100 2],2);
            assertEqual([0; 0; 0; 1], spmat([99 2; 100 2; 101 2; 102 2]));
            assertEqual(CircularSparseMat([4 1],1,[4 1]), spmat(99:102,2));

            assertEqual([0; 0; 0; 1], spmat([99 2; 100 2; 101 2; 102 2]));
            assertEqual(CircularSparseMat([],[],[1 2]), spmat(1,:));
            assertEqual(CircularSparseMat([],[],[1 2],2), spmat2(1,:));
            assertEqual(CircularSparseMat([],[],[100 1],1), spmat(:,1));
            assertEqual(CircularSparseMat([2 1],1,[100 1],1), spmat(:,2));
            assertEqual(CircularSparseMat([],[],[100 1]), spmat(1:100,1));
            assertEqual(CircularSparseMat([2 1],1,[100 1]), spmat(1:100,2));

            assertEqual([0; 0; 0; 1], spmat([99 2; 100 2; 101 2; 102 2]));
            assertEqual(CircularSparseMat([4 1],1,[4 1]), spmat(99:102,2));

            %assignment
            spmat = CircularSparseMat([],[],[100 2],1);
            spmat([102 2])=1;
            assertEqual(CircularSparseMat([2 2],1,[100 2],1), spmat);

            spmat = CircularSparseMat([],[],[100 2],1);
            spmat(102,2)=1;
            assertEqual(CircularSparseMat([2 2],1,[100 2],1), spmat);

            spmat = CircularSparseMat([],[],[100 2],1);
            spmat(100,3)=1;
            assertEqual(CircularSparseMat([100 3],1,[100 3],1), spmat);
        end
        
        function testDeserializationFromStruct(~)
            import edu.stanford.covert.util.CircularSparseMat;

            assertEqual(...
                CircularSparseMat([2 2], 1, [100 2], 1),...
                CircularSparseMat(struct(...
                    'subs', [2 2], 'vals', 1, 'siz', [100 2], 'circularDims', 1)));
        end
    end
end
