% SparseMat test cases
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/5/2010
classdef SparseMat_Test < TestCase
    properties
        inefficientWarning
    end
    
    %constructor
    methods
        function this = SparseMat_Test(name)
            this = this@TestCase(name);
        end
    end
    
    methods
        function setUp(this)
            this.inefficientWarning = warning('query', 'SparseMat:inefficient');
            warning('off', 'SparseMat:inefficient');
        end
        
        function tearDown(this)
            warning(this.inefficientWarning.state, 'SparseMat:inefficient');
        end
    end
    

    %tests
    methods
        function testConstructor(~)
            import edu.stanford.covert.util.SparseMat;

            %1 argument (matrix)
            mat = zeros(4,3,2);
            mat(4, 3, 1)=1;
            mat(3, 2, 2)=2;
            mat(1, 1, 2)=3;
            mat(3, 1, 1)=4;

            spmat = SparseMat(mat);
            [subs, vals]=find(spmat);
            assertEqual([4 3 2], size(spmat));
            assertEqual([3 1 1; 4 3 1; 1 1 2; 3 2 2], subs);
            assertEqual([4;1;3;2], vals);

            spmat = SparseMat(double(mat));
            assertEqual('double', valueClass(spmat));

            spmat = SparseMat(int32(mat));
            assertEqual('int32', valueClass(spmat));

            spmat = SparseMat(sparse(double(mat(:,:,1))));
            assertEqual('double', valueClass(spmat));
            
            spmat = SparseMat([],false(0,1), [1 1]);
            assertEqual('logical', valueClass(spmat));
            
            spmat = SparseMat([],false(0,0), [1 1]);
            assertEqual('logical', valueClass(spmat));
            
            spmat = SparseMat(false(0,1));
            assertEqual('logical', valueClass(spmat));

            %1 argument (struct)
            assertEqual(spmat, SparseMat(struct(...
                'subs', spmat.subs, 'vals', spmat.vals, 'siz', spmat.siz)));

            %3 arguments
            spmat = SparseMat([3 1 2; 1 2 1; 2 2 3],[1;2;3],[4 2 3]);
            [subs, vals]=find(spmat);
            assertEqual([1 2 1; 3 1 2; 2 2 3], subs);
            assertEqual([2;1;3], vals);
            assertEqual([4 2 3], size(spmat));
        end

        function testSubscriptReference(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(10,10);
            mat(10, 1) = 1;
            spmat = SparseMat(mat);

            assertEqual(zeros(0,1), spmat(zeros(0,2)))
            assertEqual(mat(sub2ind(size(mat), 10, 1)), spmat([10 1]))
            assertEqual(mat(sub2ind(size(mat), [10 10 10]', [1 2 3]')), spmat([10 1; 10 2; 10 3]))
            assertEqual(mat(sub2ind(size(mat), [10 10 10 2]', [2 1 1 1]')), spmat([10 2; 10 1; 10 1; 2 1]))

            assertEqual(SparseMat(mat([],:)), spmat([],:))
            assertEqual(SparseMat(mat([],[])), spmat([],[]))
            assertEqual(SparseMat(mat(1:2,:)), spmat(1:2,:))
            assertEqual(SparseMat(mat(1:10,:)), spmat(1:10,:))
            assertEqual(SparseMat(mat([10 10 9 8],:)), spmat([10 10 9 8], :))
            assertEqual(SparseMat(mat([10 10 9 8], [1 2 1])), spmat([10 10 9 8], [1 2 1]))
            assertEqual(SparseMat(mat(10, :)), spmat(10,:))
            assertEqual(SparseMat(mat(:, 1)), spmat(:,1))
            assertEqual(SparseMat(mat(:, 2)), spmat(:,2))
            assertEqual(SparseMat(mat(1:2,:,1)), spmat(1:2,:,1))
            
            mat = false(10, 10);
            mat(10, 1) = 1;
            spmat = SparseMat(mat);
            assertEqual(SparseMat(mat(1:2, :, 1)), spmat(1:2, :, 1))
            
            mat = zeros(10, 10);
            mat(1:2, 1:2) = [1 2; 3 4];
            spmat = SparseMat(mat);
            assertEqual(SparseMat(mat(1:2, 1:2)), spmat(1:2, 1:2))
        end

        function testSubscriptAssignment(~)
            import edu.stanford.covert.util.SparseMat;

            warning('off','SparseMat:invalidAssignment');

            spmat=SparseMat([],[],[10 10]);
            spmat([])=1;
            assertEqual(SparseMat([],[],[10 10]), spmat)

            spmat=SparseMat([],[],[10 10]);
            spmat(zeros(0,2))=1;
            assertEqual(SparseMat([],[],[10 10]), spmat)

            spmat=SparseMat([],[],[10 10]);
            spmat([1 1])=1;
            assertEqual(SparseMat([1 1],1,[10 10]), spmat)

            spmat=SparseMat([],[],[10 10]);
            spmat([1 1; 1 1])=1;
            assertEqual(SparseMat([1 1],1,[10 10]), spmat)

            spmat=SparseMat([],[],[10 10]);
            spmat([1 1; 1 1])=[1;2];
            assertEqual(SparseMat([1 1],2,[10 10]), spmat)

            spmat=SparseMat([],[],[10 10]);
            spmat(1,1)=1;
            assertEqual(SparseMat([1 1],1,[10 10]), spmat)

            spmat=SparseMat([],[],[10 10]);
            spmat(1:10,1:2)=[(1:10)' zeros(10,1)];
            assertEqual(SparseMat([(1:10)' ones(10,1)],(1:10)',[10 10]), spmat)

            spmat=SparseMat([],[],[10 10]);
            spmat(2:8,1:2)=SparseMat([2 2],3,[7 2]);
            assertEqual(SparseMat([3 2],3,[10 10]), spmat)

            spmat=SparseMat([],[],[10 10]);
            spmat([2 2],[1 1])=[1 2; 3 4];
            assertEqual(SparseMat([2 1],4,[10 10]), spmat)

            spmat=SparseMat([],[],[10 10]);
            spmat([],1)=1;
            assertEqual(SparseMat([],[],[10 10]), spmat)

            spmat=SparseMat([],[],[10 10]);
            spmat([],[])=1;
            assertEqual(SparseMat([],[],[10 10]), spmat)

            spmat=SparseMat([],[],[10 2]);
            exception = [];
            try
                spmat(1:2,:)=ones(2,3);
            catch exception
            end
            if isempty(exception)
                throw(MException('SparseMat_Test:error', 'Syntax should not be allowed'));
            end
            assertEqual(SparseMat([],[],[10 2]), spmat)
            
            spmat=SparseMat([],[],[10 10]);
            spmat(3:5,3) = [6;7;8];
            assertEqual(SparseMat([3 3; 4 3; 5 3],[6;7;8],[10 10]), spmat)
            
            spmat=SparseMat([],[],[10 10]);
            spmat(1:2,1:2) = 1;
            assertEqual(SparseMat([1 1; 2 1; 1 2; 2 2],[1;1;1;1], [10 10]), spmat)
            
            spmat=SparseMat([],[],[10 10]);
            spmat(1:2,1:2) = [1 2;3 4];
            assertEqual(SparseMat([1 1; 2 1; 1 2; 2 2],[1;3;2;4], [10 10]), spmat)
        end

        function testSubscriptReferenceAssignment(~)
            import edu.stanford.covert.util.SparseMat;

            %example 1
            st_x=SparseMat([],[],[10 10 1]);

            st_x(10,5,1)=1;
            st_x(5,7,1)=3;
            st_x(8,6,1)=4;

            x=st_x(:,:,:);

            assertEqual(ones(10,10),double(eq(st_x, st_x(:,:,:))))
            assertEqual(ones(10,10),double(eq(st_x, st_x(:,:,1))))

            %example 2
            fl_x=zeros([10 10 1]);
            st_x=SparseMat([],[],[10 10 1]);
            assertEqual(size(fl_x(1,:,:)), size(st_x(1,:,:)))
            assertEqual(size(fl_x(:,1,:)), size(st_x(:,1,:)))

            %example 3
            fl_x=zeros([10 10 1]);
            st_x=SparseMat([],[],[10 10 1]);

            fl_x(3,2)=1;
            st_x(3,2)=1;

            assertEqual(size(fl_x(3,:,:)), size(st_x(3,:,:)))
            assertEqual(size(fl_x(:,2,:)), size(st_x(:,2,:)))

            %example 4
            st_x=SparseMat([],[],[500 1 1]);
            exception = [];
            try
                st_x(1:50);
            catch exception
            end
            if isempty(exception)
                throw(MException('SparseMat_Test:error', 'Linear indexing should not be valid'));
            end

            exception = [];
            try
                st_x((1:50)');
            catch exception
            end
            if isempty(exception)
                throw(MException('SparseMat_Test:error', 'Linear indexing should not be valid'));
            end

            %example 5
            st_x=SparseMat([],[],[500 1 1]);
            exception = [];
            try
                st_x(4,5,6);
            catch exception
            end
            if isempty(exception)
                throw(MException('SparseMat_Test:error', 'Should not be able to reference unitialized dimensions'));
            end

            exception = [];
            try
                st_x([4 5 6]);
            catch exception
            end
            if isempty(exception)
                throw(MException('SparseMat_Test:error', 'Should not be able to reference unitialized dimensions'));
            end

            %example 6
            st_x=SparseMat([],[],[500 2 1]);
            assertEqual([1 1],size(st_x(4,2,1,1)));
            assertEqual([1 1],size(st_x(4,2,1,:,1)));
            assertEqual([1 2],size(st_x(4,:,1,:,1)));
            assertEqual([500 2],size(st_x(:,:,1,:,1)));

            %example 7
            st_x=SparseMat([],[],[500 2 1]);
            fl_x=st_x([(1:50)' ones(50,1)]);
            assertEqual('edu.stanford.covert.util.SparseMat',class(st_x(1:50,1)));
            assertEqual('double',class(fl_x));
            assertEqual([50 1], size(fl_x));

            %example 8
            st_x=SparseMat([],[],[500 2 1]);
            fl_x=zeros(500,2);
            assertEqual(size(fl_x([],:)), size(st_x([],:)));
            assertEqual(size(fl_x([],2)), size(st_x([],2)));
            assertEqual([0 1], size(st_x(zeros(0,2))));

            %example 9
            st_x=SparseMat([],[],[3 2]);
            fl_x=zeros(3,2);

            st_x(1,1)=1;
            fl_x(1,1)=1;

            assertEqual(fl_x([1;1],1), double(st_x([1;1],1)));
            assertEqual(fl_x([1;1],:), double(st_x([1;1],:)));
            assertEqual(fl_x([1;1;2],:), double(st_x([1;1;2],:)));
            assertEqual(fl_x([1;2;1],:), double(st_x([1;2;1],:)));
            assertEqual(fl_x([1;2;1],2), double(st_x([1;2;1],2)));
            assertEqual(fl_x([1;2;2;1;3],2), double(st_x([1;2;3;1;3],2)));
            assertEqual(fl_x([1;1;2;2],:), double(st_x([1;1;2;2],:)));
            assertEqual(fl_x(sub2ind([3 2],[1;1],[1;1])), double(st_x([1 1; 1 1])));
            assertEqual(fl_x(sub2ind([3 2],[1;1;2],[1;1;1])), double(st_x([1 1; 1 1; 2 1])));

            st_x(2,2)=3;
            fl_x(2,2)=3;

            assertEqual(fl_x([1;1],1), double(st_x([1;1],1)));
            assertEqual(fl_x([1;1],:), double(st_x([1;1],:)));
            assertEqual(fl_x([1;1;2],:), double(st_x([1;1;2],:)));
            assertEqual(fl_x([1;2;1],:), double(st_x([1;2;1],:)));
            assertEqual(fl_x([1;2;1],2), double(st_x([1;2;1],2)));
            assertEqual(fl_x([1;2;2;1;3],2), double(st_x([1;2;2;1;3],2)));
            assertEqual(fl_x([1;1;2;2],:), double(st_x([1;1;2;2],:)));
            assertEqual(fl_x([1;2;2;1;3],:), double(st_x([1;2;2;1;3],:)));
            assertEqual(fl_x(sub2ind([3 2],[1;1],[1;1])), double(st_x([1 1; 1 1])));
            assertEqual(fl_x(sub2ind([3 2],[1;1;2],[1;1;1])), double(st_x([1 1; 1 1; 2 1])));
        end

        function testSubscriptAssignmentReference(~)
            import edu.stanford.covert.util.SparseMat;

            warning('off','SparseMat:invalidAssignment');

            %example 1
            st_x=SparseMat([],[],[10 10 1]);
            st_x(8,6,:)=4;
            assertEqual([10 10],size(st_x))

            %example 2
            st_x=SparseMat([],[],[500 1 1]);
            exception = [];
            try
                st_x(1:50)=1;
            catch exception
            end
            if isempty(exception)
                throw(MException('SparseMat_Test:error', 'Invalid syntax'));
            end

            %example 3
            st_x=SparseMat([],[],[500 1 1]);
            exception = [];
            try
                st_x((1:50)')=1;
            catch exception
            end
            if isempty(exception)
                throw(MException('SparseMat_Test:error', 'Linear indexing should not be valid'));
            end

            %example 4
            st_x=SparseMat([],[],[100 1 1]);
            st_y=SparseMat(ones(10,1));
            exception = [];
            try
                st_x(1:10)=st_y;
            catch exception
            end
            if isempty(exception)
                throw(mexception('SparseMat_test:error', 'linear indexing should not be valid'));
            end

            st_z=SparseMat([],[],[100 1 1]);
            exception = [];
            try
                st_z(1:10)=1;
            catch exception
            end
            if isempty(exception)
                throw(MException('SparseMat_Test:error', 'Invalid syntax'));
            end

            %example 5
            st_x=SparseMat([],[],[100 1 1]);
            st_y=SparseMat(ones(10,1));
            exception = [];
            try
                st_x((1:10)')=st_y;
            catch exception
            end
            if isempty(exception)
                throw(mexception('SparseMat_test:error', 'linear indexing should not be valid'));
            end

            %example 6
            st_x=SparseMat([],[],[100 1 1]);

            st_x([(1:10)' ones(10,1)])=1;
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x([(1:10)' ones(10,1)])=SparseMat(ones(10,1));
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x([(1:10)' ones(10,1)])=ones(10,1);
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x([(1:10)' ones(10,1)])=1;
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x([(1:10)' ones(10,1)])=sparse(ones(10,1));
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x(1:10,1)=SparseMat(ones(10,1));
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x(1:10,1)=ones(10,1);
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x(1:10,1)=1;
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x(1:10,1)=sparse(ones(10,1));
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            %example 7
            st_x=SparseMat([],[],[100 1 1]);
            st_x([(1:10)' ones(10,1)])=SparseMat(ones(10,1));
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x=SparseMat([],[],[100 1 1]);
            st_x([(1:10)' ones(10,1)])=ones(10,1);
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x=SparseMat([],[],[100 1 1]);
            st_x([(1:10)' ones(10,1)])=1;
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x=SparseMat([],[],[100 1 1]);
            st_x([(1:10)' ones(10,1)])=sparse(ones(10,1));
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x=SparseMat([],[],[100 1 1]);
            st_x(1:10,1)=SparseMat(ones(10,1));
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x=SparseMat([],[],[100 1 1]);
            st_x(1:10,1)=ones(10,1);
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x=SparseMat([],[],[100 1 1]);
            st_x(1:10,1)=1;
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            st_x=SparseMat([],[],[100 1 1]);
            st_x(1:10,1)=sparse(ones(10,1));
            assertEqual(ones(10,1), double(st_x(1:10,1)));
            assertEqual([100 1],size(st_x));

            %example 8
            st_x=SparseMat([],[],[100 1 1]);

            st_x([(1:4)' ones(4,1)])=SparseMat([1 0 1 0]');
            assertEqual([1 0 1 0]', double(st_x(1:4,1)));
            assertEqual([100 1],size(st_x));

            st_x([(1:4)' ones(4,1)])=[1 0 1 0]';
            assertEqual([1 0 1 0]', double(st_x(1:4,1)));
            assertEqual([100 1],size(st_x));

            st_x([(1:4)' ones(4,1)])=sparse([1 0 1 0]');
            assertEqual([1 0 1 0]', double(st_x(1:4,1)));
            assertEqual([100 1],size(st_x));

            st_x(1:4,1)=SparseMat([1 0 1 0]');
            assertEqual([1 0 1 0]', double(st_x(1:4,1)));
            assertEqual([100 1],size(st_x));

            st_x(1:4,1)=[1 0 1 0]';
            assertEqual([1 0 1 0]', double(st_x(1:4,1)));
            assertEqual([100 1],size(st_x));

            st_x(1:4,1)=sparse([1 0 1 0]');
            assertEqual([1 0 1 0]', double(st_x(1:4,1)));
            assertEqual([100 1],size(st_x));

            %example 9
            st_x=SparseMat([],[],[100 1 1]);
            st_x([(1:4)' ones(4,1)])=SparseMat([1 0 1 0]');
            assertEqual([1 0 1 0]', double(st_x(1:4,1)));
            assertEqual([100 1],size(st_x));

            st_x=SparseMat([],[],[100 1 1]);
            st_x([(1:4)' ones(4,1)])=[1 0 1 0]';
            assertEqual([1 0 1 0]', double(st_x(1:4,1)));
            assertEqual([100 1],size(st_x));

            st_x=SparseMat([],[],[100 1 1]);
            st_x([(1:4)' ones(4,1)])=sparse([1 0 1 0]');
            assertEqual([1 0 1 0]', double(st_x(1:4,1)));
            assertEqual([100 1],size(st_x));

            st_x=SparseMat([],[],[100 1 1]);
            st_x(1:4,1)=SparseMat([1 0 1 0]');
            assertEqual([1 0 1 0]', double(st_x(1:4,1)));
            assertEqual([100 1],size(st_x));

            st_x=SparseMat([],[],[100 1 1]);
            st_x(1:4,1)=[1 0 1 0]';
            assertEqual([1 0 1 0]', double(st_x(1:4,1)));
            assertEqual([100 1],size(st_x));

            st_x=SparseMat([],[],[100 1 1]);
            st_x(1:4,1)=sparse([1 0 1 0]');
            assertEqual([1 0 1 0]', double(st_x(1:4,1)));
            assertEqual([100 1],size(st_x));

            %example 10
            st_x=SparseMat([],[],[100 2]);
            st_x(:,:)=1;
            assertEqual(ones(100,2), double(st_x));

            st_x=SparseMat([],[],[4 2]);
            st_x(:,:)=1;
            st_x(1:2,1)=SparseMat([],[],[2 1]);
            assertEqual([zeros(2,1) ones(2,1);ones(2,2)], double(st_x));

            st_x=SparseMat([],[],[4 2]);
            st_x(:,:)=1;
            st_x(1:2,:)=0;
            assertEqual([zeros(2,2);ones(2,2)], double(st_x));

            st_x=SparseMat([],[],[4 2]);
            st_x(:,:)=1;
            st_x(1:2,1)=0;
            assertEqual([zeros(2,1) ones(2,1);ones(2,2)], double(st_x));

            st_x=SparseMat([],[],[100 2]);
            st_x(:,:)=1;
            st_x(:,:)=0;
            assertEqual(SparseMat([],[],[100 2]), st_x);

            %example 11
            st_x=SparseMat([],[],[100 2]);
            st_y=SparseMat([],[],[100 2]);
            st_x(1:10,:)=1;
            st_y(1:10,:)=1;
            st_x([],1)=2;
            assertEqual(st_y, st_x);

            %example 12
            st_x=SparseMat([],[],[100 2]);
            st_y=SparseMat([],[],[100 2]);
            st_x(1:10,:)=1;
            st_y(1:10,:)=1;
            st_x([],:)=2;
            assertEqual(st_y, st_x);

            %example 13
            st_x=SparseMat([],[],[100 2]);
            st_y=SparseMat([],[],[100 2]);
            st_x(1:10,:)=1;
            st_y(1:10,:)=1;
            st_x([],:)=zeros(0,2);
            assertEqual(st_y, st_x);

            %example 14
            st_x=SparseMat([],[],[100 2]);
            st_y=SparseMat([],[],[100 2]);
            st_x(1:10,:)=1;
            st_y(1:10,:)=1;
            st_x([],:)=SparseMat;
            assertEqual(st_y, st_x);

            %example 15
            st_x=SparseMat([],[],[100 2]);
            st_y=SparseMat([],[],[100 2]);
            st_x(1:10,:)=1;
            st_y(1:10,:)=1;

            exception = [];
            try
                st_x(1,:)=SparseMat;
            catch exception
            end
            if isempty(exception)
                throw(MException('SparseMat_test:error', 'RHS cannot be null with non-null subtensor'));
            end
            assertEqual(st_y, st_x);

            %example 16
            st_x=SparseMat([],[],[100 2]);
            st_y=SparseMat([],[],[100 2]);
            st_x(1:10,:)=1;
            st_y(1:10,:)=1;
            st_x(zeros(0,2))=2;
            assertEqual(st_y, st_x);

            %example 17
            st_x=SparseMat([],[],[100 2]);
            st_y=SparseMat([],[],[100 2]);
            st_x(1:10,:)=1;
            st_y(1:10,:)=1;
            st_x(zeros(0,2))=zeros(0,1);
            assertEqual(st_y, st_x);

            %example 18
            st_x=SparseMat([],[],[100 2]);
            st_y=SparseMat([],[],[100 2]);
            st_x(1:10,:)=1;
            st_y(1:10,:)=1;

            exception=[];
            try
                st_x(ones(1,2))=[];
            catch exception
            end
            if isempty(exception)
                throw(mexception('SparseMat_test:error', 'RHS cannot be null with non-null subscripts'));
            end
            assertEqual(st_y, st_x);
        end

        function testConsistentBehavior(~)
            import edu.stanford.covert.util.SparseMat;

            %example 1
            fl_x=zeros([2 2 1]);
            st_x=SparseMat([],[],[2 2 1]);

            assertEqual(size(fl_x), size(st_x));
            assertEqual(size(fl_x,5),size(st_x,5));

            %example 2
            sp_x=sparse(zeros([10 10]));
            st_x=SparseMat([],[],[10 10]);

            sp_x(10,5)=1;
            sp_x(5,7)=3;
            sp_x(12,6)=4;

            st_x(10,5)=1;
            st_x(5,7)=3;
            st_x(12,6)=4;

            [i,j]=find(sp_x);
            assertEqual([i j], find(st_x));

            %example 3
            fl_x=zeros([10 10 10]);
            st_x=SparseMat([],[],[10 10 10]);

            fl_x(10,5,3)=1;
            fl_x(5,7,2)=3;
            fl_x(12,6,1)=4;
            fl_x(10,4,3)=1;
            fl_x(5,8,2)=3;
            fl_x(12,5,1)=4;

            st_x(10,5,3)=1;
            st_x(5,7,2)=3;
            st_x(12,6,1)=4;
            st_x(10,4,3)=1;
            st_x(5,8,2)=3;
            st_x(12,5,1)=4;

            [i,j,k]=ind2sub(size(fl_x),find(fl_x));
            assertEqual([i j k], find(st_x));

            [~,vals]=find(st_x);
            assertEqual(fl_x(find(fl_x)),vals);
        end

        function testInternalConsistency(~)
            import edu.stanford.covert.util.SparseMat;

            %example 1
            st_x=SparseMat([],[],[100 1 1]);
            st_y=SparseMat(ones(10,1));
            st_x(1:10,1,1)=st_y;

            subs=find(st_x);
            assertEqual(length(size(st_x)), size(subs,2));

            %example 2
            st_x=SparseMat([],[],[100 1 1]);
            st_y=SparseMat(ones(10,1));
            st_x(1:10,1,1,1)=st_y;

            subs=find(st_x);
            assertEqual(length(size(st_x)), size(subs,2));

            %example 3
            st_x=SparseMat([],[],[500 1 1]);
            st_x(1:3,1,1)=SparseMat(ones(3,1));

            st_y=SparseMat([],[],[3 1 1]);
            st_y(1,1,1)=3;
            st_y(2,1,1)=4;
            st_y(3,1,1)=5;
            st_x(4:6,1,1)=st_y;

            subs=find(st_x);
            assertEqual(length(size(st_x)), size(subs,2));
        end

        function testEnd(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            assertEqual(SparseMat(mat(2:end,:,:,:)), spmat(2:end,:,:,:));
            assertEqual(SparseMat(mat(:,:,end,:)), spmat(:,:,end,:));
            assertEqual(SparseMat(mat(:,:,:,end-1)), spmat(:,:,:,end-1));

            exception = [];
            try
                assertEqual(SparseMat(mat(1:end)), spmat(1:end));
            catch exception
            end
            if isempty(exception)
                throw(MException('SparseMat:unsupportedSubsref','SparseMat doesn''t support linear indexing'));
            end
        end

        function testSub2ind(~)
            import edu.stanford.covert.util.SparseMat;

            inds=sub2ind([100 2 3], ...
                reshape(repmat((1:100)',[1 2 3]),[],1),...
                reshape(repmat(1:2,[100 1 3]),[],1),...
                reshape(repmat(permute(1:3,[1 3 2]),[100 2 1]),[],1));

            spmat = SparseMat([],[],[100 2 3]);
            subs = [...
                reshape(repmat((1:100)',[1 2 3]),[],1) ...
                reshape(repmat(1:2,[100 1 3]),[],1) ...
                reshape(repmat(permute(1:3,[1 3 2]),[100 2 1]),[],1)];
            assertEqual(inds, sub2ind(spmat,subs));
        end

        function testInd2sub(~)
            import edu.stanford.covert.util.SparseMat;

            [i,j,k]=ind2sub([100 2 3], (1:100*2*3)');

            spmat = SparseMat([],[],[100 2 3]);
            assertEqual([i j k], ind2sub(spmat,(1:100*2*3)'));
        end

        function testSize(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            %isempty
            assertEqual(false, isempty(spmat));
            assertEqual(true, isempty(SparseMat));

            %length
            assertEqual(4, length(spmat));

            %size
            assertEqual([4 1 3 2], size(spmat));
            assertEqual(4, size(spmat,1));
            assertEqual(1, size(spmat,5));

            %numel
            assertEqual(4*1*3*2, numel(spmat));

            %ndims
            assertEqual(4, ndims(spmat));

            %nnz
            assertEqual(7, nnz(spmat));

            %permute
            assertEqual(SparseMat(permute(mat,[2 4 1 3])), permute(spmat, [2 4 1 3]));
            assertEqual(SparseMat(permute(mat,[2 5 4 1 3])), permute(spmat, [2 5 4 1 3]));

            %transpose
            mat2 = permute(sum(mat,4),[1 3 2]);
            assertEqual(SparseMat(mat2'), transpose(SparseMat(mat2)));
            assertEqual(SparseMat(mat2'), SparseMat(mat2)');

            %ctranspose
            mat2 = permute(sum(mat,4),[1 3 2])+sqrt(-1);
            assertEqual(SparseMat(ctranspose(mat2)), ctranspose(SparseMat(mat2)));

            %reshape
            assertEqual(SparseMat(reshape(mat,[3 2 4])), reshape(SparseMat(mat),[3 2 4]));
            assertEqual(SparseMat(reshape(zeros(0,2),[0 3])), reshape(SparseMat(zeros(0,2)),[0 3]));

            %repmat
            assertEqual(spmat, repmat(spmat, 1, 1));
            assertEqual(spmat, repmat(spmat, [1 1 1]));
            assertEqual(spmat, repmat(spmat, [1 1 1 1]));
            assertEqual([spmat;spmat], repmat(spmat, 2, 1));
            assertEqual(SparseMat(repmat(mat,[2 3 4 1 4 3])), repmat(spmat, [2 3 4 1 4 3]));

            %squeeze
            assertEqual(SparseMat(squeeze(mat)), squeeze(spmat));
        end

        function testRepeatMatrix(~)
            A = edu.stanford.covert.util.SparseMat(1);
            assertEqual([1 3], size(repmat(A, 1, 3)));
            assertEqual([1 3], size(repmat(A, [1 3])));
            assertEqual([0 2], size(repmat(A, [0 2])));
            assertEqual([1 1], size(repmat(A, [1 1 1])));
            assertEqual([0 1], size(repmat(A, [0 1 1])));
            assertEqual([1 0], size(repmat(A, [1 0 1])));
            assertEqual([1 1 0], size(repmat(A, [1 1 0])));
        end
        
        function testFind(~)
            import edu.stanford.covert.util.SparseMat;

            %example 1
            mat = zeros(4,3);
            mat(4, 3)=1;
            mat(3, 2)=2;
            mat(1, 1)=3;
            mat(3, 1)=4;
            spmat = SparseMat(mat);

            [i,j,k]=find(mat);
            [subs,vals]=find(spmat);

            assertEqual([i j], subs);
            assertEqual(k, vals);

            %example 2
            mat = zeros(4,3,'int32');
            mat(4, 3)=1;
            mat(3, 2)=2;
            mat(1, 1)=3;
            mat(3, 1)=4;
            spmat = SparseMat(mat);

            [i,j,k]=find(mat);
            [subs,vals]=find(spmat);

            assertEqual([i j], subs);
            assertEqual(k, vals);
        end

        function testRandomlySelectZeros(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,3);
            mat(4, 3)=1;
            mat(3, 2)=2;
            mat(1, 1)=3;
            mat(3, 1)=4;
            spmat = SparseMat(mat);

            %example 1
            randomlySelectZeros(spmat,0.2);
            assertEqual([2 2], size(randomlySelectZeros(spmat,2)));

            %example 2
            randStream = edu.stanford.covert.util.RandStream('mcg16807');
            randStream.reset(1);
            randomlySelectZeros(spmat, 0.2, randStream);
            assertEqual([2 2], size(randomlySelectZeros(spmat, 2, randStream)));
        end

        %concatenation
        function testConcatenation(~)
            import edu.stanford.covert.util.SparseMat;

            %Vertcat_nullIdentity
            assertEqual(...
                SparseMat,...
                vertcat(SparseMat));

            %Vertcat_identity
            t = SparseMat([1 2; 3 3], [2;1], [3 4]);
            assertEqual(t, vertcat(t));

            %Vertcat_threeTensors
            u = SparseMat([1 2; 3 3], [2;1], [3 4]);
            v = SparseMat([6 3; 4 2; 1 4], [3;4;-1], [7 4]);

            assertEqual([...
                0 2 0 0;
                0 0 0 0;
                0 0 1 0;
                0 0 0 -1;
                0 0 0 0;
                0 0 0 0;
                0 4 0 0;
                0 0 0 0;
                0 0 3 0;
                0 0 0 0;
                0 2 0 0;
                0 0 0 0;
                0 0 1 0],...
                double(vertcat(u, v, u)));

            assertEqual([...
                0 2 0 0;
                0 0 0 0;
                0 0 1 0;
                0 0 0 -1;
                0 0 0 0;
                0 0 0 0;
                0 4 0 0;
                0 0 0 0;
                0 0 3 0;
                0 0 0 0;
                0 2 0 0;
                0 0 0 0;
                0 0 1 0],...
                double([u; v; u]));

            %Vertcat_nullFirst
            t = SparseMat([1 2; 3 3], [2;1], [3 4]);
            assertEqual(...
                SparseMat([5 2; 7 3], [2;1], [7 4]),...
                vertcat(SparseMat([],[],[4 4]), t));

            %Vertcat_nullLast
            t = SparseMat([1 2; 3 3], [2;1], [3 4]);
            assertEqual(...
                SparseMat([1 2; 3 3], [2;1], [7 4]),...
                vertcat(t, SparseMat([],[],[4 4])));

            %Vertcat_zeros
            assertEqual(...
                SparseMat([],[],[7 2]),...
                vertcat(...
                SparseMat([],[],[3 2]), SparseMat([],[],[4 2])));

            %Vertcat_inconsistentDimensions
            assertExceptionThrown(...
                @()vertcat(...
                SparseMat([1 2; 3 3], [2;1], [3 3]),...
                SparseMat([1 2], 3, [3 2])),...
                'SparseMat:cat');

            %Vertcat_threeDimensions
            assertEqual(...
                SparseMat([1 2 1; 3 3 1; 4 1 2], [2;1;3], [4 3 2]),...
                vertcat(...
                SparseMat([1 2 1; 3 3 1], [2;1], [3 3 2]),...
                SparseMat([1 1 2], 3, [1 3 2])));

            %Vertcat_threeInconsistentDimensions
            assertExceptionThrown(...
                @()vertcat(...
                SparseMat([1 2 1; 3 3 1], [2;1], [3 3 1]),...
                SparseMat([1 1 2], 3, [1 3 2])),...
                'SparseMat:cat');

            %Horzcat_nullIdentity
            assertEqual(...
                SparseMat,...
                horzcat(SparseMat));

            %Horzcat_identity
            t = SparseMat([1 2; 3 3], [2;1], [3 4]);
            assertEqual(t, horzcat(t));

            %Horzcat_threeTensors
            u = SparseMat([2 1; 3 3], [2;1], [4 3]);
            v = SparseMat([3 6; 2 4; 4 1], [3;4;-1], [4 7]);

            assertEqual([
                0  0  0  0  0  0  0  0  0  0  0  0  0;
                2  0  0  0  0  0  4  0  0  0  2  0  0;
                0  0  1  0  0  0  0  0  3  0  0  0  1;
                0  0  0 -1  0  0  0  0  0  0  0  0  0],...
                double(horzcat(u, v, u)));

            assertEqual([
                0  0  0  0  0  0  0  0  0  0  0  0  0;
                2  0  0  0  0  0  4  0  0  0  2  0  0;
                0  0  1  0  0  0  0  0  3  0  0  0  1;
                0  0  0 -1  0  0  0  0  0  0  0  0  0],...
                double([u, v, u]));

            %Horzcat_nullFirst
            t = SparseMat([1 2; 3 3], [2;1], [3 4]);
            assertEqual(...
                SparseMat([1 7; 3 8], [2;1], [3 9]),...
                horzcat(SparseMat([],[],[3 5]), t));

            %Horzcat_nullLast
            t = SparseMat([1 2; 3 3], [2;1], [3 4]);
            assertEqual(...
                SparseMat([1 2; 3 3], [2;1], [3 9]),...
                horzcat(t, SparseMat([],[],[3 5])));

            %Horzcat_zeros
            assertEqual(...
                SparseMat([],[],[2 7]),...
                horzcat(...
                SparseMat([],[],[2 3]), SparseMat([],[],[2 4])));

            %Horzcat_inconsistentDimensions
            assertExceptionThrown(...
                @()horzcat(...
                SparseMat([1 2; 3 3], [2;1], [3 3]),...
                SparseMat([1 2], 3, [2 3])),...
                'SparseMat:cat');

            %Horzcat_threeDimensions
            assertEqual(...
                full(SparseMat([1 2 1; 3 3 1; 2 4 1], [2;1;3], [3 4 2])),...
                full(horzcat(...
                SparseMat([1 2 1; 3 3 1], [2;1], [3 3 2]),...
                SparseMat([2 1 1], 3, [3 1 2]))));

            %Horzcat_threeInconsistentDimensions
            assertExceptionThrown(...
                @()horzcat(...
                SparseMat([1 2 1; 3 3 1], [2;1], [3 3 2]),...
                SparseMat([2 1 1], 3, [3 1 1])),...
                'SparseMat:cat');

            %padarray
            x = reshape(1:15,5,3);
            spmat = SparseMat(x);

            assertEqual(SparseMat(padarray(x,3)),padarray(spmat,3));
            assertEqual(SparseMat(padarray(x,[0 2])),padarray(spmat,[0 2]));
            assertEqual(SparseMat(padarray(x,[1 0 2])),padarray(spmat,[1 0 2]));
            assertEqual(SparseMat(padarray(x,[1 0 0])),padarray(spmat,[1 0 0]));
            assertEqual(SparseMat(padarray(x,[0 1 1],1)),padarray(spmat,[0 1 1],1));
            assertEqual(SparseMat(padarray(x,[0 1 1],2,'pre')),padarray(spmat,[0 1 1],2,'pre'));
            assertEqual(SparseMat(padarray(x,[0 1 1],2,'post')),padarray(spmat,[0 1 1],2,'post'));
            assertEqual(SparseMat(padarray(x,[0 1 1],2,'both')),padarray(spmat,[0 1 1],2,'both'));

            assertEqual(SparseMat(padarray(x,[2 2],'circular','pre')),padarray(spmat,[2 2],'circular','pre'));
            assertEqual(SparseMat(padarray(x,[2 2],'replicate','pre')),padarray(spmat,[2 2],'replicate','pre'));
            assertEqual(SparseMat(padarray(x,[2 2],'symmetric','pre')),padarray(spmat,[2 2],'symmetric','pre'));
            assertEqual(SparseMat(padarray(x,[2 2],'circular','post')),padarray(spmat,[2 2],'circular','post'));
            assertEqual(SparseMat(padarray(x,[2 2],'replicate','post')),padarray(spmat,[2 2],'replicate','post'));
            assertEqual(SparseMat(padarray(x,[2 2],'symmetric','post')),padarray(spmat,[2 2],'symmetric','post'));
            assertEqual(SparseMat(padarray(x,[2 2],'circular','both')),padarray(spmat,[2 2],'circular','both'));
            assertEqual(SparseMat(padarray(x,[2 2],'replicate','both')),padarray(spmat,[2 2],'replicate','both'));
            assertEqual(SparseMat(padarray(x,[2 2],'symmetric','both')),padarray(spmat,[2 2],'symmetric','both'));
        end

        function testAndOrXorNot(test)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            %mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            spmat2 = SparseMat(mat2);

            %not
            assertEqual(SparseMat(~mat), ~spmat);
            assertEqual(SparseMat(~mat2), ~spmat2);

            %or
            assertEqual(SparseMat(mat | 0), spmat | 0);
            assertEqual(SparseMat(mat | 2), spmat | 2);
            assertEqual(SparseMat(mat | Inf), spmat | Inf);
            assertEqual(SparseMat(mat | -Inf), spmat | -Inf);
            assertEqual(SparseMat(mat | mat2), spmat | spmat2);
            assertEqual(SparseMat(mat | mat2), spmat | mat2);

            assertEqual(SparseMat(0 | mat), 0 | spmat);
            assertEqual(SparseMat(2 | mat), 2 | spmat);
            assertEqual(SparseMat(Inf | mat), Inf | spmat);
            assertEqual(SparseMat(-Inf | mat), -Inf | spmat);
            assertEqual(SparseMat(mat2 | mat), spmat2 | spmat);
            assertEqual(SparseMat(mat2 | mat), mat2 | spmat);

            %and
            assertEqual(SparseMat(mat & 0), spmat & 0);
            assertEqual(SparseMat(mat & 2), spmat & 2);
            assertEqual(SparseMat(mat & Inf), spmat & Inf);
            assertEqual(SparseMat(mat & -Inf), spmat & -Inf);
            assertEqual(SparseMat(mat & mat2), spmat & spmat2);
            assertEqual(SparseMat(mat & mat2), spmat & mat2);

            assertEqual(SparseMat(0 & mat), 0 & spmat);
            assertEqual(SparseMat(2 & mat), 2 & spmat);
            assertEqual(SparseMat(Inf & mat), Inf & spmat);
            assertEqual(SparseMat(-Inf & mat), -Inf & spmat);
            assertEqual(SparseMat(mat2 & mat), spmat2 & spmat);
            assertEqual(SparseMat(mat2 & mat), mat2 & spmat);

            %xor
            assertEqual(SparseMat(xor(mat,0)), xor(spmat,0));
            assertEqual(SparseMat(xor(mat,2)), xor(spmat,2));
            assertEqual(SparseMat(xor(mat,Inf)), xor(spmat,Inf));
            assertEqual(SparseMat(xor(mat,-Inf)), xor(spmat,-Inf));
            assertEqual(SparseMat(xor(mat,mat2)), xor(spmat,spmat2));
            assertEqual(SparseMat(xor(mat,mat2)), xor(spmat,mat2));

            assertEqual(SparseMat(xor(0,mat)), xor(0,spmat));
            assertEqual(SparseMat(xor(2,mat)), xor(2,spmat));
            assertEqual(SparseMat(xor(Inf,mat)), xor(Inf,spmat));
            assertEqual(SparseMat(xor(-Inf,mat)), xor(-Inf,spmat));
            assertEqual(SparseMat(xor(mat2,mat)), xor(spmat2,spmat));
            assertEqual(SparseMat(xor(mat2,mat)), xor(mat2,spmat));

            %nor
            assertEqual(SparseMat(~(mat | 0)), nor(spmat,0));
            assertEqual(SparseMat(~(mat | 2)), nor(spmat,2));
            assertEqual(SparseMat(~(mat | Inf)), nor(spmat,Inf));
            assertEqual(SparseMat(~(mat | -Inf)), nor(spmat,-Inf));
            assertEqual(SparseMat(~(mat | mat2)), nor(spmat,spmat2));
            assertEqual(SparseMat(~(mat | mat2)), nor(spmat,mat2));

            assertEqual(SparseMat(~(0 | mat)), nor(0,spmat));
            assertEqual(SparseMat(~(2 | mat)), nor(2,spmat));
            assertEqual(SparseMat(~(Inf | mat)), nor(Inf,spmat));
            assertEqual(SparseMat(~(-Inf | mat)), nor(-Inf,spmat));
            assertEqual(SparseMat(~(mat2 | mat)), nor(spmat2,spmat));
            assertEqual(SparseMat(~(mat2 | mat)), nor(mat2,spmat));

            %nand
            assertEqual(SparseMat(~(mat & 0)), nand(spmat,0));
            assertEqual(SparseMat(~(mat & 2)), nand(spmat,2));
            assertEqual(SparseMat(~(mat & Inf)), nand(spmat,Inf));
            assertEqual(SparseMat(~(mat & -Inf)), nand(spmat,-Inf));
            assertEqual(SparseMat(~(mat & mat2)), nand(spmat,spmat2));
            assertEqual(SparseMat(~(mat & mat2)), nand(spmat,mat2));

            assertEqual(SparseMat(~(0 & mat)), nand(0,spmat));
            assertEqual(SparseMat(~(2 & mat)), nand(2,spmat));
            assertEqual(SparseMat(~(Inf & mat)), nand(Inf,spmat));
            assertEqual(SparseMat(~(-Inf & mat)), nand(-Inf,spmat));
            assertEqual(SparseMat(~(mat2 & mat)), nand(spmat2,spmat));
            assertEqual(SparseMat(~(mat2 & mat)), nand(mat2,spmat));
        end

        function testPlusMinus(test)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            spmat2 = SparseMat(mat2);

            %unary
            assertEqual(SparseMat(+mat), +spmat);
            assertEqual(SparseMat(-mat), -spmat);

            %addition
            assertEqual(SparseMat(mat+0), spmat+0);
            assertEqual(SparseMat(mat+2), spmat+2);
            assertEqual(SparseMat(mat+NaN), spmat+NaN);
            assertEqual(SparseMat(mat+Inf), spmat+Inf);
            assertEqual(SparseMat(mat+-Inf), spmat+-Inf);
            assertEqual(SparseMat(mat+mat2), spmat+spmat2);
            assertEqual(SparseMat(mat+mat2), spmat+mat2);

            assertEqual(SparseMat(0+mat), 0+spmat);
            assertEqual(SparseMat(2+mat), 2+spmat);
            assertEqual(SparseMat(NaN+mat), NaN+spmat);
            assertEqual(SparseMat(Inf+mat), Inf+spmat);
            assertEqual(SparseMat(-Inf+mat), -Inf+spmat);
            assertEqual(SparseMat(mat2+mat), spmat2+spmat);
            assertEqual(SparseMat(mat2+mat), mat2+spmat);

            %subtraction
            assertEqual(SparseMat(mat-0), spmat-0);
            assertEqual(SparseMat(mat-2), spmat-2);
            assertEqual(SparseMat(mat-NaN), spmat-NaN);
            assertEqual(SparseMat(mat-Inf), spmat-Inf);
            assertEqual(SparseMat(mat--Inf), spmat--Inf);
            assertEqual(SparseMat(mat-mat2), spmat-spmat2);
            assertEqual(SparseMat(mat-mat2), spmat-mat2);

            assertEqual(SparseMat(0-mat), 0-spmat);
            assertEqual(SparseMat(2-mat), 2-spmat);
            assertEqual(SparseMat(NaN-mat), NaN-spmat);
            assertEqual(SparseMat(Inf-mat), Inf-spmat);
            assertEqual(SparseMat(-Inf-mat), -Inf-spmat);
            assertEqual(SparseMat(mat2-mat), spmat2-spmat);
            assertEqual(SparseMat(mat2-mat), mat2-spmat);
        end

        function testMultiplyDivideExponentiate(test)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            spmat2 = SparseMat(mat2);

            mat3 = [
                2 Inf 0;
                0 2 -3;
                0 0 NaN];
            mat4 = [
                4 0 -7;
                -Inf 3 0
                NaN 0 0];
            spmat3 = SparseMat(mat3);
            spmat4 = SparseMat(mat4);

            %element-wise multiplication
            assertEqual(SparseMat(mat.*0), spmat.*0);
            assertEqual(SparseMat(mat.*2), spmat.*2);
            assertEqual(SparseMat(mat.*NaN), spmat.*NaN);
            assertEqual(SparseMat(mat.*Inf), spmat.*Inf);
            assertEqual(SparseMat(mat.*-Inf), spmat.*-Inf);
            assertEqual(SparseMat(mat.*mat2), spmat.*spmat2);
            assertEqual(SparseMat(mat.*mat2), spmat.*mat2);

            assertEqual(SparseMat(0.*mat), 0.*spmat);
            assertEqual(SparseMat(2.*mat), 2.*spmat);
            assertEqual(SparseMat(NaN.*mat), NaN.*spmat);
            assertEqual(SparseMat(Inf.*mat), Inf.*spmat);
            assertEqual(SparseMat(-Inf.*mat), -Inf.*spmat);
            assertEqual(SparseMat(mat2.*mat), spmat2.*spmat);
            assertEqual(SparseMat(mat2.*mat), mat2.*spmat);

            %matrix multiplication
            assertEqual(SparseMat(mat*0), spmat*0);
            assertEqual(SparseMat(mat*2), spmat*2);
            assertEqual(SparseMat(mat*NaN), spmat*NaN);
            assertEqual(SparseMat(mat*Inf), spmat*Inf);
            assertEqual(SparseMat(mat*-Inf), spmat*-Inf);
            assertEqual(SparseMat(mat3*mat4), spmat3*spmat4);
            assertEqual(SparseMat(mat3*mat4), spmat3*mat4);

            assertEqual(SparseMat(0*mat), 0*spmat);
            assertEqual(SparseMat(2*mat), 2*spmat);
            assertEqual(SparseMat(NaN*mat), NaN*spmat);
            assertEqual(SparseMat(Inf*mat), Inf*spmat);
            assertEqual(SparseMat(-Inf*mat), -Inf*spmat);
            assertEqual(SparseMat(mat4*mat3), spmat4*spmat3);
            assertEqual(SparseMat(mat4*mat3), mat4*spmat3);

            %element-wise right division
            assertEqual(SparseMat(mat./0), spmat./0);
            assertEqual(SparseMat(mat./2), spmat./2);
            assertEqual(SparseMat(mat./NaN), spmat./NaN);
            assertEqual(SparseMat(mat./Inf), spmat./Inf);
            assertEqual(SparseMat(mat./-Inf), spmat./-Inf);
            assertEqual(SparseMat(mat./mat2), spmat./spmat2);
            assertEqual(SparseMat(mat./mat2), spmat./mat2);

            assertEqual(SparseMat(0./mat), 0./spmat);
            assertEqual(SparseMat(2./mat), 2./spmat);
            assertEqual(SparseMat(NaN./mat), NaN./spmat);
            assertEqual(SparseMat(Inf./mat), Inf./spmat);
            assertEqual(SparseMat(-Inf./mat), -Inf./spmat);
            assertEqual(SparseMat(mat2./mat), spmat2./spmat);
            assertEqual(SparseMat(mat2./mat), mat2./spmat);

            %matrix right division
            assertEqual(SparseMat(mat/0), spmat/0);
            assertEqual(SparseMat(mat/2), spmat/2);
            assertEqual(SparseMat(mat/NaN), spmat/NaN);
            assertEqual(SparseMat(mat/Inf), spmat/Inf);
            assertEqual(SparseMat(mat/-Inf), spmat/-Inf);
            % assertEqual(SparseMat(mat/mat2), spmat/spmat2);
            % assertEqual(SparseMat(mat/mat2), spmat/mat2);

            assertEqual(SparseMat(0/mat), 0/spmat);
            assertEqual(SparseMat(2/mat), 2/spmat);
            assertEqual(SparseMat(NaN/mat), NaN/spmat);
            assertEqual(SparseMat(Inf/mat), Inf/spmat);
            assertEqual(SparseMat(-Inf/mat), -Inf/spmat);
            % assertEqual(SparseMat(mat2/mat), spmat2/spmat);
            % assertEqual(SparseMat(mat2/mat), mat2/spmat);

            %element-wise left division
            assertEqual(SparseMat(mat.\0), spmat.\0);
            assertEqual(SparseMat(mat.\2), spmat.\2);
            assertEqual(SparseMat(mat.\NaN), spmat.\NaN);
            assertEqual(SparseMat(mat.\Inf), spmat.\Inf);
            assertEqual(SparseMat(mat.\-Inf), spmat.\-Inf);
            assertEqual(SparseMat(mat.\mat2), spmat.\spmat2);
            assertEqual(SparseMat(mat.\mat2), spmat.\mat2);

            assertEqual(SparseMat(0.\mat), 0.\spmat);
            assertEqual(SparseMat(2.\mat), 2.\spmat);
            assertEqual(SparseMat(NaN.\mat), NaN.\spmat);
            assertEqual(SparseMat(Inf.\mat), Inf.\spmat);
            assertEqual(SparseMat(-Inf.\mat), -Inf.\spmat);
            assertEqual(SparseMat(mat2.\mat), spmat2.\spmat);
            assertEqual(SparseMat(mat2.\mat), mat2.\spmat);

            %matrix left division
            assertEqual(SparseMat(mat\0), spmat\0);
            assertEqual(SparseMat(mat\2), spmat\2);
            assertEqual(SparseMat(mat\NaN), spmat\NaN);
            assertEqual(SparseMat(mat\Inf), spmat\Inf);
            assertEqual(SparseMat(mat\-Inf), spmat\-Inf);
            % assertEqual(SparseMat(mat\mat2), spmat\spmat2);
            % assertEqual(SparseMat(mat\mat2), spmat\mat2);

            assertEqual(SparseMat(0\mat), 0\spmat);
            assertEqual(SparseMat(2\mat), 2\spmat);
            assertEqual(SparseMat(NaN\mat), NaN\spmat);
            assertEqual(SparseMat(Inf\mat), Inf\spmat);
            assertEqual(SparseMat(-Inf\mat), -Inf\spmat);
            % assertEqual(SparseMat(mat2\mat), spmat2\spmat);
            % assertEqual(SparseMat(mat2\mat), mat2\spmat);

            %element-wise exponentiation
            assertEqual(SparseMat(mat.^0), spmat.^0);
            assertEqual(SparseMat(mat.^2), spmat.^2);
            assertEqual(SparseMat(mat.^NaN), spmat.^NaN);
            assertEqual(SparseMat(mat.^Inf), spmat.^Inf);
            assertEqual(SparseMat(mat.^-Inf), spmat.^-Inf);
            assertEqual(SparseMat(mat.^mat2), spmat.^spmat2);
            assertEqual(SparseMat(mat.^mat2), spmat.^mat2);

            assertEqual(SparseMat(0.^mat), 0.^spmat);
            assertEqual(SparseMat(2.^mat), 2.^spmat);
            assertEqual(SparseMat(NaN.^mat), NaN.^spmat);
            assertEqual(SparseMat(Inf.^mat), Inf.^spmat);
            assertEqual(SparseMat(-Inf.^mat), -Inf.^spmat);
            assertEqual(SparseMat(mat2.^mat), spmat2.^spmat);
            assertEqual(SparseMat(mat2.^mat), mat2.^spmat);

            %matrix exponentiation
            % assertEqual(SparseMat(mat^0), spmat^0);
            % assertEqual(SparseMat(mat^2), spmat^2);
            % assertEqual(SparseMat(mat^NaN), spmat^NaN);
            % assertEqual(SparseMat(mat^Inf), spmat^Inf);
            % assertEqual(SparseMat(mat^-Inf), spmat^-Inf);
            % assertEqual(SparseMat(mat^mat2), spmat^spmat2);
            % assertEqual(SparseMat(mat^mat2), spmat^mat2);

            % assertEqual(SparseMat(0^mat), 0^spmat);
            % assertEqual(SparseMat(2^mat), 2^spmat);
            % assertEqual(SparseMat(NaN^mat), NaN^spmat);
            % assertEqual(SparseMat(Inf^mat), Inf^spmat);
            % assertEqual(SparseMat(-Inf^mat), -Inf^spmat);
            % assertEqual(SparseMat(mat2^mat), spmat2^spmat);
            % assertEqual(SparseMat(mat2^mat), mat2^spmat);
        end

        function testEquality(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            spmat2 = SparseMat(mat2);

            assertEqual(SparseMat(mat==0), spmat==0);
            assertEqual(SparseMat(mat==2), spmat==2);
            assertEqual(SparseMat(mat==NaN), spmat==NaN);
            assertEqual(SparseMat(mat==Inf), spmat==Inf);
            assertEqual(SparseMat(mat==-Inf), spmat==-Inf);
            assertEqual(SparseMat(mat==mat2), spmat==spmat2);
            assertEqual(SparseMat(mat==mat2), spmat==mat2);

            assertEqual(SparseMat(0==mat), 0==spmat);
            assertEqual(SparseMat(2==mat), 2==spmat);
            assertEqual(SparseMat(NaN==mat), NaN==spmat);
            assertEqual(SparseMat(Inf==mat), Inf==spmat);
            assertEqual(SparseMat(-Inf==mat), -Inf==spmat);
            assertEqual(SparseMat(mat2==mat), spmat2==spmat);
            assertEqual(SparseMat(mat2==mat), mat2==spmat);
        end

        function testInequality(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            spmat2 = SparseMat(mat2);

            assertEqual(SparseMat(mat~=0), spmat~=0);
            assertEqual(SparseMat(mat~=2), spmat~=2);
            assertEqual(SparseMat(mat~=NaN), spmat~=NaN);
            assertEqual(SparseMat(mat~=Inf), spmat~=Inf);
            assertEqual(SparseMat(mat~=-Inf), spmat~=-Inf);
            assertEqual(SparseMat(mat~=mat2), spmat~=spmat2);
            assertEqual(SparseMat(mat~=mat2), spmat~=mat2);

            assertEqual(SparseMat(0~=mat), 0~=spmat);
            assertEqual(SparseMat(2~=mat), 2~=spmat);
            assertEqual(SparseMat(NaN~=mat), NaN~=spmat);
            assertEqual(SparseMat(Inf~=mat), Inf~=spmat);
            assertEqual(SparseMat(-Inf~=mat), -Inf~=spmat);
            assertEqual(SparseMat(mat2~=mat), spmat2~=spmat);
            assertEqual(SparseMat(mat2~=mat), mat2~=spmat);
        end

        function testGreaterThan(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            spmat2 = SparseMat(mat2);

            assertEqual(SparseMat(mat>0), spmat>0);
            assertEqual(SparseMat(mat>2), spmat>2);
            assertEqual(SparseMat(mat>NaN), spmat>NaN);
            assertEqual(SparseMat(mat>Inf), spmat>Inf);
            assertEqual(SparseMat(mat>-Inf), spmat>-Inf);
            assertEqual(SparseMat(mat>mat2), spmat>spmat2);
            assertEqual(SparseMat(mat>mat2), spmat>mat2);

            assertEqual(SparseMat(0>mat), 0>spmat);
            assertEqual(SparseMat(2>mat), 2>spmat);
            assertEqual(SparseMat(NaN>mat), NaN>spmat);
            assertEqual(SparseMat(Inf>mat), Inf>spmat);
            assertEqual(SparseMat(-Inf>mat), -Inf>spmat);
            assertEqual(SparseMat(mat2>mat), spmat2>spmat);
            assertEqual(SparseMat(mat2>mat), mat2>spmat);
        end

        function testGreaterThanOrEqualTo(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            spmat2 = SparseMat(mat2);

            assertEqual(SparseMat(mat>=0), spmat>=0);
            assertEqual(SparseMat(mat>=2), spmat>=2);
            assertEqual(SparseMat(mat>=NaN), spmat>=NaN);
            assertEqual(SparseMat(mat>=Inf), spmat>=Inf);
            assertEqual(SparseMat(mat>=-Inf), spmat>=-Inf);
            assertEqual(SparseMat(mat>=mat2), spmat>=spmat2);
            assertEqual(SparseMat(mat>=mat2), spmat>=mat2);

            assertEqual(SparseMat(0>=mat), 0>=spmat);
            assertEqual(SparseMat(2>=mat), 2>=spmat);
            assertEqual(SparseMat(NaN>=mat), NaN>=spmat);
            assertEqual(SparseMat(Inf>=mat), Inf>=spmat);
            assertEqual(SparseMat(-Inf>=mat), -Inf>=spmat);
            assertEqual(SparseMat(mat2>=mat), spmat2>=spmat);
            assertEqual(SparseMat(mat2>=mat), mat2>=spmat);
        end

        function testLessThan(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            spmat2 = SparseMat(mat2);

            assertEqual(SparseMat(mat<0), spmat<0);
            assertEqual(SparseMat(mat<2), spmat<2);
            assertEqual(SparseMat(mat<NaN), spmat<NaN);
            assertEqual(SparseMat(mat<Inf), spmat<Inf);
            assertEqual(SparseMat(mat<-Inf), spmat<-Inf);
            assertEqual(SparseMat(mat<mat2), spmat<spmat2);
            assertEqual(SparseMat(mat<mat2), spmat<mat2);

            assertEqual(SparseMat(0<mat), 0<spmat);
            assertEqual(SparseMat(2<mat), 2<spmat);
            assertEqual(SparseMat(NaN<mat), NaN<spmat);
            assertEqual(SparseMat(Inf<mat), Inf<spmat);
            assertEqual(SparseMat(-Inf<mat), -Inf<spmat);
            assertEqual(SparseMat(mat2<mat), spmat2<spmat);
            assertEqual(SparseMat(mat2<mat), mat2<spmat);
        end

        function testLessThanOrEqualTo(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            spmat2 = SparseMat(mat2);

            assertEqual(SparseMat(mat<=0), spmat<=0);
            assertEqual(SparseMat(mat<=2), spmat<=2);
            assertEqual(SparseMat(mat<=NaN), spmat<=NaN);
            assertEqual(SparseMat(mat<=Inf), spmat<=Inf);
            assertEqual(SparseMat(mat<=-Inf), spmat<=-Inf);
            assertEqual(SparseMat(mat<=mat2), spmat<=spmat2);
            assertEqual(SparseMat(mat<=mat2), spmat<=mat2);

            assertEqual(SparseMat(0<=mat), 0<=spmat);
            assertEqual(SparseMat(2<=mat), 2<=spmat);
            assertEqual(SparseMat(NaN<=mat), NaN<=spmat);
            assertEqual(SparseMat(Inf<=mat), Inf<=spmat);
            assertEqual(SparseMat(-Inf<=mat), -Inf<=spmat);
            assertEqual(SparseMat(mat2<=mat), spmat2<=spmat);
            assertEqual(SparseMat(mat2<=mat), mat2<=spmat);
        end

        function testScalarFunctions(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            %abs
            assertEqual(SparseMat(abs(mat)), abs(spmat));

            %sign
            assertEqual(SparseMat(sign(mat)), sign(spmat));

            %sqrt
            assertEqual(SparseMat(sqrt(mat)), sqrt(spmat));

            %realsqrt
            assertEqual(SparseMat(realsqrt(abs(mat))), realsqrt(abs(spmat)));

            %exp
            assertEqual(SparseMat(exp(mat)), exp(spmat));

            %expm1
            assertEqual(SparseMat(expm1(mat)), expm1(spmat));

            %log
            assertEqual(SparseMat(log(mat)), log(spmat));

            %log1p
            assertEqual(SparseMat(log1p(mat)), log1p(spmat));

            %log2
            assertEqual(SparseMat(log2(mat)), log2(spmat));

            %log10
            assertEqual(SparseMat(log10(mat)), log10(spmat));

            %reallog
            assertEqual(SparseMat(reallog(abs(mat))), reallog(abs(spmat)));
        end

        function testComplexNumbers(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3+sqrt(-1);
            mat(3, 1, 1, 1)=4-sqrt(-1);
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            %abs
            assertEqual(SparseMat(abs(mat)), abs(spmat));

            %sign
            assertEqual(SparseMat(sign(mat)), sign(spmat));

            %real
            assertEqual(SparseMat(real(mat)), real(spmat));

            %imag
            assertEqual(SparseMat(imag(mat)), imag(spmat));

            %angle
            assertEqual(SparseMat(angle(mat)), angle(spmat));

            %conj
            assertEqual(SparseMat(conj(mat)), conj(spmat));

            %isreal
            assertEqual(SparseMat(isreal(mat)), isreal(spmat));
        end

        function testSetOperators(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            mat2=ones(4,1);
            spmat2=SparseMat(mat2);

            %unique
            assertEqual(SparseMat(unique(mat)), unique(spmat));
            assertEqual(SparseMat(unique(mat2)), unique(spmat2));
            assertEqual(SparseMat(unique(zeros(0,1))), unique(SparseMat));
            assertEqual(SparseMat(unique(zeros(1,0))), unique(SparseMat));
        end

        function testAnyAll(~)
            import edu.stanford.covert.util.SparseMat;

            fl_x=zeros(3,2);
            fl_x(1,3)=2;
            fl_x(2,:)=-1;

            st_x=SparseMat(fl_x);

            %all
            assertEqual(SparseMat(all(fl_x)), all(st_x));
            assertEqual(SparseMat(all(fl_x,1)), all(st_x,1));
            assertEqual(SparseMat(all(fl_x,2)), all(st_x,2));
            assertEqual(SparseMat(all(fl_x,3)), all(st_x,3));

            %any
            assertEqual(SparseMat(any(fl_x)), any(st_x));
            assertEqual(SparseMat(any(fl_x,1)), any(st_x,1));
            assertEqual(SparseMat(any(fl_x,2)), any(st_x,2));
            assertEqual(SparseMat(any(fl_x,3)), any(st_x,3));
        end

        function testReduction(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=-1;
            mat2(2, 1, 3, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(1, 1, 2, 1)=4;
            mat2(3, 1, 2, 2)=Inf;
            mat2(4, 1, 2, 1)=NaN;
            mat2(3, 1, 2, 1)=-Inf;
            spmat2=SparseMat(mat2);

            %min
            assertEqual(SparseMat(min(mat,0)), min(spmat,0));
            assertEqual(SparseMat(min(mat,2)), min(spmat,2));
            assertEqual(SparseMat(min(mat,-2)), min(spmat,-2));
            assertEqual(SparseMat(min(mat,NaN)), min(spmat,NaN));
            assertEqual(SparseMat(min(mat,Inf)), min(spmat,Inf));
            assertEqual(SparseMat(min(mat,-Inf)), min(spmat,-Inf));

            assertEqual(SparseMat(min(0,mat)), min(0,spmat));
            assertEqual(SparseMat(min(2,mat)), min(2,spmat));
            assertEqual(SparseMat(min(-2,mat)), min(-2,spmat));
            assertEqual(SparseMat(min(NaN,mat)), min(NaN,spmat));
            assertEqual(SparseMat(min(Inf,mat)), min(Inf,spmat));
            assertEqual(SparseMat(min(-Inf,mat)), min(-Inf, spmat));

            assertEqual(SparseMat(min(mat,mat2)), min(spmat,spmat2));
            assertEqual(SparseMat(min(mat,mat2)), min(spmat,mat2));
            assertEqual(SparseMat(min(mat2,mat)), min(spmat2,spmat));
            assertEqual(SparseMat(min(mat,mat2)), min(mat2,spmat));

            assertEqual(SparseMat(min(mat)), min(spmat));
            assertEqual(SparseMat(min(mat,[],1)), min(spmat,[],1));
            assertEqual(SparseMat(min(mat,[],2)), min(spmat,[],2));
            assertEqual(SparseMat(min(mat,[],3)), min(spmat,[],3));
            assertEqual(SparseMat(min(mat,[],4)), min(spmat,[],4));
            assertEqual(SparseMat(min(mat,[],5)), min(spmat,[],5));

            %max
            assertEqual(SparseMat(max(mat,0)), max(spmat,0));
            assertEqual(SparseMat(max(mat,2)), max(spmat,2));
            assertEqual(SparseMat(max(mat,-2)), max(spmat,-2));
            assertEqual(SparseMat(max(mat,NaN)), max(spmat,NaN));
            assertEqual(SparseMat(max(mat,Inf)), max(spmat,Inf));
            assertEqual(SparseMat(max(mat,-Inf)), max(spmat,-Inf));

            assertEqual(SparseMat(max(0,mat)), max(0,spmat));
            assertEqual(SparseMat(max(2,mat)), max(2,spmat));
            assertEqual(SparseMat(max(-2,mat)), max(-2,spmat));
            assertEqual(SparseMat(max(NaN,mat)), max(NaN,spmat));
            assertEqual(SparseMat(max(Inf,mat)), max(Inf,spmat));
            assertEqual(SparseMat(max(-Inf,mat)), max(-Inf, spmat));

            assertEqual(SparseMat(max(mat,mat2)), max(spmat,spmat2));
            assertEqual(SparseMat(max(mat,mat2)), max(spmat,mat2));
            assertEqual(SparseMat(max(mat2,mat)), max(spmat2,spmat));
            assertEqual(SparseMat(max(mat,mat2)), max(mat2,spmat));

            assertEqual(SparseMat(max(mat)), max(spmat));
            assertEqual(SparseMat(max(mat,[],1)), max(spmat,[],1));
            assertEqual(SparseMat(max(mat,[],2)), max(spmat,[],2));
            assertEqual(SparseMat(max(mat,[],3)), max(spmat,[],3));
            assertEqual(SparseMat(max(mat,[],4)), max(spmat,[],4));
            assertEqual(SparseMat(max(mat,[],5)), max(spmat,[],5));

            %range
            assertEqual(SparseMat(range(mat)), range(spmat));
            assertEqual(SparseMat(range(mat,1)), range(spmat,1));
            assertEqual(SparseMat(range(mat,2)), range(spmat,2));
            assertEqual(SparseMat(range(mat,3)), range(spmat,3));
            assertEqual(SparseMat(range(mat,4)), range(spmat,4));
            assertEqual(SparseMat(range(mat,5)), range(spmat,5));

            %sum
            assertEqual(SparseMat(sum(mat)), sum(spmat));
            assertEqual(SparseMat(sum(mat,1)), sum(spmat,1));
            assertEqual(SparseMat(sum(mat,2)), sum(spmat,2));
            assertEqual(SparseMat(sum(mat,3)), sum(spmat,3));
            assertEqual(SparseMat(sum(mat,4)), sum(spmat,4));
            assertEqual(SparseMat(sum(mat,5)), sum(spmat,5));
            assertEqual(SparseMat(sum(true(5, 1))), sum(SparseMat(true(5, 1))));

            %collapse
            assertEqual(sum(sum(sum(sum(mat)))), collapse(spmat));
            assertEqual(squeeze(sum(mat,2)), collapse(spmat,2));
            assertEqual(squeeze(sum(sum(mat,1),3)), collapse(spmat,[1 3]));
            assertEqual(squeeze(sum(sum(sum(mat,3),4),2)), collapse(spmat,[2 3 4]));
            assertEqual(squeeze(sum(sum(sum(mat,2),3),4)), collapse(spmat,[2 3 4]));
            assertEqual(squeeze(sum(sum(sum(mat,1),3),4)), collapse(spmat,[1 3 4]));

            %mean
            assertEqual(SparseMat(mean(mat)), mean(spmat));
            assertEqual(SparseMat(mean(mat,1)), mean(spmat,1));
            assertEqual(SparseMat(mean(mat,2)), mean(spmat,2));
            assertEqual(SparseMat(mean(mat,3)), mean(spmat,3));
            assertEqual(SparseMat(mean(mat,4)), mean(spmat,4));
            assertEqual(SparseMat(mean(mat,5)), mean(spmat,5));

            %var
            assertElementsAlmostEqual(var(mat), double(var(spmat)));
            assertElementsAlmostEqual(var(mat,1), double(var(spmat,1)));

            % assertElementsAlmostEqual(var(mat,0,1), double(var(spmat,0,1)));
            % assertElementsAlmostEqual(var(mat,0,2), double(var(spmat,0,2)));
            % assertElementsAlmostEqual(var(mat,0,3), double(var(spmat,0,3)));
            % assertElementsAlmostEqual(var(mat,0,4), double(var(spmat,0,4)));
            % assertElementsAlmostEqual(var(mat,0,5), double(var(spmat,0,5)));

            assertElementsAlmostEqual(var(mat,1,1), double(var(spmat,1,1)));
            assertElementsAlmostEqual(var(mat,1,2), double(var(spmat,1,2)));
            assertElementsAlmostEqual(var(mat,1,3), double(var(spmat,1,3)));
            assertElementsAlmostEqual(var(mat,1,4), double(var(spmat,1,4)));
            assertElementsAlmostEqual(var(mat,1,5), double(var(spmat,1,5)));

            %std
            assertEqual(SparseMat(std(mat)), std(spmat));
            assertEqual(SparseMat(std(mat,1)), std(spmat,1));

            % assertEqual(SparseMat(std(mat,0,1)), std(spmat,0,1));
            % assertEqual(SparseMat(std(mat,0,2)), std(spmat,0,2));
            % assertEqual(SparseMat(std(mat,0,3)), std(spmat,0,3));
            % assertEqual(SparseMat(std(mat,0,4)), std(spmat,0,4));
            % assertEqual(SparseMat(std(mat,0,5)), std(spmat,0,5));

            assertEqual(SparseMat(std(mat,1,1)), std(spmat,1,1));
            assertEqual(SparseMat(std(mat,1,2)), std(spmat,1,2));
            assertEqual(SparseMat(std(mat,1,3)), std(spmat,1,3));
            assertEqual(SparseMat(std(mat,1,4)), std(spmat,1,4));
            assertEqual(SparseMat(std(mat,1,5)), std(spmat,1,5));

            %norm
            assertEqual(SparseMat(norm(reshape(mat,[],1),1)), norm(spmat,1));
            assertEqual(SparseMat(norm(reshape(mat,[],1),2)), norm(spmat,2));
            assertEqual(SparseMat(norm(reshape(mat,[],1),Inf)), norm(spmat,Inf));
            assertEqual(SparseMat(norm(reshape(mat,[],1),-Inf)), norm(spmat,-Inf));
            assertEqual(SparseMat(norm(reshape(mat,[],1),'fro')), norm(spmat,'fro'));
        end

        function testSum(~)
            import edu.stanford.covert.util.SparseMat;

            x = SparseMat([],[],[3 1]);
            assertEqual(x, sum(x,3));
            assertEqual(x, sum(x,2));
            assertEqual(SparseMat([],[],[1 1]), sum(x,1));

            x = SparseMat([],[],[1 3]);
            assertEqual(x, sum(x,3));
            assertEqual(x, sum(x,1));
            assertEqual(SparseMat([],[],[1 1]), sum(x,2));
        end

        function testIss(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            assertEqual(SparseMat(isnan(mat)), isnan(spmat));
            assertEqual(SparseMat(isinf(mat)), isinf(spmat));
            assertEqual(SparseMat(isfinite(mat)), isfinite(spmat));
        end

        %type casting
        function testTypeCasting(~)
            import edu.stanford.covert.util.SparseMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            spmat = SparseMat(mat);

            mat2 = sparse(permute(sum(mat,4),[1 3 2]));
            spmat2 = SparseMat(mat2);

            %full
            assertEqual(mat, full(spmat));

            %sparse
            assertEqual(mat2, sparse(spmat2));

            %cast
            assertEqual(cast(mat,'double'), cast(spmat,'double'));
            assertEqual(cast(mat,'int32'), cast(spmat,'int32'));

            warning('off','MATLAB:nonIntegerTruncatedInConversionToChar')
            assertEqual(char(mat), char(spmat));
            warning('on','MATLAB:nonIntegerTruncatedInConversionToChar')

            exception=[];
            try
                assertEqual(logical(mat), logical(spmat));
            catch exception
            end
            if isempty(exception)
                throw(MException('SparseMat:typeConversion','NaN''s cannot be converted to logicals'));
            end
            assertEqual(int8(mat), int8(spmat));
            assertEqual(int16(mat), int16(spmat));
            assertEqual(int32(mat), int32(spmat));
            assertEqual(int64(mat), int64(spmat));
            assertEqual(uint8(mat), uint8(spmat));
            assertEqual(uint16(mat), uint16(spmat));
            assertEqual(uint32(mat), uint32(spmat));
            assertEqual(uint64(mat), uint64(spmat));
            assertEqual(single(mat), single(spmat));
            assertEqual(double(mat), double(spmat));
        end
        
        function testSort_subs(~)
            import edu.stanford.covert.util.SparseMat;
            
            for i = 1:1000
                A = randi(100, [6 3]);
                
                assertEqual(sortrows(A, size(A, 2):-1:1), SparseMat.sort_subs(A, max(A, [], 1), size(A, 2):-1:1));
                assertEqual(sortrows(A, 1:size(A, 2)), SparseMat.sort_subs(A, max(A, [], 1), 1:size(A, 2)));
                assertEqual(sortrows(A, size(A, 2):-1:1), SparseMat.sort_subs(A, max(A, [], 1)));
                assertEqual(sortrows(A, [3 1 2]), SparseMat.sort_subs(A, max(A, [], 1), [3 1 2]));
                assertEqual(sortrows(A, [2 1 3]), SparseMat.sort_subs(A, max(A, [], 1), [2 1 3]));
                assertEqual(sortrows(A, [2 3 1]), SparseMat.sort_subs(A, max(A, [], 1), [2 3 1]));
                assertEqual(sortrows(A, [1 3 2]), SparseMat.sort_subs(A, max(A, [], 1), [1 3 2]));
                
                assertEqual(sortrows(A, [1 2]), SparseMat.sort_subs(A, max(A, [], 1), [1 2]));
                assertEqual(sortrows(A, [1 3]), SparseMat.sort_subs(A, max(A, [], 1), [1 3]));
                assertEqual(sortrows(A, [2 1]), SparseMat.sort_subs(A, max(A, [], 1), [2 1]));
                assertEqual(sortrows(A, [2 3]), SparseMat.sort_subs(A, max(A, [], 1), [2 3]));
                assertEqual(sortrows(A, [3 1]), SparseMat.sort_subs(A, max(A, [], 1), [3 1]));
                assertEqual(sortrows(A, [3 2]), SparseMat.sort_subs(A, max(A, [], 1), [3 2]));
                
                assertEqual(sortrows(A, 1), SparseMat.sort_subs(A, max(A, [], 1), 1));
                assertEqual(sortrows(A, 2), SparseMat.sort_subs(A, max(A, [], 1), 2));
                assertEqual(sortrows(A, 3), SparseMat.sort_subs(A, max(A, [], 1), 3));
            end
        end
        
        function testDiskIO(~)
            import edu.stanford.covert.util.SparseMat;
            
            for i = 1:10
                A = randi(100, [20 20 20]);
                
                S = SparseMat(A);
                
                SparseMat.toDisk(S,'spmat.bin');
                assertEqual(S(4:8,12:14,18:20),SparseMat.fromDisk('spmat.bin', [4 8], [12 14], [18 20]));
                assertEqual(S(5:10,5:10,5:10),SparseMat.fromDisk('spmat.bin', [5 10], [5 10], [5 10]));
                assertEqual(S(8:8,14:14,20:20),SparseMat.fromDisk('spmat.bin', [8 8], [14 14], [20 20]));
                delete('spmat.bin');
            end
        end
    end
end
