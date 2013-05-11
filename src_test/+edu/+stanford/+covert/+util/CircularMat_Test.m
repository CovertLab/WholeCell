% CircularMat test cases
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/5/2010
classdef CircularMat_Test < TestCase
    %constructor
    methods
        function this = CircularMat_Test(name)
            this = this@TestCase(name);
        end
    end

    %tests
    methods
        function testConstructor(this)
            import edu.stanford.covert.util.CircularMat;

            %0 arguments
            cmat = CircularMat();

            %1 argument
            mat = zeros(4,3,2);
            mat(4, 3, 1)=1;
            mat(3, 2, 2)=2;
            mat(1, 1, 2)=3;
            mat(3, 1, 1)=4;

            cmat = CircularMat(mat);

            %2 arguments
            mat = zeros(4,3,2);
            mat(4, 3, 1)=1;
            mat(3, 2, 2)=2;
            mat(1, 1, 2)=3;
            mat(3, 1, 1)=4;

            cmat = CircularMat(mat,1);
        end

        function testSubscriptReferenceAssignment(this)
            import edu.stanford.covert.util.CircularMat;

            %reference
            mat = zeros(100,2);
            mat(2,2)=1;
            cmat = CircularMat(mat,1);
            cmat2 = CircularMat(mat,2);
            assertEqual(CircularMat([0; 0; 0; 1]), cmat(99:102,2));
            assertEqual(CircularMat([0;0;0;1]), cmat(99:102,2));

            assertEqual(CircularMat(zeros(1,2)), cmat(1,:));
            assertEqual(CircularMat(zeros(1,2),2), cmat2(1,:));
            assertEqual(CircularMat(zeros(100,1),1), cmat(:,1));
            assertEqual(CircularMat([0;1;zeros(98,1)],1), cmat(:,2));
            assertEqual(CircularMat(zeros(100,1)), cmat(1:100,1));
            assertEqual(CircularMat([0;1;zeros(98,1)]), cmat(1:100,2));

            %assignment
            cmat = CircularMat(zeros(100,2),1);
            cmat(102,2)=1;
            assertEqual(CircularMat([0 0; 0 1; zeros(98,2)],1), cmat);

            cmat = CircularMat(zeros(100,2),1);
            cmat(100,3)=1;
            mat = zeros(100,2);
            mat(100,3)=1;
            assertEqual(CircularMat(mat,1), cmat);
        end

        function testEnd(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            assertEqual(CircularMat(mat(2:end,:,:,:)), cmat(2:end,:,:,:));
            assertEqual(CircularMat(mat(:,:,end,:)), cmat(:,:,end,:));
            assertEqual(CircularMat(mat(:,:,:,end-1)), cmat(:,:,:,end-1));
            assertEqual(CircularMat(mat(1:end)), cmat(1:end));
        end

        function testSize(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            %isempty
            assertEqual(false, isempty(cmat));
            assertEqual(true, isempty(CircularMat));

            %length
            assertEqual(4, length(cmat));

            %size
            assertEqual([4 1 3 2], size(cmat));
            assertEqual(4, size(cmat,1));
            assertEqual(1, size(cmat,5));

            %numel
            assertEqual(4*1*3*2, numel(cmat));

            %ndims
            assertEqual(4, ndims(cmat));

            %permute
            assertEqual(CircularMat(permute(mat,[2 4 1 3])), permute(cmat, [2 4 1 3]));
            assertEqual(CircularMat(permute(mat,[2 5 4 1 3])), permute(cmat, [2 5 4 1 3]));

            %transpose
            mat2 = permute(sum(mat,4),[1 3 2]);
            assertEqual(CircularMat(mat2'), transpose(CircularMat(mat2)));
            assertEqual(CircularMat(mat2'), CircularMat(mat2)');

            %ctranspose
            mat2 = permute(sum(mat,4),[1 3 2])+sqrt(-1);
            assertEqual(CircularMat(ctranspose(mat2)), ctranspose(CircularMat(mat2)));

            %reshape
            assertEqual(CircularMat(reshape(mat,[3 2 4])), reshape(CircularMat(mat),[3 2 4]));
            assertEqual(CircularMat(reshape(zeros(0,2),[0 3])), reshape(CircularMat(zeros(0,2)),[0 3]));

            %squeeze
            assertEqual(CircularMat(squeeze(mat)), squeeze(cmat));
        end

        function testFind(this)
            import edu.stanford.covert.util.CircularMat;

            %example 1
            mat = zeros(4,3);
            mat(4, 3)=1;
            mat(3, 2)=2;
            mat(1, 1)=3;
            mat(3, 1)=4;
            cmat = CircularMat(mat);

            [i1,j1,k1]=find(mat);
            [i2,j2,k2]=find(cmat);

            assertEqual(i1, i2);
            assertEqual(j1, j2);
            assertEqual(k1, k2);

            %example 2
            mat = zeros(4,3,'int32');
            mat(4, 3)=1;
            mat(3, 2)=2;
            mat(1, 1)=3;
            mat(3, 1)=4;
            cmat = CircularMat(mat);

            [i1,j1,k1]=find(mat);
            [i2,j2,k2]=find(cmat);

            assertEqual(i1, i2);
            assertEqual(j1, j2);
            assertEqual(k1, k2);
        end

        %concatenation
        function testConcatenation(this)
            import edu.stanford.covert.util.CircularMat;

            %Vertcat_nullIdentity
            assertEqual(...
                CircularMat,...
                vertcat(CircularMat));

            %Vertcat_identity
            t = CircularMat([1 2; 3 3], 1);
            assertEqual(t, vertcat(t));

            %Vertcat_threeTensors
            u = CircularMat([1 2; 3 3], 1);
            v = CircularMat([6 3; 4 2; 1 4], 1);

            assertEqual(CircularMat([...
                1 2;
                3 3;
                6 3;
                4 2;
                1 4;
                1 2;
                3 3]),...
                [u; v; u]);

            %Vertcat_nullFirst
            t = CircularMat([1 2; 3 3], 1);
            assertEqual(...
                CircularMat([1 2; 3 3]),...
                vertcat(CircularMat(), t));

            %Vertcat_nullLast
            t = CircularMat([1 2; 3 3], 1);
            assertEqual(...
                CircularMat([1 2; 3 3]),...
                vertcat(t, CircularMat()));

            %Vertcat_zeros
            assertEqual(...
                CircularMat(zeros(7,2)),...
                vertcat(...
                CircularMat(zeros(3,2)), CircularMat(zeros(4,2))));

            %Horzcat_nullIdentity
            assertEqual(...
                CircularMat,...
                horzcat(CircularMat));

            %Horzcat_identity
            t = CircularMat([1 2; 3 3], 1);
            assertEqual(t, horzcat(t));


            %Horzcat_nullFirst
            t = CircularMat([1 2; 3 3], 1);
            assertEqual(...
                CircularMat([zeros(2,3) [1 2; 3 3]]),...
                horzcat(CircularMat(zeros(2,3)), t));

            %Horzcat_nullLast
            t = CircularMat([1 2; 3 3], 1);
            assertEqual(...
                CircularMat([[1 2; 3 3] zeros(2,3)],1),...
                horzcat(t, CircularMat(zeros(2,3))));

            %Horzcat_zeros
            assertEqual(...
                CircularMat(zeros([2 7])),...
                horzcat(...
                CircularMat(zeros([2 3])), CircularMat(zeros([2 4]))));

            %padarray
            x = reshape(1:15,5,3);
            cmat = CircularMat(x);

            assertEqual(CircularMat(padarray(x,3)),padarray(cmat,3));
            assertEqual(CircularMat(padarray(x,[0 2])),padarray(cmat,[0 2]));
            assertEqual(CircularMat(padarray(x,[1 0 2])),padarray(cmat,[1 0 2]));
            assertEqual(CircularMat(padarray(x,[1 0 0])),padarray(cmat,[1 0 0]));
            assertEqual(CircularMat(padarray(x,[0 1 1],1)),padarray(cmat,[0 1 1],1));
            assertEqual(CircularMat(padarray(x,[0 1 1],2,'pre')),padarray(cmat,[0 1 1],2,'pre'));
            assertEqual(CircularMat(padarray(x,[0 1 1],2,'post')),padarray(cmat,[0 1 1],2,'post'));
            assertEqual(CircularMat(padarray(x,[0 1 1],2,'both')),padarray(cmat,[0 1 1],2,'both'));

            assertEqual(CircularMat(padarray(x,[2 2],'circular','pre')),padarray(cmat,[2 2],'circular','pre'));
            assertEqual(CircularMat(padarray(x,[2 2],'replicate','pre')),padarray(cmat,[2 2],'replicate','pre'));
            assertEqual(CircularMat(padarray(x,[2 2],'symmetric','pre')),padarray(cmat,[2 2],'symmetric','pre'));
            assertEqual(CircularMat(padarray(x,[2 2],'circular','post')),padarray(cmat,[2 2],'circular','post'));
            assertEqual(CircularMat(padarray(x,[2 2],'replicate','post')),padarray(cmat,[2 2],'replicate','post'));
            assertEqual(CircularMat(padarray(x,[2 2],'symmetric','post')),padarray(cmat,[2 2],'symmetric','post'));
            assertEqual(CircularMat(padarray(x,[2 2],'circular','both')),padarray(cmat,[2 2],'circular','both'));
            assertEqual(CircularMat(padarray(x,[2 2],'replicate','both')),padarray(cmat,[2 2],'replicate','both'));
            assertEqual(CircularMat(padarray(x,[2 2],'symmetric','both')),padarray(cmat,[2 2],'symmetric','both'));
        end

        function testAndOrXorNot(test)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            %mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            cmat2 = CircularMat(mat2);

            %not
            assertEqual(CircularMat(~mat), ~cmat);
            assertEqual(CircularMat(~mat2), ~cmat2);

            %or
            assertEqual(CircularMat(mat | 0), cmat | 0);
            assertEqual(CircularMat(mat | 2), cmat | 2);
            assertEqual(CircularMat(mat | Inf), cmat | Inf);
            assertEqual(CircularMat(mat | -Inf), cmat | -Inf);
            assertEqual(CircularMat(mat | mat2), cmat | cmat2);
            assertEqual(CircularMat(mat | mat2), cmat | mat2);

            assertEqual(CircularMat(0 | mat), 0 | cmat);
            assertEqual(CircularMat(2 | mat), 2 | cmat);
            assertEqual(CircularMat(Inf | mat), Inf | cmat);
            assertEqual(CircularMat(-Inf | mat), -Inf | cmat);
            assertEqual(CircularMat(mat2 | mat), cmat2 | cmat);
            assertEqual(CircularMat(mat2 | mat), mat2 | cmat);

            %and
            assertEqual(CircularMat(mat & 0), cmat & 0);
            assertEqual(CircularMat(mat & 2), cmat & 2);
            assertEqual(CircularMat(mat & Inf), cmat & Inf);
            assertEqual(CircularMat(mat & -Inf), cmat & -Inf);
            assertEqual(CircularMat(mat & mat2), cmat & cmat2);
            assertEqual(CircularMat(mat & mat2), cmat & mat2);

            assertEqual(CircularMat(0 & mat), 0 & cmat);
            assertEqual(CircularMat(2 & mat), 2 & cmat);
            assertEqual(CircularMat(Inf & mat), Inf & cmat);
            assertEqual(CircularMat(-Inf & mat), -Inf & cmat);
            assertEqual(CircularMat(mat2 & mat), cmat2 & cmat);
            assertEqual(CircularMat(mat2 & mat), mat2 & cmat);

            %xor
            assertEqual(CircularMat(xor(mat,0)), xor(cmat,0));
            assertEqual(CircularMat(xor(mat,2)), xor(cmat,2));
            assertEqual(CircularMat(xor(mat,Inf)), xor(cmat,Inf));
            assertEqual(CircularMat(xor(mat,-Inf)), xor(cmat,-Inf));
            assertEqual(CircularMat(xor(mat,mat2)), xor(cmat,cmat2));
            assertEqual(CircularMat(xor(mat,mat2)), xor(cmat,mat2));

            assertEqual(CircularMat(xor(0,mat)), xor(0,cmat));
            assertEqual(CircularMat(xor(2,mat)), xor(2,cmat));
            assertEqual(CircularMat(xor(Inf,mat)), xor(Inf,cmat));
            assertEqual(CircularMat(xor(-Inf,mat)), xor(-Inf,cmat));
            assertEqual(CircularMat(xor(mat2,mat)), xor(cmat2,cmat));
            assertEqual(CircularMat(xor(mat2,mat)), xor(mat2,cmat));
        end

        function testPlusMinus(test)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            cmat2 = CircularMat(mat2);

            %unary
            assertEqual(CircularMat(+mat), +cmat);
            assertEqual(CircularMat(-mat), -cmat);

            %addition
            assertEqual(CircularMat(mat+0), cmat+0);
            assertEqual(CircularMat(mat+2), cmat+2);
            assertEqual(CircularMat(mat+NaN), cmat+NaN);
            assertEqual(CircularMat(mat+Inf), cmat+Inf);
            assertEqual(CircularMat(mat+-Inf), cmat+-Inf);
            assertEqual(CircularMat(mat+mat2), cmat+cmat2);
            assertEqual(CircularMat(mat+mat2), cmat+mat2);

            assertEqual(CircularMat(0+mat), 0+cmat);
            assertEqual(CircularMat(2+mat), 2+cmat);
            assertEqual(CircularMat(NaN+mat), NaN+cmat);
            assertEqual(CircularMat(Inf+mat), Inf+cmat);
            assertEqual(CircularMat(-Inf+mat), -Inf+cmat);
            assertEqual(CircularMat(mat2+mat), cmat2+cmat);
            assertEqual(CircularMat(mat2+mat), mat2+cmat);

            %subtraction
            assertEqual(CircularMat(mat-0), cmat-0);
            assertEqual(CircularMat(mat-2), cmat-2);
            assertEqual(CircularMat(mat-NaN), cmat-NaN);
            assertEqual(CircularMat(mat-Inf), cmat-Inf);
            assertEqual(CircularMat(mat--Inf), cmat--Inf);
            assertEqual(CircularMat(mat-mat2), cmat-cmat2);
            assertEqual(CircularMat(mat-mat2), cmat-mat2);

            assertEqual(CircularMat(0-mat), 0-cmat);
            assertEqual(CircularMat(2-mat), 2-cmat);
            assertEqual(CircularMat(NaN-mat), NaN-cmat);
            assertEqual(CircularMat(Inf-mat), Inf-cmat);
            assertEqual(CircularMat(-Inf-mat), -Inf-cmat);
            assertEqual(CircularMat(mat2-mat), cmat2-cmat);
            assertEqual(CircularMat(mat2-mat), mat2-cmat);
        end

        function testMultiplyDivideExponentiate(test)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            cmat2 = CircularMat(mat2);

            mat3 = [
                2 Inf 0;
                0 2 -3;
                0 0 NaN];
            mat4 = [
                4 0 -7;
                -Inf 3 0
                NaN 0 0];
            cmat3 = CircularMat(mat3);
            cmat4 = CircularMat(mat4);

            %element-wise multiplication
            assertEqual(CircularMat(mat.*0), cmat.*0);
            assertEqual(CircularMat(mat.*2), cmat.*2);
            assertEqual(CircularMat(mat.*NaN), cmat.*NaN);
            assertEqual(CircularMat(mat.*Inf), cmat.*Inf);
            assertEqual(CircularMat(mat.*-Inf), cmat.*-Inf);
            assertEqual(CircularMat(mat.*mat2), cmat.*cmat2);
            assertEqual(CircularMat(mat.*mat2), cmat.*mat2);

            assertEqual(CircularMat(0.*mat), 0.*cmat);
            assertEqual(CircularMat(2.*mat), 2.*cmat);
            assertEqual(CircularMat(NaN.*mat), NaN.*cmat);
            assertEqual(CircularMat(Inf.*mat), Inf.*cmat);
            assertEqual(CircularMat(-Inf.*mat), -Inf.*cmat);
            assertEqual(CircularMat(mat2.*mat), cmat2.*cmat);
            assertEqual(CircularMat(mat2.*mat), mat2.*cmat);

            %matrix multiplication
            assertEqual(CircularMat(mat*0), cmat*0);
            assertEqual(CircularMat(mat*2), cmat*2);
            assertEqual(CircularMat(mat*NaN), cmat*NaN);
            assertEqual(CircularMat(mat*Inf), cmat*Inf);
            assertEqual(CircularMat(mat*-Inf), cmat*-Inf);
            assertEqual(CircularMat(mat3*mat4), cmat3*cmat4);
            assertEqual(CircularMat(mat3*mat4), cmat3*mat4);

            assertEqual(CircularMat(0*mat), 0*cmat);
            assertEqual(CircularMat(2*mat), 2*cmat);
            assertEqual(CircularMat(NaN*mat), NaN*cmat);
            assertEqual(CircularMat(Inf*mat), Inf*cmat);
            assertEqual(CircularMat(-Inf*mat), -Inf*cmat);
            assertEqual(CircularMat(mat4*mat3), cmat4*cmat3);
            assertEqual(CircularMat(mat4*mat3), mat4*cmat3);

            %element-wise right division
            assertEqual(CircularMat(mat./0), cmat./0);
            assertEqual(CircularMat(mat./2), cmat./2);
            assertEqual(CircularMat(mat./NaN), cmat./NaN);
            assertEqual(CircularMat(mat./Inf), cmat./Inf);
            assertEqual(CircularMat(mat./-Inf), cmat./-Inf);
            assertEqual(CircularMat(mat./mat2), cmat./cmat2);
            assertEqual(CircularMat(mat./mat2), cmat./mat2);

            assertEqual(CircularMat(0./mat), 0./cmat);
            assertEqual(CircularMat(2./mat), 2./cmat);
            assertEqual(CircularMat(NaN./mat), NaN./cmat);
            assertEqual(CircularMat(Inf./mat), Inf./cmat);
            assertEqual(CircularMat(-Inf./mat), -Inf./cmat);
            assertEqual(CircularMat(mat2./mat), cmat2./cmat);
            assertEqual(CircularMat(mat2./mat), mat2./cmat);

            %matrix right division
            % assertEqual(CircularMat(mat/0), cmat/0);
            % assertEqual(CircularMat(mat/2), cmat/2);
            % assertEqual(CircularMat(mat/NaN), cmat/NaN);
            % assertEqual(CircularMat(mat/Inf), cmat/Inf);
            % assertEqual(CircularMat(mat/-Inf), cmat/-Inf);
            % assertEqual(CircularMat(mat/mat2), cmat/cmat2);
            % assertEqual(CircularMat(mat/mat2), cmat/mat2);

            % assertEqual(CircularMat(0/mat), 0/cmat);
            % assertEqual(CircularMat(2/mat), 2/cmat);
            % assertEqual(CircularMat(NaN/mat), NaN/cmat);
            % assertEqual(CircularMat(Inf/mat), Inf/cmat);
            % assertEqual(CircularMat(-Inf/mat), -Inf/cmat);
            % assertEqual(CircularMat(mat2/mat), cmat2/cmat);
            % assertEqual(CircularMat(mat2/mat), mat2/cmat);

            %element-wise left division
            assertEqual(CircularMat(mat.\0), cmat.\0);
            assertEqual(CircularMat(mat.\2), cmat.\2);
            assertEqual(CircularMat(mat.\NaN), cmat.\NaN);
            assertEqual(CircularMat(mat.\Inf), cmat.\Inf);
            assertEqual(CircularMat(mat.\-Inf), cmat.\-Inf);
            assertEqual(CircularMat(mat.\mat2), cmat.\cmat2);
            assertEqual(CircularMat(mat.\mat2), cmat.\mat2);

            assertEqual(CircularMat(0.\mat), 0.\cmat);
            assertEqual(CircularMat(2.\mat), 2.\cmat);
            assertEqual(CircularMat(NaN.\mat), NaN.\cmat);
            assertEqual(CircularMat(Inf.\mat), Inf.\cmat);
            assertEqual(CircularMat(-Inf.\mat), -Inf.\cmat);
            assertEqual(CircularMat(mat2.\mat), cmat2.\cmat);
            assertEqual(CircularMat(mat2.\mat), mat2.\cmat);

            %matrix left division
            % assertEqual(CircularMat(mat\0), cmat\0);
            % assertEqual(CircularMat(mat\2), cmat\2);
            % assertEqual(CircularMat(mat\NaN), cmat\NaN);
            % assertEqual(CircularMat(mat\Inf), cmat\Inf);
            % assertEqual(CircularMat(mat\-Inf), cmat\-Inf);
            % assertEqual(CircularMat(mat\mat2), cmat\cmat2);
            % assertEqual(CircularMat(mat\mat2), cmat\mat2);

            % assertEqual(CircularMat(0\mat), 0\cmat);
            % assertEqual(CircularMat(2\mat), 2\cmat);
            % assertEqual(CircularMat(NaN\mat), NaN\cmat);
            % assertEqual(CircularMat(Inf\mat), Inf\cmat);
            % assertEqual(CircularMat(-Inf\mat), -Inf\cmat);
            % assertEqual(CircularMat(mat2\mat), cmat2\cmat);
            % assertEqual(CircularMat(mat2\mat), mat2\cmat);

            %element-wise exponentiation
            assertEqual(CircularMat(mat.^0), cmat.^0);
            assertEqual(CircularMat(mat.^2), cmat.^2);
            assertEqual(CircularMat(mat.^NaN), cmat.^NaN);
            assertEqual(CircularMat(mat.^Inf), cmat.^Inf);
            assertEqual(CircularMat(mat.^-Inf), cmat.^-Inf);
            assertEqual(CircularMat(mat.^mat2), cmat.^cmat2);
            assertEqual(CircularMat(mat.^mat2), cmat.^mat2);

            assertEqual(CircularMat(0.^mat), 0.^cmat);
            assertEqual(CircularMat(2.^mat), 2.^cmat);
            assertEqual(CircularMat(NaN.^mat), NaN.^cmat);
            assertEqual(CircularMat(Inf.^mat), Inf.^cmat);
            assertEqual(CircularMat(-Inf.^mat), -Inf.^cmat);
            assertEqual(CircularMat(mat2.^mat), cmat2.^cmat);
            assertEqual(CircularMat(mat2.^mat), mat2.^cmat);

            %matrix exponentiation
            mat = reshape(1:9,3,3);
            mat2 = reshape(10:18,3,3);
            cmat = CircularMat(mat);
            cmat2 = CircularMat(mat2);

            assertEqual(CircularMat(mat^0), cmat^0);
            assertEqual(CircularMat(mat^2), cmat^2);
            assertEqual(CircularMat(mat^NaN), cmat^NaN);
            assertEqual(CircularMat(mat^Inf), cmat^Inf);
            assertEqual(CircularMat(mat^-Inf), cmat^-Inf);
            % assertEqual(CircularMat(mat^mat2), cmat^cmat2);
            % assertEqual(CircularMat(mat^mat2), cmat^mat2);

            assertEqual(CircularMat(0^mat), 0^cmat);
            assertEqual(CircularMat(2^mat), 2^cmat);
            assertEqual(CircularMat(NaN^mat), NaN^cmat);
            assertEqual(CircularMat(Inf^mat), Inf^cmat);
            assertEqual(CircularMat(-Inf^mat), -Inf^cmat);
            % assertEqual(CircularMat(mat2^mat), cmat2^cmat);
            % assertEqual(CircularMat(mat2^mat), mat2^cmat);
        end

        function testEquality(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            cmat2 = CircularMat(mat2);

            assertEqual(CircularMat(mat==0), cmat==0);
            assertEqual(CircularMat(mat==2), cmat==2);
            assertEqual(CircularMat(mat==NaN), cmat==NaN);
            assertEqual(CircularMat(mat==Inf), cmat==Inf);
            assertEqual(CircularMat(mat==-Inf), cmat==-Inf);
            assertEqual(CircularMat(mat==mat2), cmat==cmat2);
            assertEqual(CircularMat(mat==mat2), cmat==mat2);

            assertEqual(CircularMat(0==mat), 0==cmat);
            assertEqual(CircularMat(2==mat), 2==cmat);
            assertEqual(CircularMat(NaN==mat), NaN==cmat);
            assertEqual(CircularMat(Inf==mat), Inf==cmat);
            assertEqual(CircularMat(-Inf==mat), -Inf==cmat);
            assertEqual(CircularMat(mat2==mat), cmat2==cmat);
            assertEqual(CircularMat(mat2==mat), mat2==cmat);
        end

        function testInequality(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            cmat2 = CircularMat(mat2);

            assertEqual(CircularMat(mat~=0), cmat~=0);
            assertEqual(CircularMat(mat~=2), cmat~=2);
            assertEqual(CircularMat(mat~=NaN), cmat~=NaN);
            assertEqual(CircularMat(mat~=Inf), cmat~=Inf);
            assertEqual(CircularMat(mat~=-Inf), cmat~=-Inf);
            assertEqual(CircularMat(mat~=mat2), cmat~=cmat2);
            assertEqual(CircularMat(mat~=mat2), cmat~=mat2);

            assertEqual(CircularMat(0~=mat), 0~=cmat);
            assertEqual(CircularMat(2~=mat), 2~=cmat);
            assertEqual(CircularMat(NaN~=mat), NaN~=cmat);
            assertEqual(CircularMat(Inf~=mat), Inf~=cmat);
            assertEqual(CircularMat(-Inf~=mat), -Inf~=cmat);
            assertEqual(CircularMat(mat2~=mat), cmat2~=cmat);
            assertEqual(CircularMat(mat2~=mat), mat2~=cmat);
        end

        function testGreaterThan(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            cmat2 = CircularMat(mat2);

            assertEqual(CircularMat(mat>0), cmat>0);
            assertEqual(CircularMat(mat>2), cmat>2);
            assertEqual(CircularMat(mat>NaN), cmat>NaN);
            assertEqual(CircularMat(mat>Inf), cmat>Inf);
            assertEqual(CircularMat(mat>-Inf), cmat>-Inf);
            assertEqual(CircularMat(mat>mat2), cmat>cmat2);
            assertEqual(CircularMat(mat>mat2), cmat>mat2);

            assertEqual(CircularMat(0>mat), 0>cmat);
            assertEqual(CircularMat(2>mat), 2>cmat);
            assertEqual(CircularMat(NaN>mat), NaN>cmat);
            assertEqual(CircularMat(Inf>mat), Inf>cmat);
            assertEqual(CircularMat(-Inf>mat), -Inf>cmat);
            assertEqual(CircularMat(mat2>mat), cmat2>cmat);
            assertEqual(CircularMat(mat2>mat), mat2>cmat);
        end

        function testGreaterThanOrEqualTo(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            cmat2 = CircularMat(mat2);

            assertEqual(CircularMat(mat>=0), cmat>=0);
            assertEqual(CircularMat(mat>=2), cmat>=2);
            assertEqual(CircularMat(mat>=NaN), cmat>=NaN);
            assertEqual(CircularMat(mat>=Inf), cmat>=Inf);
            assertEqual(CircularMat(mat>=-Inf), cmat>=-Inf);
            assertEqual(CircularMat(mat>=mat2), cmat>=cmat2);
            assertEqual(CircularMat(mat>=mat2), cmat>=mat2);

            assertEqual(CircularMat(0>=mat), 0>=cmat);
            assertEqual(CircularMat(2>=mat), 2>=cmat);
            assertEqual(CircularMat(NaN>=mat), NaN>=cmat);
            assertEqual(CircularMat(Inf>=mat), Inf>=cmat);
            assertEqual(CircularMat(-Inf>=mat), -Inf>=cmat);
            assertEqual(CircularMat(mat2>=mat), cmat2>=cmat);
            assertEqual(CircularMat(mat2>=mat), mat2>=cmat);
        end

        function testLessThan(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            cmat2 = CircularMat(mat2);

            assertEqual(CircularMat(mat<0), cmat<0);
            assertEqual(CircularMat(mat<2), cmat<2);
            assertEqual(CircularMat(mat<NaN), cmat<NaN);
            assertEqual(CircularMat(mat<Inf), cmat<Inf);
            assertEqual(CircularMat(mat<-Inf), cmat<-Inf);
            assertEqual(CircularMat(mat<mat2), cmat<cmat2);
            assertEqual(CircularMat(mat<mat2), cmat<mat2);

            assertEqual(CircularMat(0<mat), 0<cmat);
            assertEqual(CircularMat(2<mat), 2<cmat);
            assertEqual(CircularMat(NaN<mat), NaN<cmat);
            assertEqual(CircularMat(Inf<mat), Inf<cmat);
            assertEqual(CircularMat(-Inf<mat), -Inf<cmat);
            assertEqual(CircularMat(mat2<mat), cmat2<cmat);
            assertEqual(CircularMat(mat2<mat), mat2<cmat);
        end

        function testLessThanOrEqualTo(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=5;
            mat2(3, 1, 2, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(3, 1, 3, 1)=4;
            cmat2 = CircularMat(mat2);

            assertEqual(CircularMat(mat<=0), cmat<=0);
            assertEqual(CircularMat(mat<=2), cmat<=2);
            assertEqual(CircularMat(mat<=NaN), cmat<=NaN);
            assertEqual(CircularMat(mat<=Inf), cmat<=Inf);
            assertEqual(CircularMat(mat<=-Inf), cmat<=-Inf);
            assertEqual(CircularMat(mat<=mat2), cmat<=cmat2);
            assertEqual(CircularMat(mat<=mat2), cmat<=mat2);

            assertEqual(CircularMat(0<=mat), 0<=cmat);
            assertEqual(CircularMat(2<=mat), 2<=cmat);
            assertEqual(CircularMat(NaN<=mat), NaN<=cmat);
            assertEqual(CircularMat(Inf<=mat), Inf<=cmat);
            assertEqual(CircularMat(-Inf<=mat), -Inf<=cmat);
            assertEqual(CircularMat(mat2<=mat), cmat2<=cmat);
            assertEqual(CircularMat(mat2<=mat), mat2<=cmat);
        end

        function testScalarFunctions(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            %abs
            assertEqual(CircularMat(abs(mat)), abs(cmat));

            %sign
            assertEqual(CircularMat(sign(mat)), sign(cmat));

            %sqrt
            assertEqual(CircularMat(sqrt(mat)), sqrt(cmat));

            %realsqrt
            assertEqual(CircularMat(realsqrt(abs(mat))), realsqrt(abs(cmat)));

            %exp
            assertEqual(CircularMat(exp(mat)), exp(cmat));

            %expm1
            assertEqual(CircularMat(expm1(mat)), expm1(cmat));

            %log
            assertEqual(CircularMat(log(mat)), log(cmat));

            %log1p
            assertEqual(CircularMat(log1p(mat)), log1p(cmat));

            %log2
            assertEqual(CircularMat(log2(mat)), log2(cmat));

            %log10
            assertEqual(CircularMat(log10(mat)), log10(cmat));

            %reallog
            assertEqual(CircularMat(reallog(abs(mat))), reallog(abs(cmat)));
        end

        function testComplexNumbers(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3+sqrt(-1);
            mat(3, 1, 1, 1)=4-sqrt(-1);
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            %abs
            assertEqual(CircularMat(abs(mat)), abs(cmat));

            %sign
            assertEqual(CircularMat(sign(mat)), sign(cmat));

            %real
            assertEqual(CircularMat(real(mat)), real(cmat));

            %imag
            assertEqual(CircularMat(imag(mat)), imag(cmat));

            %angle
            assertEqual(CircularMat(angle(mat)), angle(cmat));

            %conj
            assertEqual(CircularMat(conj(mat)), conj(cmat));

            %isreal
            assertEqual(isreal(mat), isreal(cmat));
        end

        function testSetOperators(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            mat2=ones(4,1);
            cmat2=CircularMat(mat2);

            %unique
            assertEqual(CircularMat(unique(mat)), unique(cmat));
            assertEqual(CircularMat(unique(mat2)), unique(cmat2));
            assertEqual(CircularMat(unique(zeros(0,1))), unique(CircularMat(zeros(0,1))));
            assertEqual(CircularMat(unique(zeros(1,0))), unique(CircularMat(zeros(1,0))));
        end

        function testAnyAll(this)
            import edu.stanford.covert.util.CircularMat;

            fl_x=zeros(3,2);
            fl_x(1,3)=2;
            fl_x(2,:)=-1;

            st_x=CircularMat(fl_x);

            %all
            assertEqual(CircularMat(all(fl_x)), all(st_x));
            assertEqual(CircularMat(all(fl_x,1)), all(st_x,1));
            assertEqual(CircularMat(all(fl_x,2)), all(st_x,2));
            assertEqual(CircularMat(all(fl_x,3)), all(st_x,3));

            %any
            assertEqual(CircularMat(any(fl_x)), any(st_x));
            assertEqual(CircularMat(any(fl_x,1)), any(st_x,1));
            assertEqual(CircularMat(any(fl_x,2)), any(st_x,2));
            assertEqual(CircularMat(any(fl_x,3)), any(st_x,3));
        end

        function testReduction(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            mat2 = zeros(4,1,3,2);
            mat2(4, 1, 3, 1)=-1;
            mat2(2, 1, 3, 2)=2;
            mat2(1, 1, 1, 2)=3;
            mat2(1, 1, 2, 1)=4;
            mat2(3, 1, 2, 2)=Inf;
            mat2(4, 1, 2, 1)=NaN;
            mat2(3, 1, 2, 1)=-Inf;
            cmat2=CircularMat(mat2);

            %min
            assertEqual(CircularMat(min(mat,0)), min(cmat,0));
            assertEqual(CircularMat(min(mat,2)), min(cmat,2));
            assertEqual(CircularMat(min(mat,-2)), min(cmat,-2));
            assertEqual(CircularMat(min(mat,NaN)), min(cmat,NaN));
            assertEqual(CircularMat(min(mat,Inf)), min(cmat,Inf));
            assertEqual(CircularMat(min(mat,-Inf)), min(cmat,-Inf));

            assertEqual(CircularMat(min(0,mat)), min(0,cmat));
            assertEqual(CircularMat(min(2,mat)), min(2,cmat));
            assertEqual(CircularMat(min(-2,mat)), min(-2,cmat));
            assertEqual(CircularMat(min(NaN,mat)), min(NaN,cmat));
            assertEqual(CircularMat(min(Inf,mat)), min(Inf,cmat));
            assertEqual(CircularMat(min(-Inf,mat)), min(-Inf, cmat));

            assertEqual(CircularMat(min(mat,mat2)), min(cmat,cmat2));
            assertEqual(CircularMat(min(mat,mat2)), min(cmat,mat2));
            assertEqual(CircularMat(min(mat2,mat)), min(cmat2,cmat));
            assertEqual(CircularMat(min(mat,mat2)), min(mat2,cmat));

            assertEqual(CircularMat(min(mat)), min(cmat));
            assertEqual(CircularMat(min(mat,[],1)), min(cmat,[],1));
            assertEqual(CircularMat(min(mat,[],2)), min(cmat,[],2));
            assertEqual(CircularMat(min(mat,[],3)), min(cmat,[],3));
            assertEqual(CircularMat(min(mat,[],4)), min(cmat,[],4));
            assertEqual(CircularMat(min(mat,[],5)), min(cmat,[],5));

            %max
            assertEqual(CircularMat(max(mat,0)), max(cmat,0));
            assertEqual(CircularMat(max(mat,2)), max(cmat,2));
            assertEqual(CircularMat(max(mat,-2)), max(cmat,-2));
            assertEqual(CircularMat(max(mat,NaN)), max(cmat,NaN));
            assertEqual(CircularMat(max(mat,Inf)), max(cmat,Inf));
            assertEqual(CircularMat(max(mat,-Inf)), max(cmat,-Inf));

            assertEqual(CircularMat(max(0,mat)), max(0,cmat));
            assertEqual(CircularMat(max(2,mat)), max(2,cmat));
            assertEqual(CircularMat(max(-2,mat)), max(-2,cmat));
            assertEqual(CircularMat(max(NaN,mat)), max(NaN,cmat));
            assertEqual(CircularMat(max(Inf,mat)), max(Inf,cmat));
            assertEqual(CircularMat(max(-Inf,mat)), max(-Inf, cmat));

            assertEqual(CircularMat(max(mat,mat2)), max(cmat,cmat2));
            assertEqual(CircularMat(max(mat,mat2)), max(cmat,mat2));
            assertEqual(CircularMat(max(mat2,mat)), max(cmat2,cmat));
            assertEqual(CircularMat(max(mat,mat2)), max(mat2,cmat));

            assertEqual(CircularMat(max(mat)), max(cmat));
            assertEqual(CircularMat(max(mat,[],1)), max(cmat,[],1));
            assertEqual(CircularMat(max(mat,[],2)), max(cmat,[],2));
            assertEqual(CircularMat(max(mat,[],3)), max(cmat,[],3));
            assertEqual(CircularMat(max(mat,[],4)), max(cmat,[],4));
            assertEqual(CircularMat(max(mat,[],5)), max(cmat,[],5));

            %range
            assertEqual(CircularMat(range(mat)), range(cmat));
            assertEqual(CircularMat(range(mat,1)), range(cmat,1));
            assertEqual(CircularMat(range(mat,2)), range(cmat,2));
            assertEqual(CircularMat(range(mat,3)), range(cmat,3));
            assertEqual(CircularMat(range(mat,4)), range(cmat,4));
            assertEqual(CircularMat(range(mat,5)), range(cmat,5));

            %sum
            assertEqual(CircularMat(sum(mat)), sum(cmat));
            assertEqual(CircularMat(sum(mat,1)), sum(cmat,1));
            assertEqual(CircularMat(sum(mat,2)), sum(cmat,2));
            assertEqual(CircularMat(sum(mat,3)), sum(cmat,3));
            assertEqual(CircularMat(sum(mat,4)), sum(cmat,4));
            assertEqual(CircularMat(sum(mat,5)), sum(cmat,5));

            %mean
            assertEqual(CircularMat(mean(mat)), mean(cmat));
            assertEqual(CircularMat(mean(mat,1)), mean(cmat,1));
            assertEqual(CircularMat(mean(mat,2)), mean(cmat,2));
            assertEqual(CircularMat(mean(mat,3)), mean(cmat,3));
            assertEqual(CircularMat(mean(mat,4)), mean(cmat,4));
            assertEqual(CircularMat(mean(mat,5)), mean(cmat,5));

            %var
            assertElementsAlmostEqual(var(mat), double(var(cmat)));
            assertElementsAlmostEqual(var(mat,1), double(var(cmat,1)));

            % assertElementsAlmostEqual(var(mat,0,1), double(var(cmat,0,1)));
            % assertElementsAlmostEqual(var(mat,0,2), double(var(cmat,0,2)));
            % assertElementsAlmostEqual(var(mat,0,3), double(var(cmat,0,3)));
            % assertElementsAlmostEqual(var(mat,0,4), double(var(cmat,0,4)));
            % assertElementsAlmostEqual(var(mat,0,5), double(var(cmat,0,5)));

            assertElementsAlmostEqual(var(mat,1,1), double(var(cmat,1,1)));
            assertElementsAlmostEqual(var(mat,1,2), double(var(cmat,1,2)));
            assertElementsAlmostEqual(var(mat,1,3), double(var(cmat,1,3)));
            assertElementsAlmostEqual(var(mat,1,4), double(var(cmat,1,4)));
            assertElementsAlmostEqual(var(mat,1,5), double(var(cmat,1,5)));

            %std
            assertEqual(CircularMat(std(mat)), std(cmat));
            assertEqual(CircularMat(std(mat,1)), std(cmat,1));

            % assertEqual(CircularMat(std(mat,0,1)), std(cmat,0,1));
            % assertEqual(CircularMat(std(mat,0,2)), std(cmat,0,2));
            % assertEqual(CircularMat(std(mat,0,3)), std(cmat,0,3));
            % assertEqual(CircularMat(std(mat,0,4)), std(cmat,0,4));
            % assertEqual(CircularMat(std(mat,0,5)), std(cmat,0,5));

            assertEqual(CircularMat(std(mat,1,1)), std(cmat,1,1));
            assertEqual(CircularMat(std(mat,1,2)), std(cmat,1,2));
            assertEqual(CircularMat(std(mat,1,3)), std(cmat,1,3));
            assertEqual(CircularMat(std(mat,1,4)), std(cmat,1,4));
            assertEqual(CircularMat(std(mat,1,5)), std(cmat,1,5));

            %norm
            mat = reshape(1:20,4,5);
            cmat = CircularMat(mat);
            assertEqual(norm(mat,1), norm(cmat,1));
            assertEqual(norm(mat,2), norm(cmat,2));
            assertEqual(norm(mat,Inf), norm(cmat,Inf));
            assertEqual(norm(reshape(mat,[],1),-Inf), norm(reshape(cmat,[],1),-Inf));
            assertEqual(norm(mat,'fro'), norm(cmat,'fro'));
        end

        function testIss(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            assertEqual(CircularMat(isnan(mat)), isnan(cmat));
            assertEqual(CircularMat(isinf(mat)), isinf(cmat));
            assertEqual(CircularMat(isfinite(mat)), isfinite(cmat));
        end

        %type casting
        function testTypeCasting(this)
            import edu.stanford.covert.util.CircularMat;

            mat = zeros(4,1,3,2);
            mat(4, 1, 3, 1)=-1;
            mat(3, 1, 2, 2)=2;
            mat(1, 1, 1, 2)=3;
            mat(3, 1, 1, 1)=4;
            mat(3, 1, 2, 1)=Inf;
            mat(1, 1, 2, 1)=NaN;
            mat(2, 1, 2, 1)=-Inf;
            cmat = CircularMat(mat);

            mat2 = sparse(permute(sum(mat,4),[1 3 2]));
            cmat2 = CircularMat(mat2);

            %sparse
            assertEqual(mat2, sparse(cmat2));

            %cast
            assertEqual(cast(mat,'double'), cast(cmat,'double'));
            assertEqual(cast(mat,'int32'), cast(cmat,'int32'));

            warning('off','MATLAB:nonIntegerTruncatedInConversionToChar')
            assertEqual(char(mat), char(cmat));
            warning('on','MATLAB:nonIntegerTruncatedInConversionToChar')

            exception=[];
            try
                assertEqual(logical(mat), logical(cmat));
            catch exception
            end
            if isempty(exception)
                throw(MException('CircularMat:typeConversion','NaN''s cannot be converted to logicals'));
            end
            assertEqual(int8(mat), int8(cmat));
            assertEqual(int16(mat), int16(cmat));
            assertEqual(int32(mat), int32(cmat));
            assertEqual(int64(mat), int64(cmat));
            assertEqual(uint8(mat), uint8(cmat));
            assertEqual(uint16(mat), uint16(cmat));
            assertEqual(uint32(mat), uint32(cmat));
            assertEqual(uint64(mat), uint64(cmat));
            assertEqual(single(mat), single(cmat));
            assertEqual(double(mat), double(cmat));
        end
    end
end