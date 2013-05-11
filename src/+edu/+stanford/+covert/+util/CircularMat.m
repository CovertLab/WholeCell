% CircularMat
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/18/2010
classdef CircularMat
    properties (SetAccess=protected)
        mat
        circularDims
    end

    methods
        function this = CircularMat(mat, circularDims)
            if ~exist('mat','var')
                mat = [];
            end

            if ~exist('circularDims','var')
                circularDims = zeros(1,0);
            end

            if ...
                    any(circularDims < 1) || ...
                    ~isequal(circularDims, ceil(real(circularDims))) || ...
                    any(circularDims > ndims(mat)) || ...
                    ~isvector(circularDims) || ...
                    numel(unique(circularDims)) < numel(circularDims)
                throw(MException('CircularMat:invalidDimensions',...
                    'circular dimensions must be unique positive real integers <= ndims(mat)'));
            end

            this.mat = mat;
            this.circularDims = reshape(sort(circularDims),1,[]);
        end
    end

    %size
    methods
        function value = isempty(this)
            value = isempty(this.mat);
        end

        function value = numel(this)
            value = numel(this.mat);
        end

        function value = length(this)
            value = length(this.mat);
        end

        function value = size(this, varargin)
            value = size(this.mat, varargin{:});
        end

        function value = ndims(this)
            value = ndims(this.mat);
        end

        function this = reshape(this, varargin)
            if nargin == 3 && numel(varargin{1})>1
                this.mat = reshape(this.mat, varargin{1});
                if ...
                        any(varargin{2} < 1) || ...
                        ~isequal(varargin{2}, ceil(real(varargin{2}))) || ...
                        any(varargin{2} > ndims(this.mat)) || ...
                        ~isvector(varargin{2}) || ...
                        numel(unique(varargin{2})) < numel(varargin{2})
                    throw(MException('CircularMat:invalidDimensions',...
                        'circular dimensions must be unique positive real integers <= ndims(mat)'));
                end

                this.circularDims = sort(reshape(varargin{2},1,[]));
            else
                this.mat = reshape(this.mat, varargin{:});
                this.circularDims = zeros(1,0);
            end
        end

        function this = transpose(this)
            this.mat=transpose(this.mat);
            this.circularDims(this.circularDims==1)=-1;
            this.circularDims(this.circularDims==2)=1;
            this.circularDims(this.circularDims==-1)=2;
            this.circularDims = sort(this.circularDims);
        end

        function this = ctranspose(this)
            this = conj(transpose(this));
        end

        function this = permute(this, dims)
            this.mat=permute(this.mat, dims);
            this.circularDims = find(ismember(dims, this.circularDims));
        end

        function this = squeeze(this)
            [~,this.circularDims] = intersect(find(size(this.mat)==1), this.circularDims);
            this.mat=squeeze(this.mat);
        end

    end

    %subscript reference, assignment
    methods
        function value = subsref(this, s)
            import edu.stanford.covert.util.CircularMat;

            %Case 1: Dot reference to properties
            if strcmp(s(1).type,'.')
                value = builtin('subsref', this, s);
                return;
            end

            if numel(s.subs)==1
                %Case 2: linear indexing
                value = CircularMat(subsref(this.mat, s));
            else
                %Case 3: subscript indexing
                keepCircularDims = false(size(this.circularDims));
                for i=1:numel(this.circularDims)
                    if this.circularDims(i) > numel(s.subs)
                        break;
                    end

                    subs = s.subs{this.circularDims(i)};

                    if isnumeric(subs)
                        subs = mod(subs-1, size(this.mat, this.circularDims(i)))+1;
                    elseif ischar(subs) && strcmp(subs,':')
                        keepCircularDims(i)=true;
                    end

                    s.subs{this.circularDims(i)} = subs;
                end

                value = CircularMat(subsref(this.mat, s), this.circularDims(1,keepCircularDims));
            end
        end

        function this = subsasgn(this, s, rhs)
            import edu.stanford.covert.util.CircularMat;

            %Case 1: Dot reference to properties
            if strcmp(s(1).type,'.')
                this = builtin('subsasgn', this, s, rhs);
                return;
            end

            if numel(s.subs)==1
                %Case 2: linear indexing
                this.mat = subsasgn(this.mat, s, rhs);
            else
                %Case 3: subscript indexing
                keepCircularDims = true(size(this.circularDims));
                for i=1:numel(this.circularDims)
                    if this.circularDims(i) > numel(s.subs)
                        break;
                    end

                    subs = s.subs{this.circularDims(i)};

                    if isnumeric(subs)
                        subs = mod(subs-1, size(this.mat, this.circularDims(i)))+1;
                        if any(subs > size(this.mat, this.circularDims(i)))
                            keepCircularDims(i)=false;
                        end
                    end

                    s.subs{this.circularDims(i)} = subs;
                end

                this.mat = subsasgn(this.mat, s, rhs);
                this.circularDims=this.circularDims(1,keepCircularDims);
            end
        end

        function inds = subsindex(this)
            inds = subsindex(this.mat);
        end

        function e = end(this, k, n)
            e = feval('end',this.mat, k, n);
        end
    end

    %concatenation
    methods
        function C = cat(dim, varargin)
            import edu.stanford.covert.util.CircularMat;

            for i=1:numel(varargin)
                if isa(varargin{i}, 'CircularMat');
                    C = varargin{i};
                    break;
                end
            end

            values = cell(size(varargin));
            for i=1:numel(varargin)
                if isa(varargin{i}, 'CircularMat');
                    values{i}=varargin{i}.mat;
                else
                    values{i}=varargin{i};
                end
            end

            C.mat = cat(dim, values{:});
            if ~isempty(C.circularDims) && numel(values)>1
                C.circularDims(:, C.circularDims==dim)=[];
            end
        end

        function C = horzcat(varargin)
            C = cat(2, varargin{:});
        end
        function C = vertcat(varargin)
            C = cat(1, varargin{:});
        end

        function this = padarray(this, varargin)
            oldsiz = size(this.mat);
            this.mat = padarray(this.mat, varargin{:});
            newsiz = size(this.mat);
            this.circularDims = intersect(this.circularDims, find(oldsiz==newsiz(1:numel(oldsiz))));
        end
    end

    %element-wise functions
    methods
        function this = abs(this)
            this.mat=abs(this.mat);
        end

        function this = sign(this)
            this.mat = sign(this.mat);
        end

        function this = sqrt(this)
            this.mat = sqrt(this.mat);
        end

        function this = realsqrt(this)
            this.mat = realsqrt(this.mat);
        end

        function this = exp(this)
            this.mat = exp(this.mat);
        end

        function this = expm1(this)
            this.mat = expm1(this.mat);
        end

        function this = log(this)
            this.mat = log(this.mat);
        end

        function this = log1p(this)
            this.mat = log1p(this.mat);
        end

        function this = log2(this)
            this.mat = log2(this.mat);
        end

        function this = log10(this)
            this.mat = log10(this.mat);
        end

        function this = reallog(this)
            this.mat = reallog(this.mat);
        end
    end

    %complex numbers
    methods
        function this = real(this)
            this.mat = real(this.mat);
        end

        function this = imag(this)
            this.mat = imag(this.mat);
        end

        function this = angle(this)
            this.mat = angle(this.mat);
        end

        function this = conj(this)
            this.mat = conj(this.mat);
        end

        function val = isreal(this)
            val = isreal(this.mat);
        end
    end

    %logic
    methods
        function this = not(this)
            this.mat = not(this.mat);
        end

        function C = or(A, B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = or(A.mat, B.mat);
                else
                    C.mat = or(A.mat, B);
                end
            else
                C = B;
                C.mat = or(A, B.mat);
            end
        end

        function C = and(A, B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = and(A.mat, B.mat);
                else
                    C.mat = and(A.mat, B);
                end
            else
                C = B;
                C.mat = and(A, B.mat);
            end
        end

        function C = xor(A, B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = xor(A.mat, B.mat);
                else
                    C.mat = xor(A.mat, B);
                end
            else
                C = B;
                C.mat = xor(A, B.mat);
            end
        end
    end

    %opertors
    methods
        function C = uplus(A)
            C = A;
        end

        function C = uminus(A)
            C = A;
            C.mat = -C.mat;
        end

        function C = plus(A, B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = plus(A.mat, B.mat);
                else
                    C.mat = plus(A.mat, B);
                end
            else
                C = B;
                C.mat = plus(A, B.mat);
            end
        end

        function C = minus(A, B)
            C = A + -B;
        end

        function C = times(A, B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = times(A.mat, B.mat);
                else
                    C.mat = times(A.mat, B);
                end
            else
                C = B;
                C.mat = times(A, B.mat);
            end
        end

        function C = mtimes(A, B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = mtimes(A.mat, B.mat);
                else
                    C.mat = mtimes(A.mat, B);
                end

                if ~any(size(B))
                    C.circularDims = zeros(1,0);
                elseif numel(B)>1
                    C.circularDims(C.circularDims==2)=[];
                end
            else
                C = B;
                C.mat = mtimes(A, B.mat);
                C.circularDims(C.circularDims==1)=[];
                C.circularDims = C.circularDims + ndims(A) - 1;

                if ~any(size(A))
                    C.circularDims = zeros(1,0);
                elseif numel(A)>1
                    C.circularDims(C.circularDims==1)=[];
                    C.circularDims = C.circularDims + ndims(A)-1;
                end
            end
        end

        function C = rdivide(A,B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = rdivide(A.mat, B.mat);
                else
                    C.mat = rdivide(A.mat, B);
                end
            else
                C = B;
                C.mat = rdivide(A, B.mat);
            end
        end

        function C = ldivide(A,B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = ldivide(A.mat, B.mat);
                else
                    C.mat = ldivide(A.mat, B);
                end
            else
                C = B;
                C.mat = ldivide(A, B.mat);
            end
        end

        %function C = mrdivide(A, B)
        %end

        %function C = mldivide(A, B)
        %end

        function C = power(A,B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = power(A.mat, B.mat);
                else
                    C.mat = power(A.mat, B);
                end
            else
                C = B;
                C.mat = power(A, B.mat);
            end
        end

        function C = mpower(A,B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = mpower(A.mat, B.mat);
                else
                    C.mat = mpower(A.mat, B);
                end
            else
                C = B;
                C.mat = mpower(A, B.mat);
            end
        end
    end


    %comparison
    methods
        function C = isequal(A,B)
            C = strcmp(class(A),class(B)) && ...
                isequal(A.circularDims, B.circularDims) && ...
                isequal(A.mat, B.mat);
        end

        function assertElementsAlmostEqual(A, B, varargin)
            import edu.stanford.covert.util.CircularMat;

            assertEqual(size(A), size(B));

            if isa(A, 'CircularMat')
                if isa(B, 'CircularMat')
                    assertElementsAlmostEqual(A.mat, B.mat, varargin{:});
                else
                    assertElementsAlmostEqual(A.mat, B, varargin{:});
                end
            else
                assertElementsAlmostEqual(A, B.mat, varargin{:});
            end
        end

        function C = eq(A, B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = eq(A.mat, B.mat);
                else
                    C.mat = eq(A.mat, B);
                end
            else
                C = B;
                C.mat = eq(A, B.mat);
            end
        end

        function C = ne(A, B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = ne(A.mat, B.mat);
                else
                    C.mat = ne(A.mat, B);
                end
            else
                C = B;
                C.mat = ne(A, B.mat);
            end
        end

        function C = le(A, B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = le(A.mat, B.mat);
                else
                    C.mat = le(A.mat, B);
                end
            else
                C = B;
                C.mat = le(A, B.mat);
            end
        end

        function C = lt(A, B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = lt(A.mat, B.mat);
                else
                    C.mat = lt(A.mat, B);
                end
            else
                C = B;
                C.mat = lt(A, B.mat);
            end
        end

        function C = ge(A, B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = ge(A.mat, B.mat);
                else
                    C.mat = ge(A.mat, B);
                end
            else
                C = B;
                C.mat = ge(A, B.mat);
            end
        end

        function C = gt(A, B)
            import edu.stanford.covert.util.CircularMat;

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = gt(A.mat, B.mat);
                else
                    C.mat = gt(A.mat, B);
                end
            else
                C = B;
                C.mat = gt(A, B.mat);
            end
        end
    end

    %set operators
    methods
        function this = unique(this)
            this.mat = unique(this.mat);
            this.circularDims = zeros(1,0);
        end

        %function C = union(A, B, varargin); end
        %function C = intersect(A, B, varargin); end
        %function C = setdiff((A, B, varargin); end
        %function C = setxor(A, B, varargin); end
        %function C = ismember(A, B, varargin); end
    end

    %reduction
    methods
        function this = all(this, varargin)
            oldsiz = size(this.mat);
            this.mat = all(this.mat, varargin{:});
            newsiz = size(this.mat);
            this.circularDims = intersect(this.circularDims, find(oldsiz(1:numel(newsiz))==newsiz));
        end

        function this = any(this, varargin)
            oldsiz = size(this.mat);
            this.mat = any(this.mat, varargin{:});
            newsiz = size(this.mat);
            this.circularDims = intersect(this.circularDims, find(oldsiz(1:numel(newsiz))==newsiz));
        end

        function C = min(A, B, dim)
            import edu.stanford.covert.util.CircularMat;

            if ~exist('B','var')
                C = A;
                C.mat = min(C.mat,[],1);
                C.circularDims(1,C.circularDims==1)=[];
                return
            end

            if isempty(B)
                C = A;
                C.mat = min(C.mat,[],dim);
                C.circularDims(1,C.circularDims==dim)=[];
                return;
            end

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = min(A.mat, B.mat);
                else
                    C.mat = min(A.mat, B);
                end
            else
                C = B;
                C.mat = min(A, B.mat);
            end
        end

        function C = max(A, B, dim)
            import edu.stanford.covert.util.CircularMat;

            if ~exist('B','var')
                C = A;
                C.mat = max(C.mat,[],1);
                C.circularDims(1,C.circularDims==1)=[];
                return
            end

            if isempty(B)
                C = A;
                C.mat = max(C.mat,[],dim);
                C.circularDims(1,C.circularDims==dim)=[];
                return;
            end

            if isa(A, 'CircularMat')
                C = A;
                if isa(B, 'CircularMat')
                    C.mat = max(A.mat, B.mat);
                else
                    C.mat = max(A.mat, B);
                end
            else
                C = B;
                C.mat = max(A, B.mat);
            end
        end

        function this = range(this, varargin)
            oldsiz = size(this.mat);
            this.mat = range(this.mat, varargin{:});
            newsiz = size(this.mat);
            this.circularDims = intersect(this.circularDims, find(oldsiz(1:numel(newsiz))==newsiz));
        end

        function this = sum(this, varargin)
            oldsiz = size(this.mat);
            this.mat = sum(this.mat, varargin{:});
            newsiz = size(this.mat);
            this.circularDims = intersect(this.circularDims, find(oldsiz(1:numel(newsiz))==newsiz));
        end

        function this = mean(this, varargin)
            oldsiz = size(this.mat);
            this.mat = mean(this.mat, varargin{:});
            newsiz = size(this.mat);
            this.circularDims = intersect(this.circularDims, find(oldsiz(1:numel(newsiz))==newsiz));
        end

        function this = var(this, varargin)
            oldsiz = size(this.mat);
            this.mat = var(this.mat, varargin{:});
            newsiz = size(this.mat);
            this.circularDims = intersect(this.circularDims, find(oldsiz(1:numel(newsiz))==newsiz));
        end

        function this = std(this, varargin)
            oldsiz = size(this.mat);
            this.mat = std(this.mat, varargin{:});
            newsiz = size(this.mat);
            this.circularDims = intersect(this.circularDims, find(oldsiz(1:numel(newsiz))==newsiz));
        end

        function val = norm(this, varargin)
            val = norm(this.mat, varargin{:});
        end
    end

    %find
    methods
        function varargout = find(this, varargin)
            switch nargout
                case 1
                    varargout{1} = find(this.mat, varargin{:});
                case 2
                    [varargout{1},varargout{2}] = find(this.mat, varargin{:});
                case 3
                    [varargout{1},varargout{2},varargout{3}] = find(this.mat, varargin{:});
                otherwise
                    [varargout{:}] = find(this.mat, varargin{:});
            end
        end
    end

    %iss
    methods
        function this = isnan(this)
            this.mat = isnan(this.mat);
        end

        function this = isinf(this)
            this.mat = isinf(this.mat);
        end

        function this = isfinite(this)
            this.mat = isfinite(this.mat);
        end
    end

    methods
        function val = matClass(this)
            val = class(this.mat);
        end

        function this = matCast(this, newclass)
            this.mat = cast(this.mat, newclass);
        end

        function val = cast(this, newclass)
            val = this.(newclass)();
        end

        function val = sparse(this)
            val = sparse(this.mat);
        end

        function val = char(this)
            val = cast(this.mat, 'char');
        end

        function val = logical(this)
            val = cast(this.mat, 'logical');
        end

        function val = int8(this)
            val = cast(this.mat, 'int8');
        end

        function val = int16(this)
            val = cast(this.mat, 'int16');
        end

        function val = int32(this)
            val = cast(this.mat, 'int32');
        end

        function val = int64(this)
            val = cast(this.mat, 'int64');
        end

        function val = uint8(this)
            val = cast(this.mat, 'uint8');
        end

        function val = uint16(this)
            val = cast(this.mat, 'uint16');
        end

        function val = uint32(this)
            val = cast(this.mat, 'uint32');
        end

        function val = uint64(this)
            val = cast(this.mat, 'uint64');
        end

        function val = single(this)
            val = cast(this.mat, 'single');
        end

        function val = double(this)
            val = cast(this.mat, 'double');
        end
    end

    %display
    methods
        function disp(this)
            display(this);
        end

        function display(this)
            siz = size(this);
            fprintf('%d',siz(1));
            fprintf('x%d',siz(2:end));
            fprintf(' CircularMat of type %s with ',class(this.mat));
            if isempty(this.circularDims)
                fprintf('no circular dimensions\n');
            else
                 fprintf('circular dimensions {%d',this.circularDims(1));
                 if numel(this.circularDims)>1
                     fprintf(', %d',this.circularDims(2:end));
                 end
                 fprintf('}\n',this.circularDims(2:end));
            end
            display(this.mat);
        end
    end
end