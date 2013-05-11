% SparseMat
%   Multidimensional sparse matrix with support for
%   - subscript reference/assignment
%   - sizing (isempty, numel, nnz, size, ndims, nnz, reshape, repmat, squeeze,
%     transpose, ctranspose, permute)
%   - concatenation
%   - logical operators
%   - relational operators
%   - element-wise algebra
%   - reduction (eg. all, any, min, max, range, sum, mean, var, std, collapse,
%     norm)
%   - find
%   - isinf, isfinite, isnan
%   - type casting
%   - abs, sign, sqrt, realsqrt, exp, expm1, log, logp1, log2, log10, reallog
%   - complex numbers (ctranspose, conf, real, imag, isreal, angle)
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/18/2010
classdef SparseMat
    properties (SetAccess = protected)
        subs
        vals
        siz
    end
    
    properties (Dependent = true, SetAccess = protected)
        inds
    end
    
    %constructor
    methods
        %three ways to instatiate:
        %(1) x = SparseMat()
        %    Creates an empty sparse matrix of size 0x2
        %
        %(2) x = SparseMat(mat)
        %    Casts matrix mat to sparse matrix
        %
        %(3) x = SparseMat(subs, vals, sz)
        %    Creates sparse matrix with size sz and non-zeros with subscripts
        %    subs and values vals
        function this = SparseMat(varargin)
            import edu.stanford.covert.util.SparseMat;
            
            switch nargin
                case 3
                    subs = varargin{1};
                    vals = varargin{2};
                    siz = varargin{3};
                    
                    if ~any(size(subs)); subs = zeros(0, size(siz, 2)); end;
                    if ~any(size(vals)); vals = cast(zeros(0, 1), class(vals)); end;
                    if numel(vals) == 1; vals = vals(ones(size(subs, 1), 1), 1); end;
                    
                    %squeeze trailing singleton dimensions
                    idx = max(2, find(siz ~= 1, 1, 'last'));
                    if ~isempty(idx)
                        if ~all(subs(:, idx+1:end) == 1)
                            throw(MException('SparseMat:invalidDimensions','Index exceeds matrix dimensions.'))
                        end
                        siz = siz(:, 1:idx);
                    end
                    if size(subs, 2) > size(siz, 2) && all(all(subs(:, size(siz, 2) + 1:end) == 1))
                        subs = subs(:, 1:size(siz,2));
                    end
                    
                    if size(subs, 2) ~= size(siz, 2)
                        throw(MException('SparseMat:error','dimensions of subs must equal length of siz'));
                    elseif size(subs, 1) ~= size(vals, 1)
                        throw(MException('SparseMat:error','numbers of subs and vals must be equal'));
                    elseif size(vals, 2) > 1
                        throw(MException('SparseMat:error','length of second dimension of vals can be at most 1'));
                    elseif any(siz < 0 | siz ~= ceil(siz))
                        throw(MException('SparseMat:error','siz must be non-negative integers'));
                    elseif any(subs < 1 | subs ~= ceil(subs))
                        throw(MException('SparseMat:error','subs must be positive integers'));
                    elseif ~isempty(subs) && any(max(subs, [], 1) > siz)
                        throw(MException('SparseMat:error','subs cannot be greater than lengths of dimensions'));
                    end
                    
                    if ~this.isunique_subs(subs, siz)
                        [~, idxs] = unique(subs, 'rows', 'last');
                        warning('SparseMat:invalidAssignment', 'subscripts are repeated, setting element to last value')
                        subs = subs(idxs, :);
                        vals = vals(idxs, :);
                    end
                    
                    this.siz  = siz;
                    this.subs = subs;
                    this.vals = vals;
                    
                    this = this.normalize(true, true, true, zeros(1, this.ndims));
                case 0
                    this.siz = [0 2];
                    this.subs = zeros(0,2);
                    this.vals = zeros(0,1);
                case 1
                    arg = varargin{1};
                    if isstruct(arg)
                        this.subs = arg.subs;
                        this.vals = arg.vals;
                        this.siz = arg.siz;
                    else
                        this.siz = size(arg);
                        
                        if isa(arg,'SparseMat')
                            [this.subs, this.vals] = find(arg);
                        else
                            [inds, ~, vals] = find(arg(:)); %#ok<*PROP>
                            this.subs = this.ind2sub(inds);
                            this.vals = cast(vals, class(arg));
                        end
                        
                        this = this.normalize(true, false, false, zeros(1, this.ndims));
                    end
                otherwise
                    throw(MException('SparseMat:error','no constructor matches the calling signature'));
            end
        end
    end
    
    %getters
    methods
        function inds = get.inds(this)
            inds = this.sub2ind(this.subs);
        end
    end
    
    %size
    methods
        function val = isempty(this)
            val = ~all(this.siz);
        end
        
        function n = numel(this)
            n = prod(this.siz);
        end
        
        function n = nnz(this)
            n = size(this.subs, 1);
        end
        
        function l = length(this)
            l = max(this.siz);
        end
        
        function siz = size(this, k)
            if ~exist('k','var')
                siz = this.siz;
            elseif k < 1 || k~=ceil(k)
                throw(MException('SparseMat:invalidDimensions','k must be a positive integer'))
            elseif k > size(this.siz,2)
                siz = 1;
            else
                siz = this.siz(k);
            end
        end
        
        function n = ndims(this)
            n = numel(this.siz);
        end
        
        function this = reshape(this, varargin)
            if numel(varargin)==0
                throw(MException('SparseMat:invalidDimensions','Not enough input arguments.'));
            elseif numel(varargin)==1
                dims=varargin{1};
            else
                idxs = cellfun(@isempty,varargin);
                if sum(idxs)>1
                    throw(MException('SparseMat:invalidDimensions','Size can only have one unknown dimension.'));
                end
                dims = zeros(1, numel(varargin));
                dims(~idxs) = cell2mat(varargin(~idxs));
                dims(idxs) = max(0, this.numel/prod(dims(~idxs)));
            end
            
            if prod(dims) ~= this.numel
                throw(MException('SparseMat:invalidDimensions','the number of elements must not change'));
            elseif any(dims<0) || ~all(isint(dims))
                throw(MException('SparseMat:invalidDimensions','dims must be non-negative integers'))
            end
            
            inds = this.inds;
            this.siz = dims;
            this.subs = this.ind2sub(inds);
            this = this.normalize(true, false, false, zeros(1,this.ndims));
        end
        
        function this = repmat(this, varargin)
            if numel(varargin) == 0
                throw(MException('SparseMat:invalidDimensions','Not enough input arguments.'));
            end
            copies = cell2mat(varargin);
            if any(copies<0) || ~all(isint(copies))
                throw(MException('SparseMat:invalidDimensions','copies must be non-negative integers'))
            end
            copies = copies(1:max(2,find(copies~=1, 1, 'last')));
            len = length(copies);
            if all(copies)
                if any(copies>1)
                    subs = this.subs;
                    siz = this.siz;
                    subs(:,end+1:len) = 1;
                    siz(end+1:len) = 1;
                    for i = 1:len
                        n = copies(i);
                        if n > 1
                            nVals = size(subs,1);
                            subs = repmat(subs, n, 1);
                            subs(nVals+1:end,i) = subs(nVals+1:end,i) + ...
                                reshape(repmat(siz(i) * (1:n-1), [nVals 1]), [], 1);
                        end
                    end
                    this.subs = subs;
                    this.vals = repmat(this.vals, prod(copies), 1);
                    this.siz(1:len) = siz(1:len) .* copies;
                    oldDims = 1 : this.ndims;
                    oldDims(copies ~= 1) = 0;
                    this = this.normalize(false, false, true, oldDims);
                end
            else
                this.subs = [];
                this.vals = [];
                this.siz(end+1:len) = 1;
                this.siz(1:len) = this.siz(1:len) .* copies;
            end
        end
        
        function this = transpose(this)
            if this.ndims > 2
                throw(MException('SparseMat:invalidDimensions','Transpose on ND array is not defined'));
            end
            this = this.permute([2 1]);
        end
        
        function this = ctranspose(this)
            this = transpose(conj(this));
        end
        
        function this = permute(this, dims)
            if ~all(diff(sort(dims)))
                throw(MException('SparseMat:invalidDimensions','ORDER cannot contain repeated permutation indices'))
            elseif any(dims < 1) || ~isequal(dims, ceil(dims))
                throw(MException('SparseMat:invalidDimensions','dims must be positive integers'))
            elseif sum(ismembc(dims, 1:this.ndims)) ~= this.ndims
                throw(MException('SparseMat:invalidDimensions','dims must contain each dimension of the matrix'))
            end
            
            oldDims = ismembc2(dims, 1:this.ndims);
            [~, idxs] = ismember(1:this.ndims, dims);
            
            subs = ones(size(this.subs,1), numel(dims));
            subs(:, idxs) = this.subs;
            
            siz = ones(1, numel(dims));
            siz(:, idxs) = this.siz;
            
            this.subs = subs;
            this.siz = siz;
            
            this = this.normalize(true, false, true, oldDims);
        end
        
        function this = squeeze(this)
            oldDims = find(this.siz~=1);
            dims = find(this.siz==1);
            this.subs(:,dims)=[];
            this.siz(:,dims)=[];
            
            if this.ndims<2
                this.subs = [this.subs ones(size(this.subs,1), 2-size(this.subs,2))];
                this.siz = [this.siz ones(1, 2-size(this.siz,2))];
                oldDims = [oldDims 0];
            end
            
            this = this.normalize(false, false, false, oldDims);
        end
        
        function inds = sub2ind(this, subs, siz)
            if ~exist('siz','var')
                siz = this.siz;
            end
            
            mult = [1 cumprod(siz(1:end-1))];
            inds = (subs - 1) * mult' + 1;
        end
        
        function subs = ind2sub(this, inds, siz)
            if ~exist('siz','var')
                siz = this.siz;
            end
            
            n = length(siz);
            subs = zeros(numel(inds), n);
            
            if isempty(inds)
                return;
            end
            
            k = [1 cumprod(siz(1:end-1))];
            inds = inds - 1;
            
            for i = n : -1 : 1
                subs(:,i) = floor(inds / k(i)) + 1;
                inds = rem(inds, k(i));
            end
        end
    end
    
    %subscript reference, assignment
    methods
        %Supports three syntaxes
        %1. Dot assignment is referred to the super class subsasgn method.
        %2. Subscripts of elements, eg.
        %   x([r11, r21, ..., rN1; r21, r22, ..., rN2; ...])=y.
        %3. Subscripts ranges, eg. x(R1,R2,...,RN)=y where Rn are ranges.
        %
        %For Cases (2) and (3) y can either be a SparseMat or a builtin MATLAB
        %type. In these cases y will be cast to the class of x.vals (one of
        %MATLAB's builtin types).
        function this = subsasgn(this, s, rhs)
            import edu.stanford.covert.util.SparseMat;
            
            %Case 1: Dot reference to properties
            if strcmp(s.type,'.') || strcmp(s.type, '{}')
                this = builtin('subsasgn', this, s, rhs);
                return;
            end
            
            rhs = cast(rhs, valueClass(this));
            
            if numel(s.subs)==1
                %Case 2: Subscripts of elements
                subs = s.subs{1};
                
                %do nothing if subscripts are empty
                if isempty(subs)
                    if numel(rhs)~=1 && ~isequal(size(rhs), [0 1])
                        throw(MException('SparseMat:invalidAssignment','Subscripted assignment dimension mismatch.'));
                    else
                        return;
                    end
                end
                
                %check that subscripts have 2 dimensions
                if ndims(subs) ~= 2 || size(subs,2)~= this.ndims
                    throw(MException('SparseMat:invalidDimensions','Subscripts must have number of columns equal to dimensions of matrix.'));
                end
                
                %check that subscripts are positive integers
                if any(any(subs<1)) || ~isequal(subs, ceil(real(subs)))
                    throw(MException('SparseMat:invalidDimensions','Subscript indices must either be real positive integers or logicals.'));
                end
                
                %check that subscripts are in range
                if any(max(subs,[],1) > this.siz)
                    throw(MException('SparseMat:invalidDimensions','Index exceeds matrix dimensions.'));
                end
                
                %check that rhs isn't empty
                if isempty(rhs)
                    throw(MException('SparseMat:invalidAssignment','Subscripted assignment dimension mismatch.'));
                end
                
                %check that rhs is vector
                if numel(rhs) > 1 && sum(size(rhs)>1)>1
                    throw(MException('SparseMat:invalidAssignment','Subscripted assignment dimension mismatch.'));
                end
                rhs = rhs(:);
                
                %check that rhs size matches that of subscripts
                if numel(rhs) > 1 && numel(rhs)~=size(subs,1)
                    throw(MException('SparseMat:invalidAssignment','Subscripted assignment dimension mismatch.'));
                end
                
                %remove duplicates
                if this.isunique_subs(subs, this.siz)
                    idxs = (1:size(subs, 1))';
                else
                    [subs, idxs] = unique(subs, 'rows', 'last');
                    warning('SparseMat:invalidAssignment','subscripts are repeated, setting element to last value')
                end
                
                if numel(rhs) > 1
                    sref = struct('type', '()', 'subs', []);
                    sref.subs = {idxs, 1};
                    rhs = subsref(rhs, sref);
                end
                
                %assign
                [tfs, idxs] = this.ismember_subs(subs, this.subs, this.siz);
                if numel(rhs) == 1
                    if isa(rhs,'SparseMat')
                        rhs = full(rhs);
                    end
                    
                    if rhs == 0
                        this.subs(idxs(tfs),:) = [];
                        this.vals(idxs(tfs),:) = [];
                    else
                        this.vals(idxs(tfs),:) = rhs;
                        
                        this.subs = [this.subs; subs(~tfs,:)];
                        this.vals = [this.vals; rhs(ones(size(this.subs,1)-size(this.vals,1), 1), 1)];
                        this = this.normalize(false, false, true, 1:this.ndims);
                    end
                elseif isa(rhs, 'SparseMat')
                    sref = struct('type', '()', 'subs', []);
                    sref.subs = {find(tfs), 1};
                    this.vals(idxs(tfs),:) = subsref(rhs, sref);
                    
                    [subIdxs, rhsIdxs] = intersect(rhs.subs(:,1), find(~tfs));
                    
                    this.subs = [this.subs; subs(subIdxs, :)];
                    this.vals = [this.vals; rhs.vals(rhsIdxs, :)];
                    this = this.normalize(false, true, true, 1:this.ndims);
                else
                    this.vals(idxs(tfs),:) = rhs(tfs,:);
                    
                    idxs = find(~tfs & rhs);
                    this.subs = [this.subs; subs(idxs, :)];
                    this.vals = [this.vals; rhs(idxs, :)];
                    this = this.normalize(false, true, true, 1:this.ndims);
                end
            else
                %Case 3: Subscripts ranges
                subs = s.subs;
                
                %Do all short circuting and error checking before extraction to
                %avoid doing the work of extraction if isn't going to be
                %necessary
                %1. check that subcripts of dimensions beyond that of the matrix
                %   are either : or 1.
                %2. check if any subs{n} is empty
                %3. check that subscripts are vectors, are within range, and are
                %   real positive integers
                %4. if fewer subscripts provided than dimensions of matrix, assume
                %   all missing subscripts are 1
                
                %compute size that RHS must be
                lhsSiz = zeros(1, numel(subs));
                for n = 1:numel(subs)
                    if ischar(subs{n}) && (subs{n} == ':')
                        if n > this.ndims
                            lhsSiz(n) = -1;
                        else
                            lhsSiz(n) = this.siz(n);
                        end
                    else
                        lhsSiz(n) = numel(subs{n});
                    end
                end
                
                %check RHS size matches that of subscripts
                rhsSiz = size(rhs);
                lhsDims = find(lhsSiz ~= 1);
                rhsDims = find(rhsSiz ~= 1);
                if ~any(rhsSiz)
                    throw(MException('SparseMat:invalidAssignment', 'RHS cannot be null.'));
                elseif numel(rhs) ~= 1
                    if numel(lhsDims) ~= numel(rhsDims)
                        throw(MException('SparseMat:invalidAssignment', 'Subscripted assignment dimension mismatch.'));
                    end
                    
                    for n = 1:numel(lhsDims)
                        if lhsSiz(lhsDims(n)) == -1
                            if lhsDims(n) > this.ndms  && rhsSiz(rhsDims(b)) > 1
                                throw(MException('SparseMat:invalidAssignment', 'Subscripted assignment dimension mismatch.'));
                            end
                            subs{lhsDims(n)} = (1:rhsSiz(rhsDims(n)))';
                            lhsSiz(lhsDims(n)) = rhsSize(rhsDims(n));
                        end
                    end
                    
                    if ~isequal(lhsSiz(lhsSiz ~= 1), rhsSiz(rhsSiz ~= 1))
                        throw(MException('SparseMat:invalidAssignment', 'Subscripted assignment dimension mismatch.'));
                    end
                else
                    for i = 1:numel(lhsSiz)
                        if lhsSiz(i) == -1;
                            subs{i} = 1;
                            lhsSiz(i) = 1;
                        end
                    end
                end
                
                if ~all(lhsSiz)
                    return;
                elseif numel(rhs) > 1
                    rhs = reshape(rhs, lhsSiz);
                end
                
                % reshape subscripts to vectors, and check subscripts are real
                % positive integers
                for n = 1:numel(subs)
                    if ischar(subs{n}) && subs{n} == ':'
                        continue;
                    end
                    
                    %reshape subs{n} to vector
                    if ~(size(subs{n}, 2) == 1 && ndims(subs{n}) == 2)
                        subs{n} = subs{n}(:);
                    end
                    
                    %check that subscripts are positive integers
                    if any(subs{n} < 1) || ~isequal(subs{n}, ceil(real(subs{n})))
                        throw(MException('SparseMat:invalidDimensions', 'Subscript indices must either be real positive integers or logicals.'));
                    end
                end
                
                %remove duplicates
                for i = 1:numel(lhsDims)
                    n = lhsDims(i);
                    
                    if ischar(subs{n}) && subs{n} == ':'
                        continue;
                    end
                    
                    if isscalar(subs{n}) || edu.stanford.covert.util.SetUtil.isunique(subs{n})
                        oldIdxs = (1:numel(subs{n}))';
                    else
                        [subs{n}, oldIdxs] = unique(subs{n}, 'last');
                        warning('SparseMat:invalidAssignment', 'subscripts are repeated, setting element to last value');
                    end
                    
                    lhsSiz(n) = numel(subs{n});
                    
                    if numel(rhs) > 1
                        m = rhsDims(i);
                        
                        if isa(rhs, 'SparseMat')
                            newIdxs = (1:size(subs{n}, 1))';
                            [tfs, idxs] = ismember(rhs.subs(:, m), oldIdxs);
                            rhs.subs(tfs, m) = newIdxs(idxs(tfs));
                            rhs.subs(~tfs, :) = [];
                            rhs.vals(~tfs, :) = [];
                        else
                            sref = struct('type', '()', 'subs', []);
                            sref.subs = num2cell(repmat(':', 1, ndims(rhs)));
                            sref.subs{m} = oldIdxs;
                            rhs = subsref(rhs, sref);
                        end
                    end
                end
                
                % if fewer subscripts provided than dimensions of matrix, assume
                % all missing subscripts are 1
                for n = numel(subs) + 1:this.ndims
                    subs{n} = ones(1, 1);
                    lhsSiz(n) = 1;
                end
                
                %compute new size
                siz = [this.siz ones(1, numel(subs) - size(this.siz, 2))];
                for n = 1:numel(subs)
                    if ischar(subs{n}) && subs{n} == ':'
                        siz(n) = lhsSiz(n);
                    else
                        siz(n) = max(siz(n), max(subs{n}));
                    end
                end
                
                if numel(rhs) == 1
                    siz(siz == -1) = 1;
                else
                    ansDims = find(siz ~= 1);
                    rhsDims = find(rhsSiz ~= 1);
                    siz(ansDims(siz(ansDims) == -1)) = rhsSiz(rhsDims(siz(ansDims) == -1));
                end
                
                % Delete what currently occupies the specified range
                deleteRange = true;
                for n = this.ndims + 1:numel(subs)
                    if isnumeric(subs{n}) && ~any(subs{n} == 1)
                        deleteRange = false;
                    end
                end
                
                if deleteRange
                    delIdxs = 1:this.nnz;
                    for n = 1:this.ndims
                        if ischar(subs{n}) && subs{n} == ':'
                            continue;
                        end
                        
                        delIdxs = delIdxs(ismembc(this.subs(delIdxs, n), sort(subs{n})));
                    end
                    
                    this.subs(delIdxs, :) = [];
                    this.vals(delIdxs, :) = [];
                end
                
                %assign
                orginalDims = [1:this.ndims zeros(1, this.ndims - size(siz, 2))];
                orginalDims(this.siz ~= siz(1:this.ndims)) = 0;
                
                if numel(rhs) == 1
                    if isa(rhs, 'SparseMat')
                        rhs = full(rhs);
                    end
                    
                    if rhs~=0
                        newsubs = zeros(prod(lhsSiz), size(lhsSiz, 2));
                        for n = 1:numel(subs)
                            subsSize = lhsSiz;
                            subsSize(n) = 1;
                            if n == 1
                                subsDims = [1 2];
                            else
                                subsDims = n:-1:1;
                            end
                            if ischar(subs{n}) && subs{n} == ':'
                                newsubs(:, n) = reshape(repmat(permute((1:siz(n))', subsDims), subsSize), [], 1);
                            else
                                newsubs(:, n) = reshape(repmat(permute(subs{n}, subsDims), subsSize), [], 1);
                            end
                        end
                        
                        this.siz = siz;
                        this.subs = [this.subs ones(size(this.subs, 1), size(lhsSiz, 2) - size(this.subs, 2));
                            newsubs];
                        this.vals = [this.vals;
                            rhs(ones(size(newsubs, 1), 1), 1)];
                    end
                else
                    if ~isa(rhs,'SparseMat')
                        rhs = SparseMat(rhs);
                    end
                    
                    ndim = max([size(rhs.subs, 2) numel(subs) numel(lhsSiz) numel(rhsSiz) numel(siz)]);
                    
                    rhs.subs = [rhs.subs ones(size(rhs.subs, 1), ndim-size(rhs.subs, 2))];
                    rhsSubs = [rhs.subs ones(size(rhs.subs, 1), ndim-size(rhs.subs, 2))];
                    subs(numel(lhsSiz) + 1:ndim) = {1};
                    lhsSiz = [lhsSiz ones(1, ndim - numel(lhsSiz))];
                    rhsSiz = [rhsSiz ones(1, ndim - numel(rhsSiz))];
                    siz = [siz ones(1, ndim - numel(siz))];
                    lhsDims = [find(lhsSiz ~= 1) find(lhsSiz == 1)];
                    rhsDims = [find(rhsSiz ~= 1) find(rhsSiz == 1)];
                    for i = 1:numel(lhsDims)
                        n = lhsDims(i);
                        m = rhsDims(i);
                        if isnumeric(subs{n})
                            rhsSubs(:,m) = subs{n}(rhs.subs(:, m));
                        else
                            tmpsubs = (1:siz(n))';
                            rhsSubs(:,m) = tmpsubs(rhs.subs(:, m));
                        end
                    end
                    
                    this.siz = siz;
                    this.subs = [this.subs ones(size(this.subs, 1), size(rhsSubs, 2) - size(this.subs, 2));
                        rhsSubs ones(size(rhsSubs, 1), size(this.subs, 2)-size(rhsSubs, 2))];
                    this.vals = [this.vals;
                        rhs.vals];
                end
                
                this = this.normalize(true, false, true, orginalDims);
            end
        end
        
        %Supports three syntaxes
        %1. Dot reference is referred to the super class subsref method.
        %2. Subscripts of elements, eg.
        %   y=x([r11, r21, ..., rN1; r21, r22, ..., rN2; ...]). This returns a
        %   matrix of type equal to that of the type of the vals property (that
        %   is a matrix of one of MATLAB's builtin types).
        %3. Subscripts ranges, eg. y=x(R1,R2,...,RN) where Rn are ranges. This
        %   returns an object of SparseMat type
        function value = subsref(this, s)
            import edu.stanford.covert.util.SparseMat;
            
            %Case 1: Dot reference to properties
            if strcmp(s.type,'.') || strcmp(s.type, '{}')
                value = builtin('subsref', this, s);
                return;
            end
            
            if numel(s.subs) == 1
                %Case 2: Subscripts of elements
                subs = s.subs{1};
                
                %check that subscripts have 2 dimensions
                if ndims(subs) ~= 2 || size(subs,2) ~= this.ndims
                    throw(MException('SparseMat:invalidDimensions','Subscripts must have 2 dimensions.'));
                end
                
                %check that subscripts are real positive integers
                validateattributes(subs, {'numeric'}, {'real', 'positive', 'integer'});
                
                %check that subscripts are in range
                if ~isempty(subs) && any(max(subs,[],1) > this.siz)
                    throw(MException('SparseMat:invalidDimensions','Index exceeds matrix dimensions.'));
                end
                
                %extract
                [tfs, idxs] = this.ismember_subs(subs, this.subs, this.siz);
                switch class(this.vals)
                    case 'logical'
                        value = false(size(subs,1),1);
                    case 'char'
                        value = char(zeros(size(subs,1),1));
                    otherwise
                        value = cast(zeros(size(subs,1),1), class(this.vals));
                end
                
                value(tfs) = this.vals(idxs(tfs),:);
            else
                %Case 3: Subscript ranges
                subs = s.subs;
                
                %Do all short circuting and error checking before extraction to
                %avoid doing the work of extraction if isn't going to be
                %necessary
                %1. check that subcripts of dimensions beyond that of the matrix
                %   are either : or 1.
                %2. check if any subs{n} is empty
                %3. check that subscripts are vectors, are within range, and are
                %   real positive integers
                %4. if fewer subscripts provided than dimensions of matrix, assume
                %   all missing subscripts are 1
                
                %check that subcripts of dimensions beyond that of the matrix
                %are either : or 1.
                for n = this.ndims+1:numel(subs)
                    if ~((ischar(subs{n}) && (subs{n} == ':')) || (isnumeric(subs{n}) && all(subs{n}(:)==1)))
                        throw(MException('SparseMat:invalidDimensions','Index exceeds matrix dimensions.'));
                    end
                end
                subs = subs(1:min(numel(subs),this.ndims));
                
                %compute size of return value
                originalDims = zeros(1, this.ndims);
                newsiz = zeros(1, numel(subs));
                for n = 1:numel(subs)
                    if ischar(subs{n}) && subs{n} == ':'
                        newsiz(n) = this.siz(n);
                        originalDims(n) = n;
                    else
                        newsiz(n) = numel(subs{n});
                    end
                end
                
                %check if any subs{n} is empty
                if ~all(newsiz)
                    value = this;
                    value.subs = zeros(0, this.ndims);
                    value.vals = cast(zeros(0, 1), class(this.vals));
                    value.siz = [newsiz ones(1, this.ndims - numel(newsiz))];
                    value = value.normalize(true, false, false, [originalDims zeros(1, this.ndims - size(originalDims,2))]);
                    return;
                end
                
                % check that subscripts are vectors, are within range, and are
                % real positive integers
                for n = 1:numel(subs)
                    if ischar(subs{n}) && subs{n} == ':'
                        continue;
                    end
                    
                    %check that subs{n} is vector
                    validateattributes(subs, {'cell'}, {'vector'});
                    subs{n} = subs{n}(:);
                    
                    %check that subscripts are real positive integers in range
                    validateattributes(subs{n}, {'numeric'}, {'real', 'positive', 'integer', '<=', this.siz(n)});
                end
                
                % if fewer subscripts provided than dimensions of matrix, assume
                % all missing subscripts are 1
                for n = numel(subs) + 1:this.ndims
                    subs{n} = ones(1, 1);
                    newsiz(n) = 1;
                end
                
                %extract in order of new lengths of dimensions
                [~, order] = sort(newsiz);
                newsubs = this.subs;
                newvals = this.vals;
                for i = 1:numel(subs)
                    n = order(i);
                    
                    if ischar(subs{n}) && subs{n} == ':'
                        continue;
                    elseif isscalar(subs{n})
                        tfs = subs{n} == newsubs(:, n);
                        newsubs = [newsubs(tfs, 1:n-1) ones(sum(tfs), 1) newsubs(tfs, n+1:end)];
                        newvals = newvals(tfs, :);
                    elseif all(diff(subs{n}(:)))
                        tmpsubs = zeros(0, this.ndims);
                        tmpvals = cast(zeros(0, 1), class(this.vals));
                        tmp = subs{n};
                        for j = 1:numel(tmp)
                            tfs = find(newsubs(:, n) == tmp(j));
                            tmpsubs = [tmpsubs; newsubs(tfs, 1:n-1) j(ones(numel(tfs), 1)) newsubs(tfs, n+1:end)]; %#ok<AGROW>
                            tmpvals = [tmpvals; newvals(tfs, :)]; %#ok<AGROW>
                        end
                        
                        newsubs = tmpsubs;
                        newvals = tmpvals;
                    else
                        [tmp, ~, idxs] = this.unique_subs(subs{n}, this.siz(n));
                        tmpsubs = zeros(0, this.ndims);
                        tmpvals = cast(zeros(0, 1), class(this.vals));
                        for j = 1:numel(tmp)
                            idxs1 = find(idxs == j);
                            m = size(idxs1);
                            idxs2 = find(tmp(j) == newsubs(:, n));
                            idxs1 = idxs1(:)';
                            idxs2 = idxs2(:);
                            idxs1 = idxs1(ones(size(idxs2)), :);
                            idxs2 = idxs2(:, ones(m));
                            idxs1 = idxs1(:);
                            idxs2 = idxs2(:);
                            tmpsubs = [tmpsubs; newsubs(idxs2, 1:n-1) idxs1 newsubs(idxs2, n+1:end)]; %#ok<AGROW>
                            tmpvals = [tmpvals; newvals(idxs2, :)]; %#ok<AGROW>
                        end
                        
                        newsubs = tmpsubs;
                        newvals = tmpvals;
                    end
                end
                
                value = this;
                value.subs = newsubs;
                value.vals = newvals;
                value.siz = newsiz;
                value = value.normalize(true, false, true, [originalDims zeros(1, this.ndims - size(originalDims,2))]);
            end
        end
        
        function inds = subsindex(this) %#ok<STOUT,MANU>
            throw(MException('SparseMat:unsupportedSyntax', 'Function not implemented'));
        end
        
        function e = end(this, k, ~)
            e = this.siz(k);
        end
        
        function this = normalize(this, reSize, reFind, reSort, ~) %last argument is oldDims and is used by CircularSparseMat
            if reSize
                if isempty(this.subs)
                    this.subs = zeros(0, this.ndims);
                    this.vals = cast(zeros(0, 1), class(this.vals));
                end
                
                %pad size
                if isempty(this.siz)
                    this.siz = [0 0];
                    if ~isempty(this.subs) || ~isempty(this.vals)
                        throw(MException('SparseMat:invalidDimensions','Index exceeds matrix dimensions.'))
                    end
                    this.subs = zeros(0, 0);
                    this.vals = cast(zeros(0, 1), class(this.vals));
                elseif size(this.siz, 2) == 1
                    this.siz = [this.siz 1];
                    this.subs = [this.subs ones(this.nnz, 1)];
                end
                
                %squeeze trailing singleton dimensions
                if this.siz(end) == 1
                    idx = max(2, find(this.siz ~= 1, 1, 'last'));
                    if ~isempty(idx)
                        if any(this.subs(:, idx + 1:end) ~= 1)
                            throw(MException('SparseMat:invalidDimensions', 'Index exceeds matrix dimensions.'))
                        end
                        this.siz = this.siz(1, 1:idx);
                        this.subs = this.subs(:, 1:idx);
                    end
                end
            end
            
            %remove zeros
            if reFind && ~isempty(this.vals)
                tfs = this.vals ~= 0;
                this.subs = this.subs(tfs, :);
                this.vals = this.vals(tfs, 1);
            end
            
            %sort rows
            if reSort && ~isempty(this.subs)
                [this.subs, order] = this.sort_subs(this.subs, this.siz);
                this.vals = this.vals(order, :);
            end
        end
    end
    
    %concatenation
    methods
        function this = cat(dim, varargin)
            import edu.stanford.covert.util.SparseMat;
            
            validateattributes(dim, {'numeric'}, {'scalar', 'real', 'positive', 'integer'});
            
            if nargin < 2
                return;
            end

            this = varargin{1};
            if ~isa(this, 'SparseMat')
                this = SparseMat(this);
            end
            
            nRows = 0;
            edges = zeros(numel(varargin) + 1, 1);
            for i = 1:numel(varargin)
                B = varargin{i};
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                    varargin{i} = B;
                end
                edges(i) = nRows + 1;
                nRows = nRows + size(B.vals, 1);
            end
            edges(end) = nRows + 1;
            
            siz  = this.siz;
            subs = zeros(nRows, max(size(this.subs, 2), dim));
            vals = zeros(nRows, 1);
            
            if dim > this.ndims
                siz = [siz ones(1, dim - size(siz, 2))];
            end
            siz(dim) = 0;
            
            for i = 1:numel(varargin)
                B = varargin{i};
                
                Bsiz  = B.siz;
                Bsubs = B.subs;
                Bvals = B.vals;
                
                if max(dim, size(Bsiz, 2)) > size(siz, 2)
                    siz = [siz ones(1, max(dim, size(Bsiz, 2)) - size(siz, 2))]; %#ok<AGROW>
                    subs = [subs ones(size(subs, 1), max(dim, size(Bsiz, 2)) - size(subs, 2))]; %#ok<AGROW>
                end
                
                if max(dim, size(siz, 2)) > size(Bsiz, 2)
                    Bsiz = [Bsiz ones(1, max(dim, size(siz, 2)) - size(Bsiz, 2))]; %#ok<AGROW>
                    Bsubs = [Bsubs ones(size(Bsubs, 1), max(dim, size(siz, 2)) - size(Bsubs, 2))]; %#ok<AGROW>
                end
                
                if ~isequal(siz([1:dim-1 dim+1:end]), Bsiz([1:dim-1 dim+1:end]))
                    throw(MException('SparseMat:cat', 'Inconsistent dimensions'));
                end
                
                Bsubs(:, dim) = Bsubs(:, dim) + siz(dim);
                
                subs(edges(i):edges(i+1) - 1, :) = Bsubs;
                vals(edges(i):edges(i+1) - 1, :) = Bvals;
                
                siz(dim) = siz(dim) + Bsiz(dim);
            end
            
            this.subs = subs;
            this.vals = vals;
            this.siz = siz;
            oldDims = 1:this.ndims;
            oldDims(dim) = 0;
            if dim < this.ndims
                this = this.normalize(true, true, true, oldDims);
            else
                this = this.normalize(true, false, false, oldDims);
            end
        end
        
        function this = horzcat(this, varargin)
            this = cat(2, this, varargin{:});
        end
        
        function this = vertcat(this, varargin)
            this = cat(1, this, varargin{:});
        end
        
        function this = padarray(this, padsize, padval, direction)
            import edu.stanford.covert.util.SparseMat;
            
            %padsize
            if any(padsize < 0) || ~isequal(padsize, ceil(real(padsize))) || ~all(isfinite(padsize)) || ~isvector(padsize)
                throw(MException('SparseMat:invalidDimensions','padsize must be vector of non-negative real integers'));
            end
            idx = find(padsize~=0,1,'last');
            if isempty(idx)
                return
            end
            padsize = [padsize(1:idx) zeros(1, this.ndims-idx)];
            if this.ndims < numel(padsize)
                this.siz = [this.siz ones(1, numel(padsize)-size(this.siz,2))];
                this.subs = [this.subs ones(this.nnz, numel(padsize)-size(this.subs,2))];
            end
            
            %padval
            if ~exist('padval','var')
                padval = 0;
            elseif ~isnumeric(padval) && ~(ischar(padval) && ismember(padval, {'circular','replicate','symmetric'}))
                throw(MException('SparseMat:invalidValue','padval must be numeric of one of the strings ''circular'', ''replicate'', or ''symmetric'''));
            end
            
            %direction
            if ~exist('direction','var')
                direction='both';
            end
            
            %pad
            switch direction
                case 'both'
                    if isnumeric(padval) && padval==0
                        this.siz = this.siz + 2*padsize;
                        this.subs = this.subs + padsize(ones(this.nnz, 1), :);
                    elseif isnumeric(padval)
                        for i=1:numel(padsize)
                            tmpsiz = this.siz;
                            tmpsiz(i)=padsize(i);
                            tmp = SparseMat([],[],tmpsiz);
                            tmp.subs = tmp.ind2sub(1:tmp.numel);
                            tmp.vals = padval(ones(tmp.numel, 1), 1);
                            
                            this = cat(i, tmp, this, tmp);
                        end
                    else
                        for i = 1:numel(padsize)
                            sref = struct('type', '()', 'subs', []);
                            tmp  = {':'};
                            sref.subs = tmp(1, ones(ndims(this), 1));
                            preref=sref;
                            postref=sref;
                            
                            switch padval
                                case 'circular'
                                    preref.subs{i} = mod(this.siz(i)-padsize(i)+(1:padsize(i))-1,this.siz(i))+1;
                                    postref.subs{i} = mod((1:padsize(i))-1,this.siz(i))+1;
                                case 'replicate'
                                    preref.subs{i} = ones(padsize(i),1);
                                    postref.subs{i} = this.siz(i * ones(padsize(i), 1));
                                case 'symmetric'
                                    preref.subs{i} = mod((padsize(i):-1:1)-1,this.siz(i))+1;
                                    postref.subs{i} = mod(this.siz(i)-(0:padsize(i)-1)-1,this.siz(i))+1;
                            end
                            
                            this = cat(i, this.subsref(preref), this, this.subsref(postref));
                        end
                    end
                case 'post'
                    if isnumeric(padval) && padval==0
                        this.siz = this.siz + padsize;
                    elseif isnumeric(padval)
                        for i=1:numel(padsize)
                            tmpsiz = this.siz;
                            tmpsiz(i)=padsize(i);
                            
                            tmp = SparseMat([],[],tmpsiz);
                            tmp.subs = tmp.ind2sub(1:tmp.numel);
                            tmp.vals = padval(ones(tmp.numel,1),1);
                            
                            this = cat(i, this, tmp);
                        end
                    else
                        for i=1:numel(padsize)
                            sref = struct('type', '()', 'subs', []);
                            tmp  = {':'};
                            sref.subs = tmp(1, ones(ndims(this), 1));
                            
                            switch padval
                                case 'circular'
                                    sref.subs{i} = mod((1:padsize(i))-1,this.siz(i))+1;
                                case 'replicate'
                                    sref.subs{i} = this.siz(i * ones(padsize(i), 1));
                                case 'symmetric'
                                    sref.subs{i} = mod(this.siz(i)-(0:padsize(i)-1)-1,this.siz(i))+1;
                            end
                            
                            this = cat(i, this, this.subsref(sref));
                        end
                    end
                case 'pre'
                    if isnumeric(padval) && padval==0
                        this.siz = this.siz + padsize;
                        this.subs = this.subs + padsize(ones(this.nnz,1), :);
                    elseif isnumeric(padval)
                        for i=1:numel(padsize)
                            tmpsiz = this.siz;
                            tmpsiz(i)=padsize(i);
                            
                            tmp = SparseMat([],[],tmpsiz);
                            tmp.subs = tmp.ind2sub(1:tmp.numel);
                            tmp.vals = padval(ones(tmp.numel,1),1);
                            
                            this = cat(i, tmp, this);
                        end
                    else
                        for i = 1:numel(padsize)
                            sref = struct('type', '()', 'subs', []);
                            tmp  = {':'};
                            sref.subs = tmp(1, ones(ndims(this), 1));
                            
                            switch padval
                                case 'circular'
                                    sref.subs{i} = mod(this.siz(i)-padsize(i)+(1:padsize(i))-1, this.siz(i)) + 1;
                                case 'replicate'
                                    sref.subs{i} = ones(1, padsize(i), 1);
                                case 'symmetric'
                                    sref.subs{i} = mod((padsize(i):-1:1)-1, this.siz(i)) + 1;
                            end
                            
                            this = cat(i, this.subsref(sref), this);
                        end
                    end
            end
        end
    end
    
    %scalar functions
    methods
        function this = abs(this)
            this.vals = abs(this.vals);
        end
        
        function this = sign(this)
            this.vals = sign(this.vals);
        end
        
        function C = sqrt(A)
            C = A;
            C.vals = sqrt(A.vals);
        end
        
        function C = realsqrt(A)
            C = A;
            C.vals = realsqrt(A.vals);
        end
        
        function C = exp(A)
            C = A;
            
            C.subs = A.ind2sub(1:A.numel);
            C.vals = ones(A.numel, 1);
            
            C.vals(A.inds) = exp(A.vals);
            C.subs(A.inds(A.vals==-Inf),:)=[];
            C.vals(A.inds(A.vals==-Inf),:)=[];
        end
        
        function C = expm1(A)
            C = A;
            
            C.vals = expm1(A.vals);
        end
        
        function C = log(A)
            C = A;
            
            C.subs = A.ind2sub(1:A.numel);
            C.vals = -Inf(A.numel, 1);
            
            C.vals(A.inds) = log(A.vals);
            C.subs(A.inds(A.vals==1),:)=[];
            C.vals(A.inds(A.vals==1),:)=[];
        end
        
        function C = log2(A)
            C = A;
            
            C.subs = A.ind2sub(1:A.numel);
            C.vals = -Inf(A.numel, 1);
            
            C.vals(A.inds) = log2(A.vals);
            C.subs(A.inds(A.vals==1),:)=[];
            C.vals(A.inds(A.vals==1),:)=[];
        end
        
        function C = log10(A)
            C = A;
            
            C.subs = A.ind2sub(1:A.numel);
            C.vals = -Inf(A.numel, 1);
            
            C.vals(A.inds) = log10(A.vals);
            C.subs(A.inds(A.vals==1),:)=[];
            C.vals(A.inds(A.vals==1),:)=[];
        end
        
        function C = log1p(A)
            C = A;
            
            C.vals = log1p(A.vals);
        end
        
        function C = reallog(A)
            C = A;
            
            C.subs = A.ind2sub(1:A.numel);
            C.vals = -Inf(A.numel, 1);
            
            C.vals(A.inds) = reallog(A.vals);
            C.subs(A.inds(A.vals==1),:)=[];
            C.vals(A.inds(A.vals==1),:)=[];
        end
        
        %         function C = sqrtm(A); end;
        %         function C = expm(A); end;
        %         function C = logm(A); end;
        %         function C = funm(A); end;
    end
    
    %complex numbers
    methods
        function this = real(this)
            this.vals = real(this.vals);
            
            this = this.normalize(false, true, false, 1:this.ndims);
        end
        
        function this = imag(this)
            this.vals = imag(this.vals);
            
            this = this.normalize(false, true, false, 1:this.ndims);
        end
        
        function this = angle(this)
            this.vals = angle(this.vals);
            
            this = this.normalize(false, true, false, 1:this.ndims);
        end
        
        function this = conj(this)
            this.vals = conj(this.vals);
        end
        
        function this = isreal(this)
            import edu.stanford.covert.util.SparseMat;
            
            if isreal(this.vals)
                this.subs = [1 1];
                this.vals = true;
            else
                this.subs = zeros(0,2);
                this.vals = false(0,1);
            end
            this.siz = [1 1];
            
            this = this.normalize(true, false, false, zeros(1, 2));
        end
    end
    
    %logic
    methods
        function this = not(this)
            warning('SparseMat:inefficient', 'Negating sparse matrices can be very inefficient. Consider refactoring your code to avoid negation.')
            this.subs = this.ind2sub(setdiff(1:this.numel, this.inds));
            this.vals = true(size(this.subs,1), 1);
        end
        
        function C = or(A, B)
            import edu.stanford.covert.util.SparseMat;
            
            if numel(A) == 1 || ~isa(A, 'SparseMat')
                tmp = B;
                B = A;
                A = tmp;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                if B
                    subs = A.ind2sub(1:A.numel);
                else
                    subs = A.subs;
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes', 'size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                subs = A.unique_subs([A.subs; B.subs], A.siz);
            end
            
            C = A;
            C.subs = subs;
            C.vals = true(size(subs,1), 1);
            C = C.normalize(true, true, true, 1:A.ndims);
        end
        
        function C = and(A, B)
            import edu.stanford.covert.util.SparseMat;
            
            if numel(A) == 1 || ~isa(A, 'SparseMat')
                tmp = B;
                B = A;
                A = tmp;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                if B
                    subs = A.subs;
                else
                    subs = zeros(0,2);
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes', 'size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                subs = A.intersect_subs(A.subs, B.subs, A.siz);
            end
            
            C = A;
            C.subs = subs;
            C.vals = true(size(subs,1), 1);
            C = C.normalize(true, true, true, 1:A.ndims);
        end
        
        function C = xor(A, B)
            import edu.stanford.covert.util.SparseMat;
            
            if numel(A) == 1 || ~isa(A, 'SparseMat')
                tmp = B;
                B = A;
                A = tmp;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                if B
                    subs = A.ind2sub(setdiff(1:A.numel, A.inds));
                else
                    subs = A.subs;
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes', 'size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                subs = setxor(A.subs, B.subs, 'rows');
            end
            
            C = A;
            C.subs = subs;
            C.vals = true(size(subs,1), 1);
            C = C.normalize(true, true, true, 1:A.ndims);
        end
        
        function C = nor(A, B)
            import edu.stanford.covert.util.SparseMat;
            
            if numel(A) == 1 || ~isa(A, 'SparseMat')
                tmp = B;
                B = A;
                A = tmp;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                if B
                    subs = zeros(0,1);
                else
                    subs = A.ind2sub(setdiff(1:A.numel, A.inds));
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes', 'size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                subs = A.ind2sub(intersect(setdiff(1:A.numel, A.inds), setdiff(1:B.numel, B.inds)));
            end
            
            C = A;
            C.subs = subs;
            C.vals = true(size(subs,1), 1);
            C = C.normalize(true, true, true, 1:A.ndims);
        end
        
        function C = nand(A, B)
            import edu.stanford.covert.util.SparseMat;
            
            if numel(A) == 1 || ~isa(A, 'SparseMat')
                tmp = B;
                B = A;
                A = tmp;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                if B
                    subs = A.ind2sub(setdiff(1:A.numel, A.inds));
                else
                    subs = A.ind2sub(1:A.numel);
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes', 'size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                subs = A.ind2sub(setdiff(1:A.numel, intersect(A.inds, B.inds)));
            end
            
            C = A;
            C.subs = subs;
            C.vals = true(size(subs,1), 1);
            C = C.normalize(true, true, true, 1:A.ndims);
        end
    end
    
    %operators
    methods
        function this = uminus(this)
            this.vals = -this.vals;
        end
        
        function this = uplus(this)
        end
        
        function C = plus(A, B)
            import edu.stanford.covert.util.SparseMat;
            
            if numel(A) == 1 || ~isa(A, 'SparseMat')
                tmp = B;
                B = A;
                A = tmp;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B, 'SparseMat')
                    B=full(B);
                end
                
                if B==0
                    subs = A.subs;
                    vals = A.vals;
                else
                    subs = [A.subs;
                        A.ind2sub(setdiff(1:A.numel,A.inds))];
                    vals = [A.vals + B;
                        B(ones(A.numel-A.nnz, 1),1)];
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes', 'size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                subs = A.subs;
                vals = A.vals;
                
                [tfs, idxs] = B.ismember_subs(B.subs, A.subs, B.siz);
                vals(idxs(tfs)) = vals(idxs(tfs)) + B.vals(tfs);
                
                subs = [subs; B.subs(~tfs,:)];
                vals = [vals; B.vals(~tfs,:)];
            end
            
            C = A;
            C.subs = subs;
            C.vals = vals;
            C = C.normalize(true, true, true, 1:A.ndims);
        end
        
        function C = minus(A,B)
            C = A + -B;
        end
        
        function C = times(A,B)
            import edu.stanford.covert.util.SparseMat;
            
            if numel(A) == 1 || ~isa(A, 'SparseMat')
                tmp = B;
                B = A;
                A = tmp;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                if isnan(B)
                    subs = A.ind2sub(1:A.numel);
                    vals = B(ones(A.numel, 1),1);
                elseif isinf(B)
                    subs = [
                        A.subs;
                        A.ind2sub(setdiff(1:A.numel, A.inds))];
                    vals = [
                        A.vals * B;
                        NaN(A.numel - A.nnz, 1)];
                else
                    subs = A.subs;
                    vals = A.vals * B;
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes', 'size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                [~, idxs1, idxs2] = intersect(A.inds, B.inds);
                subs = A.subs(idxs1, :);
                vals = A.vals(idxs1, :).*B.vals(idxs2, :);
                
                tfs1 = A.ismember_subs(A.subs, B.subs, A.siz);
                tfs2 = B.ismember_subs(B.subs, A.subs, B.siz);
                subs = [subs;
                    A.subs(~tfs1, :);
                    B.subs(~tfs2, :)];
                vals = [vals;
                    A.vals(~tfs1, :) * 0;
                    B.vals(~tfs2, :) * 0];
            end
            
            C = A;
            C.subs = subs;
            C.vals = vals;
            C = C.normalize(true, true, true, 1:A.ndims);
        end
        
        function C = rdivide(A, B)
            import edu.stanford.covert.util.SparseMat;
            
            if isa(B,'SparseMat')
                subs = [
                    B.subs;
                    B.ind2sub(setdiff(1:B.numel,B.inds))];
                vals = [
                    1./B.vals;
                    Inf(B.numel-B.nnz, 1)];                
                B = SparseMat(subs, vals, B.siz);
            else
                B = 1./B;
            end
            C = A .* B;
        end
        
        function C = ldivide(A, B)
            C = B ./ A;
        end
        
        function C = mtimes(A, B)
            if numel(A)==1 || numel(B)==1
                C = A.*B;
            elseif ndims(A)<=2 && ndims(B)<=2
                C = tprod(A, B, 2, 1);
            else
                throw(MException('SparseMat:invalidDimensions','Inputs must be 2-D, or at least one input must be scalar'));
            end
        end
        
        function C = tprod(A, B, Aid, Bid)
            import edu.stanford.covert.util.SparseMat;
            
            if any(Aid<1) || ~isequal(Aid, ceil(real(Aid))) || any(Aid > ndims(A)) || ~all(diff(sort(Aid)))
                throw(MException('SparseMat:invalidDimensions','inner dimensions must be unique positive real integers'));
            end
            if any(Bid<1) || ~isequal(Bid, ceil(real(Bid))) || any(Bid > ndims(B)) || ~all(diff(sort(Bid)))
                throw(MException('SparseMat:invalidDimensions','inner dimensions must be unique positive real integers'));
            end
            if numel(Aid) ~= numel(Bid)
                throw(MException('SparseMat:invalidDimensions','inner dimensions must have same size'));
            end
            
            %convert to SparseMat type
            if ~isa(A, 'SparseMat')
                A = SparseMat(A);
            end
            
            if ~isa(B, 'SparseMat')
                B = SparseMat(B);
            end
            
            %compute size of result
            Aod = setdiff(1:A.ndims, Aid);
            Bod = setdiff(1:B.ndims, Bid);
            siz = [size(A, Aod) size(B, Bod)];
            
            %compute tensor product
            Anansubs = A.subs(isnan(A.vals),:);
            Bnansubs = B.subs(isnan(B.vals),:);
            A.subs(isnan(A.vals),:)=[];
            B.subs(isnan(B.vals),:)=[];
            A.vals(isnan(A.vals),:)=[];
            B.vals(isnan(B.vals),:)=[];
            
            Ainfsubs = A.subs(isinf(A.vals),:);
            Binfsubs = B.subs(isinf(B.vals),:);
            Azrosubs = A.ind2sub(setdiff(1:A.numel, A.inds));
            Bzrosubs = B.ind2sub(setdiff(1:B.numel, B.inds));
            
            subs = [repmat(A.subs, B.nnz, 1) B.sort_subs(repmat(B.subs, A.nnz, 1), B.siz)
                repmat(Ainfsubs, size(Bzrosubs,1), 1) B.sort_subs(repmat(Bzrosubs, size(Ainfsubs,1), 1), B.siz);
                A.sort_subs(repmat(Azrosubs, size(Binfsubs,1), 1), A.siz) repmat(Binfsubs, size(Azrosubs,1), 1)];
            vals = [reshape(A.vals*B.vals',[],1);
                NaN(size(Ainfsubs,1)*size(Bzrosubs,1),1)
                NaN(size(Binfsubs,1)*size(Azrosubs,1),1)];
            
            badIdxs = find(subs(:,Aid)~=subs(:,A.ndims + Bid));
            subs(badIdxs, :)=[];
            vals(badIdxs, :)=[];
            subs(:, [Aid A.ndims + Bid])=[];
            
            subs = [subs;
                repmat(Anansubs(:,Aod),prod(B.siz(Bod)),1) B.ind2sub(sort(repmat(1:prod(B.siz(Bod)), size(Anansubs,1), 1)), B.siz(Bod));
                A.ind2sub(sort(repmat(1:prod(A.siz(Aod)), size(Bnansubs,1), 1)), A.siz(Aod)) repmat(Bnansubs(:,Bod),prod(A.siz(Aod)),1)];
            vals = [vals;
                NaN(size(Anansubs,1) * prod(B.siz(Bod)),1);
                NaN(size(Bnansubs,1) * prod(A.siz(Aod)),1)];
            
            [subs, order] = A.sort_subs(subs, siz);
            sums = vals(order,:);
            for i=2:size(subs,1)
                if isequal(subs(i-1,:),subs(i,:))
                    sums(i) = sums(i-1)+sums(i);
                    sums(i-1)=0;
                end
            end
            
            idxs = find(sums);
            subs = subs(idxs,:);
            vals = sums(idxs,:);
            
            subs(vals==0,:)=[];
            vals(vals==0,:)=[];
            
            %return as SparseMat
            if ~isa(A,'SparseMat') || (isa(B,'SparseMat') && ~strcmp(class(B),'edu.stanford.covert.util.SparseMat'))
                C = B;
            else
                C = A;
            end
            C.subs = subs;
            C.vals = vals;
            C.siz = siz;
            C = C.normalize(true, true, true, zeros(1, C.ndims));
        end
        
        function C = mrdivide(A, B)
            if numel(A)==1 || numel(B)==1
                C = A./B;
            else
                throw(MException('SparseMat:undefined','operation not yet implemented'));
            end
        end
        
        function C = mldivide(A, B)
            if numel(A)==1 || numel(B)==1
                C = A.\B;
            else
                throw(MException('SparseMat:undefined','operation not yet implemented'));
            end
        end
        
        function C = power(A, B)
            import edu.stanford.covert.util.SparseMat;
            
            if numel(A) == 1
                %A is scalar
                
                if isa(A,'SparseMat')
                    A=full(A);
                end
                
                if isnan(A)
                    subs = B.ind2sub(1:B.numel);
                    vals = A(ones(B.numel, 1),1);
                elseif isinf(A)
                    subs = [
                        B.subs;
                        B.ind2sub(setdiff(1:B.numel, B.inds))];
                    vals = [
                        A.^ B.vals;
                        repmat(A.^0, B.numel-B.nnz, 1)];
                else
                    subs = B.ind2sub(1:B.numel);
                    vals = ones(B.numel,1);
                    vals(B.inds)=A.^B.vals;
                end
            elseif numel(B) == 1
                %B is a scalar
                
                if isa(B,'SparseMat')
                    B=full(B);
                end
                
                if isnan(B)
                    subs = A.ind2sub(1:A.numel);
                    vals = B(ones(A.numel, 1), 1);
                elseif isinf(B)
                    subs = [
                        A.subs;
                        A.ind2sub(setdiff(1:A.numel, A.inds))];
                    vals = [
                        A.vals .^ B;
                        repmat(0.^B, A.numel-A.nnz, 1)];
                elseif B==0
                    subs = A.ind2sub(1:A.numel);
                    vals = ones(A.numel,1);
                    vals(A.inds)=A.vals.^0;
                else
                    subs = A.subs;
                    vals = A.vals .^ B;
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes','size RHS must equal that of LHS'))
                end
                
                if ~isa(A, 'SparseMat')
                    A = SparseMat(A);
                end
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                [~,idxs1,idxs2]=intersect(A.inds, B.inds);
                subs = A.subs(idxs1,:);
                vals = A.vals(idxs1,:).^B.vals(idxs2,:);
                
                tfs1 = A.ismember_subs(A.subs, B.subs, A.siz);
                tfs2 = B.ismember_subs(B.subs, A.subs, B.siz);
                subs = [subs;
                    A.subs(~tfs1,:);
                    B.subs(~tfs2,:)];
                vals = [vals;
                    A.vals(~tfs1,:).^0;
                    0.^B.vals(~tfs2,:)];
                
                subs = [subs;
                    A.ind2sub(setdiff(1:A.numel,[A.inds;B.inds]))];
                vals = [vals;
                    ones(size(subs,1)-size(vals,1),1)];
            end
            
            subs(vals==0,:)=[];
            vals(vals==0,:)=[];
            
            if ~isa(A,'SparseMat') || (isa(B,'SparseMat') && ~strcmp(class(B),'edu.stanford.covert.util.SparseMat'))
                C = B;
            else
                C = A;
            end
            C.subs = subs;
            C.vals = vals;
            C = C.normalize(true, true, true, 1:C.ndims);
        end
        
        function C = mpower(A, B)
            throw(MException('SparseMat:undefined','operation not yet implemented'));
        end
    end
    
    %comparison
    methods
        function val = isequal(A, B)
            import edu.stanford.covert.util.SparseMat;
            
            val = isa(A, 'SparseMat') && ...
                isequal(class(A), class(B)) && ...
                isequal(A.siz, B.siz) && ...
                isequal(A.subs, B.subs) && ...
                isequal(A.vals, B.vals);
        end
        
        function assertElementsAlmostEqual(A, B, varargin)
            import edu.stanford.covert.util.SparseMat;
            
            assertEqual(size(A), size(B));
            
            if ~isa(A, 'SparseMat')
                A = SparseMat(A);
            end
            if ~isa(B, 'SparseMat')
                B = SparseMat(B);
            end
            
            assertEqual(A.subs, B.subs);
            assertElementsAlmostEqual(A.vals, B.vals, varargin{:});
        end
        
        function C = eq(A, B)
            import edu.stanford.covert.util.SparseMat;
            
            if numel(A) == 1 || ~isa(A, 'SparseMat')
                tmp = B;
                B = A;
                A = tmp;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                if B == 0
                    subs = A.ind2sub(setdiff(1:A.numel, A.inds));
                else
                    subs = A.subs(A.vals == B, :);
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes', 'size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                %non-zeros
                [~, idxs] = intersect([A.inds A.vals], [B.inds B.vals], 'rows');
                subs = A.subs(idxs,:);
                
                %zeros
                subs = [subs;
                    A.ind2sub(setdiff(1:A.numel,[A.inds; B.inds]))];
            end
            
            C = A;
            C.subs = subs;
            C.vals = true(size(subs,1), 1);
            C = C.normalize(true, true, true, 1:A.ndims);
        end
        
        function C = ne(A,B)
            import edu.stanford.covert.util.SparseMat;
            
            if numel(A) == 1 || ~isa(A, 'SparseMat')
                tmp = B;
                B = A;
                A = tmp;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B,'SparseMat')
                    B = full(B);
                end
                
                if B == 0
                    subs = A.subs;
                else
                    subs = [
                        A.ind2sub(setdiff(1:A.numel, A.inds));
                        A.subs(A.vals ~= B, :)];
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes', 'size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                [~, idxs1, idxs2] = intersect(A.inds, B.inds);
                
                tfs1 = A.ismember_subs(A.subs, B.subs, A.siz);
                tfs2 = B.ismember_subs(B.subs, A.subs, B.siz);
                
                subs = [
                    A.subs(idxs1(A.vals(idxs1, :) ~= B.vals(idxs2, :)), :);
                    A.subs(~tfs1, :);
                    B.subs(~tfs2, :)];
            end
            
            C = A;
            C.subs = subs;
            C.vals = true(size(subs,1), 1);
            C = C.normalize(true, true, true, 1:A.ndims);
        end
        
        function C = ge(A,B)
            import edu.stanford.covert.util.SparseMat;
            
            if ~isa(A, 'SparseMat')
                C = B <= A;
                return;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                if isnan(B)
                    subs = zeros(0,2);
                elseif B > 0
                    subs = A.subs(A.vals >= B, :);
                else
                    subs = [
                        A.subs(A.vals >= B, :);
                        A.ind2sub(setdiff(1:A.numel, A.inds))];
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes', 'size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                [~, idxs1, idxs2] = intersect(A.inds, B.inds);
                subs = A.subs(idxs1(A.vals(idxs1, :) >= B.vals(idxs2, :)), :);
                
                tfs1 = A.ismember_subs(A.subs, B.subs, A.siz);
                tfs2 = B.ismember_subs(B.subs, A.subs, B.siz);
                idxs1 = find(~tfs1);
                idxs2 = find(~tfs2);
                
                subs = [subs;
                    A.subs(idxs1(A.vals(idxs1, :) >= 0), :);
                    B.subs(idxs2(B.vals(idxs2, :) <= 0), :)
                    A.ind2sub(intersect(setdiff(1:A.numel, A.inds), setdiff(1:B.numel, B.inds)))];
            end
            
            C = A;
            C.subs = subs;
            C.vals = true(size(subs,1), 1);
            C = C.normalize(true, true, true, 1:A.ndims);
        end
        
        function C = gt(A,B)
            import edu.stanford.covert.util.SparseMat;
            
            if ~isa(A,'SparseMat')
                C = B < A;
                return;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                if isnan(B)
                    subs = zeros(0,2);
                elseif B>=0
                    subs = A.subs(A.vals > B, :);
                else
                    subs = [
                        A.subs(A.vals > B, :);
                        A.ind2sub(setdiff(1:A.numel, A.inds))];
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes','size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                [~, idxs1, idxs2] = intersect(A.inds, B.inds);
                subs = A.subs(idxs1(A.vals(idxs1,:) > B.vals(idxs2,:)),:);
                
                tfs1 = A.ismember_subs(A.subs, B.subs, A.siz);
                tfs2 = B.ismember_subs(B.subs, A.subs, B.siz);
                idxs1 = find(~tfs1);
                idxs2 = find(~tfs2);
                
                subs = [subs;
                    A.subs(idxs1(A.vals(idxs1,:)>0),:);
                    B.subs(idxs2(B.vals(idxs2,:)<0),:)];
            end
            
            C = A;
            C.subs = subs;
            C.vals = true(size(subs,1), 1);
            C = C.normalize(true, true, true, 1:A.ndims);
        end
        
        function C = le(A,B)
            import edu.stanford.covert.util.SparseMat;
            
            if ~isa(A,'SparseMat')
                C = B >= A;
                return;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                if isnan(B)
                    subs = zeros(0,2);
                elseif B < 0
                    subs = A.subs(A.vals <= B, :);
                else
                    subs = [
                        A.subs(A.vals <= B, :);
                        A.ind2sub(setdiff(1:A.numel, A.inds))];
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes','size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                [~, idxs1, idxs2] = intersect(A.inds, B.inds);
                subs = A.subs(idxs1(A.vals(idxs1,:) <= B.vals(idxs2,:)),:);
                
                tfs1 = A.ismember_subs(A.subs, B.subs, A.siz);
                tfs2 = B.ismember_subs(B.subs, A.subs, B.siz);
                idxs1 = find(~tfs1);
                idxs2 = find(~tfs2);
                
                subs = [subs;
                    A.subs(idxs1(A.vals(idxs1,:)<=0),:);
                    B.subs(idxs2(B.vals(idxs2,:)>=0),:);
                    A.ind2sub(intersect(setdiff(1:A.numel,A.inds), setdiff(1:B.numel, B.inds)))];
            end
            
            C = A;
            C.subs = subs;
            C.vals = true(size(subs,1), 1);
            C = C.normalize(true, true, true, 1:A.ndims);
        end
        
        function C = lt(A,B)
            import edu.stanford.covert.util.SparseMat;
            
            if ~isa(A,'SparseMat')
                C = B > A;
                return;
            end
            
            %B is scalar
            if numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                if isnan(B)
                    subs = zeros(0,2);
                elseif B <= 0
                    subs = A.subs(A.vals < B, :);
                else
                    subs = [
                        A.subs(A.vals < B, :);
                        A.ind2sub(setdiff(1:A.numel, A.inds))];
                end
            else
                if ~isequal(size(A), size(B))
                    throw(MException('SparseMat:invalidSizes', 'size RHS must equal that of LHS'))
                end
                
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                [~, idxs1, idxs2] = intersect(A.inds, B.inds);
                subs = A.subs(idxs1(A.vals(idxs1,:) < B.vals(idxs2,:)),:);
                
                tfs1 = A.ismember_subs(A.subs, B.subs, A.siz);
                tfs2 = B.ismember_subs(B.subs, A.subs, B.siz);
                idxs1 = find(~tfs1);
                idxs2 = find(~tfs2);
                
                subs = [subs;
                    A.subs(idxs1(A.vals(idxs1,:)<0),:);
                    B.subs(idxs2(B.vals(idxs2,:)>0),:)];
            end
            
            C = A;
            C.subs = subs;
            C.vals = true(size(subs,1), 1);
            C = C.normalize(true, true, true, 1:A.ndims);
        end
    end
    
    %set operators
    methods
        function val = unique(this)
            import edu.stanford.covert.util.SparseMat;
            
            vals = unique(this.vals);
            
            if this.nnz < this.numel
                vals = sort([vals;0]);
            end
            
            subs = [(1:size(vals,1))' ones(size(vals))];
            
            if ~any(this.siz)
                siz = [0 0];
            elseif this.ndims==2 && this.siz(1)==1 && ~isempty(vals)
                siz = [1 size(vals,1)];
            else
                siz = [size(vals,1) 1];
            end
            
            val = this;
            val.subs = subs;
            val.vals = vals;
            val.siz = siz;
            val = val.normalize(true, true, true, zeros(1, 2));
        end
        
        %function C = union(A, B, varargin); end
        %function C = intersect(A, B, varargin); end
        %function C = setdiff((A, B, varargin); end
        %function C = setxor(A, B, varargin); end
        %function C = ismember(A, B, varargin); end
    end
    
    %reduction
    methods
        function B = all(A, dim)
            import edu.stanford.covert.util.SparseMat;
            
            if ~exist('dim', 'var')
                dim = 1;
            end
            
            if dim < 1 || dim ~= ceil(dim)
                throw(MException('SparseMat:invalidDimension', 'invalid dimension'));
            elseif dim>size(A.siz, 2)
                B = A;
                B.vals = true(size(A.vals));
                B = B.normalize(true, true, true, 1:A.ndims);
                return;
            end
            
            siz = A.siz;
            siz(dim) = 1;
            
            subs = A.sort_subs(A.subs(:,[1:dim-1 dim+1:end]), A.siz(1,[1:dim-1 dim+1:end]));
            
            if ~A.isunique_subs(subs, A.siz([1:dim-1 dim+1:end]))
                [~, first, ~]     = unique(subs, 'rows', 'first');
                [subs, last, ind] = unique(subs, 'rows', 'last');
                subs = subs(ind(last-first+1 == A.siz(dim)), :);
            end
            subs = [subs(:, 1:dim-1) ones(size(subs, 1), 1) subs(:, dim:end)];
            
            vals = true(size(subs, 1), 1);
            
            B = A;
            B.subs = subs;
            B.vals = vals;
            B.siz = siz;
            oldDims = 1:A.ndims;
            oldDims(dim) = 0;
            B = B.normalize(true, true, true, oldDims);
        end
        
        function B = any(A, dim)
            import edu.stanford.covert.util.SparseMat;
            
            if ~exist('dim', 'var')
                dim = 1;
            end
            
            if dim < 1 || dim ~= ceil(dim)
                throw(MException('SparseMat:invalidDimension','invalid dimension'));
            elseif dim > A.ndims
                B = A;
                B.vals = true(size(A.vals));
                B = B.normalize(false, false, false, 1:A.ndims);
                return;
            end
            
            siz = A.siz;
            siz(dim) = 1;
            
            subs = A.unique_subs(A.subs(:, [1:dim-1 dim+1:end]), A.siz(1, [1:dim-1 dim+1:end]));
            subs = [subs(:, 1:dim-1) ones(size(subs, 1), 1) subs(:, dim:end)];
            
            vals = true(size(subs, 1), 1);
            
            B = A;
            B.subs = subs;
            B.vals = vals;
            B.siz = siz;
            oldDims = 1:A.ndims;
            oldDims(dim) = 0;
            B = B.normalize(true, true, true, oldDims);
        end
        
        function C = min(A, B, dim)
            import edu.stanford.covert.util.SparseMat;
            
            if ~exist('B', 'var')
                B = [];
                dim = 1;
            end
            
            if numel(A) == 1 || ~isa(A, 'SparseMat')
                tmp = B;
                B = A;
                A = tmp;
            end
            
            if exist('dim', 'var')
                if ~isempty(B)
                    throw(MException('SparseMat:unsupportedSyntax', 'MIN with two matrices to compare and a working dimension is not supported.'))
                end
                
                if dim < 1 || dim ~= ceil(dim)
                    throw(MException('SparseMat:invalidDimension', 'invalid dimension'));
                elseif dim > size(A.siz, 2)
                    C = A;
                    return
                else
                    subsNaNVals = sortrows([A.vals ~isnan(A.vals) A.subs(:,[1:dim-1 dim+1:end])], A.ndims+1:-1:1);
                    subsVals = subsNaNVals(:, [1 3:end]);
                    mins = true(A.nnz, 1);
                    nnzs = zeros(A.nnz, 1);
                    mins(1) = true;
                    nnzs(1) = 1;
                    for i = 2:A.nnz
                        if isequal(subsVals(i-1, 2:end),subsVals(i, 2:end))
                            nnzs(i) = nnzs(i-1);
                            nnzs(i-1) = 0;
                            
                            if isnan(subsVals(i-1,1))
                                mins(i-1) = false;
                            else
                                mins(i) = false;
                            end
                        end
                        nnzs(i) = nnzs(i) + 1;
                    end
                    
                    mins = find(mins);
                    [~, ~, nnzs] = find(nnzs);
                    mins(nnzs < A.siz(dim) & (subsVals(mins, 1) > 0 | isnan(subsVals(mins, 1)))) = [];
                    
                    C = A;
                    C.vals = subsVals(mins, 1);
                    oldDims = 1:A.ndims;
                    if dim == A.ndims && dim > 2
                        C.subs = subsVals(mins, 2:dim);
                        C.siz(dim) = [];
                        oldDims(dim) = [];
                    else
                        C.subs = [subsVals(mins,2:dim) ones(length(mins), 1) subsVals(mins, dim+1:end)];
                        C.siz(dim) = 1;
                        oldDims(dim) = 0;
                    end
                    
                    C = C.normalize(true, false, false, oldDims);
                end
            elseif numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                C = A;
                if isnan(B)
                elseif B > 0
                    C.vals(A.vals > B | isnan(A.vals)) = B;
                elseif B == 0
                    C.subs(A.vals > B | isnan(A.vals), :) = [];
                    C.vals(A.vals > B | isnan(A.vals), :) = [];
                else
                    C.subs = A.ind2sub(1:A.numel);
                    C.vals = B(ones(A.numel, 1), 1);
                    C.vals(A.inds) = min(A.vals, B);
                end
            elseif ~isequal(size(A), size(B))
                throw(MException('SparseMat:invalidInput', 'Matrix dimensions must agree.'));
            else
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                [subs, idxs1, idxs2] = A.intersect_subs(A.subs, B.subs, A.siz);
                vals = min(A.vals(idxs1, :), B.vals(idxs2, :));
                
                tfs1 = A.ismember_subs(A.subs, B.subs, A.siz);
                tfs2 = B.ismember_subs(B.subs, A.subs, B.siz);
                idxs1 = find(~tfs1);
                idxs2 = find(~tfs2);
                idxs1 = idxs1(A.vals(idxs1) < 0);
                idxs2 = idxs2(B.vals(idxs2) < 0);
                
                subs = [subs;
                    A.subs(idxs1, :);
                    B.subs(idxs2, :)];
                vals = [vals;
                    A.vals(idxs1, :);
                    B.vals(idxs2, :)];
                
                C = A;
                C.subs = subs;
                C.vals = vals;
                C = C.normalize(true, true, true, 1:A.ndims);
            end
        end
        
        function C = max(A, B, dim)
            import edu.stanford.covert.util.SparseMat;
            
            if ~exist('B', 'var')
                B = [];
                dim = 1;
            end
            
            if numel(A) == 1 || ~isa(A, 'SparseMat')
                tmp = B;
                B = A;
                A = tmp;
            end
            
            if exist('dim', 'var')
                if ~isempty(B)
                    throw(MException('SparseMat:unsupportedSyntax', 'MIN with two matrices to compare and a working dimension is not supported.'))
                end
                
                if dim < 1 || dim ~= ceil(dim)
                    throw(MException('SparseMat:invalidDimension', 'invalid dimension'));
                elseif dim > size(A.siz, 2)
                    C = A;
                    return
                else
                    subsNaNVals = sortrows([A.vals ~isnan(A.vals) A.subs(:,[1:dim-1 dim+1:end])], A.ndims+1:-1:1);
                    subsVals = subsNaNVals(:, [1 3:end]);
                    maxs = false(A.nnz, 1);
                    nnzs = zeros(A.nnz, 1);
                    maxs(1) = true;
                    nnzs(1) = 1;
                    for i = 2:A.nnz
                        if isequal(subsVals(i-1, 2:end), subsVals(i, 2:end))
                            nnzs(i) = nnzs(i-1);
                            nnzs(i-1) = 0;
                            
                            maxs(i) = true;
                            maxs(i-1) = false;
                        else
                            maxs(i) = true;
                        end
                        nnzs(i) = nnzs(i) + 1;
                    end
                    
                    maxs = find(maxs);
                    [~, ~, nnzs] = find(nnzs);
                    maxs(nnzs < A.siz(dim) & (subsVals(maxs, 1) < 0 | isnan(subsVals(maxs, 1)))) = [];
                    
                    C = A;
                    C.vals = subsVals(maxs, 1);
                    oldDims = 1:A.ndims;
                    if dim == A.ndims && dim > 2
                        C.subs = subsVals(maxs, 2:dim);
                        C.siz(dim) = [];
                        oldDims(dim) = [];
                    else
                        C.subs = [subsVals(maxs, 2:dim) ones(length(maxs), 1) subsVals(maxs, dim + 1:end)];
                        C.siz(dim) = 1;
                        oldDims(dim) = 0;
                    end
                    
                    C = C.normalize(true, false, false, oldDims);
                end
            elseif numel(B) == 1
                if isa(B, 'SparseMat')
                    B = full(B);
                end
                
                C = A;
                if isnan(B)
                elseif B < 0
                    C.vals(A.vals < B | isnan(A.vals)) = B;
                elseif B == 0
                    C.subs(A.vals < B | isnan(A.vals), :) = [];
                    C.vals(A.vals < B | isnan(A.vals), :) = [];
                else
                    C.subs = A.ind2sub(1:A.numel);
                    C.vals = B(ones(A.numel, 1),1);
                    C.vals(A.inds) = max(A.vals, B);
                end
            elseif ~isequal(size(A), size(B))
                throw(MException('SparseMat:invalidInput','Matrix dimensions must agree.'));
            else
                if ~isa(B, 'SparseMat')
                    B = SparseMat(B);
                end
                
                [subs, idxs1, idxs2] = A.intersect_subs(A.subs, B.subs, A.siz);
                vals = max(A.vals(idxs1,:), B.vals(idxs2, :));
                
                tfs1 = A.ismember_subs(A.subs, B.subs, A.siz);
                tfs2 = B.ismember_subs(B.subs, A.subs, B.siz);
                idxs1 = find(~tfs1);
                idxs2 = find(~tfs2);
                idxs1 = idxs1(A.vals(idxs1) > 0);
                idxs2 = idxs2(B.vals(idxs2) > 0);
                
                subs = [subs;
                    A.subs(idxs1, :);
                    B.subs(idxs2, :)];
                vals = [vals;
                    A.vals(idxs1, :);
                    B.vals(idxs2, :)];
                
                C = A;
                C.subs = subs;
                C.vals = vals;
                C = C.normalize(true, true, true, 1:A.ndims);
            end
        end
        
        function C = range(A, dim)
            import edu.stanford.covert.util.SparseMat;
            
            if ~exist('dim', 'var')
                dim = 1;
            end
            
            if dim < 1 || dim ~= ceil(dim)
                throw(MException('SparseMat:invalidDimension', 'invalid dimension'));
            elseif dim > A.ndims
                C = A;
                C.subs(isfinite(C.vals), :) = [];
                C.vals = NaN(size(C.subs, 1), 1);
                return;
            end
            
            subsNaNVals = sortrows([A.vals ~isnan(A.vals) A.subs(:,[1:dim-1 dim+1:end])], A.ndims+1:-1:1);
            subsVals = subsNaNVals(:, [1 3:end]);
            mins = true(A.nnz, 1);
            maxs = false(A.nnz, 1);
            nnzs = zeros(A.nnz, 1);
            mins(1) = true;
            maxs(1) = true;
            nnzs(1) = 1;
            for i = 2:A.nnz
                if isequal(subsVals(i-1, 2:end),subsVals(i, 2:end))
                    nnzs(i) = nnzs(i-1);
                    nnzs(i-1) = 0;
                    
                    if isnan(subsVals(i-1,1))
                        mins(i-1) = false;
                    else
                        mins(i) = false;
                    end
                    
                    maxs(i) = true;
                    maxs(i-1) = false;
                else
                    maxs(i) = true;
                end
                nnzs(i) = nnzs(i) + 1;
            end
            
            [~, ~, nnzs] = find(nnzs);
            mins = find(mins);
            maxs = find(maxs);
            mins(nnzs < A.siz(dim) & (subsVals(mins, 1) > 0 | isnan(subsVals(mins, 1)))) = [];
            maxs(nnzs < A.siz(dim) & (subsVals(maxs, 1) < 0 | isnan(subsVals(maxs, 1)))) = [];
            
            minSubs = subsVals(mins, 2:end);
            maxSubs = subsVals(maxs, 2:end);
            minVals = subsVals(mins, 1);
            maxVals = subsVals(maxs, 1);
            
            [subs, idxs1, idxs2] = intersect(maxSubs, minSubs, 'rows');
            vals = maxVals(idxs1, :) - minVals(idxs2, :);
            subs(vals == 0, :) = [];
            vals(vals == 0, :) = [];
            
            subs = [subs;
                minSubs(setdiff(1:size(minSubs, 1), idxs2), :);
                maxSubs(setdiff(1:size(maxSubs, 1), idxs1), :)];
            vals = [vals;
                -minVals(setdiff(1:size(minSubs, 1), idxs2), :);
                maxVals(setdiff(1:size(maxSubs, 1), idxs1), :)];
            
            [subs, order] = sortrows(subs, A.ndims-1:-1:1);
            vals = vals(order, :);
            
            C = A;
            C.vals = vals;
            oldDims = 1:A.ndims;
            if dim == A.ndims && dim > 2
                C.subs = subs;
                C.siz(dim) = [];
                oldDims(dim) = [];
            else
                C.subs = [subs(:, 1:dim-1) ones(size(subs, 1), 1) subs(:, dim:end)];
                C.siz(dim) = 1;
                oldDims(dim) = 0;
            end
            
            C = C.normalize(true, false, false, oldDims);
        end
        
        function C = sum(A, dim)
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.util.SparseMat;
            
            if ~exist('dim', 'var')
                dim = 1;
            end
            
            if dim < 1 || dim ~= ceil(dim)
                throw(MException('SparseMat:invalidDimension', 'invalid dimension'));
            elseif dim > A.ndims
                C = A;
                return;
            end
            
            subs = A.subs;
            subs(:, dim) = 1;
            siz = A.siz;
            siz(dim) = 1;
            
            [subs, ~, J] = A.unique_subs(subs, siz);
            vals = accumarray(J, A.vals, [size(subs, 1) 1]);
            
            if isa(A, 'edu.stanford.covert.util.CircularSparseMat')
                oldDims = 1:A.ndims;
                oldDims(dim) = [];
                C = CircularSparseMat(subs, vals, siz, oldDims);
            else
                C = SparseMat(subs, vals, siz);
            end
        end
        
        function C = mean(A, dim)
            if ~exist('dim', 'var')
                dim = 1;
            end
            
            if dim < 1 || dim ~= ceil(dim)
                throw(MException('SparseMat:invalidDimension', 'invalid dimension'));
            elseif dim > A.ndims
                C = A;
                return;
            end
            
            C = sum(A, dim);
            C.vals = C.vals / A.siz(dim);
        end
        
        function C = var(A, flag, dim)
            import edu.stanford.covert.util.SparseMat;
            
            if ~exist('flag', 'var') || flag == 0
                norm = 2;
            else
                norm = 0;
            end
            
            if ~exist('dim', 'var')
                dim = 1;
            elseif flag == 0 && dim > 1
                throw(MException('SparseMat:unsupportedSyntax', 'unsupported syntax'))
            end
            
            if dim < 1 || dim ~= ceil(dim)
                throw(MException('SparseMat:invalidDimension', 'invalid dimension'));
            elseif dim > A.ndims
                C = A;
                C.subs(isfinite(C.vals), :) = [];
                C.vals = NaN(size(C.subs, 1), 1);
                return;
            end
            
            [subs, order] = A.sort_subs(A.subs(:, [1:dim-1 dim+1:end]), A.siz(1, [1:dim-1 dim+1:end]));
            
            sums = A.vals(order,:);
            sumSquares = A.vals(order, :) .^ 2;
            for i=2:A.nnz
                if isequal(subs(i-1, :), subs(i, :))
                    sums(i) = sums(i-1) + sums(i);
                    sums(i - 1) = 0;
                    
                    sumSquares(i) = sumSquares(i - 1) + sumSquares(i);
                    sumSquares(i - 1) = 0;
                end
            end
            
            idxs = find(sumSquares);
            subs = subs(idxs, :);
            vals = sumSquares(idxs, :) / (A.siz(dim) - norm) - (sums(idxs, :) / (A.siz(dim) - norm)) .^ 2;
            
            subs(vals == 0, :) = [];
            vals(vals == 0, :) = [];
            
            C = A;
            C.vals = vals;
            oldDims = 1:A.ndims;
            if dim == A.ndims && dim > 2
                C.subs = subs;
                C.siz(dim) = [];
                oldDims(dim) = [];
            else
                C.subs = [subs(:, 1:dim-1) ones(size(subs, 1), 1) subs(:, dim:end)];
                C.siz(dim) = 1;
                oldDims(dim) = 0;
            end
            
            C = C.normalize(true, false, false, oldDims);
        end
        
        function C = std(A, flag, dim)
            if ~exist('flag', 'var')
                flag = 0;
            end
            
            if ~exist('dim', 'var')
                dim = 1;
            end
            
            C = var(A, flag, dim);
            C.vals(~isfinite(C.vals), :) = NaN;
            
            C.subs(C.vals < 0, :) = [];
            C.vals(C.vals < 0, :) = [];
            
            C.vals = sqrt(C.vals);
        end
        
        %         function value = prod(this); end
        %         function value = cumsum(this); end
        %         function value = diff(this); end
        %         function value = cumprod(this); end
        %         function value = mode(this); end
        %         function value = median(this); end
        %         function value = cov(this); end
        
        % Similar to sum, but returns a collapsed object of builtin MATLAB type.
        % Also if dim is negative, it sums over all but the specified
        % dimensions.
        function C = collapse(A, dim)
            if ~exist('dim', 'var')
                C = sum(A.vals);
                return;
            end
            
            validateattributes(dim, {'numeric'}, {'real', 'integer'})
            validateattributes(abs(dim), {'numeric'}, {'positive', '<=' A.ndims});
            if ~all(diff(sort(dim)))
                throw(MException('SparseMat:invalidDimension', 'invalid dimension'));
            end
            
            if all(dim < 0)
                dim = setdiff(1:A.ndims, -dim);
            elseif any(dim < 0)
                throw(MException('SparseMat:invalidDimension', 'invalid dimension'));
            end
            
            if sum(ismembc(dim, 1:A.ndims)) == A.ndims
                C = sum(A.vals);
                return
            end
            
            C = A;
            for i = 1:length(dim)
                C = sum(C, dim(i));
            end
            
            C = double(squeeze(C));
        end
        
        function C = norm(A, option)
            import edu.stanford.covert.util.SparseMat;
            
            if ~exist('option', 'var')
                option = 2;
            end
            
            vals = norm(A.vals, option);
            
            if A.numel > A.nnz && ~isnan(vals)
                if option == Inf
                    vals = max(0, vals);
                elseif option == -Inf
                    vals = min(0, vals);
                end
            end
            
            C = A;
            C.subs = ones(numel(vals), 2);
            C.vals = vals;
            C.siz = [1 1];
            C = C.normalize(false, true, false, zeros(1, 2));
        end
    end
    
    %find
    methods
        function [subs, vals] = find(this)
            subs = this.subs;
            vals = this.vals;
        end
        
        %if n>=1, randomly select n zero entries of a sparse tensor
        %if n< 1, randomly select each entry with probability n, return only zero entries
        function subs = randomlySelectZeros(this, n, randStream)
            if n<1
                if exist('randStream','var')
                    inds = randStream.randi(this.numel, [randStream.random('poiss', n*this.numel) 1]);
                else
                    inds = randi(this.numel, [random('poiss', n*this.numel) 1]);
                end
                subs = this.ind2sub(setdiff(inds, this.inds));
            else
                n = min(this.numel-this.nnz, n);
                selectedInds = [];
                while numel(selectedInds)<n
                    if exist('randStream','var')
                        testInds = randStream.randi(this.numel, [n-numel(selectedInds) 1]);
                    else
                        testInds = randi(this.numel, [n-numel(selectedInds) 1]);
                    end
                    selectedInds = [selectedInds; testInds]; %#ok<AGROW>
                end
                
                subs = this.ind2sub(selectedInds);
            end
        end
    end
    
    %iss
    methods
        function this = isnan(this)
            this.subs = this.subs(isnan(this.vals), :);
            this.vals = true(size(this.subs,1),1);
        end
        
        function this = isinf(this)
            this.subs = this.subs(isinf(this.vals), :);
            this.vals = true(size(this.subs,1),1);
        end
        
        function this = isfinite(this)
            this.subs = this.sort_subs([
                this.ind2sub(setdiff(1:this.numel, this.inds));
                this.subs(isfinite(this.vals), :)], ...
                this.siz);
            this.vals = true(size(this.subs,1),1);
        end
    end
    
    %type casting
    methods
        function type = valueClass(this)
            type = class(this.vals);
        end
        
        function val = full(this)
            val = cast(this, valueClass(this));
        end
        
        function val = sparse(this)
            if this.ndims > 2
                throw(MException('SparseMat:invalidDimensions','cannot convert to MATLAB builtin sparse with greater than 2 dimensions'))
            elseif ~isa(this.vals,'double')
                warning('SparseMat:invalidType','converting from SparseMat of type ''%s'' to MATLAB builtin sparse of type ''double''', valueClass(this))
            end
            val = sparse(this.subs(:,1), this.subs(:,2), this.vals,this.siz(1),this.siz(2));
        end
        
        function val = cast(this, newclass)
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.util.SparseMat;
            
            switch newclass
                case 'sptensor',           val = sptensor(this);
                case 'CircularSparseMat',  val = CircularSparseMat(this);
                case 'SparseMat'        ,  val = SparseMat(this);
                case 'logical'
                    val = false(this.siz);
                    val(this.inds) = this.vals;
                otherwise
                    val = zeros(this.siz, newclass);
                    val(this.inds) = this.vals;
            end
        end
        
        function this = valueCast(this, newclass)
            this.vals = cast(this.vals, newclass);
        end
        
        function val = char(this)
            val = char(zeros(this.siz));
            val(this.inds) = this.vals;
        end
        
        function val = logical(this)
            val = false(this.siz);
            val(this.inds) = this.vals;
        end
        
        function val = int8(this)
            val = cast(this, 'int8');
        end
        
        function val = int16(this)
            val = cast(this, 'int16');
        end
        
        function val = int32(this)
            val = cast(this, 'int32');
        end
        
        function val = int64(this)
            val = cast(this, 'int64');
        end
        
        function val = uint8(this)
            val = cast(this, 'uint8');
        end
        
        function val = uint16(this)
            val = cast(this, 'uint16');
        end
        
        function val = uint32(this)
            val = cast(this, 'uint32');
        end
        
        function val = uint64(this)
            val = cast(this, 'uint64');
        end
        
        function val = single(this)
            val = cast(this, 'single');
        end
        
        function val = double(this)
            val = cast(this, 'double');
        end
        
        function val = sptensor(this)
            val = sptensor(this.subs, this.vals, this.siz);
        end
    end
    
    methods (Static)
        function [subs2, I, J] = unique_subs(subs, siz)
            if isempty(subs)
                I = zeros(0, 1);
                J = zeros(0, 1);
                subs2 = subs;
                return;
            end
            
            [inds, order] = sort( (subs-1)*[1 cumprod(siz(1:end-1))]' );
            subs2 = subs(order, :);
            subs2 = subs2([diff(inds) ~= 0; true], :);
            
            if nargout >= 2
                I = order([diff(inds) ~= 0; true]);
            end
            if nargout >= 3
                [~, J] = edu.stanford.covert.util.SparseMat.ismember_subs(subs, subs2, siz);
            end
        end
        
        function tf = isunique_subs(subs, siz)
            if isempty(subs)
                tf = true;
                return;
            end
            tf = all(diff(sort(  (subs-1)*[1 cumprod(siz(1:end-1))]'  )));
        end
        
        function [subs, order] = sort_subs(subs, siz, colSortOrder)
            if nargin < 3
                colSortOrder = numel(siz):-1:1;
            end
            tmp = zeros(size(siz));
            tmp(colSortOrder(end)) = 1;
            tmp(colSortOrder(end-1:-1:1)) = cumprod(siz(colSortOrder(end:-1:2)));
            
            %optimized for memory
            if nargin == 3
                [~, order] = sort((subs-1)*tmp');
                subs = subs(order, :);
                return;
            end
            
            %optimized for run time
            inds = (subs-1)*tmp';            
            if  ~issorted(inds)
                [~, order] = sort(inds);
                subs = subs(order, :);
            elseif nargout > 1
                order = (1:size(subs, 1))';
            end
        end
        
        function [tfs, idxs] = ismember_subs(subs1, subs2, siz)
            [tfs, idxs] = ismember(  (subs1-1)*[1 cumprod(siz(1:end-1))]',  (subs2-1)*[1 cumprod(siz(1:end-1))]'  );
        end
        
        function [C,IA,IB] = intersect_subs(subs1, subs2, siz)
            [~,IA,IB] = intersect(  (subs1-1)*[1 cumprod(siz(1:end-1))]',  (subs2-1)*[1 cumprod(siz(1:end-1))]'  );
            C = subs1(IA, :);
        end
        
        function val = issorted_subs(subs, siz)
            val = issorted((subs-1)*[1 cumprod(siz(1:end-1))]');
        end
    end
    
    %mapped to disk
    methods (Static)
        function toDisk(spmat, fileName)
            import edu.stanford.covert.util.SparseMat;
            
            fid = fopen(fileName, 'w+');
            [subs vals] = find(spmat);
            c = 0;
            c = c + 8*fwrite(fid,size(subs,1),'uint64'); % number of elements
            c = c + 8*fwrite(fid,size(subs,2),'uint64'); % number of dimensions
            fwrite(fid,zeros(1,128-c));
            matrix = SparseMat.sort_subs([subs vals],[size(spmat) NaN],[1 2 3]);
            %             matrix = sortrows([subs vals],1);
            fwrite(fid,matrix(:, 1),'uint32');
            fwrite(fid,matrix(:, 2),'uint8');
            fwrite(fid,matrix(:, 3),'uint32');
            fwrite(fid,matrix(:, 4),'single');
            fclose(fid);
        end
        
        function ret = fromDisk(fileName, dim1, dim2, dim3)
            import edu.stanford.covert.util.SparseMat;
            fid = fopen(fileName, 'r');
            out = fread(fid,2,'uint64');
            nElem = out(1);
            nDim = out(2);
            fclose(fid);
            m = memmapfile(fileName,'Offset',128*1,'Format', ...
                {'uint32', [nElem 1], 'c1'
                'uint8', [nElem 1], 'c2'
                'uint32', [nElem 1], 'c3'
                'single', [nElem 1], 'c4'
                }...
                );
            dim1L = dim1(1);
            dim1U = dim1(2);
            dim2L = dim2(1);
            dim2U = dim2(2);
            dim3L = dim3(1);
            dim3U = dim3(2);
            
            subs = [];
            vals = [];
            % TODO: handle trivial case of dim1L == dim1U (is this actually
            % trivial?)
            L1 = SparseMat.binSearch(m, dim1L, 1, nElem, 1, true);
            U1 = SparseMat.binSearch(m, dim1U, 1, nElem, 1, false);
            
            absLowerInd_d1 = L1.indL1;
            absUpperInd_d1 = U1.indU1;
            
            dim1_lb = absLowerInd_d1;
            while (dim1_lb <= absUpperInd_d1)
                tmp1 = SparseMat.binSearch(m, m.Data(1).c1(dim1_lb), dim1_lb, absUpperInd_d1, 1, false);
                dim1_ub = tmp1.indU1;
                
                % ==========
                % Second Dimension (second column)
                % TODO: handle trivial case of dim1_lb == dim1_ub
                L2 = SparseMat.binSearch(m, dim2L, dim1_lb, dim1_ub, 2, true);
                U2 = SparseMat.binSearch(m, dim2U, dim1_lb, dim1_ub, 2, false);
                
                absLowerInd_d2 = L2.indL2;
                absUpperInd_d2 = U2.indU2;
                
                dim2_lb = absLowerInd_d2;
                while(dim2_lb <= absUpperInd_d2)
                    tmp2 = SparseMat.binSearch(m, m.Data(1).c2(dim2_lb), dim2_lb, absUpperInd_d2, 2, false);
                    dim2_ub = tmp2.indU2;
                    
                    % =====
                    % Third Dimension (third column)
                    % TODO: handle trivial case of dim2_lb == dim2_ub
                    L3 = SparseMat.binSearch(m, dim3L, dim2_lb, dim2_ub, 3, true);
                    U3 = SparseMat.binSearch(m, dim3U, dim2_lb, dim2_ub, 3, false);
                    
                    absLowerInd_d3 = L3.indL3;
                    absUpperInd_d3 = U3.indU3;
                    
                    dim3_lb = absLowerInd_d3;
                    while(dim3_lb <= absUpperInd_d3)
                        tmp3 = SparseMat.binSearch(m, m.Data(1).c3(dim3_lb), dim3_lb, absUpperInd_d3, 3, false);
                        dim3_ub = tmp3.indU3;
                        %%% Do real stuff here!%%%
                        if(tmp3.found)
                            subs = [subs; ...
                                double(m.Data(1).c1(dim3_lb)) - dim1L + 1, ...
                                double(m.Data(1).c2(dim3_lb)) - dim2L + 1, ...
                                double(m.Data(1).c3(dim3_lb)) - dim3L + 1]; %#ok<AGROW>
                            vals = [vals; double(m.Data(1).c4(dim3_lb))]; %#ok<AGROW>
                        end
                        dim3_lb = dim3_ub + 1;
                    end
                    
                    % =====
                    dim2_lb = dim2_ub + 1;
                end
                
                % ==========
                dim1_lb = dim1_ub + 1;
            end
%             ret = struct;
%             
%             ret.subs = subs;
%             ret.vals = vals;
            ret = SparseMat(subs, vals, [dim1U - dim1L + 1 dim2U - dim2L + 1 dim3U - dim3L + 1]);
        end
        
        function ret = binSearch(m, value, low, high, col, lowerend)
            plow = low;
            phigh = high;
            pguess = floor(mean([plow,phigh]));
            colVec = m.Data(1).(sprintf('c%d',col));
            if(lowerend)
                prev = max(pguess-1,1);
                while phigh > plow && (colVec(pguess) ~= value || colVec(prev) == value) 
                    if(colVec(pguess) > value || colVec(prev) == value)
                        phigh = pguess - 1;
                    else
                        plow = pguess + 1;
                    end
                    pguess = floor(mean([plow,phigh]));
                    pguess = max(pguess, low);
                    pguess = min(pguess, high);
                    prev = max(pguess-1, 1);
                end
            else
                next = min(pguess+1,high);
                while phigh > plow && (colVec(pguess) ~= value || colVec(next) == value)
                    if(colVec(pguess) < value || colVec(next) == value)
                        plow = pguess + 1;
                    else
                        phigh = pguess - 1;
                    end
                    pguess = floor(mean([plow,phigh]));
                    pguess = max(pguess, low);
                    pguess = min(pguess, high);
                    next = min(pguess+1, high);
                end
            end
           
            ret = struct;
            if colVec(pguess) == value
%             if phigh > plow
                ret.found = true;
            else
                ret.found = false;
            end
            if(lowerend)
                c = 'L';
            else
                c = 'U';
            end
            ret.(sprintf('ind%c%d',c,col)) = pguess;
        end
    end
    
    %display
    methods
        function disp(this)
            this.display();
        end
        
        function display(this)
            fprintf('%d',this.siz(1));
            fprintf('x%d',this.siz(2:end));
            fprintf(' %s of type %s with %d non-zeros\n', class(this), valueClass(this), this.nnz);
            for i=1:this.nnz
                fprintf('(%d',this.subs(i,1));
                fprintf(',%d',this.subs(i,2:end));
                fprintf(')\t\t%.4f\n',this.vals(i));
            end
        end
    end
end