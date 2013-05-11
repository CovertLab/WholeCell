% CircularSparseMat
%   Multidimensional circular sparse matrix. Builds on SparseMat adding circular
%   wrapping of indices
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/18/2010
classdef CircularSparseMat < edu.stanford.covert.util.SparseMat
    properties (SetAccess=protected)
        circularDims
    end
    
    %constructor
    methods
        %five ways to instatiate:
        %(1) x = CircularSparseMat()
        %    Creates an empty sparse matrix of size 0x2 with no dimensions
        %    wrapping
        %
        %(2) x = CircularSparseMat(mat)
        %    Casts matrix mat to sparse matrix with no dimensions wrapping
        %
        %(3) x = CircularSparseMat(mat, circularDims)
        %    Casts matrix mat to sparse matrix with selected dimensions wrapping
        %
        %(4) x = CircularSparseMat(subs, vals, sz)
        %    Creates sparse matrix with size sz and non-zeros with subscripts
        %    subs and values vals with no dimensions wrapping
        %
        %(5) x = CircularSparseMat(subs, vals, sz, circularDims)
        %    Creates sparse matrix with size sz and non-zeros with subscripts
        %    subs and values vals with selected dimensions wrapping
        function this = CircularSparseMat(varargin)
            import edu.stanford.covert.util.CircularSparseMat;
            
            switch nargin
                case 4, supervarargin = varargin(1:3);
                    circDims = varargin{end};
                    
                    subs = varargin{1};
                    siz  = varargin{3};
                    if ~isempty(subs)
                        for i = 1:numel(circDims)
                            subs(:, circDims(i)) = mod(subs(:, circDims(i)) - 1, siz(circDims(i))) + 1;
                        end
                        supervarargin{1} = subs;
                    end
                case 0,
                case 1, supervarargin = varargin(1);
                    circDims = varargin{1}.circularDims;
                case 2, supervarargin = varargin(1);
                    circDims = varargin{end};
                case 3, supervarargin = varargin(1:3);
                    circDims = zeros(1,0);
                otherwise
                    throw(MException('CircularSparseMat:error','no constructor matches the calling signature'));
            end
            
            this = this@edu.stanford.covert.util.SparseMat(supervarargin{:});
            
            validateattributes(circDims, {'numeric'}, {'positive', 'real', 'integer', '<=', this.ndims});
            if numel(circDims) > 1 && ~all(diff(sort(circDims)))
                throw(MException('CircularSparseMat:invalidDimensions', 'circularDims must be positive real intergers in the range 1..ndims(this)'));
            end
            
            if numel(circDims) == 1
                this.circularDims = circDims;
            else
                this.circularDims = sort(circDims);
            end
        end
        
        function this = normalize(this, reSize, reFind, reSort, oldDims)
            this = this.normalize@edu.stanford.covert.util.SparseMat(reSize, reFind, reSort, oldDims);
            
            %inlined is ismember optimized for small number of circular dimensions
            if isscalar(this.circularDims)
                this.circularDims = find(oldDims == this.circularDims);
            else
                tfs = false(size(oldDims));
                for i = 1:numel(oldDims)
                    tfs(i) = any(this.circularDims == oldDims(i));
                end
                this.circularDims = find(tfs);
            end
        end
        
        function value = isDimCircular(this, dim)
            value = any(this.circularDims == dim);
        end
        
        function tf = isequal(A, B)
            tf = isequal@edu.stanford.covert.util.SparseMat(A, B) && ...
                isequal(A.circularDims, B.circularDims);
        end
        
        function C = tprod(A, B, Aid, Bid)
            Aod = setdiff(1:A.ndims, Aid);
            Bod = setdiff(1:B.ndims, Bid);
            
            circDims = [];
            if isa('CircularSparseMat',A)
                [tfs, idxs] = ismember(A.circularDims, Aod);
                circDims = idxs(tfs);
            end
            if isa('CircularSparseMat',A)
                [tfs, idxs] = ismember(B.circularDims, Bod);
                circDims = [circDims idxs(tfs)+numel(Aod)];
            end
            
            C = tprod@edu.stanford.covert.util.SparseMat(A, B);
            C.circularDims = circDims;
        end
        
        function this = subsasgn(this, s, rhs)
            %Case 1: Dot reference to properties
            if strcmp(s.type,'.') || strcmp(s.type, '{}')
                this = this.subsasgn@edu.stanford.covert.util.SparseMat(s, rhs);
                return;
            end
            
            if numel(s.subs)==1
                %Case 2: Subscripts of elements
                subs = s.subs{1};
                
                if size(subs,2) ~= this.ndims
                    throw(MException('SparseMat:invalidDimensions','Subscripts must have number of columns equal to dimensions of matrix.'));
                end
                for i=1:numel(this.circularDims)
                    subs(:,this.circularDims(i)) = mod(subs(:,this.circularDims(i))-1, this.siz(this.circularDims(i)))+1;
                end
                
                s.subs{1}=subs;
            else
                %Case 3: Subscripts ranges
                for i=1:numel(this.circularDims)
                    if this.circularDims(i) > numel(s.subs)
                        break;
                    end
                    
                    subs = s.subs{this.circularDims(i)};
                    
                    if isnumeric(subs)
                        subs = mod(subs-1, this.siz(this.circularDims(i)))+1;
                    end
                    
                    s.subs{this.circularDims(i)} = subs;
                end
            end
            
            this = this.subsasgn@edu.stanford.covert.util.SparseMat(s, rhs);
        end
        
        function value = subsref(this, s)
            %Case 1: Dot reference to properties
            if strcmp(s.type,'.') || strcmp(s.type, '{}')
                value = this.subsref@edu.stanford.covert.util.SparseMat(s);
                return;
            end
            
            if numel(s.subs)==1
                %Case 2: Subscripts of elements
                subs = s.subs{1};
                
                if size(subs,2) ~= this.ndims
                    throw(MException('SparseMat:invalidDimensions','Subscripts must have number of columns equal to dimensions of matrix.'));
                end
                for i=1:numel(this.circularDims)
                    subs(:,this.circularDims(i)) = mod(subs(:,this.circularDims(i))-1, this.siz(this.circularDims(i)))+1;
                end
                
                s.subs{1}=subs;
            else
                %Case 3: Subscripts ranges
                for i=1:numel(this.circularDims)
                    if this.circularDims(i) > numel(s.subs)
                        break;
                    end
                    
                    subs = s.subs{this.circularDims(i)};
                    
                    if isnumeric(subs)
                        subs = mod(subs-1, this.siz(this.circularDims(i)))+1;
                    end
                    
                    s.subs{this.circularDims(i)} = subs;
                end
            end
            
            value = this.subsref@edu.stanford.covert.util.SparseMat(s);
        end
        
        function display(this)
            if isempty(this.circularDims)
                circDims = 'no circular dimensions';
            else
                circDims = sprintf('circular dimensions {%d', this.circularDims(1));
                if numel(this.circularDims)>1
                    circDims = [circDims sprintf(',%d', this.circularDims(2:end))];
                end
                circDims = [circDims '}'];
            end
            
            fprintf('%d',this.siz(1));
            fprintf('x%d',this.siz(2:end));
            fprintf(' %s of type %s with %s and %d non-zeros\n', ...
                class(this), valueClass(this), circDims, this.nnz);
            for i=1:this.nnz
                fprintf('(%d',this.subs(i,1));
                fprintf(',%d',this.subs(i,2:end));
                fprintf(')\t\t%.4f\n',this.vals(i));
            end
        end
    end
end