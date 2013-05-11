function d = blockdiag(varargin)
%BLOCKDIAG Build a block diagonal matrix of the input arguments.
%
%   BLOCKDIAG(A, B, C, ...) returns a block diagonal matrix which has the
%   matrices A, B, C, ... on the main diagonal.
%
%   BLOCKDIAG(A, B, C, ... , 'sparse') returns the same, but as a sparse
%   matrix.
%
%   See also DIAG, BLKDIAG, HORZCAT, VERTCAT.

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-07-20 12:55:00 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, Inf, nargin))

   % should we return a sparse matrix?
   want_sparse = strcmp(varargin{end}, 'sparse');
   if want_sparse
      varargin = varargin(1:end-1);
   end

   nargsin = length(varargin);
   error(nargchk(1, Inf, nargsin))

   % calculate the size of the output matrix
   rows = 0;            % number of rows in output matrix
   cols = 0;            % number of columns in output matrix
   nnz  = 0;            % number of non-zero elements in output matrix
   for k = 1:nargsin
      [r, c, d] = size(varargin{k});
      if d > 1
         error('Matrices can not have more than 2 dimensions.');
      end
      rows = rows + r;
      cols = cols + c;
      if want_sparse
         nnz = nnz + length(find(varargin{k}));
         varargin{k} = sparse(varargin{k});
      end
   end

   % initialize output matrix
   if want_sparse
      d = sparse([], [], [], rows, cols, nnz);
      if nnz == 0, return, end
   else
      d = zeros(rows, cols);
   end

   % fill the input matrices into the output matrix
   i = 0;
   j = 0;
   for k = 1:nargsin
      [r, c] = size(varargin{k});
      d(i+1:i+r , j+1:j+c) = varargin{k};
      i = i + r;
      j = j + c;
   end
