function B = repelm(A,varargin)
%REPELM Replicate elements in an array.
%
%   REPELM(A,M,N) and REPELM(A,[M N]) replicates each element in the matrix
%   A so it becomes an M-by-N matrix.
%
%   REPELM(A,M,N,P,...) and REPELM(A,[M N P ...]) replicates each element in
%   A so it becomes an M-by-N-by-P-by-... array.  A can be N-D.
%
%   REPELM(A,M,N) when A is a matrix is the same as KRON(A,ONES(M,N)), but
%   the former is much faster since it does no multiplications.  The former
%   also works when A is not of class double.
%
%   Example:
%       repelm(magic(2),2,3)
%       repelm(NaN,2,3)
%
%   See also REPMAT, MESHGRID, NDGRID.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:21:10 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

if nargin < 2
   error( 'Not enough input arguments.' );
elseif nargin == 2
   if length(varargin{1}) == 1          % REPELM(A,M)
      siz = [ varargin{1} varargin{1} ];
   else                                 % REPELM(A,[M N P ...])
      siz = varargin{1};
   end
else                                    % REPELM(A,M,N,P,...)
   siz = [ varargin{:} ];
end

if length(A) == 1
   nelems = prod(siz);
   if nelems > 0
      % Since B doesn't exist, the first statement creates a B with
      % the right size and type.  Then use scalar expansion to
      % fill the array.  Finally reshape to the specified size.
      B(nelems) = A;
      B(:) = A;
      B = reshape(B,siz);
   else
      B = A(ones(siz));
   end
elseif ndims(A)==2 & length(siz)==2
   [m,n] = size(A);
   mind = 1:m;
   nind = 1:n;
   mind = mind(ones(1,siz(1)),:);
   nind = nind(ones(1,siz(2)),:);
   B = A(mind,nind);
else

%   Asiz = size(A);
%   Asiz = [Asiz ones(1,length(siz)-length(Asiz))];
%   siz = [siz ones(1,length(Asiz)-length(siz))];
%   for i = 1:length(Asiz)
%      ind = 1:Asiz(i);
%      subs{i} = ind(ones(1,siz(i)),:);
%   end
%   B = A(subs{:});

   Asiz = size(A);
   Adim = length(Asiz);
   Rdim = length(siz);
   Bdim = max(Adim, Rdim);
   Asiz  = [Asiz ones(1,Bdim-Adim)];
   subs = {':'};
   subs = subs(ones(1,Bdim));
   for i = 1:Rdim
      if siz(i) > 1
         ind = 1:Asiz(i);
         subs{i} = ind(ones(1,siz(i)),:);
      end
   end
   B = A(subs{:});

end
