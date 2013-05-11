function B = xrepmat(A,varargin)
%XREPMAT Replicate and tile an array.
%
%   XREPMAT(A,M,N) and XREPMAT(A,[M N]) replicates and tiles the matrix
%   A to produce a M-by-N block matrix.
%
%   XREPMAT(A,(M,N,P,...) and XREPMAT(A,[M N P ...]) tiles the array A
%   to produce a M-by-N-by-P-by-... block array.  A can be N-D.
%
%   XREPMAT(A,M,N) when A is a scalar is commonly used to produce an
%   M-by-N matrix filled with A's value.  This can be much faster than
%   A*ONES(M,N) when M and/or N are large.
%
%   Example:
%       xrepmat(magic(2),2,3)
%       xrepmat(NaN,2,3)
%
%   See also REPMAT, MESHGRID, NDGRID.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:50:27 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

if nargin < 2
   error( 'Not enough input arguments.' );
elseif nargin == 2
   if length(varargin{1}) == 1          % XREPMAT(A,M)
      Rsiz = [ varargin{1} varargin{1} ];
   else                                 % XREPMAT(A,[M N P ...])
      Rsiz = varargin{1};
   end
else                                    % XREPMAT(A,[M N P ...])
   Rsiz = [ varargin{:} ];
end

if length(A) == 1
   nelems = prod(Rsiz);
   if nelems > 0
      % Since B doesn't exist, the first statement creates a B with the
      % right size and type.  Then use scalar expansion to fill the
      % array.  Finally reshape to the specified size.
      B(nelems) = A;
      B(:) = A;
      B = reshape(B,Rsiz);
   else
      B = A(ones(Rsiz));
   end
else
   Asiz = size(A);
   Adim = length(Asiz);
   Rdim = length(Rsiz);
   if ( Adim == 2 ) & ( Rdim == 2 )
      mind = (1:Asiz(1))';
      nind = (1:Asiz(2))';
      mind = mind(:,ones(1,Rsiz(1)));
      nind = nind(:,ones(1,Rsiz(2)));
      B = A(mind,nind);
   else
      Bdim = max(Adim,Rdim);
      Asiz = [ Asiz ones(1,Bdim-Adim) ];
      Rsiz = [ Rsiz ones(1,Bdim-Rdim) ];

      % This might be `clever', but this is usually slower than the
      % for-loop below.
%      v = reshape( reshape( 1:2*Bdim, Bdim, 2 )', 1, 2*Bdim );
%      A = A(:);
%      B = reshape( A(:,ones(1,prod(Rsiz))), [ Asiz Rsiz ] );
%      B = reshape( permute( B, v ), Asiz.*Rsiz );

      for i = 1:Bdim
         ind = (1:Asiz(i))';
         subs{i} = ind(:,ones(1,Rsiz(i)));
      end
      B = A(subs{:});

   end
end
