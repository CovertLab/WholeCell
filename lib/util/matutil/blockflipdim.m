function y = blockflipdim(x, dim, k)
%BLOCKFLIPDIM Flip each block in an array along specified dimension.
%
%   BLOCKFLIPDIM(X, DIM, K) returns X with dimension DIM flipped in
%   blocks of length K.  If S is the size vector of X, then each block
%   will have size [ S(1) S(2) ... S(DIM)/K ... ].
%
%   BLOCKFLIPDIM(X, DIM) returns X flipped along dimension DIM.  It is
%   equivalent to BLOCKFLIPDIM(X, DIM, SIZE(X,DIM))
%
%   BLOCKFLIPDIM(X, 1, K) is equilvalent to BLOCKFLIPUD(X, K), and
%   BLOCKFLIPDIM(X, 2, K) is equilvalent to BLOCKFLIPLR(X, K).
%
%   See also FLIPDIM, BLOCKFLIPUD, BLOCKFLIPLR.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:40 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(2, 3, nargin));

   s    = size(x);
   sdim = size(x,dim);          % S(DIM) fails if DIM > NDIMS(X)

   if nargin == 2
      k = sdim;
   else
      if rem(sdim, k)
         error('K doesn''t evenly divide the length along dimension DIM.');
      end
   end

   if k == 1

      % Flipping along a singleton dimension changes nothing.
      y = x;

   else

      % Dimension permutation vector.
      p = 1:ndims(x);
      p([ 1 dim ]) = [ dim 1 ];

      % Move DIMth dimension up front.
      y = permute(x, p);

      % Reshape if K differs from length along dimension DIM.
      if k ~= sdim
         y = reshape(y, [ k prod(s)/k ]);
      end

      % Flip blocks, reshape back and restore dimension order.
      y = permute(reshape(y(k:-1:1,:), s(p)), p);

   end
