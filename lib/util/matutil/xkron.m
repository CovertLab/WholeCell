function Z = xkron(X,Y)
%XKRON Kronecker tensor product.
%
%   XKRON(X,Y), where X and Y are matrices, is the Kronecker tensor product
%   of X and Y.  The result is a large matrix formed by taking all possible
%   products between the elements of X and those of Y.
%
%   For example, if X is 2-by-3, and Y is any N-dimensional array, then
%   KRON(X,Y) is
%
%     [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
%       X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]
%
%   If either X or Y is sparse, and neither X nor Y has dimension
%   greater than 2, then only nonzero elements are multiplied in the
%   computation, and the result is sparse.
%
%   XKRON works also when X and/or Y are N-dimensional arrays.  There does
%   not seem to be a unique definition of an N-dimensional Kronecker tensor
%   product, so I simply picked the one I found most logical.
%
%   If X is 3-by-4-by-2, and Y is any N-dimensional array, then
%
%     cat( 3, [ X(1,1,1)*Y  X(1,2,1)*Y  X(1,3,1)*Y  X(1,4,1)*Y ;  ...
%               X(2,1,1)*Y  X(2,2,1)*Y  X(2,3,1)*Y  X(2,4,1)*Y ;  ...
%               X(3,1,1)*Y  X(3,2,1)*Y  X(3,3,1)*Y  X(3,4,1)*Y ], ...
%             [ X(1,1,2)*Y  X(1,2,2)*Y  X(1,3,2)*Y  X(1,4,2)*Y ;  ...
%               X(2,1,2)*Y  X(2,2,2)*Y  X(2,3,2)*Y  X(2,4,2)*Y ;  ...
%               X(3,1,2)*Y  X(3,2,2)*Y  X(3,3,2)*Y  X(3,4,2)*Y ] )

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:50:29 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   xs = size(X);
   ys = size(Y);
   xd = length(xs);
   yd = length(ys);

   if ( issparse(X) | issparse(Y) ) & ( xd <= 2 ) & ( yd <= 2 )

      % When at least one input is sparse, and neither argument has
      % dimension greater than 2, the result is sparse.

      mx = xs(1); nx = xs(2);
      my = ys(1); ny = ys(2);

      [ix, jx, sx] = find(X);
      [iy, jy, sy] = find(Y);
      ix = ix(:); jx = jx(:); sx = sx(:);
      iy = iy(:); jy = jy(:); sy = sy(:);
      kx = ones(size(sx));
      ky = ones(size(sy));
      t = my*(ix-1)';
      ik = t(ky,:)+iy(:,kx);
      t = ny*(jx-1)';
      jk = t(ky,:)+jy(:,kx);
      Z = sparse(ik, jk, sy*sx.', mx*my, nx*ny);

   else

      dz = max(xd, yd);
      xs = [ xs ones(1,dz-xd) ];
      ys = [ ys ones(1,dz-yd) ];

      v = reshape( reshape( 1:2*dz, dz, 2 )', 1, 2*dz );
      Z = reshape( Y(:)*X(:).', [ ys xs ] );
      Z = reshape( permute( Z, v ), xs.*ys );

   end
