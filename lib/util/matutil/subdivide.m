function y = subdivide(x, p, q)
%SUBDIVIDE Subdivide a matrix into several submatrices.
%
%   Y = SUBDIVIDE(X, P, Q) where X is an M-by-N matrix and P and Q are
%   integers that are factors of P and Q respectively, returns a
%   P-by-Q-by-M/P-by-N/Q matrix so that Y(:,:,I,J) is the (I,J)th
%   submatrix of X.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:50:33 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(3, 3, nargin));
   if ndims(x) ~= 2
      error('First argument must be a 2D matrix.');
   end

   [ m, n ] = size(x);
   if rem(m, p) | rem(n, q)
      error(sprintf([ 'Can''t create %d-by-%d submatrices from' ...
                      ' a %d-by-%d matrix.' ], p, q, m, n));
   end

   y = reshape(x, [ p m/p q n/q ]);
   y = permute(y, [ 1 3 2 4 ]);
