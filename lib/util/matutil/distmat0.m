function d = distmat0(x, y)
%DISTMAT0 Distance matrix (no for-loops).
%
%   D = DISTMAT0(X, Y) returns the distance matrix with all distances
%   between the points represented by the rows of X and Y.
%
%   DISTMAT0(X) is equivalent to DISTMAT0(X, X), but the former computes the
%   distance matrix faster.
%
%   Distance is Euclidean.
%
%   The calculation is done with no for-loops.
%
%   See also DISTMAT1, DISTMAT2, DISTMAT3.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:29 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 2, nargin));

   if nargin == 1                       % DISTMAT0(X)

      if ndims(x) ~= 2
         error('Input must be a matrix.');
      end

      m = size(x, 1);
      [i, j] = find(triu(ones(m), 1));  % trick to get indices

      % Calculate the distances on the upper triangular part and then let
      % the lower triangular part be mirror image of the upper triangular
      % part.
      d = zeros(m, m);                  % initialise output matrix
      d( i + m*(j-1) ) = sqrt(sum(abs(x(i,:)-x(j,:)).^2, 2));
      d( j + m*(i-1) ) = d( i + m*(j-1) );

   else                                  % DISTMAT0(X, Y)

      if ndims(x) ~= 2 | ndims(y) ~= 2
         error('Input must be two matrices.');
      end

      [mx, nx] = size(x);
      [my, ny] = size(y);

      if nx ~= ny
         error('Both matrices must have the same number of columns.');
      end

      m = mx;                   % number of rows in distance matrix
      n = my;                   % number of columns in distance matrix
      p = nx;                   % dimension of each point

      % This is fast and reasonably memory efficient (I think).
      x = permute(x, [ 1 3 2 ]);
      y = permute(y, [ 3 1 2 ]);
      d = sqrt(sum((x(:,ones(1,n),:) - y(ones(1,m),:,:)).^2, 3));

   end
