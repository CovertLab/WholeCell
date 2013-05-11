function d = distmat3(x, y)
%DISTMAT3 Distance matrix (three nested for-loops).
%
%   D = DISTMAT3(X, Y) returns the distance matrix with all distances
%   between the points represented by the rows of X and Y.
%
%   DISTMAT3(X) is equivalent to DISTMAT3(X, X), but the former computes the
%   distance matrix faster.
%
%   Distance is Euclidean.
%
%   The calculation is done with three for-loops.
%
%   See also DISTMAT0, DISTMAT1, DISTMAT2.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:22 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 2, nargin));

   if nargin == 1               % DISTMAT3(X)

      if ndims(x) ~= 2
         error('Input must be a matrix.');
      end

      [m, p] = size(x);
      d = zeros(m, m);          % initialise output matrix

      for i = 1:m-1
         for j = i+1:m
            ssq = 0;            % sum of squares
            for k = 1:p
               ssq = ssq + abs(x(i,k) - x(j,k))^2;
            end
            d(i,j) = sqrt(ssq);
            d(j,i) = d(i,j);
         end
      end

   else                         % DISTMAT3(X)

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
      d = zeros(m, n);          % initialise output matrix

      for i = 1:m
         for j = 1:n
            ssq = 0;         % Sum of squares.
            for k = 1:p
               ssq = ssq + abs(x(i,k) - y(j,k))^2;
            end
            d(i,j) = sqrt(ssq);
         end
      end

   end
