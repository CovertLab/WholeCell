function d = distmat2( x, y )
%DISTMAT2 Distance matrix (two nested for-loops).
%
%   D = DISTMAT2(X, Y) returns the distance matrix with all distances
%   between the points represented by the rows of X and Y.
%
%   DISTMAT1(X) is equivalent to DISTMAT2(X, X), but the former computes the
%   distance matrix faster.
%
%   Distance is Euclidean.
%
%   The calculation is done with two for-loops.
%
%   See also DISTMAT0, DISTMAT1, DISTMAT3.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:24 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 2, nargin));

   if nargin == 1               % DISTMAT2(X)

      if ndims(x) ~= 2
         error('Input must be a matrix.');
      end

      m = size(x, 1);
      d = zeros(m, m);          % initialise output matrix

      for i = 1:m-1
         for j = i+1:m
            diff = x(i,:) - x(j,:);        % difference
            d(i,j) = sqrt(diff*diff');     % distance
            d(j,i) = d(i,j);
         end
      end

   else                            % DISTMAT2(X, Y)

      if ndims(x) ~= 2 | ndims(y) ~= 2
         error('Input must be two matrices.');
      end

      [mx, nx] = size(x);
      [my, ny] = size(y);

      if nx ~= ny
         error('Both matrices must have the same number of columns.');
      end

      m = mx;                      % number of rows in distance matrix
      n = my;                      % number of columns in distance matrix
      d = zeros(m, n);             % initialise output matrix

      for i = 1:m
         for j = 1:n
            diff   = x(i,:) - y(j,:);      % difference
            d(i,j) = sqrt(diff*diff');     % distance
         end
      end

   end
