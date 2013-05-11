function d = distmat1(x, y)
%DISTMAT1 Distance matrix (one for-loop).
%
%   D = DISTMAT1(X, Y) returns the distance matrix with all distances
%   between the points represented by the rows of X and Y.
%
%   DISTMAT1(X) is equivalent to DISTMAT1(X, X), but the former computes the
%   distance matrix faster.
%
%   Distance is Euclidean.
%
%   The calculation is done with one for-loop.
%
%   See also DISTMAT0, DISTMAT2, DISTMAT3.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:28 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 2, nargin));

   if nargin == 1                       % DISTMAT1(X)

      if ndims(x) ~= 2
         error('Input must be a matrix.');
      end

      m = size(x, 1);
      d = zeros(m, m);                  % initialise output matrix

      for i = 1:m-1
         xi = x(i,:);
         diff = xi(ones(m-i,1),:) - x(i+1:m,:);    % difference
         dist = sqrt(sum(abs(diff).^2, 2));        % distance
         d(i+1:m,i) = dist;
         d(i,i+1:m) = dist.';
      end

   else                                 % DISTMAT1(X, Y)

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

      % The for-loop is applied on the shortest dimension, the other is
      % vectorized.
      if m < n
         idx = ones(1, n);
         for i = 1:m
            xi     = x(i,:);
            d(i,:) = sqrt(sum(abs(y - xi(idx,:)).^2, 2)).';
         end
      else
         idx = ones(1, m);
         for j = 1:n
            yj     = y(j,:);
            d(:,j) = sqrt(sum(abs(x - yj(idx,:)).^2, 2));
         end
      end

   end
