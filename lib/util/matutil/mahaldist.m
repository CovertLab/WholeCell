function D = mahaldist(X, Y, W)
%MAHALDIST Mahalanobis distance.
%
%   D = MAHALDIST(X, Y, W) computes a distance matrix D where the element
%   D(i,j) is the Mahalanobis distance between the points X(i,:) and Y(j,:)
%   using W as the weight matrix so that
%
%      D(i,j) = (X(i,:) - Y(j,:)) * W * (X(i,:) - Y(j,:))';
%
%   If W is the identity matrix, then D is the squared Euclidean distances.
%
%   MAHALDIST(X, [], W) is the same as MAHALDIST(X, X, W), but the former is
%   computed more efficiently.
%
%   MAHALDIST([], Y, W) is the same as MAHALDIST(Y, Y, W), but the former is
%   computed more efficiently.
%
%   D = MAHALDIST(X, Y) is the same as MAHALDIST(X, MEAN(Y,1), INV(COV(Y))),
%   except that the former is computed more efficiently.  In this case, D(i)
%   is the Mahalanobis distance from the point X(i,:) to the centroid of set
%   Y (i.e., the mean of the rows in Y) with respect to the variance in the
%   set Y (i.e., the variance matrix of the rows in Y).
%
%   MAHALDIST(X) is the same as MAHALDIST(X, X), but the former is computed
%   more efficiently.
%
%   NOTE: The Mahalanobis is sometimes defined as the square root of the
%   values returned by MAHALDIST.
%
%   NOTE: The way this function interprets the input argument has changed
%   since earlier versions.  Specifically, the two-argument call now assumes
%   the arguments are swapped compared to earlier releases.
%
%   When X and Y are real, then MAHAL(X, Y) is the same as MAHALDIST(X, Y)
%   and PDIST(X, 'mahal') returns the lower triangular part of
%   SQRT(MAHALDIST(X, [], INV(COV(X)))).  Both MAHAL and PDIST are in the
%   Statistics Toolbox.
%
%   See also MAHAL (Statistics Toolbox), PDIST (Statistics Toolbox).

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-06-19 12:43:59 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;
   error(nargchk(1, 3, nargsin));

   if nargsin <= 2

      %
      % MAHALDIST(X), MAHALDIST(X, Y)
      %
      % No weight matrix specified, so use the inverse of the unbiased
      % variance matrix.
      %

      if ndims(X) ~= 2
         error('X must be a 2D array.');
      end

      [mx, nx] = size(X);

      if nargsin == 1

         M = sum(X, 1)/mx;              % centroid (mean)
         Xc = X - M(ones(mx,1),:);      % subtract centroid of X
         W = (Xc' * Xc)/(mx - 1);       % variance matrix

      else

         if ndims(Y) ~= 2
            error('Y must be a 2D array.');
         end

         [my, ny] = size(Y);

         if nx ~= ny
            error('X and Y must have the same number of columns.');
         end

         M = sum(Y, 1)/my;              % centroid (mean)
         Yc = Y - M(ones(my,1),:);      % subtract centroid of Y
         W = (Yc' * Yc)/(my - 1);       % variance matrix

         Xc = X - M(ones(mx,1),:);      % subtract centroid of Y
      end

      % The call to REAL here is only to remove ``numerical noise''.
      D = real(sum((Xc / W) .* conj(Xc), 2));   % Mahalanobis distances

   else

      %
      % MAHALDIST(X, Y, W), MAHALDIST(X, [], W), MAHALDIST([], Y, W)
      %

      if ndims(X) ~= 2 | ndims(Y) ~= 2 | ndims(W) ~= 2
         error('X, Y, and W must be 2D arrays.');
      end

      [mx, nx] = size(X);
      [my, ny] = size(Y);
      [mw, nw] = size(W);

      if mw ~= nw
         error('W must be a square matrix.');
      end

      if isempty(X) | isempty(Y)
         if isempty(X)                  % MAHALDIST(X, [], W)
            if ny ~= nw
               error('Y and W must have the same number of columns.');
            end
            m = my;
            X = Y;
         else                           % MAHALDIST([], Y, W)
            if nx ~= nw
               error('X and W must have the same number of columns.');
            end
            m  = mx;
         end

         % The loop below is a semi-vectorized version of the code
         %
         %   D = zeros(m, m);
         %   for j = 1 : m
         %      for i = 1 : m
         %         Dij = X(i,:) - X(j,:);
         %         D(i,j) = Dij * W * Dij';
         %      end
         %   end

         D = zeros(m, m);
         for i = 1 : m-1
            Xi = X(i,:);
            Xi = Xi(ones(m-i,1),:);
            Dc = X(i+1:m,:) - Xi;
            D(i+1:m,i) = real(sum((Dc * W) .* conj(Dc), 2));
            D(i,i+1:m) = D(i+1:m,i).';
         end

      else                              % MAHALDIST([], Y, W)

         if nx ~= ny | nx ~= nw
            error('X, Y, and W must have the same number of columns.');
         end

         D = zeros(mx, my);

         % The loops below are semi-vectorized versions of the code
         %
         %   D = zeros(mx, my);
         %   for j = 1 : my
         %      for i = 1 : mx
         %         Dij = X(i,:) - Y(j,:);
         %         D(i,j) = Dij * W * Dij';
         %      end
         %   end

         % loop over the shorter of the two dimensions
         if mx <= my

            for i = 1 : mx
               Xi = X(i,:);
               Xi = Xi(ones(my,1),:);
               Dc = Xi - Y;
               D(i,:) = reshape(real(sum((Dc * W) .* conj(Dc), 2)), [1 my]);
            end

         else

            for j = 1 : my
               Yj = Y(j,:);
               Yj = Yj(ones(mx,1),:);
               Dc = X - Yj;
               D(:,j) = real(sum((Dc * W) .* conj(Dc), 2));
            end

         end
      end

   end
