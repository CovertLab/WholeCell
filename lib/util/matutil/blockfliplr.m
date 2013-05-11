function y = blockfliplr(x, k)
%BLOCKFLIPLR Flip each block in an array in left/right direction.
%
%   BLOCKFLIPLR(X, K) returns X with the columns flipped in left/right
%   direction in blockes of length K.  If X is M-by-N*K-by-P-..., then
%   each block will have size M-by-K-by-P-...
%
%   BLOCKFLIPLR(X) returns X flipped in left/right direction.  It is
%   equivalent to BLOCKFLIPLR(X, SIZE(X,2)).
%
%   See also FLIPLR.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:38 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 2, nargin));

   s = size(x);
   if (nargin == 1) | (k == s(1))
      % Number of rows is not specified, or specified row number matches
      % number of rows in X, so flip the whole thing upside-down.
      y = reshape(x(:,s(2):-1:1,:), s)
   else
      if k == 1
         % Flipping along a singleton dimension changes nothing.
         y = x;
      else
         if rem(s(2), k)
            error('K doesn''t evenly divide the number of columns in X.');
         end
         y = reshape(x, [ s(1) k prod(s(2:end))/k ]);
         y = reshape(y(:,k:-1:1,:), s);
      end
   end
