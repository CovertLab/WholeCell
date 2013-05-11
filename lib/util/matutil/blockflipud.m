function y = blockflipud(x, k)
%BLOCKFLIPUD Flip each block in an array in up/down direction.
%
%   BLOCKFLIPUD(X, K) returns X with the rows flipped in up/down
%   direction in blockes of length K.  If X is M*K-by-N-by-P-..., then
%   each block will have size K-by-N-by-P-...
%
%   BLOCKFLIPUD(X) returns X flipped in up/down direction.  It is
%   equivalent to BLOCKFLIPUD(X, SIZE(X,1)).
%
%   See also FLIPUD.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:36 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 2, nargin));

   s = size(x);
   if (nargin == 1) | (k == s(1))
      % Number of rows is not specified, or specified row number matches
      % number of rows in X, so flip the whole thing upside-down.
      y = reshape(x(s(1):-1:1,:), s);
   else
      if k == 1
         % Flipping along a singleton dimension changes nothing.
         y = x;
      else
         if rem(s(1), k)
            error('K doesn''t evenly divide the number of rows in X.');
         end
         % Flip blocks upside-down, K rows at a time.
         y = reshape(x, [ k prod(s)/k ]);
         y = reshape(y(k:-1:1,:), s);
      end
   end
