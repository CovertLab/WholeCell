function y = asvec(x, dim, order)
%ASVEC  Return elements of an array as a vector.
%
%   ASCOL(X) returns the elements of X as a column vector, running through
%   the elements of X in column major order.
%
%   ASVEC(X, DIM) returns the elements of X as a vector along dimension DIM.
%   For instance, ASVEC(X, 1) returns the elements of X as a column vector
%   and ASVEC(X, 2) returns the elements as a row vector.
%
%   ASVEC(X, DIM, ORDER) will run through the dimensions of X in the order
%   specified in ORDER.  Unspecified dimensions will be run through in
%   increasing order after the specified dimensions.  For instance, if X is
%   4-D, then ASVEC(X, DIM, [2 3]) = ASVEC(X, DIM, [2 3 1 4]).  Row major
%   order is obtained with ASVEC(X, DIM, 2).
%
%   ASVEC(X, DIM) is equivalent to ASVEC(X, DIM, 1:NDIMS(X)).
%
%   See also ASCOL, ASROW, PERMUTE, IPERMUTE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-07-07 22:57:10 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;
   error(nargchk(1, 3, nargsin));

   % if dimension order is specified, then reorder the dimension of `x'
   if nargsin == 3
      nd = max(ndims(x), max(order));   % largest dimension to deal with
      rest = ones(1, nd);               % one `1' for each dimension
      rest(order) = 0;                  % remove dimension in `order'
      pv = [order find(rest)];          % dimension permutation vector
      x = permute(x, pv);               % permute dimensions in `x'
   end

   % default dimension for output vector is `1' (i.e., column vector)
   if nargsin < 2
      dim = 1;
   end

   sy = ones(1, max(2, dim));           % initialize size vector for `y'
   sy(dim) = prod(size(x));             % all elements along dimension `dim'
   y = reshape(x, sy);                  % convert `x' into a vector
