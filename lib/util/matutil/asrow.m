function y = asrow(x, varargin)
%ASROW  Return elements of an array as a row vector.
%
%   ASROW(X) returns the elements of X as a row vector, running through the
%   elements X in column major order.
%
%   ASROW(X, ORDER) will run through the dimensions of X in the order
%   specified in ORDER.  Unspecified dimensions will be run through in
%   increasing order after the specified dimensions.  For instance, if X is
%   4-D, then ASROW(X, [2 3]) = ASROW(X, [2 3 1 4]).
%
%   For example, if
%
%      X = [4 6 8   then  ASROW(X)    = [4 5 6 7 8 9]
%           5 7 9]        ASROW(X, 2) = [4 6 8 5 7 9]
%
%   where ASROW(X) = ASROW(X, 1) = ASROW(X, [1 2]) and
%   ASROW(X, 2) = ASROW(X, [2 1]).
%
%   See also ASCOL, ASVEC, RESHAPE, PERMUTE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-07-07 22:57:10 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;
   error(nargchk(1, 2, nargsin));

   y = asvec(x, 2, varargin{:});
