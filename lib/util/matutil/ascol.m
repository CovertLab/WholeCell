function y = ascol(x, varargin)
%ASCOL  Return elements of an array as a column vector.
%
%   ASCOL(X) returns the elements of X as a column vector, running through
%   the elements X in column major order.
%
%   ASCOL(X, ORDER) will run through the dimensions of X in the order
%   specified in ORDER.  Unspecified dimensions will be run through in
%   increasing order after the specified dimensions.  For instance, if X is
%   4-D, then ASCOL(X, [2 3]) = ASCOL(X, [2 3 1 4]).
%
%   For example, if
%
%      X = [4 6 8   then  ASCOL(X) = [4   and  ASCOL(X, 2) = [4
%           5 7 9]                    5                       6
%                                     6                       8
%                                     7                       5
%                                     8                       7
%                                     9]                      9]
%
%   where ASCOL(X) = ASCOL(X, 1) = ASCOL(X, [1 2]) and
%   ASCOL(X, 2) = ASCOL(X, [2 1]).
%
%   See also ASROW, ASVEC, RESHAPE, PERMUTE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-07-07 22:57:30 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;
   error(nargchk(1, 2, nargsin));

   y = asvec(x, 1, varargin{:});
