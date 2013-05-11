function y = acosplus1(x)
%ACOSPLUS1 Inverse cosine of argument plus one.
%
%   ACOSPLUS1(X) is ACOS(X+1) calculated in a way that is numerically better
%   when X is close to one.
%
%   This illustrates the difference
%
%      x = -1e-15*(0:0.01:1);
%      plot( x, acos(x+1), 'b-', x, acosplus1(x), 'r-' );
%
%   See also COSMINUS1.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 14:52:43 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(1, 1, nargin));

   y = 2 * asin(sqrt(-(x / 2)));
