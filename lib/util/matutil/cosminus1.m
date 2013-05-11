function y = cosminus1(x)
%COSMINUS1 Cosine minus one.
%
%   COSMINUS1(X) is COS(X)-1 calculated in a way that is numerically better
%   when X is close zero, or more generally, when X is close to 2*pi*K where K
%   is an integer.
%
%   This illustrates the difference
%
%      x = 4e-8*(-1:0.01:1);
%      plot( x, cos(x)-1, 'b-', x, cosminus1(x), 'r-' );
%
%   See also ACOSPLUS1.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 15:17:49 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(1, 1, nargin));

   y = -2 * sin(x / 2).^2;
