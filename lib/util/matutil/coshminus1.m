function y = coshminus1(x)
%COSHMINUS1 Hyperbolic cosine minus one.
%
%   COSHMINUS1(X) is COSH(X)-1 calculated in a way that is numerically better
%   when X is close zero.
%
%   This illustrates the difference
%
%      x = 5e-8*(-1:0.01:1);
%      plot(x, cosh(x)-1, 'b-', x, coshminus1(x), 'r-');
%
%   See also ACOSHPLUS1.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 14:57:06 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(1, 1, nargin));

   y = 2 * sinh(x / 2).^2;
