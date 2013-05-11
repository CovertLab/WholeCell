function y = acoshplus1(x)
%ACOSHPLUS1 Inverse hyperbolic cosine of argument plus one.
%
%   ACOSHPLUS1(X) is ACOSH(X+1) calculated in a way that is numerically better
%   when X is close to one.
%
%   This illustrates the difference
%
%      x = 2e-15*(0:0.01:1);
%      plot(x, acosh(x+1), 'b-', x, acoshplus1(x), 'r-');
%
%   See also COSHMINUS1.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 14:54:17 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(1, 1, nargin));

   y = 2 * asinh(sqrt(x / 2));
