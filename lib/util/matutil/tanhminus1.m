function y = tanhminus1(x)
%TANHMINUS1 Hyperbolic tangent minus one.
%
%   TANHMINUS1(X) is TANH(X)-1 calculated in a way that is numerically better
%   when X is large and positive.
%
%   This example illustrates the difference
%
%      x = 17.5:0.01:20.5;
%      plot(x, tanh(x)-1, 'y-', x, tanhminus1(x), 'r-');

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 15:04:13 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(1, 1, nargin));

   y = -2 ./ (1 + exp(2 * x));
