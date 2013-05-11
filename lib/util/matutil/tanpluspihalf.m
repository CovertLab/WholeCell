function y = tanpluspihalf(x)
%TANPLUSPIHALF Tangent of argument plus pi/2.
%
%   TANPLUSPIHALF(X) is TAN(X+PI/2) calculated in a way that is numerically
%   better when X is close to -PI/2.
%
%   See also ATANMINUSPIHALF.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 15:04:23 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(1, 1, nargin));

   y = zeros(size(x));

   c = -pi / 2;

   k = x < c;
   y(k) = tan(x(k) + pi/2);

   k = x > c;
   y(k) = -1 ./ tan(x(k));
