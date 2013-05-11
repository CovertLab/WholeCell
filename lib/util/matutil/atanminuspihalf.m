function y = atanminuspihalf(x)
%ATANMINUSPIHALF Inverse tangent minus pi/2.
%
%   ATANMINUSPIHALF(X) is ATAN(X)-PI/2 calculated in a way that is numerically
%   better when X is large and positive.
%
%   See also TANPLUSPIHALF.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-13 14:55:25 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(1, 1, nargin));

   y = zeros(size(x));

   k = x <= 0;
   y(k) = atan(x(k)) - (pi / 2);

   k = x > 0;
   y(k) = atan(-1 ./ x(k));
