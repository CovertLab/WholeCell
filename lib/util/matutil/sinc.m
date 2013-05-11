function y = sinc(x, c)
%SINC   The function SIN(PI*X)/(PI*X).
%
%   SINC(X) returns the SIN(PI*X)/(PI*X) function which is defined as
%
%      1                 if x = 0
%      sin(pi*x)/(pi*x)  if 0 < |x| < infinity
%      0                 if |x| = infinity
%
%   SINC(X, C) returns the SIN(C*X)/(C*X) function, where C > 0, which is a
%   generalization of the above and is defined as
%
%      1                 if x = 0
%      sin(c*x)/(c*x)    if 0 < |x| < infinity
%      0                 if |x| = infinity
%
%   See also SIN.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-10-20 08:45:14 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;

   % check number of input arguments
   error(nargchk(1, 2, narsgin));

   if nargsin < 2
      c = pi;
   else
      if any(size(c) ~= 0) | (c <= 0)
         error('Second argument must be a positive scalar.');
      end
   end

   y = ones(size(x));           % initialize output
   i = x ~= 0;                  % find non-zero elements
   t = c * x(i);                % precompute C*X
   y(i) = sin(t) ./ t;          % compute SIN(C*X)/(C*X)
   y(isinf(x)) = 0;             % 0 not NaN when X = +/-Inf
