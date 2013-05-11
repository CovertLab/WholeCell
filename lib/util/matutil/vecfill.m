function y = vecfill(x, val)
%VECFILL Fill last non-zero value into elements with zero value.
%
%   Y = VECFILL(X) fills the last non-zero value of the elements of X into
%   the following elements of X with zero value.  For instance,
%   VECFILL([ 0 4 0 0 0 3 0 2 0 0 ]) returns [ 0 4 4 4 4 3 3 2 2 2].
%
%   VECFILL(X, VAL) does the filling based on the value VAL rather than
%   zero.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:50:31 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % check number of input arguments
   error(nargchk(1, 2, nargin));

   % assign default value to VAL if not given otherwise check it
   if nargin < 2
      val = 0;
   end
   if any(size(val) ~= 1)
      error('VAL must be a scalar.');
   end

   % find elements different from VAL
   if isnan(val)
      idx = ~isnan(x);
   else
      idx = x ~= val;
   end

   % now fill the values
   idx(1) = 1;                  % required if x(1) = val
   fillvals = x(idx);
   y = fillvals(cumsum(idx));
